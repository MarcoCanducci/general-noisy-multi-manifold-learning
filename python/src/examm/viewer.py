"""Standalone interactive viewer (one HTML file) for one or several fitted hierarchies.

    examm.viewer.export_viewer("viewer.html", X, {"run A": H_a, "run B": H_b}, names=feature_names)

The page needs no server: open it in a browser. It shows the points in 3-D (rotate, zoom, pick
the three features), the graphs of the structures, a tree to switch structures, parents or leaves
on and off, and an outlier selector working on any of the per-point scores, by percentile or by
an explicit threshold, over all shown points or within each structure.
"""
from __future__ import annotations

import base64
import json
from pathlib import Path

import numpy as np

from ._util import as_points

__all__ = ["export_viewer", "viewer_payload"]

SCORES = [
    # key, label, low values are outliers?
    ("logdens", "log density (global model)", True),
    ("dist", "orthogonal distance from the structure", False),
    ("d2", "orthogonal Mahalanobis distance (squared)", False),
    ("pval", "model tail probability of the orthogonal offset", True),
    ("rank", "regional rank of the orthogonal distance", False),
    ("bg", "background probability", False),
]


def _b64(a, dtype) -> str:
    return base64.b64encode(np.ascontiguousarray(a, dtype=dtype).tobytes()).decode("ascii")


def _graph_layout(means: np.ndarray, edges: np.ndarray, dim: int):
    """Planar layout of an abstract graph and, for chains, the arc length of every node.

    Uses a force layout (networkx) when available, otherwise the two leading principal
    components of the node positions."""
    n = len(means)
    arc = np.zeros(0)
    if dim == 1 and n >= 2:
        try:
            from .crawling import ManifoldGraph
            from .viz import chain_order

            order, s = chain_order(ManifoldGraph(means, edges, 1))
            if len(order) == n:
                arc = np.zeros(n)
                arc[order] = s
        except Exception:  # noqa: BLE001 - layout is best effort
            arc = np.zeros(0)
    try:
        import networkx as nx

        G = nx.Graph()
        G.add_nodes_from(range(n))
        G.add_edges_from(np.asarray(edges).tolist())
        pos = nx.kamada_kawai_layout(G) if n <= 1200 else nx.spring_layout(G, seed=0, iterations=60)
        xy = np.array([pos[k] for k in range(n)], dtype=float)
    except Exception:  # noqa: BLE001
        c = means - means.mean(0)
        xy = c @ np.linalg.svd(c, full_matrices=False)[2][:2].T if n > 2 else np.c_[np.arange(n), np.zeros(n)]
    return xy, arc


def _version_payload(H, X, hops: int, min_region: int) -> dict:
    res = H.regional_deviation(X, hops=hops, min_region=min_region)
    leaves = H.leaves()
    structures, node_id = [], {}

    def walk(node, parent_id):
        for c in node.children:
            sid = len(structures)
            node_id[id(c)] = sid
            g = c.model.graph
            xy, arc = _graph_layout(c.model.means_, g.edges, int(g.dim))
            structures.append({"id": sid, "label": c.label, "parent": parent_id, "level": c.level, "dim": int(g.dim),
                               "kind": "group" if c.children else "leaf", "weight": float(c.weight),
                               "nodes": _b64(c.model.means_, np.float32), "n_nodes": int(g.n_nodes),
                               "edges": _b64(g.edges, np.int32), "graph_of": sid,
                               "layout": _b64(xy, np.float32), "arc": _b64(arc, np.float32) if len(arc) else ""})
            walk(c, sid)
            if c.children and c.pi_residual > 0:
                rid = len(structures)
                structures.append({"id": rid, "label": f"{c.label} (residual)", "parent": sid, "level": c.level + 1,
                                   "dim": int(g.dim), "kind": "residual", "weight": float(c.weight * c.pi_residual),
                                   "nodes": "", "n_nodes": 0, "edges": "", "graph_of": sid,
                                   "layout": "", "arc": ""})
                node_id[("res", id(c))] = rid

    walk(H.root_, -1)
    # map leaf index (into H.leaves()) to structure id; background leaves never win in regional_deviation
    lookup = np.full(len(leaves), -1)
    for k, (n, kind, _) in enumerate(leaves):
        if kind != "background":
            lookup[k] = node_id[("res", id(n))] if kind == "residual" else node_id[id(n)]
    sid = lookup[res["leaf"]]
    counts = np.bincount(sid[sid >= 0], minlength=len(structures))
    for s in structures:
        s["points"] = int(counts[s["id"]])
    return {
        "structures": structures,
        "sid": _b64(sid, np.int16),
        "node": _b64(res["node"], np.int32),
        "toward": _b64(res["toward"], np.int32),
        "frac": _b64(res["fraction"], np.float32),
        "scores": {"logdens": _b64(H.score_samples(X), np.float32), "dist": _b64(res["distance"], np.float32),
                   "d2": _b64(res["d2"], np.float32), "pval": _b64(res["p_value"], np.float32),
                   "rank": _b64(res["regional_rank"], np.float32), "bg": _b64(res["background_probability"], np.float32)},
        "summary": H.summary(),
    }


def _extra_payload(extra, n: int) -> list:
    out = []
    for name, values in (extra or {}).items():
        a = np.asarray(values)
        if a.shape != (n,):
            raise ValueError(f"extra column {name!r} must hold one value per point")
        if np.issubdtype(a.dtype, np.number) or a.dtype == bool:
            out.append({"name": str(name), "kind": "num", "values": _b64(a.astype(float), np.float32), "cats": []})
        else:
            cats, codes = np.unique(a.astype(str), return_inverse=True)
            out.append({"name": str(name), "kind": "cat", "values": _b64(codes, np.int32), "cats": cats.tolist()})
    return out


def viewer_payload(X, versions: dict, names=None, hops: int = 2, min_region: int = 30, extra=None) -> dict:
    """The data embedded in the viewer page (JSON-serialisable)."""
    X = as_points(X)
    n, d = X.shape
    names = list(names) if names is not None else [f"x{k + 1}" for k in range(d)]
    if len(names) != d:
        raise ValueError("names must have one entry per feature")
    if not versions:
        raise ValueError("at least one version is required")
    return {"n": n, "d": d, "names": names, "X": _b64(X, np.float32), "extra": _extra_payload(extra, n),
            "scores": [{"key": k, "label": l, "low": low} for k, l, low in SCORES],
            "versions": [{"name": str(name), **_version_payload(H, X, hops, min_region)} for name, H in versions.items()]}


def export_viewer(path, X, versions, names=None, title: str = "examm viewer", plotlyjs: str = "inline",
                  hops: int = 2, min_region: int = 30, extra=None) -> Path:
    """Write the standalone viewer page.

    Parameters
    ----------
    path : file to write (``.html``)
    X : (n, d) data
    versions : dict name -> fitted :class:`examm.HierarchicalAGTM` (or a single hierarchy)
    names : feature names
    plotlyjs : ``"inline"`` embeds plotly.js in the page (works offline, about 5 MB);
        ``"cdn"`` loads it from the plotly CDN.
    hops, min_region : passed to :meth:`HierarchicalAGTM.regional_deviation`
    extra : dict name -> (n,) array, optional
        Additional per-point labels (numeric, or text/categorical) offered for colouring the
        points and the graph nodes, e.g. physical parameters. More can be added later in the
        page itself by loading a CSV file.
    """
    if not isinstance(versions, dict):
        versions = {"model": versions}
    payload = viewer_payload(X, versions, names, hops, min_region, extra)
    if plotlyjs == "inline":
        from plotly.offline import get_plotlyjs

        lib = "<script>" + get_plotlyjs() + "</script>"
    elif plotlyjs == "cdn":
        lib = '<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>'
    else:
        raise ValueError("plotlyjs must be 'inline' or 'cdn'")
    html = (_TEMPLATE.replace("/*__DATA__*/", json.dumps(payload, separators=(",", ":")))
            .replace("<!--__PLOTLY__-->", lib).replace("__TITLE__", title))
    path = Path(path)
    path.write_text(html, encoding="utf-8")
    return path


_TEMPLATE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
:root{--ink:#0b0b0b;--muted:#52514e;--line:#dddcd6;--bg:#fcfcfb;--panel:#f4f3ef;--accent:#2a78d6;--out:#e34948}
*{box-sizing:border-box}html,body{height:100%;margin:0}
body{font:13px/1.4 system-ui,-apple-system,Segoe UI,Roboto,sans-serif;color:var(--ink);background:var(--bg);display:flex}
#side{width:350px;min-width:350px;height:100%;overflow:auto;background:var(--panel);border-right:1px solid var(--line);padding:12px 14px}
#main{flex:1;min-width:0;height:100%;display:flex;flex-direction:column;overflow:hidden}
#plot{flex:1 1 0;min-height:0;height:0;overflow:hidden;position:relative}
#status{flex:none;height:30px;line-height:18px;background:var(--bg);position:relative;z-index:2;padding:6px 12px;border-top:1px solid var(--line);color:var(--muted);font-size:12px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
h1{font-size:15px;margin:0 0 8px}h2{font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:var(--muted);margin:16px 0 6px}
label{display:block;margin:6px 0 2px;color:var(--muted);font-size:12px}
select,input[type=number]{width:100%;padding:4px 6px;border:1px solid var(--line);border-radius:4px;background:#fff;font:inherit}
.row{display:flex;gap:6px;align-items:center}.row>*{flex:1}.row.tight>*{flex:none}
button{font:inherit;padding:4px 8px;border:1px solid var(--line);border-radius:4px;background:#fff;cursor:pointer}
button:hover{border-color:var(--accent)}
#tabs{flex:none;display:flex;gap:4px;padding:8px 12px 0;border-bottom:1px solid var(--line);background:var(--panel)}
.tab{border:1px solid var(--line);border-bottom:none;border-radius:6px 6px 0 0;background:#e9e8e3;padding:6px 14px}
.tab.active{background:var(--bg);font-weight:600;position:relative;top:1px}
#dimtog{display:flex;flex-wrap:wrap;gap:10px;margin:0 0 6px}
#dimtog label{display:flex;align-items:center;gap:4px;margin:0;color:var(--ink)}
#catlegend{display:flex;flex-wrap:wrap;gap:4px 10px;margin-top:4px;font-size:11px}#catlegend span{display:flex;align-items:center;gap:4px}
.viewonly{display:none}body.v3d .only3d{display:block}body.vgraph .onlygraph{display:block}
#tree{border:1px solid var(--line);border-radius:4px;background:#fff;max-height:300px;overflow:auto;padding:4px 0}
.node{display:flex;align-items:center;gap:5px;padding:1px 6px;white-space:nowrap;cursor:default}
.node:hover{background:#f1f5fb}.node input{margin:0}
.sw{width:10px;height:10px;border-radius:2px;flex:none}.cnt{color:var(--muted);font-size:11px;margin-left:auto;padding-left:8px}
.kind{color:var(--muted);font-size:11px}
.check{display:flex;align-items:center;gap:6px;margin:4px 0;color:var(--ink)}.check input{margin:0}
#outinfo{margin-top:6px;padding:6px 8px;background:#fff;border:1px solid var(--line);border-radius:4px;font-size:12px}
#outinfo b{color:var(--out)}
details{margin-top:10px}summary{cursor:pointer;color:var(--muted);font-size:12px}
pre{font-size:10.5px;white-space:pre;overflow:auto;max-height:240px;background:#fff;border:1px solid var(--line);border-radius:4px;padding:6px}
.hint{color:var(--muted);font-size:11px;margin:4px 0 0}
</style><!--__PLOTLY__--></head>
<body>
<div id="side">
  <h1>__TITLE__</h1>
  <label for="version">Version</label><select id="version"></select>
  <h2>View</h2>
  <div class="row"><div><label for="ax0">x</label><select id="ax0"></select></div><div><label for="ax1">y</label><select id="ax1"></select></div><div><label for="ax2">z</label><select id="ax2"></select></div></div>
  <label for="colour">Colour points by</label><select id="colour"></select>
  <div class="row"><div><label for="psize">Point size</label><input id="psize" type="range" min="1" max="6" step="0.5" value="2"></div>
  <div><label for="popac">Opacity</label><input id="popac" type="range" min="0.1" max="1" step="0.05" value="0.7"></div></div>
  <label class="check"><input type="checkbox" id="showpts" checked> points</label>
  <label class="check"><input type="checkbox" id="showgraphs" checked> graphs of the selected structures</label>
  <div class="viewonly onlygraph">
  <h2>Graph view</h2>
  <label for="gstruct">Structure</label><select id="gstruct"></select>
  <label for="gcolour">Colour nodes by (surfaces)</label><select id="gcolour"></select>
  <label for="gy">Vertical axis (curves)</label><select id="gy"></select>
  <label class="check"><input type="checkbox" id="gpts" checked> points of the structure</label>
  <p class="hint">Surfaces and higher: the abstract graph laid out in the plane; each point sits at its position on the graph (with a little jitter). Curves: each point at its position along the curve (horizontal) against the outlier score (vertical); grey lines mark the nodes.</p>
  </div>
  <h2>Structures</h2>
  <div id="dimtog"></div>
  <div class="row tight" style="margin-bottom:6px"><button id="all">All</button><button id="none">None</button><button id="leavesonly">Leaves only</button><button id="top">Residuals only</button></div>
  <div id="tree"></div>
  <p class="hint">Ticking a parent ticks everything under it. "(residual)" holds the points of a parent that none of its children explains; its graph is the parent's.</p>
  <h2>Extra labels</h2>
  <input type="file" id="xfile" accept=".csv,.tsv,.txt">
  <p class="hint">Load a CSV with a header line to colour by your own per-point values (numbers or categories): either one line per point in the original order, or with a column named <b>row</b> holding the point index. The file stays in your browser.</p>
  <div id="xinfo" class="hint"></div><div id="catlegend"></div>
  <h2>Outliers</h2>
  <label for="score">Score</label><select id="score"></select>
  <label class="check"><input type="radio" name="mode" value="pct" checked> by percentile</label>
  <div class="row"><input id="pct" type="range" min="0.05" max="20" step="0.05" value="1"><input id="pctnum" type="number" min="0.01" max="50" step="0.05" value="1" style="max-width:90px"></div>
  <label class="check"><input type="radio" name="mode" value="thr"> by threshold</label>
  <div class="row"><input id="thr" type="range" step="any"><input id="thrnum" type="number" step="any" style="max-width:110px"></div>
  <label for="scope">Applied</label><select id="scope"><option value="shown">over all shown points</option><option value="each">within each shown structure</option></select>
  <label class="check"><input type="checkbox" id="hl" checked> highlight outliers</label>
  <label class="check"><input type="checkbox" id="only"> show outliers only</label>
  <div id="outinfo"></div>
  <div class="row tight" style="margin-top:6px"><button id="csv">Download outliers (CSV)</button><button id="csvall">Download shown points (CSV)</button></div>
  <details><summary>Model summary of this version</summary><pre id="summary"></pre></details>
</div>
<div id="main"><div id="tabs"><button class="tab active" data-view="3d">3-D view</button><button class="tab" data-view="graph">Graph view (one structure, unfolded)</button></div>
<div id="plot"></div><div id="status"></div></div>
<script>
const DATA = /*__DATA__*/;
const PAL = ["#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4","#008300","#6250d6","#e34948"];
const DIMCOL = ["#b9b8b2","#e34948","#2a78d6","#1baf7a","#eda100","#6250d6","#e87ba4","#008300","#52514e"];
function dec(b64, T){ if(!b64) return new T(0); const s=atob(b64), u=new Uint8Array(s.length); for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i); return new T(u.buffer); }
function colourOf(i){ if(i<PAL.length) return PAL[i]; const h=(i*137.508)%360; return `hsl(${h.toFixed(0)},62%,${i%2?44:56}%)`; }
const N=DATA.n, D=DATA.d, X=dec(DATA.X, Float32Array);
const V = DATA.versions.map(v => ({name:v.name, summary:v.summary, sid:dec(v.sid,Int16Array), node:dec(v.node,Int32Array),
  toward:dec(v.toward,Int32Array), frac:dec(v.frac,Float32Array),
  scores:Object.fromEntries(Object.entries(v.scores).map(([k,b])=>[k,dec(b,Float32Array)])),
  S:v.structures.map(s=>({...s, nodes:dec(s.nodes,Float32Array), edges:dec(s.edges,Int32Array), layout:dec(s.layout,Float32Array), arc:dec(s.arc,Float32Array), on:s.kind!=="group", colour:colourOf(s.id)}))}));
let view = "3d";
const EXTRA = (DATA.extra||[]).map(e=>({name:e.name, kind:e.kind, cats:e.cats||[], values: e.kind==="num" ? dec(e.values,Float32Array) : dec(e.values,Int32Array)}));
function extraOf(key){ return key && key.startsWith("x:") ? EXTRA[+key.slice(2)] : null; }
function refreshExtraOptions(){
  const keep=(id,items)=>{ const el=$(id), old=el.value; fill(el, items); if(items.some(([v])=>v===old)) el.value=old; };
  const xs=EXTRA.map((e,i)=>["x:"+i, e.name+(e.kind==="cat"?" (categories)":"")]), xn=xs.filter((_,i)=>EXTRA[i].kind==="num");
  keep("colour", [["structure","structure"],["dim","dimension of the structure"]].concat(DATA.scores.map(s=>[s.key,s.label]), xs));
  keep("gcolour", [["count","points per node"],["median","median of the outlier score"],["nout","outliers per node"]].concat(DATA.names.map((n,i)=>["pc:"+i,"mean "+n]), xn.map(([v,t])=>[v,"mean "+t])));
  keep("gy", [["score","the outlier score"]].concat(DATA.names.map((n,i)=>["pc:"+i,n]), xn));
}
function pointValue(key){ // per-point numeric accessor for "pc:k" or a numeric extra
  if(key.startsWith("pc:")){ const k=+key.slice(3); return i=>X[i*D+k]; }
  const e=extraOf(key); return e && e.kind==="num" ? (i=>e.values[i]) : null; }
function parseLine(line, delim){ const out=[]; let cur="", q=false;
  for(let i=0;i<line.length;i++){ const c=line[i]; if(q){ if(c==='"'){ if(line[i+1]==='"'){cur+='"';i++;} else q=false; } else cur+=c; } else if(c==='"') q=true; else if(c===delim){ out.push(cur); cur=""; } else cur+=c; }
  out.push(cur); return out.map(t=>t.trim()); }
function loadExtraText(text, fname){
  const lines=text.split(/\r?\n/).filter(l=>l.trim().length); if(lines.length<2) throw new Error("the file has no data lines");
  const delim = lines[0].includes("\t") ? "\t" : lines[0].includes(";") && !lines[0].includes(",") ? ";" : ",";
  const head=parseLine(lines[0],delim), rows=lines.slice(1).map(l=>parseLine(l,delim));
  const ri=head.findIndex(h=>h.toLowerCase()==="row"); let index;
  if(ri>=0){ index=rows.map(r=>parseInt(r[ri],10)); if(index.some(k=>!(k>=0&&k<N))) throw new Error("the row column holds values outside 0.."+(N-1)); }
  else { if(rows.length!==N) throw new Error(`the file has ${rows.length} data lines but there are ${N} points; add a "row" column with the point index`); index=rows.map((_,k)=>k); }
  const added=[], skipped=[];
  head.forEach((name,c)=>{ if(c===ri||!name) return;
    const raw=new Array(N).fill(""); rows.forEach((r,k)=>{ raw[index[k]] = r[c]===undefined ? "" : r[c]; });
    const filled=raw.filter(t=>t!==""&&t.toLowerCase()!=="nan"); if(!filled.length){ skipped.push(name+" (empty)"); return; }
    let e;
    if(filled.every(t=>Number.isFinite(Number(t)))){ const vals=new Float32Array(N); for(let i=0;i<N;i++){ const t=raw[i]; vals[i] = (t===""||t.toLowerCase()==="nan") ? NaN : Number(t); } e={name,kind:"num",cats:[],values:vals}; }
    else { const cats=[...new Set(filled)].sort(); if(cats.length>60){ skipped.push(name+" (more than 60 categories)"); return; }
      const code=new Map(cats.map((t,k)=>[t,k])), vals=new Int32Array(N); for(let i=0;i<N;i++) vals[i]=code.has(raw[i])?code.get(raw[i]):-1; e={name,kind:"cat",cats,values:vals}; }
    const old=EXTRA.findIndex(x=>x.name===name); if(old>=0) EXTRA[old]=e; else EXTRA.push(e); added.push(name); });
  refreshExtraOptions();
  $("xinfo").textContent = `${fname}: ${rows.length.toLocaleString()} lines, added ${added.length? added.join(", "):"nothing"}`+(skipped.length? "; skipped "+skipped.join(", "):"");
  if(added.length){ const last=EXTRA.findIndex(x=>x.name===added[0]); $("colour").value="x:"+last; }
  draw();
}
const $ = id => document.getElementById(id);
let cur = 0, lastOut = null;
function fill(sel, items, val){ sel.innerHTML=""; items.forEach(([v,t])=>{const o=document.createElement("option");o.value=v;o.textContent=t;sel.appendChild(o);}); if(val!==undefined) sel.value=val; }
fill($("version"), V.map((v,i)=>[i,v.name]));
for(let a=0;a<3;a++) fill($("ax"+a), DATA.names.map((n,i)=>[i,n]), Math.min(a,D-1));

fill($("score"), DATA.scores.map(s=>[s.key,s.label]));
function children(v, id){ return v.S.filter(s=>s.parent===id); }
function setSubtree(v, s, on){ s.on = on; children(v,s.id).forEach(c=>setSubtree(v,c,on)); }
function buildTree(){
  const v=V[cur], box=$("tree"); box.innerHTML=""; syncParents(v);
  const add=(s,depth)=>{
    const row=document.createElement("div"); row.className="node"; row.style.paddingLeft=(6+depth*16)+"px";
    const cb=document.createElement("input"); cb.type="checkbox"; cb.checked=s.on; cb.dataset.id=s.id;
    if(s.kind==="group"){ const st=leafStates(v,s); cb.indeterminate = st.some(x=>x) && !st.every(x=>x); }
    cb.addEventListener("change",()=>{ setSubtree(v,s,cb.checked); syncParents(v); buildTree(); draw(); });
    const sw=document.createElement("span"); sw.className="sw"; sw.style.background = s.kind==="group" ? "transparent" : s.colour; if(s.kind==="group") sw.style.border="1px solid var(--muted)";
    const lab=document.createElement("span"); lab.textContent=s.label;
    const kd=document.createElement("span"); kd.className="kind"; kd.textContent = s.dim+"-D"+(s.kind==="group"?" · parent":"");
    const cnt=document.createElement("span"); cnt.className="cnt"; cnt.textContent = s.kind==="group" ? subtreePoints(v,s).toLocaleString() : s.points.toLocaleString();
    row.append(cb,sw,lab,kd,cnt); box.appendChild(row);
    children(v,s.id).forEach(c=>add(c,depth+1));
  };
  children(v,-1).forEach(s=>add(s,0));
  $("summary").textContent=v.summary;
}
function leafStates(v,s){ return s.kind==="group" ? children(v,s.id).flatMap(c=>leafStates(v,c)) : [s.on]; }
function subtreePoints(v,s){ return s.kind==="group" ? children(v,s.id).reduce((a,c)=>a+subtreePoints(v,c),0) : s.points; }
function syncParents(v){ const rec=s=>{ if(s.kind!=="group") return s.on; const k=children(v,s.id).map(rec); s.on=k.some(x=>x); return s.on; }; children(v,-1).forEach(rec); }
function quantile(arr, q){ if(!arr.length) return NaN; const a=Float64Array.from(arr).sort(); const p=(a.length-1)*q, lo=Math.floor(p), hi=Math.ceil(p); return a[lo]+(a[hi]-a[lo])*(p-lo); }
function scoreMeta(){ return DATA.scores.find(s=>s.key===$("score").value); }
function setThrRange(){
  const v=V[cur], s=v.scores[$("score").value]; let lo=Infinity, hi=-Infinity;
  for(let i=0;i<N;i++){ const x=s[i]; if(Number.isFinite(x)){ if(x<lo)lo=x; if(x>hi)hi=x; } }
  if(!Number.isFinite(lo)){lo=0;hi=1;}
  const vals=[]; for(let i=0;i<N;i+=7){ if(Number.isFinite(s[i])) vals.push(s[i]); }
  const start = scoreMeta().low ? quantile(vals,0.01) : quantile(vals,0.99);
  $("thr").min=lo; $("thr").max=hi; $("thr").value=start; $("thrnum").value=(+start).toPrecision(4);
}
function compute(only){
  // only: if given, the selection is restricted to that one structure (graph view)
  const v=V[cur], on=v.S.map(s=> only===undefined ? (s.on && s.kind!=="group") : s.id===only), sid=v.sid, sc=v.scores[$("score").value], low=scoreMeta().low;
  const mode=document.querySelector("input[name=mode]:checked").value, scope=$("scope").value;
  const shown=new Uint8Array(N); let nshown=0;
  for(let i=0;i<N;i++){ const s=sid[i]; if(s>=0 && on[s]){shown[i]=1;nshown++;} }
  const out=new Uint8Array(N); let nout=0, info=""; const cuts=new Map();
  if(mode==="pct"){
    const q=Math.min(Math.max(+$("pctnum").value,0),100)/100;
    const groups = new Map();
    for(let i=0;i<N;i++){ if(!shown[i]||!Number.isFinite(sc[i])) continue; const key = scope==="each" ? sid[i] : -1; if(!groups.has(key)) groups.set(key,[]); groups.get(key).push(i); }
    const thr=[];
    groups.forEach((idx,key)=>{ const t=quantile(idx.map(i=>sc[i]), low?q:1-q); thr.push(t); cuts.set(key,t); idx.forEach(i=>{ if(low? sc[i]<t : sc[i]>t){out[i]=1;nout++;} }); });
    info = scope==="each" ? `${groups.size} structures, each with its own cut` : `cut at ${thr.length?thr[0].toPrecision(4):"–"}`;
  } else {
    const t=+$("thrnum").value; cuts.set(-1,t);
    for(let i=0;i<N;i++){ if(shown[i] && Number.isFinite(sc[i]) && (low? sc[i]<t : sc[i]>t)){out[i]=1;nout++;} }
    info=`${low?"score below":"score above"} ${t.toPrecision(4)}`;
  }
  return {shown,out,nshown,nout,info,low,cuts};
}
function seqColour(vals){ // single-hue ramp, robust range
  const f=vals.filter(Number.isFinite); const lo=quantile(f,0.02), hi=quantile(f,0.98); return {cmin:lo,cmax:hi}; }
function draw(){ buildDimToggles(); catLegend(); if(view==="graph") drawGraph(); else draw3d(); }
function catLegend(){ const e=extraOf($("colour").value), box=$("catlegend"); box.innerHTML="";
  if(!e||e.kind!=="cat"||view==="graph") return;
  e.cats.forEach((t,k)=>{ const sp=document.createElement("span"), sw=document.createElement("span"); sw.className="sw"; sw.style.background=colourOf(k); sp.append(sw, document.createTextNode(t)); box.appendChild(sp); }); }
function buildDimToggles(){
  const v=V[cur], box=$("dimtog"); box.innerHTML="";
  const dims=[...new Set(v.S.filter(s=>s.kind!=="group").map(s=>s.dim))].sort((a,b)=>a-b);
  dims.forEach(dm=>{
    const mine=v.S.filter(s=>s.kind!=="group" && s.dim===dm), st=mine.map(s=>s.on);
    const lab=document.createElement("label"), cb=document.createElement("input"); cb.type="checkbox";
    cb.checked=st.every(x=>x); cb.indeterminate=st.some(x=>x)&&!st.every(x=>x);
    cb.addEventListener("change",()=>{ mine.forEach(s=>{s.on=cb.checked;}); syncParents(v); buildTree(); draw(); });
    const n=mine.reduce((a,s)=>a+s.points,0);
    lab.append(cb, document.createTextNode(`${dm}-D`), Object.assign(document.createElement("span"),{className:"kind",textContent:`(${mine.length}, ${n.toLocaleString()} pts)`}));
    box.appendChild(lab); });
}
function fillGraphSelect(){
  const v=V[cur], opts=v.S.filter(s=>s.kind!=="group" && s.points>0).sort((a,b)=>b.points-a.points).map(s=>[s.id, `${s.label} · ${s.dim}-D · ${s.points.toLocaleString()} points`]);
  fill($("gstruct"), opts);
}
function hash01(i,k){ let x=((i+1)*2654435761 + k*40503)>>>0; x^=x>>>15; x=Math.imul(x,2246822519)>>>0; x^=x>>>13; return x/4294967296; }
function median(a){ if(!a.length) return NaN; const b=Float64Array.from(a).sort(); const m=b.length>>1; return b.length%2? b[m] : 0.5*(b[m-1]+b[m]); }
function drawGraph(){
  const v=V[cur], sidSel=+$("gstruct").value, s=v.S[sidSel];
  if(!s){ Plotly.react("plot",[],{}); return; }
  const g=v.S[s.graph_of], K=g.n_nodes, L=g.layout, E=g.edges, P=g.nodes;
  const r=compute(sidSel); lastOut=r;
  const key=$("score").value, sc=v.scores[key], meta=scoreMeta(), hl=$("hl").checked, only=$("only").checked;
  const idx=[]; for(let i=0;i<N;i++) if(v.sid[i]===sidSel) idx.push(i);
  const traces=[], ht="row %{customdata[0]}<br>score %{customdata[1]:.4g}<extra></extra>";
  const is1d = g.dim===1 && g.arc.length===K;
  const cutKey = (document.querySelector("input[name=mode]:checked").value==="pct" && $("scope").value==="each") ? sidSel : -1, cut = r.cuts.has(cutKey) ? r.cuts.get(cutKey) : NaN;
  let layout;
  if(is1d){
    const elen=(a,b)=>{ let q=0; for(let k=0;k<D;k++){ const t=P[a*D+k]-P[b*D+k]; q+=t*t; } return Math.sqrt(q); };
    const px=[],py=[],cd=[],ox=[],oy=[],ocd=[];
    const yk=$("gy").value, yf = yk==="score" ? (i=>sc[i]) : pointValue(yk), ylab = yk==="score" ? meta.label : (yk.startsWith("pc:") ? DATA.names[+yk.slice(3)] : extraOf(yk).name);
    idx.forEach(i=>{ const a=v.node[i], b=v.toward[i]; const pos=g.arc[a]+Math.sign(g.arc[b]-g.arc[a])*v.frac[i]*elen(a,b), y=yf(i);
      if(!Number.isFinite(y)) return;
      if(r.out[i]){ ox.push(pos);oy.push(y);ocd.push([i,sc[i]]); if(hl||only) return; } if(only) return; px.push(pos);py.push(y);cd.push([i,sc[i]]); });
    if($("gpts").checked) traces.push({type:"scattergl",mode:"markers",x:px,y:py,customdata:cd,hovertemplate:ht,marker:{size:5,color:s.colour,opacity:0.7},name:"points",showlegend:false});
    if(ox.length) traces.push({type:"scattergl",mode:"markers",x:ox,y:oy,customdata:ocd,hovertemplate:ht,marker:{size:9,color:"#e34948",symbol:"diamond",line:{color:"#fff",width:0.5}},name:"outliers",showlegend:false});
    const shapes=[]; for(let k=0;k<K;k++) shapes.push({type:"line",xref:"x",yref:"paper",x0:g.arc[k],x1:g.arc[k],y0:0,y1:1,line:{color:"#dddcd6",width:1}});
    if(Number.isFinite(cut) && yk==="score") shapes.push({type:"line",xref:"paper",yref:"y",x0:0,x1:1,y0:cut,y1:cut,line:{color:"#e34948",width:1,dash:"dash"}});
    layout={shapes, xaxis:{title:{text:"position along the curve (arc length)"},gridcolor:"#eeede8",zeroline:false}, yaxis:{title:{text:ylab},gridcolor:"#eeede8",zeroline:false}};
  } else {
    const n=E.length/2, ex=new Array(3*n), ey=new Array(3*n); let med=[];
    for(let e=0;e<n;e++){ const a=E[2*e], b=E[2*e+1]; ex[3*e]=L[2*a];ex[3*e+1]=L[2*b];ex[3*e+2]=null; ey[3*e]=L[2*a+1];ey[3*e+1]=L[2*b+1];ey[3*e+2]=null; med.push(Math.hypot(L[2*a]-L[2*b],L[2*a+1]-L[2*b+1])); }
    const jit=0.3*(median(med)||0.01);
    traces.push({type:"scattergl",mode:"lines",x:ex,y:ey,line:{color:"#c9c8c2",width:1},hoverinfo:"skip",showlegend:false});
    const cnt=new Float64Array(K), nout=new Float64Array(K), lists=Array.from({length:K},()=>[]), cm=$("gcolour").value;
    const pcSum=new Float64Array(K), pcN=new Float64Array(K), pf=pointValue(cm);
    const px=[],py=[],cd=[],ox=[],oy=[],ocd=[];
    idx.forEach(i=>{ const a=v.node[i], b=v.toward[i], f=v.frac[i]; cnt[a]++; if(r.out[i]) nout[a]++; if(Number.isFinite(sc[i])) lists[a].push(sc[i]); if(pf){ const t=pf(i); if(Number.isFinite(t)){ pcSum[a]+=t; pcN[a]++; } }
      const ang=2*Math.PI*hash01(i,1), rad=jit*Math.sqrt(hash01(i,2));
      const x=L[2*a]+f*(L[2*b]-L[2*a])+rad*Math.cos(ang), y=L[2*a+1]+f*(L[2*b+1]-L[2*a+1])+rad*Math.sin(ang);
      if(r.out[i]){ ox.push(x);oy.push(y);ocd.push([i,sc[i]]); if(hl||only) return; } if(only) return; px.push(x);py.push(y);cd.push([i,sc[i]]); });
    if($("gpts").checked && px.length) traces.push({type:"scattergl",mode:"markers",x:px,y:py,customdata:cd,hovertemplate:ht,marker:{size:3,color:"#b9b8b2",opacity:0.6},name:"points",showlegend:false});
    const val=new Array(K);
    for(let k=0;k<K;k++){ val[k] = cm==="count" ? cnt[k] : cm==="nout" ? nout[k] : cm==="median" ? median(lists[k]) : (pcN[k]? pcSum[k]/pcN[k] : NaN); }
    const title = cm==="count" ? "points per node" : cm==="nout" ? "outliers per node" : cm==="median" ? "median "+meta.label : "mean "+(cm.startsWith("pc:") ? DATA.names[+cm.slice(3)] : extraOf(cm).name);
    const nx=[],ny=[],ncd=[]; for(let k=0;k<K;k++){ nx.push(L[2*k]); ny.push(L[2*k+1]); ncd.push([k,cnt[k],nout[k],val[k]]); }
    traces.push({type:"scattergl",mode:"markers",x:nx,y:ny,customdata:ncd,hovertemplate:"node %{customdata[0]}<br>%{customdata[1]} points, %{customdata[2]} outliers<br>value %{customdata[3]:.4g}<extra></extra>",
      marker:{size:8,color:val,colorscale:"Blues",reversescale:true,line:{color:"#fff",width:0.5},colorbar:{title:{text:title,side:"right"},thickness:12,len:0.6}},name:"nodes",showlegend:false});
    if(ox.length) traces.push({type:"scattergl",mode:"markers",x:ox,y:oy,customdata:ocd,hovertemplate:ht,marker:{size:9,color:"#e34948",symbol:"diamond",line:{color:"#fff",width:0.5}},name:"outliers",showlegend:false});
    layout={xaxis:{visible:false}, yaxis:{visible:false,scaleanchor:"x"}};
  }
  Plotly.react("plot", traces, Object.assign({margin:{l:is1d?60:10,r:10,t:10,b:is1d?50:10},paper_bgcolor:"#fcfcfb",plot_bgcolor:"#ffffff",uirevision:"g"+cur+"-"+sidSel,dragmode:"pan"},layout), {displaylogo:false,responsive:true,scrollZoom:true});
  $("outinfo").innerHTML = `<b>${r.nout.toLocaleString()}</b> outliers among ${r.nshown.toLocaleString()} points of this structure (${r.nshown? (100*r.nout/r.nshown).toFixed(2):"0"} %)<br><span class="kind">${meta.low?"low":"high"} values are outliers · ${Number.isFinite(cut)? "cut at "+cut.toPrecision(4) : r.info}</span>`;
  $("status").textContent = `${v.name} · structure ${s.label}: ${s.dim}-D, ${K} nodes, ${idx.length.toLocaleString()} points · drag to pan, scroll to zoom`;
}
function setView(w){ view=w; document.body.className = w==="graph" ? "vgraph" : "v3d"; document.querySelectorAll(".tab").forEach(t=>t.classList.toggle("active", t.dataset.view===w)); Plotly.purge("plot"); draw(); }
function draw3d(){
  const v=V[cur], ax=[+$("ax0").value,+$("ax1").value,+$("ax2").value], cm=$("colour").value;

  const r=compute(); lastOut=r;
  const only=$("only").checked, hl=$("hl").checked, sc=v.scores[$("score").value], xe=extraOf(cm);
  const px=[],py=[],pz=[],pc=[],cd=[], ox=[],oy=[],oz=[],ocd=[];
  for(let i=0;i<N;i++){
    if(!r.shown[i]) continue; const b=i*D, isOut=r.out[i]===1;
    if(isOut && (hl||only)){ ox.push(X[b+ax[0]]);oy.push(X[b+ax[1]]);oz.push(X[b+ax[2]]); ocd.push([i,v.S[v.sid[i]].label,sc[i]]); if(hl) continue; }
    if(only) continue;
    px.push(X[b+ax[0]]);py.push(X[b+ax[1]]);pz.push(X[b+ax[2]]); cd.push([i,v.S[v.sid[i]].label,sc[i]]);
    if(cm==="structure") pc.push(v.S[v.sid[i]].colour); else if(cm==="dim") pc.push(DIMCOL[Math.min(v.S[v.sid[i]].dim,8)]);
    else if(xe) pc.push(xe.kind==="num" ? xe.values[i] : (xe.values[i]<0 ? "#b9b8b2" : colourOf(xe.values[i]))); else pc.push(v.scores[cm][i]);
  }
  const traces=[];
  const ht="row %{customdata[0]} · %{customdata[1]}<br>score %{customdata[2]:.4g}<extra></extra>";
  if($("showpts").checked && px.length){
    const marker={size:+$("psize").value, opacity:+$("popac").value, color:pc, line:{width:0}};
    if(cm!=="structure" && cm!=="dim" && !(xe && xe.kind==="cat")){ const rg=seqColour(pc); Object.assign(marker,{colorscale:"Blues",reversescale:true,cmin:rg.cmin,cmax:rg.cmax,
      colorbar:{title:{text: xe ? xe.name : DATA.scores.find(s=>s.key===cm).label,side:"right"},thickness:12,len:0.6}}); }
    traces.push({type:"scatter3d",mode:"markers",x:px,y:py,z:pz,customdata:cd,hovertemplate:ht,marker,name:"points",showlegend:false});
  }
  if($("showgraphs").checked){
    const done=new Set();
    v.S.forEach(s=>{ if(!s.on||s.kind==="group") return; const g=v.S[s.graph_of]; if(done.has(g.id)||!g.n_nodes) return; done.add(g.id);
      const E=g.edges, P=g.nodes, n=E.length/2, gx=new Array(3*n), gy=new Array(3*n), gz=new Array(3*n);
      for(let e=0;e<n;e++){ const a=E[2*e]*D, b=E[2*e+1]*D; gx[3*e]=P[a+ax[0]];gx[3*e+1]=P[b+ax[0]];gx[3*e+2]=null; gy[3*e]=P[a+ax[1]];gy[3*e+1]=P[b+ax[1]];gy[3*e+2]=null; gz[3*e]=P[a+ax[2]];gz[3*e+1]=P[b+ax[2]];gz[3*e+2]=null; }
      traces.push({type:"scatter3d",mode:"lines",x:gx,y:gy,z:gz,line:{color: cm==="structure"? (s.kind==="residual"? s.colour : g.colour) : "#0b0b0b", width:3},hoverinfo:"skip",name:"graph "+g.label,showlegend:false}); });
  }
  if(ox.length) traces.push({type:"scatter3d",mode:"markers",x:ox,y:oy,z:oz,customdata:ocd,hovertemplate:ht,marker:{size:+$("psize").value+2.5,color:"#e34948",symbol:"diamond",line:{color:"#ffffff",width:0.5}},name:"outliers",showlegend:false});
  const axis=t=>({title:{text:t},backgroundcolor:"#ffffff",gridcolor:"#dddcd6",zerolinecolor:"#dddcd6",color:"#52514e"});
  Plotly.react("plot", traces, {margin:{l:0,r:0,t:0,b:0},paper_bgcolor:"#fcfcfb",uirevision:"keep",
    scene:{aspectmode:"data",xaxis:axis(DATA.names[ax[0]]),yaxis:axis(DATA.names[ax[1]]),zaxis:axis(DATA.names[ax[2]]),uirevision:"keep"}},
    {displaylogo:false,responsive:true});
  $("outinfo").innerHTML = `<b>${r.nout.toLocaleString()}</b> outliers among ${r.nshown.toLocaleString()} shown points (${r.nshown? (100*r.nout/r.nshown).toFixed(2):"0"} %)<br><span class="kind">${scoreMeta().low?"low":"high"} values are outliers · ${r.info}</span>`;
  $("status").textContent = `${v.name} · ${r.nshown.toLocaleString()} of ${N.toLocaleString()} points shown · drag to rotate, scroll to zoom, right-drag to pan`;
}
function download(onlyOut){
  const v=V[cur], r=lastOut||compute(); const keys=DATA.scores.map(s=>s.key);
  const rows=["row,structure,dim,node,outlier,"+keys.join(",")+","+DATA.names.join(",")];
  for(let i=0;i<N;i++){ if(!r.shown[i]||(onlyOut&&!r.out[i])) continue; const s=v.S[v.sid[i]];
    const xs=[]; for(let k=0;k<D;k++) xs.push(X[i*D+k].toPrecision(6));
    rows.push([i,'"'+s.label+'"',s.dim,v.node[i],r.out[i],...keys.map(k=>v.scores[k][i]),...xs].join(",")); }
  const a=document.createElement("a"); a.href=URL.createObjectURL(new Blob([rows.join("\n")],{type:"text/csv"})); a.download=(onlyOut?"outliers":"points")+".csv"; a.click(); URL.revokeObjectURL(a.href);
}
function setAll(fn){ const v=V[cur]; v.S.forEach(s=>{ s.on=fn(s,v); }); syncParents(v); buildTree(); draw(); }
$("all").onclick=()=>setAll(()=>true); $("none").onclick=()=>setAll(()=>false);
$("leavesonly").onclick=()=>setAll(s=>s.kind==="leaf");
$("top").onclick=()=>setAll(s=>s.kind==="residual");
$("version").onchange=()=>{ cur=+$("version").value; buildTree(); fillGraphSelect(); setThrRange(); draw(); };
document.querySelectorAll(".tab").forEach(t=>t.addEventListener("click",()=>setView(t.dataset.view)));
refreshExtraOptions();
["gstruct","gcolour","gpts","gy"].forEach(id=>$(id).addEventListener("input",draw));
$("xfile").addEventListener("change",ev=>{ const f=ev.target.files[0]; if(!f) return; const rd=new FileReader();
  rd.onload=()=>{ try{ loadExtraText(String(rd.result), f.name); }catch(err){ $("xinfo").textContent="Could not use "+f.name+": "+err.message; } };
  rd.readAsText(f); });
["ax0","ax1","ax2","colour","psize","popac","showpts","showgraphs","scope","hl","only"].forEach(id=>$(id).addEventListener("input",draw));
$("score").onchange=()=>{ setThrRange(); draw(); };
$("pct").oninput=()=>{ $("pctnum").value=$("pct").value; document.querySelector("input[name=mode][value=pct]").checked=true; draw(); };
$("pctnum").oninput=()=>{ $("pct").value=$("pctnum").value; document.querySelector("input[name=mode][value=pct]").checked=true; draw(); };
$("thr").oninput=()=>{ $("thrnum").value=(+$("thr").value).toPrecision(4); document.querySelector("input[name=mode][value=thr]").checked=true; draw(); };
$("thrnum").oninput=()=>{ $("thr").value=$("thrnum").value; document.querySelector("input[name=mode][value=thr]").checked=true; draw(); };
document.querySelectorAll("input[name=mode]").forEach(e=>e.addEventListener("change",draw));
$("csv").onclick=()=>download(true); $("csvall").onclick=()=>download(false);
document.body.className="v3d"; buildTree(); fillGraphSelect(); setThrRange(); draw();
</script></body></html>
"""
