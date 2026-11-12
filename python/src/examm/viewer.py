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


def _version_payload(H, X, hops: int, min_region: int) -> dict:
    res = H.regional_deviation(X, hops=hops, min_region=min_region)
    leaves = H.leaves()
    structures, node_id = [], {}

    def walk(node, parent_id):
        for c in node.children:
            sid = len(structures)
            node_id[id(c)] = sid
            g = c.model.graph
            structures.append({"id": sid, "label": c.label, "parent": parent_id, "level": c.level, "dim": int(g.dim),
                               "kind": "group" if c.children else "leaf", "weight": float(c.weight),
                               "nodes": _b64(c.model.means_, np.float32), "n_nodes": int(g.n_nodes),
                               "edges": _b64(g.edges, np.int32), "graph_of": sid})
            walk(c, sid)
            if c.children and c.pi_residual > 0:
                rid = len(structures)
                structures.append({"id": rid, "label": f"{c.label} (residual)", "parent": sid, "level": c.level + 1,
                                   "dim": int(g.dim), "kind": "residual", "weight": float(c.weight * c.pi_residual),
                                   "nodes": "", "n_nodes": 0, "edges": "", "graph_of": sid})
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
        "scores": {"logdens": _b64(H.score_samples(X), np.float32), "dist": _b64(res["distance"], np.float32),
                   "d2": _b64(res["d2"], np.float32), "pval": _b64(res["p_value"], np.float32),
                   "rank": _b64(res["regional_rank"], np.float32), "bg": _b64(res["background_probability"], np.float32)},
        "summary": H.summary(),
    }


def viewer_payload(X, versions: dict, names=None, hops: int = 2, min_region: int = 30) -> dict:
    """The data embedded in the viewer page (JSON-serialisable)."""
    X = as_points(X)
    n, d = X.shape
    names = list(names) if names is not None else [f"x{k + 1}" for k in range(d)]
    if len(names) != d:
        raise ValueError("names must have one entry per feature")
    if not versions:
        raise ValueError("at least one version is required")
    return {"n": n, "d": d, "names": names, "X": _b64(X, np.float32),
            "scores": [{"key": k, "label": l, "low": low} for k, l, low in SCORES],
            "versions": [{"name": str(name), **_version_payload(H, X, hops, min_region)} for name, H in versions.items()]}


def export_viewer(path, X, versions, names=None, title: str = "examm viewer", plotlyjs: str = "inline",
                  hops: int = 2, min_region: int = 30) -> Path:
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
    """
    if not isinstance(versions, dict):
        versions = {"model": versions}
    payload = viewer_payload(X, versions, names, hops, min_region)
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
#plot{flex:1 1 0;min-height:0;height:0}
#status{flex:none;padding:6px 12px;border-top:1px solid var(--line);color:var(--muted);font-size:12px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
h1{font-size:15px;margin:0 0 8px}h2{font-size:11px;text-transform:uppercase;letter-spacing:.06em;color:var(--muted);margin:16px 0 6px}
label{display:block;margin:6px 0 2px;color:var(--muted);font-size:12px}
select,input[type=number]{width:100%;padding:4px 6px;border:1px solid var(--line);border-radius:4px;background:#fff;font:inherit}
.row{display:flex;gap:6px;align-items:center}.row>*{flex:1}.row.tight>*{flex:none}
button{font:inherit;padding:4px 8px;border:1px solid var(--line);border-radius:4px;background:#fff;cursor:pointer}
button:hover{border-color:var(--accent)}
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
  <h2>Structures</h2>
  <div class="row tight" style="margin-bottom:6px"><button id="all">All</button><button id="none">None</button><button id="leavesonly">Leaves only</button><button id="top">Residuals only</button></div>
  <div id="tree"></div>
  <p class="hint">Ticking a parent ticks everything under it. "(residual)" holds the points of a parent that none of its children explains; its graph is the parent's.</p>
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
<div id="main"><div id="plot"></div><div id="status"></div></div>
<script>
const DATA = /*__DATA__*/;
const PAL = ["#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4","#008300","#6250d6","#e34948"];
const DIMCOL = ["#b9b8b2","#e34948","#2a78d6","#1baf7a","#eda100","#6250d6","#e87ba4","#008300","#52514e"];
function dec(b64, T){ if(!b64) return new T(0); const s=atob(b64), u=new Uint8Array(s.length); for(let i=0;i<s.length;i++)u[i]=s.charCodeAt(i); return new T(u.buffer); }
function colourOf(i){ if(i<PAL.length) return PAL[i]; const h=(i*137.508)%360; return `hsl(${h.toFixed(0)},62%,${i%2?44:56}%)`; }
const N=DATA.n, D=DATA.d, X=dec(DATA.X, Float32Array);
const V = DATA.versions.map(v => ({name:v.name, summary:v.summary, sid:dec(v.sid,Int16Array), node:dec(v.node,Int32Array),
  scores:Object.fromEntries(Object.entries(v.scores).map(([k,b])=>[k,dec(b,Float32Array)])),
  S:v.structures.map(s=>({...s, nodes:dec(s.nodes,Float32Array), edges:dec(s.edges,Int32Array), on:s.kind!=="group", colour:colourOf(s.id)}))}));
const $ = id => document.getElementById(id);
let cur = 0, lastOut = null;
function fill(sel, items, val){ sel.innerHTML=""; items.forEach(([v,t])=>{const o=document.createElement("option");o.value=v;o.textContent=t;sel.appendChild(o);}); if(val!==undefined) sel.value=val; }
fill($("version"), V.map((v,i)=>[i,v.name]));
for(let a=0;a<3;a++) fill($("ax"+a), DATA.names.map((n,i)=>[i,n]), Math.min(a,D-1));
fill($("colour"), [["structure","structure"],["dim","dimension of the structure"]].concat(DATA.scores.map(s=>[s.key,s.label])));
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
function compute(){
  const v=V[cur], on=v.S.map(s=>s.on && s.kind!=="group"), sid=v.sid, sc=v.scores[$("score").value], low=scoreMeta().low;
  const mode=document.querySelector("input[name=mode]:checked").value, scope=$("scope").value;
  const shown=new Uint8Array(N); let nshown=0;
  for(let i=0;i<N;i++){ const s=sid[i]; if(s>=0 && on[s]){shown[i]=1;nshown++;} }
  const out=new Uint8Array(N); let nout=0, info="";
  if(mode==="pct"){
    const q=Math.min(Math.max(+$("pctnum").value,0),100)/100;
    const groups = new Map();
    for(let i=0;i<N;i++){ if(!shown[i]||!Number.isFinite(sc[i])) continue; const key = scope==="each" ? sid[i] : -1; if(!groups.has(key)) groups.set(key,[]); groups.get(key).push(i); }
    const thr=[];
    groups.forEach((idx,key)=>{ const t=quantile(idx.map(i=>sc[i]), low?q:1-q); thr.push(t); idx.forEach(i=>{ if(low? sc[i]<t : sc[i]>t){out[i]=1;nout++;} }); });
    info = scope==="each" ? `${groups.size} structures, each with its own cut` : `cut at ${thr.length?thr[0].toPrecision(4):"–"}`;
  } else {
    const t=+$("thrnum").value;
    for(let i=0;i<N;i++){ if(shown[i] && Number.isFinite(sc[i]) && (low? sc[i]<t : sc[i]>t)){out[i]=1;nout++;} }
    info=`${low?"score below":"score above"} ${t.toPrecision(4)}`;
  }
  return {shown,out,nshown,nout,info,low};
}
function seqColour(vals){ // single-hue ramp, robust range
  const f=vals.filter(Number.isFinite); const lo=quantile(f,0.02), hi=quantile(f,0.98); return {cmin:lo,cmax:hi}; }
function draw(){
  const v=V[cur], ax=[+$("ax0").value,+$("ax1").value,+$("ax2").value], cm=$("colour").value;
  const r=compute(); lastOut=r;
  const only=$("only").checked, hl=$("hl").checked, sc=v.scores[$("score").value];
  const px=[],py=[],pz=[],pc=[],cd=[], ox=[],oy=[],oz=[],ocd=[];
  for(let i=0;i<N;i++){
    if(!r.shown[i]) continue; const b=i*D, isOut=r.out[i]===1;
    if(isOut && (hl||only)){ ox.push(X[b+ax[0]]);oy.push(X[b+ax[1]]);oz.push(X[b+ax[2]]); ocd.push([i,v.S[v.sid[i]].label,sc[i]]); if(hl) continue; }
    if(only) continue;
    px.push(X[b+ax[0]]);py.push(X[b+ax[1]]);pz.push(X[b+ax[2]]); cd.push([i,v.S[v.sid[i]].label,sc[i]]);
    if(cm==="structure") pc.push(v.S[v.sid[i]].colour); else if(cm==="dim") pc.push(DIMCOL[Math.min(v.S[v.sid[i]].dim,8)]); else pc.push(v.scores[cm][i]);
  }
  const traces=[];
  const ht="row %{customdata[0]} · %{customdata[1]}<br>score %{customdata[2]:.4g}<extra></extra>";
  if($("showpts").checked && px.length){
    const marker={size:+$("psize").value, opacity:+$("popac").value, color:pc, line:{width:0}};
    if(cm!=="structure" && cm!=="dim"){ const rg=seqColour(pc); Object.assign(marker,{colorscale:"Blues",reversescale:true,cmin:rg.cmin,cmax:rg.cmax,
      colorbar:{title:{text:DATA.scores.find(s=>s.key===cm).label,side:"right"},thickness:12,len:0.6}}); }
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
$("version").onchange=()=>{ cur=+$("version").value; buildTree(); setThrRange(); draw(); };
["ax0","ax1","ax2","colour","psize","popac","showpts","showgraphs","scope","hl","only"].forEach(id=>$(id).addEventListener("input",draw));
$("score").onchange=()=>{ setThrRange(); draw(); };
$("pct").oninput=()=>{ $("pctnum").value=$("pct").value; document.querySelector("input[name=mode][value=pct]").checked=true; draw(); };
$("pctnum").oninput=()=>{ $("pct").value=$("pctnum").value; document.querySelector("input[name=mode][value=pct]").checked=true; draw(); };
$("thr").oninput=()=>{ $("thrnum").value=(+$("thr").value).toPrecision(4); document.querySelector("input[name=mode][value=thr]").checked=true; draw(); };
$("thrnum").oninput=()=>{ $("thr").value=$("thrnum").value; document.querySelector("input[name=mode][value=thr]").checked=true; draw(); };
document.querySelectorAll("input[name=mode]").forEach(e=>e.addEventListener("change",draw));
$("csv").onclick=()=>download(true); $("csvall").onclick=()=>download(false);
buildTree(); setThrRange(); draw();
</script></body></html>
"""
