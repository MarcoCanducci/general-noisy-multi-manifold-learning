"""Renderer-independent 3-D scene description and its two browser backends.

``plotly`` draws the scene with WebGL in the browser (interactive, nothing rendered on the
server, works everywhere). ``pyvista`` renders the scene with VTK on the server and shows the
image, with sliders for the camera; it needs an off-screen OpenGL library (a display, EGL or
OSMesa) and is meant for very large scenes or publication-quality views.
"""
from __future__ import annotations

import subprocess
import sys
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np

__all__ = ["Scene", "pyvista_available", "render", "to_xyz"]


def to_xyz(P) -> np.ndarray:
    """Pad 2-column coordinates with z = 0."""
    P = np.asarray(P, float)
    if P.shape[1] == 2:
        P = np.c_[P, np.zeros(len(P))]
    return np.ascontiguousarray(P[:, :3])


@dataclass
class Scene:
    """Ordered list of layers: ("points", name, xyz, colors, size, opacity) and
    ("lines", name, nodes, edges, color, width). ``colors`` is a hex string or an (n, 3) uint8
    array; line ``color`` is a hex string."""

    layers: list = field(default_factory=list)
    labels: tuple = ("x", "y", "z")
    flat: bool = False

    def points(self, name, xyz, color="#2a78d6", size=2.0, opacity=1.0):
        xyz = to_xyz(xyz)
        if len(xyz):
            self.layers.append(("points", name, xyz, color, float(size), float(opacity)))
        return self

    def lines(self, name, nodes, edges, color="#0b0b0b", width=2.0):
        nodes, edges = to_xyz(nodes), np.asarray(edges, int).reshape(-1, 2)
        if len(edges):
            self.layers.append(("lines", name, nodes, edges, color, float(width)))
        return self

    def sphere(self, name, centre, radius, color="#e34948"):
        """Three great circles marking a ball of the given radius."""
        t = np.linspace(0, 2 * np.pi, 49)
        c = np.asarray(to_xyz(np.atleast_2d(centre))[0])
        z = np.zeros_like(t)
        rings = [np.c_[np.cos(t), np.sin(t), z], np.c_[np.cos(t), z, np.sin(t)], np.c_[z, np.cos(t), np.sin(t)]]
        if self.flat:
            rings = rings[:1]
        nodes = np.vstack([c + radius * r for r in rings])
        edges = np.vstack([np.c_[np.arange(48), np.arange(1, 49)] + 49 * k for k in range(len(rings))])
        return self.lines(name, nodes, edges, color, 3.0)


@lru_cache(maxsize=1)
def pyvista_available() -> bool:
    """True if PyVista can render off-screen here (probed in a subprocess, because a missing
    OpenGL library crashes the interpreter instead of raising)."""
    code = ("import numpy as np, pyvista as pv\n"
            "pl = pv.Plotter(off_screen=True)\n"
            "pl.add_points(pv.PolyData(np.random.rand(10, 3)))\n"
            "pl.render()\n")
    try:
        r = subprocess.run([sys.executable, "-c", code], capture_output=True, timeout=60)
        return r.returncode == 0
    except Exception:
        return False


def _hex_to_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) for i in (0, 2, 4))


def _plotly(scene: Scene, height: int):
    import plotly.graph_objects as go

    traces = []
    for layer in scene.layers:
        if layer[0] == "points":
            _, name, xyz, color, size, opacity = layer
            if not isinstance(color, str):
                color = [f"rgb({r},{g},{b})" for r, g, b in np.asarray(color)]
            traces.append(go.Scatter3d(x=xyz[:, 0], y=xyz[:, 1], z=xyz[:, 2], mode="markers", name=name,
                                       marker=dict(size=size, color=color, opacity=opacity, line=dict(width=0)),
                                       hoverinfo="skip"))
        else:
            _, name, nodes, edges, color, width = layer
            seg = np.full((len(edges) * 3, 3), np.nan)
            seg[0::3], seg[1::3] = nodes[edges[:, 0]], nodes[edges[:, 1]]
            traces.append(go.Scatter3d(x=seg[:, 0], y=seg[:, 1], z=seg[:, 2], mode="lines", name=name,
                                       line=dict(color=color, width=width), hoverinfo="skip",
                                       connectgaps=False))
    fig = go.Figure(traces)
    axis = dict(backgroundcolor="white", gridcolor="#dddcd6", zerolinecolor="#dddcd6", color="#52514e")
    cam = dict(eye=dict(x=0, y=0, z=2.2), up=dict(x=0, y=1, z=0)) if scene.flat else None
    fig.update_layout(height=height, margin=dict(l=0, r=0, t=0, b=0), paper_bgcolor="white",
                      legend=dict(orientation="h", y=0.99, x=0.01, font=dict(size=11), itemsizing="constant"),
                      scene=dict(aspectmode="data", camera=cam,
                                 xaxis=dict(title=scene.labels[0], **axis), yaxis=dict(title=scene.labels[1], **axis),
                                 zaxis=dict(title=scene.labels[2] if len(scene.labels) > 2 else "", **axis)),
                      uirevision="keep")
    return fig


def _pyvista(scene: Scene):
    import pyvista as pv

    pl = pv.Plotter(off_screen=True)
    pl.set_background("white")
    for layer in scene.layers:
        if layer[0] == "points":
            _, name, xyz, color, size, opacity = layer
            cloud = pv.PolyData(xyz)
            if isinstance(color, str):
                pl.add_points(cloud, color=color, point_size=size + 1, opacity=opacity)
            else:
                cloud["rgb"] = np.asarray(color, np.uint8)
                pl.add_points(cloud, scalars="rgb", rgb=True, point_size=size + 1, opacity=opacity)
        else:
            _, name, nodes, edges, color, width = layer
            mesh = pv.PolyData(nodes)
            mesh.lines = np.c_[np.full(len(edges), 2), edges].ravel()
            pl.add_mesh(mesh, color=color, line_width=width)
    if scene.flat:
        pl.view_xy()
    else:
        pl.view_isometric()
    return pl


def render(scene: Scene, backend: str = "plotly", height: int = 640, pane=None, camera=(35.0, 25.0, 1.0)):
    """Return ``(panel pane, backend used)`` for the scene. An existing pane of the same backend
    can be passed as ``pane`` to be updated in place. ``camera`` = (azimuth, elevation, zoom) is
    used by the pyvista backend."""
    import panel as pn

    if backend == "auto":
        backend = "plotly"
    if backend == "pyvista":
        from io import BytesIO

        from PIL import Image

        pl = _pyvista(scene)
        pl.reset_camera()
        if not scene.flat:
            pl.camera.azimuth = camera[0] - 45.0
            pl.camera.elevation = camera[1] - 35.264
        pl.camera.zoom(camera[2])
        img = pl.screenshot(return_img=True, window_size=(1100, height))
        pl.close()
        buf = BytesIO()
        Image.fromarray(img).save(buf, format="PNG")
        if isinstance(pane, pn.pane.PNG):
            pane.object = buf.getvalue()
            return pane, backend
        return pn.pane.PNG(buf.getvalue(), sizing_mode="scale_width", max_width=1100), backend
    fig = _plotly(scene, height)
    if isinstance(pane, pn.pane.Plotly):
        pane.object = fig
        return pane, "plotly"
    return pn.pane.Plotly(fig, height=height, sizing_mode="stretch_width", config={"displaylogo": False}), "plotly"
