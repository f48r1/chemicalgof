# [ ] latest update module: any improvements ?

import math
import random
from collections.abc import Callable

import networkx as nx

def force_layout(
    graph: nx.Graph,
    sizes: dict,                    # node -> (width, height)
    *,
    iterations: int = 300,
    k: float | None = None,         # ideal spring length (auto if None)
    gravity: float = 0.05,          # pull toward centre of mass
    repulsion_strength: float = 1.2,
    overlap_padding: float = 8.0,   # extra gap to maintain between boxes
    initial_positions: dict | None = None,
    seed: int | None = None,
    progress_cb: Callable[[int, float], None] | None = None,
) -> dict:
    """
    Compute non-overlapping positions for rectangle nodes in *graph*.

    Parameters
    ----------
    graph : nx.Graph
        The graph whose nodes need to be laid out.
    sizes : dict
        Mapping ``node -> (width, height)`` for every node.
    iterations : int
        Maximum number of simulation steps.
    k : float, optional
        Natural spring length. Defaults to ``sqrt(area / n)`` where *area*
        is estimated from the total node area.
    gravity : float
        Strength of the centre-gravity force that prevents the layout from
        drifting off to infinity.
    repulsion_strength : float
        Multiplier for the node-to-node repulsion force.
    overlap_padding : float
        Minimum pixel gap enforced between any two rectangles.
    initial_positions : dict, optional
        Seed positions ``node -> (x, y)`` (top-left corner).  Random if omitted.
    seed : int, optional
        Random seed for reproducibility.
    progress_cb : callable, optional
        Called as ``progress_cb(iteration, temperature)`` after each step.

    Returns
    -------
    dict
        ``node -> (x, y)`` giving the **centre** of each rectangle.
    """
    if len(graph) == 0:
        return {}

    rng = random.Random(seed)
    nodes = list(graph.nodes())
    n = len(nodes)

    # --- default spring length based on average node size ----------------
    avg_w = sum(sizes[v][0] for v in nodes) / n
    avg_h = sum(sizes[v][1] for v in nodes) / n
    if k is None:
        # pack area * 4 so nodes have breathing room
        total_area = sum(sizes[v][0] * sizes[v][1] for v in nodes)
        k = math.sqrt(total_area * 4.0 / max(n, 1))

    # --- initial positions (centres) -------------------------------------
    spread = k * math.sqrt(n)
    if initial_positions:
        pos = {
            v: list(initial_positions[v])
            for v in nodes
        }
    else:
        pos = {
            v: [rng.uniform(-spread, spread), rng.uniform(-spread, spread)]
            for v in nodes
        }

    # --- cooling schedule ------------------------------------------------
    # temperature starts high (free movement) and cools to near zero
    t = spread * 0.5
    t_min = k * 0.01
    cooling = (t_min / t) ** (1.0 / max(iterations, 1))

    # --- main loop -------------------------------------------------------
    for iteration in range(iterations):
        disp = {v: [0.0, 0.0] for v in nodes}

        # ----- repulsion (rectangle-aware) -------------------------------
        _apply_repulsion(nodes, pos, sizes, disp,
                         k, repulsion_strength, overlap_padding)

        # ----- attraction along edges ------------------------------------
        _apply_attraction(graph, nodes, pos, sizes, disp, k)

        # ----- gravity toward barycentre ---------------------------------
        cx = sum(pos[v][0] for v in nodes) / n
        cy = sum(pos[v][1] for v in nodes) / n
        for v in nodes:
            disp[v][0] -= gravity * (pos[v][0] - cx)
            disp[v][1] -= gravity * (pos[v][1] - cy)

        # ----- cap displacement by temperature & apply -------------------
        for v in nodes:
            dx, dy = disp[v]
            length = math.hypot(dx, dy) or 1e-9
            scale = min(length, t) / length
            pos[v][0] += dx * scale
            pos[v][1] += dy * scale

        t = max(t * cooling, t_min)

        if progress_cb:
            progress_cb(iteration, t)

    # --- post-pass: hard-push any remaining overlaps ---------------------
    _resolve_overlaps(nodes, pos, sizes, overlap_padding, passes=20)

    return {v: tuple(pos[v]) for v in nodes}


# ---------------------------------------------------------------------------
# Force helpers
# ---------------------------------------------------------------------------

def _rect_gap(ax, ay, aw, ah, bx, by, bw, bh):
    """
    Return (dx, dy, gap) where (dx, dy) is the vector from B's centre to
    A's centre and *gap* is the signed minimum separation between the two
    rectangles (negative means overlap).
    """
    dx = ax - bx
    dy = ay - by
    # half-extents
    ahw, ahh = aw / 2.0, ah / 2.0
    bhw, bhh = bw / 2.0, bh / 2.0

    # overlap amount on each axis (positive = overlap)
    ox = (ahw + bhw) - abs(dx)
    oy = (ahh + bhh) - abs(dy)
    gap = -min(ox, oy)          # negative when overlapping
    return dx, dy, gap


def _apply_repulsion(nodes, pos, sizes, disp,
                     k, strength, padding):
    """O(n²) rectangle repulsion – fast enough up to ~500 nodes."""
    k2 = k * k * strength
    for i, u in enumerate(nodes):
        ux, uy = pos[u]
        uw, uh = sizes[u]
        for v in nodes[i + 1:]:
            vx, vy = pos[v]
            vw, vh = sizes[v]

            dx, dy, gap = _rect_gap(ux, uy, uw, uh, vx, vy, vw, vh)
            dist = math.hypot(dx, dy) or 1e-9

            if gap < padding:
                # strong repulsion when close / overlapping
                needed = padding - gap          # how much to push
                force = k2 / dist + needed * 2.0
            else:
                # mild long-range repulsion (Fruchterman-Reingold style)
                force = k2 / dist

            fx = (dx / dist) * force
            fy = (dy / dist) * force
            disp[u][0] += fx
            disp[u][1] += fy
            disp[v][0] -= fx
            disp[v][1] -= fy


def _apply_attraction(graph, nodes, pos, sizes, disp, k):
    """Spring attraction: anchored at the closest border point of each rect."""
    for u, v in graph.to_undirected().edges():
        if u not in pos or v not in pos:
            continue
        ux, uy = pos[u]
        vx, vy = pos[v]
        uw, uh = sizes[u]
        vw, vh = sizes[v]

        # anchor points on the rectangle borders facing each other
        ax, ay = _border_point(ux, uy, uw, uh, vx, vy)
        bx, by = _border_point(vx, vy, vw, vh, ux, uy)

        dx = bx - ax
        dy = by - ay
        dist = math.hypot(dx, dy) or 1e-9

        # Hooke's law: F = dist² / k  (Fruchterman-Reingold attraction)
        force = dist * dist / k
        fx = (dx / dist) * force
        fy = (dy / dist) * force

        disp[u][0] += fx
        disp[u][1] += fy
        disp[v][0] -= fx
        disp[v][1] -= fy


def _border_point(cx, cy, w, h, tx, ty):
    """
    Return the point on the border of the rectangle (cx,cy,w,h) that lies
    on the line segment from the rectangle centre toward target (tx,ty).
    """
    dx = tx - cx
    dy = ty - cy
    if dx == 0 and dy == 0:
        return cx, cy
    hw, hh = w / 2.0, h / 2.0
    # clip to rectangle border
    if abs(dy) * hw <= abs(dx) * hh:
        # intersects left or right side
        scale = hw / abs(dx)
    else:
        # intersects top or bottom side
        scale = hh / abs(dy)
    return cx + dx * scale, cy + dy * scale


# ---------------------------------------------------------------------------
# Post-process overlap resolver
# ---------------------------------------------------------------------------

def _resolve_overlaps(nodes, pos, sizes, padding, passes=10):
    """
    Iterative hard-push pass: push overlapping rectangles apart until
    no two rectangles overlap (plus *padding*).  Much faster than the main
    simulation for the final pixel-perfect cleanup.
    """
    for _ in range(passes):
        moved = False
        for i, u in enumerate(nodes):
            ux, uy = pos[u]
            uw, uh = sizes[u]
            for v in nodes[i + 1:]:
                vx, vy = pos[v]
                vw, vh = sizes[v]

                # how much do they overlap?
                ox = (uw + vw) / 2.0 + padding - abs(ux - vx)
                oy = (uh + vh) / 2.0 + padding - abs(uy - vy)

                if ox > 0 and oy > 0:
                    # push along the axis of minimum overlap
                    if ox < oy:
                        push = ox / 2.0 + 0.5
                        sign = 1 if ux >= vx else -1
                        pos[u][0] += sign * push
                        pos[v][0] -= sign * push
                    else:
                        push = oy / 2.0 + 0.5
                        sign = 1 if uy >= vy else -1
                        pos[u][1] += sign * push
                        pos[v][1] -= sign * push
                    moved = True
        if not moved:
            break