from ..gof import DiGraphFrags, FragNode
from .layout import force_layout
from .nodes import drawNode

def _rect_border_point(center, size, direction):
    import numpy as np
    w, h = size
    dx, dy = direction

    tx = (w / 2) / abs(dx) if dx != 0 else np.inf
    ty = (h / 2) / abs(dy) if dy != 0 else np.inf

    t = min(tx, ty)

    return center + direction * t

def drawGoF(
    graph : DiGraphFrags,
    random_seed : int | None = None,
    vert_or_horiz:str='horiz',
    dpi=100,
    custom_positions : dict[FragNode,tuple[float,float]] | None = None,
    custom_node_images : dict[FragNode,tuple[float,float]] | None = None,
):
    
    import numpy as np
    from matplotlib import pyplot as plt

    nodes :list[FragNode] = list(graph.nodes)
    edges:list[tuple[FragNode, FragNode]] = list(graph.to_undirected().edges)

    if custom_node_images is not None:
        node_image = custom_node_images
    else:
        node_image = {node:drawNode(node) for node in nodes}

    node_size = {node:img.size for node,img in node_image.items()}

    if custom_positions is not None:
        positions = custom_positions
    else:
        positions = force_layout(graph, node_size, k=60, overlap_padding=1., seed=random_seed)

    pos_arr = np.array([ tuple(positions[node]) for node in nodes])
    sizes_arr = np.array([ tuple(node_size[node]) for node in nodes ])

    xy_min = np.min(pos_arr - sizes_arr / 2, axis=0)
    xy_max = np.max(pos_arr + sizes_arr / 2, axis=0)

    width, height = xy_max - xy_min

    if (vert_or_horiz == 'horiz' and height > width) or (vert_or_horiz == 'vert' and width > height):
        positions = {node:(y,x) for node,(x,y) in positions.items()}
        width, height = height, width
        xy_min = xy_min[::-1]
        xy_max = xy_max[::-1]

    fig = plt.figure(figsize=(width/dpi, height/dpi), dpi=dpi)
    ax:plt.Axes = fig.add_axes([0, 0, 1, 1])  # occupa tutta la figura

    ax.axis("off")
    ax.set_xlim(xy_min[0], xy_max[0])
    ax.set_ylim(xy_min[1], xy_max[1])
    # ax.set_aspect("equal")

    for node in nodes:
        cx, cy = positions[node]
        w_, h_ = node_size[node]
        hw, hh = w_/2, h_/2

        ax.imshow(node_image[node],
                    extent=(cx-hw, cx+hw, cy-hh, cy+hh),
                    zorder=3)

    circle_radius = 10
    circle_size = 3.14 * circle_radius**2
    for src_and_dst in edges:
        xy_nodes = np.array([ positions[node] for node in src_and_dst ])

        diff = xy_nodes[1] - xy_nodes[0]
        dist = np.linalg.norm(diff)
        direction = diff / dist

        border_padded_nodes = np.array([
            _rect_border_point(positions[node], node_size[node], direction * factor_direction ) + direction * factor_direction * circle_radius
            for factor_direction,node in zip([1,-1], src_and_dst)
        ])

        ax.plot(*border_padded_nodes.T, '-', zorder=1, c='black')

        for idx_node in range(2):
            src = src_and_dst[idx_node]
            if src.fragment.num_connector <= 1:
                continue

            dst = src_and_dst[(idx_node - 1) * (-1)]

            edge_data = graph[src][dst]
            label = str(edge_data['aB'])
            if edge_data['stereo'] is not None:
                label+=edge_data['stereo']

            center_circle = border_padded_nodes[idx_node]

            ax.scatter(*center_circle, fc='white', marker='o', s=circle_size, zorder=4, ec='black')
            ax.text(*center_circle, label, ha='center', va='center', zorder=5, fontdict={'fontsize':8})


    plt.close(fig)
    return fig
