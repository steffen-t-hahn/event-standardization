"""."""
import numpy as np
import matplotlib.pyplot as plt

from itertools import product
from functools import reduce


def get_transformation(mat, azimuth, cycle=6, max_angle=2 * np.pi,
                       sp_trafo=np.eye(2), inv_sp_trafo=None):
    """Compute the transformation matrix for the given azimuth.

    Example implementation of computing the right transformation matrix
    for a given cyclic transfromation.

    Args:
        mat: basic cyclic transformation matrix (2, 2)
        azimuth: event angle

    Kwargs:
        cycle: (mat ** (cycle + 1)) = mat
        max_angle: maximum angle we bin the azimuth
        sp_trafo: transformation into integral coordinates of the grid
        inv_sp_trafo: inverse sp_trafo

    """
    bins = np.linspace(np.deg2rad(30), max_angle - np.deg2rad(30),
                       num=cycle)

    # 0 is left of bins[0] => therefore -1
    n = np.digitize(azimuth, bins=bins) - 1
    n = n if n >= 0 else cycle - 1

    trafo = reduce(lambda l, r: l @ r, [np.eye(2)] + [mat] * n)

    if inv_sp_trafo is None:
        inv_sp_trafo = np.linalg.inv(sp_trafo)

    azimuth_new = np.arccos(((inv_sp_trafo @ trafo) @ sp_trafo @
                             [np.cos(azimuth), np.sin(azimuth)])[0])

    return trafo, azimuth_new


def hex_grid(ax, coordinates, radius=1. / np.sqrt(3),
             polygon_kwargs={'edgecolor': 'k', 'joinstyle': 'round'},
             threshold=1e-9):
    """Plot a hexagonal grid on ax."""
    from matplotlib.patches import RegularPolygon

    # (0) basic variables
    X_u, Y_u = np.array((1, 0)), np.array((-0.5, np.sqrt(3) / 2))

    x_l = 0.5 * np.array([-1, -1, 0, 1, 1, 0])
    y_l = 0.5 / np.sqrt(3) * np.array([-1, 1, 2, 1, -1, -2])

    for c in coordinates:
        x_c, y_c = c[0] * X_u + c[1] * Y_u

        ax.fill(x_l + x_c, y_l + y_c,
                facecolor='red', **polygon_kwargs)

    for c in product(np.arange(-10, 7), np.arange(-6, 10)):
        x_c, y_c = c[0] * X_u + c[1] * Y_u

        if x_c**2 + y_c**2 > 36:
            continue

        hex_patch = RegularPolygon(
            (x_c, y_c), numVertices=6, radius=radius,
            orientation=np.radians(0),
            facecolor='none', lw=0.2, **polygon_kwargs)
        ax.add_patch(hex_patch)


def plot_transformation(old_coordinates, new_coodinates,
                        old_azimuth, new_azimuth):
    """."""
    fig, axes = plt.subplots(nrows=1, ncols=2, dpi=300)

    axes[0].set_ylabel('y')
    for ax, coords, azimuth in zip(axes.ravel(),
                                   [old_coordinates, new_coodinates],
                                   [old_azimuth, new_azimuth]):
        plt.setp(ax, aspect='equal', xlim=(-6, 6), ylim=(-6, 6),
                 xticks=[], yticks=[], xlabel='x')
        hex_grid(ax, coordinates=coords)

        ax.arrow(+0.5 * np.cos(azimuth) + 4.8, +0.5 * np.sin(azimuth) - 4.8,
                 -1 * np.cos(azimuth), -1 * np.sin(azimuth),
                 head_width=0.4, head_length=0.3)

        # write angle phi
        ax.text(0.1, 0.92,
                "${:0.1f}^\\circ$".format(np.rad2deg(azimuth)),
                ha='left', va='top',
                transform=ax.transAxes,
                bbox=dict(facecolor='white', alpha=0.9))

    fig.tight_layout()
    plt.show()

    # fig.savefig('output/transformation_check.png')


if __name__ == '__main__':
    # (1) define mock data
    relativ_position = [[0, 0], [1, 0], [1, 1], [1, 2], [3, 3],
                        [1, 3], [2, 2], [-1, 0], [0, -1], [-2, 2]]
    relativ_position = np.array(relativ_position)

    # (2) compute right transformation matrix
    azimuth_real = np.deg2rad(157)
    sp_trafo = np.array([[1, 1. / np.sqrt(3)],
                         [0, 2. / np.sqrt(3)]])

    # (2.1) compute rotation matrix
    trafo_rot, azimuth_rot = \
        get_transformation(np.array([[0, 1], [-1, 1]]), azimuth_real,
                           sp_trafo=sp_trafo)

    # (2.2) compute reflection matrix
    trafo_reflect, azimuth_final = \
        get_transformation(np.array([[0, 1], [1, 0]]), azimuth_rot,
                           cycle=2, max_angle=np.deg2rad(90),
                           sp_trafo=sp_trafo)

    # (3) compute and plot final position
    final_positon = (trafo_reflect @ (trafo_rot @ relativ_position.T)).T

    plot_transformation(relativ_position, final_positon,
                        azimuth_real, azimuth_final)
