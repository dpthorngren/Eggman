import numpy as np

import eggman


def plot_biellipse(bell, res=200):
    from matplotlib import pyplot as plt
    _, ax = plt.subplots(figsize=(12, 8))

    # Principle axes
    x, y, _ = bell["position"]
    axes = eggman.utils.biellipsoid_axes(bell)
    for i in range(4):
        color = ['red', 'green', 'blue', "orange"][i]
        ax.plot([x, x + axes[0, i]], [y, y + axes[1, i]], zorder=1, color=color, alpha=.3)

    # Limb outline
    xyz = biellipse_outline(bell, res)
    forward, _ = position_masks(bell, xyz)
    for mask, color in zip([forward, ~forward], ['red', 'green']):
        x = np.ma.masked_where(mask, xyz[0])
        y = np.ma.masked_where(mask, xyz[1])
        ax.plot(x, y, color=color, zorder=10)

    # Meridians and equator
    for i in range(3):
        xyz = biellipsoid_meridians(bell, i, res)
        _, visible = position_masks(bell, xyz)
        # x = np.ma.masked_where(~visible, xyz[0])
        # y = np.ma.masked_where(~visible, xyz[1])
        # ax.plot(x, y, ':', color='black', zorder=10, alpha=.1)
        x = np.ma.masked_where(visible, xyz[0])
        y = np.ma.masked_where(visible, xyz[1])
        ax.plot(x, y, color='black', zorder=10, alpha=.5)


def position_masks(bell, positions):
    '''For positions on the biellipsoid surface, returns whether they are visible
    and whether they are on the front or back half of the biellipse.'''
    positions = positions - bell['position'][:, None]
    dir_forward = np.array(bell['rot'][::3])
    forward = dir_forward.dot(positions) * dir_forward.dot([0., 0., -1]) >= 0.
    visible = np.zeros_like(forward)
    f_limb = np.cross(bell['f1'], bell['f2'])
    f_limb /= np.linalg.norm(f_limb) * np.sign(f_limb[2])
    b_limb = np.cross(bell['b1'], bell['b2'])
    b_limb /= np.linalg.norm(b_limb) * np.sign(b_limb[2])
    visible[forward] = (np.dot(f_limb, positions) <= 0.)[forward]
    visible[~forward] = (np.dot(b_limb, positions) <= 0.)[~forward]
    return forward, visible


def biellipse_outline(bell, res=200) -> np.ndarray:
    forward = np.array(bell['rot'][::3])
    t = np.linspace(0, 2 * np.pi, res)[:, None]
    e_forward = np.cos(t) * np.array(bell['f1']) + np.sin(t) * np.array(bell['f2'])
    mask = e_forward @ forward > 0
    e_forward = e_forward[mask]
    e_backward = np.cos(t) * np.array(bell['b1']) + np.sin(t) * np.array(bell['b2'])
    mask = e_backward @ forward > 0
    e_backward = e_backward[~mask]
    result = np.vstack([e_forward, e_backward])
    result = sorted(result, key=lambda x: np.arctan2(x[1], x[0]))
    result = np.vstack([result, result[0]])
    return result.T + np.array(bell['position'])[:, None]


def biellipsoid_axes(bell):
    '''Get the principle axes from a biellipsoid, for visualization and testing purposes.'''
    rot = np.array(bell['rot']).reshape(3, 3)
    radii = bell['radii']
    forward = rot @ np.diag([radii[0], radii[2], radii[3]])
    back = rot @ np.diag(radii[1:])
    back[:, 0] *= -1
    return np.column_stack([forward[:, 0], back])


def biellipsoid_meridians(bell, axis=0, res=200):
    axes = biellipsoid_axes(bell)
    forward = axes[:, 0]
    t = np.linspace(0, 2 * np.pi, res)[:, None]
    if axis == 0:
        v = np.cos(t) * axes[:, 2] + np.sin(t) * axes[:, 3]
        return v.T + bell['position'][:, None]
    elif axis == 1:
        u = np.cos(t) * axes[:, 0] + np.sin(t) * axes[:, 3]
        v = np.cos(t) * axes[:, 1] + np.sin(t) * axes[:, 3]
    elif axis == 2:
        u = np.cos(t) * axes[:, 0] + np.sin(t) * axes[:, 2]
        v = np.cos(t) * axes[:, 1] + np.sin(t) * axes[:, 2]
    else:
        raise ValueError(f"Axis must be 0, 1, or 2, not {axis}.")
    u = u[np.dot(u, forward) > 0]
    v = v[np.dot(v, forward) < 0]
    result = np.vstack([u, v])
    result = sorted(result, key=lambda x: np.arctan2(x[1], x[0]))
    result = np.vstack([result, result[0]])
    return np.array(result).T + bell['position'][:, None]


def _decompose_position_(x, y, z):
    '''Determines radius, sin and cosine of the azimuthal angle, and sin and cosine of
    the inclination angle consistent with the given position. These are only orbital
    elements if the orbit is circular. There are multiple solutions if z=0, assumes
    i=90 in that case. Used mainly for creating valid test cases.'''
    rad = np.sqrt(x*x + y*y + z*z)
    st = x / rad
    ct = np.sign(z) * np.sqrt(1 - st*st)
    ci = z / (rad*ct) if ct != 0 else 0.
    si = np.sqrt(1 - ci*ci)
    return rad, ct, st, ci, si
