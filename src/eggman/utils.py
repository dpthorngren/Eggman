import numpy as np

from .cy_eggman import BiellipsoidWrap, OrbitWrap


def rotation_matrix(angle, axis, degrees=True):
    '''Creates a matrix representing a rotation counterclockwise around axis by angle.'''
    if degrees:
        angle *= np.pi / 180.
    ci = np.cos(angle)
    si = np.sin(angle)
    rot = np.array([[1., 0., 0.], [0, ci, -si], [0, si, ci]])
    rot = np.roll(rot, axis, axis=0)
    rot = np.roll(rot, axis, axis=1)
    return rot


def plot_transit(
        r_forward,
        r_back,
        r_pole,
        r_side,
        theta,
        phi,
        gamma,
        t,
        t0,
        period,
        semimajor,
        inclination,
        eccen=0,
        lon_periapse=0,
        rotate_with_orbit=True,
        outline_args=dict(),
        area_args=dict(),
        meridian_args=dict(),
):
    from matplotlib import pyplot as plt
    _, ax = plt.subplots(figsize=(12, 8))

    orb = OrbitWrap(period, t0, semimajor, eccen, inclination, lon_periapse)
    assert (r_forward > 0) and (r_back > 0) and (r_side > 0), "Radii must be positive."
    ci = orb.cos_inc if rotate_with_orbit else -2.
    x, y, z = orb.get_position(t)
    bell = BiellipsoidWrap(x, y, z, theta, phi, gamma, r_forward, r_back, r_pole, r_side, ci)
    print(bell.position)
    # bell.plot_outline()
    bell.plot_area(**area_args)
    # bell.plot_meridians()

    # TODO: Better source plotting
    x = np.linspace(-1, 1, 200)
    y = np.sqrt(1 - x**2)
    zorder = -100 if z > 0 else 100
    plt.fill_between(x, -y, y, zorder=zorder, color="silver")
