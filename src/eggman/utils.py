import numpy as np

from .cy_eggman import BiellipsoidWrap, LightSourceWrap, OrbitWrap


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
        limb_type,
        limb,
        eccen=0,
        lon_periapse=0,
        rotate_with_orbit=True,
        area_args=dict(),
        meridian_args=dict(visible_only=True),
        source_args=dict(),
):
    from matplotlib import pyplot as plt
    _, ax = plt.subplots(figsize=(12, 8))

    orb = OrbitWrap(period, t0, semimajor, eccen, inclination, lon_periapse)
    assert (r_forward > 0) and (r_back > 0) and (r_side > 0), "Radii must be positive."
    ci = orb.cos_inc if rotate_with_orbit else -2.
    x, y, z = orb.get_position(t)

    # Draw the planet
    # TODO: Handle catwoman emulation
    bell = BiellipsoidWrap(x, y, z, theta, phi, gamma, r_forward, r_back, r_pole, r_side, ci)
    bell.plot_area(**area_args)
    bell.plot_meridians(**meridian_args)

    # Draw the parent star
    source = LightSourceWrap(limb_type, limb)
    source_args.setdefault('zorder', -100 if z > 0 else 100)
    source.plot_brightness(pcm_args=source_args)
