import numpy as np


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
