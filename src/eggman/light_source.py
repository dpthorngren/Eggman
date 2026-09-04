import typing
from enum import Enum

import cython
import numpy as np
# Cython cimports (LSP cannot process)
from cython.cimports.eggman import cy_eggman as cye  # type: ignore

if not cython.compiled and typing.TYPE_CHECKING:
    # Allow editor LSP to find existing compiled code...
    SourceType = Enum(
        "SourceType",
        ["NoEmission", "Lambertian", "QuadraticLimb", "NonLinearLimb", "DayNight", "EmissionMap"])
    from eggman import Shape
else:  # ...but provide the real cimports during compilation
    from cython.cimports.eggman.cy_eggman import Shape, SourceType  # type: ignore


@cython.cclass
class LightSource:
    csource: cye.LightSource

    _source_type_names_ = [
        'no_emission', 'lambertian', 'quadratic_limb', 'nonlinear_limb', 'day_night', 'emission_map'
    ]
    _source_n_params_ = [0, 1, 3, 5, 2, 2]
    _source_enum_ = [
        SourceType.NoEmission,
        cye.SourceType.Lambertian,
        cye.SourceType.QuadraticLimb,
        cye.SourceType.NonLinearLimb,
        cye.SourceType.DayNight,
        cye.SourceType.EmissionMap,
    ]

    @property
    def source_type(self):
        return self._source_type_names_[self.csource.stype]

    @property
    def source_type_code(self):
        return self.csource.stype

    @property
    def source_params(self):
        return np.array(self.csource.params)[:self._source_n_params_[self.source_type_code]]

    @property
    def limb_norm(self):
        return self.csource.limb_norm

    @property
    def n(self):
        return self.csource.get_n()

    @property
    def m(self):
        return self.csource.get_m()

    @property
    def map_size(self):
        return self.csource.get_map_size()

    def __init__(self, source_type: str, source_params: list[float]):
        source_type = source_type.lower().strip()
        assert source_type in self._source_type_names_, \
            f"Error: unknown source_type {source_type}, must be one of {self._source_type_names_}"

        # Check the number of parameters given matches the required number
        source_code: cye.SourceType = self._source_enum_[self._source_type_names_.index(
            source_type)]
        n_params = len(source_params)
        expect_n_params = self._source_n_params_[source_code]
        if source_type in ['quadratic_limb', 'nonlinear_limb'] and n_params == expect_n_params - 1:
            source_params.insert(0, 1 / np.pi)
            n_params += 1
        assert len(source_params) == expect_n_params, \
            f"Error: wrong number of source parameters; should be {expect_n_params},"\
            f"was {n_params}."

        # Create a c-readable form of the source parameters
        source_params_c: cye.Array_SourceParams = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
        for i, p in enumerate(source_params):
            source_params_c[i] = p
        self.csource = cye.CLightSource(source_code, source_params_c)

    def get_brightness(self, x, y, bell: Shape, ring_front: cython.bint = True):
        cython.declare(i=int, j=int)

        types = [float, np.float64, np.float32, np.int32, np.int64, int]
        if type(x) in types and type(y) in types:
            if bell.r_up == 0:
                loc = bell.raycast(x, y)[1]
                loc = bell.sphere_to_world(loc)
                if ring_front != loc[2] > 0.:
                    return 0.
            return self.csource.get_brightness(x, y, bell.cshape)
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        assert x.ndim == 1 and y.ndim == 1

        output = np.zeros((x.shape[0], y.shape[0]), dtype=np.double)
        x_view: cye.Array1d_f64 = x
        y_view: cye.Array1d_f64 = y
        output_view: cye.Array2d_f64 = output

        for i in range(x_view.shape[0]):
            for j in range(y_view.shape[0]):
                if bell.r_up == 0:
                    loc = bell.raycast(x_view[i], y_view[j])[1]
                    loc = bell.sphere_to_world(loc)
                    if ring_front != loc[2] > 0.:
                        output_view[i, j] = 0.
                        continue
                output_view[i, j] = self.csource.get_brightness(x_view[i], y_view[j], bell.cshape)
        return output

    def get_brightness_sphere(self, x, y):
        cython.declare(i=int, j=int)
        types = [float, np.float64, np.float32, np.int32, np.int64, int]
        if type(x) in types and type(y) in types:
            return self.csource.get_brightness_sphere(x, y)
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        assert x.ndim == 1 and y.ndim == 1

        output = np.zeros((len(x), len(y)), dtype=np.double)
        x_view: cye.Array1d_f64 = x
        y_view: cye.Array1d_f64 = y
        output_view: cye.Array2d_f64 = output

        for i in range(x_view.shape[0]):
            for j in range(y_view.shape[0]):
                output_view[i, j] = self.csource.get_brightness_sphere(x_view[i], y_view[j])
        return output

    def get_emission_point(self, i: int) -> float:
        return self.csource.get_emission_point(i)

    def set_emission_point(self, i: int, value: float) -> None:
        self.csource.set_emission_point(i, value)

    def get_emission_location(self, i: int) -> np.ndarray:
        cython.declare(loc=cye.Vec3)
        loc = self.csource.get_emission_location(i)
        return np.array([loc.x, loc.y, loc.z])

    def interp_emission(self, x: float, y: float, z: float) -> float:
        cython.declare(loc=cye.Vec3)
        loc = cye.Vec3(x, y, z)
        return self.csource.interp_emission(loc)

    def plot_brightness(
            self, shp: Shape, res=400, axis=None, ring_front: cython.bint = True, **args):
        from matplotlib import pyplot as plt
        from matplotlib.colors import ListedColormap

        xlim = shp.x_bounds()
        ylim = shp.y_bounds()
        x = np.linspace(xlim[0], xlim[1], res)
        y = np.linspace(ylim[0], ylim[1], res)

        if self.source_type_code != 0:
            brightness = self.get_brightness(x, y, shp, ring_front)
            brightness = np.ma.masked_equal(brightness, 0.)
        else:
            brightness = shp.line_intersects(x, y, ring_front)
            brightness = np.ma.masked_equal(brightness, 0.)
            brightness *= 0.

        # Setting color works badly with masked array inputs, swap for a cmap
        if 'color' in args.keys():
            args['cmap'] = ListedColormap([args['color']] * 2)
            del args['color']
        # Other plotting defaults
        args.setdefault('cmap', 'autumn')
        args.setdefault('shading', 'gouraud')
        args.setdefault('snap', True)
        if axis is None:
            axis = plt.gca()
        axis.pcolormesh(x, y, brightness.T, **args)
