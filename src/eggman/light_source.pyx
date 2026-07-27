
cdef class LightSource:
    _source_type_names_ = ['none', 'lambertian', 'quadratic_limb', 'nonlinear_limb', 'day_night']
    _source_n_params_ = [0, 1, 3, 5, 3]
    _source_enum_ = [SourceType.None, SourceType.Lambertian, SourceType.QuadraticLimb, SourceType.NonLinearLimb, SourceType.DayNight]

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

    def __init__(self, str source_type, list[double] source_params):
        source_type = source_type.lower().strip()
        assert source_type in self._source_type_names_, f"Error: source_type {source_type} not recognized, must be one of {self._source_type_names_}"
        cdef SourceType source_code = self._source_enum_[self._source_type_names_.index(source_type)]
        if source_type in ['quadratic_limb', 'nonlinear_limb'] and len(source_params) == self._source_n_params_[source_code] - 1:
            source_params.insert(0, 1/np.pi)
        assert len(source_params) == self._source_n_params_[source_code], f"Error: wrong number of source parameters; should be {self._source_n_params_[source_code]}, was {len(source_params)}."
        cdef double[MAX_SOURCE_PARAMS] source_params_c = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]
        for i, p in enumerate(source_params):
            source_params_c[i] = p
        self.csource = CLightSource(source_code, source_params_c)

    def get_brightness(self, x, y, Shape bell):
        types = [float, np.float64, np.float32, np.int32, np.int64, int]
        if type(x) in types and type(y) in types:
            return self.csource.get_brightness_sphere(x, y)
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        assert x.ndim == 1 and y.ndim == 1
        xlen = len(x)
        ylen = len(y)

        output = np.zeros((xlen, ylen), dtype=np.double)
        cdef double[:] x_view = x
        cdef double[:] y_view = y
        cdef double[:, :] output_view = output

        cdef int i, j
        for i in range(len(x)):
            for j in range(len(y)):
                output_view[i, j] = self.csource.get_brightness_sphere(x_view[i], y_view[j])
        return output

    def get_brightness_sphere(self, x, y):
        types = [float, np.float64, np.float32, np.int32, np.int64, int]
        if type(x) in types and type(y) in types:
            return self.csource.get_brightness_sphere(x, y)
        x = np.atleast_1d(x)
        y = np.atleast_1d(y)
        assert x.ndim == 1 and y.ndim == 1
        xlen = len(x)
        ylen = len(y)

        output = np.zeros((xlen, ylen), dtype=np.double)
        cdef double[:] x_view = x
        cdef double[:] y_view = y
        cdef double[:, :] output_view = output

        cdef int i, j
        for i in range(len(x)):
            for j in range(len(y)):
                output_view[i, j] = self.csource.get_brightness_sphere(x_view[i], y_view[j])
        return output

    def plot_brightness(self, Shape bell, res=400, pcm_args=dict()):
        from matplotlib import pyplot as plt

        xlim = bell.x_bounds()
        ylim = bell.y_bounds()
        x = np.linspace(*xlim, res)
        y = np.linspace(*ylim, res)

        brightness = self.get_brightness(x, y, bell)
        brightness = np.ma.masked_equal(brightness, 0.)

        pcm_args.setdefault('cmap', 'autumn')
        plt.pcolormesh(x, y, brightness.T, **pcm_args)
