
cdef class LightSourceWrap:
    @property
    def limb_type(self):
        return ['quadratic', 'nonlinear'][self.source.source_type]

    @property
    def limb_code(self):
        return self.source.source_type

    @property
    def limb_params(self):
        return np.array([self.source.limb])

    @property
    def normalization(self):
        return self.source.normalization

    def __init__(self, str limb_type, list[double] limb_params):
        cdef int limb_code = -1
        limb_type = limb_type.lower().strip()
        if limb_type == "quadratic":
            assert len(limb_params) == 2, "Error: Quadratic limb darkening takes exactly two parameters."
            limb_code = 0
            limb_params = [limb_params[0], limb_params[1], 0., 0.]
        elif limb_type == "nonlinear":
            assert len(limb_params) == 4, "Error: Nonlinear limb darkening takes exactly two parameters."
            limb_code = 1
        else:
            raise ValueError("limb_type not recognized, must be one of quadratic, nonlinear.")
        self.source = LightSource(limb_code, limb_params[0], limb_params[1], limb_params[2], limb_params[3])

    def get_brightness(self, x, y):
        types = [float, np.float64, np.float32, np.int32, np.int64, int]
        if type(x) in types and type(y) in types:
            return self.source.get_brightness(x, y)
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
                output_view[i, j] = self.source.get_brightness(x_view[i], y_view[j])
        return output

    def plot_brightness(source, res=400, pcm_args=dict()):
        from matplotlib import pyplot as plt

        # TODO: Allow biellipsoid input for positioning
        x = np.linspace(-1, 1, res)
        y = np.linspace(-1, 1, res)

        brightness = source.get_brightness(x, y)
        brightness = np.ma.masked_equal(brightness, 0.)

        pcm_args.setdefault('cmap', 'autumn')
        plt.pcolormesh(x, y, brightness, **pcm_args)
