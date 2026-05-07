
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

    def get_brightness(self, double x, double y):
        return self.source.get_brightness(x, y)
