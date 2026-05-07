
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

    def __init__(lelf, str limb_type, list[double] limb_params):
        cdef int limb_code = -1
        limbType = limbType.lower().strip()
        if limbType == "quadratic":
            assert len(limb) == 2, "Error: Quadratic limb darkening takes exactly two parameters."
            limb_code = 0
            limb = [limb[0], limb[1], 0., 0.]
        elif limbType == "nonlinear":
            assert len(limb) == 4, "Error: Nonlinear limb darkening takes exactly two parameters."
            limb_code = 1
        else:
            raise ValueError("limbType not recognized, must be one of quadratic, nonlinear.")
        self.source = LightSource(limb_code, limb[0], limb[1], limb[2], limb[3])

    def get_brightness(double x, double y):
        return self.source.get_brightness(x, y)
