import unittest
import numpy as np
import eggman

class EggmanTester(unittest.TestCase):
    def testIsTransiting(self):
        # Ensure we correctly can tell transit from non-transit
        limb = np.array([.3, .2, 0., 0., 0.])
        baseArgs = {'t':0., 't0': 0., 'period':10., 'semimajor':10.,
                    'inclination':89., 'limbType':0, 'limb':limb}
        result = eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=5.))
        self.assertEqual(result, 1.)
        result = eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=3.))
        self.assertEqual(result, 1.)
        result = eggman.asymmetricTransit(1.1, 1.1, 1.1, **dict(baseArgs, t=2.49, semimajor=1.))
        self.assertLess(result, 1.)
        result = eggman.asymmetricTransit(1.1, 1.1, 1.1, **dict(baseArgs, t=2.51, semimajor=1.))
        self.assertEqual(result, 1.)
        result = eggman.asymmetricTransit(.01, .015, .011, **dict(baseArgs, t=.01))
        self.assertLess(result, 1.)
        result = eggman.asymmetricTransit(.01, .015, .011, **dict(baseArgs, t=-10.01))
        self.assertLess(result, 1.)
        result = eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=-25.01))
        self.assertEqual(result, 1.)
        result = eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=-.01, inclination=80.))
        self.assertEqual(result, 1.)

    def testAsymmetricAnalytic(self):
        # Test cases where the transit depth is known analytically
        limb = np.array([0., 0., 0., 0., 0.])
        baseArgs = {'t':0., 't0': 0., 'period':10., 'semimajor':10.,
                    'inclination':89., 'limbType':0, 'limb':limb}
        result = eggman.asymmetricTransit(.1, .1, .1, **baseArgs)
        self.assertAlmostEqual(result, 1-.1*.1, delta=1e-6)
        result = eggman.asymmetricTransit(.1, .2, -1, **baseArgs)
        self.assertAlmostEqual(result, 1-(.1*.1+.2*.2)/2., delta=1e-6)
        result = eggman.asymmetricTransit(.12, .19, .15, **baseArgs)
        self.assertAlmostEqual(result, 1-(.12*.15+.19*.15)/2., delta=1e-6)
        result = eggman.asymmetricTransit(.01, .01, .02, **baseArgs)
        self.assertAlmostEqual(result, 1-(.01*.02+.01*.02)/2., delta=1e-6)

    def testAsymmetricNumerical(self):
        limb = np.array([.1, .3, 0., 0., 0.])
        baseArgs = {'t':0.001, 't0': 0., 'period':1., 'semimajor':15.,
                    'inclination':90., 'limbType':0, 'limb':limb}
        # Test cases where the transit depth was calculated from other codes.
        result = eggman.asymmetricTransit(.11, .1, -1, **baseArgs)
        self.assertAlmostEqual(result, 0.9879552837718931)
        result = eggman.asymmetricTransit(.11, .1, -1, **dict(baseArgs, t=0.01))
        self.assertAlmostEqual(result, 0.9923392211352866)

if __name__ == "__main__":
    unittest.main()
