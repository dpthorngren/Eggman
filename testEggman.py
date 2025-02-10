import unittest
import numpy as np
import eggman

class EggmanTester(unittest.TestCase):

    def testGrazing(self):
        # Grazing transits can cause numerical issues, let's double-check that doesn't happen.
        # Grazing on ingress
        baseArgs = {'t':np.array([-np.arcsin(1.1/10)]), 't0': 0., 'period':2*np.pi, 'semimajor':10.,
                    'inclination':90., 'limbType':'quadratic', 'limb':[.3, .2]}
        try:
            eggman.asymmetricTransit(.1+1e-3, .11, .1, **dict(baseArgs))
            eggman.asymmetricTransit(.1+1e-6, .09, -1, **dict(baseArgs))
            eggman.asymmetricTransit(.1+1e-12, .1, .09, **dict(baseArgs))
        except:
            self.fail("Error raised on ingress grazing transit.")
        # Just barely non-grazing
        try:
            result1 = eggman.asymmetricTransit(.1-1e-3, .11, .1, **dict(baseArgs))
            result2 = eggman.asymmetricTransit(.1-1e-6, .09, -1, **dict(baseArgs))
            result3 = eggman.asymmetricTransit(.1-1e-12, .1, .09, **dict(baseArgs))
        except:
            self.fail("Error raised on ingress barely non-grazing transit.")
        self.assertEqual(result1, 1.)
        self.assertEqual(result2, 1.)
        self.assertEqual(result3, 1.)

        # Grazing on egress
        baseArgs['t'] = -baseArgs['t']
        try:
            eggman.asymmetricTransit(.11, .1+1e-3, .1, **dict(baseArgs))
            eggman.asymmetricTransit(.09, .1+1e-6, -1, **dict(baseArgs))
            eggman.asymmetricTransit(.1, .1+1e-12, .09, **dict(baseArgs))
        except:
            self.fail("Error raised on egress grazing transit.")
        eggman.asymmetricTransit(.1+1e-6, .1, .1, **dict(baseArgs, t=np.array([5.])))

        # Grazing on top (true grazing transit)
        baseArgs['inclination'] = 90. - abs(baseArgs['t'][0])*180/np.pi
        baseArgs['t'] = np.array([0.])
        try:
            eggman.asymmetricTransit(.11, .1, .1+1e-3, **dict(baseArgs))
            eggman.asymmetricTransit(.09, .1, .1+1e-6, **dict(baseArgs))
            eggman.asymmetricTransit(.1, .09, .1+1e-12, **dict(baseArgs))
        except:
            self.fail("Error raised on top-side (t=mid-transit, i<90) grazing transit.")

        # Barely non-grazing on top
        try:
            result1 = eggman.asymmetricTransit(.11, .1, .1-1e-3, **dict(baseArgs))
            result2 = eggman.asymmetricTransit(.09, .1, .1-1e-6, **dict(baseArgs))
            result3 = eggman.asymmetricTransit(.1, .09, .1-1e-12, **dict(baseArgs))
        except:
            self.fail("Error raised on top-side (t=mid-transit, i<90) barely non-transit.")
        self.assertEqual(result1, 1.)
        self.assertEqual(result2, 1.)
        self.assertEqual(result3, 1.)


    def testIsTransiting(self):
        # Ensure we correctly can tell transit from non-transit
        baseArgs = {'t':np.array([0.]), 't0': 0., 'period':10., 'semimajor':10.,
                    'inclination':89., 'limbType':'quadratic', 'limb':[.3, .2]}
        result = eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=np.array([5.])))
        self.assertEqual(result[0], 1.)
        result = eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=np.array([3.])))
        self.assertEqual(result[0], 1.)
        result = eggman.asymmetricTransit(1.1, 1.1, 1.1, **dict(baseArgs, t=np.array([2.49]), semimajor=1.))
        self.assertLess(result[0], 1.)
        result = eggman.asymmetricTransit(1.1, 1.1, 1.1, **dict(baseArgs, t=np.array([2.51]), semimajor=1.))
        self.assertEqual(result[0], 1.)
        result = eggman.asymmetricTransit(.01, .015, .011, **dict(baseArgs, t=np.array([.01])))
        self.assertLess(result[0], 1.)
        result = eggman.asymmetricTransit(.01, .015, .011, **dict(baseArgs, t=np.array([-10.01])))
        self.assertLess(result[0], 1.)
        result = eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=np.array([-25.01])))
        self.assertEqual(result[0], 1.)
        result = eggman.asymmetricTransit(.1, .1, .1, **dict(baseArgs, t=np.linspace(-10,10,1000), inclination=80.))
        self.assertTrue(all(result == 1.))

    def testAsymmetricAnalytic(self):
        # Test cases where the transit depth is known analytically
        baseArgs = {'t':np.array([0.]), 't0': 0., 'period':10., 'semimajor':10.,
                    'inclination':89., 'limbType':'quadratic', 'limb':[0.,0.]}
        result = eggman.asymmetricTransit(.1, .1, .1, **baseArgs)
        self.assertAlmostEqual(result, 1-.1*.1, delta=1e-6)
        result = eggman.asymmetricTransit(.1, .2, -1, **baseArgs)
        self.assertAlmostEqual(result, 1-(.1*.1+.2*.2)/2., delta=1e-6)
        result = eggman.asymmetricTransit(.12, .19, .15, **baseArgs)
        self.assertAlmostEqual(result, 1-(.12*.15+.19*.15)/2., delta=1e-6)
        result = eggman.asymmetricTransit(.01, .01, .02, **baseArgs)
        self.assertAlmostEqual(result, 1-(.01*.02+.01*.02)/2., delta=1e-6)

    def testAsymmetricNumerical(self):
        baseArgs = {'t':np.array([0.001, .01]), 't0': 0., 'period':1., 'semimajor':15.,
                    'inclination':90., 'limbType':'quadratic', 'limb':[.1, .3]}
        # Test cases where the transit depth was calculated from other codes.
        result = eggman.asymmetricTransit(.11, .1, -1, **baseArgs)
        self.assertAlmostEqual(result[0], 0.9879552837718931)
        self.assertAlmostEqual(result[1], 0.9923392211352866)

if __name__ == "__main__":
    unittest.main()
