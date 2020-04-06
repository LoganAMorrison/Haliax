import unittest
from lanre.utils import besselk0, besselk0e
from lanre.utils import besselk1, besselk1e
from lanre.utils import besselk2, besselk2e
from lanre.utils import besselk3, besselk3e
from lanre.utils import besselkn, besselkne
from scipy.special import k0, k1, k0e, k1e, kn, kve


class TestBessel(unittest.TestCase):
    def test_besselk0(self):
        xs = [0.1, 1.0, 10.0]
        for x in xs:
            sp = k0(x)
            hal = besselk0(x)
            relerr = abs((sp - hal) / sp)
            self.assertLessEqual(relerr, 0.05)

    def test_besselk0e(self):
        xs = [0.1, 1.0, 10.0]
        for x in xs:
            sp = k0e(x)
            hal = besselk0e(x)
            relerr = abs((sp - hal) / sp)
            self.assertLessEqual(relerr, 0.05)

    def test_besselk1(self):
        xs = [0.1, 1.0, 10.0]
        for x in xs:
            sp = k1(x)
            hal = besselk1(x)
            relerr = abs((sp - hal) / sp)
            self.assertLessEqual(relerr, 0.05)

    def test_besselk1e(self):
        xs = [0.1, 1.0, 10.0]
        for x in xs:
            sp = k1e(x)
            hal = besselk1e(x)
            relerr = abs((sp - hal) / sp)
            self.assertLessEqual(relerr, 0.05)

    def test_besselk2(self):
        xs = [0.1, 1.0, 10.0]
        for x in xs:
            sp = kn(2, x)
            hal = besselk2(x)
            relerr = abs((sp - hal) / sp)
            self.assertLessEqual(relerr, 0.05)

    def test_besselk2e(self):
        xs = [0.1, 1.0, 10.0]
        for x in xs:
            sp = kve(2, x)
            hal = besselk2e(x)
            relerr = abs((sp - hal) / sp)
            self.assertLessEqual(relerr, 0.05)

    def test_besselk3(self):
        xs = [0.1, 1.0, 10.0]
        for x in xs:
            sp = kn(3, x)
            hal = besselk3(x)
            relerr = abs((sp - hal) / sp)
            self.assertLessEqual(relerr, 0.05)

    def test_besselk3e(self):
        xs = [0.1, 1.0, 10.0]
        for x in xs:
            sp = kve(3, x)
            hal = besselk3e(x)
            relerr = abs((sp - hal) / sp)
            self.assertLessEqual(relerr, 0.05)

    def test_besselk10(self):
        xs = [0.1, 1.0, 10.0]
        for x in xs:
            sp = kn(10, x)
            hal = besselkn(10, x)
            relerr = abs((sp - hal) / sp)
            self.assertLessEqual(relerr, 0.05)

    def test_besselk10e(self):
        xs = [0.1, 1.0, 10.0]
        for x in xs:
            sp = kve(10, x)
            hal = besselkne(10, x)
            relerr = abs((sp - hal) / sp)
            self.assertLessEqual(relerr, 0.05)


if __name__ == '__main__':
    unittest.main()
