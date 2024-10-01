import unittest
import logging
from retmap import density


class TestDensityCase(unittest.TestCase):
    def test_ozawa(self):
        result = density.ozawa(1, 1, 1, 1)
        self.assertTrue(isinstance(result, (float, int)))

    def test_hauer(self):
        result = density.hauer(1, 1, 1, 1)
        self.assertTrue(isinstance(result, (float, int)))

    def test_empirical(self):
        result = density.empirical(1, 1, 1)
        self.assertTrue(isinstance(result, (float, int)))

    def test_extrapolation(self):
        result = density.extrapolation(1000, "local", "Ar")
        self.assertTrue(isinstance(result, (float, int)))

        result = density.extrapolation(100, "local", "Ar")
        self.assertTrue(isinstance(result, (float, int)))


if __name__ == '__main__':
    unittest.main()
