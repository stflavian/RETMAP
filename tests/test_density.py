import unittest
import logging
from adsorpyon import density


class TestDensityCase(unittest.TestCase):
    def test_ozawa(self):
        self.assertEqual(
            first=density.ozawa(10, 10, 1, 0.5),
            second=1,
            msg="Incorrect density at boiling temperature from Ozawa model")

        self.assertEqual(
            first=density.ozawa(10, 10, 0, 0.5),
            second=0,
            msg="Density from Ozawa is not 0 when boiling density is 0")

    def test_hauer(self):
        self.assertEqual(
            first=density.hauer(10, 10, 1, 0.5),
            second=1,
            msg="Incorrect density at boiling temperature from Hauer model")

        self.assertEqual(
            first=density.hauer(10, 10, 0, 0.5),
            second=0,
            msg="Density from Hauer is not 0 when boiling density is 0")


if __name__ == '__main__':
    unittest.main()
