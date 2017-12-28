"""
A file to test Macaulay's formulation.
"""
import sys
sys.path.insert(0, '../src/')

import unittest
import sympy as sym

from macaulay import MacaulayResultant

c, d = sym.symbols("a, b")
x, y = sym.symbols("x, y")

p = sym.lambdify((x, y), c * x + y)
q = sym.lambdify((x, y), x + d * y)
polynomials = [p, q]
variables = [x, y]

macaulay = MacaulayResultant(polynomials=polynomials, variables=variables)

class TestMacaulayResultant(unittest.TestCase):

    def test_dixon_resultant_init(self):
        """Test init method of MacaulayResultant."""
        a = sym.IndexedBase("alpha")

        self.assertEqual(macaulay.polynomials, polynomials)
        self.assertEqual(macaulay.variables, variables)
        self.assertEqual(macaulay.n, 2)
        self.assertEqual(macaulay.degrees, [1, 1])
        self.assertEqual(macaulay.degree_m, 1)
        self.assertEqual(macaulay.monomials_size, 2)

    def test_get_max_degrees(self):
        max_degrees = macaulay.get_max_degrees()

        self.assertIsInstance(max_degrees, list)
        self.assertEqual(len(max_degrees), 2)

    def test_get_polynomial_degree(self):
        self.assertEqual(macaulay.get_polynomial_degree(polynomials[0]), 1)
        self.assertEqual(macaulay.get_polynomial_degree(polynomials[1]), 1)

    def test_get_degree_m(self):
        self.assertEqual(macaulay.get_degree_m(), 1)

    def test_get_size(self):
        self.assertEqual(macaulay.get_size(), 2)