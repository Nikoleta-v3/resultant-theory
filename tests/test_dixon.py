"""
A file to test Dixon's matrix.
"""
import sys
sys.path.insert(0, '../src/')

import unittest
import sympy as sym

from dixon import DixonResultant

class TestDixonResultant(unittest.TestCase):

    def test_dixon_resultant_init(self):
        """Test init method of DixonResultant."""
        c, d = sym.symbols("a, b")
        x, y = sym.symbols("x, y")

        p = sym.lambdify((x, y), c * x + y)
        q = sym.lambdify((x, y), x + d * y)

        polynomials = [p, q]
        variables = [x, y]

        dixon = DixonResultant(polynomials=polynomials, variables=variables)

        self.assertEqual(dixon.polynomials, polynomials)
        self.assertEqual(dixon.variables, variables)
        self.assertEqual(dixon.n, 2)
        self.assertEqual(dixon.m, 2)
        self.assertEqual(len(dixon.dummy_variables), 2)
        self.assertEqual(dixon.max_degree, 2)