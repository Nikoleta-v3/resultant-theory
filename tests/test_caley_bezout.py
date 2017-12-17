"""
A file to test Caley-Bezout's matrix.
"""
import sys
sys.path.insert(0, '../src/')

import unittest
import sympy as sym

from caley_bezout import caley_bezout_matrix


class TestCaleyBezout(unittest.TestCase):

    def test_caley_bezout_matrix_generic_case(self):
        """Test Caley-Bezout's matrix for a generic case."""
        b_2, b_1, b_0 = sym.symbols("b_2, b_1, b_0")
        x = sym.symbols("x")

        p = sym.lambdify(x, b_2 * x ** 2 + b_1 * x + b_0)
        q = sym.lambdify(x, sym.diff(p(x), x))

        matrix = caley_bezout_matrix(p, q, x)
        self.assertEqual(matrix, sym.Matrix([[2 * b_2 ** 2, b_1 * b_2], 
                                             [b_1 * b_2, -2 * b_0 * b_2 + b_1 ** 2]]))

        self.assertEqual(matrix.det().factor(), -b_2 ** 2 * (4 * b_0 * b_2 - b_1 ** 2))
    
    def test_caley_bezout_matrix_example(self):
        """Test Caley-Bezout's matrix for a numerical example."""
        x = sym.symbols("x")

        p = sym.lambdify(x, x ** 2 - 5 * x + 6)
        q = sym.lambdify(x, x ** 2 - 3 * x + 2)

        matrix = caley_bezout_matrix(p, q, x)
        self.assertEqual(matrix, sym.Matrix([[2, -4], 
                                             [-4, 8]]))

        self.assertEqual(matrix.det().factor(), 0)