"""
A file to test Sylvester's matrix.
"""
import sys
sys.path.insert(0, '../src/')

import unittest
import sympy as sym

from sylvesters import sylvester_matrix


class TestSylvester(unittest.TestCase):

    def test_sylvester_matrix_generic_case(self):
        """Test Sylvester's matrix for a generic case."""
        a = sym.IndexedBase("a")
        b = sym.IndexedBase("b")
        x = sym.symbols("x")

        p = a[1] * x + a[0]
        q = b[2] * x ** 2 + b[1] * x + b[0]

        matrix = sylvester_matrix(p, q, x)
        self.assertEqual(matrix, sym.Matrix([[a[1], a[0], 0],
                                             [0, a[1], a[0]],
                                             [b[2], b[1], b[0]]]))
        self.assertEqual(matrix.det(),
                         a[0] ** 2 * b[2] - a[0] * a[1] * b[1] + a[1] ** 2 * b[0])

    def test_sylvester_matrix_example(self):
        """Test Sylvester's matrix for a numerical example."""
        x = sym.symbols("x")

        p = x ** 2 - 5 * x + 6
        q = x ** 2 - 3 * x + 2

        matrix = sylvester_matrix(p, q, x)
        self.assertEqual(matrix, sym.Matrix([[1, -5, 6, 0],
                                             [0, 1, -5, 6],
                                             [1, -3, 2, 0],
                                             [0, 1, -3, 2]]))
        self.assertEqual(matrix.det(), 0)

    def test_sylvester_matrix_example_two(self):
        """Test Sylvester's matrix for a numerical example."""
        x = sym.symbols("x")

        p = x ** 3 + 1
        q = x + 1

        matrix = sylvester_matrix(p, q, x)
        self.assertEqual(matrix.det(), 0)
