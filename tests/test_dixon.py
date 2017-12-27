"""
A file to test Dixon's matrix.
"""
import sys
sys.path.insert(0, '../src/')

import unittest
import sympy as sym

from dixon import DixonResultant

c, d = sym.symbols("a, b")
x, y = sym.symbols("x, y")

p = sym.lambdify((x, y), c * x + y)
q = sym.lambdify((x, y), x + d * y)
polynomials = [p, q]
variables = [x, y]

dixon = DixonResultant(polynomials=polynomials, variables=variables)

class TestDixonResultant(unittest.TestCase):

    def test_dixon_resultant_init(self):
        """Test init method of DixonResultant."""
        a = sym.IndexedBase("alpha")

        self.assertEqual(dixon.polynomials, polynomials)
        self.assertEqual(dixon.variables, variables)
        self.assertEqual(dixon.n, 2)
        self.assertEqual(dixon.m, 2)
        self.assertEqual(dixon.dummy_variables, [a[0], a[1]])
        self.assertEqual(dixon.max_degrees, [1, 1])

    def test_get_dummy_variables(self):

        dummy_variables = dixon.get_dummy_variables()

        self.assertIsInstance(dummy_variables, list)
        self.assertEqual(len(dummy_variables), 2)

    def test_get_max_degrees(self):

        max_degrees = dixon.get_max_degrees()

        self.assertIsInstance(max_degrees, list)
        self.assertEqual(len(max_degrees), 2)

    def test_get_dixon_polynomial_numerical(self):
        """Test Dixon's polynomial for a numerical example."""

        x, y = sym.symbols('x, y')
        a = sym.IndexedBase("alpha")

        p = sym.lambdify((x, y), x + y)
        q = sym.lambdify((x, y), x ** 2 + y **3)
        h = sym.lambdify((x, y), x ** 2 + y)

        dixon = DixonResultant([p, q, h], [x, y])
        polynomial = -x * y ** 2 * a[0] - x * y ** 2 * a[1] - x * y * a[0] * a[1] \
        - x * y * a[1] ** 2 - x * a[0] * a[1] ** 2 + x * a[0] - y ** 2 * a[0] * a[1] \
        + y ** 2 * a[1] - y * a[0] * a[1] ** 2 + y * a[1] ** 2

        self.assertEqual(dixon.get_dixon_polynomial().factor(), polynomial)

    def test_get_coefficients_of_alpha_numerical(self):
        """Test Dixon's coefficients of a_1,...,a_n products for a numerical example."""

        x, y = sym.symbols('x, y')
        a = sym.IndexedBase("alpha")

        p = sym.lambdify((x, y), x + y)
        q = sym.lambdify((x, y), x ** 2 + y **3)
        h = sym.lambdify((x, y), x ** 2 + y)

        dixon = DixonResultant([p, q, h], [x, y])
        coefficients = [-x - y, -x * y - y ** 2, -x * y ** 2 + x, -x  * y + y,
                        -x * y ** 2 + y ** 2]
        polynomial = dixon.get_dixon_polynomial()

        self.assertEqual(dixon.get_coefficients_of_alpha(polynomial), coefficients)

    def test_get_dixon_matrix(self):
        """Test Dixon's coefficients of a_1,...,a_n products for a numerical example."""

        x, y = sym.symbols('x, y')
        a = sym.IndexedBase("alpha")

        p = sym.lambdify((x, y), x + y)
        q = sym.lambdify((x, y), x ** 2 + y **3)
        h = sym.lambdify((x, y), x ** 2 + y)

        dixon = DixonResultant([p, q, h], [x, y])
        polynomial = dixon.get_dixon_polynomial()

        self.assertEqual(dixon.get_dixon_matrix(polynomial).det(), 0)
