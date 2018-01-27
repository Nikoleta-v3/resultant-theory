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

    def test_macaulay_example(self):
        """Tests the Macaulay formulation. This example is from:https://dl.acm.org/citation.cfm?id=550525"""
        
        x, y, z = sym.symbols('x, y, z')
        a_1_1, a_1_2, a_1_3, a_2_2, a_2_3, a_3_3 = sym.symbols('a_1_1, a_1_2, a_1_3, a_2_2, a_2_3, a_3_3')
        b_1_1, b_1_2, b_1_3, b_2_2, b_2_3, b_3_3 = sym.symbols('b_1_1, b_1_2, b_1_3, b_2_2, b_2_3, b_3_3')
        c_1, c_2, c_3 = sym.symbols('c_1, c_2, c_3')

        f_1 = sym.lambdify((x, y, z), a_1_1 * x ** 2 + a_1_2 * x * y + a_1_3 * x * z \
                   + a_2_2 * y ** 2 + a_2_3 * y * z + a_3_3 * z ** 2)
        f_2 = sym.lambdify((x, y, z), b_1_1 * x ** 2 + b_1_2 * x * y + b_1_3 * x * z \
                   + b_2_2 * y ** 2 + b_2_3 * y * z + b_3_3 * z ** 2)
        f_3 = sym.lambdify((x, y, z), c_1 * x + c_2 * y + c_3 * z)

        mac = MacaulayResultant([f_1, f_2, f_3], [x, y, z])

        self.assertEqual(mac.degrees, [2, 2, 1])
        self.assertEqual(mac.degree_m, 3)

        mac.get_monomials_set()
        self.assertEqual(mac.monomial_set, [x ** 3,x ** 2 * y, x ** 2 * z,
                                            x * y ** 2, x * y * z, x * z ** 2,
                                            y ** 3, y ** 2 *z, y * z ** 2,
                                            z ** 3])
        self.assertEqual(mac.monomials_size, 10)
        self.assertEqual(mac.get_row_coefficients(), [[x, y, z], [x, y, z],
                                                      [x * y, x * z, y * z, z ** 2]])

        matrix = mac.get_matrix()
        self.assertEqual(matrix.shape, (mac.monomials_size, mac.monomials_size))
        self.assertEqual(mac.get_submatrix(matrix), sym.Matrix([[a_1_1, a_2_2],
                                                                [b_1_1, b_2_2]]))

    def test_macaulay_example(self):
        """Tests the Macaulay formulation. This example is from:http://isc.tamu.edu/resources/preprints/1996/1996-02.pdf"""
        
        x, y, z = sym.symbols('x, y, z')
        a_0, a_1, a_2 = sym.symbols('a_0, a_1, a_2')
        b_0, b_1, b_2 = sym.symbols('b_0, b_1, b_2')
        c_0, c_1, c_2,c_3, c_4 = sym.symbols('c_0, c_1, c_2, c_3, c_4')

        f = sym.lambdify((x, y, z), a_0 * y -  a_1 * x + a_2 * z)
        g = sym.lambdify((x, y, z), b_1 * x ** 2 + b_0 * y ** 2 - b_2 * z ** 2)
        h = sym.lambdify((x, y, z), c_0 * y - c_1 * x ** 3 + c_2 * x ** 2 * z - c_3 * x * z ** 2 + c_4 * z ** 3)

        mac = MacaulayResultant([f, g, h], [x, y, z])

        self.assertEqual(mac.degrees, [1, 2, 3])
        self.assertEqual(mac.degree_m, 4)

        mac.get_monomials_set()

        self.assertEqual(mac.monomials_size, 15)
        self.assertEqual(len(mac.get_row_coefficients()), mac.n)

        matrix = mac.get_matrix()
        self.assertEqual(matrix.shape, (mac.monomials_size, mac.monomials_size))
        self.assertEqual(mac.get_submatrix(matrix), sym.Matrix([[-a_1, a_0, a_2, 0],
                                                                [0, -a_1, 0, 0],
                                                                [0, 0, -a_1, 0],
                                                                [0, 0, 0, -a_1]]))