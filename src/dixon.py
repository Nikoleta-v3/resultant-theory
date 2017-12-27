"""
A class for creating the Dixon's matrix and calculating the Dixon's
resultant. Literature: https://dl.acm.org/citation.cfm?doid=190347.190372.

Dixon's resultant is a generalisation of the Cayley-Bezout resultant.
"""

import sympy as sym
import numpy as np
import functools

from sympy.polys.monomials import itermonomials

class DixonResultant():
    """
    A class for retrieving the Dixon's resultant of a multivariate system.
    """

    def __init__(self, polynomials, variables):
        """
        A class that takes two lists, a list of polynomials and list of
        variables. Returns the Dixon matrix of the multivariate system.

        Parameters
        ----------
        variables: list
            A list of all n variables
        polynomials : list of sympy polynomials
            A  list of m n-degree polynomials
        """
        self.polynomials = polynomials
        self.variables = variables

        self.n = len(self.variables)
        self.m = len(self.polynomials)

        self.dummy_variables = self.get_dummy_variables()
        self.max_degrees = self.get_max_degrees()

    def get_dummy_variables(self):
        """
        Returns
        -------

        dummy_variables: list
            A list of n alpha variables. These are the replacing variables
        """
        a = sym.IndexedBase("alpha")
        dummy_variables = [a[i] for i in range(self.n)]

        return dummy_variables

    def get_max_degrees(self):
        """
        Returns
        -------

        degrees: list
            A list of the d_max of each variable. The max degree is the
            max(degree(p_1, x_i), ..., degree(p_m, x_i))
        """
        max_degrees = []
        for v in self.variables:
            max_degrees.append(max([sym.Poly(f(*self.variables)).degree(v) for f in self.polynomials]))
        return max_degrees

    def get_dixon_polynomial(self):
        """
        Returns
        -------

        dixon_polynomial: sympy polynomial
            A sympy polynomial which formulated by Dixon's formulation.
            Dixon's polynomial is calculated as:

            delta = Delta(A) / ((x_1 - a_1) ... (x_n - a_n)) where,

            A =  |p_1(x_1,... x_n), ..., p_n(x_1,... x_n)|
                 |p_1(a_1,... x_n), ..., p_n(a_1,... x_n)|
                 |...             , ...,              ...|
                 |p_1(a_1,... a_n), ..., p_n(a_1,... a_n)|
        """
        if self.m != (self.n + 1):
            raise Exception('Method invalid for given combination.')

        # first row
        rows = [[poly(*self.variables) for poly in self.polynomials]]

        temp = [*self.variables]
        iterator = iter(self.dummy_variables)

        for idx in (idx for idx in range(self.n)):
            temp[idx] = next(iterator)
            rows.append([poly(*temp) for poly in self.polynomials])

        A = sym.Matrix(rows)
        product_of_differences = functools.reduce(lambda x, y: x * y,
                                                  [a - b for a, b in zip(self.variables, self.dummy_variables)])
        dixon_polynomial = (A.det() / product_of_differences).factor()
        return sym.Poly(dixon_polynomial, *self.dummy_variables)

    def get_coefficients_of_alpha(self, polynomial):
        """
        Returns
        --------
        coefficients: list
            A list of coefficients (in x_i, ..., x_n terms) of the power products
            a_1, ..., _n in Dixon's polynomial
        """
        coefficients = []
        for powers in polynomial.monoms():
            monomial = functools.reduce(lambda i, j: i * j,
                                    [a ** b for a, b in zip(self.dummy_variables, powers)])
            coefficients.append(polynomial.coeff_monomial(monomial))

        return coefficients
    
    def get_dixon_matrix(self, polynomial):
        """
        Construct the Dixon matrix from the coefficients of polynomial \alpha. Each coefficient is
        viewed as a polynomial of x and y.
        """
        coefficients = self.get_coefficients_of_alpha(polynomial)
        size = len(polynomial.monoms())
        monomials = list(itermonomials(self.variables, sum(self.max_degrees)))

        array = np.array([[sym.Poly(c, *self.variables).coeff_monomial(m) for m in monomials]
                                                         for c in coefficients])

        dixon_matrix = sym.Matrix(np.delete(array, np.nonzero((array==0).sum(axis=0) == size), 
                                                                        axis=1))
        return dixon_matrix
