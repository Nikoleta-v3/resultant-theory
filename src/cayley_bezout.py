"""
A function for creating the Cayley-Bezout's matrix and calculating the Bezout's
resultant. 

The Bezout formulation has gone over different generalizations. The most common
one is the Cayley that is used after for the Dixon formulation as well.

Literature: https://dl.acm.org/citation.cfm?doid=190347.190372.
"""

import sympy as sym


def cayley_bezout_matrix(p, q, x):
    """
    A function that takes two univariate non zero polynomials and returns the
    Caley-Bezout formulation matrix $n\timesn$, where $n$ is the maximum degree
    of the two polynomials.

    Parameters
    ----------
    p: sympy lambdify function
        A non zero polynomial.
    q: sympy lambdify function
        A non zero polynomial.
    x: sympy symbol
        Variable for which we are solving.
    """
    a = sym.symbols('a')
    degree = max(sym.Poly(p(x)).degree(), sym.Poly(q(x)).degree())

    matrix = sym.Matrix([[p(x), q(x)], [p(a), q(a)]])

    bezout_polynomial = (matrix.det() / (x - a)).factor().collect(a)
    coefficients = sym.Poly(bezout_polynomial, a).all_coeffs()

    monomials = [x ** power for power in range(degree)]

    bezout_matrix = sym.Matrix([[sym.Poly(row, x).coeff_monomial(m) for m in monomials]
                                for row in coefficients])

    return bezout_matrix

