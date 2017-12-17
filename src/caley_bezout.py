"""
A function for creating the Caley-Bezout's matrix and calculating the Bezout's
resultant. 

The Bezout formulation has gone over different generalizations. The most common
one is the Cayley that is used after for the Dixon formulation as well.

Literature: https://dl.acm.org/citation.cfm?id=550525.
"""

import sympy as sym


def caley_bezout_matrix(p, q, x):
    """
    A function that takes two non zero polynomials and returns the Caley-Bezout
    formulation matrix $n\timesn$, where $n$ is the maximum degree of the two
    polynomials.

    Parameters
    ----------
    p: sympy lambdify function
        A non zero polynomial.
    q: sympy lambdify function
        A non zero polynomial.
    x: sympy symbol
        Variable for which we are solving.
    """
    a = sym.symbols('alpha')
    bezout_matrix = sym.Matrix([[p(x), q(x)], [p(a), q(a)]])

    bezout_polynomial = (bezout_matrix.det() / (x - a)).factor().collect(a)
    coefficients = sym.Poly(bezout_polynomial, a).coeffs()

    matrix = sym.Matrix([sym.Poly(row, x).coeffs() for row in coefficients])

    return matrix