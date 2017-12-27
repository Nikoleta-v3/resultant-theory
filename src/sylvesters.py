"""
A function for creating the Sylvester's matrix and calculating the Sylvester's
resultant. Literature: http://comet.lehman.cuny.edu/vpan/pdf/DEKPalg.pdf
"""

import sympy as sym


def sylvester_matrix(p, q, x):
    """
    A function that takes two non zero polynomials, of degree m and n respectively
    and returns the Sylvester's matrix $(m+n)\times(m+n)$.

    Parameters
    ----------
    p: sympy expression
        A non zero polynomial of degree $m$.
    q: sympy expression
        A non zero polynomial of degree $n$.
    x: sympy symbol
        Variable for which we are solving.
    """
    p_polynomial = sym.Poly(p, x)
    q_polynomial = sym.Poly(q, x)

    m_degree, n_degree = p_polynomial.degree(), q_polynomial.degree()

    p_coefficients = p_polynomial.all_coeffs()
    q_coefficients = q_polynomial.all_coeffs()

    matrix_size = m_degree + n_degree

    matrix = []
    for row in [p_coefficients, q_coefficients]:
        matrix.append(row)
        while len(row) < matrix_size:
            row = [0] + row
            matrix.append(row)

    for row in matrix:
        while len(row) < matrix_size:
            row.append(0)

    return sym.Matrix(matrix)

