"""
A class for creating the Macaulay's matrix and calculating the Macaulay's
resultant. Literature:
"""

import sympy as sym

from sympy.polys.monomials import itermonomials
from sympy.functions.combinatorial.factorials import binomial

class Macaulay():
    """
    A class for calculating the Macaulay resultant.
    """
    def __init__(self, polynomials, variables):
        """
        Parameters
        ----------          
        n: integer
            Number of variables and polynomials.
        variables: list
            A list of all n variables
        polynomials : list of sympy polynomials
            A list of n ndegree popynomials
        degrees: list
            A list of max degree of each polynomial
        degree_m: integer
            The degree_m (as referenced in literature)
        mononomials_size: integer
            The size of the set containg all the possible monomials of our variables of degree d_m
        """
        self.polynomials = polynomials
        self.variables = variables
        self.n = len(variables)
        self.degrees = [self.get_polynomial_degree(poly) for poly in self.polynomials]
        self.degree_m = self.get_degree_m()
        self.mononomials_size = self.get_size()


    def get_polynomial_degree(self, poly):
        return sym.Poly(poly(*self.variables)).degree()

    def get_degree_m(self):
        return 1 + sum([d - 1 for d in self.degrees])

    def get_size(self):
        return binomial(self.degree_m + self.n - 1, self.n - 1)

    def get_monomials_of_certain_degree(self, degree):
        return list(itermonomials(self.variables, degree) -
                    itermonomials(self.variables, degree - 1))

    def get_monomials_set(self):
        monomials = self.get_monomials_of_certain_degree(self.degree_m)
        self.monomials = monomials

    def get_row_coefficients(self):
        row_coeff = []
        divisable = []
        for i in range(self.n):
            if i == 0:
                row_coeff.append(self.get_monomials_of_certain_degree(self.degree_m - self.degrees[i]))

            else:
                divisable.append(self.variables[i - 1] ** self.degrees[i - 1])
                poss_rows = self.get_monomials_of_certain_degree(self.degree_m - self.degrees[i])
                for div in divisable:
                    for p in poss_rows:
                        if p % div == 0:
                            poss_rows.remove(p)
                row_coeff.append(poss_rows)   
        return row_coeff

    def get_matrix(self):
        rows = []
        row_coeff = self.get_row_coefficients()
        for i in range(self.n):
            for multiplier in row_coeff[i]:
                coeffs = []
                poly = sym.Poly(self.polynomials[i](*self.variables) * multiplier, *self.variables)

                for mono in self.monomial_set:
                    coeffs.append(poly.coeff_monomial(mono))
                rows.append(coeffs)
        return sym.Matrix(rows)

    def get_reduced_nonreduced(self):

        divisible = []
        for m in self.monomial_set:
            temp = []
            for i, v in enumerate(self.variables):
                temp.append(sym.Poly(m, v).degree() >= self.degrees[i])
            divisible.append(temp)

        reduced = [i for i, r in enumerate(divisible) if sum(r) < 2]
        non_reduced = [i for i, r in enumerate(divisible) if sum(r) >= 2]

        return reduced, non_reduced

    def get_submatrix(self, matrix):
        reduced, non_reduced = self.get_reduced_nonreduced()
        
        ais = list([self.polynomials[i](*self.variables).coeff(self.variables[i] ** self.degrees[i]) 
                    for i in range(self.n)])

        reduced_matrix = matrix[:, reduced]
        keep = []
        for row in range(reduced_matrix.rows):
            check = [ai in reduced_matrix[row, :] for ai in ais]
            if True not in check:
                keep.append(row)

        return matrix[keep, non_reduced]