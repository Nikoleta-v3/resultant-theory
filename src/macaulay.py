"""
A class for creating the Macaulay's matrix and calculating the Macaulay's
resultant. Literature:
"""

import sympy as sym

from sympy.polys.monomials import itermonomials
from sympy.polys.orderings import monomial_key
from sympy.functions.combinatorial.factorials import binomial

class MacaulayResultant():
    """
    A class for calculating the Macaulay resultant.
    """
    def __init__(self, polynomials, variables):
        """
        Parameters
        ----------          
        A class that takes two lists, a list of n polynomials and list of n
        variables. Returns the Maucalay's matrices of the multivariate system.

        Parameters
        ----------
        variables: list
            A list of all n variables
        polynomials : list of sympy polynomials
            A  list of m n-degree polynomials
        """
        self.polynomials = polynomials
        self.variables = variables
        self.n = len(variables)

        self.degrees = self.get_max_degrees()
        self.degree_m = self.get_degree_m()
        self.monomials_size = self.get_size()

    def get_max_degrees(self):
        """
        Returns
        -------
        degrees: list
            A list of the d_max of each polynomial
        """
        degrees = [self.get_polynomial_degree(poly) for poly in self.polynomials]
        return degrees

    def get_polynomial_degree(self, poly):
        """
        Returns
        -------
        degree: int
            The degree of a polynomial
        """
        return max(sym.Poly(poly(*self.variables)).degree_list())

    def get_degree_m(self):
        """
        Returns
        -------
        degree_m: int
            The degree_m is calculated as  1 + \sum_1 ^ n (d_i - 1), where
            d_i is the degree of the i polynomial
        """
        return 1 + sum([d - 1 for d in self.degrees])

    def get_size(self):
        """
        Returns
        -------
        size: int
            The size of set T. Set T is the set of all possible monomials of
            the n variables for degree equal to the degree_m 
        """
        return binomial(self.degree_m + self.n - 1, self.n - 1)

    def get_monomials_of_certain_degree(self, degree):
        """
        Returns
        -------
        monomials: list
            A list of monomials of a certain degree. Sympy returns up to a degree
        """
        monomials = list(itermonomials(self.variables, degree) -
                         itermonomials(self.variables, degree - 1))

        return sorted(monomials, key=monomial_key('lex', self.variables))[::-1]

    def get_monomials_set(self):
        """
        Returns
        -------
        self.monomial_set: set
            The set T. Set of all possible monomials of degree degree_m
        """
        monomial_set = self.get_monomials_of_certain_degree(self.degree_m)
        self.monomial_set = monomial_set

    def get_row_coefficients(self):
        """
        Returns
        -------
        row_coefficients: list
            The row coefficients of Macaulay's matrix
        """
        row_coefficients = []
        divisible = []
        for i in range(self.n):
            if i == 0:
                row_coefficients.append(self.get_monomials_of_certain_degree(self.degree_m - self.degrees[i]))
            else:
                divisible.append(self.variables[i - 1] ** self.degrees[i - 1])
                degree = self.degree_m - self.degrees[i]
                if degree == 0:
                    poss_rows = [1]
                else:
                    poss_rows = self.get_monomials_of_certain_degree(self.degree_m - self.degrees[i])
                for div in divisible:
                    for p in poss_rows:
                        if sym.fraction((p / div).expand())[1] == 1 :
                            poss_rows = [item for item in poss_rows if item != p]
                row_coefficients.append(poss_rows)
        return row_coefficients

    def get_matrix(self):
        """
        Returns
        -------
        macaulay_matrix: sym Matrix
            The Macaulay's matrix
        """
        rows = []
        row_coefficients = self.get_row_coefficients()
        for i in range(self.n):
            for multiplier in row_coefficients[i]:
                coefficients = []
                poly = sym.Poly(self.polynomials[i](*self.variables) * multiplier, *self.variables)

                for mono in self.monomial_set:
                    coefficients.append(poly.coeff_monomial(mono))
                rows.append(coefficients)

        macaulay_matrix = sym.Matrix(rows)
        return macaulay_matrix

    def get_reduced_nonreduced(self):
        """
        Returns
        -------
        reduced: list
            A list of the reduced monomials
        non_reduced: list
            A list of the monomials that are not reduced
        """
        divisible = []
        for m in self.monomial_set:
            temp = []
            for i, v in enumerate(self.variables):
                temp.append(sym.Poly(m, v).degree() >= self.degrees[i])
            divisible.append(temp)

        non_reduced = []
        minus = 1
        while not non_reduced:
            reduced = [i for i, r in enumerate(divisible) if sum(r) < self.n - minus]
            non_reduced = [i for i, r in enumerate(divisible) if sum(r) >= self.n - minus]
            minus += 1

        return reduced, non_reduced 

    def get_submatrix(self, matrix):
        """
        Returns
        -------
        macaulay_submatrix: sym Matrix
            The Macaulay's matrix
        """
        reduced, non_reduced = self.get_reduced_nonreduced()

        reduction_set = [v ** self.degrees[i] for i, v in enumerate(self.variables)]

        ais = list([self.polynomials[i](*self.variables).coeff(reduction_set[i])
                    for i in range(self.n)])

        reduced_matrix = matrix[:, reduced]
        keep = []
        for row in range(reduced_matrix.rows):
            check = [ai in reduced_matrix[row, :] for ai in ais]
            if True not in check:
                keep.append(row)

        return matrix[keep, non_reduced]
