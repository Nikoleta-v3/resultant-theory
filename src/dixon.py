"""
A class for creating the Dixon's matrix and calculating the Dixon's
resultant. Literature:
"""

import sympy as sym

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
        n: integer
            Number of variables
        m: integer
            Number of $n$-degree polynomials $(n + 1)$.
        variables: list
            A list of all $n$ variables
        polynomials : list of sympy polynomials
            A  list of $m$ $n$-degree polynomials
        degree: integer
        """

        self.polynomials = polynomials
        self.variables = variables

        self.n = len(self.variables)
        self.m = len(self.polynomials)

        a = sym.IndexedBase("alpha")
        self.dummy_variables = [a[i] for i in range(self.n)]
        
        self.max_degree = self.degree()

    def degree(self):
        degree = []
        for v in self.variables:
            degree.append(max([sym.Poly(f(*self.variables)).degree(v) for f in self.polynomials]))
        return sum(degree)

    def get_dixon_polynomial(self):
        
        if self.m != (self.n + 1):
            raise Exception('Method invalid for given combination.')
        
        # first row
        rows = [[poly(*self.variables) for poly in self.polynomials]]
        
        temp = [*self.variables]
        iterator = iter(self.dummy_variables)
        
        for idx in (idx for idx in range(self.n)):
            temp[idx] = next(iterator)
            rows.append([poly(*temp) for poly in self.polynomials])
 
        matrix = sym.Matrix(rows)
        product_of_differences = functools.reduce(lambda x, y: x * y, 
                                                  [a - b for a, b in zip(self.variables, self.dummy_variables)])
        return sym.Poly((matrix.det() / product_of_differences).factor(), *self.dummy_variables)

    
    def get_coeff_of_alpha(self, polynomial):
        """
        Returns the coefficients of terms x_i...x_n in \delta when viewed as a polynomial
        of \alpha_i...\alpha_n.
        """
        coeffs = []
        for monoms_powers in polynomial.monoms():
            mono = functools.reduce(lambda i, j: i * j,
                                    [a ** b for a, b in zip(self.dummy_variables, monoms_powers)])
            coeffs.append(polynomial.coeff_monomial(mono))

        return coeffs
    
    def dixon_matrix(self, coeffs, size):
        """
        Construct the Dixon matrix from the coefficients of polynomial \alpha. Each coefficient is
        viewed as a polynomial of x and y.
        """  
        mononomials = list(itermonomials(variables, self.max_degree))
        array = np.array([[sym.Poly(c, *variables).coeff_monomial(m) for m in mononomials] for c in coeffs])
    
        return sym.Matrix(np.delete(array, np.nonzero((array==0).sum(axis=0) == size), axis=1))