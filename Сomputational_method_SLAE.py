import numpy as np
from fractions import Fraction
import copy
from math import sqrt

def subtract(x, y, coeff):
    return list(map(lambda a,b: a - coeff * b, x, y))

#Gauss method with numpy arrays and float values, 
#gives not an exact result, but close to the style of Python

def GAUSS(matrix, b, x):
    N = len(b)
    for i, mtrxstr in enumerate(matrix):
        for j,val in enumerate(matrix[i+1:]):
            coeff = float(mtrxstr[i] / val[i])
            b[i+j+1] = b[i+j+1] * coeff - b[i]
            newstr = coeff * val - mtrxstr
            matrix[i+j+1] = np.array(newstr, dtype='f')
    for i, val in enumerate(matrix[::-1]):
        sum = 0.
        for j in range(N-i,N):
            sum += val[j] * x[j]
        x[N-i-1] = (b[N-i-1] - sum) / val[N-i-1]

    return x


def Determinant(matrix):
    N = len(matrix)
    for j in range(N):
        for i in range(j+1,N):
            coeff = matrix[i][j] / matrix[j][j]
            for k in range(j,N):
                matrix[i][k] -= coeff * matrix[j][k]
            #matrix[i] = subtract(matrix[i],matrix[j], coeff)
    res = Fraction(1,1)
    for i in range(N):
        res *= matrix[i][i]
    return res



#Gauss method with Rational numbers, gives exact result
def Gauss(matrix, b, x): 
    matr = copy.deepcopy(matrix)
    newb = copy.deepcopy(b)
    N = len(matrix)
    for j in range(N):
        for i in range(j+1,N):
            coeff = matr[i][j] / matr[j][j]
            for k in range(j,N):
                matr[i][k] = matr[i][k] - coeff * matr[j][k]
            #matrix[i] = subtract(matrix[i],matrix[j], coeff)           #"for k" in one string
            newb[i] = newb[i] - coeff * newb[j]
    for i in range(N-1,-1,-1):
        sum = Fraction(0,1)
        for j in range(i+1,N):
            sum = sum + matr[i][j] * x[j]
        x[i] = (newb[i] - sum) / matr[i][i]      
    return x

def Kramer(matrix, b , x):
    matr = copy.deepcopy(matrix)
    N = len(matr)
    for i in range(N):
        new_matrix = list(matr)
        for j in range(N):
            new_matrix[j][i] = b[j]
        x[i] = Determinant(new_matrix) / Determinant(matr)
    return x


def LU(matrix, L, U):
    N = len(matrix)
    for i in range(N):
        U[i][i] = 1
        L[i][0] = matrix[i][0]
        U[0][i] = matrix[0][i] / L[0][0]
    for i in range(1, N):
        for j in range(1, N):
            if i >= j:
                sum = Fraction(0,1)
                for k in range(j):
                    sum += L[i][k] * U[k][j]
                L[i][j] = matrix[i][j] - sum
            else:
                sum = Fraction(0,1)
                for k in range(i):
                    sum += L[i][k] * U[k][j]
                U[i][j] = (matrix[i][j] - sum) / L[i][i]
    return L, U

def LUsolution(matrix, b, L, U):
    N = len(matrix)
    y = []
    x = []
    L, U = LU(matrix, L, U)

    # SOLUTION L*y = b

    y.append( b[0] / L[0][0])
    for i in range(1, N):
        sum = 0
        for j in range(i):
            sum  += L[i][j] * y[j]
        y.append((b[i] - sum) / L[i][i])

    # SOLUTION U*x = y
    
    x.append(y[N-1])
    for i in range(1,N):
        sum = 0
        for j in range(i):
            sum += U[N - i - 1][N - j - 1] * x[j]
        x.append((y[N - i - 1] - sum) / U[ N - i - 1][ N - i - 1])
    x = x[::-1]
    return x


def SeidelMethod(matrix,b, eps):
    N = len(matrix)
    x = []
    for i in range(N):
        x.append(Fraction(0,1))
    converge = False
    while not converge:
        x_new = copy.deepcopy(x)
        for i in range(N):
            s1 = sum(matrix[i][j] * x_new[j] for j in range(i))
            s2 = sum(matrix[i][j] * x[j] for j in range(i + 1, N))
            x_new[i] = (b[i] - s1 - s2) / matrix[i][i]           
        converge = sqrt(sum((x_new[i] - x[i])**2 for i in range(N))) <= eps
        x = x_new
    for i in range(N):
        x[i] = round(x[i].numerator / x[i].denominator , len(str(eps.denominator)))
    return x


def simpleMethod(matrix, b, eps):
        prev_x = []
        N = len(matrix)
        for i in range(N):
            prev_x.append(Fraction(0,1))
        while(True):
            curr_x = []
            for i in range(N):
                curr_x.append(b[i])
                for j in range(N):
                    if i != j:
                        curr_x[i] -= matrix[i][j] * prev_x[j]
                curr_x[i] /= matrix[i][i]
            error = Fraction(0,1)
            for i in range(N):
                error += abs(curr_x[i] - prev_x[i])
            if error < eps:
                break
            prev_x = curr_x
        return prev_x
            
    
    


'''x = np.zeros(3)      
lst = np.array([[2,4,5],
       [3,-1,2],    
       [-4,1,1]], dtype='f')
b = np.array([11,4,-2], dtype='f')
GAUSS(lst, b, x)
print(lst)
print(b)
print(x)'''

L = [[Fraction(0,1),Fraction(0,1),Fraction(0,1)],
     [Fraction(0,1),Fraction(0,1),Fraction(0,1)],
     [Fraction(0,1),Fraction(0,1),Fraction(0,1)]]

U = [[Fraction(0,1),Fraction(0,1),Fraction(0,1)],
     [Fraction(0,1),Fraction(0,1),Fraction(0,1)],
     [Fraction(0,1),Fraction(0,1),Fraction(0,1)]]


matrix = [[Fraction(10,1),Fraction(1,1),Fraction(1,1)],
          [Fraction(2,1),Fraction(10,1),Fraction(1,1)],
          [Fraction(2,1),Fraction(2,1),Fraction(10,1)]]
b1 = [Fraction(12,1),Fraction(13,1),Fraction(14,1)]
x1 = [Fraction(0,1),Fraction(0,1),Fraction(0,1)]

x = []

res = Gauss(matrix, b1, x1)
print("Gauss method")
print(res)

res1 = Kramer(matrix, b1, x1)
print('\n' + "Kramer method")
print(res1)
#print(matrix)

x = LUsolution(matrix, b1, L, U)
print('\n' + "L matrix")
print(L)
print('\n' + "U matrix")
print(U)
print('\n'+"LU solution")
print(x)

x = SeidelMethod(matrix, b1, Fraction(1,100))
print('\n' + "Soidel solution")
print(x)

x = simpleMethod(matrix, b1, Fraction(1,100))
print('\n' + "Simple Method")
print(x)
#print(Determinant(matrix))