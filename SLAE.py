import numpy as np
from fractions import Fraction
import copy
from math import sqrt
import time

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
    start_time = time.time()
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
    Gauss_time = time.time() - start_time
    return x, Gauss_time

def Kramer(matrix, b , x):
    matr = copy.deepcopy(matrix)
    N = len(matr)
    start_time = time.time()
    for i in range(N):
        new_matrix = list(matr)
        for j in range(N):
            new_matrix[j][i] = b[j]
        x[i] = Determinant(new_matrix) / Determinant(matr)
    Kramer_time = time.time() - start_time
    return x, Kramer_time


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
    start_time = time.time()
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
    LU_time = time.time() - start_time
    return x, LU_time


def SeidelMethod(matrix,b, eps):
    if not check_matrix(matrix):
        print("The matrix does not satisfy the condition of diagonal predominance of elements. Enter the correct matrix")
        return
    N = len(matrix)
    x = []
    start_time = time.time()
    for i in range(N):
        x.append(Fraction(0,1))
    converge = False
    while not converge:
        x_new = copy.deepcopy(x)
        for i in range(N):
            s1 = sum(matrix[i][j] * x_new[j] for j in range(i))
            s2 = sum(matrix[i][j] * x[j] for j in range(i + 1, N))
            x_new[i] = (b[i] - s1 - s2) / matrix[i][i]           
        #converge = sqrt(sum((x_new[i] - x[i])**2 for i in range(N))) <= eps
        error = Fraction(0,1)
        for i in range(N):
                error += abs(x_new[i] - x[i])
        if error < eps:
            converge = True
        x = x_new
    #for i in range(N):
    #    x[i] = round(x[i].numerator / x[i].denominator , len(str(eps.denominator)))
    Seidel_time = time.time() - start_time
    return x, Seidel_time


def simpleMethod(matrix, b, eps):
    if not check_matrix(matrix):
        print("The matrix does not satisfy the condition of diagonal predominance of elements. Enter the correct matrix")
        return
    start_time = time.time()
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
    Simple_Method_time = time.time() - start_time
    return prev_x, Simple_Method_time

def check_matrix(matrix): 
    for i, line in enumerate(matrix): 
        s = Fraction(0,1)
        for j, elem in enumerate(line):
            if i != j:
                s += abs(line[j])
        if abs(line[i]) < (s - abs(line[i])): 
            return False 
    return True
            
    
def main(matrix, b):
    for i, line in enumerate(matrix):
        b[i] = Fraction(float(b[i]))
        for j, val in enumerate(line):
            matrix[i][j] = Fraction(float(val))
    res = []
    x = []
    L = []
    U = []
    N = len(matrix)
    for i in range(N):
        x.append(Fraction(0,1))
        temp = []
        for j in range(N):
            temp.append(Fraction(0,1))
        L.append(temp)
        U.append(temp)

    res, gauss_time = Gauss(matrix, b, x)
    print("Gauss method-------" + str(gauss_time))
    print(res)

    res, kramer_time = Kramer(matrix, b, x)
    print('\n' + "Kramer method------" + str(kramer_time))
    print(res)

    x, lu_time = LUsolution(matrix, b, L, U)
    print('\n' + "L matrix")
    print(L)
    print('\n' + "U matrix")
    print(U)
    print('\n'+"LU solution--------" + str(lu_time))
    print(x)

    x, seidel_time = SeidelMethod(matrix, b, Fraction(1,100))
    print('\n' + "Seidel solution------" + str(seidel_time))
    print(x)

    x, simple_time = simpleMethod(matrix, b, Fraction(1,100))
    print('\n' + "Simple Method-------" + str(simple_time))
    print(x)
