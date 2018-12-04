import numpy as np
from fractions import Fraction
import copy
from math import sqrt
import time
import random

def subtract(x, y, coeff):
    return list(map(lambda a,b: a - coeff * b, x, y))


def Determinant(matrix):
    N = len(matrix)
    for j in range(N):
        for i in range(j+1,N):
            coeff = matrix[i][j] / matrix[j][j]
            for k in range(j,N):
                matrix[i][k] -= coeff * matrix[j][k]
            #matrix[i] = subtract(matrix[i],matrix[j], coeff)
    #res = Fraction(1,1)
    res = 1.0
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
        #sum = Fraction(0,1)
        sum = 0.0
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
        #new_matrix = list(matr)
        new_matrix = copy.deepcopy(matr)
        for j in range(N):
            new_matrix[j][i] = b[j]
        x[i] = Determinant(new_matrix) / Determinant(matr)
        #x[i] = np.linalg.det(new_matrix) / np.linalg.det(matr)
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
                #sum = Fraction(0,1)
                sum = 0.0
                for k in range(j):
                    sum += L[i][k] * U[k][j]
                L[i][j] = matrix[i][j] - sum
            else:
                #sum = Fraction(0,1)
                sum = 0.0
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
        #x.append(Fraction(0,1))
        x.append(0.0)
    converge = False
    while not converge:
        x_new = copy.deepcopy(x)
        for i in range(N):
            s1 = sum(matrix[i][j] * x_new[j] for j in range(i))
            s2 = sum(matrix[i][j] * x[j] for j in range(i + 1, N))
            x_new[i] = (b[i] - s1 - s2) / matrix[i][i]           
        #converge = sqrt(sum((x_new[i] - x[i])**2 for i in range(N))) <= eps
        #error = Fraction(0,1)
        error = 0.0
        for i in range(N):
                error += abs(x_new[i] - x[i])
        if error < eps:
            converge = True
        x = x_new
    for i in range(N):
        x[i] = round(x[i], len(str(eps))-2)
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
        #prev_x.append(Fraction(0,1))
        prev_x.append(0.0)
    while(True):
        curr_x = []
        for i in range(N):
            curr_x.append(b[i])
            for j in range(N):
                if i != j:
                    curr_x[i] -= matrix[i][j] * prev_x[j]
            curr_x[i] /= matrix[i][i]
        #error = Fraction(0,1)
        error = 0.0
        for i in range(N):
            error += abs(curr_x[i] - prev_x[i])
        if error < eps:
            break
        prev_x = curr_x
    Simple_Method_time = time.time() - start_time
    return prev_x, Simple_Method_time

def check_matrix(matrix): 
    for i, line in enumerate(matrix): 
        #s = Fraction(0,1)
        s = 0.0
        for j, elem in enumerate(line):
            if i != j:
                s += abs(line[j])
        if abs(line[i]) < (s - abs(line[i])): 
            return False 
    return True

def test_matrix(N):
    matrix = []
    b = []
    for i in range(N):
        line = []
        #b.append(Fraction(random.randint(1,20)))
        b.append(random.randint(1,20))
        for j in range(N):
            if i != j:
                #line.append(Fraction(random.randint(1,20)))
                line.append(random.randint(1,20))
            else:
                #line.append(Fraction(N * 20 + random.randint(10,30)))
                line.append(N * 20 + random.randint(10,30))
        matrix.append(line)
    return matrix, b
            
    
def main(matrix, b):
    for i, line in enumerate(matrix):
        #b[i] = Fraction(float(b[i]))
        b[i] = float(b[i])
        for j, val in enumerate(line):
            #matrix[i][j] = Fraction(float(val))
            matrix[i][j] = float(val)
    res = []
    x = []
    L = []
    U = []
    time_methods = []
    N = len(matrix)
    for i in range(N):
        #x.append(Fraction(0,1))
        x.append(0.0)
        temp = []
        for j in range(N):
            #temp.append(Fraction(0,1))
            temp.append(0.0)
        L.append(temp)
        U.append(temp)
    
    avg_gauss_time = 0; avg_kramer_time = 0; avg_lu_time = 0;
    avg_seidel_time = 0; avg_simple_time = 0;
    for i in range(100):
        #N = 10
        #x = []
        #L = []
        #U = []
        #for i in range(N):
            #x.append(Fraction(0,1))
        #    x.append(0.0)
        #for i in range(N):
        #    tmp = []
        #    for j in range(N):
        #        #tmp.append( Fraction(0,1))
        #        tmp.append(0.0)
        #    L.append(tmp)
        #    U.append(tmp)

        #matrix, b = test_matrix(N)
        res, gauss_time = Gauss (matrix, b, x)
        avg_gauss_time += gauss_time
        #print("Gauss method-------" + str(gauss_time))

        #matr = copy.deepcopy(matrix)
        #b1 = copy.deepcopy(b)
        #for i,line in enumerate(matr):
        #    b1[i] = b1[i].numerator / b1[i].denominator
        #    for j,val in enumerate(line):
        #        matr[i][j] = val.numerator / val.denominator

        res, kramer_time = Kramer(matrix, b, x)
        avg_kramer_time += kramer_time
        #print('\n' + "Kramer method------" + str(kramer_time))

        x, lu_time = LUsolution(matrix, b, L, U)
        avg_lu_time += lu_time
        #print('\n'+"LU solution--------" + str(lu_time))

        x, seidel_time = SeidelMethod(matrix, b, 0.001)#Fraction(1,100))
        avg_seidel_time += seidel_time
        #print('\n' + "Seidel solution------" + str(seidel_time))

        x, simple_time = simpleMethod(matrix, b, 0.001)#Fraction(1,100))
        avg_simple_time += simple_time
        #print('\n' + "Simple Method-------" + str(simple_time))

    print("Total Gauss method time after 100 iterations = " + str(avg_gauss_time))
    print("Total Kramer method time after 100 iterations = " + str(avg_kramer_time))
    print("Total LU method time after 100 iterations = " + str(avg_lu_time))
    print("Total Seidel method time after 100 iterations = " + str(avg_seidel_time))
    print("Total Simple method time after 100 iterations = " + str(avg_simple_time) + '\n\n')

    time_methods.append(avg_gauss_time)
    time_methods.append(avg_kramer_time)
    time_methods.append(avg_lu_time)
    time_methods.append(avg_seidel_time)
    time_methods.append(avg_simple_time)

    print("Average Gauss method time after 100 iterations = " + str(avg_gauss_time / 100))
    print("Average Kramer method time after 100 iterations = " + str(avg_kramer_time / 100))
    print("Average LU method time after 100 iterations = " + str(avg_lu_time / 100))
    print("Average Seidel method time after 100 iterations = " + str(avg_seidel_time / 100))
    print("Average Simple method time after 100 iterations = " + str(avg_simple_time / 100))

    for i, line in enumerate(matrix):
        b[i] = Fraction(float(b[i]))
        for j, val in enumerate(line):
            matrix[i][j] = Fraction(float(val))
    #res = GAUSS_Rational(matrix, b, x)
    res,t = SeidelMethod(matrix, b, 0.001)
    return time_methods, res
#main()