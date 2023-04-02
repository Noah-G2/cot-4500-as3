import numpy as np

def f(x,y):
    return x- (y*y)

def EulerMethod(xLower,y0,xUpper,n):
    
    h = (xUpper-xLower)/n
    
    for i in range(n):
        slope = f(xLower, y0)
        yn = y0 + h * slope
        y0 = yn
        xLower = xLower+h
    
    print("%.5f" %yn)
    print()
    
def RungeKuttaMethod(x0,y0,xn,n):
    
    h = (xn-x0)/n
    
    for i in range(n):
        k1 = h * (f(x0, y0))
        k2 = h * (f((x0+h/2), (y0+k1/2)))
        k3 = h * (f((x0+h/2), (y0+k2/2)))
        k4 = h * (f((x0+h), (y0+k3)))
        k = (k1+2*k2+2*k3+k4)/6
        yn = y0 + k
        y0 = yn
        x0 = x0+h
    
    print("%.5f" %yn)
    print()
    
def GaussianElimination(Matrix):
    for i in range(len(Matrix)):
        max_row = i
        for j in range(i+1, len(Matrix)):
            if abs(Matrix[j][i]) > abs(Matrix[max_row][i]):
                max_row = j
        Matrix[i], Matrix[max_row] = Matrix[max_row], Matrix[i]
        pivot = Matrix[i][i]
        for j in range(i, len(Matrix[i])):
            Matrix[i][j] /= pivot
        for j in range(len(Matrix)):
            if j != i:
                factor = Matrix[j][i]
                for k in range(i, len(Matrix[i])):
                    Matrix[j][k] -= factor * Matrix[i][k]
    x = [row[-1] for row in Matrix]
    
    x = [int(o) for o in x]
    print(x)
    
    print()
    
def lu_decomposition(matrix):
    n = len(matrix)
    L = np.zeros((n, n))
    U = np.zeros((n, n))
    for i in range(n):
        L[i][i] = 1
        for j in range(i, n):
            U[i][j] = matrix[i][j]
            for k in range(i):
                U[i][j] -= L[i][k] * U[k][j]
        for j in range(i + 1, n):
            L[j][i] = matrix[j][i]
            for k in range(i):
                L[j][i] -= L[j][k] * U[k][i]
            L[j][i] /= U[i][i]
    
    return L, U

def determinant(matrix):
    L, U = lu_decomposition(matrix)
    det = np.prod(np.diag(U))
    return det    

def isDiagonallyDominantMatrix(A):
    for i, row in enumerate(A):
        s = sum(abs(v) for j, v in enumerate(row) if i != j)
        if s > abs(row[i]):
            return False
    return True

if __name__ == "__main__":

    xLower = 0
    xUpper = 2
    InitialPoint = 1
    Iteration = 10

    EulerMethod(xLower,InitialPoint,xUpper,Iteration)

    RungeKuttaMethod(xLower,InitialPoint,xUpper,Iteration)

    Matrix3 = [[2, -1, 1, 6],
             [1, 3, 1, 0],
         [-1, 5, 4, -3]]

    GaussianElimination(Matrix3)

    Matrix4 = np.array([[1, 1, 0, 3],
              [2, 1, -1, 1],
              [3, -1, -1, 2],
              [-1, 2, 3, -1]])

    det = determinant(Matrix4)

    L, U = lu_decomposition(Matrix4)

    Matrix5 = [[9,0,5,2,1],
     [3,9,1,2,1],
     [4,2,3,12,2],
     [3,2,4,0,8]]

    print()
    print("%.5f" %det)
    print()
    print(L)
    print()
    print(U)
    print()
    print(isDiagonallyDominantMatrix(Matrix5))
    print()

    arr = np.array([[2,2,1],[2,3,0],[1,0,2]])
    res = np.all(np.linalg.eigvals(arr) > 0)

    print(res)
