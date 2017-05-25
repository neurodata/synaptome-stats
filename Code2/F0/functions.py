import numpy as np

def buff( x, bf ):
    return [x - bf, x + bf]

inRange = lambda x, y, bf: bf <= x[0] <= (y[0]-1-bf) and\
    bf <= x[1] <= (y[1]-1-bf) and bf <= x[2] <= (y[2]-1-bf)

def get11Bool(loc , big, bf = 5 ):
    xyz = np.around(loc, decimals=0).astype(int)
    y = [xyz[1] - bf, xyz[1] + bf + 1]
    z = [xyz[2] - bf, xyz[2] + bf + 1]
    #
    Bool = min(x) >= 0 and min(y) >= 0 and min(z) >=0 and\
        max(x) < big.shape[0] and max(y) < big.shape[1] and max(z) < big.shape[2]
    #
    if Bool:
        cube = big[(xyz[0] - bf):(xyz[0] + bf + 1):1,
                   (xyz[1] - bf):(xyz[1] + bf + 1):1,
                   (xyz[2] - bf):(xyz[2] + bf + 1):1]
        return cube

def OLDgetCube(loc, big, bf = 5):
    xyz = np.around(loc, decimals=0).astype(int)
    #
    cube = big[(xyz[0] - bf):(xyz[0] + bf + 1):1,
               (xyz[1] - bf):(xyz[1] + bf + 1):1,
               (xyz[2] - bf):(xyz[2] + bf + 1):1]
    return cube

def getCube(cube, bf, at):
    sub = cube[[slice(i - bf, i + bf + 1) for i in at]]
    return(sub)

def getCubeCenter(cube, bf):
    ata = [(i-1)/2 for i in cube.shape]
    sub = cube[[slice(i - bf, i + bf + 1) for i in ata]]
    return(sub)

def cube2vec(cube):
    vec = np.concatenate(np.concatenate(cube))
    return(vec)


def distMat3(bf):
    A = np.reshape(
        np.array([np.sqrt((i-(bf+1))**2 + (j-(bf+1))**2 + (k-(bf+1))**2)\
        for i in range(1, 2*bf+2)\
        for j in range(1, 2*bf+2)\
        for k in range(1, 2*bf+2)]),(2*bf+1, 2*bf+1, 2*bf+1))
    A[np.where(A == 0)] = 1
    return(A)


def distMat2(bf):
    A = np.reshape(np.array([np.sqrt((i-(bf+1))**2 + (j-(bf+1))**2)
    for i in range(1, 2*bf+2)
    for j in range(1, 2*bf+2)]), (2*bf+1, 2*bf+1))
    A[np.where(A == 0)] = 1
    return(A)

def wm(cube):
    M = np.where(cube == cube.max())
    ind = [[M[0][i], M[1][i], M[2][i]] for i in range(len(M))]
    return(ind)

def f0(cube, bf = 5):
    c11 = getCubeCenter(cube, bf = 5)
    d = distMat3(bf = 5)
    M = np.transpose(np.where(c11 == c11.max()))[0]
    F0 = np.sum(c11)
    #F1 = np.sum(c11/(d**2))
    F = [F0]
    return(F)

def Fall(cube):
    c11 = getCubeCenter(cube, bf = 5)
    d = distMat3(bf = 5)
    M = np.transpose(np.where(c11 == c11.max()))[0]
    F0 = np.sum(c11)
    F1 = np.sum(c11/(d**2))
    if F0 > 0:
        F2 = np.sum((c11*d)/F0)
        F3 = np.sum((c11*(d**2))/F0)
        F4 = d[np.where(c11 == c11.max())][0]
        F5 = np.sum(getCube(cube, 2, M))
        F = [F0, F1, F2, F3, F4, F5]
    else:
        F = [F0]
    return(F)

