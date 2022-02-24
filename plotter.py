from matplotlib import projections
import matplotlib.pyplot as plt
import numpy as np
from skimage import measure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def readFile(filename):
    f = open(filename)
    
    firstLine = f.readlines()[0].split(',')
    f.close()
    m = int(firstLine[0])
    n = int(firstLine[1])
    p = int(firstLine[2])

    f = open(filename)
    phi = np.zeros((m,n,p))
    lines = f.readlines()[1:]
    count = 0
    for k in range(p):
        for j in range(n):
            for i in range(m):
                phi[i,j,k] = float(lines[count].split(',')[3])
                count += 1
    return phi, m, n, p

def getSurface(volume, m, n, p, level=0, plot=True):
    verts, faces, normals, values = measure.marching_cubes(volume, level)
    size = max(m, n, p)
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.set_box_aspect([1,1,1])
        mesh = Poly3DCollection(verts[faces]/size)
        mesh.set_edgecolor('k')
        ax.add_collection3d(mesh)
        plt.tight_layout()
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_zlim(0, 1)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.show()

def printScatter(filename):
    f = open(filename)

    x = []
    y = []
    z = []

    count = 0
    for line in f.readlines():
        if count == 0:
            count += 1
            continue
        count += 1
        parsedLine = line.split(',')
        x.append(float(parsedLine[0]))
        y.append(float(parsedLine[1]))
        z.append(float(parsedLine[2]))

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(x,y,z)
    plt.xlabel('x')
    plt.ylabel('y')

    plt.show()

def main():
    phi, m, n, p = readFile('0.000000.txt')
    getSurface(phi, m, n, p, 0, True)
    phi, m, n, p = readFile('1.000000.txt')
    getSurface(phi, m, n, p, 0, True)
    phi, m, n, p = readFile('2.000000.txt')
    getSurface(phi, m, n, p, 0, True)
    phi, m, n, p = readFile('3.000000.txt')
    getSurface(phi, m, n, p, 0, True)
    phi, m, n, p = readFile('4.000000.txt')
    getSurface(phi, m, n, p, 0, True)
    # printScatter('file.txt')

if __name__=='__main__':
    main()