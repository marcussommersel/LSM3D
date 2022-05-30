from matplotlib import projections
import matplotlib.pyplot as plt
import numpy as np
from skimage import measure # scikit-image library
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# reads a .txt-file and returns the signed distance field and the number of nodes in each direction
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

# takes a signed distance field and the number of nodes in each direction
# uses the marching cubes algorithm to find the zero-contour and plots this contour
def getSurface(volume, m, n, p, level=0, plot=True, filename='fig'):
    verts, faces, normals, values = measure.marching_cubes(volume, level) # marching cubes algorithm
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
        ax.view_init(elev=0., azim=0)
        ax.set_xlabel('x', fontsize=14, style='italic')
        ax.set_ylabel('y', fontsize=14, style='italic')
        ax.set_zlabel('z', fontsize=14, style='italic')
        plt.savefig(filename + '.pdf', dpi=900, format='pdf',bbox_inches='tight')
        plt.close()

# plots the particles from the particle level set method
def plotParticle(filename):
    f = open(filename + '.txt')

    x = []
    y = []
    z = []

    for line in f.readlines():
        parsedLine = line.split(',')
        x.append(float(parsedLine[0]))
        y.append(float(parsedLine[1]))
        z.append(float(parsedLine[2]))

    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(x,y,z)
    plt.xlabel('x', fontsize=14, style='italic')
    plt.ylabel('y', fontsize=14, style='italic')
    plt.savefig(filename + '.pdf', dpi=900, format='pdf',bbox_inches='tight')
    plt.close()

# main function for plotting
def main():
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman']
    path = 'figures/'

    f = open(path + 'plotTimes.txt')
    # plots signed distance field for all time steps recorded in plotTimes.txt
    for line in f.readlines():
        phi, m, n, p = readFile(path + line[:-1] + '.txt')
        getSurface(phi, m, n, p, 0, True, path + line[:-1])

    g = open(path + 'plotTimesParticle.txt')
    # plots particles for all time steps recorded in plotTimesParticle.txt
    for line in g.readlines():
        plotParticle(line[:-1])

if __name__=='__main__':
    main()
