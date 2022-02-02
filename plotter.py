from matplotlib import projections
import matplotlib.pyplot as plt


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
    printScatter('file.txt')

if __name__=='__main__':
    main()