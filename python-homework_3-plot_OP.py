import sys
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np


def load_csv(csv):

    my_data = np.genfromtxt(csv, delimiter=',', skip_header=1)
    return my_data


def slice_by_ax(atoms, ax, start, end):

    new_atoms = {}
    if ax == 'x':
        mask = np.logical_and(atoms[:, 2] > start, atoms[:, 2] < end)
        new_atoms = atoms[mask]
    elif ax == 'y':
        mask = np.logical_and(atoms[:, 3] > start, atoms[:, 3] < end)
        new_atoms = atoms[mask]
    elif ax == 'z':
        mask = np.logical_and(atoms[:, 4] > start, atoms[:, 4] < end)
        new_atoms = atoms[mask]

    return new_atoms


def plot_countour(atoms, path, ax):

    if ax == 'x':
        x, y, z = atoms[:, 7], atoms[:, 8], atoms[:, 9]
        x_label, y_label = 'NORMALIZED Y', 'NORMALIZED Z'
    elif ax == 'y':
        x, y, z = atoms[:, 6], atoms[:, 8], atoms[:, 9]
        x_label, y_label = 'NORMALIZED X', 'NORMALIZED Z'
    elif ax == 'z':
        x, y, z = atoms[:, 6], atoms[:, 7], atoms[:, 9]
        x_label, y_label = 'NORMALIZED X', 'NORMALIZED Y'

    
    bounds = np.arange(0, 1.01, 0.1)

    triang = tri.Triangulation(x, y)

    fig1, ax1 = plt.subplots()

    ax1.set_aspect('equal')
    tcf = ax1.tricontourf(triang, z, 10, cmap='rainbow', levels=bounds)

    cbar = fig1.colorbar(tcf, cax=None, ax=None, ticks=bounds, fraction=0.025)
    cbar.ax.set_ylabel('Order Parameter')

    ax1.tricontour(triang, z, 10, colors='k', levels=bounds, linewidths=0.5)
    ax1.set_ylim(-1, 1)
    ax1.set_xlim(-1, 1)

    # Set x, y labels
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    
    plt.show()


def main():
    # Load arguments
    # Example command: python3 plot_OP.py dump.al_GB-210-n1_494_31680_879_atoms_list.csv z 10 15
    fp = sys.argv[1]
    ax = sys.argv[2]
    slice_start = float(sys.argv[3])
    slice_stop = float(sys.argv[4])

    # Run plot
    atoms = load_csv(fp)
    atoms = slice_by_ax(atoms, ax, slice_start, slice_stop)
    plot_countour(atoms, fp, ax)


if __name__ == "__main__":
    main()
