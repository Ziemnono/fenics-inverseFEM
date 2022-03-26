#!/usr/bin/env python
# description     :File used after incompressible-neo-hooke.py to plot the results.
# authors         :Arnaud Mazier & Jack Hale
# contact         :mazier.arnaud@gmail.com
# date            :10/2021
# python_version  :3.8.10
# paper          :
# ==============================================================================

# Imports
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

# Global variables
DATA_PATH = "results.npy"


# Class declarations

# Function declarations
def read_data(data_file):
    with open(data_file, 'rb') as f:
        matrix = np.load(f, allow_pickle=True)
        return matrix[0]*1e-3, matrix[1]*1e-5, matrix[2]*1e5

def plot_results_3d(x, y, z):
    fontsize=20
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    # cmhot = plt.get_cmap("viridis")
    c = np.abs(z)
    ax.set_zlim(0, 1);

    ax.set_xlabel('$\mu$ [$10^{3}$ Pa]', fontsize=fontsize)
    ax.set_ylabel('$\lambda$ [$10^{5}$ Pa]', fontsize=fontsize)
    ax.set_zlabel('MAE [$10^{-5}$ m]', fontsize=fontsize)
    ax.xaxis.labelpad=fontsize
    ax.yaxis.labelpad=fontsize
    ax.zaxis.labelpad=fontsize


    ax.tick_params(axis='both', which='major', labelsize=fontsize)
    scatter = ax.scatter(x, y, z, c=c, cmap="viridis", vmin = 0, vmax = 1.5)
    fig.colorbar( scatter, shrink=1)
    plt.show()


def main():
    mu, lmbda, error = read_data(DATA_PATH)
    plot_results_3d(mu, lmbda, error)


# Main body
if __name__ == '__main__':
    main()
