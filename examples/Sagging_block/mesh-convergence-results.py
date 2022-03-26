#!/usr/bin/env python
# description     :File used after mesh-convergence-study.py to plot the results
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
DATA_PATH = "mesh_analysis.npy"


# Class declarations

# Function declarations
def read_data(data_file):
    with open(data_file, 'rb') as f:
        matrix = np.load(f, allow_pickle=True)
        return matrix[0], matrix[1]

def plot_results(x, y):
    fontsize = 20
    # plt.plot(x, y, "-o", color="red", markersize=20, linewidth=5)
    plt.semilogy(x, y, "-o", color="red", markersize=fontsize, linewidth=5)
    plt.ylabel("Mean Absolute Error (MAE) [m]", fontsize=fontsize)
    plt.xlabel("Number of elements (tetrahedron)", fontsize=fontsize)
    xticks = np.arange(0, 45001, 5000)
    # # xticks[0] = 2000
    plt.xticks(xticks)
    plt.xlim([0, 45001])
    # plt.yticks(np.arange(10e-7, 10e-5, 10e-6))
    plt.ylim([3*10e-9, 3*10e-6])
    plt.tick_params(axis='both', which='major', labelsize=fontsize)
    plt.grid()
    plt.show()


def main():
    nb_elements, error = read_data(DATA_PATH)
    plot_results(nb_elements, error)


# Main body
if __name__ == '__main__':
    main()
