# Groundwater Modeling Coding Assignment #2
# Jim Finnegan
# 1D Transport Equation
# comparison of methods

from methods import finite_difference, finite_element, analytical
from matplotlib import pyplot as plt
import numpy as np

# compare results for the following D-values
D = [0.1, 1]
for i in D:
    # calculate results
    c_fd = finite_difference(i, 1)
    c_fe = finite_element(i, 1)
    c_an = analytical(i)
    c_an = np.array(c_an, dtype=np.float)

    # compare results
    for j in range(1, 5):
        fe_diff = np.subtract(c_fe[(j * 10), :], c_an[(j * 10), :])
        # fe_percent = fe_diff/c_an
        fd_diff = np.subtract(c_fd[(j * 10), :], c_an[(j * 10), :])
        # fd_percent = fd_diff / c_an[(j*10), :]

        x = np.linspace(0, 200, num=101)
        plt.plot(x, fe_diff, fd_diff)
        title_string = 'Solution comparison - difference from analytical\n' \
                       + 'D = ' + str(i) + ', R = 1' + ', t = ' + str(j * 100) + ' days'
        plt.title(title_string)
        plt.xlabel('distance (m)')
        plt.ylabel('difference in C/C0')
        plt.legend(['FE', 'FD'])
        if i == 0.1:
            filename = 'comparison_' + 'D_1_t' + str(j * 100) + '.png'
        else:
            filename = 'comparison_' + 'D1_t' + str(j * 100) + '.png'
        plt.savefig(filename)
        # plt.show()
        plt.clf()

    # fd_diff = [(c_fd[10, :] - c_an[10, :]), (c_fd[20, :] - c_an[20, :]), (c_fd[30, :] - c_an[30, :]),
    #            (c_fd[40, :] - c_an[40, :])]
    # fe_diff = [(c_fe[10, :] - c_an[10, :]), (c_fe[20, :] - c_an[20, :]), (c_fe[30, :] - c_an[30, :]),
    #            (c_fe[40, :] - c_an[40, :])]
    # fd_diff = np.vstack(fe_diff)
    # fe_diff = np.vstack(fe_diff)
