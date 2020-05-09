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
    fe_diff = np.subtract(c_fe[8, :], c_an[8, :])
    fd_diff = np.subtract(c_fd[8, :], c_an[8, :])

    x = np.linspace(0, 200, num=101)
    plt.plot(x, fe_diff, fd_diff)
    title_string = 'Solution comparison - difference from analytical\n' \
                       + 'D = ' + str(i) + ', R = 1' + ', t = 400 days'
    plt.title(title_string)
    plt.xlabel('distance (m)')
    plt.ylabel('difference in C/C0')
    plt.legend(['FE', 'FD'])
    plt.show()


