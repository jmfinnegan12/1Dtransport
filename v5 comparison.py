# Groundwater Modeling Coding Assignment #2
# Jim Finnegan
# 1D Transport Equation
# comparison for v = 0.5 m/d and t = 200 d

from methods import finite_difference, finite_element, analytical_vfive
from matplotlib import pyplot as plt
import numpy as np

c_fd = finite_difference(1, 1)
c_fe = finite_element(1, 1)
c_an = analytical_vfive(1)
c_an = np.array(c_an, dtype=np.float)

fe_diff = np.subtract(c_fe[20, :], c_an)
fd_diff = np.subtract(c_fd[20, :], c_an)

x = np.linspace(0, 200, num=101)
plt.plot(x, fe_diff, fd_diff)
title_string = 'Solution comparison - difference from analytical\n' \
               + 'D = 1, R = 1, t = 200 days'
plt.title(title_string)
plt.xlabel('distance (m)')
plt.ylabel('difference in C/C0')
plt.legend(['FE', 'FD'])
plt.show()

