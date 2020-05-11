# Groundwater Modeling Coding Assignment #2
# Jim Finnegan
# 1D Transport Equation
# comparison for v = 0.5 m/d and t = 200 d

from methods import finite_difference, finite_element, analytical_vfive
from matplotlib import pyplot as plt
import numpy as np

c_fd = finite_difference(0.1, 1)
c_fe = finite_element(0.1, 1)
c_an = analytical_vfive(0.1)
c_an = np.array(c_an, dtype=np.float)


