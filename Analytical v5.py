# Groundwater Modeling Coding Assignment #2
# Jim Finnegan
# 1D Transport Equation
# Analytical Solution

import numpy as np
from math import exp, sqrt
from scipy.special import erfc
from matplotlib import pyplot as plt

# initial conditions
# for R = 1
v = 0.5
D = float(input('Enter D (m^2/d): '))
L, dx = 200, 2
dist = np.linspace(2, L, num=int(L/dx))
dist = [int(x) for x in dist]
# print(dist)
t = 200

# calculate C using analytical solution for all x>0
C = np.zeros(len(dist)+1)
for x in range(len(dist)):
    try:
        C[x] = (1/2)*(exp(v*dist[x]/D)*erfc((dist[x]+v*t)/(2*sqrt(D*t)))+erfc((dist[x]-v*t)/(2*sqrt(D*t))))
    except OverflowError:
        C[x] = 0    # set C = 0 if math overflow error

# plot
dist.insert(0, 0)
plt.plot(dist, C)

title_string = 'Analytical solution\n' + 'D = ' + str(D) + ', R = 1'
plt.title(title_string)
plt.xlabel('distance (m)')
plt.ylabel('C/C0')
plt.legend(['0 days', '100 days', '200 days', '300 days', '400 days'])
plt.show()
