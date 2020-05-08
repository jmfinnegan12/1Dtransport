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
v = 0.1
D = float(input('Enter D (m^2/d): '))
L, dx = 200, 2
dist = np.linspace(2, L, num=int(L/dx))
dist = [int(x) for x in dist]
# print(dist)
t, dt = 400, 10
time = np.linspace(10, t, num=int(t/dt))
time = [int(t) for t in time]
# print(time)

# initialize grid for C
# x is distance (one column is 2 ft), y is time (one row is 10 days)
C = np.zeros((len(time)+1, len(dist)+1))
C[:, 0] = 1     # boundary condition: C/C0 = 1 at x=0

# calculate C using analytical solution for all x>0
for x in range(len(dist)):
    for t in range(len(time)):
        C[t+1][x+1] = (1/2)*(exp(v*dist[x]/D)*erfc((dist[x]+v*time[t])/(2*sqrt(D*time[t])))+erfc((dist[x]-v*time[t])/(2*sqrt(D*time[t]))))

# plot
dist.insert(0, 0)
plt.plot(dist, C[0, :])
plt.plot(dist, C[10, :])
plt.plot(dist, C[20, :])
plt.plot(dist, C[30, :])
plt.plot(dist, C[40, :])

title_string = 'Analytical solution\n' + 'D = ' + str(D) + ', R = 1'
plt.title(title_string)
plt.xlabel('distance (m)')
plt.ylabel('C/C0')
plt.legend(['0 days', '100 days', '200 days', '300 days', '400 days'])
plt.show()
