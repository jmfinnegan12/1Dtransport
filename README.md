# 1Dtransport
## EOS-288 Groundwater Modeling, Tufts University

This repository consists of 1D solute transport models using python developed as part of the coursework for Tufts EOS-288, Groundwater Modleing. The project considers a 1-D solute transport problem of contaminant migration in a porous medium, described by the following equation:

![alt text](https://github.com/jmfinnegan12/1Dtransport/blob/master/readme%20photos/PDE.PNG)

subject to initial condition: C(x, 0) = 0 
and boundary conditions: C(0, t) = C0 and C(inf, t) = 0

where	C is	the	solute	concentration,	D is	the	hydrodynamic	dispersion	coefficient,	v is	the	average	linear	velocity	
and	R is	the	retardation	coefficient.	The	analytical	solution	to	this	mathematical	problem	is	based	after Ogata	
and	Banks	(1961):

![alt text](https://github.com/jmfinnegan12/1Dtransport/blob/master/readme%20photos/Analytical.PNG)

The solution is modeled using both the Finite Difference Method (FDM) and the Finite Element Method (FEM) and compared to the analytical solution. 

Source: EOS-288 class notes and homework assignment text

The code is summarized below:


### Finite Difference Method

The finite difference solution is derived by approximating each of the terms in the 1D flow and transport PDE with the Crank-Nicholson centered difference scheme:

![alt text](https://github.com/jmfinnegan12/1Dtransport/blob/master/readme%20photos/C-N%20scheme.PNG)

The solution can be simplified to a tridiagonal system:

![alt text](https://github.com/jmfinnegan12/1Dtransport/blob/master/readme%20photos/tridiag.PNG)

The system is solved directly at each time step by using the numpy.dot and the numpy.linalg.solve functions

The code accepts user inputs for values of D and R

![alt text](https://github.com/jmfinnegan12/1Dtransport/blob/master/Final%20Plots/FD%20R1%20D_1.png)

### Finite Element Method

The finite element solution employs the Galerkin method by defining the operator:

![alt text](https://github.com/jmfinnegan12/1Dtransport/blob/master/readme%20photos/operator.PNG)

with a linear trial function:

![alt text](https://github.com/jmfinnegan12/1Dtransport/blob/master/readme%20photos/trial.PNG)


The element stiffness and storage matrices are defined as follows:

![alt text](https://github.com/jmfinnegan12/1Dtransport/blob/master/readme%20photos/stiffness-storage.PNG)

The global stiffness and storage matrices are square, tridaigonal matrices that have dimension equal to the number of elements in the mesh. Assembling the global stiffness and storage matrices is a simple matter for a 1D problem.

The final finite element equation can be simplified to:

![alt text](https://github.com/jmfinnegan12/1Dtransport/blob/master/readme%20photos/FE%20equation.PNG)

This equation is solved directly using the numpy.linalg.solve function, with the boundary conditions being reset to a 'dummy equation' to keep a concentration ratio of 1 at x = 0 for each time step

The code accepts user inputs for values of D and R

![alt text](https://github.com/jmfinnegan12/1Dtransport/blob/master/Final%20Plots/FE%20R1%20D_1.png)

### Analytical Solution 

The analytical solution uses the complimentary error function (scipy.special.erfc) assuming a retardation factor (R) of 1 to solve the analytical solution equation shown above for every time step and distance element used in the numerical models. 

The code accepts user inputs for values of D and R

![alt text](https://github.com/jmfinnegan12/1Dtransport/blob/master/Final%20Plots/A%20R1%20D_1.png)

### Comparison

The comparison code iterates over all scenarios specified in the problem statement to produce and save plots of the difference in C/C0 between the analytical solution and both numerical methods. The value of C/C0 decreases to zero at distance, so a percentage error was impossible to calculate without encountering a divide by zero condition. Calculating a more complex error value such as relative percent difference was beyond the scope of this assignment

![alt text](https://github.com/jmfinnegan12/1Dtransport/blob/master/Comparison%20Plots/comparison_D_1_t400.png)
