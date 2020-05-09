# 1Dtransport
## EOS-288 Groundwater Modeling, Tufts University

This repository consists of 1D solute transport models using python developed as part of the coursework for Tufts EOS-288, Groundwater Modleing. The project considers a 1-D solute transport problem of contaminant migration in a porous medium, described by the following equation:

![alt text](https://github.com/jmfinnegan12/1Dtransport/blob/master/PDE.PNG)

subject to initial condition: C(x, 0) = 0 
and boundary conditions: C(0, t) = C0 and C(inf, t) = 0

where	C is	the	solute	concentration,	D is	the	hydrodynamic	dispersion	coefficient,	v is	the	average	linear	velocity	
and	R is	the	retardation	coefficient.	The	analytical	solution	to	this	mathematical	problem	is	based	after Ogata	
and	Banks	(1961):

![alt text](https://github.com/jmfinnegan12/1Dtransport/blob/master/Analytical.PNG)

The solution is modeled using both the Finite Difference Method (FDM) and the Finite Element Method (FEM) and compared to the analytical solution. 

Source: EOS-288 class notes and homework assignment text

The code is summarized below:


### Finite Difference Code

### Finite Element Code


