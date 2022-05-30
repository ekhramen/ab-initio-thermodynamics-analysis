#Ab initio thermodynamics analysis

This script can be executed if the required values are passed in the dictinary of the scripts. 
The dictionary implyes the names of the structures (keys) and the values of the energies calucalted in CP2K (values). However, any eenrgy unit can be chosen with appropriate modifications to the script.

## Installation

```
git clone https://github.com/ekhramen/ab-initio-thermodynamics-analysis.git
```

This scripts was created and executed using Spyder 5.0.5 and Python 3.9.6.
More information on Spyder could be found here: https://www.spyder-ide.org/

## Instructions

The script could be used if the following input parameters are present:
	1. The equilibrium describing the reaction of formation the bimetallic species in zeolites. 
	In out case, we have employes the follwoing equation:  
	(2*p - q + 2 - 3*n - 2*m) / 4 * O2 +  n * AlOH/zeol + m * CuO + (1 - n)/zeol ->
         -> CumAlnOpHq/zeol + (q + n -2) / 2 * H2O
	Therefore, you have to change the equation in the script and do the appropriate adjustments for calculating the Gibss free energy.
	2. The calculated energies of the all components of the equilibrium computed at the same level of theory in the same energy units.
	3. Include the tabulated data on the values of enthaply and entropy of the gaseous species present in the equilibrium. We included the data on H2O.txt and O2.txt 
	The thermochemical data on other gaseous species could be found here: https://webbook.nist.gov/chemistry/fluid/
	4. Decide on the range of chemical potentials that you are interested in (They might have to correspond to the particular T and P values)


As an output, depending on your choice, it can produce the 3D or 2D diagrams containing the dependency of the Gibbs free enregy change on the chemical potentials of water or/(and) oxygen.

