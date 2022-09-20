# AB INITIO THERMODYNAMICS ANALYSIS (aiTA_stability_analysis.py)

This script can be executed if the required values are passed in the dictionary of the scripts. 
The dictionary contains the names of the structures (keys) and the respective energy values (values). 


## Installation

```
git clone https://github.com/ekhramen/ab-initio-thermodynamics-analysis.git
```

* [structures]: This folder contains the structures, whose stabilities have been evaluated in the aiTA scripts.
* [scripts]: This folder contains the scripts used to run aiTA.
This script was created and executed using Spyder 5.0.5 and Python 3.9.6.
More information on Spyder could be found here: https://www.spyder-ide.org/

## Instructions

The script could be used if the following input parameters are present:

	1. The equilibrium describes the reaction of the formation of the bimetallic species in zeolites. 
In our case, we have employed the following equation: 
        ```
	(2p - q + 2 - 3n - 2m)/4 (molecular_oxygen) + n(AlOH/zeol) + m(CuO) + (1 - n)/zeol -> CumAlnOpHq/zeol + (q + n -2)/2 (water)
        ```

        2. Make the changes to the equation in the script for calculating the Gibbs free energy formation of your structures.

	3. Make sure that the energies of all components of the equilibrium are computed at the same level of theory and are present in the same energy units.

	4. Include the tabulated data on the values of enthalpy and entropy of the gaseous species present in the equilibrium. We included the data on H2O.txt and O2.txt. 
The thermochemical data on other gaseous species could be found here:
        ```
        https://webbook.nist.gov/chemistry/fluid/
        ```
	
        5. Decide on the range of chemical potentials that you are interested in (These values generally correspond to the particular T and P values)


As an output, depending on your choice, it can produce the 3D or 2D diagrams containing the dependency of the Gibbs free energy change on the chemical potentials of water or/(and) oxygen.

