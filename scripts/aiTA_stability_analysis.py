
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  8 15:29:44 2021

@author: ekhramenkova
"""

### Import required packages ###
import numpy as np
import math
np.set_printoptions(suppress=True)
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import seaborn as sns

#This dictionary uses the energies of the strucutres calculated in CP2K 6.2  
ref = {
"CH4"     :     -24.06553,
"R"       :      8.3144598, #J/K*mol
"P_0"     :      0.986923169314,
 "P_1"    :      101325, #Paskal (1atm = 101325 Pa) kg*m^(-1)*s^(-2)
"p"       :      3.14159265359,
"k_b"     :      1.38064852E-23, # kg*m^2/(s^2*K)
"h"       :      6.62607015E-34, #kg*m^2/s
"H2/MOR" :     -1733.253350871367047,
"H2O"    :     -17.220233011652674,
"O2"    :       -31.928779841938734,
"CuO"    :  np.array([-256.514242347907896,-256.604542815571165,-256.605290843659930,-256.606036036061596,-256.606433574234472])/4,
#"Cu2AlO3"    :    [-1878.7080953662, -1878.7064409702], # Cu2AlO3 m2 m4 multiplicities
#"Cu2AlO3" :    [-1878.72375324063023, 2, 1, 3, 0], #optimized cluster 234, multiplicity 2
"Cu3O3"   :    [-1924.612981, 3, 0, 3, 0],  ### the lowest model from the JPC letters CumAlnOpHq
"Cu2AlO3" :    [-1878.7080953662, 2, 1, 3, 0], #optimized cluster 234, multiplicity 2
"Cu2AlO4H":    [-1895.341404, 2, 1, 4, 1],  # Cu2AlO4H from '82' min configuration
"AlO2H3":      [-1768.175892392781634, 1, 2, 3], # energie, x = Al, y = O, z = H
"Al2O2"  :     [-1768.535512663725058, 2, 2, 0],
"AlOH"    :     [-1750.855320466132525, 1, 1, 1],
"Al2O3H2":     [-1785.826027086544173, 2, 3, 2],
"Al2O4H4" :    [-1803.101121642846238, 2, 4, 4],
"Al3O4H"  :    [-1803.45277047845456, 3, 4, 1],
"Al3O5H3" :    [-1820.7297111536829, 3, 5, 3],
"Al4O5"    :    [-1821.029194080930893, 4, 5, 0],
"Al4O6H2" :    [-1838.30959578877354, 4, 6, 2],
"119" : [-1895.323621,	-1895.332919, -1895.316322], # Cu2AlO4H m1 m3 m5 multiplicities
"283" : [-1895.338171,	-1895.311485, -1895.307079],  #Cu2AlO4H m1 m3 m5 multiplicities
"403" : [-1895.308247, 	-1895.311469, -1895.307122], #Cu2AlO4H m1 m3 m5 multiplicities
"497" : [-1895.313289,	-1895.313686, -1895.308527], #Cu2AlO4H m1 m3 m5 multiplicities
"675" : [-1895.329183,	-1895.328071, -1895.323865], #Cu2AlO4H m1 m3 m5 multiplicities
"82"  : [-1895.325136,	-1895.337373, -1895.341404], #Cu2AlO4H m1 m3 m5 multiplicities
"336" : [-1895.290082,	-1895.290587, -1895.294806], #Cu2AlO4H m1 m3 m5 multiplicities
"448" : [-1895.320813,	-1895.318501, -1895.313271],  #Cu2AlO4H m1 m3 m5 multiplicities
"633" : [-1895.317184,	-1895.322503, -1895.32258]#Cu2AlO4H m1 m3 m5 mul
}

################# AITD (Gibbs free energy vs CH4 chemical potential)#####################
#AITD (Gibbs free energy as a function of O2 and H2O chemicl potentials)
        # (2*p - q + 2 - 3*n - 2*m) / 4 * O2 +  n * AlOH/zeol + m * CuO + (1 - n)/zeol ->
        # -> CumAlnOpHq/zeol + (q + n -2) / 2 * H2O
        # 2*CuO + AlOH/MOR + 0.5*O2  -> Cu2AlO4H/MOR 
        # 2*CuO + AlOH/MOR + 0.25*O2  -> Cu2AlO3/MOR + 0.5*H2O 
class AITD:
    # Define the change in water chemical potential
    X = np.arange(-2.0, -0.1, 0.15) 
    
    # upload the tabulated data for H2O
    list_H2O = np.loadtxt("list_H2O.txt")
    H_H2O_0K, H_H2O_700, S_H2O_700 = list_H2O[0,2], list_H2O[8,2], list_H2O[8,1]
    
    # upload the tabulated data for O2
    list_O2 = np.loadtxt("list_O2.txt")
    H_O2_0K, H_O2_700, S_O2_700 = list_O2[0,2], list_O2[8,2], list_O2[8,1]

    def __init__(self,  reference, cluster):
        """
        

        Parameters
        ----------
        reference : key of the dictionary
            This is the dictionary with my data which inclused the name od the structure and respective energies
        cluster : key from the dictionary
            This is the key from the dictionary which corresponds to a particular structure

        Returns
        -------
        None.

        """
        
        
        self.reference = reference
        self.cluster = cluster
    
    def calculate_E(self): # E, x, y, z
        # (2*p - q + 2 - 3*n - 2*m) / 4 * O2 +  n * AlOH/zeol + m * CuO + (1 - n)/zeol ->
        # -> CumAlnOpHq/zeol + (q + n -2) / 2 * H2O
        DE = 27.211 * ( self.cluster[0] - ( ( self.cluster[4] + self.cluster[2] - 2 ) / 2 ) * self.reference['H2O'] \
        - ( ( 2*self.cluster[3] - self.cluster[4] + 2 - 3*self.cluster[2] - 2*self.cluster[1] ) / 4 ) * self.reference['O2'] \
        - self.cluster[2] * self.reference['AlOH'][0] \
        - ( 1 - self.cluster[2] ) * self.reference['H2/MOR'] - self.cluster[1] * self.reference['CuO'][4] )
        
        return DE
              
    def calculate_3d_G(self, DE):
        """
        

        Parameters
        ----------
        DE : float
            the reaction energy from the previous function for a particular structure. Value is given in eV.

        Returns
        -------
        DH : float
            the Gibbs free energy change in a range of water and oxygen chemical potentials.

        """
        
        X, Y = np.meshgrid(np.arange(-2.0, -0.1, 0.15) , np.arange(-1.0, -0.01, 0.077) )
        
        # add the change in the water and oxygen chemical potentials
        DG =  DE - ( ( self.cluster[4] + self.cluster[2] - 2 ) / 2 ) * X \
        - ( ( 2*self.cluster[3] - self.cluster[4] + 2 - 3*self.cluster[2] - 2*self.cluster[1] ) / 2 ) * Y 

        return DG
        
    def print_3d_diagram(self, list_DG, list_names, list_colors):
        """
        

        Parameters
        ----------
        list_DG : list of arrays
            the Gibbs free energy change 2D-arrays for a range of structures..
        list_names : list of names
            names of the structures included in the analysis.
        list_colors : list of colors
            the colors of your choice.

        Returns
        -------
        None.

        """
        
        #Create the figure
        fig = plt.figure()
        #Add the axes
        ax = fig.gca(projection='3d')

        X, Y = np.meshgrid(np.arange(-2.0, -0.1, 0.15) , np.arange(-1.0, -0.01, 0.077) )
        #Add the axes
        #Plotting the planes using 2d mesh and 
        #Iterating trough the list of structures to keep them together in one figure
        for DG, name, color in zip(list_DG, list_names, list_colors):
            surf = ax.plot_surface(Y, X, DG, color=color, alpha = 1, label= name)
            surf._facecolors2d = surf._facecolor3d
            surf._edgecolors2d = surf._edgecolor3d
            ax.set_xlabel(r'oxygen chem. pot., eV')        
            ax.set_ylabel(r'water chem. pot., eV')
            ax.set_zlabel('$ΔG, eV$')
            ax.set_title(r'$Phase$ $diagram$')
            plt.rcParams['svg.fonttype'] = 'none'
            ax.view_init(10, 10)
            lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0, prop={'size': 10})
        
        # Save 3d graph
        plt.savefig("Phase_diagram_Cu2AlO3_vs_Cu2AlO4H_3d.svg", dpi=300)
        plt.show()
        
    def calculate_2d_G(self, DE):
        # (2*p - q + 2 - 3*n - 2*m) / 4 * O2 +  n * AlOH/zeol + m * CuO + (1 - n)/zeol ->
        # -> CumAlnOpHq/zeol + (q + n -2) / 2 * H2O
        # Calculate the water chemical potential change at the following conditions:
        # 700 K and 1 atm for X definition (water chemical potential)
        X = ( ( ( AITD.H_H2O_700 - AITD.H_H2O_0K ) * 1000 ) - ( 700 * AITD.S_H2O_700 ) + \
                       self.reference["R"] * 700 *\
                math.log( 1 / self.reference["P_0"] ) ) / 96485 # eV of water chemical potentail
        # the oxygen chemical potential is defined in the following range:
        Y = np.arange(-1.0, -0.01, 0.077) # eV of oxygen chemical potential
        DH =  DE - ( ( self.cluster[4] + self.cluster[2] - 2 ) / 2 ) * X \
        - ( ( 2*self.cluster[3] - self.cluster[4] + 2 - 3*self.cluster[2] - 2*self.cluster[1] ) / 2 ) * Y 
        return DH
        
    def print_2d_diagram(self, list_DG, list_names, list_colors):
        sns.set_style('whitegrid')
        Y = np.arange(-1.0, -0.01, 0.077) # eV of oxygen chemical potential
        # calculate the oxygen chemical potential change at 700 K and 1 atm
        m_O2_700_1 = 0.5 * ( ( ( AITD.H_O2_700 - AITD.H_O2_0K ) * 1000 ) - ( 700 * AITD.S_O2_700 ) + \
               self.reference["R"] * 700 *\
               math.log( 1 / self.reference["P_0"] ) ) / 96485 # eV of water chemical potentail
        for DG, name, color in zip(list_DG, list_names, list_colors):
            plt.plot(Y, DG, color=color, label=name)
            plt.xlabel(r'oxygen chem. pot., eV')
            plt.ylabel(r'$ΔG, eV$')
            plt.rcParams['svg.fonttype'] = 'none'
            plt.xlim(-1, 0)
            plt.legend()
        # Save 2d graph
        plt.savefig("Phase_diagram_Cu2AlO3_vs_Cu2AlO4H_2d.svg", dpi=300)
        # add the value of the oxygen chemical potential corresponding to 700 K and 1 atm

        plt.axvline(x = m_O2_700_1)
        plt.show()

# These structures are included for testing and can be replaced with structures of your choice.
if __name__ == "__main__":
    ##calculate the Gibbs free enegry change (DG) in eV 
    ### for 3D diagram ###
    aitd_1 = AITD(ref, ref['Cu2AlO4H'])
    free_energies_1 = aitd_1.calculate_3d_G(aitd_1.calculate_E())
    
    aitd_2 = AITD(ref, ref['Cu2AlO3'])
    free_energies_2 = aitd_2.calculate_3d_G(aitd_2.calculate_E())
    
    ### Plot all the planes together ###
    aitd_1.print_3d_diagram([free_energies_1, free_energies_2],\
                         ['Cu2AlO4H', 'Cu2AlO3'], ['lawngreen', 'magenta'])
        
    ##calculate the Gibbs free enegry change (DG) in eV 
    ### for 2D diagram ###
    free_energies_0_5_O2_1 = aitd_1.calculate_2d_G(aitd_1.calculate_E())
    free_energies_0_5_O2_2 = aitd_2.calculate_2d_G(aitd_2.calculate_E())
    ### Plot all the line together in 2D plot ###
    aitd_2.print_2d_diagram([free_energies_0_5_O2_1, free_energies_0_5_O2_2], \
                            ['Cu2AlO4H', 'Cu2AlO3'], ['lawngreen', 'magenta'])
    























