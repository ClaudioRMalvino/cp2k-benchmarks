######################################################################
----------------------------------------------------------------------
                                README                               
----------------------------------------------------------------------

Parameter sets for neural network potentials (NNPs) for bulk water 
based on BLYP and RPBE density-functionals with and without van der 
Waals (vdW) corrections employing the Grimme D3 method. 

These files represent supplementary material for the manuscript

'How van der Waals interactions determine the properties of water'
by
Tobias Morawietz, Andreas Singraber, Christoph Dellago, and 
Joerg Behler

Electronic Mail:
tobias.morawietz@theochem.ruhr-uni-bochum.de
joerg.behler@theochem.ruhr-uni-bochum.de

######################################################################
 Hints for usage:
######################################################################

NNPs have a highly flexible non-physical functional form which can 
adapt freely in order to reproduce reference potential-energy surfaces
with very high numerical accuracy. However, due to their flexibility, 
NNPs have only limited extrapolation capabilities and their 
performance for situations that are very different from the 
configurations in the reference data set needs to be carefully 
checked to avoid unreliable predictions.

Such extrapolation situations can and should be detected by comparing 
the symmetry function values of the given configuration with the 
corresponding minimum and maximum values of each symmetry function
in the reference data set (included in the scaling.data file).
Results obtained from simulations in which the symmetry function 
are outside the range of the values of the training set for the
respective function should be treated with caution and are very 
likely to be unreliable. This is no fundamental limitation of the NNP 
method but a consequence of the specific data set that has been used 
to determine the parameters of the NNP (bulk water in a certain 
energy-volume range in the present case). 

Even in cases where the symmetry function values do not indicate any 
extrapolations, the NNPs could yield incorrect results. It is 
therefore of great importance to only apply the NNPs to situations 
that are not too far away from the reference configurations.

The NNPs presented here are based on periodic configurations of 
liquid and crystalline water (ice Ih, XI, IX, II, XIV, XV, VIII, and 
X), sampled from geometry optimizations and molecular dynamics 
simulations with maximum temperatures of 400 K. The minimum and 
maximum volumes (in angstrom^3/molecule) of the configurations 
employed in the development of the different NNPs is given in the 
following table:

-----------------------------
 NNP	   : V min  : V max
-----------------------------
 BLYP      :   7.7  :  42.8
 BLYP-vdW  :   7.7  :  36.9
 RPBE	   :   7.6  :  53.6
 RPBE-vdW  :   7.6  :  39.9
-----------------------------

The training structures have been in part derived from MD simulations
under various conditions. While many different atomic configurations 
are therefore included, the dissociation of water molecules has not 
been enforced artificially beyond strongly distorted water 
configurations that occur spontaneously in ab initio MD simulations 
under the specified conditions (although the NNP method is in 
principle able to handle the dissociation and recombination of water
molecules). Such structures are therefore not well-represented in the 
training sets underlying the current potentials and attempts to study 
proton transport or the simulation of systems involving solvated 
protons or hydroxide ions is strongly discouraged with the current 
parameter set. The extension of the potential to include such 
structures is possible and currently in progress.

######################################################################
 NNP Parameters:
######################################################################

For each NNP (BLYP, BLYP-vdW, RPBE, RPBE-vdW), three files are 
provided: scaling.data, weights.001.data, and weights.008.data.

------------------
 scaling.data     
------------------
Contains minimum, maximum and average value for each symmetry 
function.

Column # 1: Element number (1 = hydrogen, 2 = oxygen)
Column # 2: Symmetry function number (counter)
Column # 3: Minimum symmetry function value
Column # 4: Maximum symmetry function value
Column # 5: Average symmetry function value

------------------
 weights.001.data 
------------------
Contains all parameters (weight parameters and bias weights) for the 
atomic neural network describing hydrogen atoms.

Column # 1: Parameter value
Column # 2: Type of parameter (a = weight parameter, b = bias weight)
Column # 3: Parameter number
Column # 4-7: Information on the role of the parameter. In case of a 
              weight parameter connecting two nodes, four numbers are 
              given specifying the source layer and node as well as 
              the target layer and node. In case of a bias weight 
              only the target layer and node are given.

------------------
 weights.008.data 
------------------
Contains all parameters (weight parameters and bias weights) for the 
atomic neural network describing oxygen atoms.

Column # 1: Parameter value
Column # 2: Type of parameter (a = weight parameter, b = bias weight)
Column # 3: Parameter number
Column # 4-7: Information on the role of the parameter. In case of a 
              weight parameter connecting two nodes, four numbers are 
              given specifying the source layer and node as well as 
              the target layer and node. In case of a bias weight 
              only the target layer and node are given.

------------------
Symmetry function specifications
------------------
Further, the definitions of the symmetry functions are part of the 
potential. Any change in the symmetry function parameters including 
the cutoff radius must be avoided. The list of symmetry function 
parameters can be found in the supplementary material of the 
manuscript in Tables S1 - S4.

######################################################################
