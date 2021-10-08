README file for gcFront

gcFront searches for reaction/gene knockouts that will growth-couple a metabolite of interest in a constraint-based model.

gcFront is available under a GNU general public license. In short, this is a copyleft license that gives you the freedom to use, modify and make/transmit copies of gcFront, as long as you extend those same freedoms to others who you supply with gcFront or any code derived from it. For a full list of the terms and conditions associated with the use of gcFront, see the 'License.txt' file.

If you find gcFront useful, please remember to cite the paper in which gcFront was developed.

________________________________________________________

Using gcFront:

To run gcFront, you must first install:
MATLAB + MATLAB's global optimisation toolbox
The COBRA toolbox

Links:
https://www.mathworks.com/products/matlab.html
https://www.mathworks.com/products/global-optimization.html
https://opencobra.github.io/cobratoolbox/stable/

You will also need a linear algebra solver. You probably don't HAVE to install one, since you probably already have the solver GLPK as part of the COBRA toolbox. However, while GLPK appears to work with gcFront, other programs have reported issues when using parallel processing with GLPK, and other solvers are typically faster. Thus, it is recommended that you install a different COBRA-compatible linear algebra solver.

Links to some compatible solvers with free academic licenses that we've tested:
https://www.gurobi.com/products/gurobi-optimizer/ 
https://www.ibm.com/uk-en/analytics/cplex-optimizer
See COBRA toolbox documentation for a list of all compatible solvers.

After these pre-requisites have been installed, gcFront can be run. See the file 'gcFront documentation.pdf' for a summary of how gcFront works, 'gcFront tutorial.pdf' for a walk-through of how to use gcFront, and 'TUTORIAL_SUCCINATE.m' for an example of code that uses gcFront to search for KOs that couple succinate in the E. Coli core model.

gcFront was developed in MATLAB 2017a with Gurobi 8.1.1 and the COBRA toolbox version 3.1. If gcFront doesn't work for you, try running it with these software versions.

________________________________________________________

Contents:

'License.txt' - the GNU general public license that this code is distributed under.

'gcFront_Documentation.pdf' - explanation of how gcFront works.

'README.txt' - the file you are currently reading.

gcFront:
	'gcFront.m' - the gcFront function. Use this function to start gcFront.
	'iML1515.mat' - the iML1515 metabolic model, as downloaded from the BiGG database.
	'e_coli_core.mat' - the E. coli core metabolic model, as downloaded from the BiGG database.
	Other files - functions that are used within gcFront.

Regen_MainFig1:
	'RUNME_FindAlgTimesAndSols.m' - compares the speed and designs found by gcFront and several other gc design software packages. Used to get the data for figure 1 in the gcFront paper.
	'fGetEssGenes.m' - function used by FindAlgTimesAndSols. Returns a list of E. coli essential genes found by Goodall et al. (2018).
	Other folders - contain original and modified code from FastPros, GCopt, the COBRA toolbox implementation of OptGene, and the OptPipe implementation of RobustKnock, so that they can be compared to gcFront by FindAlgTimesAndSols.

Regen_SuppFig2:
	'RUNME_SUPP_FIG_2.m' - the script used to plot product and coupling strength for two-reaction knockout pairs in e_coli_core. Used to create the plot in Supplementary figure 2 in the gcFront paper.
	Other files - functions used to carry out calculations in SUPP_FIG_2.

Tutorial_SuccSynthStrains:
	'gcFront_tutorial.pdf' - Text that gives a worked example on how to use gcFront to find growth coupled designs.
	'TUTORIAL_SUCCINATE.m' - Code that searches for growth-coupled designs in the E. coli core model. 
	'e_coli_core.mat' - the E. coli core metabolic model, as downloaded from the BiGG database.