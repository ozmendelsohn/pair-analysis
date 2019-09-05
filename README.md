# pair-analysis
In this repository, I have uploaded my method of detecting the stacking fault in a binary alloys system in MD (Molacure Dynamics) simulations.

This method is base on counting the number of asymmetric bond for each atom (A-B & B-A) and assigning it as a "Pair" parameter. And by analyzing this the distribution of the pair parameter i the nighest neighbors I am able to detect the different stacking fault in the superlattice.

My specific implementation is about detecting the different stacking fault in the L1_{2} superlattice system. By applying this method, I was able to identify three different kinds of stacking faults in a compression test of rounded edges Wulff shape nanoparticles from an MD simulation.

