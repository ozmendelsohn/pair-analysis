# pair-analysis
In this repository, I have uploaded my method of detecting the stacking fault in a binary alloys system in MD (Molecular Dynamics) simulations.

This method is based on counting the number of asymmetric bond for each atom (A-B & B-A) and assigning it as a "Pair" parameter. And by analyzing this the distribution of the pair parameter of the nearest neighbors, I was able to detect the different stacking fault in the superlattice.

My specific implementation is about detecting the different stacking fault in the L1<sub>2</sub> superlattice system. By applying this method, I was able to identify three different kinds of stacking faults in a compression test of rounded edges Wulff shape nanoparticles from an MD simulation.

The distribution of the pair parameter for each stacking fault type for B<sub>3</sub>A L1<sub>2</sub> alloy:

|      Fault Type      | Structure Type | Particle type | Pair Parameter | Distribution of the Nighbor Pair Parameter [Particle type: Pair Parameter - Count of Each Type] |
|:------ | :------: | :------: | :------: | :------ |
| 1: Anti-Phase Boundaries (ABP)| FCC (1) | A | 11 | A:11 - 1 ; B:3 - 2 ; B:4 - 9 |
| 1: Anti-Phase Boundaries (ABP)| FCC (1) | B | 3 | A:11 - 2 ; A:12 - 2 ; B:3 - 1 ; B:4 - 8 |
| 1: Anti-Phase Boundariess (ABP)| FCC (1) | B | 4 | A:11 - 3 ; A:12 - 2 ; B:3 - 3 ; B:4 - 5 |
| 2: Complex Stacking Faults (CSF)| HCP (2) | A | 11 | - |
| 2: Complex Stacking Faults (CSF)| HCP (2) | B | 3 | A:11 - >1; B:4 - >4 |
| 2: Complex Stacking Faults (CSF)| HCP (2) | B | 4 | A:11 - >1; B:3 - >2; B:4 - >2 |
| 3: Super-Intrinsic Stacking Faults (SISF)| HCP(2) | A | 12 | - | 
| 3: Super-Intrinsic Stacking Faultss (SISF)| HCP(2) | B | 4 | A:12 - >1; B:4 - >4|
