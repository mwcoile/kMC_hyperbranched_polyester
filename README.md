# Overview

This code was developed to model the reversible polymerization of a lactone that results in a hyperbranched polyester. A publication is forthcoming describing the results. In the manuscript, we refer to monomer and chain units 1-5. In the code, the following letters are used instead: M = cEp, 1 = cE, 2 = Emp, 3 = E, 4 = Em, 5 = Ep.

To run the code from the command line using a bash shell, first run "make" to compile the file into "program", then specify reactions you would like to simulate in this system in input.txt, and then run the code using ./program. Alternatively, reactions can be specified in "inputs.cpp" but this will require recompiling the model.

About the algorithm:

This code implements the kinetic Monte Carlo algorithm. Four random numbers are used at each time step. First, r<sub>1</sub> is used to compute the time step between reaction events. It is generated on the interval (0,1]. From interarrival times in a Poisson distribution (independent events occuring at some average rate), the time step is calculated as (1/totalRate)*log(1/r<sub>1</sub>), and zero is excluded from the r<sub>1</sub> interval because that would lead to an infinite time step.

Next, r<sub>2</sub> is used to calculate the reaction event to take place. It is generated on the range (0,1].
Example of criteria: sumRates<r2*totalRate
so if reaction channel 1 has a rate of 2 and reaction channel 2 has a rate of 1 and reaction channel 3 has a rate is 4 then total rate is 7

In other words, when sumRates>=r2*totalRate we choose it

So,
* if r<sub>2</sub>*totalRate is (0,2], we choose reaction channel 1

* if r<sub>2</sub>*totalRate is (2,3], we choose reaction channel 2

* if r<sub>2</sub>*totalRate is (3,7], we choose reaction channel 3

Then, r<sub>3</sub> and r<sub>4</sub> (r<sub>4</sub> is only needed if a bimolecular reaction channel is chosen) are used to determine the specific moieties to react. These are then multiplied by the total number of reactants N to generate a selector on the range [0,N) where N is the number of moieties participating in that reaction channel. Then, [0,1), selects the first moiety, [1,2) selects the second, up to [N-1,N) selecting the Nth moiety.

# How is polymer detail stored?
* Three moiety trackers, each implemented as C++ maps that map a moiety name to the number present
    + Track overall moieties present
    + Track moieties present on a particular chain
    + Track moieties present on each monomer within chain
* Molecular weight trackers 
    + Overall number and weight average molecular weight
    + Molecular weight of each polymer chain
    + Weight of each monomer in the system including those in polymer chains 
* Structure of each chain stored in an adjacency list where monomers are represented as vertices and bonds are represented as edges (utilizes C++ Boost Graph Library)
* Each chain object includes adjacency list structure, molecular weight, moieties present, and is stored in vector of chains termed “chainPool”
* Link tracker: the identity of the two chain units participating in the reaction as well as the position of the link. Analogous to the moiety trackers, this is tracked at the chain level and the overall level.

# Using the model

As described above, this code is designed to track moieties on particular monomers. If the reactivity of a particular group, say, OH, depends on position, e.g., position 1, 2, or 3, these may be labeled as such, i.e., OH1, OH2, OH3, and their reactions may be specified separately with different kinetic parameters in the input file.

Sometimes, though, the reactivity of a particular functional group or bond may depend not only on its position but also on the state (i.e., reacted or unreacted) of the other functional groups on the monomer. In this case, each monomeric chain unit may be assigned a single functional group that serves as its name, e.g., cEp, and then subsequent reactions transform the functional group into a new one. In the convention of paragraph 1, a typical step growth reaction would be specified as "A + B -> ab, linking". In this convention, a step growth reaction might look like "AA + BB -> Aa + bB, linking".

* Model options
Specifying the model options is an important step before running the model. Default options are included in the C++ code, but can be edited. These decide how to deal with cyclization reactions as well as which analysis functions to run -- i.e., what results the user would like to have returned

trackLoopPairs decides whether the code should bother explicitly updating the number of 
trackFullStructureAsGraph; // track the full explicit topology, or only functional groups and weights of each polymer chain
Tirrell: for linear p;olymers only -- compare to Tirrell sequence length theory (linear polymers only)
Changxia: write out concentrations to file to compare with Changxia's experimental data
printWeightFractions: if true, print wt fractions for comparison with most probable distribution. Assumes all monomers have the same molecular weight
eG: if true, write end group analysis to file
trackAllLinks: if true, track not only the edge name but also the name of each monomer (stored in a singular functional group on each monomer) participating in a particular bond. For this to work, the
loopingSeparateChannel: are looping reactions considered to be separate reaction channels, or should they just be lumped into the bimolecular reactions? If true that looping is considered a separate reaction channel, then loopPairs should be tracked (particuarly below the gel point)--otherwise the code doesn't know what to select... 
};

