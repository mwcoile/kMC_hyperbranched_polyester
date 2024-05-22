#pragma once
#include<iomanip>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<math.h>
#include<array>
#include<limits>
#include "boost/multi_array.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/graph/topological_sort.hpp"
#include "boost/graph/copy.hpp"
#include<boost/graph/connected_components.hpp>
#include<random>
#include<map>
#include<algorithm>
#include<chrono>
#include<type_traits>

struct MolWtVariables {
    // https://stackoverflow.com/questions/5914422/proper-way-to-initialize-c-structs
    double Mn; // number average molecular weight
    double Mw; // weight average molecular weight
    double PDI; // polydispersity index
    bool sameChain; // intrachain reaction or interchain reaction?
    bool isNewChain; // new chain formed from two monomers?
    bool monomerAdded; // does the reaction consist of a single monomer added to an existing chain?
    // for clarification on the following 3 variables, see Hiemenz Polymer Chemistry Section 1.7.1
    double sumNi; // sumN_i, i.e., total number of polymer molecules
    double sumMiNi; // Mi is mass of polymer chain i. sumMiNi then is the total polymer mass in the system
    double sumMi2Ni; // second moment of weight distribution, i.e., sum over all chains of mass of chain i squared
    double Mi_1; // product chain 1 weight
    double Mi_2; // product chain 2 weight
    double Mi_1_old; // for exchange reactions, reactant chain 1 weight
    double Mi_2_old; // for exchange reactions, reactant chain (or monomer) 2 weight
    double weightCutoff; // at what weight is something considered to be a polymer? (i.e., when will is show up in a GPC trace)
    std::string rxnType; // what type of reaction is taking place

};

MolWtVariables molecularWeightExchangeRxn(MolWtVariables mwv);
MolWtVariables molecularWeightBreakingRxn( MolWtVariables mwv);
MolWtVariables molecularWeightLinkingRxn(MolWtVariables mwv);
MolWtVariables molecularWeight(MolWtVariables mwv);
