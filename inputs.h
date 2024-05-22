#pragma once
#include<map>
#include<string>
#include<vector>
#include"chainUpdate.h"


std::map< std::string, std::map< std::string, double>> getRateConstants();
void readInLHHWLinking(std::vector<Reaction> &reactionList);
void readInLHHWBreaking(std::vector<Reaction> &reactionList);


std::vector<double> convertStringVectortoDoubleVector(std::vector<std::string>& stringVector);
std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str);
Reaction getNextRxn(std::string line);

struct InitialConditions {
    double T;
    SystemVariables cSS0;
    std::vector<Reaction> myReactions;
    unsigned long simulationSize;
    std::map<std::string,signed long> numMoieties0; // initial total number of functional groups of each type
    double totalMonomerConcentration;
    std::vector<double> monomerWeights;
};

InitialConditions readInitialConditionsFromInput();
std::map< std::string, std::map < std::string, double > > getCapK(std::string strdate,std::string runID);

