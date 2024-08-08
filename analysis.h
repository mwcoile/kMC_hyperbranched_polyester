#pragma once
#include"chainUpdate.h"
#include<vector>


void print_graph(std::vector<Chain> chainPool,std::string TimePoint);

void recordChangxiaConcentrations(std::map<std::string,signed long> numMoieties, 
std::map<std::string,signed long> numMoieties0,double over_x,double simTime,
std::ofstream &ChangxiaOutput, double V, double const Na, double Mn, double Mw);

void TirrellTheory(std::vector<Chain> &chainPool,std::map<std::string,signed long> &numMoieties,
std::map<std::string,signed long> &numMoieties0, std::map<std::string,signed long> &numMonomers,
std::map<std::string,Monomer> &monomers,std::ofstream &TirrellOutput, double simTime, double Mn, 
double Mw );

//void printwfs(std::vector<Chain> chainPool,std::string strdate, std::string filename, 
////std::vector<double> monomerWeights,std::map<std::string,signed long> numMonomers,
//signed long simulationSize);
void printwfs(std::vector<Chain> chainPool,std::string strdate, std::string filename, std::vector<double> monomerWeights,std::map<std::string,signed long> numMonomers,signed long simulationSize);

void printEndGrps(const SystemVariables &cSS,std::ofstream &eGOutput, double simTime);
void printExecutions(std::vector<long> reactionsExecutedTally,std::vector<Reaction> reactionList,std::string strdate,std::string runID);
void calculateDPrw(const SystemVariables &cSS,std::ofstream &conv);
std::vector<std::string> enumerateLinks(std::vector<Reaction> reactionList, std::ofstream &linksFile);
void printAllLinksToFile(std::map<std::string, signed long> links, std::ofstream &linksFile, std::vector<std::string> allPossibleLinks,double simTime,double V, double const Na);
