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
#include"chainUpdate.h"

bool test1(std::map<std::string,long> numMoieties, std::vector<Chain> chainPool, std::map<std::string,std::vector<std::string>> funGroups, 
std::map<std::string, signed long> numMonomers, std::map<std::string,Monomer> monomers);

bool test2(std::vector<Chain> chainPool);

bool test3(SystemVariables cSS,ModelOptions myOptions);

bool test4(SystemVariables cSS,ModelOptions myOptions);

void testCode(SystemVariables currentSystemState);
void countLinks(SystemVariables & cSS);

// void print_map(std::string comment, const std::map<std::string, signed long>& m);

