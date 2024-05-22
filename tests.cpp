#include"tests.h"


bool test1(std::map<std::string,long> numMoieties, std::vector<Chain> chainPool, std::map<std::string,std::vector<std::string>> funGroups, 
std::map<std::string, signed long> numMonomers, std::map<std::string,Monomer> monomers) {
    // maybe we should check here that all functional groups are adding up to what we expect..
    // this should all go in a function that could be used as one of many unit tests later...
    std::map<std::string, long> sumChainFunCounters;
    std::map<std::string, long> sumChainUnitFunCounters;
    std::map<std::string, long> sumMonomerFunCounters;
    // note: sumMonomerFunCounters+sumChainFunCounters == numMoieties[particular fun]
    // note 2: sumChainFunCounters == sumChainUnitFunCounters
    // compare for all elements in numMoieties
    for (auto iter=numMoieties.begin(); iter!=numMoieties.end(); ++iter) {
        // add up all elements in chainPool
        for (int i=0;i<chainPool.size();++i) {
            if (chainPool[i].funs.count(iter->first) ) {
                sumChainFunCounters[iter->first]+=chainPool[i].funs[iter->first];
            }
            for (auto vxs=boost::vertices(chainPool[i].polymer); vxs.first!=vxs.second; ++vxs.first) {
                if (chainPool[i].polymer[*vxs.first].funs.count(iter->first)) {
                    sumChainUnitFunCounters[iter->first]+=chainPool[i].polymer[*vxs.first].funs[iter->first]; //vxs.first->funs[iter->first];
                }
            }
        }

        // add up monomer functional groups
        sumMonomerFunCounters[iter->first] = monomerFunGrpCount(funGroups, iter->first, numMonomers, monomers);

        // perform check
        if (numMoieties[iter->first]!=sumChainUnitFunCounters[iter->first]+sumMonomerFunCounters[iter->first]) {
            return false;
            
        }
        if (sumChainUnitFunCounters[iter->first] != sumChainFunCounters[iter->first]) {
            return false;
            //throw std::runtime_error("sumChainUnitFunCounters!-sumChainFunCounters");
        }
    }
    return true;
}

bool test2(std::vector<Chain> chainPool) {
    // for verification, check that all polymers are still attached to one another (this should be a unit test of the code overall, tbh)
    for (int i=0;i<chainPool.size();++i) {
        if (!checkConnectivity(chainPool[i].polymer).isSingleComponent) {
            // only need to update fun group counts and molecular weight trackers
            return false;
        }
    }
    return true;
}

bool test3(SystemVariables cSS,ModelOptions myOptions) {
    // check to see whether the contents of each chain matches the tracker info in chain.links

    // if (myOptions.trackFullStructureAsGraph == false) {
    //     std::cout << "not tracking full structure --unable to complete test3" << std::endl;
    //     return true;
    // }
    // if (myOptions.trackAllLinks == false) {
    //     std::cout << "not tracking all links -- unable to complete test3" << std::endl;
    // }
    
    // // loop through every chain in chainPool
    // for (int i=0;i<cSS.chainPool.size();++i) {
    //     std::map< std::string, int> checkLinkTracker; 
    //     Graph polymer=cSS.chainPool[i].polymer;
    //     // iterate through all edges in polymer...
    //     for (auto edgs = boost::edges(polymer); edgs.first!=edgs.second; ++edgs.first) {
    //         std::string source_name;
    //         std::string target_name;
    //         Vertex src = source(*edgs.first,polymer);
    //         Vertex tgt = target(*edgs.first,polymer);
    //         // get names of vertices
    //         for (auto j=polymer[src].funs.begin();j!=polymer[src].funs.end();j++) {
    //             if (j->second != 0) {
    //                 source_name=j->first;
    //                 break;
    //             }
    //         }
    //         for (auto j=polymer[tgt].funs.begin();j!=polymer[tgt].funs.end();j++) {
    //             if (j->second != 0) {
    //                 target_name=j->first;
    //                 break;
    //             }
    //         }
    //         // get name of edge
    //         std::string edge_name=polymer[*edgs.first].name;
    //         checkLinkTracker[source_name+"-"+edge_name+"-"+target_name]++;
    //     }
    //     for (auto k=checkLinkTracker.begin();k!=checkLinkTracker.end();++k) {
    //         if (k->second != cSS.chainPool[i].linkTotals[k->first]) {
    //             throw std::runtime_error("helpppp");
    //             return false;
    //         }
    //     }
    // }

    return true;
}

bool test4(SystemVariables cSS,ModelOptions myOptions) {

    // confirm that loopPairs are being tracked correctly
    if (myOptions.trackLoopPairs==false) {
        return true; // no need to check that loop pairs are being tracked correctly
    }
    std::map < std::string, std::map < std::string , signed long long > > checkLoopPairs;
    for (int i=0;i<cSS.chainPool.size();++i) {
        // add up loop pairs on this chain
        for (int j=0;j<cSS.loopingRxns.size();++j) {
            std::string fun1=cSS.loopingRxns[j].reactants[0];
            std::string fun2=cSS.loopingRxns[j].reactants[1];
            checkLoopPairs[fun1][fun2] += cSS.chainPool[i].funs[fun1] * cSS.chainPool[i].funs[fun2];
            checkLoopPairs[fun2][fun1] = checkLoopPairs[fun1][fun2];
        }

    }
    for (std::map < std::string, std::map < std::string , signed long long > >::iterator it = checkLoopPairs.begin(); it!= checkLoopPairs.end();++it) {
        for (auto it2 = (it->second).begin(); it2!=(it->second).end();++it2) {
            if (it2->second!=cSS.loopPairs[it->first][it2->first] ) {
                throw std::runtime_error("check tracking of looping pairs");
                return false;
            }
        }
    }


    std::map < std::string , int > testMe;
    for (auto iter = testMe.begin(); iter!=testMe.end();++iter) {
        int foo = iter->second;
    }
    return true;
}

// test function to confirm that the number of allophanate moieties is 1/2 the number of allophanate edges in chainPool

// test function to confirm that number of alkenes is equal to number of amines and CO2 when no amine side reactions present

// test function to confirm that alkenes + amines + isocyanates + carbamates*2 + allophanates*3 + alcohols == total initial number of functional groups

void testCode(SystemVariables currentSystemState) {
    // this should get moved to more formal unit testing or a separate c++ file
    // also all parameters should be passed here, but I plan to rework this at some point
    // this section of code is just to ensure everything is working correctly
    // call this function at the start of every kMC loop to ensure that everything looks as expected--
    // useful for tracking down origins of bugs

    if (!test1(currentSystemState.numMoieties,currentSystemState.chainPool, currentSystemState.funGroups, currentSystemState.numMonomers, currentSystemState.monomers) ) {
        throw std::runtime_error("functional group counters don't add up!");
    }
    // if (!test2(currentSystemState.chainPool)) {
    //     throw std::runtime_error("all monomers within any given polymer chain should be connected at the start of every kMC iteration");
    // }
    return;
}

// void print_map(std::string comment, const std::map<std::string, signed long>& m)
// {
//     // this function is just for debugging, to give the ability to print a map to std::cout
//     std::cout << comment;
//     for (const auto& [key, value] : m) {
//         std::cout << key << " = " << value << "; ";
//     }
//     std::cout << "\n";
//     return;
// }

void countLinks(SystemVariables & cSS) {
    // for (auto it=cSS.links.begin(); it !=cSS.links.end(); ++it) {
    //     std::string linkName = it->first;
    //     int chainLinkCountersSum=0;
    //     int polymerLinkCountersSum=0;
    //     for (int i=0;i<cSS.chainPool.size();++i) {
    //         if (cSS.chainPool[i].linkTotals.count(linkName)) {
    //             chainLinkCountersSum+=cSS.chainPool[i].linkTotals[linkName]; 
    //         }
    //         for (auto edgs=edges(cSS.chainPool[i].polymer); edgs.first!=edgs.second; ++edgs.first) {
    //             auto polymer = cSS.chainPool[i].polymer;
    //             Vertex src=source(*edgs.first,polymer);
    //             Vertex tgt=target(*edgs.first,polymer);
    //             // determine edge name
    //             std::string srcFun="";
    //             std::string tgtFun="";
    //             // this code is specific to the Changxia case -- in the future just use monomer names instead of assuming there is only one functional group per monomer?
    //             for (auto j=polymer[src].funs.begin();j!=polymer[src].funs.end();++j) {
    //                 if (j->second != 0) {
    //                     srcFun=j->first;
    //                     break;
    //                 }
    //             }
    //             for (auto j=polymer[tgt].funs.begin();j!=polymer[tgt].funs.end();++j) {
    //                 if (j->second != 0) {
    //                     tgtFun=j->first;
    //                     break;
    //                 }
    //             }
    //             if (srcFun=="" || tgtFun=="") {
    //                 throw std::runtime_error("there is a problem. most likely stems from an issue with edge counting. ");
    //             }
    //             std::string currentLink=srcFun+'-'+cSS.chainPool[i].nmp[*edgs.first]+'-'+tgtFun;
    //             if (currentLink==linkName) {
    //                 polymerLinkCountersSum++;
    //             }
    //         }
    //     }
    //     if (chainLinkCountersSum!=cSS.links[linkName] || polymerLinkCountersSum!=cSS.links[linkName]) {
    //         throw std::runtime_error("woah smth wrong");
    //     }
    //}
    std::cout << "this function is broken" << std::endl;

    //     if (linkCounter>whichLink) {
    //         linkCounter-=cSS.chainPool[i].linkTotals[fullLinkNameToBust];

    //         // examine every edge in cSS.chainPool[i]... need to obtain edge list
    //         // https://www.boost.org/doc/libs/1_79_0/libs/graph/doc/EdgeListGraph.html
    //         auto polymer=cSS.chainPool[i].polymer;
    //         auto edgs=edges(polymer);
    //         while (edgs.first!=edgs.second) {
    //             Vertex src=source(*edgs.first,polymer);
    //             Vertex tgt=target(*edgs.first,polymer);
    //             // determine edge name
    //             std::string srcFun="";
    //             std::string tgtFun="";
    //             // this code is specific to the Changxia case -- in the future just use monomer names instead of assuming there is only one functional group per monomer?
    //             for (auto j=polymer[src].funs.begin();j!=polymer[src].funs.end();++j) {
    //                 if (j->second != 0) {
    //                     srcFun=j->first;
    //                     break;
    //                 }
    //             }
    //             for (auto j=polymer[tgt].funs.begin();j!=polymer[tgt].funs.end();++j) {
    //                 if (j->second != 0) {
    //                     tgtFun=j->first;
    //                     break;
    //                 }
    //             }
    //             if (srcFun=="" || tgtFun=="") {
    //                 throw std::runtime_error("there is a problem. most likely stems from an issue with edge counting. ");
    //             }
    //             std::string currentLink=srcFun+'-'+cSS.chainPool[i].nmp[*edgs.first]+'-'+tgtFun;
    //             if (currentLink==fullLinkNameToBust) {
    //                 linkCounter++;
    //             }
    //             if (linkCounter>whichLink) {
    //                 myEdge.chainIndex=i;
    //                 myEdge.src=src;
    //                 myEdge.tgt=tgt;
    //                 breakOutOfLoop=true;
    //                 break;
    //             }
    //             edgs.first++;
    //         }
    //         if (breakOutOfLoop) {
    //             break;
    //         }
    //     }
    // }
    // if (breakOutOfLoop==false) {
    //     throw std::runtime_error("problem! edge to break not found. probably stems from edge counting issue");
    // }
}
