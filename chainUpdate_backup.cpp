#include"chainUpdate.h"

long monomerFunGrpCount(std::map<std::string,std::vector<std::string>> funGroups, std::string grp, 
std::map<std::string, signed long> numMonomers, std::map<std::string,Monomer> monomers) {
    long countGrp=0; // total number of this functional group found on monomers

    if (funGroups.count( grp )) {
        for (int i=0; i<funGroups[grp].size();i++) {
            // for each monomer that contains the group, 
            // check the number of monomers remaining and the number of functional groups per monomer
            std::string monomerWithFunGroup=funGroups[grp][i];
            int funGroupsOnMonomer=monomers[monomerWithFunGroup].funs[grp];
            int funGroupsRemainingAttachedToMonomerI=numMonomers[monomerWithFunGroup]*funGroupsOnMonomer;
            countGrp+=funGroupsRemainingAttachedToMonomerI;
        }
    }
    return countGrp;
}

Connectivity checkConnectivity(Graph d) {
    // this function checks whether the input graph is connected, in which case it returns true,
    // or whether it contains two or more fragments, in which case it returns false.

    Connectivity myConnect;
    UndirectedGraph uG; // create undirected graph to check for connectivity
    copy_graph(d, uG); // there may be a more efficient way of creating an undirected graph from a directed graph
    
    myConnect.component.resize(boost::num_vertices (uG)); // declare vector that will label connected pieces ("components") of the graph with different integers (starting with 0)
    int numComponents = boost::connected_components (uG, &myConnect.component[0]); // determine what topological fragment each vertex in the polymer chain belongs to
    if (numComponents>2) {
        throw std::runtime_error( "Oh boy. Shouldn't have more than two fragments after one reverse reaction...");
    }
    if (numComponents==1) {
        myConnect.isSingleComponent=true; // if everything is still connected, return true boolean
        return myConnect;
    }
    else {
        // if the polymer is no longer connected, return false
        myConnect.isSingleComponent=false; // if everything is not still connected, return false boolean
        return myConnect;
    }

}

EdgeParticipants linearSearchForMatchingEdge(std::vector<Chain> chainPool,double whichLink, std::string linkName,bool trackFullStructureAsGraph) {
    EdgeParticipants selectedEdge{};
    int linkCounter=0; // just count up number of links explored of type linkOriginal
    //int whichChain=-1; 
    bool foundLink=false;
    for (int i=0;i<chainPool.size();++i) {
        // check for presence of link on each chain
        linkCounter+=chainPool[i].funs[linkName];
        if (linkCounter>whichLink) {
            //whichChain=i;
            selectedEdge.chainIndex=i; // chain containing link to undergo exchange reaction
            if (trackFullStructureAsGraph==false) {
                return selectedEdge;
            }
            linkCounter-=chainPool[i].funs[linkName];
            auto polymer=chainPool[i].polymer; // shorthand
            auto edgs=edges(polymer); // obtain start and end edges for the polymer
            while (edgs.first!=edgs.second) {
                // loop through all edges in the polymer; if the edge is of the right type (linkOriginal), increment link counter
                if (chainPool[i].nmp[*edgs.first]==linkName) {
                    // std::string debuggin = chainPool[i].polymer[*edgs.first].name;
                    linkCounter++;
                }
                // else if (chainPool[i].nmp[*edgs.first]!="allophanate" && chainPool[i].nmp[*edgs.first]!="NCO-OH") {
                //     throw std::runtime_error("woah this does not make sense!");
                // }
                if (linkCounter>whichLink) {
                    // the edge has been located
                    selectedEdge.src=source(*edgs.first,chainPool[i].polymer);
                    selectedEdge.tgt=target(*edgs.first,chainPool[i].polymer);
                    //selectedEdge.chainIndex=whichChain;
                    if (selectedEdge.src>10000 || selectedEdge.tgt>10000) {
                        throw std::runtime_error("check that edge selected makes sense");
                    }
                    foundLink=true;
                    break;
                } 
                edgs.first++;
            }
            if (foundLink!=true) {
                throw std::runtime_error("something wrong!");
            }
            break;
        }
    }
    return selectedEdge;
}

GroupLocation linearSearchForParticipatingGroup(const SystemVariables &cSS,const double &randomNumber, std::string groupName, bool trackFullStructureAsGraph) {
    // find and return a particular group to react
    GroupLocation chosenGroup;
    double whichGrp = randomNumber*cSS.numMoieties.at(groupName);
    int grpCounter=0;
    // first find all monomers that contain the groups
    std::map<std::string,std::vector<std::string>>::const_iterator monomersWithFunGrp = cSS.funGroups.find(groupName);
    if (monomersWithFunGrp!=cSS.funGroups.end()) {
        for (int i=0; i<monomersWithFunGrp->second.size();i++) { 
            // for each monomer that contains the group, 
            // check the number of monomers remaining and the number of functional groups per monomer
            std::string monomerWithFunGroup=cSS.funGroups.at(groupName)[i];
            int funGroupsOnMonomer=cSS.monomers.at(monomerWithFunGroup).funs.at(groupName);
            int funGroupsRemainingAttachedToMonomerI=cSS.numMonomers.at(monomerWithFunGroup)*funGroupsOnMonomer;
            grpCounter +=funGroupsRemainingAttachedToMonomerI;
            // check whether whichReactant has selected this monomer type to react
            if (grpCounter>whichGrp) {
                // If chosen to react, store monomer type selected
                chosenGroup.monomerSelected=monomerWithFunGroup; 
                break;
            }
        }
    }
    if (chosenGroup.monomerSelected=="") {
        for (int i=0; i<cSS.chainPool.size();i++) {
            // loop through each Graph object... and search for 
            std::map<std::string,signed long>::const_iterator funsCount=cSS.chainPool[i].funs.find(groupName);

            if (funsCount!=cSS.chainPool[i].funs.end()) {
                grpCounter+=funsCount->second;
            }

            // Note other options for dealing with const map, e.g. 

            // option 1: same amount of time as find when I benchmarked it
            // if (cSS.chainPool[i].funs.count( groupName )) {
            //     grpCounter+=cSS.chainPool[i].funs.at(groupName);
            // } 
            
            // option 2: considerably slower than find when I benchmarked it
            // try {
            //     grpCounter+=cSS.chainPool[i].funs.at(groupName);
            // }
            // catch (std::out_of_range) {
            //     continue;
            // }

            // check whether this is the chain to react
            if (grpCounter>whichGrp) {
                // if the moiety to react is found on this chain, store the chain number to react
                chosenGroup.whichChain=i;
                
                if (trackFullStructureAsGraph) {
                    // loop through chain and find specific monomer to react...
                    grpCounter-=cSS.chainPool[i].funs.at(groupName);
                    
                    auto vpair = boost::vertices(cSS.chainPool[i].polymer);
                    for (auto iter = vpair.first; iter!=vpair.second; iter++) {
                        // add the funs on the given monomer to the counter
                        std::map<std::string,signed long>::const_iterator funsOnMonomerCount=cSS.chainPool[i].polymer[*iter].funs.find(groupName);
                        if (funsOnMonomerCount!=cSS.chainPool[i].polymer[*iter].funs.end()) {
                            grpCounter+=funsOnMonomerCount->second;
                        }
                        // check if this monomer is the selected one
                        if (grpCounter>whichGrp) {
                            // found node to react, store its location!
                            chosenGroup.whichVertex=*iter;
                            break;
                        }
                    }
                }
                break;
            }
        }
    }
    if (grpCounter<=whichGrp) {
        std::cout << cSS.seed << "  "  << cSS.executionNumber << "  "<< std::endl;
        throw std::runtime_error ( " darn, unable to locate group...check counters");
    }
    return chosenGroup;
} 

const Edge myNullEdge = {};

std::string getChainUnitName(Monomer myMonomer) {
    // this code is specific to the Changxia case -- in the future just use monomer names instead of assuming there is only one functional group per monomer?
    std::string chainUnitName="";
    for (auto j=myMonomer.funs.begin();j!=myMonomer.funs.end();++j) {
        if (j->second != 0) {
            chainUnitName=j->first;
            break;
        }
    }
    if (chainUnitName=="") {
        throw std::runtime_error("unexpected behavior!");
    }
    return chainUnitName;
}
std::string getLinkName(Graph polymer, Edge edg) {

    // this is more complicated than it needs to be -- should use monomer names rather than assuming 1 fun group per monomer unit. alternatively, could store full link name as edge name. i'm not really sure why i haven't done this second one (i think mainly bc I think the former is a better and more robust solution, but I haven't gotten the time to implement this change as it is rather significant at this point.)

    Vertex tgt = target( edg, polymer );
    Vertex src = source( edg, polymer);
    // determine edge name
    std::string srcFun="";
    std::string tgtFun="";
    srcFun = getChainUnitName(polymer[src]);
    tgtFun = getChainUnitName(polymer[tgt]);
    std::string linkName = srcFun + '-' + polymer[edg].name + '-' + tgtFun;
    return linkName;
}

template < typename edgeIter>
void updateLinks(std::map<std::string, signed long> &links,Chain &myChain,const std::string &oldVertName, const std::string &newVertName, std::pair<edgeIter,edgeIter> edgs, const Edge &newEdge) { 
    // for reactions in which a reaction can occur between two chain units that are bonded together, and the reaction depends not only on the bond but also the identity of the chain unit (e.g., the
    // identity of the other bonds on each chain unit, this function updates the number of links between particular monomers following a reaction)
    bool isOutEdges=std::is_same<edgeIter,boost::graph_traits<Graph>::out_edge_iterator>::value; // check the directionality of the edges that we're examining
    // if out edges, source remains constant
    if (isOutEdges) {
        while (edgs.first!=edgs.second) {
            // if the edge is the same as the new edge, ignore it
            if (*edgs.first==newEdge) {
                edgs.first++;
                continue;
            }
            std::string linkType=myChain.polymer[*edgs.first].name; //cSS.chainPool[chainSelected].nmp[*out_edgs.first]; //get linkage name. OLD: namemaps[chainSelected][*out_edgs.first];
            std::string linkstgt="";
            int linkstgtVertex=target(*edgs.first,myChain.polymer);
            
            // find linkage name of linkstgt vertex. FUTURE - might just use the "monomer name" to figure this out
            // for (auto i=myChain.polymer[linkstgtVertex].funs.begin();i!=myChain.polymer[linkstgtVertex].funs.end();i++) {
            //     if (i->second != 0) {
            //         linkstgt=i->first;
            //         break;
            //     }
            // }
            // if (linkstgt=="") {
            //     throw std::runtime_error("unexpected behavior");
            // }
            linkstgt=getChainUnitName(myChain.polymer[linkstgtVertex]);
            links[oldVertName+'-'+linkType+'-'+linkstgt]--;
            myChain.linkTotals[oldVertName+'-'+linkType+'-'+linkstgt]--;
            links[newVertName+'-'+linkType+'-'+linkstgt]++;
            myChain.linkTotals[newVertName+'-'+linkType+'-'+linkstgt]++;
            if (links[oldVertName+'-'+linkType+'-'+linkstgt]<0 || myChain.linkTotals[oldVertName+'-'+linkType+'-'+linkstgt]<0){
                throw std::runtime_error("darnit!!");
            }
            edgs.first++;
        }
    }
    // if in edges, target remains constant
    else {
        while (edgs.first!=edgs.second) {
            // if the edge is the same as the new edge, ignore it
            if (*edgs.first==newEdge) {
                edgs.first++;
                continue;
            }
            std::string linkType=myChain.polymer[*edgs.first].name; // previously: cSS.chainPool[chainSelected].nmp[*in_edgs.first]; 
            std::string linkssrc="";
            int linkssrcVertex=source(*edgs.first,myChain.polymer);
            
            // find name of chain unit that is the link's src
            for (auto i=myChain.polymer[linkssrcVertex].funs.begin();i!=myChain.polymer[linkssrcVertex].funs.end();i++) {
                if (i->second != 0) {
                    linkssrc=i->first;
                    break;
                }
            }
            if (linkssrc!=""){
                links[linkssrc+'-'+linkType+'-'+oldVertName]--;
                links[linkssrc+'-'+linkType+'-'+newVertName]++;
                myChain.linkTotals[linkssrc+'-'+linkType+'-'+oldVertName]--;
                myChain.linkTotals[linkssrc+'-'+linkType+'-'+newVertName]++;
            }
            if (links[linkssrc+'-'+linkType+'-'+oldVertName]<0 || myChain.linkTotals[linkssrc+'-'+linkType+'-'+oldVertName]<0){
                throw std::runtime_error("darnit!!");
            }
            edgs.first++;
        }
    }
    return;
}

void updateAdjacentLinks(std::map<std::string, signed long> &links,Chain &myChain, Reaction rxn, EdgeParticipants modifiedEdge,const ModelOptions myOptions,const Edge &newEdge = myNullEdge) {
    if (myOptions.trackAllLinks==false) {
        return;
    }
    // declare the following for readability
    Vertex tgt=modifiedEdge.tgt;
    Vertex src=modifiedEdge.src;
    std::string reactant1 = rxn.reactants[0];
    std::string reactant2 = rxn.reactants[1];
    std::string product1 = rxn.products[0];
    std::string product2 = rxn.products[1];

    // adjust links tracker by links OUT of tgt
    auto out_edgs = boost::out_edges(tgt,myChain.polymer);
    std::string oldVertName=reactant2;
    std::string newVertName=product2;
    updateLinks(links,myChain, oldVertName, newVertName, out_edgs, newEdge);

    // adjust links tracker by links IN to tgt
    auto in_edgs=boost::in_edges(tgt,myChain.polymer);
    updateLinks(links,myChain, oldVertName, newVertName, in_edgs, newEdge);

    // adjust links tracker by links OUT of src
    out_edgs=boost::out_edges(src,myChain.polymer);
    oldVertName=reactant1;
    newVertName=product1;
    updateLinks(links,myChain, oldVertName, newVertName, out_edgs, newEdge);

    // adjust links tracker by links IN to src
    in_edgs=boost::in_edges(src,myChain.polymer);
    updateLinks(links,myChain, oldVertName, newVertName, in_edgs, newEdge);
}

void updateLinksAdjacentToOneVertex(std::map<std::string, signed long> &links,Chain &myChain, Reaction rxn, Vertex modifiedVertex,const ModelOptions myOptions) {
    if (myOptions.trackAllLinks==false) {
        return;
    } 
    // declare the following for readability
    std::string oldVertexName = rxn.reactants[0];
    std::string newVertexName = rxn.products[0];

    // adjust links tracker by links OUT of modifiedVertex
    auto out_edgs=boost::out_edges(modifiedVertex,myChain.polymer);
    updateLinks(links,myChain, oldVertexName, newVertexName, out_edgs, myNullEdge);

    // adjust links tracker by links IN to modifiedVertex
    auto in_edgs=boost::in_edges(modifiedVertex,myChain.polymer);
    updateLinks(links,myChain, oldVertexName, newVertexName, in_edgs, myNullEdge);
}

bool isCatalystBoundForR1R2R3Linking(std::map<std::string, signed long> myfuns) {
    return (myfuns["E1*"]>0 || myfuns["Em*p"]>0 || myfuns["E3*"]>0);
}

std::vector<PossibleLoopingRxn> enumeratePossibleLoops(SystemVariables &cSS) {

    // create a vector of candidate looping pairs identified. Looping Rxn is < pair identity | topological distance | index of chain in chainPool | std::pair vertices (if less than dmax) >
    std::vector<PossibleLoopingRxn> candidateLoopingPairs;

    for (int p = 0; p<cSS.chainPool.size(); ++p) {
        // no sense in interrogating chains that are too short to loop
        int DP_p = boost::num_vertices(cSS.chainPool[p].polymer); // degree of polymerization of chain p
        if (DP_p<cSS.minLoopSize) {
            continue;
        }
        // store counts of Em*p, E1* and E3* identified within distance dmax of cE (v1)
        std::map<std::string,int> numCatAlcoholFound = {{"Em*p",0},{"E1*",0},{"E3*",0}};
        Vertex v1 = -1; // store location of cE. If v1 is -1, then no cE was found

        if (cSS.chainPool[p].funs["cE"] == 1 && isCatalystBoundForR1R2R3Linking(cSS.chainPool[p].funs)) {
            // find cE as the starting vertex
            for (auto q = vertices(cSS.chainPool[p].polymer); q.first!=q.second; ++q.first) {
                if ("cE" == getChainUnitName(cSS.chainPool[p].polymer[*q.first]) ) {
                    v1 = *q.first;
                    break;
                }
            }
            if (v1 == -1) {
                throw std::runtime_error("something wrong with counting of cE's...");
            }
            // have v1, now need to find E1*, Em*p, E3* within distance d of v1. 
            int dmax = cSS.dmax; // dmax is the maximum topological distance from v1 to search for E1*, Em*p, E3*
            // example: for a linear chain of length 5, cE (i.e., v1) is distance 0 from itself and the other end of the chain is distance d 4 from cE (i.e., v1)

            // Create a vector of length dmax, where each entry is a vector of vertices
            std::vector< std::vector< Vertex> > connectedVerticesToCheck(dmax+1, std::vector<Vertex>(0));
            connectedVerticesToCheck[0].push_back(v1); // push v1 onto the vector of vertices to check

            for (int d=1; d<=dmax; ++d ) {
                // loop over all vertices in connectedVerticesToCheck[d-1]
                for (int v=0; v<connectedVerticesToCheck[d-1].size(); ++v) {
                    // get all the vertices connected to the vth vertex
                    auto in_edgs = boost::in_edges(connectedVerticesToCheck[d-1][v],cSS.chainPool[p].polymer);
                    while (in_edgs.first!=in_edgs.second) { 
                        Vertex src = source(*in_edgs.first,cSS.chainPool[p].polymer);
                        // push this vertex onto the vector of vertices to check, provided that d<dmax
                        if (d!=dmax) {
                            connectedVerticesToCheck[d].push_back(src);
                        }
                        std::string unitName = getChainUnitName(cSS.chainPool[p].polymer[src]);
                        if (unitName == "cE") {
                            std::runtime_error("you've done something wrong!");
                        }

                        if (unitName == "E1*" || unitName == "E3*") {
                            // increment the count of E1* or E3* identified within distance dmax of cE (v1)
                            numCatAlcoholFound[unitName]+=1;
                            if (d-1>cSS.minLoopSize) {
                                // if the topological distance is greater than the minimum loop size, then it is a candidate looping pair
                                candidateLoopingPairs.push_back({unitName,d,p, std::make_pair(v1,src)});
                            }
                        }

                        if (unitName == "Em*p") {
                            // when there is a Em*p that can loop with cE, then there are two possible looping pairs -- m and p.
                            
                            // increment the count of Em*p identified within distance dmax of cE (v1) by 1, since 1 Em*p was found
                            numCatAlcoholFound[unitName]+=1;
                            if (d-1>cSS.minLoopSize) {
                                // if the topological distance is greater than the minimum loop size, then it is a candidate looping pair
                                candidateLoopingPairs.push_back({"Em*p_m",d,p, std::make_pair(v1,src)}); // account for the case of the m hydroxyl reacting
                                candidateLoopingPairs.push_back({"Em*p_p",d,p, std::make_pair(v1,src)}); // account for the case of the p hydroxyl reacting
                            }
                        }

                        ++in_edgs.first; // increment the iterator
                    }
                } 
                
            }
            // in the future, add something like the below lines to consider looping pairs that are not within distance dmax of cE (v1). For now, exclude these larger loops from forming
            // // calculate number of candidate looping pairs a distance d > dmax from cE (v1)
            // for (auto numCatAlcoholFound_it = numCatAlcoholFound.begin(); numCatAlcoholFound_it!=numCatAlcoholFound.end(); ++numCatAlcoholFound_it) {
            //     while (cSS.chainPool[p].funs[numCatAlcoholFound_it->first] - numCatAlcoholFound_it->second > 0) {
            //         // note in the line below that dmax+1 is used to indicate that the topological distance is greater than dmax. Also note that the second vertex is -1 to indicate that it is not known
            //         candidateLoopingPairs.push_back({numCatAlcoholFound_it->first,dmax+1,p, std::make_pair(v1,-1)});
            //         numCatAlcoholFound_it->second+=1;
            //     }
            // }
        }
    }
    return candidateLoopingPairs;
}

void linkingRxn(SystemVariables &cSS, Reaction rxn,double r3,double r4, ModelOptions myOptions,std::map<std::string,DoubleSec>& timers) {
    // this function carries out the polymerization linking reaction that forms a bond between a monomer
    // containing reactant 1 and a monomer containing reactant 2.
    std::string reactant1 = rxn.reactants[0];
    std::string reactant2 = rxn.reactants[1];
    std::string productLinkage;
    if (myOptions.trackAllLinks == false ) {
        productLinkage = rxn.products[0]; 
    }
    else if (myOptions.trackAllLinks == true && rxn.linkFormed != "") {
        productLinkage = rxn.linkFormed;
    }
    else {
        throw std::runtime_error("If we are tracking links formed, the type of link formed should be specified...");
    }
    
    
    // step 1: identify which two specific species containing these moeities to react together!
    // auto start = Time::now(); // start clock
    GroupLocation selectedReactant1 = linearSearchForParticipatingGroup(cSS,r3,reactant1,myOptions.trackFullStructureAsGraph);
    GroupLocation selectedReactant2 = linearSearchForParticipatingGroup(cSS,r4,reactant2,myOptions.trackFullStructureAsGraph);
    // auto stop = Time::now(); // stop clock
    // timers["linearSearchForParticipatingGroup"] += (stop - start);
    // use the following for shorthand so that lines of code don't get too long
    int whichChain1=selectedReactant1.whichChain; // index of selected group in chainPool
    int whichChain2=selectedReactant2.whichChain; // index of selected group in chainPool
    Vertex whichVertex1 = selectedReactant1.whichVertex; // vertex/monomer containing group in chainPool
    Vertex whichVertex2 = selectedReactant2.whichVertex; // vertex/monomer containing group in chainPool
    std::string monomerSelected1=selectedReactant1.monomerSelected; // if not in chainPool, type of monomer to react
    std::string monomerSelected2=selectedReactant2.monomerSelected; // if not in chainPool, type of monomer to react
    
    // if two monomers are chosen, react them and form a new chain
    if (whichChain1==-1 && whichChain2==-1) {
        // auto startMonomerLinkingTime = Time::now();

        // create the new polymer
        Chain newChain;
        Graph newPolymer; 
        boost::property_map<Graph, std::string EdgeData::*>::type namemap=get(&EdgeData::name,newPolymer); // following example https://stackoverflow.com/questions/3527313/using-bundled-properties-as-the-weight-map-in-dijkstra-shortest-paths
        if (myOptions.trackFullStructureAsGraph) {
            // combine the two monomers into one chain
            Vertex v1 = add_vertex(newPolymer);
            newPolymer[v1]=cSS.monomers[monomerSelected1];
            Vertex v2 = add_vertex(newPolymer);
            newPolymer[v2]=cSS.monomers[monomerSelected2];
            // add edge and react functional groups
            newPolymer[v1].funs[reactant1]--;
            newPolymer[v2].funs[reactant2]--;
            // track product formation on the *monomer* fun group counters
            newPolymer[v1].funs[rxn.products[0]]++;
            // if there is a second group, increment that accordingly too
            if (rxn.products.size()==2) {
                newPolymer[v2].funs[rxn.products[1]]++;
            }         
            Edge newEdge = add_edge(v1,v2,newPolymer).first;
            namemap[newEdge]=productLinkage;
            
            // record the new link... specific for changxia's system
            if (rxn.products.size()==2 && myOptions.trackAllLinks) {
                std::string product1=rxn.products[0];
                std::string product2=rxn.products[1];
                if (cSS.links.count(product1+'-'+productLinkage+'-'+product2) ) {
                    cSS.links[product1+'-'+productLinkage+'-'+product2]++; // increment   
                }
                else {
                    cSS.links[product1+'-'+productLinkage+'-'+product2]=1;
                }
            } 
        }
        // update linkages in EdgeData (make linkage map for updating EdgeData later)
        
        // decrease reacting functional group *overall* counters and monomer counts
        cSS.numMoieties[reactant1]--;
        cSS.numMoieties[reactant2]--;
        
        // increment *overall* counters of product moieties
        cSS.numMoieties[rxn.products[0]]++;
        if (rxn.products.size()==2) {
            cSS.numMoieties[rxn.products[1]]++;
        }
        // decrement monomer counts
        cSS.numMonomers[monomerSelected1]--;
        cSS.numMonomers[monomerSelected2]--;
        if (monomerSelected1=="" || monomerSelected2=="") {
            throw std::runtime_error( "Ruh roh! This isn't right" );
        }
        
        // insert structure for the new chain
        newChain.polymer =newPolymer;

        /* Update functional groups counters on the *chain*
        Functional groups on the new chain are the sum of the functional groups 
        on the two monomers, minus the two reacted funs, plus any products that form. 
        Note that if the linkage is assigned as a functional group to more than one 
        monomer then this line no longer holds true.
        */
        newChain.funs=cSS.monomers[monomerSelected1].funs;
        if (rxn.products.size()==2) {
            // specifically for Changxia's system
            std::string product1=rxn.products[0];
            std::string product2=rxn.products[1];
            newChain.linkTotals[product1+'-'+productLinkage+'-'+product2]=1; 
        }
        for(std::map<std::string,signed long>::iterator it = cSS.monomers[monomerSelected2].funs.begin(); it != cSS.monomers[monomerSelected2].funs.end(); ++it) newChain.funs[it->first] += it->second;
        newChain.funs[reactant1]--;
        newChain.funs[reactant2]--;
        newChain.funs[rxn.products[0]]++;
        if (rxn.products.size()==2) {
            newChain.funs[rxn.products[1]]++;
        }

        // update looping pairs
        if (myOptions.trackLoopPairs) {
            for (int i=0;i<cSS.loopingRxns.size();++i) {
                std::string fun1=cSS.loopingRxns[i].reactants[0];
                std::string fun2=cSS.loopingRxns[i].reactants[1];
                if (newChain.funs.count(fun1) && newChain.funs.count(fun2)) {
                    cSS.loopPairs[fun1][fun2]+=newChain.funs[fun1]*newChain.funs[fun2];
                    cSS.loopPairs[fun2][fun1]=cSS.loopPairs[fun1][fun2];
                }
            }
        }
        // update molecular weight
        cSS.mwv.Mi_1=cSS.monomers[monomerSelected1].weight;
        cSS.mwv.Mi_2=cSS.monomers[monomerSelected2].weight;
        cSS.mwv.isNewChain=true;

        // weight of the new chain is the sum of the monomers (minus any products that may have been released)
        newChain.weight=cSS.monomers[monomerSelected1].weight+cSS.monomers[monomerSelected2].weight;

        // finally, add chain to chainPool
        newChain.nmp=namemap;
        cSS.chainPool.push_back(newChain);
        // auto stopMonomerLinkingTime = Time::now();
        // timers["monomer-monomerLinking"]+=stopMonomerLinkingTime-startMonomerLinkingTime;
    }

    // if one monomer and one chain are chosen to react
    else if ((whichChain1==-1 && whichChain2!=-1) || (whichChain1!=-1 && whichChain2==-1)) {
        // auto startMonomerChainLinkingTime = Time::now();
        // these lines specifically apply to changxia's system
        Vertex src;
        Vertex tgt;
        
        // figure out which is the monomer and which is the chain
        int chainSelected;
        Vertex vertexSelected; //boost::range_detail::integer_iterator<unsigned long>
        std::string monomerSelected;
        std::string monomerFun;
        std::string chainFun;
        bool chainIsWhichVertex1; // keep track of which vertex corresponded to the monomer, so that products can be assigned appropriately
        if (whichChain1!=-1) {
            // if the chain is selected for reactant1, track the index of the chain and vertex
            chainSelected=whichChain1;
            
            if (myOptions.trackFullStructureAsGraph) {
                vertexSelected=whichVertex1;
                chainIsWhichVertex1=true;
            }
            // if the functional group reactant1 is on a chain, that means reactant2 must be on a monomer
            monomerSelected=monomerSelected2; // monomer type selected to react
            // track which functional groups to decrement 
            chainFun=reactant1;
            monomerFun=reactant2;
            
        }
        else if (whichChain2!=-1) {
            // if the chain is selected for reactant2, track the index of the chain and vertex
            chainSelected=whichChain2;
            if (myOptions.trackFullStructureAsGraph) {
                vertexSelected=whichVertex2;
                chainIsWhichVertex1=false;
            }
            // if the functional group reactant2 is on a chain, that means reactant1 must be on a monomer
            monomerSelected=monomerSelected1; // monomer type selected to react
            // track which functional groups to decrement 
            chainFun=reactant2;
            monomerFun=reactant1;
        }
        else {
            // if neither one is the chain to react, then something went wrong
            throw std::runtime_error( "Ruh roh! Check the code's logic" );
        }

        // update looping pairs (part 1)
        if (myOptions.trackLoopPairs) {
            for (int i=0;i<cSS.loopingRxns.size();++i) {
                std::string fun1=cSS.loopingRxns[i].reactants[0];
                std::string fun2=cSS.loopingRxns[i].reactants[1];
                if (cSS.chainPool[chainSelected].funs.count(fun1) && cSS.chainPool[chainSelected].funs.count(fun2)) {
                    cSS.loopPairs[fun1][fun2]-=cSS.chainPool[chainSelected].funs[fun1]*cSS.chainPool[chainSelected].funs[fun2];
                    cSS.loopPairs[fun2][fun1]=cSS.loopPairs[fun1][fun2];
                }
            }
        }
        
        Edge newEdge; // declare it; assign it later in an if statement because directionality matters; need it later when counting links in Changxia's case
        if (myOptions.trackFullStructureAsGraph) {
            // execute the reaction by (1) adding the vertex to the polymer graph and (2) adding the edge
            Vertex v1 = add_vertex(cSS.chainPool[chainSelected].polymer); // add vertex to the respective polymer graph structure
            cSS.chainPool[chainSelected].polymer[v1]=cSS.monomers[monomerSelected]; // populate vertex with the monomer chosen to be added
            
            if (monomerSelected==""){
                throw std::runtime_error( "Ruh roh! This isn't right" );
            }
            if (chainIsWhichVertex1) {
                newEdge = add_edge(vertexSelected,v1,cSS.chainPool[chainSelected].polymer).first; // add the edge
                // save src and target vertex information
                src=vertexSelected;
                tgt=v1;

            }
            else {
                newEdge = add_edge(v1, vertexSelected,cSS.chainPool[chainSelected].polymer).first; // add the edge
                tgt=vertexSelected;
                src=v1;
            }
            cSS.chainPool[chainSelected].nmp[newEdge]=productLinkage; // update name of new linkage (OLD: namemaps[chainSelected][newEdge]=productLinkage; )

            // update moiety count on the individual *monomers*
            cSS.chainPool[chainSelected].polymer[v1].funs[monomerFun]--;
            cSS.chainPool[chainSelected].polymer[vertexSelected].funs[chainFun]--;
            // update product count on *monomers* corresponding to appropriate original functional groups
            if (chainIsWhichVertex1) {
                cSS.chainPool[chainSelected].polymer[vertexSelected].funs[rxn.products[0]]++;
                if (rxn.products.size()==2) {
                    cSS.chainPool[chainSelected].polymer[v1].funs[rxn.products[1]]++;
                }
            }
            else if (!chainIsWhichVertex1) {
                cSS.chainPool[chainSelected].polymer[v1].funs[rxn.products[0]]++;
                if (rxn.products.size()==2) {
                    cSS.chainPool[chainSelected].polymer[vertexSelected].funs[rxn.products[1]]++;
                }
            }
            
        }

        cSS.numMonomers[monomerSelected]--; // update the number of monomers tracker
        if (monomerSelected=="") {
            throw std::runtime_error( "Ruh roh! This isn't right" );
        }

        // update moiety count on the individual *chain* by part (1) and part (2) below
        // part (1): update functionals that have changed due to additional monomer
        for(auto it = cSS.monomers[monomerSelected].funs.begin(); it != cSS.monomers[monomerSelected].funs.end(); ++it) cSS.chainPool[chainSelected].funs[it->first] += it->second; 
        // part (2): update moiety count due to reaction
        cSS.chainPool[chainSelected].funs[chainFun]--;
        cSS.chainPool[chainSelected].funs[monomerFun]--;
        cSS.chainPool[chainSelected].funs[rxn.products[0]]++;
        if (rxn.products.size()==2) {
            cSS.chainPool[chainSelected].funs[rxn.products[1]]++;
        }

        // update number average and weight average molecular weights
        cSS.mwv.Mi_1=cSS.monomers[monomerSelected].weight;
        cSS.mwv.Mi_2=cSS.chainPool[chainSelected].weight; // important: this only applies BEFORE the tracker is updated
        cSS.mwv.monomerAdded=true;

        // update chain weight
        cSS.chainPool[chainSelected].weight+=cSS.monomers[monomerSelected].weight;

        // update *overall* moiety count
        cSS.numMoieties[reactant1]--; // reactant 1
        cSS.numMoieties[reactant2]--; // reactant 2
        cSS.numMoieties[rxn.products[0]]++; // product 1
        if (rxn.products.size()==2) {
            cSS.numMoieties[rxn.products[1]]++;  // product 2
        }

        // update looping pairs (part 2)
        if (myOptions.trackLoopPairs) {

            for (int i=0;i<cSS.loopingRxns.size();++i) {
                std::string fun1=cSS.loopingRxns[i].reactants[0];
                std::string fun2=cSS.loopingRxns[i].reactants[1];
                if (cSS.chainPool[chainSelected].funs.count(fun1) && cSS.chainPool[chainSelected].funs.count(fun2)) {
                    cSS.loopPairs[fun1][fun2]+=cSS.chainPool[chainSelected].funs[fun1]*cSS.chainPool[chainSelected].funs[fun2];
                    cSS.loopPairs[fun2][fun1]=cSS.loopPairs[fun1][fun2];
                }
            }
        }

        // need to be tracking full structure as graph...
        // update links tracker with the new number of linkages of each type
        if (rxn.products.size()==2 && myOptions.trackFullStructureAsGraph && myOptions.trackAllLinks) {
            std::string product1=rxn.products[0];
            std::string product2=rxn.products[1];

            // first, increment links tracker by new link produced 
            if (cSS.links.count(product1+'-'+productLinkage+'-'+product2) ) {
                cSS.links[product1+'-'+productLinkage+'-'+product2]++; // increment 
            }
            else {
                cSS.links[product1+'-'+productLinkage+'-'+product2]=1;
            }
            // also, increment chain fun groups tracker by new link produced 
            if (cSS.chainPool[chainSelected].linkTotals.count(product1+'-'+productLinkage+'-'+product2) ) {
                cSS.chainPool[chainSelected].linkTotals[product1+'-'+productLinkage+'-'+product2]++; // increment 
            }
            else {
                cSS.chainPool[chainSelected].linkTotals[product1+'-'+productLinkage+'-'+product2]=1;
            }
        
            // update the linkage pairs affected by the new chain units   (first must identify the old linkage, must have tracked src and target.. )
            EdgeParticipants modifiedEdge = {src,tgt,-1};
            updateAdjacentLinks(cSS.links,cSS.chainPool[chainSelected], rxn, modifiedEdge,myOptions, newEdge);
        }

        // auto stopMonomerChainLinkingTime = Time::now();
        // timers["monomer-chainLinking"]+=stopMonomerChainLinkingTime-startMonomerChainLinkingTime;
    }
    // else, if two chains are chosen to react, execute the following code
    else if (whichChain1!=-1 && whichChain2!=-1) {
        // auto startChainChainLinkingTime = Time::now();
        Vertex src; // specifically for changxia's system, store source and target vertex indices within whichChain1
        Vertex tgt; // specifically for changxia's system, store source and target vertex indices within whichChain1
        // Check if the two chains are the same 
        bool isLoop = false;
        if (whichChain1==whichChain2) {
            // Note: loop formation is a unimolecular reaction. If loop formation were allowed here, the rate of loop formation would be dependent upon system size--in a larger system, the same would be less likely to serendipitously choose two of the same chain
            if (myOptions.loopingSeparateChannel==true) {
            //    std::cout << "Same Chain Chosen for Not Looping Rxn! Skip execution # " << cSS.executionNumber << "\n";
                cSS.mwv.sameChain=true; // don't want to update molecular weight if no reaction was executed!
                cSS.rejectedSteps++;
                return;
            }
            isLoop=true;
        }

        // For speed, copy into the chain with the smaller chain number, but for now just copy into whichChain1
        if (whichChain1!=whichChain2 || myOptions.loopingSeparateChannel==false) {
            if (whichChain1==whichChain2 && myOptions.trackFullStructureAsGraph) {
                // no one-membered loops
                if (whichVertex1==whichVertex2) {
                    // need to get out and not allow any loop to form
                    return;
                }
                // no two-membered loops
                if (boost::edge(whichVertex1,whichVertex2,cSS.chainPool[whichChain1].polymer).second || 
                boost::edge(whichVertex2,whichVertex1,cSS.chainPool[whichChain1].polymer).second) {
                    std::cout << "edge already exists" << std::endl;
                    return;
                }
                // no three-membered loops
                LoopParticipants pairToLoop;
                pairToLoop.src=whichVertex1;
                pairToLoop.tgt=whichVertex2;
                pairToLoop.chainIndex=whichChain1;
                if (isThreeMemberedRing(pairToLoop, cSS.chainPool[pairToLoop.chainIndex].polymer)) {
                    std::cout << " this would form a 3 membered ring -- 3 membered ring formation illegal" << std::endl;
                    return;
                    // prohibit 3 membered ring formation. this is particularly important when trying to have simultaneous breaking reactions, as the updateLinks doesn't work with such code...
                }
            }
            if (myOptions.trackLoopPairs) {
                // // update looping pairs (part 1 for whichChain1)
                for (int i=0;i<cSS.loopingRxns.size();++i) {
                    std::string fun1=cSS.loopingRxns[i].reactants[0];
                    std::string fun2=cSS.loopingRxns[i].reactants[1];
                    if (cSS.chainPool[whichChain1].funs.count(fun1) && cSS.chainPool[whichChain1].funs.count(fun2)) {
                        cSS.loopPairs[fun1][fun2]-=cSS.chainPool[whichChain1].funs[fun1]*cSS.chainPool[whichChain1].funs[fun2];
                        cSS.loopPairs[fun2][fun1]=cSS.loopPairs[fun1][fun2];
                    }
                }

                // update looping pairs (part 1 for whichChain2)
                for (int i=0;i<cSS.loopingRxns.size();++i) {
                    std::string fun1=cSS.loopingRxns[i].reactants[0];
                    std::string fun2=cSS.loopingRxns[i].reactants[1];
                    if (cSS.chainPool[whichChain2].funs.count(fun1) && cSS.chainPool[whichChain2].funs.count(fun2)) {
                        cSS.loopPairs[fun1][fun2]-=cSS.chainPool[whichChain2].funs[fun1]*cSS.chainPool[whichChain2].funs[fun2];
                        cSS.loopPairs[fun2][fun1]=cSS.loopPairs[fun1][fun2];
                    }
                }
            }
            Edge newEdge; // assign in the next section, but will need value again when updating edge counts in the linkage tracker
            if (myOptions.trackFullStructureAsGraph) {
                // first mark the two vertices that are to be linked
                cSS.chainPool[whichChain1].polymer[whichVertex1].chain1=true; // reactant 1
                cSS.chainPool[whichChain2].polymer[whichVertex2].chain2=true; // reactant 2
                
                // then copy the graph
                if (whichChain1!=whichChain2) {
                    // no need to copy graph if a loop is formed
                    copy_graph(cSS.chainPool[whichChain2].polymer,cSS.chainPool[whichChain1].polymer);
                }
                
                // then add the appropriate edge (this might be hard to do)

                // loop through until chain 1 and chain 2 are found
                auto vpair = boost::vertices(cSS.chainPool[whichChain1].polymer);
                for (auto iter = vpair.first; iter!=vpair.second; iter++) {
                    //boost::range_detail::integer_iterator<unsigned long> whichVertex;
                    if (cSS.chainPool[whichChain1].polymer[*iter].chain1==true) {
                        whichVertex1=*iter;
                    }
                    if (cSS.chainPool[whichChain1].polymer[*iter].chain2==true) {
                        whichVertex2=*iter;
                    }
                }
                src=whichVertex1; // store source vertex to update linkage counts later, specifically applies to Changxia's system
                tgt=whichVertex2; // store target vertex to update linkage counts later, specifically applies to Changxia's system
                newEdge = add_edge(whichVertex1,whichVertex2,cSS.chainPool[whichChain1].polymer).first;
                cSS.chainPool[whichChain1].nmp[newEdge]=productLinkage;
                
                cSS.chainPool[whichChain1].polymer[whichVertex2].chain2 =false; // reactant 2
                cSS.chainPool[whichChain1].polymer[whichVertex1].chain1=false; // reactant 1
                // update counters on monomers
                cSS.chainPool[whichChain1].polymer[whichVertex1].funs[reactant1]--; // reactant 1
                cSS.chainPool[whichChain1].polymer[whichVertex2].funs[reactant2]--; // reactant 2
                cSS.chainPool[whichChain1].polymer[whichVertex1].funs[rxn.products[0]]++; // product 1
                if (rxn.products.size()==2) {
                    cSS.chainPool[whichChain1].polymer[whichVertex2].funs[rxn.products[1]]++; // product 2
                }
            }
            // then update the counters in the new, combined chain. Funs available are the sum of the two minus the two that reacted... Note, if a new group formed, it should be added to the fun group list.
            if (isLoop==false) {
                for(auto it = cSS.chainPool[whichChain2].funs.begin(); it != cSS.chainPool[whichChain2].funs.end(); ++it) {
                    cSS.chainPool[whichChain1].funs[it->first] += it->second; 
                }
            }
            
            cSS.chainPool[whichChain1].funs[reactant1]--;
            cSS.chainPool[whichChain1].funs[reactant2]--;
            cSS.chainPool[whichChain1].funs[rxn.products[0]]++; 
            if (rxn.products.size()==2) {
                cSS.chainPool[whichChain1].funs[rxn.products[1]]++;
            }

            // update looping pairs (part 2)
            if (myOptions.trackLoopPairs) {
                if (isLoop) {
                    // tracking looping pairs only currently works when looping is considered as a separate reaction channel
                    throw std::runtime_error("this situation has not been handled");
                }
                for (int i=0;i<cSS.loopingRxns.size();++i) {
                    std::string fun1=cSS.loopingRxns[i].reactants[0];
                    std::string fun2=cSS.loopingRxns[i].reactants[1];
                    cSS.loopPairs[fun1][fun2]+=cSS.chainPool[whichChain1].funs[fun1]*cSS.chainPool[whichChain1].funs[fun2];
                    cSS.loopPairs[fun2][fun1]=cSS.loopPairs[fun1][fun2];
                }
            }

            
            // update number average and weight average molecular weight
            cSS.mwv.sameChain=isLoop;
            cSS.mwv.Mi_1=cSS.chainPool[whichChain1].weight;
            cSS.mwv.Mi_2=cSS.chainPool[whichChain2].weight;

            // update weight of the newly formed chain (sum of the two original chains)
            if (isLoop==false) {
                cSS.chainPool[whichChain1].weight+=cSS.chainPool[whichChain2].weight;
            }
            
            // update links tracker with the new number of linkages of each type
            if (myOptions.trackAllLinks && rxn.products.size()==2 && myOptions.trackFullStructureAsGraph) {
                std::string product1=rxn.products[0];
                std::string product2=rxn.products[1];

                // first, increment links tracker by new link produced 
                if (cSS.links.count(product1+'-'+productLinkage+'-'+product2) ) {
                    cSS.links[product1+'-'+productLinkage+'-'+product2]++; // increment 
                    cSS.chainPool[whichChain1].linkTotals[product1+'-'+productLinkage+'-'+product2]++; // increment 

                }
                else {
                    cSS.links[product1+'-'+productLinkage+'-'+product2]=1;
                    cSS.chainPool[whichChain1].linkTotals[product1+'-'+productLinkage+'-'+product2]=1;
                }
                for (auto quiter = cSS.chainPool[whichChain2].linkTotals.begin(); quiter!=cSS.chainPool[whichChain2].linkTotals.end(); ++quiter) {
                    cSS.chainPool[whichChain1].linkTotals[quiter->first] += quiter->second;
                }
                // update the linkage pairs affected by the new chain units   (first must identify the old linkage, must have tracked src and target.. )  

                EdgeParticipants modifiedEdge={src,tgt,-1}; // track the src and tgt vertices of the modified edge
                updateAdjacentLinks(cSS.links,cSS.chainPool[whichChain1], rxn, modifiedEdge,myOptions,newEdge);
            }

            // delete the old chain and property maps
            if (isLoop==false) {
                cSS.chainPool.erase(cSS.chainPool.begin()+whichChain2);
            }
            // update all the system-wide counters of moieties, etc to indicate that the reaction has taken place
            cSS.numMoieties[reactant1]--;
            cSS.numMoieties[reactant2]--;
            cSS.numMoieties[rxn.products[0]]++; 
            if (rxn.products.size()==2) {
                cSS.numMoieties[rxn.products[1]]++;
            }

        }
        // auto stopChainChainLinkingTime = Time::now();
        // timers["chain-chainLinking"]+=stopChainChainLinkingTime-startChainChainLinkingTime;
    }
    return;
}

void exchangeRxn(SystemVariables &cSS, Reaction rxn,double r3,double r4, ModelOptions myOptions,std::map<std::string,DoubleSec>& timers) {
    // this code executes an exchange reaction such as ab + C -> ac + B.
    // cSS = currentSystemState, i.e., the current values of all tracked quantities, polymers, etc.
    // note that ab + O -> ob + A would also be possible... the question is which group are we replacing, src or tgt?
    if (rxn.reactants.size()!=2 || rxn.products.size()!=2) {
        throw std::runtime_error("for an exchange reaction to work, there need to be two reactants and two products");
    }
    std::string exchangedGroup="tgt";
    std::string linkOriginal = rxn.reactants[0];
    std::string freeGrp = rxn.reactants[1];
    std::string linkProduct = rxn.products[0];
    std::string freeGrpProduct = rxn.products[1];
    // in an exchange reaction input file, let's always have the first reactant be the link and the second reactant be the free group
    // need to identify free grp on chain or on own
    // same with linkOriginal
    // check which chain linkOriginal is on... simple linear search through chainPool
    double whichLink=r3*cSS.numMoieties[linkOriginal];

    EdgeParticipants selectedEdge = linearSearchForMatchingEdge(cSS.chainPool,whichLink,linkOriginal,myOptions.trackFullStructureAsGraph);

    // next, identify the free group that will be participating in the exchange reaction
    double whichFreeGrp=r4*cSS.numMoieties[freeGrp];
    int findFreeGrp=0;
    std::string monomerSelected = ""; // if a monomer is selected, record which type of monomer has been selected
    // first, check if any monomers contain the groups
    if (cSS.funGroups.count( freeGrp )) {
        for (int i=0; i<cSS.funGroups[freeGrp].size();i++) {
            // for each monomer that contains the group, 
            // check the number of monomers remaining and the number of functional groups per monomer
            std::string monomerWithFunGroup=cSS.funGroups[freeGrp][i];
            int funGroupsOnMonomer=cSS.monomers[monomerWithFunGroup].funs[freeGrp];
            int funGroupsRemainingAttachedToMonomerI=cSS.numMonomers[monomerWithFunGroup]*funGroupsOnMonomer;
            findFreeGrp+=funGroupsRemainingAttachedToMonomerI;
            // check whether whichReactant has selected this monomer type to react
            if (findFreeGrp>whichFreeGrp) {
                // If chosen to react, store monomer type selected
                monomerSelected=monomerWithFunGroup; 
                // add goto statement here instead.
                break;
            }
        }
    }

    long findFreeGrpLong=findFreeGrp;
    if (findFreeGrpLong != monomerFunGrpCount(cSS.funGroups, freeGrp, cSS.numMonomers, cSS.monomers) ) {
        throw std::runtime_error("helpp");
    }
    // if the monomers didn't contain the free functional group, will need to search the chain_pool
    // simple linear search 
    int chainContainingFreeGroup=-1;
    boost::graph_traits < Graph >::vertex_descriptor vFreeGrp; // store vertex that contains the free group

    if (selectedEdge.src>10000) {
        throw std::runtime_error("this should nothappen");
    }

    if (monomerSelected=="") {
        for (int i=0; i<cSS.chainPool.size();i++) {
            // loop through each Graph object... and search for
            if (cSS.chainPool[i].funs.count(freeGrp)) {
                // interestingly , without the if /count statement, this led to erasing the stored "edge"
                findFreeGrp+=cSS.chainPool[i].funs[freeGrp];
            }

            // check whether this is the chain to react
            if (findFreeGrp>whichFreeGrp) {
                // if the moiety to react is found on this chain, store the chain number to react
                chainContainingFreeGroup=i;
                if (selectedEdge.src>10000) {
                    throw std::runtime_error("this should nothappen");
                }
                if (myOptions.trackFullStructureAsGraph) {
                    // loop through chain and find specific monomer to react...
                    findFreeGrp-=cSS.chainPool[i].funs[freeGrp];
                    
                    auto vpair = boost::vertices(cSS.chainPool[i].polymer);
                    bool findFreeGrpSuccessful=false;

                    for (auto iter = vpair.first; iter!=vpair.second; iter++) {
                        if (selectedEdge.src>10000) {
                            throw std::runtime_error("this should nothappen");
                        }
                        
                        // add the funs on the given monomer to the counter
                        if (cSS.chainPool[i].polymer[*iter].funs.count(freeGrp)) {
                            // interestingly , without the if /count statement, this led to erasing the "edge"
                            findFreeGrp+=cSS.chainPool[i].polymer[*iter].funs[freeGrp]; 
                        }
                        // check if this monomer is the selected one
                        if (findFreeGrp>whichFreeGrp) {
                            // found node to react, store its location!
                            vFreeGrp=*iter;
                            findFreeGrpSuccessful=true;
                            break;
                        }
                    }
                    if (!findFreeGrpSuccessful) {
                        throw std::runtime_error("failed to find vFreeGrp");
                        // if this happens, it means that there are a different number of freeGrps in the chain+monomer counters than in numMoieties[freeGrp] or in chain vs monomer counter (actually only the latter)
                    }
                }
                break;
            }
            
        }
        if (findFreeGrp<=whichFreeGrp) {
            throw std::runtime_error("free grp not found!!");
        }
    }    

    // now we have identified the groups to react and we must carry out the actual reaction
    if (selectedEdge.src>10000 || (vFreeGrp>1000 && monomerSelected=="") ) {
        throw std::runtime_error("this should nothappen");    
    }

    // // two new polymers formed (well possibly a monomer instead of a polymer, but let's just store it as a polymer...

    int whichChain=selectedEdge.chainIndex; // index of chain to bust, just for ease of referencing here

    // if monomer free group + link to exchange
    if (monomerSelected!="") {
        boost::remove_edge(selectedEdge.src,selectedEdge.tgt,cSS.chainPool[whichChain].polymer);     // bust the edge

        // add the new monomer to the chain
        vFreeGrp = add_vertex(cSS.chainPool[whichChain].polymer);
        cSS.chainPool[whichChain].polymer[vFreeGrp]=cSS.monomers[monomerSelected];

        // add functional groups from new monomer
        auto iter=cSS.chainPool[whichChain].polymer[vFreeGrp].funs.begin();
        while (iter != cSS.chainPool[whichChain].polymer[vFreeGrp].funs.end() ) {
            cSS.chainPool[whichChain].funs[iter->first]+=iter->second;
            iter++;
        }
        // record old molecular weight of polymer
        cSS.mwv.Mi_2_old=cSS.chainPool[whichChain].weight;
        cSS.mwv.Mi_1_old=cSS.monomers[monomerSelected].weight; // could just as easily be written chainPool[whichChain].polymer[vFreeGrp].weight
        // update molecular weight of polymer
        cSS.chainPool[whichChain].weight+=cSS.chainPool[whichChain].polymer[vFreeGrp].weight;

    }
    else {
        // else free group and link are both on polymers
        // the two polymers make two new polymers
        if (whichChain==chainContainingFreeGroup) {
           // std::cout << "disallow intrachain exchange rxns *for now* " << std::endl;
            cSS.mwv.sameChain=true; // molecular weight or whatever should not be affected--we are rejecting this step.
            return; 
        }
        boost::remove_edge(selectedEdge.src,selectedEdge.tgt,cSS.chainPool[whichChain].polymer); // bust the edge

        // copy chainContainingFreeGroup into whichChain
        // first need to mark "vFreeGrp" to update... let's use the chain1 bool... maybe i should call them "marker1" and "marker 2"
        cSS.chainPool[chainContainingFreeGroup].polymer[vFreeGrp].chain1=true;
        int whichChainSize= boost::num_vertices( cSS.chainPool[whichChain].polymer );
        copy_graph(cSS.chainPool[chainContainingFreeGroup].polymer,cSS.chainPool[whichChain].polymer);
        cSS.chainPool[chainContainingFreeGroup].polymer[vFreeGrp].chain1=false; // mark as false in the original chain

        // update vFreeGrp
        if (cSS.chainPool[whichChain].polymer[whichChainSize+vFreeGrp].chain1) {
            vFreeGrp=whichChainSize+vFreeGrp;
            cSS.chainPool[whichChain].polymer[vFreeGrp].chain1=false; // reset to false in the copied chain
        }
        else {
            // loop through till the new vFreeGrp is found
            // again a linear search through the polymer to find...
            auto vpair = boost::vertices(cSS.chainPool[whichChain].polymer);
            bool updateSuccessful=false;
            for (auto iter = vpair.first; iter!=vpair.second; iter++) {
                if ( cSS.chainPool[whichChain].polymer[*iter].chain1) {
                    vFreeGrp=*iter;
                    updateSuccessful=true;
                    break;
                }
            }
            cSS.chainPool[whichChain].polymer[vFreeGrp].chain1=false; // reset to false in the copied chain
            if (!updateSuccessful) {
                throw std::runtime_error("vFreeGrp was not updated successfully after a chain+chain exchange rxn");
            }
        }


        // add functional groups from chainContainingFreePolymer
        auto iter=cSS.chainPool[chainContainingFreeGroup].funs.begin();
        while (iter != cSS.chainPool[chainContainingFreeGroup].funs.end()) {
            cSS.chainPool[whichChain].funs[iter->first]+=iter->second;
            iter++;
        }

        // record old molecular weight of polymer
        cSS.mwv.Mi_2_old=cSS.chainPool[whichChain].weight;
        cSS.mwv.Mi_1_old=cSS.chainPool[chainContainingFreeGroup].weight;
        // update molecular weight of polymer
        cSS.chainPool[whichChain].weight+=cSS.chainPool[chainContainingFreeGroup].weight;

    }
    // add the new edge
    Edge newEdge = add_edge(selectedEdge.src,vFreeGrp,cSS.chainPool[whichChain].polymer).first;
    cSS.chainPool[whichChain].nmp[newEdge]=linkProduct;

    // update functional group counts on monomer
    (cSS.chainPool[whichChain].polymer[selectedEdge.tgt].funs[freeGrpProduct])++;
    (cSS.chainPool[whichChain].polymer[selectedEdge.src].funs[linkProduct])++;
    (cSS.chainPool[whichChain].polymer[selectedEdge.src].funs[linkOriginal])--;
    (cSS.chainPool[whichChain].polymer[vFreeGrp].funs[freeGrp])--;
    if (cSS.chainPool[whichChain].polymer[selectedEdge.src].funs[linkOriginal] <0 || cSS.chainPool[whichChain].polymer[vFreeGrp].funs[freeGrp] <0) {
        throw std::runtime_error("this shouldn't happen");
    }
    // update functional group counts on polymer
    (cSS.chainPool[whichChain].funs[freeGrpProduct])++;
    (cSS.chainPool[whichChain].funs[linkProduct])++;
    (cSS.chainPool[whichChain].funs[linkOriginal])--;
    (cSS.chainPool[whichChain].funs[freeGrp])--;
    if (cSS.chainPool[whichChain].funs[linkOriginal]<0) {
        throw std::runtime_error("this shouldn't happen");
    }
    
    // check if tgt and src are still connected (possible if loops form) 
    Connectivity myConnectivity=checkConnectivity(cSS.chainPool[whichChain].polymer);
    if (myConnectivity.isSingleComponent) {
        // only need to update fun group counts and molecular weight trackers
        throw std::runtime_error("for now , simultanous looping and exchange reactions are not allowed! functionality to be added at later date");
    }
    else {
        // break chains apart
        Chain newFragment=cSS.chainPool[whichChain];
        // create undirected graph to check for connectivity
        UndirectedGraph uG;
        copy_graph(cSS.chainPool[whichChain].polymer, uG); // there may be a more efficient way of creating an undirected graph from a directed graph

        // now delete all the ones *not connected* to component 0 from the *old* graph
        // delete all the ones *connected* to component 0 from the *new* graph
        Graph d=cSS.chainPool[whichChain].polymer;
        int countDeletedEdges=0; // count deleted edges (corresponding to component 1) on old graph (component 0)
        int countDeletedEdgesNewGraph=0; // count deleted edges (corresponding to component 0) on new graph (component 1)
        int totalVertices=boost::num_vertices(cSS.chainPool[whichChain].polymer); // start and end vertices
        for (int i = 0; i < totalVertices; ++i) {

            if (myConnectivity.component[i] != 0) {
                int index = i-countDeletedEdges; // calculate index in the old chain
                // update chain's functional group count at the end...
                boost::clear_vertex(i-countDeletedEdges,cSS.chainPool[whichChain].polymer);
                boost::remove_vertex(i-countDeletedEdges,cSS.chainPool[whichChain].polymer);
                countDeletedEdges++;
            }
            else {
                int index=i-countDeletedEdgesNewGraph; // calculate index of the new chain
                // calculate the effect of deleting the particular monomer
                for (auto itfun = d[index].funs.begin(); itfun != d[index].funs.end(); itfun++ ) {
                    newFragment.funs[itfun->first]-=itfun->second;
                } 
                // update chain weight
                newFragment.weight-=d[index].weight;
                boost::clear_vertex(i-countDeletedEdgesNewGraph,d); // delete all edges associated with the vertex
                boost::remove_vertex(i-countDeletedEdgesNewGraph,d);
                countDeletedEdgesNewGraph++;
            }
        }
        
        // fix the fun group count on the old chain
        newFragment.polymer=d;
        for (auto it=cSS.chainPool[whichChain].funs.begin(); it!=cSS.chainPool[whichChain].funs.end();it++) {
            // funs on the old fragment are equal to the original number of funs minus the ones on the new fragment..
            it->second=it->second-newFragment.funs[it->first];
        }
        // fix the weight on the old chain
        (cSS.chainPool[whichChain].weight)-=newFragment.weight;
        
        // put newFragment in the chainPool
        if (monomerSelected!="") {
            cSS.chainPool.push_back(newFragment); // if there was only 1 chain reacting at the beginning, push back new fragment
            // update namemap
            cSS.chainPool.back().nmp=get(&EdgeData::name,newFragment.polymer); // following example https://stackoverflow.com/questions/3527313/using-bundled-properties-as-the-weight-map-in-dijkstra-shortest-paths
        }
        else {
            cSS.chainPool[chainContainingFreeGroup]=newFragment; // if two chains were reacted, place "new fragment" in the chainContainingFreeGroup
        }
        
        // update the molecular weight
        cSS.mwv.Mi_1=newFragment.weight;
        cSS.mwv.Mi_2=cSS.chainPool[whichChain].weight;

        // update numMonomers
        if (monomerSelected!="") {
            cSS.numMonomers[monomerSelected]--;
        }
            
    } 

    // update number of moieties available
    cSS.numMoieties[linkProduct]++;
    cSS.numMoieties[freeGrpProduct]++;
    cSS.numMoieties[freeGrp]--;
    cSS.numMoieties[linkOriginal]--;
    return;
}

LoopParticipants linearSearchForLoopingPair(SystemVariables &cSS,double randomPair, Reaction rxn, ModelOptions myOptions) {

    
    // identify the chain to loop
    LoopParticipants selectedPair;
    std::string reactant1 = rxn.reactants[0];
    std::string reactant2 = rxn.reactants[1];

    signed long long findPairToLoop=0;
    signed long long chainToLoop=-1; 
    for (int i=0; i<cSS.chainPool.size();i++) {
            // loop through each Graph object... and search for
            findPairToLoop+=cSS.chainPool[i].funs[reactant2]*cSS.chainPool[i].funs[reactant1];

            // check whether this is the chain to react
            if (findPairToLoop>randomPair) {
                // if the moiety to react is found on this chain, store the chain number to react
                selectedPair.chainIndex=i;
                break;
            }
    }

    // if we're tracking the full structure, also need to return the specific vertices to react
    if (myOptions.trackFullStructureAsGraph ) {
        findPairToLoop-=cSS.chainPool[selectedPair.chainIndex].funs[reactant2]*
        cSS.chainPool[selectedPair.chainIndex].funs[reactant1]; // decrement by number of looping pairs on chain
        bool foundPair=false;
        // a big question here will be whether there is any restriction on minimum loop size... i.e. can two functional groups on the same monomer react? 
        // I think maybe we should disallow that, but the problem is that the counting of pairs on a chain as funs1*funs2 doesn't work if one monomer has functional groups of both types
        // might need to store a matrix on each chain of "allowedLoopPairs"... 
        Graph polymer = cSS.chainPool[selectedPair.chainIndex].polymer;
        
        for (auto vxs=boost::vertices(polymer); vxs.first!=vxs.second; ++vxs.first) {
            int numFuns = polymer[*vxs.first].funs[reactant1];
            if (numFuns!=0) {
                
                for (auto vxs2=boost::vertices(polymer); vxs2.first!=vxs2.second; ++vxs2.first) {
                    // find all pairs
                    int numFuns2=0;
                    if ( polymer[*vxs2.first].funs.count(reactant2) ) {
                        numFuns2=polymer[*vxs2.first].funs[reactant2];
                    }
                    findPairToLoop+=numFuns*numFuns2;
                    if (findPairToLoop > randomPair) {
                        selectedPair.src = *vxs.first;
                        selectedPair.tgt = *vxs2.first;
                        foundPair=true;
                        break;
                    }
                }
            }
            if (foundPair) {
                break;
            }
        }
        if (foundPair==false) {
            throw std::runtime_error("did not correctly identify looping pair... ");
        }
    }
    return selectedPair;
}

bool isThreeMemberedRing(LoopParticipants pairToLoop,Graph polymer) {
    // this function takes as input two vertices that an edge could be placed between
    // as well as the polymer that contains these vertices. This function then returns
    // a bool specifying whether a 3 membered ring would be formed by adding this edge. 
    // note that this will also return true if a 2 membered ring would be formed.

    bool isThreeMembered=false;

    // loop edge hasn't yet been added yet

    // start with out edges from src
    auto out_edges_src = boost::out_edges(pairToLoop.src,polymer);
    while (out_edges_src.first!=out_edges_src.second) {
        // for each target of these out edges, see if they have either in or out edges connecting them to pairToLoop.tgt
        Vertex srctgt = target(*out_edges_src.first, polymer);
        
        // check out edges
        auto otedgs = boost::out_edges(srctgt, polymer); 
        while (otedgs.first!=otedgs.second) {
            Vertex srctgttgt = target(*otedgs.first,polymer);
            if (srctgttgt == pairToLoop.tgt) {
                return true;
            }
            otedgs.first++;
        }

        // check in edges
        auto inedgs = boost::in_edges(srctgt,polymer);
        while (inedgs.first!=inedgs.second ) {
            Vertex srctgtsrc = source(*inedgs.first,polymer);
            if (srctgtsrc == pairToLoop.tgt) {
                return true;
            }
            inedgs.first++;
        }
        out_edges_src.first++;
    }

    // then check in edges from src
    auto in_edges_src = boost::in_edges(pairToLoop.src,polymer);
    while (in_edges_src.first!=in_edges_src.second) {
        Vertex srcsrc = source(*in_edges_src.first,polymer);

        // check out edges
        auto otedgs = boost::out_edges(srcsrc, polymer); 
        while (otedgs.first!=otedgs.second) {
            Vertex srcsrctgt = target(*otedgs.first,polymer);
            if (srcsrctgt == pairToLoop.tgt) {
                return true;
            }
            otedgs.first++;
        }

        // check in edges
        auto inedgs = boost::in_edges(srcsrc,polymer);
        while (inedgs.first!=inedgs.second ) {
            Vertex srcsrcsrc = source(*inedgs.first,polymer);
            if (srcsrcsrc == pairToLoop.tgt) {
                return true;
            }
            inedgs.first++;
        }

        in_edges_src.first++;
    }
    return false;
}

// double radiusOfGyration(int DP, double weight) {
//     // calculate radius of gyration of a polymer 

//     // From Hiemenz Polymer Chemistry page 232: 
//     // "A useful rule of thumb for polymers is that Rg is about 100 A when 
//     // M = 1E5 g/mol. This number can be used to estimate Rg for any other M 
//     // by recalling the proportionality of Rg to M^1/2
//     return sqrt( (weight * DP) / 1.0E5 ) * 100; // in Angstroms
// }

double meanEndToEndDistance(long topologicalDistance, double persistenceLength) {
    // Hiemenz Polymer Chemistry 2nd edition Equation 6.4.7 (page 227)
    // <h^2>0 = L*lk, where L is the contour length and lk is the kuhn length
    double lk = persistenceLength/2.0; // kuhn length
    double l = 5.52851; // edge length in Angstroms; i.e. the length of a monomer
    // see file "calculatingMonomerLengthsChangxia.gjf" and notes on Feb 3.
    double L = topologicalDistance * l; // contour length
    return sqrt(L*lk);
}
double polymerChainVolumeInLiters(PossibleLoopingRxn candidateLoop, int dmax) {
    // calculate the volume of a polymer chain in liters based upon topological distance
    //double Rg = radiusOfGyration( boost::num_vertices(polymer) , monomerWeight); // in Angstroms (Hiemenz Polymer Chemistry rule of thumb chapter 6)
    
    double persistenceLength = 100; // in Angstroms. 
    // Persistence length taken from entry for Poly(p-phenylene terephthalamide) in Methane sulfonic acid solvent
    // in Table 6.2 from Hiemenz Polymer Chemistry 2nd edition.

    if (candidateLoop.d > dmax) {
        return 0;
    }
    else {
        double radius = meanEndToEndDistance(candidateLoop.d, persistenceLength);
        return (4.0/3.0) * 3.1415 * pow(radius,3) / 1e+33; // volume of a sphere, 1e33 to convert between cubic angstroms and L 
}
    }
    
// LoopParticipants SearchForLoopingPair(SystemVariables &cSS,double whichPair, Reaction rxn, ModelOptions myOptions) {

//     // this function identifies the pair to loop

//     std::string reactant1 = rxn.reactants[0];
//     std::string reactant2 = rxn.reactants[1];

//     double cumulativeProbability = 0;
//     double Rg;
//     double Vp;
//     LoopParticipants theseParticipants;
//     theseParticipants.chainIndex=-1; // a negative entry here allows determination of whether a candidate loop was properly located

//     for (int p = 0; p<cSS.chainPool.size(); ++p) {
//         // want to exclude all loops size 3 and smaller
//         int DP_p = boost::num_vertices(cSS.chainPool[p].polymer);
//         double Vp;
//         if (DP_p>3) {
//             signed long pairs = cSS.chainPool[p].funs[reactant1]*cSS.chainPool[p].funs[reactant2];
//             if (pairs!=0) {
//                 // calculate Rg for the polymer
//                 //Rg = sqrt(1.0E5 / ( cSS.monomers[0].weight * DP_p) ) * 100; // in Angstroms
//                 // Vp = (4.0/3.0) * 3.1415 * pow(Rg,3) / 1e+33; // 1e33 to convert between cubic angstroms and L 
//                 Vp = polymerChainVolumeInLiters(cSS.chainPool[p].polymer, cSS.monomers["cEp"].weight); 
//                 cumulativeProbability += pairs / Vp;
//             }
//             if (cumulativeProbability>whichPair) {
//                 theseParticipants.chainIndex = p;
//                 cumulativeProbability -= pairs / Vp;
//             }
//         }
//         // now go through and figure out which specific monomers should be looped...
//         if (theseParticipants.chainIndex == p) {

//             // we might simplify this since we KNOW that there is only 1 cE...

//             // the cE is always the zeroth reactant. 
//             if (reactant1!="cE") {
//                 throw std::runtime_error( " looping has only been developed for Changxia Shi's hyperbranched polyester system!");
//             }
//             // find the monomers
//             bool foundtgt = false;
//             bool foundsrc = false;
//             for (auto vxs = vertices(cSS.chainPool[p].polymer); vxs.first!=vxs.second; ++vxs.first) {
                
//                 std::string chainUnitName = getChainUnitName(cSS.chainPool[p].polymer[*vxs.first]);

//                 if (foundtgt==false) {
//                     if ( foundtgt == false && chainUnitName == reactant2 ) {
//                         cumulativeProbability += 1.0 / Vp;
//                     }

//                     if ( cumulativeProbability > whichPair ) {
//                         theseParticipants.tgt = *vxs.first;
//                         foundtgt=true;
//                     }

//                 }
                
//                 if ( chainUnitName == "cE" ) {
//                     theseParticipants.src = *vxs.first;
//                     foundsrc=true;
//                 }
//                 if (foundtgt && foundsrc) {
//                     break;
//                 }
//             }
//             if (foundtgt==false || foundsrc == false) {
//                 throw std::runtime_error("unable to locate vertices correctly..."); 
//             }

//             break;
//         }
//     }

    
//     return theseParticipants;
// }

void loopingRxn(SystemVariables &cSS, PossibleLoopingRxn lr, ModelOptions myOptions, std::vector<Reaction> reactionList) {
    // right now, do not allow loops > dmax. This should change in the future, but for a first pass, just allow loops of size minLoopSize to dmax...
    // in this case, the vertices have ALREADY been identified, so we just need to connect them.

    int reactionNo = -1;

    if (lr.catAlcoholName=="Em*p_p") {
        // this is the looping analog of R6, 2f
        // Note: this is reaction 71
        reactionNo = 71;

    }
    if (lr.catAlcoholName=="Em*p_m") {
        // this is the looping analog of R8, 2f
        // Note: this is reaction 73
        reactionNo = 73;
    }
    if (lr.catAlcoholName=="E1*") {
        
        // this is the looping analog of R5, 2f
        // Note: this is reaction 70
        reactionNo = 70;    
    }
    if (lr.catAlcoholName=="E3*") {
        // this is the looping analog of R7, 2f
        // Note: this is reaction 72
        reactionNo = 72;
    }

    if (reactionNo==-1) {
        throw std::runtime_error("unexpected catAlcoholName encountered...");
    }
    if (myOptions.trackFullStructureAsGraph == false || myOptions.trackAllLinks == false) { 
        throw std::runtime_error("looping reactions are only implemented when full structure and all links are tracked...");
    }
    // could also throw an error if we track links but our linking + looping reactions only have one product...

    // check to make sure the list of reactants and products are long enough. Same for th reactionList
    if (reactionList[reactionNo].reactants.size()!=2 || reactionList[reactionNo].products.size()!=2) {
        throw std::runtime_error("looping reactions are only implemented when there are exactly 2 reactants and 2 products...");
    }
    if (reactionList.size()<reactionNo+1) {
        throw std::runtime_error("reactionList is not long enough to contain the looping reaction...");
    }

    // get the reactants and products
    std::string reactant1 = reactionList[reactionNo].reactants[0];
    std::string reactant2 = reactionList[reactionNo].reactants[1];
    std::string product2 = reactionList[reactionNo].products[1];
    std::string product1 = reactionList[reactionNo].products[0];

    // add the edge
    Edge e = boost::add_edge(lr.vertices.first,lr.vertices.second,cSS.chainPool[lr.p].polymer).first;;

    // set edge name
    cSS.chainPool[lr.p].polymer[e].name=reactionList[reactionNo].linkFormed;

    // could run some tests here to make sure loop is not too small. However, topological distance should have already been determined in the function that identified the candidate pair,
    // and the pair should have been rejected if the loop was too small.

    // if (boost::edge(lr.vertices.first,lr.vertices.second,cSS.chainPool[lr.p].polymer).second || 
    //     boost::edge(lr.vertices.second,lr.vertices.first,cSS.chainPool[lr.p].polymer).second) {
    //     // the proposed edge is a two membered loop!
    //     //std::cout << "edge arleady exists" << std::endl;
    //     cSS.catBindingTally=0; // allow a few more catalyst binding/unbinding steps to take place such that the catalyst chooses to sit somewhere that the formation of a two membered loop is less likely
    //     return;
    // }
    // else if (isThreeMemberedRing(pairToLoop, cSS.chainPool[pairToLoop.chainIndex].polymer)) {
    //     //
    //     //std::cout << " this would form a 3 membered ring -- 3 membered ring formation illegal" << std::endl;
    //     cSS.catBindingTally=0; // allow a few more catalyst binding/unbinding steps to take place such that the catalyst chooses to sit somewhere that the formation of a three membered loop is less likely
    //     return;
    //     // prohibit 3 membered ring formation. this is particularly important when trying to have simultaneous breaking reactions, as the updateLinks doesn't work with such code...
    // }

    //cSS.chainPool[lr.p].nmp[e]=reactionList[reactionNo].linkFormed; // update name of new linkage

    // update overall trackers
    cSS.numMoieties[reactant1]--;
    // cSS.numMoieties["*"]--; // catalyst
    // throw std::runtime_error("EDIT THE ABOVE LINE");
    cSS.numMoieties[reactant2]--;
    cSS.numMoieties[product1]++;
    cSS.numMoieties[product2]++;

    // update chain specific trackers
    cSS.chainPool[lr.p].funs[reactant1]--;
    cSS.chainPool[lr.p].funs[reactant2]--;
    cSS.chainPool[lr.p].funs[product1]++;
    cSS.chainPool[lr.p].funs[product2]++;

    // update specific chain unit trackers 
    cSS.chainPool[lr.p].polymer[lr.vertices.first].funs[reactant1]--;
    cSS.chainPool[lr.p].polymer[lr.vertices.first].funs[product1]++;
    cSS.chainPool[lr.p].polymer[lr.vertices.second].funs[reactant2]--;
    cSS.chainPool[lr.p].polymer[lr.vertices.second].funs[product2]++;
        
    // update links trackers
    // first, increment links tracker by new link produced 
    if (cSS.links.count(product1+'-'+reactionList[reactionNo].linkFormed+'-'+product2) ) {
        cSS.links[product1+'-'+reactionList[reactionNo].linkFormed+'-'+product2]++; // increment 
    }
    else {
        cSS.links[product1+'-'+reactionList[reactionNo].linkFormed+'-'+product2]=1;
    }
    // also, increment chain fun groups tracker by new link produced 
    if (cSS.chainPool[lr.p].linkTotals.count(product1+'-'+reactionList[reactionNo].linkFormed+'-'+product2) ) {
        cSS.chainPool[lr.p].linkTotals[product1+'-'+reactionList[reactionNo].linkFormed+'-'+product2]++; // increment 
    }
    else {
        cSS.chainPool[lr.p].linkTotals[product1+'-'+reactionList[reactionNo].linkFormed+'-'+product2]=1;
    }

    EdgeParticipants modifiedEdge;
    modifiedEdge.src = source(e,cSS.chainPool[lr.p].polymer);
    modifiedEdge.tgt = target(e,cSS.chainPool[lr.p].polymer);
    modifiedEdge.chainIndex=-1;
    updateAdjacentLinks(cSS.links,cSS.chainPool[lr.p], reactionList[reactionNo], modifiedEdge,myOptions,e);
    
    cSS.mwv.sameChain=true; // note that a loop was formed when calling the molecular weight function
    
}
// void loopingRxn(SystemVariables &cSS, Reaction rxn,double r3,double r4, ModelOptions myOptions) {
//     // set reactants and products
//     std::string reactant1 = rxn.reactants[0];
//     std::string reactant2 = rxn.reactants[1];
//     std::string product1 = rxn.products[0];
//     std::string product2;
//     if (rxn.products.size()==2) {
//         product2 = rxn.products[1];
//     }
    
//     double whichPair;
//     // if (!myOptions.trackLoopPairs) {
//     //     throw std::runtime_error("for a looping reaction, need to know number of candidate looping pairs... but they aren't being tracked");
//     // }
//     throw std::runtime_error("need to fix this before running!");
//     whichPair = 0; // r3*rxn.currentTotalLoopingProbability; // could divide out by the whole k thing so that we're just left with pairs / Vp.. yea i did this now

//     LoopParticipants pairToLoop = SearchForLoopingPair(cSS,whichPair,rxn,myOptions);
    
//     // got to add the edge and increment 
//     Edge e;
//     if (myOptions.trackFullStructureAsGraph) {

//         if (boost::edge(pairToLoop.src,pairToLoop.tgt,cSS.chainPool[pairToLoop.chainIndex].polymer).second || 
//         boost::edge(pairToLoop.tgt,pairToLoop.src,cSS.chainPool[pairToLoop.chainIndex].polymer).second) {
//             // the proposed edge is a two membered loop!
//             //std::cout << "edge arleady exists" << std::endl;
//             cSS.catBindingTally=0; // allow a few more catalyst binding/unbinding steps to take place such that the catalyst chooses to sit somewhere that the formation of a two membered loop is less likely
//             return;
//         }
//         else if (isThreeMemberedRing(pairToLoop, cSS.chainPool[pairToLoop.chainIndex].polymer)) {
//             //
//             //std::cout << " this would form a 3 membered ring -- 3 membered ring formation illegal" << std::endl;
//             cSS.catBindingTally=0; // allow a few more catalyst binding/unbinding steps to take place such that the catalyst chooses to sit somewhere that the formation of a three membered loop is less likely
//             return;
//             // prohibit 3 membered ring formation. this is particularly important when trying to have simultaneous breaking reactions, as the updateLinks doesn't work with such code...
//         }
//         else {
//             //std::cout << " finally, allowed ring here :) " << std::endl;
//             e = boost::add_edge(pairToLoop.src,pairToLoop.tgt,cSS.chainPool[pairToLoop.chainIndex].polymer).first;
//         }
//         cSS.chainPool[pairToLoop.chainIndex].nmp[e]=rxn.linkFormed; // update name of new linkage
//     }

//     // update overall trackers
//     cSS.numMoieties[reactant1]--;
//     // cSS.numMoieties["*"]--; // catalyst
//     // throw std::runtime_error("EDIT THE ABOVE LINE");
//     cSS.numMoieties[reactant2]--;
//     cSS.numMoieties[product1]++;

//     // update chain specific trackers
//     cSS.chainPool[pairToLoop.chainIndex].funs[reactant1]--;
//     cSS.chainPool[pairToLoop.chainIndex].funs[reactant2]--;
//     cSS.chainPool[pairToLoop.chainIndex].funs[product1]++;

//     // update specific chain unit trackers if the full topology is being tracked
//     if (myOptions.trackFullStructureAsGraph) {
//         cSS.chainPool[pairToLoop.chainIndex].polymer[pairToLoop.src].funs[reactant1]--;
//         cSS.chainPool[pairToLoop.chainIndex].polymer[pairToLoop.src].funs[product1]++;
//         cSS.chainPool[pairToLoop.chainIndex].polymer[pairToLoop.tgt].funs[reactant2]--;
//     }
    
//     if (rxn.products.size()==2) {
//         cSS.numMoieties[product2]++;
//         cSS.chainPool[pairToLoop.chainIndex].funs[product2]++;
//         if (myOptions.trackFullStructureAsGraph) {
//             cSS.chainPool[pairToLoop.chainIndex].polymer[pairToLoop.tgt].funs[product2]++;
//         }
        
//     }

//     bool updatedCurrentReactionAlready = false;

//     // need to update number of looping pairs... but now all on that chain are affected... can split this into another function
//     // loop through every looping pair in the chain.. but this is strange bc i don't know which pairs are looping pairs... use reaction list...
//     // loop through reaction list and check for looping reactions
//     for (int i=0; i<cSS.loopingRxns.size();++i) {  
//         std::string fun1=cSS.loopingRxns[i].reactants[0];
//         std::string fun2=cSS.loopingRxns[i].reactants[1];
//         // check if either of these was a participant in the above reaction
//         if (fun1==reactant1 || fun1 ==reactant2 || fun2 ==reactant1 || fun2 == reactant2 ) {
//             if (fun1==reactant1 && fun2 == reactant2) {
//                 if (updatedCurrentReactionAlready) {
//                     throw std::runtime_error("something went wrong"); // note: if two looping reactions involving 
//                     // the same two functional groups are present, throw an error. the code cannot yet handle this
//                 }
//                 // update looping pairs for this specific reaction
//                 if (fun1!=fun2) {
//                     cSS.loopPairs[reactant1][reactant2]+=cSS.chainPool[pairToLoop.chainIndex].funs[reactant2]*
//                     cSS.chainPool[pairToLoop.chainIndex].funs[reactant1]-(cSS.chainPool[pairToLoop.chainIndex].funs[reactant2]+1)*
//                     (cSS.chainPool[pairToLoop.chainIndex].funs[reactant1]+1);
//                     cSS.loopPairs[reactant2][reactant1]=cSS.loopPairs[reactant1][reactant2];
//                 }
//                 if (fun1==fun2) {
//                     // things are more interesting in this case. see derivation.
//                     throw std::runtime_error("I haven't gotten around to testing this scenario yet");
//                     signed long long n = cSS.chainPool[pairToLoop.chainIndex].funs[reactant1];

//                     cSS.loopPairs[reactant1][reactant2] -= n;
                    
//                 }
//                 updatedCurrentReactionAlready=true;
//             }
//             else {
//                 // suppose just one of these is involved
//                 std::string involvedGroup;
//                 std::string uninvolvedGroup; 
//                 if (fun1==reactant1 || fun1 == reactant2) {
//                     involvedGroup=fun1;
//                     uninvolvedGroup=fun2;
                    
//                 } 
//                 else if (fun2==reactant1 || fun2 == reactant2) {
//                     involvedGroup = fun2;
//                     uninvolvedGroup=fun1;
//                 }
//                 else {
//                     throw std::runtime_error(" this should never happen"); // one of the first two statements should be true;
//                 }
//                 // go ahead and update the number of candidate looppairs...
//                 cSS.loopPairs[fun1][fun2]-=cSS.chainPool[pairToLoop.chainIndex].funs[uninvolvedGroup];
//             }
            
//         }
//     }
//     if (myOptions.trackFullStructureAsGraph && myOptions.trackAllLinks) {
//         // update links. I think links should only be updated if a bool says we realy want to track them
//         // also, should throw an error if we track links but our linking + looping reactions only have one product...
//         // first, increment links tracker by new link produced 
//         if (cSS.links.count(product1+'-'+rxn.linkFormed+'-'+product2) ) {
//             cSS.links[product1+'-'+rxn.linkFormed+'-'+product2]++; // increment 
//         }
//         else {
//             cSS.links[product1+'-'+rxn.linkFormed+'-'+product2]=1;
//         }
//         // also, increment chain fun groups tracker by new link produced 
//         if (cSS.chainPool[pairToLoop.chainIndex].linkTotals.count(product1+'-'+rxn.linkFormed+'-'+product2) ) {
//             cSS.chainPool[pairToLoop.chainIndex].linkTotals[product1+'-'+rxn.linkFormed+'-'+product2]++; // increment 
//         }
//         else {
//             cSS.chainPool[pairToLoop.chainIndex].linkTotals[product1+'-'+rxn.linkFormed+'-'+product2]=1;
//         }

//         EdgeParticipants modifiedEdge;
//         modifiedEdge.src = source(e,cSS.chainPool[pairToLoop.chainIndex].polymer);
//         modifiedEdge.tgt = target(e,cSS.chainPool[pairToLoop.chainIndex].polymer);
//         modifiedEdge.chainIndex=-1;
//         updateAdjacentLinks(cSS.links,cSS.chainPool[pairToLoop.chainIndex], rxn, modifiedEdge,myOptions,e);
//     }
    
//     cSS.mwv.sameChain=true; // note that a loop was formed when calling the molecular weight function
    
//     // update reaction type trackers
//     cSS.catBindingTally=0;
    
//     return;
// }


void allophanateRxn(SystemVariables &cSS, Reaction rxn,double r3,double r4, ModelOptions myOptions,std::map<std::string,DoubleSec>& timers) { 
    // first need to identify the linkage... it should also be stored as a functional group, so I should be able to search through functional groups
    // linkage

    // reactant 0: urethane (or other link) = linkToBust
    // reactant 1: isocyanate (or other free group)
    // product 0: allophanate (or other new link formed) = linkFormed
    // new link added FROM src of urethane TO reactant1 (free grp)

    std::string oldLink = rxn.linkToBust;
    std::string newLink = rxn.linkFormed;
    std::string freeGrp = rxn.reactants[1];

    EdgeParticipants oldLinkDeets = linearSearchForMatchingEdge(cSS.chainPool,r3*cSS.numMoieties[oldLink], oldLink,myOptions.trackFullStructureAsGraph);
    GroupLocation freeGroupLoc = linearSearchForParticipatingGroup(cSS,r4, freeGrp, myOptions.trackFullStructureAsGraph);

    // if the two reactions are serendipitously on the same chain, reject for now
    if (oldLinkDeets.chainIndex==freeGroupLoc.whichChain) {
        //std::cout << "for now, rejecting a looping allophanate-like reaction" << std::endl;
        return;
    }
    if (myOptions.trackAllLinks) {
        throw std::runtime_error("Don't yet have ability to do an allophanateRxn and track links simultaneously");
    }
    if (myOptions.trackLoopPairs) {
        throw std::runtime_error("Don't yet have ability to track looping pairs and do allophanate reactions simultaneously");
    }

    // here is code for if two chains react... need the case of 1 chain and 1 monomer as well...
    if (freeGroupLoc.whichChain!=-1) {   
        if (myOptions.trackFullStructureAsGraph) {
            // after the old link (e.g., urethane) and free group (e.g., isocyanate) have been found, it is time to (1) change the name of the old link to "allophanate" (2)
            // include pointers from one link to the other so that one can tell which two links are a package deal (3) copy freeGrp chain into the oldLink chain
            Edge oldEdge = boost::edge(oldLinkDeets.src,oldLinkDeets.tgt,cSS.chainPool[oldLinkDeets.chainIndex].polymer).first;
            cSS.chainPool[oldLinkDeets.chainIndex].polymer[oldEdge].name = newLink; // rename the edge ... would like to introduce pointers next...

            // now create the new edge... let's keep the source of the old edge as the one to be the source of the newly created 2nd edge?
            Vertex src = source(oldEdge,cSS.chainPool[oldLinkDeets.chainIndex].polymer);
            cSS.chainPool[oldLinkDeets.chainIndex].polymer[src].chain1=true; // mark the src vertex in case it changes
            Vertex tgt = freeGroupLoc.whichVertex;
            cSS.chainPool[freeGroupLoc.whichChain].polymer[tgt].chain2=true;

            // copy chain 2 into chain 1
            copy_graph(cSS.chainPool[freeGroupLoc.whichChain].polymer,cSS.chainPool[oldLinkDeets.chainIndex].polymer); // copy polymer containing isocyanate into polymer containing urethane
            // need to add up counters here as well (this needs to be done regardless of whether tracking full structure as a graph)

            // check whether src and tgt are still src and tgt...
            if (cSS.chainPool[oldLinkDeets.chainIndex].polymer[src].chain1==false) {
                // go find src by looping thru vertices...
                for (auto vxs = vertices(cSS.chainPool[oldLinkDeets.chainIndex].polymer);vxs.first!=vxs.second;vxs.first++) {
                    if (cSS.chainPool[oldLinkDeets.chainIndex].polymer[*vxs.first].chain1) {
                        src = *vxs.first;
                        break;
                    }
                }
            }
            cSS.chainPool[oldLinkDeets.chainIndex].polymer[src].chain1=false; // reset chain1 label back to false
            // now do the same for tgt
            if (cSS.chainPool[oldLinkDeets.chainIndex].polymer[tgt].chain2==false) {
                // go find tgt by looping thru vertices...
                for (auto vxs = vertices(cSS.chainPool[oldLinkDeets.chainIndex].polymer);vxs.first!=vxs.second;vxs.first++) {
                    if (cSS.chainPool[oldLinkDeets.chainIndex].polymer[*vxs.first].chain2) {
                        tgt = *vxs.first;
                        break;
                    }
                }
            }
            cSS.chainPool[oldLinkDeets.chainIndex].polymer[tgt].chain2=false; // reset chain2 label back to false

            // now add an edge between src and tgt
            Edge newEdge = boost::add_edge(src,tgt,cSS.chainPool[oldLinkDeets.chainIndex].polymer).first;
            cSS.chainPool[oldLinkDeets.chainIndex].polymer[newEdge].name=newLink; // name the edge

            // update trackers on monomers...
            cSS.chainPool[oldLinkDeets.chainIndex].polymer[src].funs[oldLink]--;
            cSS.chainPool[oldLinkDeets.chainIndex].polymer[src].funs[newLink]++;
            cSS.chainPool[oldLinkDeets.chainIndex].polymer[tgt].funs[freeGrp]--;

        }
        // regardless of whether we're tracking the full structure as a graph, need to update moiety counters
        cSS.numMoieties[oldLink]--;
        cSS.numMoieties[freeGrp]--;
        cSS.numMoieties[newLink]++;
        // funs on new chain should be sum of the old two chains, with the following modification...
        cSS.chainPool[oldLinkDeets.chainIndex].funs[oldLink]--;
        cSS.chainPool[oldLinkDeets.chainIndex].funs[newLink]++;
        cSS.chainPool[oldLinkDeets.chainIndex].funs[freeGrp]--;
        // after the modifications above, funs on new chain are the sum of the old two chains
        for(auto it = cSS.chainPool[freeGroupLoc.whichChain].funs.begin(); it != cSS.chainPool[freeGroupLoc.whichChain].funs.end(); ++it) {
            cSS.chainPool[oldLinkDeets.chainIndex].funs[it->first] += it->second; 
        }

        // if tracking the full structure, need to update functional group counts on specific vertices, but this should have been done above.

        // figure out molecular weight stuff here

        // update number average and weight average molecular weight
        cSS.mwv.Mi_1=cSS.chainPool[oldLinkDeets.chainIndex].weight;
        cSS.mwv.Mi_2=cSS.chainPool[freeGroupLoc.whichChain].weight;

        // update weight of remaining chain
        cSS.chainPool[oldLinkDeets.chainIndex].weight+=cSS.chainPool[freeGroupLoc.whichChain].weight;
        
        // delete the old chain
        cSS.chainPool.erase(cSS.chainPool.begin()+freeGroupLoc.whichChain);
    }
    else {
        // in this case, it's a monomer plus a chain
        if ( freeGroupLoc.whichChain!=-1 ) {
            throw std::runtime_error("unexpected behavior!");
        }
        Vertex tgt = add_vertex(cSS.chainPool[oldLinkDeets.chainIndex].polymer);
        cSS.chainPool[oldLinkDeets.chainIndex].polymer[tgt] = cSS.monomers[freeGroupLoc.monomerSelected];

        // go ahead and add the edge
        Edge newEdge = add_edge(oldLinkDeets.src,tgt,cSS.chainPool[oldLinkDeets.chainIndex].polymer).first;
        cSS.chainPool[oldLinkDeets.chainIndex].polymer[newEdge].name = newLink; // name the new edge
        Edge oldEdge = boost::edge(oldLinkDeets.src,oldLinkDeets.tgt,cSS.chainPool[oldLinkDeets.chainIndex].polymer).first;
        cSS.chainPool[oldLinkDeets.chainIndex].polymer[oldEdge].name = newLink; // rename the old edge ... would like to introduce pointers next...
        
        // update trackers on monomers...
        cSS.chainPool[oldLinkDeets.chainIndex].polymer[oldLinkDeets.src].funs[oldLink]--;
        cSS.chainPool[oldLinkDeets.chainIndex].polymer[oldLinkDeets.src].funs[newLink]++;
        cSS.chainPool[oldLinkDeets.chainIndex].polymer[tgt].funs[freeGrp]--;

        // update number average and weight average molecular weight
        cSS.mwv.Mi_1=cSS.chainPool[oldLinkDeets.chainIndex].weight;
        cSS.mwv.Mi_2=cSS.monomers[freeGroupLoc.monomerSelected].weight;

        // update weight of polymer chain
        cSS.chainPool[oldLinkDeets.chainIndex].weight += cSS.monomers[freeGroupLoc.monomerSelected].weight;

        // regardless of whether we're tracking the full structure as a graph, need to update moiety counters
        cSS.numMoieties[oldLink]--;
        cSS.numMoieties[freeGrp]--;
        cSS.numMoieties[newLink]++;
        // funs on new chain should be sum of the old two chains, with the following modification...
        cSS.chainPool[oldLinkDeets.chainIndex].funs[oldLink]--;
        cSS.chainPool[oldLinkDeets.chainIndex].funs[newLink]++;
        cSS.chainPool[oldLinkDeets.chainIndex].funs[freeGrp]--;
        // after the modifications above, funs on new chain are the sum of those on the old chain and on the added monomer
        for(auto it = cSS.monomers[freeGroupLoc.monomerSelected].funs.begin(); it != cSS.monomers[freeGroupLoc.monomerSelected].funs.end(); ++it) {
            cSS.chainPool[oldLinkDeets.chainIndex].funs[it->first] += it->second; 
        }
        // update number of monomers remaining
        cSS.numMonomers[freeGroupLoc.monomerSelected]--;

    }

    // no links figured out 
    // no loop pairs figured out

    return;
}

void separateFragments(std::vector<Chain> &chainPool,long chainIndex,Connectivity myConnect) {
    // This function breaks a chain into two separate pieces, if there is no longer a topological 
    // connection between the monomers with the deleted edge. 
    
    // TODO notes:
    // the value returned could be the index of the new chain, but for now, the new chain 
    // is always the last chain in the chainPool

    // also ideally the source *should* probably always be put on the original chain and target 
    // always on the new chain, but I haven't yet implemented something like that
    
    // there is nothing for this function to do if the chain to be separated into two fragments is all connected
    if (myConnect.isSingleComponent) {
        throw std::runtime_error("yoo this function should only be called when myConnect.isSingleComponent is true");
    }
    Chain newFragment=chainPool[chainIndex]; // save the chain as a "new fragment" -- delete pieces from this to create the new chain
    Graph d=newFragment.polymer; // save the graph with the edge removed as a new copy named d

    // now delete all the ones *not connected* to component 0 from the *old* graph
    // delete all the ones *connected* to component 0 from the *new* graph
    int countDeletedEdges=0; // count deleted edges (corresponding to component 1) on old graph (component 0)
    int countDeletedEdgesNewGraph=0; // count deleted edges (corresponding to component 0) on new graph (component 1)
    int totalVertices=boost::num_vertices(chainPool[chainIndex].polymer); // start and end vertices
    for (int i = 0; i < totalVertices; ++i) {

        if (myConnect.component[i] != 0) {
            int index = i-countDeletedEdges; // calculate index in the old chain
            // update chain's functional group count at the end...
            boost::clear_vertex(i-countDeletedEdges,chainPool[chainIndex].polymer);
            boost::remove_vertex(i-countDeletedEdges,chainPool[chainIndex].polymer);
            countDeletedEdges++;
        }
        else {
            int index=i-countDeletedEdgesNewGraph; // calculate index of the new chain
            // calculate the effect of deleting the particular monomer
            for (auto itfun = d[index].funs.begin(); itfun != d[index].funs.end(); itfun++ ) {
                newFragment.funs[itfun->first]-=itfun->second;
            } 
            // update chain weight
            newFragment.weight-=d[index].weight;
            boost::clear_vertex(i-countDeletedEdgesNewGraph,d); // delete all edges associated with the vertex
            boost::remove_vertex(i-countDeletedEdgesNewGraph,d);
            countDeletedEdgesNewGraph++;
        }
    }

    if (!myConnect.isSingleComponent) {
        // fix the fun group count on the old chain
        newFragment.polymer=d;
        for (auto it=chainPool[chainIndex].funs.begin(); it!=chainPool[chainIndex].funs.end();it++) {
            // funs on the old fragment are equal to the original number of funs minus the ones on the new fragment..
            it->second=it->second-newFragment.funs[it->first];
        }
        // fix the weight on the old chain
        chainPool[chainIndex].weight-=newFragment.weight; // note: this is the weight assuming no weight was lost and weight was divided between chains according to the original monomer weights
        // give newFragment a namemap
        newFragment.nmp=get(&EdgeData::name,newFragment.polymer); // following example https://stackoverflow.com/questions/3527313/using-bundled-properties-as-the-weight-map-in-dijkstra-shortest-paths
        // put newFragment in the chainPool
        chainPool.push_back(newFragment);
    }
    return;
}

void cleanScissionRxn(SystemVariables &cSS, Reaction rxn,double r3,double r4, ModelOptions myOptions,std::map<std::string,DoubleSec>& timers) {
    
    // For now, need to be tracking full structure to do this correctly...
    if (!myOptions.trackFullStructureAsGraph){
        // In the future, groups could be randomly distributed between the two new polymer fragments
        throw std::runtime_error("Need to track full structure!");
    }
    if (rxn.products.size()!=2) {
        // for now, breaking is of type C -> A + B... the code is not equipped to handle situations where more or less than 2 products result
        throw std::runtime_error("For now, need all breaking reactions to result in two products!");
    }
    if (myOptions.trackAllLinks) {
        // if the scission reactions are link-dependent, the function adjacentBreakingRxn
        // should be used. Link tracking could be implemented here, but we have not
        // yet encountered a scenario where this is needed, and so it hasn't been
        // developed or tested
        throw std::runtime_error("code not yet capable of handling this scenario");
    }
    // STEP 0: identify which linkage to break. One frustrating step is whether to treat it as a "clean" break or not--do any functional groups get left behind? For now, treat as "clean" break
    double r_break = r3; // in range [0,1)
    double whichLinkage = r_break*cSS.numMoieties[rxn.reactants[0]]; // in range [0, # of linkages of specified type)
    std::string linkage = rxn.reactants[0]; // store the identity of the linkage to break in this reaction, which should be stored as a functional group on the source monomer
    std::string product1=rxn.products[0]; // store the identity of the produced functional group on the source monomer
    std::string product2=rxn.products[1]; // store the identity of the produced functional group on the target monomer

    EdgeParticipants selectedEdge=linearSearchForMatchingEdge(cSS.chainPool,whichLinkage, linkage,myOptions.trackFullStructureAsGraph);
    int whichChain=selectedEdge.chainIndex; // chain number within cSS.chainPool to react. 
    Vertex tgt = selectedEdge.tgt; // target vertex of the edge representing the bond to be broken
    Vertex src = selectedEdge.src; // source vertex of the edge representing the bond to be broken

    Chain newFragment; // new chain fragment resulting from break 
    Graph d; // graph to copy the topology of the new fragment into
                    
    // STEP 4 adjust functional group counters on both monomers participating in the linkage
    // note: the linkage should always be registered as a functional group on the *source* monomer
    cSS.chainPool[whichChain].polymer[src].funs[linkage]--; // decrement the linkage functional group count on the source monomer 
    // add respective products to src and target
    cSS.chainPool[whichChain].polymer[src].funs[product1]++; // increment by the functional group produced on the source monomer
    
    cSS.chainPool[whichChain].polymer[tgt].funs[product2]++; // increment by the functional group produced on the target monomer

    // STEP 5 adjust functional group counters on the chain being broken apart
    // 5.1 - add products to the chain counters
    cSS.chainPool[whichChain].funs[product1]++; // increment by the functional group produced on the source monomer
    cSS.chainPool[whichChain].funs[product2]++; // increment by the functional group produced on the target monomer
    // 5.2 - remove linkage on the chain itself
    if (!cSS.chainPool[whichChain].funs.count(linkage)) {
        // we have an issue of there are no linkages of the type selected on the chain
        throw std::runtime_error( "Oh boy.");
    }
    cSS.chainPool[whichChain].funs[linkage]--; // decrement by the linkage group
    
    // STEP 6 delete the edge!
    boost::remove_edge(src,tgt,cSS.chainPool[whichChain].polymer);

    // STEP 7 break this chain into two separate pieces, if there is no longer a topological connection between the monomers with the deleted edge. 
    
    // check whether the polymer is still a single entity
    Connectivity myConnect=checkConnectivity(cSS.chainPool[selectedEdge.chainIndex].polymer); 
    // if the polymer is not a single entity
    if (!myConnect.isSingleComponent) {
        // then separate the chain into two fragments
        separateFragments(cSS.chainPool,selectedEdge.chainIndex,myConnect);
        // and update the molecular weight accordingly
        cSS.mwv.Mi_1=cSS.chainPool[selectedEdge.chainIndex].weight;
        cSS.mwv.Mi_2=cSS.chainPool.back().weight;
    }
    
    // STEP 8 update overall moiety counts
    cSS.numMoieties[product1]++;
    cSS.numMoieties[product2]++;
    cSS.numMoieties[linkage]--;

    // end execution of breaking reaction
    return;
}


std::pair<bool,std::string> isMonomer(Chain myChain) {
    // this function returns true if myChain is a monomer. Otherwise, it returns false.
    // the weight of this code is currently hard-coded for changxia's case
    std::pair<bool,std::string> result = {false,""};
    double DP = myChain.weight/142.0;
    if (DP>1.1) {
        return result;
    }
    else if (myChain.weight>141.99 && myChain.weight < 142.01) {
        // find monomer name
        std::string monomerName;
        for (auto j=myChain.polymer[0].funs.begin();j!=myChain.polymer[0].funs.end();++j) {
            if (j->second != 0) {
                result.second=j->first;
                break;
            }
        }
        result.first=true;
    }
    else {
        throw std::runtime_error("help me plz (maybe running code for which hard-coded lactone values are no longer appropriate)");
    }
    return result;
}

// void catBindingRxn(SystemVariables &cSS, Reaction rxn,double r3,double r4,ModelOptions myOptions) {
//     // this is written in the context of Changixa's system

//     // step 1: identify the group 
//     GroupLocation myGroup = linearSearchForParticipatingGroup(cSS,r3,rxn.reactants[0],myOptions.trackFullStructureAsGraph);
//     if (myGroup.whichChain==-1) {
//         throw std::runtime_error("for changxia's system, chosen group for cat binding should not be on monomer... ");
//     }

//     // step 2: swap the group
//     // on chain
//     cSS.chainPool[myGroup.whichChain].funs[rxn.reactants[0]]--;
//     cSS.chainPool[myGroup.whichChain].funs[rxn.products[0]]++;
//     // on monomer
//     if (myOptions.trackFullStructureAsGraph) {
//         cSS.chainPool[myGroup.whichChain].polymer[myGroup.whichVertex].funs[rxn.reactants[0]]--;
//         cSS.chainPool[myGroup.whichChain].polymer[myGroup.whichVertex].funs[rxn.products[0]]++;
//     }
//     // overall
//     cSS.numMoieties[rxn.reactants[0]]--;
//     cSS.numMoieties[rxn.products[0]]++;

//     // update adjacent links
//     updateLinksAdjacentToOneVertex(cSS.links,cSS.chainPool[myGroup.whichChain], rxn, myGroup.whichVertex,myOptions);
    
//     // note that a catBindingRxn has occurred!
//     if (rxn.catalystEffect=="bound") {
//         cSS.numMoieties["*"]--;
//     }
//     else if (rxn.catalystEffect=="released") {
//         cSS.numMoieties["*"]++;
//     }
//     cSS.catBindingTally++;

//     // i think that's it!
// }

void adjacentBreakingRxn(SystemVariables &cSS, Reaction rxn,double r3,double r4,ModelOptions myOptions) {

    // this is the code to execute for reactions in which the correct link to cleave depends on the functional
    // groups present on both monomers attached to the link. This is used to simulate a hyperbranched lactone system

    // first, perform a linear search in the chainPool to find the correct link
    std::string reactant1=rxn.reactants[0];
    std::string reactant2=rxn.reactants[1];
    std::string product1=rxn.products[0];
    std::string product2=rxn.products[1];
    std::string fullLinkNameToBust=reactant1+'-'+rxn.linkToBust+'-'+reactant2;
    double whichLink=r3*cSS.links[fullLinkNameToBust]; // generate random number on range [0,# of links). (including 0, not including # of links)
    int linkCounter=0;
    int polymerToBust=-1;
    //boost::graph_traits<Graph>::edge_iterator edgeToBust;
    Vertex srcBust=-1;
    Vertex tgtBust=-1;
    bool breakOutOfLoop=false;
    for (int i=0;i<cSS.chainPool.size();++i) {
        if (cSS.chainPool[i].linkTotals.count(fullLinkNameToBust)) {
            linkCounter+=cSS.chainPool[i].linkTotals[fullLinkNameToBust];
        }
        if (linkCounter>whichLink) {
            linkCounter-=cSS.chainPool[i].linkTotals[fullLinkNameToBust];

            // examine every edge in cSS.chainPool[i]... need to obtain edge list
            // https://www.boost.org/doc/libs/1_79_0/libs/graph/doc/EdgeListGraph.html
            auto polymer=cSS.chainPool[i].polymer;
            auto edgs=edges(polymer);
            while (edgs.first!=edgs.second) {
                Vertex src=source(*edgs.first,polymer);
                Vertex tgt=target(*edgs.first,polymer);
                // determine edge name
                std::string srcFun="";
                std::string tgtFun="";
                // this code is specific to the Changxia case -- in the future just use monomer names instead of assuming there is only one functional group per monomer?
                for (auto j=polymer[src].funs.begin();j!=polymer[src].funs.end();++j) {
                    if (j->second != 0) {
                        srcFun=j->first;
                        break;
                    }
                }
                for (auto j=polymer[tgt].funs.begin();j!=polymer[tgt].funs.end();++j) {
                    if (j->second != 0) {
                        tgtFun=j->first;
                        break;
                    }
                }
                if (srcFun=="" || tgtFun=="") {
                    throw std::runtime_error("there is a problem. most likely stems from an issue with edge counting. ");
                }
                std::string currentLink=srcFun+'-'+cSS.chainPool[i].nmp[*edgs.first]+'-'+tgtFun;
                if (currentLink==fullLinkNameToBust) {
                    linkCounter++;
                }
                if (linkCounter>whichLink) {
                    polymerToBust=i;
                    srcBust=src;
                    tgtBust=tgt;
                    breakOutOfLoop=true;
                    break;
                }
                edgs.first++;
            }
            if (breakOutOfLoop) {
                break;
            }
        }
    }
    if (breakOutOfLoop==false) {
        throw std::runtime_error("problem! edge to break not found. probably stems from edge counting issue");
    }
    // let's go ahead and break the edge, go through the same procedure as when there's one product. 

    // STEP 4 adjust functional group counters on both monomers participating in the linkage
    cSS.chainPool[polymerToBust].polymer[srcBust].funs[reactant1]--;
    cSS.chainPool[polymerToBust].polymer[srcBust].funs[product1]++;
    cSS.chainPool[polymerToBust].polymer[tgtBust].funs[reactant2]--;
    cSS.chainPool[polymerToBust].polymer[tgtBust].funs[product2]++;

    // STEP 5 adjust functional group counters on the chain being broken apart
    // 5.1 - add products to the chain counters
    cSS.chainPool[polymerToBust].funs[product1]++;
    cSS.chainPool[polymerToBust].funs[product2]++;
    // 5.2 - remove reactants from the chain itself
    cSS.chainPool[polymerToBust].funs[reactant1]--;
    cSS.chainPool[polymerToBust].funs[reactant2]--;
    // if tracking Changxia-like linkages on chain, remove it here
    

    // STEP 6 delete the edge!
    boost::remove_edge(srcBust,tgtBust,cSS.chainPool[polymerToBust].polymer);

    int linkCat=0;
    for (auto it=cSS.links.begin(); it!=cSS.links.end(); ++it) {
        if ((it->first).find('*') != std::string::npos) {
            linkCat+=it->second;
        }
    }
    // 6.1 update link counts
    // first, decrement links tracker by link busted 
    cSS.links[reactant1+'-'+rxn.linkToBust+'-'+reactant2]--; 
    cSS.chainPool[polymerToBust].linkTotals[reactant1+'-'+rxn.linkToBust+'-'+reactant2]--;

    // update the linkage pairs affected by the new chain units   (first must identify the old linkage, must have tracked src and target.. )
    EdgeParticipants modifiedEdge={srcBust,tgtBust,-1}; // later, might consider using linear search function for matching edge here...
    updateAdjacentLinks(cSS.links,cSS.chainPool[polymerToBust], rxn,modifiedEdge,myOptions);

    // STEP 7 break this chain into two separate pieces, if there is no longer a topological connection between the monomers with the deleted edge. 
    bool edgeRemoved=false; // has an edge successfully been removed?
    bool Fragmentation=true; // is the chain still intact after removing the selected edge?
    Chain newFragment=cSS.chainPool[polymerToBust]; // new chain fragment resulting from break
    Graph d=newFragment.polymer; // graph to copy the topology of the new fragment into (copy all then delete extraneous vertices and edges)
    // create undirected graph to check for connectivity
    UndirectedGraph uG;
    copy_graph(d, uG); // there may be a more efficient way of creating an undirected graph from a directed graph
    std::vector<int> component (boost::num_vertices (uG)); // declare vector that will label connected pieces ("components") of the graph with different integers (starting with 0)
    auto num_components = boost::connected_components (uG, &component[0]); // determine what topological fragment each vertex in the polymer chain belongs to
    // count the number of items in each fragment
    int countComp1=0;
    int countComp2=0;
    int monomerIndex=-1;
    int possibleMonomerIndex1=-1;
    int possibleMonomerIndex2=-1;
    for (int i = 0; i<component.size(); ++i) {
        if (component[i]==0) {
            countComp1++;
            if (countComp1==1) {
                possibleMonomerIndex1=i;
            }
            
        }
        else {
            countComp2++;
            if (countComp2==1) {
                possibleMonomerIndex2=i;
            }
        }
    }
    if (num_components>2) {
        throw std::runtime_error( "Oh boy. Shouldn't have more than two fragments after one reverse reaction...");
    }
    if (num_components==1) {
        // if everything is still connected, no need to proceed further
        Fragmentation=false;
    }
    else if (countComp1==1) {
        // it's just a monomer that is coming off! Just erase the one...
        // delete the one vertex from the old graph!
        boost::clear_vertex(possibleMonomerIndex1,cSS.chainPool[polymerToBust].polymer);
        boost::remove_vertex(possibleMonomerIndex1,cSS.chainPool[polymerToBust].polymer);

        // literally just update functional units, no links to update bc monomer and already adjusted by the deleted edge...
        cSS.chainPool[polymerToBust].funs["cEp"]--;
        // adjust molecular weight of chain
        cSS.chainPool[polymerToBust].weight-=cSS.monomers["cEp"].weight;
        // increment monomer count
        cSS.numMonomers["cEp"]++;
        // set Mi1 and Mi2...
        cSS.mwv.Mi_1=cSS.monomers["cEp"].weight;
        cSS.mwv.Mi_2=cSS.chainPool[polymerToBust].weight;
    }
    else if (countComp2==1) {
        // it's just a monomer that is coming off! Just erase the one...
        // delete the one vertex from the old graph!
        boost::clear_vertex(possibleMonomerIndex2,cSS.chainPool[polymerToBust].polymer);
        boost::remove_vertex(possibleMonomerIndex2,cSS.chainPool[polymerToBust].polymer);

        // literally just update functional units, no links to update bc monomer and already adjusted by the deleted edge...
        cSS.chainPool[polymerToBust].funs["cEp"]--;
        
        // adjust molecular weight of chain
        cSS.chainPool[polymerToBust].weight-=cSS.monomers["cEp"].weight;
        // increment monomer count
        cSS.numMonomers["cEp"]++;
        // set Mi1 and Mi2...
        cSS.mwv.Mi_1=cSS.monomers["cEp"].weight;
        cSS.mwv.Mi_2=cSS.chainPool[polymerToBust].weight;
    }
    else if (Fragmentation) {
        // now delete all the ones *not connected* to component 0 from the *old* graph
        // delete all the ones *connected* to component 0 from the *new* graph
        int countDeletedEdges=0; // count deleted edges (corresponding to component 1) on old graph (component 0)
        int countDeletedEdgesNewGraph=0; // count deleted edges (corresponding to component 0) on new graph (component 1)
        int totalVertices=boost::num_vertices(cSS.chainPool[polymerToBust].polymer); // start and end vertices
        for (int i = 0; i < totalVertices; ++i) {

            if (component[i] != 0) {
                int index = i-countDeletedEdges; // calculate index in the old chain
                // update chain's functional group count at the end...
                // update links counter if that's being done here

                boost::clear_vertex(i-countDeletedEdges,cSS.chainPool[polymerToBust].polymer);
                boost::remove_vertex(i-countDeletedEdges,cSS.chainPool[polymerToBust].polymer);
                countDeletedEdges++;
            }
            else {
                int index=i-countDeletedEdgesNewGraph; // calculate index of the new chain
                // calculate the effect of deleting the particular monomer on functional group counts
                for (auto itfun = d[index].funs.begin(); itfun != d[index].funs.end(); itfun++ ) {
                    newFragment.funs[itfun->first]-=itfun->second;
                } 
                // calculate the effect of deleting the particular monomer on link counts
                
                if (myOptions.trackAllLinks) {
                    // check out edges first
                    for (auto out_edgs = boost::out_edges(index, d); out_edgs.first!=out_edgs.second; ++out_edgs.first) {
                        std::string linkName = getLinkName(d,*out_edgs.first);
                        newFragment.linkTotals[linkName]--;
                    }

                    //check in edges -- is this redundant ?
                    for (auto in_edgs = boost::in_edges(index, d); in_edgs.first!=in_edgs.second; ++in_edgs.first) {
                        // make sure that this is accurately describing src and tgt in the getLink function
                        std::string linkName = getLinkName(d,*in_edgs.first);
                        newFragment.linkTotals[linkName]--;
                    }
                }

                // update chain weight
                newFragment.weight-=d[index].weight;
                boost::clear_vertex(i-countDeletedEdgesNewGraph,d); // delete all edges associated with the vertex
                boost::remove_vertex(i-countDeletedEdgesNewGraph,d);
                countDeletedEdgesNewGraph++;
            }
        }
        //std::cout << cSS.executionNumber << std::endl;
        // STEP 8 update molecular weight information, property maps, and old chain functional group counters (calculated by subtracting new chain functional groups)
        // fix the fun group count on the old chain
        newFragment.polymer=d;
        for (auto it=cSS.chainPool[polymerToBust].funs.begin(); it!=cSS.chainPool[polymerToBust].funs.end();it++) {
            // funs on the old fragment are equal to the original number of funs minus the ones on the new fragment..
            it->second=it->second-newFragment.funs[it->first];
        }
        // fix the weight on the old chain
        cSS.chainPool[polymerToBust].weight-=newFragment.weight;

        // fix the links on the old chain
        if (myOptions.trackAllLinks) {
            // links on the old fragment are equal to the original number minus the ones on the new fragment
            // (make sure to also adjust by the link that is broken)

            // cycle through all the links
            for (auto it=cSS.chainPool[polymerToBust].linkTotals.begin();it!=cSS.chainPool[polymerToBust].linkTotals.end(); ++it) {
                cSS.chainPool[polymerToBust].linkTotals[it->first]=it->second-newFragment.linkTotals[it->first];
            }
        }
        // update namemap on newFragment
        cSS.chainPool[polymerToBust].nmp=get(&EdgeData::name,newFragment.polymer); // following example https://stackoverflow.com/questions/3527313/using-bundled-properties-as-the-weight-map-in-dijkstra-shortest-paths
        // put newFragment in the cSS.chainPool (or add to numMonomers, if appropriate)
        std::pair<bool,std::string> checkMonomer = isMonomer(newFragment);
        if (checkMonomer.first==false) {      
            cSS.chainPool.push_back(newFragment); 
        }
        else {
            // grab monomer name
            cSS.numMonomers[checkMonomer.second]++;
        }

        // update the molecular weight
        cSS.mwv.Mi_1=newFragment.weight;
        cSS.mwv.Mi_2=cSS.chainPool[polymerToBust].weight;

        // if the old fragment is a monomer, remove it from the chainPool and put it back into the monomer pile
        checkMonomer = isMonomer(cSS.chainPool[polymerToBust]);
        if (checkMonomer.first==true) {
            throw std::runtime_error("wait what this should not happen");
            // cSS.numMonomers[checkMonomer.second]++;
            // cSS.chainPool.erase(cSS.chainPool.begin()+polymerToBust);
        }

    }
    
    // STEP 9 update overall moiety counts
    cSS.numMoieties[product1]++;
    cSS.numMoieties[product2]++;
    cSS.numMoieties[reactant1]--;
    cSS.numMoieties[reactant2]--;

    // // STEP 10 update catalyst count
    // if (rxn.catalystEffect=="released" || rxn.linkToBust=="[eplinkNotFormed]*") {
    //     // if (rxn.catalystEffect=="released") {
    //     //     // it's not odd if a BREAKING reaction RELEASES catalyst
    //     //     // in Changxia's system, a linking reaction should only be a BINDING step.
    //     //     throw std::runtime_error("surprised by this");
    //     // } if statement should probably consider binding of catalyst instead of link type lol
    //     cSS.numMoieties["*"]++;
    //     cSS.catBindingTally++;

    // }
    // else {
    //     cSS.catBindingTally=0;
    // }

    // possibly here send newFragment / old chain to numMonomers if it's just a cEp
    return;
}

void carbamateDecompositionRxn(SystemVariables &cSS, Reaction rxn,double r3,double r4,ModelOptions myOptions,std::map<std::string,DoubleSec>& timers) {
    // this function encodes the reaction in which a carbamate group decomposes into an alkene, an amine, and CO2.
    // assume reaction in input file is written as NCO-OH -> NH2 + C=C + CO2(g)
    // what was originally NCO (src) loses (16+12) in molecular weight
    // what was originally OH (tgt) loses 16 in molecular weight

    // for now, assume concentrations don't change, and just pop CO2 into a "gas evolved" tracker
    if (rxn.gasEvolved=="") {
        throw std::runtime_error("need to specify gas evolved to carry out a decomposition reaction emitting a gaseous molecule");
        return;
    }
    // for now, function not capable of considering the case in which the full structure is not tracked as a graph...
    if (myOptions.trackFullStructureAsGraph==false) {
        throw std::runtime_error("not yet implemented");
    }
    // step 1: find the link to break:
    if (rxn.linkToBust!=rxn.reactants[0]) {
        throw std::runtime_error("input function has not been written as expected");
    }
    EdgeParticipants chosenEdge = linearSearchForMatchingEdge(cSS.chainPool,r3*cSS.numMoieties[rxn.reactants[0]],rxn.reactants[0],myOptions.trackFullStructureAsGraph);

    // next, carry out the reaction
    // remove the doggone edge...
    boost::remove_edge(chosenEdge.src,chosenEdge.tgt,cSS.chainPool[chosenEdge.chainIndex].polymer);

    // check for connectivity and then split into separate chains...
    Connectivity myConnect = checkConnectivity(cSS.chainPool[chosenEdge.chainIndex].polymer);
    // update molecular weight trackers BEFORE calling separateFragements because 
    // afterwards it's going to be difficult to figure out which is which
    // assume that src is the isocyanate and tgt is the alcohol...
    cSS.chainPool[chosenEdge.chainIndex].polymer[chosenEdge.tgt].weight-=16; // update weight of alcohol
    cSS.chainPool[chosenEdge.chainIndex].polymer[chosenEdge.src].weight-=(16+12); // update weight of isocyanate
    cSS.chainPool[chosenEdge.chainIndex].weight-=(16+16+12); // update weight on chain tracker
    cSS.chainPool[chosenEdge.chainIndex].polymer[chosenEdge.tgt].funs[rxn.products[1]]++; // update tgt with products (alcohol -> alkene)
    cSS.chainPool[chosenEdge.chainIndex].polymer[chosenEdge.src].funs[rxn.products[0]]++; // update src with products (isocyanate -> amine)
    cSS.chainPool[chosenEdge.chainIndex].polymer[chosenEdge.src].funs[rxn.reactants[0]]--; // update src by subtracting reactants (isocyanate -> amine)

    // update chain counts too
    cSS.chainPool[chosenEdge.chainIndex].funs[rxn.products[0]]++;
    cSS.chainPool[chosenEdge.chainIndex].funs[rxn.products[1]]++;
    cSS.chainPool[chosenEdge.chainIndex].funs[rxn.reactants[0]]--;
        
    if (myOptions.trackFullStructureAsGraph && !myConnect.isSingleComponent) {
        separateFragments(cSS.chainPool,chosenEdge.chainIndex,myConnect);
    }

    // if (myOptions.trackFullStructureAsGraph==false) {
    //     // assume the chain scission is random
    // }
    cSS.moleculesOfGasEvolved[rxn.gasEvolved]++; // update the amount of gas present

    // update overall moiety counts
    cSS.numMoieties[rxn.reactants[0]]--;
    cSS.numMoieties[rxn.products[0]]++;
    cSS.numMoieties[rxn.products[1]]++;

    return;
}

EdgeParticipants findALink( SystemVariables &cSS, Reaction rxn,double r3,double r4 ) {
    // this function is different than the function that finds just an edge, because it is finding a *LINK* 
    // which in my notation means that the identity of the adjacent monomers (functional groups present on them) 
    // matter too
    // first, perform a linear search in the chainPool to find the correct link
    std::string reactant1=rxn.reactants[0];
    std::string reactant2=rxn.reactants[1];

    std::string fullLinkNameToBust=reactant1+'-'+rxn.linkToBust+'-'+reactant2;
    double whichLink=r3*cSS.links[fullLinkNameToBust]; // generate random number on range [0,# of links). (including 0, not including # of links)
    int linkCounter=0;
    EdgeParticipants myEdge;

    if (cSS.executionNumber==94446) {
        int breakhere=3;
    }
    bool breakOutOfLoop=false;
    for (int i=0;i<cSS.chainPool.size();++i) {
        if (cSS.chainPool[i].linkTotals.count(fullLinkNameToBust)) {
            linkCounter+=cSS.chainPool[i].linkTotals[fullLinkNameToBust]; 
        }
        if (linkCounter>whichLink) {
            linkCounter-=cSS.chainPool[i].linkTotals[fullLinkNameToBust];

            // examine every edge in cSS.chainPool[i]... need to obtain edge list
            // https://www.boost.org/doc/libs/1_79_0/libs/graph/doc/EdgeListGraph.html
            auto polymer=cSS.chainPool[i].polymer;
            auto edgs=edges(polymer);
            while (edgs.first!=edgs.second) {
                Vertex src=source(*edgs.first,polymer);
                Vertex tgt=target(*edgs.first,polymer);
                // determine edge name
                std::string srcFun="";
                std::string tgtFun="";
                // this code is specific to the Changxia case -- in the future just use monomer names instead of assuming there is only one functional group per monomer?
                for (auto j=polymer[src].funs.begin();j!=polymer[src].funs.end();++j) {
                    if (j->second != 0) {
                        srcFun=j->first;
                        break;
                    }
                }
                for (auto j=polymer[tgt].funs.begin();j!=polymer[tgt].funs.end();++j) {
                    if (j->second != 0) {
                        tgtFun=j->first;
                        break;
                    }
                }
                if (srcFun=="" || tgtFun=="") {
                    throw std::runtime_error("there is a problem. most likely stems from an issue with edge counting. ");
                }
                std::string currentLink=srcFun+'-'+cSS.chainPool[i].nmp[*edgs.first]+'-'+tgtFun;
                if (currentLink==fullLinkNameToBust) {
                    linkCounter++;
                }
                if (linkCounter>whichLink) {
                    myEdge.chainIndex=i;
                    myEdge.src=src;
                    myEdge.tgt=tgt;
                    breakOutOfLoop=true;
                    break;
                }
                edgs.first++;
            }
            if (breakOutOfLoop) {
                break;
            }
        }
    }
    if (breakOutOfLoop==false) {
        throw std::runtime_error("problem! edge to break not found. probably stems from edge counting issue");
    }
    return myEdge;
}
void transformLinkTypeRxn(SystemVariables &cSS, Reaction rxn,double r3,double r4, ModelOptions myOptions,std::map<std::string,DoubleSec>& timers) {
    
    // Expected behavior:
    // rxn.linkToBust -> rxn.linkFormed

    EdgeParticipants transformThisLink = findALink(cSS,rxn,r3,r4 );

    if (rxn.catalystEffect!="released" && rxn.catalystEffect!="bound") {

        // link is between src, tgt in reactants. They'll be the same src, tgt in products if rxn does not involve release or binding of a catalyst
        // this function can be used to change the link type (due to internal reaction), or to notate the binding of a catalyst

        // first step: identify the specific link to react. This can be done with a separate function, should've been used in the reverse reaction

        // access that specific edge using known src and tgt and transform it
        Edge transformThisEdge = boost::edge(transformThisLink.src,transformThisLink.tgt,cSS.chainPool[transformThisLink.chainIndex].polymer).first;
        cSS.chainPool[transformThisLink.chainIndex].polymer[transformThisEdge].name=rxn.linkFormed;

        // update tracker on chain
        std::string oldLink = rxn.reactants[0]+'-'+rxn.linkToBust+'-'+rxn.reactants[1];
        std::string newLink = rxn.reactants[0]+'-'+rxn.linkFormed+'-'+rxn.reactants[1];
        cSS.chainPool[transformThisLink.chainIndex].linkTotals[oldLink]--;
        cSS.chainPool[transformThisLink.chainIndex].linkTotals[newLink]++;

        // update tracker in system
        cSS.links[oldLink]--;
        cSS.links[newLink]++;

        // No need to update adjacent links counters because changing just the linkType and neither vertex changes no other edges. So we should be good now with this transformation`
        //updateAdjacentLinks(cSS.links,cSS.chainPool[transformThisLink.chainIndex], rxn, transformThisLink,myOptions,newEdge);
        
        // reset the catalyst binding tally if we have a reaction take place 
        // that is neither a binding nor an unbinding reaction
        // cSS.catBindingTally=0;
        }
    // update catalyst effect if relevant
    else {
        
        // if catalyst is bound or released, then the src and tgt linkages change...
        // access that specific edge using known src and tgt and transform it
        Edge transformThisEdge = boost::edge(transformThisLink.src, transformThisLink.tgt,cSS.chainPool[transformThisLink.chainIndex].polymer).first;
        cSS.chainPool[transformThisLink.chainIndex].polymer[transformThisEdge].name=rxn.linkFormed;

        // update tracker on chain
        std::string oldLink = rxn.reactants[0]+'-'+rxn.linkToBust+'-'+rxn.reactants[1];
        std::string newLink = rxn.products[0]+'-'+rxn.linkFormed+'-'+rxn.products[1];
        cSS.chainPool[transformThisLink.chainIndex].linkTotals[oldLink]--;
        cSS.chainPool[transformThisLink.chainIndex].linkTotals[newLink]++;

        // update tracker in system
        cSS.links[oldLink]--;
        cSS.links[newLink]++;

        // update moieties here
        
        // *overall* 
        cSS.numMoieties[rxn.reactants[0]]--;
        cSS.numMoieties[rxn.reactants[1]]--;
        cSS.numMoieties[rxn.products[0]]++;
        cSS.numMoieties[rxn.products[1]]++;

        // on *chain*
        cSS.chainPool[transformThisLink.chainIndex].funs[rxn.reactants[0]]--;
        cSS.chainPool[transformThisLink.chainIndex].funs[rxn.reactants[1]]--;
        cSS.chainPool[transformThisLink.chainIndex].funs[rxn.products[0]]++;
        cSS.chainPool[transformThisLink.chainIndex].funs[rxn.products[1]]++;

        // on *monomers*
        cSS.chainPool[transformThisLink.chainIndex].polymer[transformThisLink.src].funs[rxn.reactants[0]]--;
        cSS.chainPool[transformThisLink.chainIndex].polymer[transformThisLink.src].funs[rxn.products[0]]++;
        cSS.chainPool[transformThisLink.chainIndex].polymer[transformThisLink.tgt].funs[rxn.reactants[1]]--;
        cSS.chainPool[transformThisLink.chainIndex].polymer[transformThisLink.tgt].funs[rxn.products[1]]++;

        // update adjacent links
        updateAdjacentLinks(cSS.links,cSS.chainPool[transformThisLink.chainIndex], rxn, transformThisLink,myOptions,transformThisEdge);

    }

    // if (rxn.catalystEffect=="released") {
    //     cSS.numMoieties["*"]++; // increment number of catalyst if catalyst is released
    //     cSS.catBindingTally++;
    // }
    // else if (rxn.catalystEffect=="bound") {
    //     cSS.numMoieties["*"]--; // decrement number of * if catalyst is bound in this step
    //     cSS.catBindingTally++;
    // }
    // else {
        
    //     cSS.catBindingTally=0;
    // }
    // we good!
    return;
}


void updateSystemState(SystemVariables &cSS, Reaction rxn,double r3,double r4,ModelOptions myOptions, std::map<std::string,DoubleSec>& timers ) {

    if (cSS.executionNumber==4470 || cSS.executionNumber==4502 || cSS.executionNumber == 4512) {
        bool stopHereToo=true;
    }
    if (rxn.reactionType!="linking") {
        cSS.rejectedSteps=0;
    }
    if (rxn.reactionType=="exchange") {
        // execute an exchange reaction
        exchangeRxn(cSS,rxn,r3,r4,myOptions,timers);
    }
    else if (rxn.reactionType=="looping") {
        // execute a looping reaction
        //loopingRxn(cSS,rxn,r3,r4,myOptions);
        throw std::runtime_error("looping reactions called separately");
    }
    else if (rxn.reactionType=="linking" && rxn.reactants.size()==2) {
        // execute a non-looping linking reaction 
        // i.e., a bimolecular reaction linking two monomers together (a typical polymerization reaction)
        // DoubleTimePoint start = Time::now(); //get the time at this instant 
        linkingRxn(cSS,rxn,r3,r4,myOptions,timers);
        // DoubleTimePoint stop = Time::now(); //get the time at this instant 
        // timers["linking"] += (stop - start);
    }
    else if (rxn.reactionType=="link-linking") {
        // execute an allophanate-like reaction
        // i.e., a reaction where a link (carbamate) reacts with a free functional group (isocyanate)
        // to form a new crosslink known as an allophanate 
        // DoubleTimePoint start = Time::now(); //get the time at this instant 
        allophanateRxn(cSS,rxn,r3,r4,myOptions,timers);
        // DoubleTimePoint stop = Time::now(); //get the time at this instant 
        // timers["allophanate"] += (stop - start);
    }
    else if (rxn.reactionType=="breaking" && rxn.reactants.size()==1) {
        // execute a clean random scission reaction
        // i.e., where a link is cleanly broken into two functional groups
        // DoubleTimePoint start = Time::now(); //get the time at this instant 
        cleanScissionRxn(cSS,rxn,r3,r4,myOptions,timers);
        // DoubleTimePoint stop = Time::now(); //get the time at this instant 
        // timers["breaking"] += (stop - start);
    }
    else if (rxn.reactionType=="breaking" && rxn.reactants.size()==2 && rxn.products.size()==2) {
        // execute a random scission reaction in which the functional groups on monomers participating 
        // in the link to be broken must be kjnown
        adjacentBreakingRxn(cSS,rxn,r3,r4,myOptions);
    }
    else if (rxn.reactionType=="carbamatedecomposition") {
        // DoubleTimePoint start = Time::now(); //get the time at this instant 
        carbamateDecompositionRxn(cSS,rxn,r3,r4,myOptions,timers);
        // DoubleTimePoint stop = Time::now(); //get the time at this instant 
        // timers["carbamateDecomposition"] += (stop - start);
    }
    else if(rxn.reactionType=="transformLinkType") {
        transformLinkTypeRxn(cSS,rxn,r3,r4,myOptions,timers);
    }
}
