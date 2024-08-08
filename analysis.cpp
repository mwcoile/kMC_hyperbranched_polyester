#include"analysis.h"

void print_graph(std::vector<Chain> chainPool,std::string TimePoint) {
    // write graph to file -- currently unused
    //
    //     // previously took input boost::property_map<Graph, std::string EdgeData::*>::type namemap
    Graph myGraph;
    for (int i=0;i<chainPool.size();i++) {
        copy_graph(chainPool[i].polymer,myGraph);
    }
    std::ofstream theGraph;
    std::string name="graphReadout" + TimePoint + ".csv";
    time_t t = time(0);   // get current time
    struct tm * now = localtime( & t );
    char date[80];
    strftime (date,80,"%Y%m%d_%I%M%S%p_%Z",now);
    std::string str(date);
    std::string path = "./Output/";
    std::string graphfilename=path+date+name;
    theGraph.open(graphfilename);

    std::ofstream theGraphNodes;
    std::string nameNodes="graphReadoutNodes.csv";
    std::string graphfilenameNodes=path+date+nameNodes;
    theGraphNodes.open(graphfilenameNodes);
    auto epair = boost::edges(myGraph);
    theGraph << "Source,Target,Name \n";
    for (auto iter = epair.first; iter!=epair.second;iter++) {
        //std::cout << source(*iter,myGraph) << "         " << target(*iter,myGraph) << "         " << namemap[*iter] << "\n";
        //theGraph << source(*iter,myGraph) << "," << target(*iter,myGraph) << "," << namemap[*iter] << "\n";
        theGraph << source(*iter,myGraph) << "," << target(*iter,myGraph) << "," << myGraph[*iter].name << "\n";
    }
    auto vpair = boost::vertices(myGraph);
    theGraphNodes << "Id,Name \n";
    nodeCount=0;
    for (auto iter = vpair.first; iter!=vpair.second;iter++) {
        // get the name of the functional group
        std::string nodeName = getChainUnitName(myGraph[*iter]);
        theGraphNodes << *iter << "," << nodeName << "\n";
        nodeCount++;
    }
    int i = 0
    while (i<cSS.numMonomers["cEp"]) {
        theGraphNodes << nodeCount << ",cEp \n";
        nodeCount++;
        i++;
    }
    // std::copy(vpair.first, vpair.second,
    // //     std::ostream_iterator<boost::adjacency_list<>::vertex_descriptor>{
//     std::cout, "\n"});

    theGraph.close();
    theGraphNodes.close();
    return;
}
// void print_graph(std::vector<Chain> chainPool) {
//     // write graph to file -- currently unused
// 
//     // previously took input boost::property_map<Graph, std::string EdgeData::*>::type namemap
// 
//     Graph myGraph;
//     for (int i=0;i<chainPool.size();i++) {
//         copy_graph(chainPool[i].polymer,myGraph);
//     }
// 
//     std::ofstream theGraph;
//     std::string name="graphReadout.csv";
//     time_t t = time(0);   // get current time
//     struct tm * now = localtime( & t );
//     char date[80];
//     strftime (date,80,"%Y%m%d_%I%M%S%p_%Z",now); 
//     std::string str(date);
//     std::string path = "./Output/";
//     std::string graphfilename=path+date+name; 
//     theGraph.open(graphfilename);
//     
//     // iterate through all vertices and print them out
//     // auto vpair = boost::vertices(myGraph);
//     // for (auto iter = vpair.first; iter!=vpair.second; iter++) {
//     //     std::cout << "vertex: " << *iter << "  " << myGraph[*iter].name << "  " 
//     //     << myGraph[*iter].volume << "  " << "\n";
//     // }
//     //std::cout << " START EDGE LIST\n " ;
//     // iterate through all the edges and print them out (note, I would prefer to do this in one step, 
//     // with a list of connected vertices following every vertex)
//     auto epair = boost::edges(myGraph);
//     //std::cout << "Edges \nSource    Target \n";
//     theGraph << "Source,Target,Name \n";
//     for (auto iter = epair.first; iter!=epair.second;iter++) {
//         //std::cout << source(*iter,myGraph) << "         " << target(*iter,myGraph) << "         " << namemap[*iter] << "\n";
//         //theGraph << source(*iter,myGraph) << "," << target(*iter,myGraph) << "," << namemap[*iter] << "\n";
//         theGraph << source(*iter,myGraph) << "," << target(*iter,myGraph) << "," << myGraph[*iter].name << "\n";
//     }
//     // std::copy(vpair.first, vpair.second,
//     //     std::ostream_iterator<boost::adjacency_list<>::vertex_descriptor>{
//     //     std::cout, "\n"});
//     theGraph.close();
//     return;
// }

std::string removeNonAlphabet(std::string line) {

    std::string temp = "";

    for (int i = 0; i < line.size(); ++i) {
        if ((line[i] >= 'a' && line[i] <= 'z') || (line[i] >= 'A' && line[i] <= 'Z')) {
            temp = temp + line[i];
        }
    }

    return temp;
}

void recordChangxiaConcentrations(std::map<std::string,signed long> numMoieties,std::map<std::string,signed long> numMoieties0,double over_x,double simTime,std::ofstream &ChangxiaOutput, double V, double const Na, double Mn, double Mw) {
    std::vector<double> calculatedConcentration;
    if (4==4) {



        std::map<std::string, signed long> catalystCorrected;

        for (auto i = numMoieties.begin(); i!=numMoieties.end(); ++i) {
            std::string moietyName = i->first;
            moietyName = removeNonAlphabet(moietyName);
            catalystCorrected[moietyName] += i->second;
        }   

        calculatedConcentration.push_back((catalystCorrected["E"]*1.0)/(V*Na));
        calculatedConcentration.push_back((catalystCorrected["Em"]*1.0)/(V*Na));
        calculatedConcentration.push_back((catalystCorrected["Emp"]*1.0)/(V*Na));
        calculatedConcentration.push_back((catalystCorrected["Ep"]*1.0)/(V*Na));
        calculatedConcentration.push_back((catalystCorrected["cE"]*1.0)/(V*Na));
        calculatedConcentration.push_back((catalystCorrected["cEp"]*1.0)/(V*Na));

        // calculatedConcentration.push_back((numMoieties["NCO"]*1.0)/(V*Na));
        // calculatedConcentration.push_back((numMoieties["OH"]*1.0)/(V*Na));
        // calculatedConcentration.push_back((numMoieties["NCO-OH"]*1.0)/(V*Na));
        // calculatedConcentration.push_back((numMoieties["allophanate"]*1.0)/(V*Na));
        // calculatedConcentration.push_back((numMoieties["alkene"]*1.0)/(V*Na));
        // calculatedConcentration.push_back((numMoieties["NH2"]*1.0)/(V*Na));
        // calculatedConcentration.push_back((numMoieties["urea"]*1.0)/(V*Na));
        // calculatedConcentration.push_back((numMoieties["cE"]*1.0)/(V*Na));
        // calculatedConcentration.push_back((numMoieties["cEp"]*1.0)/(V*Na));


        //         int experimentalE=chainPool.size()+numMoieties["E"];
        //         double DB=(2.0*numMoieties["E"])/(numMoieties["Em"]+numMoieties["Ep"]+2*numMoieties["E"]);
        //         ChangxiaOutput << simTime << "  " << over_x << "   " << 
        //         calculatedConcentration[0] << "   " << calculatedConcentration[1] << "   " << calculatedConcentration[2] 
        //         << "   " << calculatedConcentration[3] << "   " << calculatedConcentration[4] << "   " << calculatedConcentration[5] << "   " << Mn << "   " << Mw << "   " << DB << "\n";

        ChangxiaOutput << simTime << "  " << over_x << "   " << 
                 calculatedConcentration[0] << "   " << calculatedConcentration[1] << "   " << calculatedConcentration[2] 
                 << "   " << calculatedConcentration[3] << "   " << calculatedConcentration[4] << "   " << calculatedConcentration[5] 
                 << "   " << Mn << "   " << Mw << "\n";
    }
    else {
        double XA=1.0-numMoieties["A"]/(1.0*numMoieties0["A"]);
        ChangxiaOutput << simTime << "   " << XA << "\n";
    }
}

void TirrellTheory(std::vector<Chain> &chainPool,std::map<std::string,signed long> &numMoieties,std::map<std::string,signed long> &numMoieties0,
std::map<std::string,signed long> &numMonomers,std::map<std::string,Monomer> &monomers,std::ofstream &TirrellOutput, double simTime, double Mn, double Mw ) {
    /* Calculate NnBB and NnCC as a function of q to 
    compare with Tirrell theoretical results. Note:
    Only applies for linear polymers. Good method of 
    confirming that code works as expected.

    Does not work as written if loops are allowed to form.
    */

    // first check that all polymers are linear
    std::map<std::string,Monomer>::iterator it_tirrell=monomers.begin();
    // loop through each entry in the functional group map
    while (it_tirrell != monomers.end()) {
        if (it_tirrell->second.funs.size() > 2) {
            throw std::runtime_error( "Ruh roh! Can't do tirrell sequence analysis with branched polymers." );
        }
        it_tirrell++;
    }

    // declare variables needed for calculating NnBB and NnCC
    int chainNumber=0;
    int monomer_number=0;
    std::string currentSequence="";
    int sequenceLength=0;
    int num_seqs_BB=0;
    int num_seqs_CC=0;
    double NnBB=0;
    double NnCC=0;
    int total_sequence_length_BB=0;
    int total_sequence_length_CC=0;
    while (chainNumber < chainPool.size()) {
        //std::cout << "edges for chain number " << chainNumber << ": \n";

        boost::graph_traits<UndirectedGraph>::out_edge_iterator e, e_end;
        Graph currentPolymerDirected = chainPool[chainNumber].polymer;
        UndirectedGraph currentPolymer;
        copy_graph(currentPolymerDirected,currentPolymer);
        boost::graph_traits<UndirectedGraph>::vertex_descriptor 
            s = vertex(0, currentPolymer);
        // for (int i=0;i<num_vertices(currentPolymer);i++) {
        //     std::cout << "incident edges for vertex " << i << "\n";
        //     for (tie(e, e_end) = out_edges(i, currentPolymer); e != e_end; ++e)
        //         std::cout << "(" << source(*e, currentPolymer) 
        //         << "," << target(*e, currentPolymer) << ")" << "\n";
        // }

        // the above is just for visualization. Now, the next step is to actually do the Tirrell analysis
        chainNumber++;
        // Start with Vertex 0 and get the edge or edges that come out from it...
        boost::graph_traits<UndirectedGraph>::vertex_descriptor connection;

        // declare variables in order to properly count sequences in the middle of the chain...
        int timesThroughLoop=0; 
        std::string firstDirection=""; // in the first direction, either vertex s or connection are BB/CC. 
        std::string secondDirection=""; // in the second direction, one of the next two (excluding s) are BB/CC
        for (tie(e, e_end) = out_edges(s, currentPolymer); e != e_end; ++e) {
            /* first time through this loop does the first direction from the starting vertex.
            second time through this loop does the second direction from the starting vertex
            */
            timesThroughLoop+=1;
            if (timesThroughLoop>2) {
                throw std::runtime_error( "Ruh roh! Nonlinear polymers!" );
            }

            // for edge e, find the vertex descriptor of the monomer connected to the original monomer
            if (source(*e, currentPolymer) != s) {
                connection=source(*e, currentPolymer);
            }
            else if (target(*e, currentPolymer) != s) {
                connection=target(*e, currentPolymer);
            }
            else {
                throw std::runtime_error( "Ruh roh!" );
            }
            // track identity of sequence in the first direction
            if (timesThroughLoop==1) {
                // in the first direction, either vertex s or connection are BB/CC. 
                if (currentPolymer[connection].name=="BB") {
                    firstDirection="BB";
                }
                else if (currentPolymer[connection].name=="CC") {
                    firstDirection="CC";
                }
                else if (currentPolymer[s].name=="BB") {
                    firstDirection="BB";
                }
                else if (currentPolymer[s].name=="CC") {
                    firstDirection="CC";
                }
            }

            // print out the vertex descriptor of the monomer connected to the original monomer
            //std::cout<<connection << "\n";
            // note: there was detail on counting  / predicates in the c++ book--I think that stuff is overkill here--but could separate out this stuff into a count function
            if (currentPolymer[connection].name=="AA") {
                // do nothing
            }
            else if (currentPolymer[connection].name=="BB") {
                total_sequence_length_BB+=1;
                // start new sequence of BB
                currentSequence="BB";
                num_seqs_BB+=1;
            }
            else if (currentPolymer[connection].name=="CC") {
                total_sequence_length_CC+=1;
                // start new sequence of CC
                currentSequence="CC";
                num_seqs_CC+=1;
            }
                
            // Now, iteratively search through monomers connected to "connection"
            bool endOfChain=false;
            boost::graph_traits<UndirectedGraph>::vertex_descriptor previousVertex=s;
            boost::graph_traits<UndirectedGraph>::vertex_descriptor currentVertex=s;

            int counter=0;
            while (!endOfChain) {
                previousVertex=currentVertex;
                currentVertex=connection;
                boost::graph_traits<UndirectedGraph>::out_edge_iterator e2, e_end2;
                // loop through edges of the next vertex (now "current") to figure out what the next vertex is in the correct direction
                endOfChain=true; // unless another valid connection is found, the end of the chain has been reached
                for (tie(e2, e_end2) = out_edges(currentVertex, currentPolymer); e2 != e_end2; ++e2) {
                    if (source(*e2, currentPolymer) != previousVertex && target(*e2, currentPolymer) != previousVertex) {
                        if (source(*e2, currentPolymer) != currentVertex) {
                            connection=source(*e2, currentPolymer); //smth funny here
                        }
                        else if (target(*e2, currentPolymer) != currentVertex) {
                            connection=target(*e2, currentPolymer);
                        }
                        else {
                            throw std::runtime_error( "Ruh roh!" );
                        }
                        endOfChain=false;

                        // do counting...
                        if (currentPolymer[connection].name=="AA") {
                            // do nothing
                        }
                        else if (currentPolymer[connection].name=="BB") {
                            total_sequence_length_BB+=1;
                            // if previous sequence was BB, just continue counting!
                            if (currentSequence=="BB") {
                                sequenceLength+=1; // I don't think this sequence length variable actually matters for calculation of NnBB (only NwBB)
                            }
                            // if previous sequence was CC, start new sequence of BB 
                            else if (currentSequence=="CC") {
                                num_seqs_BB+=1;
                                sequenceLength=0; // this only matters for NwBB and NwCC, I think.
                                currentSequence="BB";
                            }        
                        }
                        else if (currentPolymer[connection].name=="CC") {
                            total_sequence_length_CC+=1;
                            if (currentSequence=="CC") {
                                sequenceLength+=1;
                            }
                            else if (currentSequence=="BB") {
                                num_seqs_CC+=1;
                                // start new sequence of CC
                                sequenceLength=0; // this only matters for NwBB and NwCC, I think.
                                currentSequence="CC";
                            }
                        }
                    }
                    // On the second side of the original monomer, after the first time through the inner loop, identify the start of the first sequence on the second side of the polymer
                    if (counter==0 && timesThroughLoop==2) {
                        // first time through, either currentVertex or connection are BB/CC.. previousVertex is s (right?)
                        if (currentPolymer[currentVertex].name=="AA" ) {
                            
                            /* Note: if currentVertex can't be assigned because connection was the 
                            end of chain, it's no big deal because currentVertex==connection.
                            */
                            secondDirection=currentPolymer[connection].name;

                        }
                        else {
                            secondDirection=currentPolymer[currentVertex].name;
                        }
                        // if the same sequence has been recorded on both sides of the original monomer, then num_seqs_BB or num_seqs_CC has been incorrectly tabulated.
                        if (firstDirection==secondDirection) {
                            if (firstDirection=="BB") {
                                num_seqs_BB-=1;
                            }
                            else if (firstDirection=="CC") {
                                num_seqs_CC-=1;
                            }
                            
                        }

                        if ((firstDirection=="AA") || (secondDirection=="AA") || 
                            (firstDirection=="") || (secondDirection=="") ) {
                                // do nothing
                                //std::cout << " end of chain reached before first/second direction could be assigned \n";
                            }
                        // check that first and second direction are both either BB or CC!
                        if (firstDirection!="BB" && firstDirection!="CC") {
                            //std::cout << " end of chain reached before first/second direction could be assigned \n";
                        }
                        if (secondDirection!="CC" && secondDirection!="BB") {
                            //std::cout << " end of chain reached before first/second direction could be assigned \n";
                        }

                    }
                    counter++;
                }
                if (endOfChain) {
                    // if all we're doing is counting BB and CC and the number of sequences (total sequence length and num seqs), nothing needs to happen here
                }

            }

        }
    }
    // output results to file
    if (num_seqs_BB != 0) {
        NnBB=(1.0)*(total_sequence_length_BB+numMonomers["BB"])/(num_seqs_BB+numMonomers["BB"]);
    }
    if (num_seqs_CC != 0) {
        NnCC=((1.0)*total_sequence_length_CC+numMonomers["CC"])/(num_seqs_CC+numMonomers["CC"]);
    }

    // calculate q1 and q2
    double q1 = (1.0*(numMoieties0["B"]-numMoieties["B"]))/(1.0*numMoieties0["B"]);
    double q2 = (1.0*(numMoieties0["C"]-numMoieties["C"]))/(1.0*numMoieties0["C"]);
    //    Tirrell << "time       q1         Nnbb       q2         Nncc\n";
    TirrellOutput << std::left << std::setw(7) << simTime << "     " << std::setw(6) << q1 << "     " << std::setw(6) << NnBB << "     " << std::setw(6) << q2 << "     " << std::setw(6) << NnCC << "  " << Mn << "   "  << Mw << "\n";
}

void printwfs(std::vector<Chain> chainPool,std::string strdate, std::string filename, std::vector<double> monomerWeights,std::map<std::string,signed long> numMonomers,signed long simulationSize) {
    // print out molecular weight distributions
    // assumes all monomers have the same molecular weight

    std::ofstream wf; // weight fractions... note that Driscoll used 1 for all his monomer weights
    //std::ofstream nf; // number fractions
    //std::ofstream wf_actual; // weight fractions using real monomer weight

    wf.open("./Output/"+strdate+"weight_fraction"+filename + ".txt");
    //nf.open(path+date+"number_fraction"+filename);
    signed long wt_fraction_size = simulationSize;

    std::vector<double> wt_fraction(wt_fraction_size,0);
    for (int i=0;i<chainPool.size();i++) {

        if (chainPool[i].weight/(1.0*monomerWeights[0]) >wt_fraction.size() ) {
            throw std::runtime_error("check how the molecular weight distributions are printed to a file");
        }
        wt_fraction[int(chainPool[i].weight/(1.0*monomerWeights[0]))]++; // if all monomers do not have the same molecular weight, need to adjust this line
        
    }
    // Do monomers count as chains of length 1? Check derivation, but I think so
    std::map<std::string,signed long>::iterator it3=numMonomers.begin();
    int totalMonomers=0;
    while (it3 != numMonomers.end()) {
        totalMonomers+=it3->second;
        it3++;
    }
    wt_fraction[1]+=totalMonomers;
    // divide every weight fraction bucket by the total inital mass and weight it by the bucket average mass
    for (int i=0;i<wt_fraction.size();i++) {
        //nf << wt_fraction[i]/simulationSize;
        wt_fraction[i]=wt_fraction[i]*i/(simulationSize*1.0);
        wf << wt_fraction[i] << std::endl;
    }
    wf.close();
    //nf.close();
}

void printEndGrps(const SystemVariables &cSS,std::ofstream &eGOutput, double simTime) {
    // crude code for calculating end group composition and writing to a file
    // assumes there are only three monomers, AA, BB, and CC with end group possibilities A, B, or C...

    // loop through all elements in chainPool and count up number of ends of each type
    // no need to restrict to linear polymers, but as written assumes only A, B, or C end groups...
    std::map<std::string,signed long> endGrpCounts; 
    for (int i =0; i< cSS.chainPool.size(); ++i) {
        // look through each chain and count up number of end groups of each type
        // check number of edges for each monomer... if there is only one, then it is an end group
        auto vxs=boost::vertices(cSS.chainPool[i].polymer);
        while (vxs.first!=vxs.second) {
            int outEdgeCount=0;
            int inEdgeCount=0;
            auto outEdgs=boost::out_edges(*vxs.first, cSS.chainPool[i].polymer);
            auto inEdgs=boost::in_edges(*vxs.first, cSS.chainPool[i].polymer);
            while (outEdgs.first!=outEdgs.second){
                outEdgeCount++;
                outEdgs.first++;
            }
            while (inEdgs.first!=inEdgs.second) {
                inEdgeCount++;
                inEdgs.first++;
            }
            if (outEdgeCount+inEdgeCount==1) {
                // now we have an End GroupTM. You can also check to be sure that we're above the weight cutoff for a polymer, or we can just plan on excluding monomers only.
                // any functional groups on this monomer that aren't edges should get added to the end group analysis... here we are restricting to A, B, and C
                if (cSS.chainPool[i].polymer[*vxs.first].funs.count("A") ) {
                    endGrpCounts["A"]+=cSS.chainPool[i].polymer[*vxs.first].funs.at("A");
                }
                if (cSS.chainPool[i].polymer[*vxs.first].funs.count("B") ) {
                    endGrpCounts["B"]+=cSS.chainPool[i].polymer[*vxs.first].funs.at("B");
                }
                if (cSS.chainPool[i].polymer[*vxs.first].funs.count("C") ) {
                    endGrpCounts["C"]+=cSS.chainPool[i].polymer[*vxs.first].funs.at("C");
                }
            }
            vxs.first++;
        }
    }

    // now we've calculated all end group counts... 
    int totalEndGrps=endGrpCounts["A"]+endGrpCounts["B"]+endGrpCounts["C"];
    double Afr=(endGrpCounts["A"]*1.0)/(1.0*totalEndGrps);
    double Bfr=(endGrpCounts["B"]*1.0)/(1.0*totalEndGrps);
    double Cfr=(endGrpCounts["C"]*1.0)/(1.0*totalEndGrps);

    eGOutput << simTime << "  " << cSS.mwv.Mn << "  " << cSS.mwv.Mw << "  " << Afr << "   " << Bfr << "   "  << Cfr << "   " << std::endl;
    return;
}

void calculateDPrw(const SystemVariables &cSS,std::ofstream &conv) {
    double max = 0;
    signed long long maxChain=0;
    for (int i=0;i<cSS.chainPool.size();++i) {
        if (cSS.chainPool[i].polymer.m_vertices.size()>max) {
            max = cSS.chainPool[i].polymer.m_vertices.size();
            maxChain=i;
        }
    }

    double DPrw=0;
    // now that we know which one is the max, we leave it out of the DPw calculation
    for (int i=0;i<cSS.chainPool.size();++i) {
        if (i!=maxChain) {
            DPrw += cSS.chainPool[i].polymer.m_vertices.size()*cSS.chainPool[i].polymer.m_vertices.size();
        }
    }
    DPrw=DPrw/(1.0*cSS.chainPool.size() ) ;

    conv << DPrw << "   ";
    return;
}

void printExecutions(std::vector<long> reactionsExecutedTally,std::vector<Reaction> reactionList,std::string strdate,std::string runID) {
    // print out the number of times each reaction was executed along with the reaction name
    std::ofstream reactionsExecutedFile;
    reactionsExecutedFile.open("./Output/" + strdate + runID + "reactionsExecuted.txt");
    for (int i=0; i<reactionList.size(); ++i) {
        reactionsExecutedFile << reactionList[i].step << "," << reactionList[i].number << " " << reactionsExecutedTally[i] << std::endl;
    }
    return;
}

    // // print reactions to a file for debugging
    // std::ofstream reactionsFile;
    // reactionsFile.open("./Output/" + strdate + runID + "reactions.txt");
    // // each line corresponds to a reaction
    // for (int i=0; i<reactionList.size(); ++i) {
    //     reactionsFile << "Reaction type: " << reactionList[i].reactionType << "  " << " reactants: ";
    //     for (int j=0; j<reactionList[i].reactants.size(); ++j) {
    //         reactionsFile << reactionList[i].reactants[j] << " ";
    //     }
    //     reactionsFile << " linkToBust: " << reactionList[i].linkToBust;
    //     reactionsFile << " linkFormed: " << reactionList[i].linkFormed;
    //     reactionsFile << " gasEvolved: " << reactionList[i].gasEvolved;
    //     reactionsFile << " catalystEffect: " << reactionList[i].catalystEffect;
    //     reactionsFile << " A: " << reactionList[i].A;
    //     reactionsFile << " products: ";
    //     for (int j=0; j<reactionList[i].products.size(); ++j) {
    //         reactionsFile << reactionList[i].products[j] << " ";
    //     }
    //     reactionsFile << std::endl;
    // }
    // reactionsFile.close();


            // // calculate what Mn *should* be
        // double totalWeight=0;
        // int chainsPastCutoff=0;
        // for (int i=0; i<cSS.chainPool.size(); ++i) {
        //     if (cSS.chainPool[i].weight > cSS.mwv.weightCutoff ) {
        //         chainsPastCutoff++;
        //         totalWeight+=cSS.chainPool[i].weight;
        //     }
        // }
        // double calculatedMn=totalWeight/(1.0*chainsPastCutoff);

        // if ( (calculatedMn-cSS.mwv.Mn) *(calculatedMn-cSS.mwv.Mn) > 0.1  ) {
        //     throw std::runtime_error("whatttt");
        // }

        // // check that the amount of catalyst always equal the initial amount of catalyst
        // int totalCat=0;
        // int linkCat=0;
        // totalCat+=cSS.numMoieties["*"];
        // for (auto it=cSS.links.begin(); it!=cSS.links.end(); ++it) {
        //     if ((it->first).find('*') != std::string::npos) {
        //         linkCat+=it->second;
        //     }
        // }
        // totalCat+=linkCat;
        // // ignore the below code... it will also always occur in a link name. actually, this checker is kinda
        // // buggy as originally written because a * in a vertex name may appear in more than one link name...
        // // so my test is no good. I might try cycling through all *edges* and ensuring that the sum of them
        // // plus vertices have 25-numMoieties * units...
        // // // catalyst can also be bound to a chain unit..
        // // totalCat+=cSS.numMoieties["E1*"];
        // // totalCat+=cSS.numMoieties["Em2*"];
        // // totalCat+=cSS.numMoieties["E3*"];
        // // totalCat+=cSS.numMoieties["Ep4*"];
        // if (totalCat!=25) {
        //     throw std::runtime_error("woah woah"); 
        // }
        // if (cSS.executionNumber%200==0) {
        //     countLinks(cSS);
        // }
        // // check that all moieties are positive
        // for (auto it=cSS.numMoieties.begin(); it!=cSS.numMoieties.end(); ++it) {
        //     if (it->second<0) {
        //         std::cout << "group " << it->first << " is negative!! Something is wrong" << "\n";
        //         std::cout << "execution number was " << cSS.executionNumber << "\n";
        //         throw std::runtime_error("watttt");
        //     }
        // }

std::vector<std::string> enumerateLinks(std::vector<Reaction> reactionList, std::ofstream &linksFile) {
    // this function returns all possible links that could be made with a given reactionList
    // it also prints the link names as headers to a file that will be used to print the link populations over time
    // the file is called "links.txt"
    linksFile << "time ";
    std::vector<std::string> allPossibleLinks;
    // loop through the reaction list and determine the product link name.
    // if the product link name is not already in allPossibleLinks, add it
    for (int i=0; i<reactionList.size(); ++i) {
        // if there is only one product, then never mind. this reaction doesn't make a link
        if (reactionList[i].products.size()<2) {
            continue;
        }
        else if (reactionList[i].products.size()==2) {
            if (reactionList[i].linkFormed!="") {
                std::string productLinkName=reactionList[i].products[0]+"-" +reactionList[i].linkFormed+"-" +reactionList[i].products[1];

                // check if productLinkName is in allPossibleLinks
                //bool alreadyInList=false;
                
                //if (!alreadyInList) {
                allPossibleLinks.push_back(productLinkName);
                linksFile << productLinkName << " ";
                //}
            }
            
        }
        else {
            throw std::runtime_error("reaction has more than two products. this is not allowed");
        }
    }
    // These are all the links that can depolymerize... now let's write to a file the links that cannot break
    allPossibleLinks.push_back("E-ep-cE");
    allPossibleLinks.push_back("Ep-ep-cE");
    allPossibleLinks.push_back("E-ep-Em");
    allPossibleLinks.push_back("Ep-ep-Em");
    allPossibleLinks.push_back("E-em-Ep");
    allPossibleLinks.push_back("Ep-em-Ep");
    allPossibleLinks.push_back("E-em-E");
    allPossibleLinks.push_back("Ep-em-E");
    allPossibleLinks.push_back("E-ep-E");
    allPossibleLinks.push_back("Ep-ep-E");

    // for (int j=0; j<allPossibleLinks.size(); ++j) {
    //     if (allPossibleLinks[j]==productLinkName) {
    //         std::cout << "unexpectd behavior in enumerateLinks" << std::endl;
    //         break;
    //     }
    // }

    auto result = allPossibleLinks; // make a copy before we start sorting the vector to check for duplicates
    // find any duplicate strings in allPossibleLinks
    std::sort(allPossibleLinks.begin(), allPossibleLinks.end());
    // if there is a duplicate string, print out a message to standard output
    for (int i=0; i<allPossibleLinks.size()-1; ++i) {
        if (allPossibleLinks[i]==allPossibleLinks[i+1]) {
            std::cout << "duplicate string in allPossibleLinks!! shouldnt happen" << std::endl;
        }
    }

    linksFile << std::endl;
    return result;
}

void printAllLinksToFile(std::map<std::string, signed long> links, std::ofstream &linksFile, std::vector<std::string> allPossibleLinks,double simTime,double V, double const Na) {
    // this function should print out all links in cSS.links to a file along with the simulation time
    // if a link is not in cSS.links, then print a 0

    linksFile << simTime << " ";
    for (int i=0; i<allPossibleLinks.size(); ++i) {
        if ( links.count(allPossibleLinks[i]) ) {
            linksFile << (links[allPossibleLinks[i]]*1.0)/(V*Na) << " ";
        }
        else {
            linksFile << 0 << " ";
        }
    }
    linksFile << std::endl;
    return;
}
