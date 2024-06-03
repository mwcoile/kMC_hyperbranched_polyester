#include"molecularWeight.h"
#include"chainUpdate.h"
#include"inputs.h"
#include"tests.h"
#include"analysis.h"


#if __cplusplus < 201103L
#error This file requires compiler and library support for the \
ISO C++ 2011 standard. This support is currently experimental, and must be \
enabled with the -std=c++11 or -std=gnu++11 compiler options.
#endif

/* 
author: Matthew Coile
purpose: kinetic model capable of tracking polymer structures during (de)polymerization
date: conceived in 2019, first working draft completed Nov 24, 2021
*/

double chainLengthDependentk(int i, int j, double k_int) {
    // i is the length of chain 1
    // j is the length of chain 2
    // k_int is the intrinsic reaction rate 
    double phi = 0.2; // 3.5E-9
    double l = 1.0;
    double Ha_ij = phi/( (sqrt(i)+sqrt(j))* (pow(i,-l)+pow(j,-l)));
    double Kij = k_int / (1.0 + Ha_ij);
    return Kij;
}

double approxChainLengthDependentk(int i) {
    // i is the length of the shorter of two chains participating in a reaction 
    // k_int is the intrinsic reaction rate 
    double phi = 1.4*4; // 3.5E-9
    double l = 1.0;
    
    double Ha_ij = phi*pow(i,l-0.5)/4.0;
    double Kij_divided_by_k_int = 1; // / (1.0 + Ha_ij);
    return Kij_divided_by_k_int; //1.0
}

chosenReactionChannel selectReactionChannel(std::vector<double> rate,std::vector<double> loopingRate, std::vector<std::map<std::string,double>> ratesByLength, double totalRate,double r2,std::vector<Reaction> reactionList) {
    // this function takes as input a vector of rates and the sum of these rates (totalRate) and chooses one of the 
    // entries in the vector of rates randomly based on their relative values. For example, for two reactions with 
    // rates = {1, 2}, an index corresponding to the second entry is twice as likely to be returned as an index 
    // corresponding to the first entry
    chosenReactionChannel chosenRxn;
    chosenRxn.isLengthDependent=false;
    int mu=0; // reaction channel selector
    double sumRates=0; // Ref: DOI: 10.1002/mats.201000041 equation 1, here implemented with all terms multiplied by total rate
    while (sumRates<r2*totalRate) { 
        if (mu<rate.size()) {
            if (mu!=4 && mu!=5 && mu!=6 && mu !=7) {
		    sumRates+=rate[mu];
	    }
	    else {
		std::string R_i = "";
	    	if (mu==4) R_i = "R5";
		else if (mu==5) R_i = "R6";
		else if (mu==6) R_i = "R7";
		else if (mu==7) R_i = "R8";
		else throw std::runtime_error("unforseen isssue!");
		for (int i=1;i<ratesByLength.size();++i) {
			sumRates+=ratesByLength[i][R_i];
			if (sumRates>=r2*totalRate) {
				chosenRxn.isLengthDependent=true;
				chosenRxn.length=i; // length of the chain containing the second moiety
				chosenRxn.chosenR_i=R_i;
				goto alldone;	
			}
		}
			
	    }
	} 
//        else {
//            // figure out which length we are corresponding to...
//            for (int i=1;i<ratesByLength.size();++i) {
//                // figure out which specific reaction we are corresponding to
//                for (auto it=ratesByLength[i].begin(); it!=ratesByLength[i].end(); ++it) {
//                    
//                    sumRates+=it->second;
//                    if (sumRates>=r2*totalRate) {
//                        chosenRxn.isLengthDependent=true;
//                        chosenRxn.length=i; // this is the length of the first or second moiety?
//                        chosenRxn.chosenR_i=it->first;
//                        if (chosenRxn.chosenR_i=="R5") {
//                            mu=4;
//                        }
//                        else if (chosenRxn.chosenR_i=="R6") {
//                            mu=5;
//                        }
//                        else if (chosenRxn.chosenR_i=="R7") {
//                            mu=6;
//                        }
//                        else if (chosenRxn.chosenR_i=="R8") {
//                            mu=7;
//                        }
//			else {
//				throw std::runtime_error("we have an unforeseen issue");
//			}
//                        //chosenRxn.reactionType[]// here set the reaction type
//                        goto alldone;
//                    }
//                }
//            }
//        }
        // else {
        //     sumRates+=loopingRate[mu-rate.size()];
        // }
	if (mu>=rate.size()) {
		throw std::runtime_error("we have a pgoblem!");
	}
        if (sumRates<r2*totalRate) mu++; 
    }
alldone:
    chosenRxn.mu=mu;
    if (mu<rate.size()) {
        chosenRxn.reactionType=reactionList[mu].reactionType;
    }
    else {
        if (chosenRxn.chosenR_i=="R5") chosenRxn.reactionType=reactionList[4].reactionType;
        else if (chosenRxn.chosenR_i=="R6") chosenRxn.reactionType=reactionList[5].reactionType;
        else if (chosenRxn.chosenR_i=="R7") chosenRxn.reactionType=reactionList[6].reactionType;
        else if (chosenRxn.chosenR_i=="R8") chosenRxn.reactionType=reactionList[7].reactionType;
	else {
		throw std::runtime_error("this should not happen !");
	}
	if (chosenRxn.reactionType!="linking") {
		throw std::runtime_error("Ah found a bug!");
	}
    }
    if (mu>=rate.size() && chosenRxn.isLengthDependent==false) {
        throw std::runtime_error("oh no did not find the reaction to execute!");
    }
    return chosenRxn; // return chosen reaction
}


int main() {

    const bool timeShit = true;

    TimePoint start;
    TimePoint startRandAndSetup;
    TimePoint stopRandAndSetup;
    TimePoint startKMC;
    TimePoint stopKMC;
    TimePoint startPropensityTimer;
    TimePoint startKMCInnerLoop;
    TimePoint stopPropensityTimer; 
    TimePoint startRxnExecution;
    TimePoint stopRxnExecution;
    TimePoint startmolwt;
    TimePoint stopmolwt;
    TimePoint startanalysis;
    TimePoint stopanalysis;

    TimePoint stopKMCInnerLoop; //get the timepoint at this instant 
    TimePoint startprint;
    
    TimePoint stopprint;
    TimePoint stop; 
    DoubleSec innerLoopTime;
    start = Time::now(); //get the timepoint at this instant 
    startRandAndSetup = Time::now(); //get the timepoint at this instant 

    std::map<std::string, DoubleSec> timers;

    if (timeShit) {
        DoubleSec zero = DoubleSec::zero(); // (Time::now()-Time::now());
        timers["analysisTime"]=zero;
        timers["molwt"]=zero;
        timers["linking"]= zero;
        timers["breaking"]=zero;
        timers["carbamateDecomposition"]=zero;
        timers["allophanate"]=zero;
        timers["propensityAndChannelSelection"]=zero;
        timers["linearSearchForParticipatingGroup"]=zero;
        timers["monomer-monomerLinking"]=zero;
        timers["monomer-chainLinking"]=zero;
        timers["chain-chainLinking"]=zero;
        // these are the ones to check Oct 9, 2023 for whether improving the split function is worth it
        timers["split"]=zero;
        timers["connectedcomponents"]=zero;
        timers["adjbreaking"]=zero;
        timers["calcgrpsByLength"]=zero;
        timers["calccns"]=zero;
        timers["calcratesByLength"]=zero;

    }
    

    // trackLoopPairs, trackFullStructureAsGraph; Tirrell; Changxia; printWeightFractions; eG, trackAllLinks, loopingSeparateChannel
    ModelOptions myOptions= {false,  // trackLoopPairs
                            true, // trackFullStructureAsGraph
                            false, // Tirrell
                            true, // Changxia
                            false, // printWeightFractions
                            false, // end group analysis
                            true, // trackAllLinks
                            true, // loopingSeparateChannel: treat looping as a separate reaction channel?
                            false}; // limitStep3
    if (myOptions.trackLoopPairs==true && myOptions.loopingSeparateChannel==false) {
        throw std::runtime_error( "inconsistent options - no point in tracking looping pairs if we aren't considering looping as a separate reaction channel");
    }

    // random number generation, slightly more sophisticated than using std::rand
    // based on example here https://stackoverflow.com/questions/288739/generate-random-numbers-uniformly-over-an-entire-range
    std::random_device rand_dev; // obtain unsigned int seed from a random device (e.g., time(0)) by calling rand_dev()
    unsigned int seed = rand_dev(); //1797740251 // currently busted with seed 961071223 // try 1386656568 //works with 1124697860; // for debugging, can set this to a constant value if desired, e.g., 1086615307; 
    std::mt19937 generator(seed); // generator(rand_dev()); 
    //std::cout << seed << std::endl; // write this out to the standard output...
    const double range_from = 0; // min value in range (inclusive)
    const double range_to = 1; // max value of range (not inclusive)
    std::uniform_real_distribution<double> distr(range_from, range_to); // in range [from, to)
    // example of how to generate a number in the range: "distr(generator)"




    // READ IN USER INPUTS FROM FILE //
    InitialConditions myConditions = readInitialConditionsFromInput();
    SystemVariables cSS = myConditions.cSS0;
    cSS.seed = seed;
    double T = myConditions.T; // temperature in K
    std::vector<Reaction> reactionList = myConditions.myReactions;
    
    std::string runID="seed_is_" + std::to_string(seed) + "simulation_size_is" + std::to_string(myConditions.simulationSize) + "dmax_is" + std::to_string(cSS.dmax); // a particular name can be appended to all filenames for a run, do so here

    // begin all output filenames with current date and time. These lines obtain the current date and time as a string suitable for this purpose
    time_t t = time(0);   // get current time
    struct tm * now = localtime( & t );
    char date[80]; // initialize place to store date
    strftime (date,80,"%Y%m%d_%H%M%S_%Z",now); // https://www.cplusplus.com/reference/ctime/strftime/ 
    std::string strdate=std::string(date);
    std::ofstream TirrellOutput;
    if (myOptions.Tirrell) {
        TirrellOutput.open("./Output/" +strdate + runID + "TirrellOutput.txt");
    }
    std::ofstream ChangxiaOutput;
    ChangxiaOutput.open("./Output/" + strdate + runID + "ChangxiaOutput.txt");

    std::ofstream eG;
    if (myOptions.eG) {
        eG.open("./Output/"+strdate+"endgrpfrs" + ".txt");
    }
    //std::ofstream // statusFile;
    // statusFile.open("./Output/"+runID+"statusFile.txt");
    // calculate simulation volume
    const double Na = 6.022E23; // avogadro's number: items per mole
    double V = (1.0*myConditions.simulationSize)/(Na*myConditions.totalMonomerConcentration); // control volume in L assuming concentration is in mol/L 
    double Ct = myConditions.totalMonomerConcentration / 500.0; // total concentration of catalyst in mol/L
    const double R = 0.008314; // ideal gas constant, kJ/(K mol)
    double simTime=0;
    
    std::vector<double> rate(reactionList.size(),0); // declare vector to store rates corresponding to each reaction
    double over_x=0; // conversion -- definition dependent upon system
    double totalRate=1; // initialize to nonzero value so that kMC loop will execute
    
    readInLHHWLinking(reactionList);
    readInLHHWBreaking(reactionList);

    if (timeShit) {
        stopRandAndSetup = Time::now(); //get the timepoint at this instant 
        timers["RandAndSetup"]=(stopRandAndSetup-startRandAndSetup);
        startKMC = Time::now(); //get the timepoint at this instant 
        innerLoopTime = startKMC - startKMC; // zero it out, but I don't feel like looking up the type
    }

    // initialize cSS.reactionsExecutedTally
    cSS.reactionsExecutedTally.resize(reactionList.size(),0);
    
    // open the output file that will print the number of each link over time
    std::ofstream linksFile;
    linksFile.open("./Output/" + strdate + runID + "links.txt");
    // open the output file that will print out the net rate analysis
    std::ofstream netRateAnalysis;
    netRateAnalysis.open("./Output/" + strdate + runID + "netRateAnalysis.txt");
    
    // determine all possible links that may form 
    std::vector<std::string> allPossibleLinks = enumerateLinks(reactionList, linksFile);


    // assume we have a known set of equilibrium constants
    std::map<std::string,std::map<std::string,double>> capK=getCapK(strdate,runID); // equilibrium constants for each reaction
    // also assume we have a known set of rate constants. Both of these will be read in from a file
    std::map<std::string,std::map<std::string,double>> rateConstants=getRateConstants(); // rate constants for each reaction

    // initialize transport_factor
    std::vector<double> transport_factor(myConditions.simulationSize+1,0);
    // for each length, calculate transport factor for each reaction
    for (long i=1; i<transport_factor.size(); ++i) {
        // for each of the 10 forward reactions -- note that we could also just do this for the step growth reactions R5-R8
        //for (int j=4;j<8;++j) {
            //std::string R_i = reactionList[j].number;
        transport_factor[i] += approxChainLengthDependentk(i);
        //}
    }
    // calculate all factors as function of chain length...
    
    // initialize numMoietiesByLength
    //std::vector<std::map<std::string,signed long>> numMoietiesByLength(myConditions.simulationSize+1,std::map<std::string,signed long>());
    //cSS.numMoietiesByLength = numMoietiesByLength;
    
    std::vector<std::map<std::string,signed long>> numMoietiesByLength(myConditions.simulationSize+1,std::map<std::string,signed long>());
        //for (int i=0;i<cSS.numMoietiesByLength;++i) {
    cSS.numMoietiesByLength = numMoietiesByLength;
    std::ofstream coverages;
    coverages.open("./Output/" + strdate + runID + "coverages.txt");
    coverages << "time  X    3R1*	2R2*	3R3*	M*	1R*	2R*	3R*	4R*	5R*	6R*	7R*	8R*	9R*	10R*	C1*	C2*	C3*	C4*	C5*	C6*	C7*	C8*	C9*	C10*	1I*	2I*	3I*	4I*	5I*	6I*	7I*	8I*	9I*	10I*	1P*	2P*	3P*	4P*	5P*	6P*	7P*	8P*	9P*	10P*	1P2*	2P2*	3P2*	4P2*	5P2*	6P2*	7P2*	8P2*	9P2*	10P2* \n";
    // coverages << "time  X  * 5-*_step1   2-*_step1 4-*_step1 M-*_step1 M-1-*_step1   M-*-M_competitive   2-b-3_step3 " <<
    //     "2-b-4-*_step3  2-a-3-*_step3   2-a-5-*_step3 2-b-3-*_step3 2-b-4-*_step3 4-b-4-*_step3 4-a-3-*_step3 4-b-1-*_step3 4-b-1-*_step3" <<
    //     "2-b-3-*_step4    2-b-4-*_step4  2-a-3-*_step4   2-a-5-*_step4  2-b-3-*_step4 2-b-4-*_step4 4-b-4-*_step4 4-a-3-*_step4 4-b-1-*_step4 4-b-1-*_step4 \n";
    std::map<std::string,bool> didIprint;
    // kMC loop!
    while (totalRate>0 && simTime<60*60*24*2 && over_x<0.95) {
        // statusFile << "start of kMC loop " << cSS.executionNumber << std::endl;
        // statusFile << std::flush;

        if (cSS.rejectedSteps>20) {
            std::cout << " breaking out of kMC loop because too many rejected steps" << std::endl;
            break;
        }
        if (cSS.rejectedSteps>0) {
            //std::cout << "rejected steps: " << cSS.rejectedSteps << std::endl;
        }

        if (timeShit) {
            startPropensityTimer = Time::now(); // start clock
            startKMCInnerLoop = Time::now(); //get the timepoint at this instant 
        }
        

        cSS.executionNumber++;

        if (cSS.executionNumber%200==0) {
            std::cout << "execution number: " << cSS.executionNumber << " sim time: " <<  simTime << " conversion: " << over_x << " Mn  "  <<  cSS.mwv.Mn << "  Mw  " << cSS.mwv.Mw << " chainPool size: "  <<  cSS.chainPool.size() << "\n";
            //std::cout << "count of Em6-ep*-Em6: " << cSS.links["Em6-ep*-Em6"] << "  count of Emp2-ep*-Em2  " << cSS.links["Emp2-ep*-Em2"] << " count of Emp9-[eplinkNotFormed]*-cE9  " << cSS.links["Emp9-[eplinkNotFormed]*-cE9"] << "  count of Emp9-ep*-cE9 "<< cSS.links["Emp9-ep*-cE9"] << "\n";
            //std::cout << "cat binding tally  " << cSS.catBindingTally << "\n";
        }
       

        // can run tests here
        //testCode(cSS);
        //test4(cSS,myOptions);
        // i think we need to write a test function here checking that links are being updated appropriately.
        // test3(cSS,myOptions);

        // // Enumerate possible looping reactions
        // std::vector<PossibleLoopingRxn> candidateLoopingReactions = enumeratePossibleLoops(cSS);

        // calculate propensity functions

        // right at the start, must calculate "normalizedCumSum." This requires knowledge
        // of "numMoietiesByLength." In the future, we will calculate update numMoietiesByLength
        // on a reaction-by-reaction basis, but for now we will update numMoietiesByLength 
        // right at the beginning
        // initialize empty vector of maps
        TimePoint startcns= Time::now(); // start clock

        // cSS.numMoietiesByLength = numMoietiesByLength;
        // for (int i=0;i<cSS.chainPool.size(); ++i) {
        //   // numMoieties on monomers doesn't matter for this... only consider chains
        //     int chainSize=boost::num_vertices(cSS.chainPool[i].polymer);
        //     for (auto it = cSS.chainPool[i].funs.begin(); it!=cSS.chainPool[i].funs.end(); ++it) {
        //         cSS.numMoietiesByLength[chainSize][it->first]+=it->second; //cSS.chainPool[i].funs[it->first]
        //     }
        // }
        //std::vector<std::map<std::string,long>> regCumSum = {cSS.maxChainSize+1,std::map<std::string,long>()};
        // create a vector of maps of type std::vector<std::map<std::string,long>> regCumSum with the length of the vector being cSS.maxChainSize+1
        std::vector<std::map<std::string,long>> regCumSum(cSS.maxChainSize+1,std::map<std::string,long>());

        // initialize normalizedCumSum with empty maps of the length of the simulation size (?does this actually need to happen)
        int myCount =0;

        //std::vector<std::map<std::string,double>> normalizedCumSum = {cSS.maxChainSize+1,std::map<std::string,double>()};

        // create a vector of maps of type std::vector<std::map<std::string,double>> normalizedCumSu with the length of the vector being cSS.maxChainSize+1
        std::vector<std::map<std::string,double>> normalizedCumSum(cSS.maxChainSize+1,std::map<std::string,double>());

        for (long i=1; i<cSS.maxChainSize+1;++i) {
            // sum up everything less than i multiplied by the corresponding diffusion factor
            for (long j=1; j<=(i-1); ++j) {
                // do this for each moiety
                for (auto it = cSS.numMoietiesByLength[j].begin(); it != cSS.numMoietiesByLength[j].end(); ++it) {
                    normalizedCumSum[i][it->first]+=it->second*transport_factor[j];
                    regCumSum[i][it->first]+=it->second;
                    //myCount++;
                    if (it->second==0) {
                        throw std::runtime_error("help");
                    }
                }
            }
        }
        //std::cout << myCount << std::endl;
        TimePoint stopcns = Time::now(); // stop clock
        timers["calccns"] +=(stopcns-startcns);

        // // just check that numMoietiesByLength always sums up to numMoieties
        // std::map<std::string,int> confirmCountsMatch;
        // for (auto iter = cSS.numMoieties.begin(); iter!=cSS.numMoieties.end(); ++iter) {
        //     
        //     for (int i=0;i<cSS.numMoietiesByLength.size();++i) {
        //         if (cSS.numMoietiesByLength[i].find(iter->first)!=cSS.numMoietiesByLength[i].end()) {
        //             confirmCountsMatch[iter->first]+=cSS.numMoietiesByLength[i].find(iter->first)->second;
        //         }
        //     }
        // }
        // for (auto iter = confirmCountsMatch.begin();iter!=confirmCountsMatch.end(); ++iter) {
        //     if (iter->second!=cSS.numMoieties[iter->first] && iter->first!="cEp") {
        //         throw std::runtime_error("something fishy here");
        //     }
        // }
        // now check that numMoietiesByLength
        // calculate free catalyst concentration...
        std::string P1 = "Emp-ep-E";
        std::string P2 = "Emp-ep-Em";
        std::string P3 = "Emp-em-E";
        std::string P4 = "Emp-em-Ep";
        std::string P5 = "Em-ep-E";
        std::string P6 = "Em-ep-Em";
        std::string P7 = "Em-em-E";
        std::string P8 = "Em-em-Ep";
        std::string P9 = "Emp-ep-cE";
        std::string P10 = "Em-ep-cE";

        // note, Ct is total catalyst concentration in mol/L
        // all numMoieties terms should be converted to mol/L
        double forwardBindingSites = capK["R1"]["Kstep1"]*cSS.numMoieties["Ep"]*1.0/(V*Na)
        +capK["R2"]["Kstep1"]*cSS.numMoieties["Emp"]*1.0/(V*Na)
        +capK["R3"]["Kstep1"]*cSS.numMoieties["Em"]*1.0/(V*Na)
        +capK["R9"]["Kstep1"]*cSS.numMoieties["cEp"]*1.0/(V*Na)
        +capK["R1"]["Kstep1"]*capK["R1"]["Kstep2"]*cSS.numMoieties["cEp"]*cSS.numMoieties["Ep"]*1.0/(V*Na)*1.0/(V*Na)
        +capK["R2"]["Kstep1"]*capK["R2"]["Kstep2"]*cSS.numMoieties["cEp"]*cSS.numMoieties["Emp"]*1.0/(V*Na)*1.0/(V*Na)
        +capK["R3"]["Kstep1"]*capK["R3"]["Kstep2"]*cSS.numMoieties["cEp"]*cSS.numMoieties["Em"]*1.0/(V*Na)*1.0/(V*Na)
        +capK["R4"]["Kstep1"]*capK["R4"]["Kstep2"]*cSS.numMoieties["cEp"]*cSS.numMoieties["Emp"]*1.0/(V*Na)*1.0/(V*Na)
        +capK["R5"]["Kstep1"]*capK["R5"]["Kstep2"]*cSS.numMoieties["cE"]*cSS.numMoieties["Ep"]*1.0/(V*Na)*1.0/(V*Na)
        +capK["R6"]["Kstep1"]*capK["R6"]["Kstep2"]*cSS.numMoieties["cE"]*cSS.numMoieties["Emp"]*1.0/(V*Na)*1.0/(V*Na)
        +capK["R7"]["Kstep1"]*capK["R7"]["Kstep2"]*cSS.numMoieties["cE"]*cSS.numMoieties["Em"]*1.0/(V*Na)*1.0/(V*Na)
        +capK["R8"]["Kstep1"]*capK["R8"]["Kstep2"]*cSS.numMoieties["cE"]*cSS.numMoieties["Emp"]*1.0/(V*Na)*1.0/(V*Na)
        +capK["R9"]["Kstep1"]*capK["R9"]["Kstep2"]*cSS.numMoieties["cEp"]*cSS.numMoieties["cEp"]*1.0/(V*Na)*1.0/(V*Na)
        +capK["R10"]["Kstep1"]*capK["R10"]["Kstep2"]*cSS.numMoieties["cE"]*cSS.numMoieties["cEp"]*1.0/(V*Na)*1.0/(V*Na);

        double otherBindingSites = capK["R9"]["Kbind"]*cSS.numMoieties["cEp"]*cSS.numMoieties["cEp"]*1.0/(V*Na)*1.0/(V*Na)
        +capK["R9"]["Kbind"]*cSS.numMoieties["cEp"]*cSS.numMoieties["Ep"]*1.0/(V*Na)*1.0/(V*Na)
        +capK["R9"]["Kbind"]*cSS.numMoieties["cEp"]*cSS.numMoieties["Em"]*1.0/(V*Na)*1.0/(V*Na)
        +2*capK["R9"]["Kbind"]*cSS.numMoieties["cEp"]*cSS.numMoieties["Emp"]*1.0/(V*Na)*1.0/(V*Na)
        +4*capK["R9"]["Kbind"]*cSS.numMoieties["Emp"]*cSS.numMoieties["Emp"]*1.0/(V*Na)*1.0/(V*Na)
        +2*capK["R9"]["Kbind"]*cSS.numMoieties["Emp"]*cSS.numMoieties["Em"]*1.0/(V*Na)*1.0/(V*Na)
        +2*capK["R9"]["Kbind"]*cSS.numMoieties["Emp"]*cSS.numMoieties["Ep"]*1.0/(V*Na)*1.0/(V*Na)
        +capK["R9"]["Kbind"]*cSS.numMoieties["Ep"]*cSS.numMoieties["Ep"]*1.0/(V*Na)*1.0/(V*Na)
        +capK["R9"]["Kbind"]*cSS.numMoieties["Ep"]*cSS.numMoieties["Em"]*1.0/(V*Na)*1.0/(V*Na)
        +capK["R9"]["Kbind"]*cSS.numMoieties["Em"]*cSS.numMoieties["Em"]*1.0/(V*Na)*1.0/(V*Na);

        double reverseBindingSites = capK["R1"]["Kprimestep5"]*capK["R1"]["Kprimestep6"]*cSS.links[P1]*1.0/(V*Na)
        +capK["R2"]["Kprimestep5"]*capK["R2"]["Kprimestep6"]*cSS.links[P2]*1.0/(V*Na)
        +capK["R3"]["Kprimestep5"]*capK["R3"]["Kprimestep6"]*cSS.links[P3]*1.0/(V*Na)
        +capK["R4"]["Kprimestep5"]*capK["R4"]["Kprimestep6"]*cSS.links[P4]*1.0/(V*Na)
        +capK["R5"]["Kprimestep5"]*capK["R5"]["Kprimestep6"]*cSS.links[P5]*1.0/(V*Na)
        +capK["R6"]["Kprimestep5"]*capK["R6"]["Kprimestep6"]*cSS.links[P6]*1.0/(V*Na)
        +capK["R7"]["Kprimestep5"]*capK["R7"]["Kprimestep6"]*cSS.links[P7]*1.0/(V*Na)
        +capK["R8"]["Kprimestep5"]*capK["R8"]["Kprimestep6"]*cSS.links[P8]*1.0/(V*Na)
        +capK["R9"]["Kprimestep5"]*capK["R9"]["Kprimestep6"]*cSS.links[P9]*1.0/(V*Na)
        +capK["R10"]["Kprimestep5"]*capK["R10"]["Kprimestep6"]*cSS.links[P10]*1.0/(V*Na)
        +capK["R1"]["Kprimestep4"]*capK["R1"]["Kprimestep5"]*capK["R1"]["Kprimestep6"]*cSS.links[P1]*1.0/(V*Na)
        +capK["R2"]["Kprimestep4"]*capK["R2"]["Kprimestep5"]*capK["R2"]["Kprimestep6"]*cSS.links[P2]*1.0/(V*Na)
        +capK["R3"]["Kprimestep4"]*capK["R3"]["Kprimestep5"]*capK["R3"]["Kprimestep6"]*cSS.links[P3]*1.0/(V*Na)
        +capK["R4"]["Kprimestep4"]*capK["R4"]["Kprimestep5"]*capK["R4"]["Kprimestep6"]*cSS.links[P4]*1.0/(V*Na)
        +capK["R5"]["Kprimestep4"]*capK["R5"]["Kprimestep5"]*capK["R5"]["Kprimestep6"]*cSS.links[P5]*1.0/(V*Na)
        +capK["R6"]["Kprimestep4"]*capK["R6"]["Kprimestep5"]*capK["R6"]["Kprimestep6"]*cSS.links[P6]*1.0/(V*Na)
        +capK["R7"]["Kprimestep4"]*capK["R7"]["Kprimestep5"]*capK["R7"]["Kprimestep6"]*cSS.links[P7]*1.0/(V*Na)
        +capK["R8"]["Kprimestep4"]*capK["R8"]["Kprimestep5"]*capK["R8"]["Kprimestep6"]*cSS.links[P8]*1.0/(V*Na)
        +capK["R9"]["Kprimestep4"]*capK["R9"]["Kprimestep5"]*capK["R9"]["Kprimestep6"]*cSS.links[P9]*1.0/(V*Na)
        +capK["R10"]["Kprimestep4"]*capK["R10"]["Kprimestep5"]*capK["R10"]["Kprimestep6"]*cSS.links[P10]*1.0/(V*Na);



        double star = Ct / (1.0 + forwardBindingSites + otherBindingSites + reverseBindingSites); 
	double denom = (1.0 + forwardBindingSites + otherBindingSites + reverseBindingSites);
        std::vector<double> ksteptimesstar={0.0,0.0,0.0,0.0};
        totalRate=0;
        std::vector<double> rate(reactionList.size(),0); // keep track of each non-looping reaction channel's rate
        //std::vector<double> loopingRates(candidateLoopingReactions.size(),0); // keep track of each looping reaction's rate

        bool allOneChain=false;

        std::vector<std::map<std::string,double>> ratesByLength(myConditions.simulationSize+1,std::map<std::string,double>()); // initialize vector of maps to store rates by chain length
	std::vector<double> ratesDebugging(20,0);
        std::map<std::string,double> lengthDependentRateSums;
        // I think it's easier to just entirely rewrite this for v5.
        // statusFile << "start of propensity calculations " << cSS.executionNumber << std::endl;
        // statusFile << std::flush;
        for (int i=0; i<reactionList.size(); ++i) {
            if (reactionList[i].reactants.size()!=2 || reactionList[i].products.size()!=2) {
                throw std::runtime_error("2 reactants and 2 products should be specified for all reactions in Changxia's system.");
            }
            if (reactionList[i].reactionType=="linking") {
                std::string R_i = reactionList[i].number;
                if (R_i == "R1" || R_i == "R2" || R_i == "R3" || R_i == "R4" || R_i == "R9" || R_i == "R10") {
                    double kf=rateConstants[R_i]["kf"];
                    double kr=rateConstants[R_i]["kr"];
                    double kfprime = capK[R_i]["Kstep1"]*capK[R_i]["Kstep2"]*kf;
                    double cfprime = kfprime/(V*Na); // kf is trimolecular for these linking reactions: catalyst+2 reactants
                    double krprime = kr * capK[R_i]["Kprimestep4"]*capK[R_i]["Kprimestep5"]*capK[R_i]["Kprimestep6"];
                    double crprime = krprime; // kr is bimolecular (catalyst + link)
                    if (cSS.chainPool.size()==1 && cSS.numMonomers["cEp"]==0) {
                        // no linking reactions possible... all linking reactions should be carried out as looping reactions
                        rate[i]=0;
                        allOneChain=true;
                        std::cout << "all one chain" << std::endl;
                        continue;
                    }
                    // calculate c parameter -- linking reaction should always involve 2 reactants and a catalyst molecule
                    if (reactionList[i].reactants.size()==2) {
                        std::string moiety1=reactionList[i].reactants[0];
                        std::string moiety2=reactionList[i].reactants[1];

                        // check for the special case for which there is ONE monomer remaining -- monomer+monomer reaction cannot occur in this case
                        if (moiety1==moiety2 && cSS.numMoieties[moiety1]==1) {
                            rate[i]=0;
                            std::cout << "special case where 1 monomer remains" << std::endl;
                        }
                        else {
                            // check that two products are specified
                            if (reactionList[i].products.size()!=2) {
                                throw std::runtime_error("Linking reaction should have two products!");
                            }
                            else if (reactionList[i].linkFormed=="") {
                                throw std::runtime_error("Linking reaction should specify a link formed!");
                            }
                            std::string P = reactionList[i].products[0]+"-"+reactionList[i].linkFormed+"-"+reactionList[i].products[1]; // this is the product link string
                            double rateFDividedByStar=cfprime*cSS.numMoieties[moiety1]*cSS.numMoieties[moiety2]; // need num of each functional group overall tracked.
                            //if (rateDividedByStar<0) {
                                // looks like the depolymerization reaction should be happening, not the polymerization reaction
                            rate[i]=rateFDividedByStar*star;
                            totalRate+=rate[i];
                            
                            //}
                            // else {
                            //     rate[i]=rateDividedByStar*star;
                            //     totalRate+=rate[i];
                            //     rate[i+10]=0;
                            // }
                        } 
                    }
                    else {
                        throw std::runtime_error("Linking reaction should be bimolecular!");
                    }
                }
                else {
		    // bool debugging = false;
		    // if (debugging) {
		    double kf=rateConstants[R_i]["kf"];
                    // double kr=rateConstants[R_i]["kr"];
                    double kfprime = capK[R_i]["Kstep1"]*capK[R_i]["Kstep2"]*kf;
                    double cfprime = kfprime/(V*Na); // kf is trimolecular for these linking reactions: catalyst+2 reactants
                    // double krprime = kr * capK[R_i]["Kprimestep3"]*capK[R_i]["Kprimestep4"];
                    //double crprime = krprime; // kr is bimolecular (catalyst + link)
                    // if (cSS.chainPool.size()==1 && cSS.numMonomers["cEp"]==0) {
                    //    // no linking reactions possible... all linking reactions should be carried out as looping reactions
                    //    ratesDebugging[i]=0;
                    //    allOneChain=true;
                    //    std::cout << "all one chain" << std::endl;
                    //    continue;
                    //}
                    // calculate c parameter -- linking reaction should always involve 2 reactants and a catalyst molecule
                    // if (reactionList[i].reactants.size()==2) {
                    //     std::string moiety1=reactionList[i].reactants[0];
                    //     std::string moiety2=reactionList[i].reactants[1];

                    //     // check for the special case for which there is ONE monomer remaining -- monomer+monomer reaction cannot occur in this case
                    //     if (moiety1==moiety2 && cSS.numMoieties[moiety1]==1) {
                    //         // ratesDebugging[i]=0;
                    //         std::cout << "special case where 1 monomer remains" << std::endl;
                    //     }
                    //     else {
                    //         // check that two products are specified
                    //         if (reactionList[i].products.size()!=2) {
                    //             throw std::runtime_error("Linking reaction should have two products!");
                    //         }
                    //         else if (reactionList[i].linkFormed=="") {
                    //             throw std::runtime_error("Linking reaction should specify a link formed!");
                    //         }
                    //         // std::string P = reactionList[i].products[0]+"-"+reactionList[i].linkFormed+"-"+reactionList[i].products[1]; // this is the product link string
                    //         // double rateFDividedByStar=cfprime*cSS.numMoieties[moiety1]*cSS.numMoieties[moiety2]; // need num of each functional group overall tracked.
                    //         //if (rateDividedByStar<0) {
                    //             // looks like the depolymerization reaction should be happening, not the polymerization reaction
                    //         // ratesDebugging[i]=rateFDividedByStar*star;
			    ksteptimesstar[i-4]=star*cfprime;
                            //totalRate+=ratesDebugging[i];

                            //}
                            // else {
                            //     rate[i]=rateDividedByStar*star;
                            //     totalRate+=rate[i];
                            //     rate[i+10]=0;
                            // }
                    //    }
                    //}
                    //else {
                    //    throw std::runtime_error("Linking reaction should be bimolecular!");
                    //}

		    //}
                    // in this case, we have chain length dependent step growth reactions
                    TimePoint startcalcratesByLength=Time::now();
                    for (long length=1; length<=cSS.maxChainSize;++length) {
                        //double kf=rateConstants[R_i]["kf"];
                        //double kr=rateConstants[R_i]["kr"];
                        //double kfprime = capK[R_i]["Kstep1"]*kf;
                        //double cfprime = kfprime/(V*Na); // kf is trimolecular for these linking reactions: catalyst+2 reactants
                        //double krprime = kr * capK[R_i]["Kprimestep3"]*capK[R_i]["Kprimestep4"];
                        //double crprime = krprime; // kr is bimolecular (catalyst + link)
                        if (cSS.chainPool.size()==1 && cSS.numMonomers["cEp"]==0) {
                            // no linking reactions possible... all linking reactions should be carried out as looping reactions
                            rate[i]=0;
                            allOneChain=true;
                            std::cout << "all one chain" << std::endl;
                            continue;
                        }
                        // calculate c parameter -- linking reaction should always involve 2 reactants and a catalyst molecule
                        if (reactionList[i].reactants.size()==2) {
                            std::string moiety1=reactionList[i].reactants[0];
                            std::string moiety2=reactionList[i].reactants[1];

                            // check for the special case for which there is ONE monomer remaining -- monomer+monomer reaction cannot occur in this case
                            if (moiety1==moiety2 && cSS.numMoieties[moiety1]==1) {
                                rate[i]=0;
                                std::cout << "special case where 1 monomer remains" << std::endl;
                            }
                            else {
                                // check that two products are specified
                                if (reactionList[i].products.size()!=2) {
                                    throw std::runtime_error("Linking reaction should have two products!");
                                }
                                else if (reactionList[i].linkFormed=="") {
                                    throw std::runtime_error("Linking reaction should specify a link formed!");
                                }
                                double rateFDividedByStar;
                                if (cSS.numMoietiesByLength[length].find(moiety2)!=cSS.numMoietiesByLength[length].end() ) {
                                    rateFDividedByStar=cfprime*(normalizedCumSum[length][moiety1]*(cSS.numMoietiesByLength[length].find(moiety2)->second) +
                                (cSS.numMoieties[moiety1]-regCumSum[length][moiety1])*(cSS.numMoietiesByLength[length].find(moiety2)->second)*transport_factor[length]); // need num of each functional group overall tracked.
                                }
                                else {
                                    rateFDividedByStar=0;
                                }
                                ratesByLength[length][R_i]=rateFDividedByStar*star;
                                lengthDependentRateSums[R_i]+=ratesByLength[length][R_i];
                                totalRate += ratesByLength[length][R_i];
                            } 
                        }
                        else {
                            throw std::runtime_error("Linking reaction should be bimolecular!");
                        }
                    }
                    TimePoint stopcalcratesByLength;
                    timers["calcratesByLength"]+= (stopcalcratesByLength - startcalcratesByLength);
                    // finished first
                    //std::cout << "yeehaw " << std::endl;
                }
            }
            else if (reactionList[i].reactionType=="breaking" && reactionList[i].reactants.size()==2) {
                std::string R_i = reactionList[i].number;
                double kf=rateConstants[R_i]["kf"];
                double kr=rateConstants[R_i]["kr"];
                double kfprime = capK[R_i]["Kstep1"]*kf;
                double cfprime = kfprime/(V*Na); // kf is trimolecular for these linking reactions: catalyst+2 reactants
                double krprime = kr * capK[R_i]["Kprimestep4"]*capK[R_i]["Kprimestep5"]*capK[R_i]["Kprimestep6"];
                double crprime = krprime; // kr is bimolecular (catalyst + link)
                // calculate c parameter -- linking reaction should always involve 2 reactants and a catalyst molecule
                if (reactionList[i].reactants.size()==2) {
                    std::string moiety1=reactionList[i].reactants[0];
                    std::string moiety2=reactionList[i].reactants[1];

                    // check for the special case for which there is ONE monomer remaining -- monomer+monomer reaction cannot occur in this case
                    if (moiety1==moiety2 && cSS.numMoieties[moiety1]==1) {
                        rate[i]=0;
                        //throw std::runtime_error("in what context does this arise??");
                        //for breaking of an Em-ep-Em bond, need two Em's! Ultimately doesn't matter tho because the rate is defined to be 
                        //proportional to the number of links with the link name P, of which there 0.
                    } // this doesn't matter for breaking??
                    else {
                        // check that two products are specified
                        if (reactionList[i-10].products.size()!=2) {
                            throw std::runtime_error("Linking reaction should have two products!");
                        }
                        else if (reactionList[i-10].linkFormed=="") {
                            throw std::runtime_error("Linking reaction should specify a link formed!");
                        }
                        std::string P = reactionList[i-10].products[0]+"-"+reactionList[i-10].linkFormed+"-"+reactionList[i-10].products[1]; // this is the product link string
                        //double rateDividedByStar= -cfprime*cSS.numMoieties[moiety1]*cSS.numMoieties[moiety2] + crprime * cSS.links[P]; // need num of each functional group overall tracked.
                        double rateRDividedByStar= crprime * cSS.links[P]; // need num of each functional group overall tracked.
                        rate[i]=rateRDividedByStar*star;
                        totalRate+=rate[i];
                        
                        // if (rateDividedByStar<=0) {
                        //     // looks like the polymerization reaction should be happening, not the depolymerization reaction
                        //     rate[i]=0;
                        // }
                        // else {
                        //     rate[i]=rateDividedByStar*star;
                        //     totalRate+=rate[i];
                        // }
                        // if (rateDividedByStar*star!= -rate[i-10]) {
                        //     //throw std::runtime_error("mayday mayday maydayy");
                        //     int asdf=3;
                        // }
                    } 
                }
            }
            else {
                throw std::runtime_error("Reaction type not expected!");
            }
            
        }

        // statusFile << "end propensity calculations " << cSS.executionNumber << std::endl;
        // statusFile << std::flush;
   	// for the chain+chain reactions, check that the sum of entries in ratesByLength is equivalent to the corresponding rate in ratesDebugging...
	// std::vector<std::string> chainChainRxns = {"R5","R6","R7","R8"};
	// for (int j=0;j<chainChainRxns.size();++j) {
	// 	double rateSum=0;
	// 	for (int i = 0; i< ratesByLength.size(); ++i){
	// 	    
	// 	    rateSum+=ratesByLength[i][chainChainRxns[j]];    
	// 	} 
	// 	if (ratesDebugging[4+j]>0 && (rateSum/ratesDebugging[4+j] > 1.000001 || rateSum/ratesDebugging[4+j] < 0.999999)) {
	// 	    throw std::runtime_error("sum does not equal ratesByLength");
	// 	}
	// 	else if (ratesDebugging[4+j]==0) {
	// 		if (rateSum>0.000001) {
	// 			throw std::runtime_error("sum does not equual ratesByLength!");
	// 		}
	// 	}
	// }
	// if total rate is 0, then break!
	if (totalRate==0) {
	    std::cout << "total rate is 0, so we are done at execution number " << cSS.executionNumber << " num cEp is " << cSS.numMoieties["cEp"] << std::endl;
	    break;
	}
        if (cSS.executionNumber%100 ==0) {
            netRateAnalysis << simTime << " ";
            netRateAnalysis << over_x << " ";
            for (int i=0;i<rate.size();++i) {
                if (i==4) netRateAnalysis << lengthDependentRateSums["R5"] << " ";
                else if (i==5) netRateAnalysis << lengthDependentRateSums["R6"] << " ";
                else if (i==6) netRateAnalysis << lengthDependentRateSums["R7"] << " ";
                else if (i==7) netRateAnalysis << lengthDependentRateSums["R8"] << " ";
                else netRateAnalysis << rate[i] << " ";
            }
            netRateAnalysis << std::endl;
            //lengthDependentRateSums
        }
	
        // // verify that the chain length dependent rates are consistent with not-chain-length dependent rates
        // //if (transport_factor.size() > 100 && transport_factor[90]==1) {
        //     // loop through each step growth reactin
        // std::vector<std::string> R_i = {"R5","R6","R7","R8"};
        // for (int i = 0; i<R_i.size();++i) {
        //     double sumRatesByLength=0;
        //     for (int j = 0; j<ratesByLength.size(); ++j) {
        //         sumRatesByLength+=ratesByLength[j][R_i[i]];
        //     }
        //     // now we compare to the regular rate... 
        //     std::string reactant11=reactionList[i+4].reactants[0];
        //     std::string reactant22=reactionList[i+4].reactants[1];
        //     double kf = rateConstants[R_i[i]]["kf"];
        //     double kfprime = capK[R_i[i]]["Kstep1"]*kf;
        //     double cfprime = kfprime/(V*Na);
        //     double rateFDividedByStar=cfprime*cSS.numMoieties[reactant11]*cSS.numMoieties[reactant22]; // need num of each functional group overall tracked.
        //     double realRate = rateFDividedByStar*star;
        //     if (std::abs(realRate-sumRatesByLength) > 0.001) {
        //         throw std::runtime_error("ok we have found an issue");
        //     }
        // }
            
       // }

        // choose timestep tau
        double r1=distr(generator); // in range [0,1)
        if (r1==0){
            // want r1 on range (0,1] but generated on range [0,1)
            r1 = 1.0;
        } 
        double tau = (1/totalRate)*log(1/r1); // calculate timestep tau
        simTime+=tau; // increment time by timestep tau

        // choose a reaction event to take place
        double r2 =  distr(generator); // in range [0,1)
        std::vector<double> loopingRates = {}; // no looping reactions in this case
        chosenReactionChannel chosenRxn = selectReactionChannel(rate,loopingRates,ratesByLength, totalRate,r2,reactionList); // mu is the reaction channel selected
        if (timeShit) {
            stopPropensityTimer = Time::now(); // start clock
            timers["propensityAndChannelSelection"] +=(stopPropensityTimer-startPropensityTimer);
            startRxnExecution = Time::now(); //get the time at this instant 
        }

        // execute reaction (i.e., update system state). This could probably be placed in a separate function "updateSystemState()" which in turn could call a function for each of these if statements, which would probably contribute to readability of the code... could even be placed in a separate c++ file...

        double r3 = distr(generator); // in range [0,1)
        double r4 = distr(generator); // in range [0,1)
        
        // record the reaction type so that the molecular weight can be updated appropriately. 
        //if (chosenRxn.mu<rate.size() ) {
        cSS.mwv.rxnType=chosenRxn.reactionType;
        //}
        //else {
        
            //cSS.mwv.rxnType="looping";
        //}
        // let's set the reaction type in "chosenRxn"
        
        
        // statusFile << "about to call updateSystemState " << cSS.executionNumber << std::endl;
        // statusFile << std::flush;
        //if (mu<rate.size()) {
        updateSystemState(cSS, reactionList[chosenRxn.mu], r3, r4, myOptions,timers,chosenRxn,reactionList,chosenRxn.length,ratesByLength[chosenRxn.length][chosenRxn.chosenR_i],transport_factor,ksteptimesstar);
            // note that this assumes no steps were rejected
            //cSS.reactionsExecutedTally[mu]++; // increment tally of this reaction type
        //}
        //else {
            // I'm just going to hard code the looping reactions here for now
            // loopingRxn(cSS, candidateLoopingReactions[mu-rate.size()], myOptions,reactionList);
            //throw std::runtime_error("mayday mayday mayday...looping reaction selected but rate of looping should be 0");
        //}
        
        if (timeShit) {
            stopRxnExecution = Time::now(); //get the time at this instant 
            timers["execution"] += stopRxnExecution-startRxnExecution;
            startmolwt = Time::now(); //get the time at this instant 
        }

        // statusFile << "done with updatesystem state, call mol weight fxn " << cSS.executionNumber << std::endl;
        // statusFile << std::flush;
        cSS.mwv=molecularWeight(cSS.mwv); 
        
        if (timeShit) {
            stopmolwt = Time::now(); //get the time at this instant 
            timers["molwt"]+=(stopmolwt-startmolwt);
            startanalysis = Time::now(); // start clock

        }

        // statusFile << "mol weight done, start writing to files if correct execution " << cSS.executionNumber << std::endl; 
        // statusFile << std::flush;
        over_x=1-((1.0*cSS.numMoieties["cEp"])/(1.0*myConditions.numMoieties0["cEp"])); // monomer conversion
        // if (cSS.executionNumber%200==0) {
        //     coverages << simTime << "   " << over_x << "   "<< 1/denom << "   "
        //     << capK["R1"]["Kstep1"]*cSS.numMoieties["Ep"]*1.0/(V*Na)/denom << "   "
        //     << capK["R2"]["Kstep1"]*cSS.numMoieties["Emp"]*1.0/(V*Na)/denom << "    "
        //     << capK["R3"]["Kstep1"]*cSS.numMoieties["Em"]*1.0/(V*Na)/denom << " "
        //     << capK["R9"]["Kstep1"]*cSS.numMoieties["cEp"]*1.0/(V*Na)/denom << "   "
        //     << capK["R10"]["Kstep1"]*cSS.numMoieties["cE"]*cSS.numMoieties["cEp"]*1.0/(V*Na)*1.0/(V*Na)/denom << "  "
        //     << capK["R9"]["Kbind"]*capK["R9"]["Kstep1"]*cSS.numMoieties["cEp"]*cSS.numMoieties["cEp"]*1.0/(V*Na)*1.0/(V*Na)/denom << "    "
        //     << capK["R1"]["Kprimestep3"]*capK["R1"]["Kprimestep4"]*cSS.links[P1]*1.0/(V*Na)/denom << "  "
        //     << capK["R2"]["Kprimestep3"]*capK["R2"]["Kprimestep4"]*cSS.links[P2]*1.0/(V*Na)/denom << "  "
        //     << capK["R3"]["Kprimestep3"]*capK["R3"]["Kprimestep4"]*cSS.links[P3]*1.0/(V*Na)/denom << "  "
        //     << capK["R4"]["Kprimestep3"]*capK["R4"]["Kprimestep4"]*cSS.links[P4]*1.0/(V*Na)/denom << "  "
        //     << capK["R5"]["Kprimestep3"]*capK["R5"]["Kprimestep4"]*cSS.links[P5]*1.0/(V*Na)/denom << "  "
        //     << capK["R6"]["Kprimestep3"]*capK["R6"]["Kprimestep4"]*cSS.links[P6]*1.0/(V*Na)/denom << "  "
        //     << capK["R7"]["Kprimestep3"]*capK["R7"]["Kprimestep4"]*cSS.links[P7]*1.0/(V*Na)/denom << "  "
        //     << capK["R8"]["Kprimestep3"]*capK["R8"]["Kprimestep4"]*cSS.links[P8]*1.0/(V*Na)/denom << "  "
        //     << capK["R9"]["Kprimestep3"]*capK["R9"]["Kprimestep4"]*cSS.links[P9]*1.0/(V*Na)/denom << "  "
        //     << capK["R10"]["Kprimestep3"]*capK["R10"]["Kprimestep4"]*cSS.links[P10]*1.0/(V*Na)/denom << "  "
        //     << capK["R1"]["Kprimestep4"]*cSS.links[P1]*1.0/(V*Na)/denom << "  "
        //     << capK["R2"]["Kprimestep4"]*cSS.links[P2]*1.0/(V*Na)/denom << "  "
        //     << capK["R3"]["Kprimestep4"]*cSS.links[P3]*1.0/(V*Na)/denom << "  "
        //     << capK["R4"]["Kprimestep4"]*cSS.links[P4]*1.0/(V*Na)/denom << "  "
        //     << capK["R5"]["Kprimestep4"]*cSS.links[P5]*1.0/(V*Na)/denom << "  "
        //     << capK["R6"]["Kprimestep4"]*cSS.links[P6]*1.0/(V*Na)/denom << "  "
        //     << capK["R7"]["Kprimestep4"]*cSS.links[P7]*1.0/(V*Na)/denom << "  "
        //     << capK["R8"]["Kprimestep4"]*cSS.links[P8]*1.0/(V*Na)/denom << "  "
        //     << capK["R9"]["Kprimestep4"]*cSS.links[P9]*1.0/(V*Na)/denom << "  "
        //     << capK["R10"]["Kprimestep4"]*cSS.links[P10]*1.0/(V*Na)/denom << "  \n";
        // }

        // if (cSS.executionNumber%200==0) {
        //     coverages << simTime << "   " << over_x << "   "<< 1/denom << "   "
        //     << capK["R1"]["Kstep1"]*cSS.numMoieties["Ep"]*1.0/(V*Na)/denom << "   "
        //     << capK["R2"]["Kstep1"]*cSS.numMoieties["Emp"]*1.0/(V*Na)/denom << "    "
        //     << capK["R3"]["Kstep1"]*cSS.numMoieties["Em"]*1.0/(V*Na)/denom << " "
        //     << capK["R9"]["Kstep1"]*cSS.numMoieties["cEp"]*1.0/(V*Na)/denom << "   "
        //     << capK["R10"]["Kstep1"]*cSS.numMoieties["cE"]*cSS.numMoieties["cEp"]*1.0/(V*Na)*1.0/(V*Na)/denom << "  "
        //     << capK["R9"]["Kbind"]*capK["R9"]["Kstep1"]*cSS.numMoieties["cEp"]*cSS.numMoieties["cEp"]*1.0/(V*Na)*1.0/(V*Na)/denom << "    "
        //     << capK["R1"]["Kprimestep3"]*capK["R1"]["Kprimestep4"]*cSS.links[P1]*1.0/(V*Na)/denom << "  "
        //     << capK["R2"]["Kprimestep3"]*capK["R2"]["Kprimestep4"]*cSS.links[P2]*1.0/(V*Na)/denom << "  "
        //     << capK["R3"]["Kprimestep3"]*capK["R3"]["Kprimestep4"]*cSS.links[P3]*1.0/(V*Na)/denom << "  "
        //     << capK["R4"]["Kprimestep3"]*capK["R4"]["Kprimestep4"]*cSS.links[P4]*1.0/(V*Na)/denom << "  "
        //     << capK["R5"]["Kprimestep3"]*capK["R5"]["Kprimestep4"]*cSS.links[P5]*1.0/(V*Na)/denom << "  "
        //     << capK["R6"]["Kprimestep3"]*capK["R6"]["Kprimestep4"]*cSS.links[P6]*1.0/(V*Na)/denom << "  "
        //     << capK["R7"]["Kprimestep3"]*capK["R7"]["Kprimestep4"]*cSS.links[P7]*1.0/(V*Na)/denom << "  "
        //     << capK["R8"]["Kprimestep3"]*capK["R8"]["Kprimestep4"]*cSS.links[P8]*1.0/(V*Na)/denom << "  "
        //     << capK["R9"]["Kprimestep3"]*capK["R9"]["Kprimestep4"]*cSS.links[P9]*1.0/(V*Na)/denom << "  "
        //     << capK["R10"]["Kprimestep3"]*capK["R10"]["Kprimestep4"]*cSS.links[P10]*1.0/(V*Na)/denom << "  "
        //     << capK["R1"]["Kprimestep4"]*cSS.links[P1]*1.0/(V*Na)/denom << "  "
        //     << capK["R2"]["Kprimestep4"]*cSS.links[P2]*1.0/(V*Na)/denom << "  "
        //     << capK["R3"]["Kprimestep4"]*cSS.links[P3]*1.0/(V*Na)/denom << "  "
        //     << capK["R4"]["Kprimestep4"]*cSS.links[P4]*1.0/(V*Na)/denom << "  "
        //     << capK["R5"]["Kprimestep4"]*cSS.links[P5]*1.0/(V*Na)/denom << "  "
        //     << capK["R6"]["Kprimestep4"]*cSS.links[P6]*1.0/(V*Na)/denom << "  "
        //     << capK["R7"]["Kprimestep4"]*cSS.links[P7]*1.0/(V*Na)/denom << "  "
        //     << capK["R8"]["Kprimestep4"]*cSS.links[P8]*1.0/(V*Na)/denom << "  "
        //     << capK["R9"]["Kprimestep4"]*cSS.links[P9]*1.0/(V*Na)/denom << "  "
        //     << capK["R10"]["Kprimestep4"]*cSS.links[P10]*1.0/(V*Na)/denom << "  \n";
        // }


        if (cSS.executionNumber%200==0) {
            coverages << simTime << "   " << over_x << "   "<< 1/denom << "   "
            << capK["R1"]["Kstep1"]*cSS.numMoieties["Ep"]*1.0/(V*Na)/denom << "  "
            << capK["R2"]["Kstep1"]*cSS.numMoieties["Emp"]*1.0/(V*Na)/denom << "  "
            << capK["R3"]["Kstep1"]*cSS.numMoieties["Em"]*1.0/(V*Na)/denom << "  "
            << capK["R9"]["Kstep1"]*cSS.numMoieties["cEp"]*1.0/(V*Na)/denom << "  "
            << capK["R1"]["Kstep1"]*capK["R1"]["Kstep2"]*cSS.numMoieties["cEp"]*cSS.numMoieties["Ep"]*1.0/(V*Na)*1.0/(V*Na)/denom << "  "
            << capK["R2"]["Kstep1"]*capK["R2"]["Kstep2"]*cSS.numMoieties["cEp"]*cSS.numMoieties["Emp"]*1.0/(V*Na)*1.0/(V*Na)/denom << "  "
            << capK["R3"]["Kstep1"]*capK["R3"]["Kstep2"]*cSS.numMoieties["cEp"]*cSS.numMoieties["Em"]*1.0/(V*Na)*1.0/(V*Na)/denom << "  "
            << capK["R4"]["Kstep1"]*capK["R4"]["Kstep2"]*cSS.numMoieties["cEp"]*cSS.numMoieties["Emp"]*1.0/(V*Na)*1.0/(V*Na) << "  "
            << capK["R5"]["Kstep1"]*capK["R5"]["Kstep2"]*cSS.numMoieties["cE"]*cSS.numMoieties["Ep"]*1.0/(V*Na)*1.0/(V*Na) << "  "
            << capK["R6"]["Kstep1"]*capK["R6"]["Kstep2"]*cSS.numMoieties["cE"]*cSS.numMoieties["Emp"]*1.0/(V*Na)*1.0/(V*Na) << "  "
            << capK["R7"]["Kstep1"]*capK["R7"]["Kstep2"]*cSS.numMoieties["cE"]*cSS.numMoieties["Em"]*1.0/(V*Na)*1.0/(V*Na) << "  "
            << capK["R8"]["Kstep1"]*capK["R8"]["Kstep2"]*cSS.numMoieties["cE"]*cSS.numMoieties["Emp"]*1.0/(V*Na)*1.0/(V*Na) << "  "
            << capK["R9"]["Kstep1"]*capK["R9"]["Kstep2"]*cSS.numMoieties["cEp"]*cSS.numMoieties["cEp"]*1.0/(V*Na)*1.0/(V*Na)/denom << "  "
            << capK["R10"]["Kstep1"]*capK["R10"]["Kstep2"]*cSS.numMoieties["cE"]*cSS.numMoieties["cEp"]*1.0/(V*Na)*1.0/(V*Na)/denom << "  "
            
            << capK["R9"]["Kbind"]*cSS.numMoieties["cEp"]*cSS.numMoieties["cEp"]*1.0/(V*Na)*1.0/(V*Na)/denom << "  "
            << capK["R9"]["Kbind"]*cSS.numMoieties["cEp"]*cSS.numMoieties["Ep"]*1.0/(V*Na)*1.0/(V*Na) <<  "  " 
            << capK["R9"]["Kbind"]*cSS.numMoieties["cEp"]*cSS.numMoieties["Em"]*1.0/(V*Na)*1.0/(V*Na) <<  "  "
            << 2*capK["R9"]["Kbind"]*cSS.numMoieties["cEp"]*cSS.numMoieties["Emp"]*1.0/(V*Na)*1.0/(V*Na) <<  "  "
            << 4*capK["R9"]["Kbind"]*cSS.numMoieties["Emp"]*cSS.numMoieties["Emp"]*1.0/(V*Na)*1.0/(V*Na) <<  "  "
            << 2*capK["R9"]["Kbind"]*cSS.numMoieties["Emp"]*cSS.numMoieties["Em"]*1.0/(V*Na)*1.0/(V*Na) <<  "  "
            << 2*capK["R9"]["Kbind"]*cSS.numMoieties["Emp"]*cSS.numMoieties["Ep"]*1.0/(V*Na)*1.0/(V*Na) <<  "  "
            << capK["R9"]["Kbind"]*cSS.numMoieties["Ep"]*cSS.numMoieties["Ep"]*1.0/(V*Na)*1.0/(V*Na) <<  "  "
            << capK["R9"]["Kbind"]*cSS.numMoieties["Ep"]*cSS.numMoieties["Em"]*1.0/(V*Na)*1.0/(V*Na) <<  "  "
            << capK["R9"]["Kbind"]*cSS.numMoieties["Em"]*cSS.numMoieties["Em"]*1.0/(V*Na)*1.0/(V*Na) <<  "  "
            << capK["R1"]["Kprimestep5"]*capK["R1"]["Kprimestep6"]*cSS.links[P1]*1.0/(V*Na)/denom << "  "
            << capK["R2"]["Kprimestep5"]*capK["R2"]["Kprimestep6"]*cSS.links[P2]*1.0/(V*Na)/denom << "  "
            << capK["R3"]["Kprimestep5"]*capK["R3"]["Kprimestep6"]*cSS.links[P3]*1.0/(V*Na)/denom << "  "
            << capK["R4"]["Kprimestep5"]*capK["R4"]["Kprimestep6"]*cSS.links[P4]*1.0/(V*Na)/denom << "  "
            << capK["R5"]["Kprimestep5"]*capK["R5"]["Kprimestep6"]*cSS.links[P5]*1.0/(V*Na)/denom << "  "
            << capK["R6"]["Kprimestep5"]*capK["R6"]["Kprimestep6"]*cSS.links[P6]*1.0/(V*Na)/denom << "  "
            << capK["R7"]["Kprimestep5"]*capK["R7"]["Kprimestep6"]*cSS.links[P7]*1.0/(V*Na)/denom << "  "
            << capK["R8"]["Kprimestep5"]*capK["R8"]["Kprimestep6"]*cSS.links[P8]*1.0/(V*Na)/denom << "  "
            << capK["R9"]["Kprimestep5"]*capK["R9"]["Kprimestep6"]*cSS.links[P9]*1.0/(V*Na)/denom << "  "
            << capK["R10"]["Kprimestep5"]*capK["R10"]["Kprimestep6"]*cSS.links[P10]*1.0/(V*Na)/denom << "  "
            << capK["R1"]["Kprimestep4"]*capK["R1"]["Kprimestep5"]*capK["R1"]["Kprimestep6"]*cSS.links[P1]*1.0/(V*Na)/denom << "  "
            << capK["R2"]["Kprimestep4"]*capK["R2"]["Kprimestep5"]*capK["R2"]["Kprimestep6"]*cSS.links[P2]*1.0/(V*Na)/denom << "  "
            << capK["R3"]["Kprimestep4"]*capK["R3"]["Kprimestep5"]*capK["R3"]["Kprimestep6"]*cSS.links[P3]*1.0/(V*Na)/denom << "  "
            << capK["R4"]["Kprimestep4"]*capK["R4"]["Kprimestep5"]*capK["R4"]["Kprimestep6"]*cSS.links[P4]*1.0/(V*Na)/denom << "  "
            << capK["R5"]["Kprimestep4"]*capK["R5"]["Kprimestep5"]*capK["R5"]["Kprimestep6"]*cSS.links[P5]*1.0/(V*Na)/denom << "  "
            << capK["R6"]["Kprimestep4"]*capK["R6"]["Kprimestep5"]*capK["R6"]["Kprimestep6"]*cSS.links[P6]*1.0/(V*Na)/denom << "  "
            << capK["R7"]["Kprimestep4"]*capK["R7"]["Kprimestep5"]*capK["R7"]["Kprimestep6"]*cSS.links[P7]*1.0/(V*Na)/denom << "  "
            << capK["R8"]["Kprimestep4"]*capK["R8"]["Kprimestep5"]*capK["R8"]["Kprimestep6"]*cSS.links[P8]*1.0/(V*Na)/denom << "  "
            << capK["R9"]["Kprimestep4"]*capK["R9"]["Kprimestep5"]*capK["R9"]["Kprimestep6"]*cSS.links[P9]*1.0/(V*Na)/denom << "  "
            << capK["R10"]["Kprimestep4"]*capK["R10"]["Kprimestep5"]*capK["R10"]["Kprimestep6"]*cSS.links[P10]*1.0/(V*Na)/denom << "  \n";
        }

        if ((cSS.executionNumber%200==0) && myOptions.Changxia) {
            recordChangxiaConcentrations(cSS.numMoieties,myConditions.numMoieties0,over_x,simTime,ChangxiaOutput,V, Na,cSS.mwv.Mn,cSS.mwv.Mw);
            printAllLinksToFile(cSS.links,linksFile,allPossibleLinks,simTime,V,Na);
            ChangxiaOutput << std::flush;
            linksFile << std::flush;
        }
        if (myOptions.Tirrell && myOptions.trackFullStructureAsGraph && (cSS.executionNumber%100==0)) { 
            TirrellTheory(cSS.chainPool,cSS.numMoieties,myConditions.numMoieties0,cSS.numMonomers,cSS.monomers,TirrellOutput,simTime,cSS.mwv.Mn,cSS.mwv.Mw);
        }
        if ((cSS.executionNumber%10==0) && myOptions.eG) {
            printEndGrps(cSS,eG,simTime);
        }
        if (simTime>60*5) {
            if (didIprint.find("5")==didIprint.end()) {
                print_graph(cSS.chainPool,runID+"5min");
                didIprint["5"]=true;
            }
        }
        
        if (simTime>60*30) {
            if (didIprint.find("30")==didIprint.end()) {
                print_graph(cSS.chainPool,runID+"30min");
                didIprint["30"]=true;
            }
        }
          
        if (simTime>60*60) {
            if (didIprint.find("60")==didIprint.end()) {
                didIprint["60"]=true;
                print_graph(cSS.chainPool,runID+"60min");
            }
        }
          
        if (simTime>60*90) {
            if (didIprint.find("90")==didIprint.end()) {
                print_graph(cSS.chainPool,runID+"90min");
                didIprint["90"]=true;
            }
        }

        if (simTime>60*120) {
            if (didIprint.find("120")==didIprint.end()) {
                print_graph(cSS.chainPool,runID+"120min");
                didIprint["120"]=true;
            }
        }

        if (simTime>60*150) {
            if (didIprint.find("150")==didIprint.end()) {
                print_graph(cSS.chainPool,runID+"150min");
                didIprint["150"]=true;
            }
        }

        if (simTime>60*210) {
            if (didIprint.find("210")==didIprint.end()) {
                print_graph(cSS.chainPool,runID+"210min");
                didIprint["210"]=true;
            }
        }

        if (simTime>60*420) {
            if (didIprint.find("420")==didIprint.end()) {
                print_graph(cSS.chainPool,runID+"420min");
                didIprint["420"]=true;
            }
        }

           
        if (timeShit) {
            stopanalysis = Time::now(); // stop clock
            timers["analysisTime"] += (stopanalysis-startanalysis);
            stopKMCInnerLoop = Time::now(); //get the timepoint at this instant 
            innerLoopTime += stopKMCInnerLoop-startKMCInnerLoop;
        }
        // statusFile << "stop writing to files: end kMC loop " << cSS.executionNumber << std::endl; 
        // statusFile << std::flush;
    }
    if (timeShit) {
        stopKMC = Time::now(); //get the timepoint at this instant 
        timers["KMCLoop"] = (stopKMC-startKMC);
        timers["KMCInner"] += (innerLoopTime);
        startprint = Time::now(); // start clock
    }

    if (myOptions.printWeightFractions) {
        printwfs(cSS.chainPool,strdate,runID,myConditions.monomerWeights,cSS.numMonomers,myConditions.simulationSize);
    }
    if (true) {
        printExecutions(cSS.reactionsExecutedTally,reactionList,strdate,runID);
    }

    // if Tirrell option is selected, need to close the Tirrell output file
    if (myOptions.Tirrell) {
        TirrellOutput.close();
    }
    ChangxiaOutput.close();
    if (myOptions.eG) {
        eG.close();
    }

    print_graph(cSS.chainPool,runID+"end");

    if (timeShit) {
        stopprint = Time::now(); // stop clock
        timers["printGraph"] = (stopprint - startprint);
        stop = Time::now(); // stop clock
        // Subtract stop and start timepoints and
        // cast it to required unit. Predefined units
        // are nanoseconds, microseconds, milliseconds,
        // seconds, minutes, hours. Use duration_cast()
        // function.
        timers["overall"] = (stop - start);
        // To get the value of duration use the count()
        // member function on the duration object
        //std::cout << timers["overall"].count() << std::endl;
    }
    
    // conv.close();
    linksFile.close();
    netRateAnalysis.close();
        // print reactions to a file for debugging
    std::ofstream reactionsFile;
    reactionsFile.open("./Output/" + strdate + runID + "reactions.txt");
    // each line corresponds to a reaction
    for (int i=0; i<reactionList.size(); ++i) {
        reactionsFile << "Reaction type: " << reactionList[i].reactionType << "  " << " reactants: ";
        for (int j=0; j<reactionList[i].reactants.size(); ++j) {
            reactionsFile << reactionList[i].reactants[j] << " ";
        }
        reactionsFile << " linkToBust: " << reactionList[i].linkToBust;
        reactionsFile << " linkFormed: " << reactionList[i].linkFormed;
        reactionsFile << " gasEvolved: " << reactionList[i].gasEvolved;
        reactionsFile << " catalystEffect: " << reactionList[i].catalystEffect;
        reactionsFile << " A: " << reactionList[i].A;
        reactionsFile << " executions: " << cSS.reactionsExecutedTally[i];
        //reactionsFile << " " cSS.rxn
        // reactionsFile << " products: ";
        // for (int j=0; j<reactionList[i].products.size(); ++j) {
        //     reactionsFile << reactionList[i].products[j] << " ";
        // }
        reactionsFile << std::endl;
    }
    // let's also print to a file times, we'll just drop it in the reactions file bc im lazy...
    std::cout << "total  time " << timers["overall"].count() << "  adj breaking  " << timers["adjBreaking"].count()  << "  connectedcomp  " << timers["connectedcomponents"].count() << " split  " << timers["split"].count() << " timers cns + numMoietiesByLength " << timers["calccns"].count() << std::endl;   
    reactionsFile.close();
    coverages.close();
    
    // statusFile << "end simulation " << cSS.executionNumber << std::endl; 
    // statusFile << std::flush;
    // statusFile.close();
    return 0;
}
