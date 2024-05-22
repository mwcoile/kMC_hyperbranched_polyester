#include"inputs.h"
// std::map< std::string, std::map < std::string, double > > getCapK() {
//     // read in equilibrium constants ("cap K") from a file
//     // the file should be formatted as follows:
//     // the first row contains the numbers of the reactions, e.g., R1, R2, R3, etc.
//     // the first column contains the step names, e.g., Kstep1, Kprimestep2, bind, etc.
//     // entries in each row, column are the equilibrium constants for the corresponding reaction and step
//     // these are stored in a map of maps, where the first key is the step name and the second key is the reaction number
//     // the value is the equilibrium constant
//     // e.g., capK["Kstep1"]["R1"] is the equilibrium constant for step 1 of reaction 1
//     // e.g., capK["bind"]["R9"] is the equilibrium constant for the alternate binding reaction of reaction 9
    
//     std::map< std::string, std::map < std::string, double > > capK;
//     std::ifstream capKfile("capK.txt");
//     std::string line;
//     std::string word;
//     std::vector<std::string> words;
//     std::vector<std::string> reactions;
//     std::vector<std::string> steps;
//     std::map<std::string,double> stepK;
//     std::string stepName;
//     std::string rxnNumber;
//     double K;
//     bool firstLine=true;
//     while (std::getline(capKfile,line)) {
//         words=getNextLineAndSplitIntoTokens(capKfile);
//         if (firstLine) {
//             for (int i=1; i<words.size(); i++) {
//                 reactions.push_back(words[i]);
//             }
//             firstLine=false;
//         } else {
//             stepName=words[0];
//             steps.push_back(stepName);
//             for (int i=1; i<words.size(); i++) {
//                 rxnNumber=reactions[i-1];
//                 K=std::stod(words[i]);
//                 stepK[rxnNumber]=K;
//             }
//             capK[stepName]=stepK;
//         }
//     }
//     return capK;
// }

// write function split to separate a string into a vector of strings given by some delimiter

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::map< std::string, std::map < std::string, double > > getCapK(std::string strdate, std::string runID) {
    // read in equilibrium constants ("cap K") from a file
    // the file should be formatted as follows:
    // the first row contains the step names, e.g., Kstep1, Kprimestep2, bind, etc.
    // the first column contains the numbers of the reactions, e.g., R1, R2, R3, etc.
    // entries in each row, column are the equilibrium constants for the corresponding step and reaction
    // these are stored in a map of maps, where the first key is the step name and the second key is the reaction number
    // the value is the equilibrium constant
    // e.g., capK["R1"]["Kstep1"] is the equilibrium constant for step 1 of reaction 1
    // e.g., capK["R9"]["bind"] is the equilibrium constant for the alternate binding reaction of reaction 9
    
    std::map< std::string, std::map < std::string, double > > capK;


    std::ifstream 
    capKfile("capK.txt");
    std::string line;
    std::string word;
    std::vector<std::string> words;
    std::vector<std::string> reactions;
    std::vector<std::string> steps;
    std::map<std::string,double> stepK;
    std::string stepName;
    std::string rxnNumber;
    double K;
    bool firstLine=true;
    while (std::getline(capKfile,line)) {
        words=split(line,'\t');
        if (firstLine) {
            for (int i=1; i<words.size(); i++) {
                reactions.push_back(words[i]);
            }
            firstLine=false;
        } else {
            stepName=words[0];
            steps.push_back(stepName);
            for (int i=1; i<words.size(); i++) {
                rxnNumber=reactions[i-1];
                K=std::stod(words[i]);
                stepK[rxnNumber]=K;
            }
            capK[stepName]=stepK;
        }
    }
    double test2=capK["R1"]["Kstep1"];
    return capK;
}

std::map< std::string, std::map< std::string, double>> getRateConstants() {
    // read rate constants in from file. File should be formatted as follows:
    // first row contains whether the reaction is a forward or reverse reaction (kf or kr)
    // first column contains the reaction numbers, e.g., R1, R2, R3, etc.
    // function returns the rate constants in a map of maps, where the first key 
    // is the reaction number and the second key is the forward or reverse rate constant (kf or kr)
    // e.g., rateConstants["R1"]["kf"] is the forward rate constant for reaction 1

    std::map< std::string, std::map< std::string, double>> rateConstants;
    std::ifstream rateConstantsFile("rateConstants.txt");
    std::string line;
    std::string word;
    std::vector<std::string> words;
    std::vector<std::string> reactions;
    std::vector<std::string> rateTypes;
    std::map<std::string,double> rateType;
    std::string rateName;
    std::string rxnNumber;
    double rate;
    bool firstLine=true;
    while (std::getline(rateConstantsFile,line)) {
        words=split(line,'\t');
        if (firstLine) {
            for (int i=1; i<words.size(); i++) {
                rateTypes.push_back(words[i]);
            }
            firstLine=false;
        } else {
            rxnNumber=words[0];
            reactions.push_back(rxnNumber);
            for (int i=1; i<words.size(); i++) {
                rateName=rateTypes[i-1];
                rate=std::stod(words[i]);
                rateType[rateName]=rate;
            }
            rateConstants[rxnNumber]=rateType;
        }
    }
    return rateConstants;
}

void readInLHHWLinking(std::vector<Reaction> &reactionList) {
    
    Reaction myReaction;
    myReaction.reactionType="linking";
    myReaction.number="R1";
    // R1: cEp + Ep <-> Emp-ep-E
    myReaction.reactants.push_back("cEp");
    myReaction.reactants.push_back("Ep");
    myReaction.products.push_back("Emp");
    myReaction.products.push_back("E");
    myReaction.linkFormed="ep";
    reactionList.push_back(myReaction);
    
    // R2: cEp + Emp <-> Emp-ep-Em
    myReaction.number="R2";
    myReaction.reactants[0]="cEp";
    myReaction.reactants[1]="Emp";
    myReaction.products[0]="Emp";
    myReaction.products[1]="Em";
    myReaction.linkFormed="ep";
    reactionList.push_back(myReaction);

    // R3: cEp + Em <-> Emp-em-E
    myReaction.number="R3";
    myReaction.reactants[0]="cEp";
    myReaction.reactants[1]="Em";
    myReaction.products[0]="Emp";
    myReaction.products[1]="E";
    myReaction.linkFormed="em";
    reactionList.push_back(myReaction);

    // R4: cEp + Emp <-> Emp-em-Ep
    myReaction.number="R4";
    myReaction.reactants[0]="cEp";
    myReaction.reactants[1]="Emp";
    myReaction.products[0]="Emp";
    myReaction.products[1]="Ep";
    myReaction.linkFormed="em";
    reactionList.push_back(myReaction);

    // R5: cE + Ep-> Em-ep-E
    myReaction.number="R5";
    myReaction.reactants[0]="cE";
    myReaction.reactants[1]="Ep";
    myReaction.products[0]="Em";
    myReaction.products[1]="E";
    myReaction.linkFormed="ep";
    reactionList.push_back(myReaction);

    // R6: cE + Emp -> Em-ep-Em
    myReaction.number="R6";
    myReaction.reactants[0]="cE";
    myReaction.reactants[1]="Emp";
    myReaction.products[0]="Em";
    myReaction.products[1]="Em";
    myReaction.linkFormed="ep";
    reactionList.push_back(myReaction);

    // R7: cE + Em -> Em-em-E
    myReaction.number="R7";
    myReaction.reactants[0]="cE";
    myReaction.reactants[1]="Em";
    myReaction.products[0]="Em";
    myReaction.products[1]="E";
    myReaction.linkFormed="em";
    reactionList.push_back(myReaction);

    // R8: cE + Emp -> Em-em-Ep
    myReaction.number="R8";
    myReaction.reactants[0]="cE";
    myReaction.reactants[1]="Emp";
    myReaction.products[0]="Em";
    myReaction.products[1]="Ep";
    myReaction.linkFormed="em";
    reactionList.push_back(myReaction);

    // R9: cEp + cEp -> Emp-ep-cE
    myReaction.number="R9";
    myReaction.reactants[0]="cEp";
    myReaction.reactants[1]="cEp";
    myReaction.products[0]="Emp";
    myReaction.products[1]="cE";
    myReaction.linkFormed="ep";
    reactionList.push_back(myReaction);

    // R10: cE + cEp -> Em-ep-cE
    myReaction.number="R10";
    myReaction.reactants[0]="cE";
    myReaction.reactants[1]="cEp";
    myReaction.products[0]="Em";
    myReaction.products[1]="cE";
    myReaction.linkFormed="ep";
    reactionList.push_back(myReaction);

}

void readInLHHWBreaking(std::vector<Reaction> &reactionList) {
    // loop through each entry in reactionList and transform it from
    // a linking reaction to a breaking reaction
    for (int i=0; i<reactionList.size(); ++i) {
        if (reactionList[i].reactionType!="linking") {
            continue;
        }
        Reaction myReaction;
        myReaction.reactionType="breaking";
        myReaction.number=reactionList[i].number;
        myReaction.reactants.push_back(reactionList[i].products[0]);
        myReaction.reactants.push_back(reactionList[i].products[1]);
        myReaction.products.push_back(reactionList[i].reactants[0]);
        myReaction.products.push_back(reactionList[i].reactants[1]);
        myReaction.linkToBust=reactionList[i].linkFormed;
        reactionList.push_back(myReaction);
    }
    return;
}


// cEp + cEp -> k_Arrhenius(6.21E+12,0) Emp9 + cE9 : linking, [eplinkNotFormed]*
// cE + cEp -> k_Arrhenius(6.21E+12,0) Em10 + cE10: linking, [eplinkNotFormed]*

// std::map<std::string,std::map<std::string,double>> grabRateConstants() {
//     // read in rate constants from a file rateConstants.txt
//     // the first row gives the the step name, e.g. "step1f"
//     // the first column gives the overall reaction name, e.g. "R1"
//     // so that entries in the rateConstants variable can be accessed like rateConstants["step1f"]["R1"]
//     std::ifstream rateConstantsFile("rateConstants.txt");
//     std::string line;
//     std::vector<std::string> lineVector;
//     std::vector<std::string> stepNames;
//     std::map<std::string,std::map<std::string,double>> rateConstants;
//     std::map<std::string,double> stepRateConstants;
//     std::string stepName;
//     std::string reactionName;
//     double rateConstant;
//     int i=0;
//     while (std::getline(rateConstantsFile,line)) {
//         // note that rate constants are delimited by tabs
//         lineVector=split(line,'\t');
//         if (i==0) {
//             // the line vector gives the step names
//             stepNames=lineVector;
//             // however some step names have \r in them so we need to remove that
//             for (int j=0;j<stepNames.size();j++) {
//                 stepNames[j].erase(std::remove(stepNames[j].begin(),stepNames[j].end(),'\r'),stepNames[j].end());
//             }
//         } else {
//             reactionName=lineVector[0];
//             for (int j=1;j<lineVector.size();j++) {
//                 stepName=stepNames[j];
//                 rateConstant=std::stod(lineVector[j]);
//                 rateConstants[stepName][reactionName]=rateConstant;
//             }
//         }
//         i++;
//     }
//     return rateConstants;
// }


void readInLactoneLoopingReactions(std::vector<Reaction> &reactionList) {
    // note which 4 reactions can be looping reactions

    // R5
    Reaction loopingR5 = reactionList[4];
    loopingR5.reactionType = "looping";
    reactionList.push_back(loopingR5);
    // R6
    Reaction loopingR6 = reactionList[5];
    loopingR6.reactionType = "looping";
    reactionList.push_back(loopingR6);
    // R7
    Reaction loopingR7 = reactionList[6];
    loopingR7.reactionType = "looping";
    reactionList.push_back(loopingR7);
    // R8
    Reaction loopingR8 = reactionList[7];
    loopingR8.reactionType = "looping";
    reactionList.push_back(loopingR8);
    if (loopingR8.reactants[1] != "Em*p") {
        throw std::runtime_error("input file not as expected.. ");
    }

}

// void readInLactoneReverseReactions(std::vector<Reaction> &reactionList) {
//     // read in the reactions needed for Changxia's system
//     Reaction reverse;

//     // REACTIONS 1-3: cEp + pOH -> Emp + XXX

//     // reaction 1 Emp-ep-Em -> cEp + Emp
//     reverse.reactionType="breaking";
//     reverse.reactants.push_back("Emp");
//     reverse.reactants.push_back("Em");
//     reverse.linkToBust="ep";
//     reverse.products.push_back("cEp");
//     reverse.products.push_back("Emp");
//     reverse.A=0; // rate constant can just be set when calculating propensity functions
//     reverse.Ea=0;
//     reactionList.push_back(reverse);
//     // reaction 2 Emp-ep-E -> cEp + Ep
//     reverse.reactants[1]="E";
//     reverse.products[1]="Ep";
//     reactionList.push_back(reverse);
//     // reaction 3 Emp-ep-cE -> cEp + cEp
//     reverse.reactants[1]="cE";
//     reverse.products[1]="cEp";
//     reactionList.push_back(reverse);

//     // REACTIONS 4-6: cE + pOH -> Em + XXX

//     // reaction 4 Em-ep-Em -> cE + Emp
//     reverse.reactants[0]="Em";
//     reverse.reactants[1]="Em";
//     reverse.products[0]="cE";
//     reverse.products[1]="Emp";
//     reactionList.push_back(reverse);
//     // reaction 5 Em-ep-E -> cE + Ep
//     reverse.reactants[1]="E";
//     reverse.products[1]="Ep";
//     reactionList.push_back(reverse);
//     // reaction 6 Em-ep-cE -> cE + cEp
//     reverse.reactants[1]="cE";
//     reverse.products[1]="cEp";
//     reactionList.push_back(reverse);

//     // REACTIONS 7-8: cEp + mOH -> Emp + XXX

//     // reaction 7 Emp-em-E -> cEp + Em
//     reverse.linkToBust="em";
//     reverse.reactants[0]="Emp";
//     reverse.reactants[1]="E";
//     reverse.products[0]="cEp";
//     reverse.products[1]="Em";
//     reactionList.push_back(reverse);

//     // reaction 8 Emp-em-Ep -> cEp + Emp
//     reverse.reactants[1]="Ep";
//     reverse.products[1]="Emp";
//     reactionList.push_back(reverse);

//     // REACTIONS 9-10 cE + mOH -> Em + XXX 

//     // reaction 9 Em-em-E -> cE + Em
//     reverse.reactants[0]="Em";
//     reverse.reactants[1]="E";
//     reverse.products[0]="cE";
//     reverse.products[1]="Em";
//     reactionList.push_back(reverse);

//     // reaction 10 Em-em-Ep -> cE + Emp
//     reverse.reactants[1]="Ep";
//     reverse.products[1]="Emp";
//     reactionList.push_back(reverse);
// }

std::vector<double> convertStringVectortoDoubleVector(std::vector<std::string>& stringVector) {
    std::vector<double> doubleVector(stringVector.size());
    for (int i=0;i<stringVector.size();i++) {
        doubleVector[i]=std::stod(stringVector[i]);
    }
return doubleVector;
}

// Function to read input data from a file, courtesy of https://stackoverflow.com/questions/1120140/how-can-i-read-and-parse-csv-files-in-c
std::vector<std::string> getNextLineAndSplitIntoTokens(std::istream& str) {
    std::vector<std::string> result;
    std::string line;
    std::getline(str,line);

    std::stringstream lineStream(line);
    std::string cell;

    while(std::getline(lineStream,cell, ','))
    {
        cell.erase(remove(cell.begin(),cell.end(),' '),cell.end()); // remove white space
        result.push_back(cell);
    }
    // This checks for a trailing comma with no data after it.
    if (!lineStream && cell.empty())
    {
        // If there was a trailing comma then add an empty element.
        result.push_back("");
    }
    return result;
}

Reaction getNextRxn(std::string line)
{
    // old input was std::istream& str 

    //std::pair<std::vector<std::string>,std::vector<std::string>> reactionInput;
    Reaction reactionInput;
    //std::vector<std::string> reactionReactants;
    //std::vector<std::string> reactionProducts;
    //std::string line;
    //std::getline(str,line);
    int arrow = line.find('->');
    int kcalc=line.find('k_');
    int closeParen=line.find(')');
    int openParen=line.find('(');
    int colon=line.find(':'); // find index of first ':' to occur in the line
    int comma=line.find(',',colon); // find the first occurrence of ',' that occurs after ':'
    std::string reactantString = line.substr(0,arrow-1);
    std::string productString = line.substr(closeParen+1,colon-closeParen-1);
    std::string arrheniusString = line.substr(openParen+1,closeParen-openParen-1);
    std::stringstream reactantStream(reactantString);
    std::stringstream productStream(productString);
    std::stringstream arrheniusStream(arrheniusString);

    std::string cell;
 
    // Reactants
    while(std::getline(reactantStream,cell, '+'))
    {
        cell.erase(remove(cell.begin(),cell.end(),' '),cell.end());
        reactionInput.reactants.push_back(cell);
    }
    // Products
    while(std::getline(productStream,cell, '+'))
    {
        cell.erase(remove(cell.begin(),cell.end(),' '),cell.end());
        reactionInput.products.push_back(cell);
    }
    // Rate Constant
    int count=0;
    while(std::getline(arrheniusStream,cell, ','))
    {
        cell.erase(remove(cell.begin(),cell.end(),' '),cell.end());
        if (count==0) {
            reactionInput.A=stod(cell);
        }
        else if (count==1) {
            reactionInput.Ea=stod(cell);
        }
        else {
             throw std::runtime_error( "Ruh roh! Check your input file." );
        }
        count++;
    }

    // Reaction type
    std::string rxnType=line.substr(colon+1,(comma-colon-1)); // grab the reaction type from the input string
    rxnType.erase(remove(rxnType.begin(),rxnType.end(),' '),rxnType.end()); // get rid of white space
    reactionInput.reactionType=rxnType; // set the reaction type 
    if (rxnType!="breaking" && rxnType!="linking" && rxnType!="exchange" && 
    rxnType!="looping" && rxnType!="release" && rxnType!="link-linking" && 
    rxnType!="carbamatedecomposition" && rxnType!="catBinding") {
        throw std::runtime_error( "Ruh roh! Check your input file. Rxn type no good" );
    }
    if (rxnType=="breaking" && reactionInput.reactants.size()!=1) {
        //throw std::runtime_error( "Ruh roh! bimolecular depolymerization reactions not yet supported!" );
        std::string linkToBust=line.substr(comma+1); // line.substr(colon+1,std::string::npos);
        linkToBust.erase(remove(linkToBust.begin(),linkToBust.end(),' '),linkToBust.end()); // get rid of white space
        reactionInput.linkToBust=linkToBust;
        // if (linkToBust.find('*') != std::string::npos) {
        //     reactionInput.catalystEffect="released";
        // }

    }
    if (rxnType=="link-linking") {
        reactionInput.linkFormed = reactionInput.products[0]; // e.g. for allophanate forming reaction, product 0 should be allophanate
        reactionInput.linkToBust = reactionInput.reactants[0]; // e.g. for allophanate forming reaction, reactant 1 should be urethane / carbamate / NCO-OH
    }
    if (rxnType=="carbamatedecomposition") {
        // expecting something like NCO-OH -> NH2 + alkene + CO2
        // first and only entry in carbamate decomposition reaction should be the carbamate / urethane linkage
        reactionInput.linkToBust = reactionInput.reactants[0];
        // the last product entry in carbamate decomposition should be the CO2 produced
        if (reactionInput.products.size()!=3) {
            throw std::runtime_error("need 3 products for the carbamate decomposition reaction");
        }
        reactionInput.gasEvolved = reactionInput.products[2]; 
    }
    if (rxnType == "linking" || rxnType == "looping") {
        // New edge name
        std::string linkFormed=line.substr(comma+1); // line.substr(colon+1,std::string::npos);
        linkFormed.erase(remove(linkFormed.begin(),linkFormed.end(),' '),linkFormed.end()); // get rid of white space
        reactionInput.linkFormed=linkFormed;
    }

    return reactionInput;
}

InitialConditions readInitialConditionsFromInput() {
    InitialConditions myConditions = {};

    std::ifstream file("input.txt");
    // std::string myInput = "Enter system temperature (K),323\n"
    // "Enter monomer names in a comma separated list,AA,BB\n"
    // "Enter concentrations of monomers above in M,1.0,1.0\n"
    // "Enter molecular weights of monomers above in g/mol,100, 100\n"
    // "Enter simulation size in terms of desired number of monomers to simulate,500\n"
    // "Enter the functional groups and the monomers on which they are present. If the same functional group is on more than one monomer add another line\n"
    // "AA, NCO, 2\n"
    // "BB, OH, 2\n"
    // "Enter reactions present that the functional groups may undergo\n"
    // "NCO + OH -> k_Arrhenius(1,0) NCO-OH: linking\n"
    // "NCO-OH -> k_Arrhenius(0.5,0) NCO + OH: breaking\n"
    // "NCO-OH + NCO -> k_Arrhenius(0.2,0) allophanate: link-linking\n"
    // "NCO-OH -> k_Arrhenius(0.2,0) NH2 + alkene + CO2: carbamate decomposition\n";
    // file << myInput << std::endl;
    // Read in temperature (K)
    // std::istringstream file(myInput);
    myConditions.T = std::stod(getNextLineAndSplitIntoTokens(file)[1]); 

    // Read in monomer names
    std::vector<std::string> speciesList=getNextLineAndSplitIntoTokens(file);
    speciesList.erase(speciesList.begin());
    // throw error if empty member of set

    // Read in concentrations 
    std::vector<std::string> speciesConcentrationsString=getNextLineAndSplitIntoTokens(file);
    speciesConcentrationsString.erase(speciesConcentrationsString.begin()); 
    if (speciesConcentrationsString.size() != speciesList.size()) {
        // throw error if any species does not have a defined concentration
        throw std::runtime_error( "Ruh roh! Check your input file." );
    }
    std::vector<double> speciesConcentrations=convertStringVectortoDoubleVector(speciesConcentrationsString);
    std::map<std::string,double> initialConcentrations;
    for (int i=0;i<speciesConcentrations.size();i++) {
        initialConcentrations[speciesList[i]]=speciesConcentrations[i];
    }

    myConditions.totalMonomerConcentration = std::accumulate(speciesConcentrations.begin(), 
    speciesConcentrations.end(), decltype(speciesConcentrations)::value_type(0));

    // Read in monomer molecular weights
    std::vector<std::string> monomerWeightsString=getNextLineAndSplitIntoTokens(file);
    monomerWeightsString.erase(monomerWeightsString.begin()); 
    if (monomerWeightsString.size() != speciesList.size()) {
        // throw error if any species does not have a defined molecular weight, or if there are more molecular weight entries than species
        throw std::runtime_error( "Ruh roh! Check your input file." );
    }
    std::vector<double> monomerWeights=convertStringVectortoDoubleVector(monomerWeightsString); // convert input strings to doubles

    // Read in simulation size 
    signed long simulationSize=strtol(getNextLineAndSplitIntoTokens(file)[1].c_str(),NULL,10);

    getNextLineAndSplitIntoTokens(file); // skip this line

    // Read in the functional groups these monomers have
    std::map<std::string,std::map<std::string,signed long>> monomerFuns; // declare map of functional groups to be initialized in the subsequent lines
    std::vector<std::string> nextLine = getNextLineAndSplitIntoTokens(file); // read in fun groups on monomers for the first fun group type

    while (nextLine[0] != "Enterreactionspresentthatthefunctionalgroupsmayundergo") {
        std::string monomer=nextLine[0];
        std::string fun=nextLine[1];
        int number=std::stoi(nextLine[2]);
        monomerFuns[monomer][fun]=number; 
        nextLine = getNextLineAndSplitIntoTokens(file);
    }

    SystemVariables cSS={}; // current system state

    // loop through fun groups and determine the monomers on which each fun group exists
    // map from fun group name to a vector of monomers on which they appear
    std::map<std::string,std::map<std::string,signed long>>::iterator it_monomerFuns=monomerFuns.begin();
    while (it_monomerFuns != monomerFuns.end()) {
        std::map<std::string,signed long>::iterator it2=it_monomerFuns->second.begin();
        std::string monomerName=it_monomerFuns->first;
        while (it2 != it_monomerFuns->second.end()) {
            std::string funGroupName = it2->first; // this is the functional group name that appears on monomer m
            // Regardless of whether the fun group is in the list already or not, the code should be the same
            cSS.funGroups[funGroupName].push_back(monomerName);
            it2++;
        }
        it_monomerFuns++;
    }
    std::vector<Reaction> reactionList;
    // Read in list of reactions -- this is the hard part of reading from the input file
    
    // figure out here how to read till end of file..
    std::string line;
    while (std::getline(file,line)) {
        reactionList.push_back(getNextRxn(line));
    }

    // create the monomers and then run the reaction
    for (int i=0;i<speciesConcentrations.size();i++) {
        cSS.numMonomers[speciesList[i]]=std::round(1.0*simulationSize*speciesConcentrations[i]/myConditions.totalMonomerConcentration);
    }

    cSS.mwv.weightCutoff=600; // declare molecular weight at which an oligomer is considered a polymer and will show up in GPC trace (600 g/mol for changxia's system)

    // Initialize list of all monomers present at start
    for (int i=0;i<speciesList.size();i++) {
        Monomer newMonomer;
        newMonomer.name=speciesList[i];
        newMonomer.funs=monomerFuns[newMonomer.name];
        newMonomer.weight=monomerWeights[i]; // molecular weight in g/mol of the monomer
        if (newMonomer.weight>cSS.mwv.weightCutoff) {
            throw std::runtime_error("check weightCutoff for Mn and Mw calculation and compare to monomer weights. \n");
        }
        newMonomer.chain1=false;
        newMonomer.chain2=false;
        cSS.monomers[speciesList[i]]=newMonomer;
    }

    std::map<std::string,std::map<std::string,signed long>>::iterator it=monomerFuns.begin();
    std::map<std::string,signed long> numMoieties0; // initial total number of functional groups of each type

    // loop through each entry in the functional group map
    while (it != monomerFuns.end()) {
        std::map<std::string,signed long>::iterator it2=it->second.begin();
        while (it2 != it->second.end()) {
            int moietiesPerMolecule=it2->second;
            double fractionMoleculesOfTotal=initialConcentrations[it->first]/myConditions.totalMonomerConcentration;
            if (numMoieties0.find(it2->first)==numMoieties0.end()) {
                // moiety not found in list

                // replaced fractionMoleculesOfTotal*simulationSize  with numMonomers[it->first]
                numMoieties0[it2->first]=cSS.numMonomers[it->first]*moietiesPerMolecule; // monomer concentration
            }
            else {
                // moiety found in list, so add additional count to existing count
                //numMoieties0[it2->first]+=fractionMoleculesOfTotal*simulationSize*moietiesPerMolecule; // monomer concentration
                numMoieties0[it2->first]+=cSS.numMonomers[it->first]*moietiesPerMolecule; // monomer concentration

            }
            it2++;
        }
        it++;
    }

    // the above is actually a potential source of a bug. functional groups present are determined by what
    // initial monomers are present. due to rounding, the calculation with the fractionMoleculesOfTotal
    // may not get them quite right. 

    cSS.numMoieties=numMoieties0; // initialize number of moieties present at start of reaction

    // initialize numMoietiesByLength (which is of form std::vector<std::map<std::string,signed long>>)
    std::vector<std::map<std::string,signed long>> numMoietiesByLength(simulationSize+1,std::map<std::string,signed long>());
    // this code is specific to Changxia's system for polymerization of monomers -- if we start with chain, this next line would need to be edited
    numMoietiesByLength[1] = numMoieties0; // initialize number of moieties present at start of reaction (should all be monomer)
    cSS.numMoietiesByLength = numMoietiesByLength;
    // pass the following back to the primary kMC loop and analysis functions
    myConditions.numMoieties0=numMoieties0;
    myConditions.simulationSize=simulationSize;
    myConditions.monomerWeights=monomerWeights;
    myConditions.cSS0=cSS;
    myConditions.myReactions=reactionList;
    return myConditions;
}