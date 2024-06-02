#pragma once
#include<string>
#include<map>
#include"molecularWeight.h"

// declare timing typedefs, following example here
// https://stackoverflow.com/questions/36751133/proper-method-of-using-stdchrono
typedef std::chrono::steady_clock Time;
typedef std::chrono::milliseconds Millisecs;
typedef std::chrono::duration<double> DoubleSec;
typedef std::chrono::time_point<Time,DoubleSec> DoubleTimePoint;
typedef std::chrono::time_point<Time> TimePoint;

struct ModelOptions {
    // consider initializing all these, e.g., bool const trackLoopPairs=false;
    // then the lines at the start of the main function are unnecessary.
    bool const trackLoopPairs;
    bool const trackFullStructureAsGraph; // track the full explicit topology, or only functional groups and weights of each polymer chain
    bool const Tirrell; // compare to Tirrell sequence length theory (linear polymers only)
    bool const Changxia; // write out values to file to compare with Changxia's experimental data
    bool const printWeightFractions; // print wt fractions for comparison with most probable distribution. Assumes all monomers have the same molecular weight
    bool const eG; // write end group analysis to file
    bool const trackAllLinks;
    bool const loopingSeparateChannel; // is looping reaction a separate reaction channel, or should it just be lumped into the bimolecular reactions?
    bool const limitStep3; // should I limit step 3 in the hyperbranched polyester system?
};
    

struct Reaction {
    std::string reactionType; // "looping" "linking" "exchange" or "breaking"
    std::string step; // "step1r" "step 1f" "step2r" "step2f" etc.
    std::string number; // "R1", "R2", etc.
    std::vector<std::string> reactants; // one or two reactants expected
    std::vector<std::string> products; // one or two products expected
    std::string linkToBust=""; // linkage being broken 
    std::string linkFormed=""; // name of linkage formed in a linking reaction
    std::string gasEvolved="";
    std::string catalystEffect=""; // "bound" "released" -- for tracking effect of reaction on total catalyst count
    double Ea; // kJ/mol
    double A; // with units relative to the type of reaction (unimolecular, bimolecular, etc), s-1 * M^(-n+1) where n is the order of reaction
};

struct EdgeData {
    std::string name; // name of edge
    // store the linkage name and the monomer pair... I'm not convinced the monomer pair needs to be stored since it's stored in "source" and "target" of the edge
    // let's just store the linkages
    //std::string linkage; // this could be a vector of multiple linkage types possibly
};

struct Monomer {
    std::string name; // name of monomer or chain unit
    std::map<std::string,signed long> funs; // list of the functional groups
    double weight; // g/mol
    bool chain1;
    bool chain2;
    //std::vector<boost::adjacency_list<boost::vecS, boost::vecS,boost::bidirectionalS,monomer,EdgeData>::edge_descriptor> linkages;
};
typedef boost::adjacency_list<boost::vecS, boost::vecS,boost::bidirectionalS,Monomer,EdgeData> Graph;
typedef boost::adjacency_list<boost::vecS, boost::vecS,boost::undirectedS,Monomer,EdgeData> UndirectedGraph;
typedef Graph::edge_descriptor Edge;
typedef Graph::vertex_descriptor Vertex;
typedef std::map<Vertex, size_t> IndexMap;


// class VertexIndexMap //maps vertex to index
// {
// public:
//     typedef boost::readable_property_map_tag category; 
//     typedef size_t  value_type;
//     typedef value_type reference;
//     typedef Vertex key_type; 

//     VertexIndexMap(const Graph& g): _g(&g) {} 

//     const Graph * _g;
// };

// namespace boost {

// template<>
// struct property_map<Graph, vertex_index_t > {
//     typedef VertexIndexMap const_type;
//     //typedef type const_type ; //-- we do not define type as "vertex_index_t" map is read-only 
// };

// VertexIndexMap get(boost::vertex_index_t, const Graph & g )
// {
//     return VertexIndexMap(g);
// }

// VertexIndexMap::value_type get(VertexIndexMap map, VertexIndexMap::key_type vertex)
// {
//     // bruh i have no idea how this works, but i guess i need a get statement that returns all the shit in a Monomer
//     //Monomer myMonomer = (*map._g)[vertex].id
//     return (*map._g)[vertex].id;
// }

// }


struct Chain {
    Graph polymer; // store topological connectivity between monomers
    std::map<std::string,signed long> funs; // list of the functional groups and the quantity that exist in a given chain... 
    double weight; // g/mol
    std::map<std::string,signed long> linkTotals; // specifically for changxia's system
    boost::property_map<Graph, std::string EdgeData::*>::type nmp; // namemap
};

struct SystemVariables {
    std::vector<Chain> chainPool; // all chains in system
    std::map<std::string,long> numMoieties; // current number of functional groups of each type within the control volume
    std::vector<std::map<std::string,signed long>> numMoietiesByLength; // track the number of moieties at each length (as you cannot have a chain or monomer of length 0, the 0th entry should always be empty)
    std::map<std::string, signed long> numMonomers; // number of each tracked monomer
    std::map<std::string,Monomer> monomers; // what monomers are present
    std::map<std::string,std::vector<std::string>> funGroups; // map of functional groups to the monomers on which they appear
    std::map<std::string, signed long> links; // map to store the number of links of each type -- really only working for Changxia's system at the moment
    // std::map<std::string,signed long long> loopingPairs; // map from pair name to number present -- also really only working for Changxia's system...
    std::map< std::string , std::map < std::string , signed long long > > loopPairs; // map from pair name to number present -- symmetric 2D map
    MolWtVariables mwv; // store current molecular weight and associated trackers
    int executionNumber=0; // track how many kMC steps have been performed
    std::vector<Reaction> loopingRxns;
    unsigned int seed;
    std::map<std::string,long> moleculesOfGasEvolved; // how much gas has been evolved over the course of the reaction?
    int rejectedSteps=0;
    int linksWithCatalyst=0;
    int minLoopSize = 4; // set minimum loop size
    int dmax = -1; // set max topological distance to explicitly consider in looping reactions
    std::vector<long> reactionsExecutedTally; // track how many times each reaction has been executed
    long maxChainSize=0;
};

long monomerFunGrpCount(std::map<std::string,std::vector<std::string>> funGroups, std::string grp, 
std::map<std::string, signed long> numMonomers, std::map<std::string,Monomer> monomers);

struct Connectivity {
    bool isSingleComponent; // is the graph a single component?
    std::vector<int> component;
};

Connectivity checkConnectivity(Graph d);

struct EdgeParticipants {
    // contains info on edge selected to partake in a reaction, which can be returned from a selecting function
    Vertex src;
    Vertex tgt;
    signed long long chainIndex;
};

typedef EdgeParticipants LoopParticipants;

struct GroupLocation {
    // contains info on group selected to react, which can be returned from a selecting function
    int whichChain=-1; // index in chainPool
    std::string monomerSelected=""; // if not in chainPool, type of monomer to react
    Vertex whichVertex; // only if tracking full structure as a graph boost::range_detail::integer_iterator<unsigned long>
};

EdgeParticipants linearSearchForMatchingEdge(std::vector<Chain> chainPool,double 
    whichLink, std::string linkName);
GroupLocation linearSearchForParticipatingGroup(const SystemVariables &cSS,
    const double &randomNumber, std::string groupName,  bool trackFullStructureAsGraph);

template < typename edgeIter>
void updateLinks(std::map<std::string, signed long> &links,Chain &myChain,const std::string &oldVertName, 
const std::string &newVertName, std::pair<edgeIter,edgeIter> edgs, const Edge &newEdge);

void updateAdjacentLinks(std::map<std::string, signed long> &links,Chain &myChain, Reaction rxn, EdgeParticipants modifiedEdge,const Edge &newEdge);

struct chosenReactionChannel {
    bool isLengthDependent;
    int mu;
    int length=0;
    std::string chosenR_i="";
    std::string reactionType="";
};

void updateSystemState(SystemVariables &cSS, Reaction rxn,double r3,double r4,ModelOptions myOptions, std::map<std::string,DoubleSec>& timers ,chosenReactionChannel chosenRxn,std::vector<Reaction> reactionList,int length, double ratesByLength_lengthRi,std::vector<double> transport_factor,std::vector<double> ksteptimesstar);

bool isThreeMemberedRing(LoopParticipants pairToLoop,Graph polymer);

struct PossibleLoopingRxn {
    std::string catAlcoholName; // name of the alcohol-catalyst complex
    long d; // topological distance between the two vertices
    long p; // index of chain in chainPool
    std::pair<Vertex,Vertex> vertices; // std::pair vertices (if less than dmax). Mark second one as -1 if not known
};

double polymerChainVolumeInLiters(PossibleLoopingRxn candidateLoop, int dmax);

std::string getChainUnitName(Monomer myMonomer);

std::vector<PossibleLoopingRxn> enumeratePossibleLoops(SystemVariables &cSS);
void loopingRxn(SystemVariables &cSS, PossibleLoopingRxn lr, ModelOptions myOptions, std::vector<Reaction> reactionList);
