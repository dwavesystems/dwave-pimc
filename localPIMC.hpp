/*
    Apache licence here
*/

/*
   Path Integral Monte Carlo code for analysis of a finite temperature transverse field Ising model,
   Defined by partition function Z = Trace[exp(-H)], and scaled Hamiltonian
   H = invTemp/2 sum_{i,j} J_{ij} \sigma^z_i \sigma^z_j + invTemp \sum_i [h_i \sigma^z_i - \Gamma\sigma^x_i]
   //NB: J is symmetric, negative J denotes a ferromagnetic interaction

   Methods exploits either single qubit Swendsen-Wang updates, or multi-qubit Swendsen-Wang updates,
   in the latter case specifically for regular and independent 1d ferromagnetic subsequences. These
   match the methods explored in A King et al. https://arxiv.org/abs/1911.03446

   Compile as: g++ -O3 -Wall -Wextra -Wno-sign-compare -std=c++11 main.cpp localPIMC.cpp -o _demo;
   To generate projected samples following 32000 sweeps from the initial condition [Figure 17 1911.03446v2]
   (convergence on 24 x 15 triangular cylinder, single qubit move) run the demo as:
   ./_demo 0 IC; %IC == -1,0,1 for counterclockwise, ordered, and clockwise initial conditions
   (convergence on 12 x 9 triangular cylinder, single qubit move) run the demo as:
   ./_demo 1 IC; %IC == -1,0,1 for counterclockwise, ordered, and clockwise initial conditions
   To generate projected samples following 32000 sweeps from the initial condition [Figure 15f 1911.03446v2]
   (convergence on square octagonal lattice, 4 qubit chain moves) run the demo as:
   ./_demo 2 IC; %IC == -1,0,1 for counterclockwise, ordered, and clockwise initial conditions
   (convergence on square octagonal lattice, 1 qubit chain moves [faster per sweep, but slower in equilibration])
   run the demo as:
   ./_demo 3 IC; %IC == -1,0,1 for counterclockwise, ordered, and clockwise initial conditions

   To reproduce other square octagonal lattice results call the localPIMC class lattice constructor
   (invTempOverJ, GammaOverJ, L,qubitsPerChain=4) and algorithm (qubitsPerChain) parameters.
   For generic models use the alternative constructor to specify the model structure.

   Code is for demonstration purposes, with full functionality and some optimizations omitted. 

   Authors: Jack Raymond, Stephen Face
   Copyright: D-Wave Systems
   License: Apache
   Last modification: March 13 2020

   See also localPIMC.cpp and main.cpp
*/
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <random>

class localPIMC {
   private:
    // Hard coding to practically convenient and statistically safe value.
    const int numTrotterSlices = 1 << 16;
    mutable std::mt19937 prng;
    int qubitsPerChain;
    int numVar;
    // symmetric coupling matrix, and associated scaled J values (excluding ferromagnetic part if chain update method,
    // including ferromagnetic part if 1 qubit method)
    std::vector<std::vector<int> > adjMat;
    std::vector<std::vector<double> > invTempJ;
    std::vector<double> invTempH;  // Scaled external longitudinal field
    double invTempGamma;           // Scaled transverse field
    double invTempJchain;          // Scaled coupling strength along chains
    int qubitsPerUpdate;           // 1 or 4 in square octagonal lattice, 1 in triangular lattice.
    void constructCouplingMatrix(int Lperiodic, double invTemp);
    void initializeWorldLines(int initialCondition, int Lperiodic, int qubitsPerChain);
    void addHToEffectiveField(std::vector<double>& effectiveFieldPerDomain, const std::vector<int>& allInterfaces,
                              double H) const;
    void addHToEffectiveField(std::vector<double>& effectiveFieldPerDomain, const std::vector<int>& componentLabels,
                              int componentOffset, const std::vector<int>& allInterfaces, double H) const;
    void addJToEffectiveField(std::vector<double>& effectiveFieldPerDomain, const std::vector<int>& allInterfaces,
                              int neighbor, double Js) const;
    void addJToEffectiveField(std::vector<double>& effectiveFieldPerDomain, const std::vector<int>& componentLabels,
                              int componentOffset, const std::vector<int>& allInterfaces, int neighbor,
                              double Js) const;
    void makeDomainGraph(int zeroChainIndex, int firstChainIndex, int sp, int chainI,
                         const std::vector<std::vector<int> >& allInterfacesEveryChain,
                         std::vector<std::vector<int> >& domainGraph) const;
    void qubitUpdate(int sp);
    void chainUpdate(int sp);
    void depthFirstDomainAssignment(const std::vector<std::vector<int> >& domainGraph,
                                    std::vector<int>& componentLabels, int componentLabel, int root) const;
    std::vector<int> makeBreakProposals() const;
    int GibbsSamplePM1(double effectiveField) const;
    double pNotJoin(int nOverlaps) const;
    void initPRNG(unsigned int seed) const;

   public:
    //Worldline specification:
    std::vector<std::vector<int> > breaks;  // Per qubit interfaces in imaginary time [0,numTrotterSlices), even.
    std::vector<int> firstSlice;            // Values of qubits on boundary spanning domain.
    localPIMC(int Lperiodic, double invTempOverJ, double GammaOverJ, int initialCondition, int qubitsPerUpdate0,
              int qubitsPerChain0, unsigned int seed);
    localPIMC(double Gamma, double Jchain, int qubitsPerUpdate0, int qubitsPerChain0,
              std::vector<std::vector<int> > adjMat0, std::vector<std::vector<double> > invTempJ0,
              std::vector<double> invTempH0, std::vector<int> classicalInitialCondition, unsigned int seed);
    std::vector<int> makeTripartiteClassification(int Lperiodic) const;
    void run(int nSweeps);
};
