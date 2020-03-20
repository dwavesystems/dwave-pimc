/*
Copyright 2020 D-Wave Systems Inc.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

/*
   Path Integral Monte Carlo code for analysis of a finite temperature transverse field Ising model,
   Defined by partition function Z = Trace[exp(-H)], and scaled Hamiltonian
   H = invTemp/2 sum_{i,j} J_{ij} \sigma^z_i \sigma^z_j + invTemp \sum_i [h_i \sigma^z_i - \Gamma\sigma^x_i]

   Methods exploit either single qubit Swendsen-Wang updates, or multi-qubit Swendsen-Wang updates,
   in the latter case specifically for regular and independent 1d ferromagnetic subsequences. These
   match the methods explored in A King et al. https://arxiv.org/abs/1911.03446

   Authors: Jack Raymond, Stephen Face
   Copyright: D-Wave Systems
   License: Apache 2
   Last modification: March 20 2020

   See also README.md localPIMC.hpp and main.cpp
*/

#include "localPIMC.hpp"
/* Public functions */

localPIMC::localPIMC(int Lperiodic, double invTempOverJ, double GammaOverJ, int initialCondition, int qubitsPerUpdate0,
                     int qubitsPerChain0, unsigned int seed) {
    // Constructor parameterized to match paper experiments:
    qubitsPerUpdate = qubitsPerUpdate0;
    qubitsPerChain = qubitsPerChain0;
    assert(qubitsPerChain >= qubitsPerUpdate);
    numVar = Lperiodic * (Lperiodic / 2 + 3) * qubitsPerChain;
    constructCouplingMatrix(Lperiodic, invTempOverJ);
    initializeWorldLines(initialCondition, Lperiodic, qubitsPerChain);
    const double Jchain = -1.8;
    invTempH = std::vector<double>(numVar, 0);
    invTempJchain = Jchain * invTempOverJ;
    invTempGamma = GammaOverJ * invTempOverJ;
    initPRNG(seed);
}
localPIMC::localPIMC(double Gamma, double Jchain, int qubitsPerUpdate0, int qubitsPerChain0,
                     std::vector<std::vector<int> > adjMat0, std::vector<std::vector<double> > invTempJ0,
                     std::vector<double> invTempH0, std::vector<int> classicalInitialCondition, unsigned int seed) {
    // Generic constructor
    qubitsPerUpdate = qubitsPerUpdate0;
    qubitsPerChain = qubitsPerChain0;
    assert(qubitsPerChain >= qubitsPerUpdate);
    numVar = adjMat0.size();
    adjMat.swap(adjMat0);
    invTempJ.swap(invTempJ0);
    invTempH.swap(invTempH0);
    invTempJchain = Jchain;
    invTempGamma = Gamma;
    breaks.resize(numVar);
    firstSlice.swap(classicalInitialCondition);
    initPRNG(seed);
};
void localPIMC::run(int nSweeps) {
    if (qubitsPerUpdate == 1) {
        std::uniform_int_distribution<> randomQubitIndex(0, numVar - 1);
        for (int sweepI = 0; sweepI < nSweeps * numVar; sweepI++) {
            qubitUpdate(randomQubitIndex(prng));
        }
    } else {
        int numChains = numVar / qubitsPerChain;
        std::uniform_int_distribution<> randomChainIndex(0, numChains - 1);
        for (int sweepI = 0; sweepI < nSweeps * numChains; sweepI++) {
            chainUpdate(randomChainIndex(prng));
        }
    }
}

/* Private functions */

void localPIMC::constructCouplingMatrix(int Lperiodic, double invTemp) {
    // Construct sparse coupling matrix for cylindrical lattices used in paper, either triangular or square
    // octagonal.
    assert(Lperiodic % 6 == 0);
    int Lopen = (Lperiodic / 2 + 3);
    assert(qubitsPerChain == 1 || qubitsPerChain == 4);
    int disp_i[] = {0, 1, 1}, disp_j[] = {1, 0, 1};
    int chainConnectionFrom[] = {(qubitsPerChain > 1 ? qubitsPerChain - 2 : 0), qubitsPerChain - 1, qubitsPerChain - 1},
        chainConnectionTo[] = {0, (qubitsPerChain > 1 ? 1 : 0), 0};
    invTempJ.resize(numVar);
    adjMat.resize(numVar);
    for (int i0 = 0; i0 < Lperiodic; i0++) {
        for (int j0 = 0; j0 < Lopen; j0++) {
            for (int couplerOrientationI = 0; couplerOrientationI < 3; couplerOrientationI++) {
                int i1 = (i0 + disp_i[couplerOrientationI]) % Lperiodic;
                int j1 = j0 + disp_j[couplerOrientationI];
                if (j1 < Lopen) {
                    int linearIndex0 = (i0 * Lopen + j0) * qubitsPerChain + chainConnectionFrom[couplerOrientationI];
                    int linearIndex1 = (i1 * Lopen + j1) * qubitsPerChain + chainConnectionTo[couplerOrientationI];
                    adjMat[linearIndex1].push_back(linearIndex0);
                    adjMat[linearIndex0].push_back(linearIndex1);

                    if ((j0 % (Lopen - 1)) || (j1 % (Lopen - 1))) {
                        invTempJ[linearIndex1].push_back(invTemp);
                        invTempJ[linearIndex0].push_back(invTemp);
                    } else {
                        // Cylindrical boundary condition:
                        invTempJ[linearIndex1].push_back(invTemp / 2);
                        invTempJ[linearIndex0].push_back(invTemp / 2);
                    }
                }
            }
        }
    }
    if (qubitsPerUpdate == 1 && qubitsPerChain > 1) {
        // Enumerate also couplings within chain if qubitUpdate instead of chainUpdate
        for (int n = 0; n < numVar; n += qubitsPerChain) {
            for (int k = 0; k < qubitsPerChain - 1; k++) {
                adjMat[n + k].push_back(n + k + 1);
                adjMat[n + k + 1].push_back(n + k);
                invTempJ[n + k].push_back(invTempJchain);
                invTempJ[n + k + 1].push_back(invTempJchain);
            }
        }
    }
}

void localPIMC::initializeWorldLines(int initialCondition, int Lperiodic, int qubitsPerChain) {
    // Initial state of Markov Chain, clockwise (1), counterclockwise (-1) or unwound (0)
    // Periodicity required in vertical direction (cylindrical lattice), dimension is multiple of 6
    assert(Lperiodic % 6 == 0 && Lperiodic >= 6);
    assert(initialCondition <= 1 && -1 <= initialCondition);
    int Lopen = 3 * (Lperiodic / 6 + 1);
    numVar = Lperiodic * Lopen * qubitsPerChain;
    // See paper description, 6 blocks, each block uses a different pseudoSpin orientation:
    int alignedMask[] = {1, 1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, 1, 1, -1, -1};  // Check paper!
    std::vector<int> tripartiteClassification = makeTripartiteClassification(Lperiodic);
    int blockSize = numVar / (6 * qubitsPerChain);
    int blockState = 0;
    int nChains = Lperiodic * Lopen;
    firstSlice.resize(numVar);
    breaks.resize(numVar);
    for (int n = 0; n < nChains; n++) {
        breaks[n].resize(0);
        assert(tripartiteClassification.size() == nChains);
        assert(blockState * 3 + tripartiteClassification[n] < 18);
        int thisSpin = alignedMask[blockState * 3 + tripartiteClassification[n]];
        for (int k = 0; k < qubitsPerChain; k++) {
            firstSlice[qubitsPerChain * n + k] = thisSpin;
        }
        if (n % blockSize == blockSize - 1) {
            // ordered (initialCondition=0) leaves pseudoSpin unchanged, otherwise rotated
            blockState = (blockState + initialCondition + 6) % 6;
        }
    }
}
void localPIMC::addHToEffectiveField(std::vector<double>& effectiveFieldPerDomain,
                                     const std::vector<int>& allInterfaces, double H) const {
    H /= numTrotterSlices;
    effectiveFieldPerDomain[0] += (numTrotterSlices + allInterfaces.front() - allInterfaces.back()) * H;
    for (int interfaceI = 1; interfaceI < allInterfaces.size(); interfaceI++) {
        effectiveFieldPerDomain[interfaceI] += (allInterfaces[interfaceI] - allInterfaces[interfaceI - 1]) * H;
    }
}

void localPIMC::addHToEffectiveField(std::vector<double>& effectiveFieldPerDomain,
                                     const std::vector<int>& componentLabels, int componentOffset,
                                     const std::vector<int>& allInterfaces, double H) const {
    H /= numTrotterSlices;
    effectiveFieldPerDomain[componentLabels[componentOffset]] +=
        (numTrotterSlices + allInterfaces.front() - allInterfaces.back()) * H;
    for (int interfaceI = 1; interfaceI < allInterfaces.size(); interfaceI++) {
        effectiveFieldPerDomain[componentLabels[componentOffset + interfaceI]] +=
            (allInterfaces[interfaceI] - allInterfaces[interfaceI - 1]) * H;
    }
}

void localPIMC::addJToEffectiveField(std::vector<double>& effectiveFieldPerDomain,
                                     const std::vector<int>& allInterfaces, int neighbor, double Js) const {
    Js /= numTrotterSlices;
    std::vector<int> allInterfacesPair(allInterfaces.size() + breaks[neighbor].size());
    std::merge(breaks[neighbor].begin(), breaks[neighbor].end(), allInterfaces.begin(), allInterfaces.end(),
               allInterfacesPair.begin());
    effectiveFieldPerDomain[0] += (numTrotterSlices + allInterfacesPair.front() - allInterfacesPair.back()) * Js;
    int interfaceI = 0;
    int interfaceAllI = 0;
    for (; interfaceAllI < allInterfacesPair.size() - 1; interfaceAllI++) {
        if (allInterfacesPair[interfaceAllI] == allInterfaces[interfaceI]) {
            if (++interfaceI == allInterfaces.size()) break;
        } else {
            Js *= -1;
        }
        effectiveFieldPerDomain[interfaceI] +=
            (allInterfacesPair[interfaceAllI + 1] - allInterfacesPair[interfaceAllI]) * Js;
    }
    for (; interfaceAllI < allInterfacesPair.size() - 1; interfaceAllI++) {
        effectiveFieldPerDomain[0] += (allInterfacesPair[interfaceAllI + 1] - allInterfacesPair[interfaceAllI]) * Js;
        Js *= -1;
    }
}
void localPIMC::addJToEffectiveField(std::vector<double>& effectiveFieldPerDomain,
                                     const std::vector<int>& componentLabels, int componentOffset,
                                     const std::vector<int>& allInterfaces, int neighbor, double Js) const {
    Js /= numTrotterSlices;
    std::vector<int> allInterfacesPair(allInterfaces.size() + breaks[neighbor].size());
    std::merge(breaks[neighbor].begin(), breaks[neighbor].end(), allInterfaces.begin(), allInterfaces.end(),
               allInterfacesPair.begin());
    effectiveFieldPerDomain[componentLabels[componentOffset]] +=
        (numTrotterSlices + allInterfacesPair.front() - allInterfacesPair.back()) * Js;
    int interfaceI = 0;
    int interfaceAllI = 0;
    for (; interfaceAllI < allInterfacesPair.size() - 1; interfaceAllI++) {
        if (allInterfacesPair[interfaceAllI] == allInterfaces[interfaceI]) {
            if (++interfaceI == allInterfaces.size()) break;
        } else {
            Js *= -1;
        }
        effectiveFieldPerDomain[componentLabels[componentOffset + interfaceI]] +=
            (allInterfacesPair[interfaceAllI + 1] - allInterfacesPair[interfaceAllI]) * Js;
    }
    for (; interfaceAllI < allInterfacesPair.size() - 1; interfaceAllI++) {
        effectiveFieldPerDomain[componentLabels[componentOffset]] +=
            (allInterfacesPair[interfaceAllI + 1] - allInterfacesPair[interfaceAllI]) * Js;
        Js *= -1;
    }
}

void localPIMC::qubitUpdate(int sp) {
    // Update 1 qubit system H = \nu_i \sigma^z_i + invTempGamma \sigma^x_i, where \nu_i are determined from
    // neighboring qubit states
    std::vector<int> breakProposals = makeBreakProposals();
    std::vector<int> allInterfaces;
    if (breakProposals.size() + breaks[sp].size() > 1) {
        allInterfaces.resize(breakProposals.size() + breaks[sp].size());
        std::merge(breakProposals.begin(), breakProposals.end(), breaks[sp].begin(), breaks[sp].end(),
                   allInterfaces.begin());
    } else
        allInterfaces.push_back(numTrotterSlices);  // Book keeping, less branching

    int nDom = allInterfaces.size();
    std::vector<double> effectiveFieldPerDomain(nDom, 0);
    addHToEffectiveField(effectiveFieldPerDomain, allInterfaces, invTempH[sp]);
    for (int neighI = 0; neighI < adjMat[sp].size(); neighI++) {
        int neighbor = adjMat[sp][neighI];
        addJToEffectiveField(effectiveFieldPerDomain, allInterfaces, neighbor,
                             firstSlice[neighbor] * invTempJ[sp][neighI]);
    }
    // Sample spin on boundary spanning domain
    int sValue = GibbsSamplePM1(effectiveFieldPerDomain[0]);
    firstSlice[sp] = sValue;
    breaks[sp].resize(0);
    // Sample states
    for (int interfaceI = 1; interfaceI < allInterfaces.size(); interfaceI++) {
        // Sample spin on subsequent domains, if flipped record interface:
        if (sValue * GibbsSamplePM1(effectiveFieldPerDomain[interfaceI]) != 1) {
            sValue *= -1;
            breaks[sp].push_back(allInterfaces[interfaceI - 1]);
        }
    }
    if (sValue != firstSlice[sp]) breaks[sp].push_back(allInterfaces.back());
}
//
void localPIMC::chainUpdate(int sp) {
    // Update multi-qubits system H = invTempJchain \sum_i \sigma^z_i \sigma^z_{i+1} + \sum_i [invTemp_i \sigma^x_i
    // + \nu_i \sigma^Z_i],
    // by Swendsen-Wang where \nu_i are determined from neighboring qubit states
    std::vector<std::vector<int> > allInterfacesEveryChain(qubitsPerChain);
    std::vector<int> nDomains(qubitsPerChain, 0);
    int nDomTotal = 0;
    for (int chainI = 0; chainI < qubitsPerChain; chainI++) {
        int qubitI = qubitsPerChain * sp + chainI;
        std::vector<int> breakProposals = makeBreakProposals();
        // Joint list of interfaces with
        if (breakProposals.size() + breaks[qubitI].size() > 1) {
            allInterfacesEveryChain[chainI].resize(breakProposals.size() + breaks[qubitI].size());
            std::merge(breakProposals.begin(), breakProposals.end(), breaks[qubitI].begin(), breaks[qubitI].end(),
                       allInterfacesEveryChain[chainI].begin());
        } else {
            allInterfacesEveryChain[chainI].push_back(numTrotterSlices);  // Book keeping, less branching later
        }
        nDomains[chainI] = nDomTotal;
        nDomTotal += allInterfacesEveryChain[chainI].size();
    }
    // Build connectivity of qubit level domains, by Swendsen-Wang rule
    std::vector<std::vector<int> > domainGraph(nDomTotal);
    for (int chainI = 0; chainI < qubitsPerChain - 1; chainI++) {
        makeDomainGraph(nDomains[chainI], nDomains[chainI + 1], sp, chainI, allInterfacesEveryChain, domainGraph);
    }
    // Find components
    int nComponents = 0;
    std::vector<int> componentLabels(nDomTotal, -1);
    for (int root = 0; root < componentLabels.size(); root++) {
        if (componentLabels[root] == -1) depthFirstDomainAssignment(domainGraph, componentLabels, nComponents++, root);
    }
    // Consolidate external fields
    std::vector<double> effectiveFieldAll(nComponents, 0);
    for (int chainI = 0; chainI < qubitsPerChain; chainI++) {
        int qubitI = qubitsPerChain * sp + chainI;
        addHToEffectiveField(effectiveFieldAll, componentLabels, nDomains[chainI], allInterfacesEveryChain[chainI],
                             invTempH[qubitI]);  // OPTIONAL FOR LATTICE CODE:
        for (int neighI = 0; neighI < adjMat[qubitI].size(); neighI++) {
            int neighbor = adjMat[qubitI][neighI];
            addJToEffectiveField(effectiveFieldAll, componentLabels, nDomains[chainI], allInterfacesEveryChain[chainI],
                                 neighbor, firstSlice[neighbor] * invTempJ[qubitI][neighI]);
        }
    }
    // sample domain spin according to effective fields
    std::vector<int> sValues(nComponents, 0);
    for (int componentLabel = 0; componentLabel < nComponents; componentLabel++) {
        sValues[componentLabel] = GibbsSamplePM1(effectiveFieldAll[componentLabel]);
    }
    // Map back assignment to each qubit and record interfaces
    for (int chainI = 0; chainI < qubitsPerChain; chainI++) {
        int qubitI = qubitsPerChain * sp + chainI;
        int dom = componentLabels[nDomains[chainI]];
        firstSlice[qubitI] = sValues[dom];
        breaks[qubitI].resize(0);
        int sValue = firstSlice[qubitI];
        for (int domainI = 1; domainI < allInterfacesEveryChain[chainI].size(); domainI++) {
            dom = componentLabels[nDomains[chainI] + domainI];
            if (sValue * sValues[dom] != 1) {
                sValue *= -1;
                breaks[qubitI].push_back(allInterfacesEveryChain[chainI][domainI - 1]);
            }
        }
        if (sValue != firstSlice[qubitI]) {
            breaks[qubitI].push_back(allInterfacesEveryChain[chainI].back());
        }
    }
}

void localPIMC::depthFirstDomainAssignment(const std::vector<std::vector<int> >& domainGraph,
                                           std::vector<int>& componentLabels, int componentLabel, int root) const {
    componentLabels[root] = componentLabel;
    for (int leaf = 0; leaf < domainGraph[root].size(); leaf++)
        if (componentLabels[domainGraph[root][leaf]] == -1)
            depthFirstDomainAssignment(domainGraph, componentLabels, componentLabel, domainGraph[root][leaf]);
}

void localPIMC::makeDomainGraph(int zeroChainIndex, int firstChainIndex, int sp, int chainI,
                                const std::vector<std::vector<int> >& allInterfacesEveryChain,
                                std::vector<std::vector<int> >& domainGraph) const {
    // If unfrustrated, attempt merge as function of interface size:
    std::uniform_real_distribution<double> probability(
        0.0, 1.0);  // Check if this is slow.. norm by max value would be good enough.
    int qubitI = sp * qubitsPerChain + chainI;
    int s0s1 = firstSlice[qubitI] * firstSlice[qubitI + 1];
    int nProposedInterfaces = allInterfacesEveryChain[chainI].size() + allInterfacesEveryChain[chainI + 1].size();
    std::vector<int> allInterfacesPair(nProposedInterfaces);
    std::merge(allInterfacesEveryChain[chainI].begin(), allInterfacesEveryChain[chainI].end(),
               allInterfacesEveryChain[chainI + 1].begin(), allInterfacesEveryChain[chainI + 1].end(),
               allInterfacesPair.begin());
    if (allInterfacesPair[0] == numTrotterSlices) {
        // No non-trivial domains, simple join:
        if (1 == s0s1 && pNotJoin(numTrotterSlices) < probability(prng)) {
            domainGraph[zeroChainIndex].push_back(firstChainIndex);
            domainGraph[firstChainIndex].push_back(zeroChainIndex);
        }
    } else {
        // At most one book keeping domain,remove:
        if (allInterfacesPair.back() == numTrotterSlices) {
            // Remove unnecessary book keeping values
            allInterfacesPair.pop_back();
        }
        // Special handling of boundary spanning domain (periodic in imaginary time):
        if (1 == s0s1 &&
            pNotJoin(numTrotterSlices - allInterfacesPair.back() + allInterfacesPair.front()) < probability(prng)) {
            domainGraph[zeroChainIndex].push_back(firstChainIndex);
            domainGraph[firstChainIndex].push_back(zeroChainIndex);
        }

        int posProposal0 = 0, posExisting0 = 0, posProposal1 = 0, posExisting1 = 0,
            chain0valid = (allInterfacesEveryChain[chainI][0] != numTrotterSlices);
        int interfaceI = 0;
        while (++interfaceI < allInterfacesPair.size()) {
            if (chain0valid && allInterfacesPair[interfaceI - 1] == allInterfacesEveryChain[chainI][posProposal0]) {
                if (allInterfacesEveryChain[chainI].size() == ++posProposal0) {
                    chain0valid = 0;
                    posProposal0 = 0;
                }
                if (breaks[qubitI].size() > posExisting0 &&
                    allInterfacesPair[interfaceI - 1] == breaks[qubitI][posExisting0]) {
                    s0s1 = s0s1 * (-1);
                    posExisting0++;  // Part of existing boundary
                }
            } else {
                if (allInterfacesEveryChain[chainI + 1].size() == ++posProposal1) {
                    posProposal1 = 0;
                }
                if (breaks[qubitI + 1].size() > posExisting1 &&
                    allInterfacesPair[interfaceI - 1] == breaks[qubitI + 1][posExisting1]) {
                    s0s1 = s0s1 * (-1);
                    posExisting1++;  // Part of existing boundary
                }
            }
            if (s0s1 == 1 &&
                pNotJoin(allInterfacesPair[interfaceI] - allInterfacesPair[interfaceI - 1]) < probability(prng)) {
                domainGraph[zeroChainIndex + posProposal0].push_back(firstChainIndex + posProposal1);
                domainGraph[firstChainIndex + posProposal1].push_back(zeroChainIndex + posProposal0);
            }
        }
    }
}

std::vector<int> localPIMC::makeBreakProposals() const {
    // Sampling of conditional distribution, exploiting continuous limit invTempGamma << numTrotterSlices:
    std::vector<int> breakProposals;
    std::uniform_real_distribution<double> probability(0.0, 1.0);
    double nInterfacesScale = numTrotterSlices / invTempGamma;
    double position = -nInterfacesScale * log(probability(prng));
    while (position < numTrotterSlices) {
        breakProposals.push_back(int(position));
        position = -nInterfacesScale * log(probability(prng)) + breakProposals.back() + 1;
    }
    return breakProposals;
}

std::vector<int> localPIMC::makeTripartiteClassification(int Lperiodic) const {
    int Lopen = 3 * (Lperiodic / 6 + 1);
    std::vector<int> tripartiteClassification(Lperiodic * Lopen);
    for (int i = 0; i < Lperiodic; i++) {
        for (int j = 0; j < Lopen; j++) {
            tripartiteClassification[i * Lopen + j] = (i + j) % 3;
        }
    }
    return tripartiteClassification;
}

int localPIMC::GibbsSamplePM1(double effectiveField) const {
    // Sample from P(s) ~ exp( - effField *s) = (1 - s tanh(effField))/2, effField is the energy associated to state
    // +1 scaled to beta=1.
    std::uniform_real_distribution<double> magThresh(-1.0, 1.0);
    return tanh(effectiveField) > magThresh(prng) ? -1 : 1;
}

double localPIMC::pNotJoin(int nOverlaps) const { return exp((2 * invTempJchain * nOverlaps) / numTrotterSlices); }

void localPIMC::initPRNG(unsigned int seed) const {
    if (seed) {
        prng = std::mt19937(seed);
    } else {
        // seed = 0 signals use of random device
        std::random_device r;
        std::seed_seq seedSeq{r(), r(), r(), r(), r(), r(), r(), r()};
        prng = std::mt19937(seedSeq);
    }
}


