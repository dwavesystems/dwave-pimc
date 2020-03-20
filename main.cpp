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

   See also README.md localPIMC.cpp and localPIMC.hpp
*/

#include <ctime>
#include <iostream>
#include "localPIMC.hpp"

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Requires two arguments: experimentNo[0,1] initialCondition[-1,0,1]\n e.g. ./_demo 1 1\n";
        return EXIT_FAILURE;
    }
    int experimentNo = atoi(argv[1]);
    int initialCondition = atoi(argv[2]);
    if (initialCondition < -1 || initialCondition > 1) {
        std::cerr << "Second argument, initialCondition must be selected from {-1,0,1}";
        return EXIT_FAILURE;
    }
    unsigned int seed = 0;  // 0 = Use random device, otherwise reproducible
    int Lperiodic = 24;     // Largest lattice
    int nSweeps = 32768;

    switch (experimentNo) {
        case 0: {
            //[Figure 17 1911.03446v2] (convergence on triangular cylinder, single qubit move)
            int qubitsPerChain = 1;
            int qubitsPerUpdate = 1;
            double invTempOverJ = 2;
            double GammaOverJ = 0.6;
            std::cout << "#24 by 15 Triangular lattice, T/J=0.5, Gamma/J=0.6, generate sample projected in "
                         "computational basis at "
                      << nSweeps << " sweeps, from classical initial condition with winding number " << initialCondition
                      << "\n";

            localPIMC localpimc(Lperiodic, invTempOverJ, GammaOverJ, initialCondition, qubitsPerUpdate, qubitsPerChain,
                                seed);

            auto tstart = std::clock();
            localpimc.run(nSweeps);
            auto tend = std::clock();
            std::cout << "#Time required:" << (tend - tstart) / (double)CLOCKS_PER_SEC << " second(s)\n";
            for (int x : localpimc.firstSlice) {
                std::cout << x << " ";
            }
            break;
        }
        case 1: {
            //[Figure 17 1911.03446v2] (convergence on triangular cylinder, single qubit move)
            Lperiodic = 12;
            int qubitsPerChain = 1;
            int qubitsPerUpdate = 1;
            double invTempOverJ = 2;
            double GammaOverJ = 0.6;
            std::cout << "#12 by 9 Triangular lattice, T/J=0.5, Gamma/J=0.6, generate sample projected in "
                         "computational basis at "
                      << nSweeps << " sweeps, from classical initial condition with winding number " << initialCondition
                      << "\n";

            localPIMC localpimc(Lperiodic, invTempOverJ, GammaOverJ, initialCondition, qubitsPerUpdate, qubitsPerChain,
                                seed);

            auto tstart = std::clock();
            localpimc.run(nSweeps);
            auto tend = std::clock();
            std::cout << "#Time required:" << (tend - tstart) / (double)CLOCKS_PER_SEC << " second(s)\n";
            for (int x : localpimc.firstSlice) {
                std::cout << x << " ";
            }
            std::cout << "\n";
            break;
        }
        case 2: {
            //[Figure 15f and others, 1911.03446v2] (convergence on square octagonal lattice cylinder, four qubit move)
            //[Various figures, 1911.03446v2] (convergence on square octagonal lattice cylinder, single qubit move)
            int qubitsPerChain = 4;
            int qubitsPerUpdate = 4;
            double invTempOverJ = 1 / 0.244;
            double GammaOverJ = 0.736;
            std::cout << "#24 by 15 Square Octagonal lattice, T/J=0.244, Gamma/J=0.736, generate sample projected in "
                         "computational basis at "
                      << nSweeps << " sweeps, from classical initial condition with winding number " << initialCondition
                      << "  (with 4-qubit spatially local moves)\n";
            localPIMC localpimc(Lperiodic, invTempOverJ, GammaOverJ, initialCondition, qubitsPerUpdate, qubitsPerChain,
                                seed);

            auto tstart = std::clock();
            localpimc.run(nSweeps);
            auto tend = std::clock();
            std::cout << "#Time required:" << (tend - tstart) / (double)CLOCKS_PER_SEC << " second(s)\n";
            for (int x : localpimc.firstSlice) {
                std::cout << x << " ";
            }
            std::cout << "\n";
            break;
        }
        case 3: {
            //[Convergence on square octagonal lattice cylinder by less efficient single qubit move)
            int qubitsPerChain = 4;
            int qubitsPerUpdate = 1;
            double invTempOverJ = 1 / 0.244;
            double GammaOverJ = 0.736;
            std::cout << "#24 by 15 Square Octagonal lattice, T/J=0.244, Gamma/J=0.736, generate sample projected in "
                         "computational basis at "
                      << nSweeps << " sweeps, from classical initial condition with winding number " << initialCondition
                      << " (with [more slowly converging] 1-qubit spatially local moves)\n";
            localPIMC localpimc(Lperiodic, invTempOverJ, GammaOverJ, initialCondition, qubitsPerUpdate, qubitsPerChain,
                                seed);

            auto tstart = std::clock();
            localpimc.run(nSweeps);
            auto tend = std::clock();
            std::cout << "#Time required:" << (tend - tstart) / (double)CLOCKS_PER_SEC << " second(s)\n";
            for (int x : localpimc.firstSlice) {
                std::cout << x << " ";
            }
            std::cout << "\n";
            break;
        }
        default:
            std::cerr << "Unknown experimentNo";
            return EXIT_FAILURE;
    }
}
