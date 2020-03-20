# dwave-pimc

Path Integral Monte Carlo code for analysis of a finite temperature transverse field Ising model, defined by partition function Z = Trace[exp(-H/[kB T])], and Hamiltonian

H = 1/2 sum_{i,j} J_{ij} \sigma^z_i \sigma^z_j +  \sum_i [h_i \sigma^z_i - \Gamma\sigma^x_i]

J are couplers, h are external longitudinal fields, and Gamma is the transverse field, T is the temperature, \sigma are Pauli Matrices, kB is the Boltzmann constant. J is symmetric, negative J denotes a ferromagnetic coupling. Code uses approximations exploiting the continuous time limit \Gamma < 2^16 kB T, (2^16 = number of Trotter Slices = large).  
 
   Methods exploit either single qubit Swendsen-Wang updates, or multi-qubit Swendsen-Wang updates,in the latter case specifically for regular and independent 1d ferromagnetic subsequences (chains). These match the methods explored in A King et al. https://arxiv.org/abs/1911.03446 . Code is for demonstration purposes with respect to this article and revised versions, with some functionality and optimizations omitted.

## Getting Started

Repository contains localPIMC.cpp, localPIMC.hpp, and main.cpp and README.md

Compile as:
```
g++ -O3 -std=c++11 main.cpp localPIMC.cpp -o _demo;
```
To generate projected samples following 32768 sweeps of operation from the initial condition (IC), for some lattice/algorithm combination (experimentNo):
```
./_demo experimentNo IC;
```
where 

IC == -1,0,1 for counterclockwise, ordered, and clockwise initial conditions respectively.

experimentNo=0: 12 x 9 triangular cylinder, single qubit move [SSE section results]), 32000 sweeps from initial condition.

experimentNo=1: 12 x 9 triangular cylinder, single qubit move [SSE section results]), 32000 sweeps from initial condition.

experimentNo=2: 24 x 15 square octagonal lattice, 4 qubit chain moves [Various figures, 1911.03446 and revised journal submission]).

experimentNo=3: 24 x 15 square octagonal lattice, single qubit chain moves [faster per sweep, but slower in equilibration, than 4-qubit method]).
 
### General usage

To  reproduce other square octagonal lattice published results call the localPIMC class lattice constructor (invTempOverJ, GammaOverJ, L, qubitsPerChain=4), specify the algorithm (qubitsPerUpdate=4), and use run() for appropriate number of sweeps.

   For generic models use the alternative constructor to specify the model structure.

## Authors


   Authors: Jack Raymond, Stephen Face

Copyright: D-Wave Systems

License: Apache

Last modification: March 20 2020
