# SAMPL8
This repository contained files and scripts used to calculate binding free energies in the SAMPL8 Host–Guest challenge

<p align="center">
  <img width="392" alt="image" src="https://github.com/HafizSaqibAli/SAMPL9/assets/88040364/4daf7b39-8393-4470-a4c7-06b138377474">
</p>

# Abstract
Free energy drives a wide range of molecular processes such as solvation, binding, chemical reactions and conformational change. 
Given the central importance of binding, a wide range of methods exist to calculate it, whether based on scoring functions, machine-learning, 
classical or electronic structure methods, alchemy, or explicit evaluation of energy and entropy. Here we present a new energy–entropy (EE) method 
to calculate the host–guest binding free energy directly from molecular dynamics (MD) simulation. Entropy is evaluated using Multiscale Cell Correlation (MCC) 
which uses force and torque covariance and contacts at two different length scales. The method is tested on a series of seven host–guest complexes 
in the SAMPL8 (Statistical Assessment of the Modeling of Proteins and Ligands) “Drugs of Abuse” Blind Challenge. The EE-MCC binding free energies are 
found to agree with experiment with an average error of 0.9 kcal mol−1. MCC makes clear the origin of the entropy changes, showing that the large loss 
of positional, orientational, and to a lesser extent conformational entropy of each binding guest is compensated for by a gain in orientational entropy 
of water released to bulk, combined with smaller decreases in vibrational entropy of the host, guest and contacting water.
