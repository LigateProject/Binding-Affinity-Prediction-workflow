Order of scripts:
0) checkProtein.sh
1) GROMACS_LiGen_integration.sh
2) solvate.sh
3) addIons.sh (only adding as many ions as are needed to neutralise the box)
3b) energyMinimisation.sh (minimal to no changes to protein and ligand structure thanks to position restraints on heavy atoms, ice-like water)
3c) equilibrationNpT.sh (equilibration with position restraints on the protein and on the ligand for unstable systems)
4a) prepareSimulations.sh
4b) runSimulations.sh
5) getFreeEnergies.sh
