# Nonlinear_DLVO_numerical

Script: Nonlinear_DLVO_numerical.m

Author: Jia-Ahn Pan

Date modified: 4/18/2023

Overview: 
This program numerically solves the full, non-linear Poisson-Bolzmann 
equations for two charged spheres in an ionic solution of equal valency (e.g. NaCl),
assuming either constant potential (CP) or constant charge (CC) at the particle surface.
It compares the computed results wit    h several approximations: (1) Derjaguin (DA) + 
linear superposition (SLA) with constant surface potential, (2) DA + linear Poisson Boltzmann (LPB) 
with constant surface potential,and (3) DA + LPB with constant surface charge. The computed
electrostatic repulsion is then combined with van der Waals attractive forces to obtain
the total interaction potential curve (DLVO theory) for both constant potential
and constant charge boundary conditions.
    The input user parameters (sphere size, surface potential, ionic 
concenctration, etc.) can be modified as desired/required. The simulation
parameters can also be modified but with caution.

Software and toolboxes required
- MATLAB
- MATLAB Partial Differential Equation Toolbox

Intermediate functions (at the end of this script)
- PoisBolzSolver_TwoSpheres_potential
- forces_two_spheres_potential
- PoisBolzSolver_TwoSpheres_charge
- forces_two_spheres_charge

Please cite the following manuscript if the software is used:
Pan, J. A.; Cho, H.; Coropceanu, I.; Wu, H.; Talapin, D. V. Stimuli-Responsive Surface Ligands for Direct Lithography of Functional Inorganic Nanomaterials. Acc. Chem. Res. 2023, 56 (17), 2286-2297. DOI: 10.1021/acs.accounts.3c00226 
