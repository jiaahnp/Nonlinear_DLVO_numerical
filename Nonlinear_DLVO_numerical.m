%{ 
Script: Nonlinear_DLVO_numerical.m
Author: Jia-Ahn Pan
Date modified: 11/29/2022

Overview: 
This program numerically solves the full, non-linear Poisson-Bolzmann 
equations for two charged spheres in an ionic solution of equal valency (e.g. NaCl),
assuming constant potential or constant charge at the particle surface.
It compares the computed results with several approximations: (1) Derjaguin (DA) + 
linear superposition (SLA), (2) DA + linear Poisson Boltzmann (LPB) with constant surface potential,
and (3) DA + LPB with constant surface charge. The computed electrostatic
repulsion is then combined with van der Waals attractive forces to obtain
the total interaction potential curve (DLVO theory) for both constant potential
and constant charge boundary conditions.
    The input user parameters (sphere size, surface potential, ionic 
concenctration, etc.) can be modified as desired/required. The simulation
parameters can also be modified but with caution.

Software and toolboxes required
- MATLAB
- MATLAB Phased Array System Toolbox
- MATLAB Partial Differential Equation Toolbox

Intermediate functions (at the end of this script)
- PoisBolzSolver_TwoSpheres_potential
- forces_two_spheres_potential
- PoisBolzSolver_TwoSpheres_charge
- forces_two_spheres_charge

%}


%% User parameters (modify freely based on system of study)
r_nm = 1.5; % radius of sphere in nm
psi_0_mV = 50; % surface potential in mV
c_0_molar = 0.001; % ionic concentration in molarity (M), e.g. c_0_molar = 1 for 1M NaCl
z = 1; %z = ionic charge of cation = ionic charge of anion (must be equal)
epsilon_r = 36.7; % relative permittivity (dielectric constant) of solvent
A = 6.22e-20; % Hamaker constant of the colloidal sphere, e.g., 6.22e-20 for CdSe
T = 298; % temperature in Kelvin

%% Constants and intermediate variables
%constants
e = 1.60217662e-19;
epsilon_0 = 8.854e-12;
N_A = 6.02214076e23;

%intermediate variables
kT = physconst('Boltzmann')*T;
c_0 = c_0_molar*N_A*1e3; %conc of ion pairs in 1/m3
r = r_nm*1e-9; %radius in meters

kappa = ((2*z^2*e^2*c_0)/(epsilon_r*epsilon_0*kT))^(0.5); % 1/debye length
r_red = r*kappa; %reduced sphere radius
psi_0_red = psi_0_mV*e/1e3/kT; %reduced surface potential

%% Simulation parameters (modify with caution)
h_closest_nm = 0.01/kappa*1e9; %closest NC surface-to-surface separation in nanometers
h_farthest_nm = 5/kappa*1e9; %farthest NC surface-to-surface separation in nanometers
%IMPORTANT: to get an accurate value of the potential (u_array_kT),
%h_farthest_nm needs to be big enough such that the force (f) is close
%to zero at this point. This can be checked in Figure 1.
h_points_nm = 50; % number of points to sample between closest and farthest surface separation
% can increase h_points_nm if want smaller data spacings
mmax = 0.05; %target maximum mesh edge length, e.g. 0.03 (smaller better but longer time)
% Note: Discontinuities/uneveness in the force curve can be observed if
% the mmax value used is too large. Reduce mmax if this is observed

h_array_nm = logspace(log10(h_closest_nm),log10(h_farthest_nm),h_points_nm)'; % array of surface distances (nm)
h_array_m = h_array_nm*1e-9;
h_array = h_array_nm.*1e-9.*kappa; %reduced surface distances

%% Estimate surface charge density 
% The surface charge density is estimated from the surface potential by 
% solving PB equation for an isolated sphere using Gauss law (differential form)

h_red = r_red*10; %use far enough distance to assume isolation
disp("Calculating surface charge density from potential...")
[f2,u,eta_thetazero,u_thetazero,charge_0_red] = forces_two_spheres_potential(h_red,r_red,psi_0_red,0.01);
charge_C_m2 = charge_0_red/e*epsilon_0*epsilon_r*kappa*kT;
charge_e_nm2 = charge_C_m2/e/(1e9)^2;

%% Electrostatics - Numerical calculation (exact)
[u_array_potential, f_array_potential] = PoisBolzSolver_TwoSpheres_potential(r_red, psi_0_red, h_array, mmax);
[u_array_charge, f_array_charge] = PoisBolzSolver_TwoSpheres_charge(r_red, psi_0_red, charge_0_red, h_array, mmax);

%Change numerical electrostatic solution to correct units
a = c_0*kT*4*pi*r*(psi_0_mV/1e3*e/kT)^2/kappa^2/kT; % new scale to units of kT
u_array_potential_kT = u_array_potential*a;
u_array_charge_kT = u_array_charge*a;

%% Plot calculated force and energy curves for numerical results
figure
hold on
plot(h_array, f_array_potential,'-diamond','Color',[0.6275,0.6431,0.8706])
plot(h_array, u_array_potential, '-diamond','Color',[0.27,0.3216,1])
plot(h_array, f_array_charge,'-square','Color',[0.9098, 0.6549, 0.6549])
plot(h_array, u_array_charge, '-square','Color',[0.9098,0.2275, 0.2275])
hold off

ttle = title({'Electrostatic force and energy between two colloidal spheres',...
    strcat('(Reduced radius = \kappar = ',num2str(r_red),'; Reduced potential = e\psi_0/kT = ',num2str(psi_0_red),...
    ')')});
ttle.FontSize = 12;
labelsize = 12;
xlabel('Surface separation (\kappa^{-1}h)', 'FontSize', labelsize)
ylabel('Interaction force or energy (reduced units)', 'FontSize', labelsize)
lgd = legend('Constant Potential (Force)','Constant Potential (Energy)','Constant Charge (Force)','Constant Charge (Energy)');
lgd.FontSize = 12;

%% Electrostatics - Analytical approximations
h_array_ana_m = (0.09e-9:0.01e-9:40e-9);
gamma = tanh(z*e*psi_0_mV/1e3/(4*kT)); 

%electrostatic double layer potential Derjaguin + linear superpos. (LSA) - constant potential
w_e = 64*pi*kT*r*c_0*gamma^2/kappa^2*exp(-kappa*h_array_ana_m);
w_e_kT = w_e./kT; %get in terms of kT

%electrostatic double layer potential HHF approx (Derjaguin + linearized PB (LPB) - constant potential) 
w_e_HHF_potential = 2*pi*epsilon_r*epsilon_0*r*(psi_0_mV/1e3)^2*(log(1+exp(-kappa*h_array_ana_m)));
w_e_HHF_potential_kT = w_e_HHF_potential/kT;

%electrostatic double layer potential HHF approx (Derjaguin + linearized PB (LPB) - constant charge) 
w_e_HHF_charge = -2*pi*epsilon_r*epsilon_0*r*(psi_0_mV/1e3)^2*(log(1-exp(-kappa*h_array_ana_m)));
w_e_HHF_charge_kT = w_e_HHF_charge/kT;


%% Van der Waals interactions (exact analytical solution)
w_v = -A./6.*(2.*r.^2./h_array_m./(4.*r+h_array_m) + 2.*r.^2./(2.*r+h_array_m).^2 + log(h_array_m.*(4.*r+h_array_m)./(2.*r+h_array_m).^2));
w_v_kT = w_v/kT; %get in terms of kT 
%% Plot electrostatic repulsion: calculated vs approximations 
figure
hold on
plot(h_array_m*1e9,u_array_potential_kT, '-diamond','Color',[0.27,0.3216,1])%, 'MarkerSize', 3)
plot(h_array_m*1e9,u_array_charge_kT, '-square','Color',[0.9098,0.2275, 0.2275])%, 'MarkerSize', 3)

plot(h_array_ana_m*1e9,w_e_kT)
plot(h_array_ana_m*1e9,w_e_HHF_potential_kT)
plot(h_array_ana_m*1e9,w_e_HHF_charge_kT)

ttle = title({'Interaction energy between two colloidal spheres',...
    strcat('(r = ',num2str(r*1e9),' nm, \psi_0 = ',num2str(psi_0_mV),...
    ' mV, \kappa^{-1} =',num2str(1e9/kappa),' nm)')});
ttle.FontSize = 12;
labelsize = 12;
xlabel('Surface separation, h (nm)', 'FontSize', labelsize)
ylabel('Interaction energy (kT)', 'FontSize', labelsize)
%xlim([0 20])
%ylim([0 7])
lgd = legend('Numerical (constant potential)','Numerical (constant charge)',...
   'Approx. - DA + LSA','Approx. - DA + LPB (constant potential)','Approx. - DA + LPB (constant charge)');
lgd.FontSize = 12;
hold off

%% Plot total DLVO (constant potential boundary)
u_dlvo_potential_kT = u_array_potential_kT+w_v_kT; %for numerical with vdw

figure
plot(h_array_m*1e9, w_v_kT,'-', h_array_m*1e9,u_array_potential_kT,'-s', h_array_m*1e9, u_dlvo_potential_kT,'-o')
ttle = title({'DLVO interaction energy (constant potential)',...
    strcat('r = ',num2str(r*1e9),' nm, \psi_0 = ',num2str(psi_0_mV),...
    ' mV, c_0 =',num2str(c_0_molar), ' M, T=', num2str(T), ' K, \kappa^{-1} =',num2str(1e9/kappa),' nm')});
ttle.FontSize = 12;
labelsize = 12;
xlabel('Surface separation, h (nm)', 'FontSize', labelsize)
ylabel('Interaction energy (kT)', 'FontSize', labelsize)
% xlim([0 20])
% ylim([-6 6])
lgd = legend('Van der Waals','Electrostatics (calculated)',...
    'Total');
lgd.FontSize = 12;
hold off

%% Plot total DLVO (constant charge boundary)
u_dlvo_charge_kT = u_array_charge_kT+w_v_kT; %for numerical with vdw

figure
plot(h_array_m*1e9, w_v_kT,'-', h_array_m*1e9,u_array_charge_kT,'-s', h_array_m*1e9, u_dlvo_charge_kT,'-o')
ttle = title({'DLVO interaction energy (constant charge)',...
    strcat('r = ',num2str(r*1e9),' nm, \sigma_0 = ', num2str(charge_C_m2), ...
    ' C/m^2, \psi^{isolated}_0\approx',num2str(psi_0_mV),...
    ' mV, c_0 =',num2str(c_0_molar), ' M, T=', num2str(T), ' K, \kappa^{-1} =',num2str(1e9/kappa),' nm')});
ttle.FontSize = 12;
labelsize = 12;
xlabel('Surface separation, h (nm)', 'FontSize', labelsize)
ylabel('Interaction energy (kT)', 'FontSize', labelsize)
% xlim([0 20])
% ylim([-6 6])
lgd = legend('Van der Waals','Electrostatics (calculated)',...
    'Total');
lgd.FontSize = 12;
hold off

%% Intermediate functions
function [u_array, f_array] = PoisBolzSolver_TwoSpheres_potential(r_red, psi_0_red, h_array, mmax)
%{
Function: PoisBolzSolver_TwoSpheres_potential

Evaluate forces at various separations and calculate potential, U

Overview: This function numerically solves the full, non-linear Poisson-Bolzmann 
equations for two charged spheres in an ionic solution of equal valency (e.g. NaCl),
assuming constant potential at the particle surface. It gives the (reduced) potential 
and (reduced) force curves as a function of the (reduced)
surface-to-surface distance between the spheres.

Input
    - r_red (reduced radius)= (radius(m))*kappa(Debye length) 
    - psi_0_red (reduced potential): psi_0_red = (Surface potential (V))*e/kT 
    - h_array (array of surface distance to be evaluated) Note: to get an
    accurate value for the potential (u_array), h_array needs to be
    evaluated to a distance far enough that the force, f is close to zero
    - mmax (mesh density)
Output
    - u_array (reduced potential between the spheres as a function of h_array)
    - f_array (reduced force between the spheres as a function of h_array)

%}

f_array = zeros(size(h_array)); %array of forces with same size as h_array

for i = 1:length(h_array) %for
    f_array(i) = forces_two_spheres_potential(h_array(i),r_red,psi_0_red,mmax);
    disp("Calculating (constant potential):" +i + " of " + length(h_array) + " steps done.")
end

% Calculate potential by integrating force from furthest to that point

h_array_flip = flipud(h_array);
f_array_flip = flipud(f_array);
u_array = cumtrapz(h_array_flip, -f_array_flip)./2; %integrate negative force and divide two (between spheres)
u_array = flipud(u_array); %flip U back to proper 

end

function [f2,u,eta_thetazero,u_thetazero,charge_reduced] = forces_two_spheres_potential(h,r,psi_0,mmax)
    
%intermediate parameters
    theta_0 = pi;
    cent_to_cent = h+2*r; %center to center distance
    u0 = psi_0; %initial guess

    % constants for cartesian to bipolar transformation
    eta_0 = acosh(1+(h/(2*r)));
    c_0 = (sqrt(cent_to_cent^2-4*r^2))/2; % or c_0 = r*sinh(eta_0)

    %parameters for pdenonlin solver: -(del)(c(del)u)+au = f
    a = 0;
    c = '-sin(y)./(cosh(x)-cos(y))'; %x is eta, y is theta
    f = strcat(num2str(c_0^2),'.*sin(y).*sinh(u)./(cosh(x)-cos(y)).^3'); 

    % Create model, modify and solve
    %create model instance
    model = createpde();

    %create square object
    R1 = [3,4,0,eta_0,eta_0,0,0,0,theta_0,theta_0]';
    gm = [R1];
    sf = 'R1';
    ns = char('R1');
    ns = ns';
    g = decsg(gm,sf,ns);

    %insert object into model
    geometryFromEdges(model,g);

    % figure
    % pdegplot(model,'EdgeLabels','on');
    % axis equal

    %boundary conditions
    applyBoundaryCondition(model,'dirichlet','Edge',2,'u',psi_0); % constant potential 2
    applyBoundaryCondition(model,'neumann','Edge',[1,3,4],'q',0, 'g',0); % rest has zero e field

    %Generate mesh -- need linear mesh for pdenonlin?
    generateMesh(model,'GeometricOrder','linear','Hmax',mmax);

    % figure
    % pdemesh(model);     
    % axis equal

    %run the numerical solver pdenonlin to get the potential u
    u = pdenonlin(model,c,a,f,'U0',u0,'Jacobian','full'); 
    
    %extract the coordinates for the potential
    eta = model.Mesh.Nodes(1,:)';
    theta = model.Mesh.Nodes(2,:)';

%     figure
%     pdeplot(model,'XYData',u,'ZData',u)
%     title('Numerical Solution');
%     xlabel('eta')
%     ylabel('theta')

    % Extract potential at median plane (eta=0) and get force between spheres

    % Extract potential at median plane (eta=0) and arrange
    pos_etazero = find(eta==0); %positions in vector for eta=0
    theta_etazero = theta(pos_etazero); %extract positions which theta = 0
    u_etazero = u(pos_etazero);
    theta_u_etazero = [theta_etazero'; u_etazero']'; %combine u and theta
    theta_u_etazero = sortrows(theta_u_etazero, 1); %sort theta and u in increasing theta
    theta_etazero = theta_u_etazero(:,1); % split up again
    u_etazero = theta_u_etazero(:,2);

    %first derivative of u at midplane
    u_etazero_der = gradient(u_etazero, theta_etazero); %take first derivative
    theta_etazero = theta_etazero(2:end); %remove first value to remove asymptotes
    u_etazero = u_etazero(2:end);
    u_etazero_der = u_etazero_der(2:end);

    % figure
    % scatter(theta_etazero, u_etazero,5, 'filled')

    % integrand for force
    f_integrand = c_0^2.*sin(theta_etazero)./((1-cos(theta_etazero)).^2)...
        .*2.*(cosh(u_etazero)-1) + sin(theta_etazero).*(u_etazero_der).^2;

    % figure
    % plot(theta_etazero,u_etazero, theta_etazero, f_integrand)

    % integrate to get reduced force , f*=-kappa^2/(nkT) * dU/dR where
    %   U=potential energy
    f = 2*pi.*trapz(theta_etazero,f_integrand); %integrate (over 0 to pi) wrt theta to get reduced force f*
    f2 = f/(2*pi*r*(psi_0)^2); %double reduced force
    

    % Calculate charge from potential gradient with Gauss law

    % Extract u as function of eta (line) with theta=0
    pos_thetazero = find(theta==0); %positions in vector for theta=0
    eta_thetazero = eta(pos_thetazero);
    u_thetazero = u(pos_thetazero);
    eta_u_thetazero = [eta_thetazero'; u_thetazero']'; %combine u and theta
    eta_u_thetazero = sortrows(eta_u_thetazero, 1); %sort theta and u in increasing theta
    eta_thetazero = eta_u_thetazero(:,1); % split up again
    u_thetazero = eta_u_thetazero(:,2);
    
    u_thetazero_der = gradient(u_thetazero, eta_thetazero); % differentiate 

    charge_reduced = (cosh(eta_0)-cos(0)) / (r*sinh(eta_0)) * u_thetazero_der(end); %calculate reduced charge

end

function [u_array, f_array] = PoisBolzSolver_TwoSpheres_charge(r_red, psi_0_red, charge_0_red, h_array, mmax)
%{
Function: PoisBolzSolver_TwoSpheres_charge

Evaluate forces at various separations and calculate potential, U

Overview: This function numerically solves the full, non-linear Poisson-Bolzmann 
equations for two charged spheres in an ionic solution of equal valency (e.g. NaCl),
assuming constant charhe at the particle surface. It gives the (reduced) potential 
and (reduced) force curves as a function of the (reduced)
surface-to-surface distance between the spheres.

Input
    - r_red (reduced radius)= (radius(m))*kappa(Debye length) 
    - psi_0_red (reduced potential): psi_0_red = (Surface potential (V))*e/kT 
    - charge_0_red (reduced charge)
    - h_array (array of surface distance to be evaluated) Note: to get an
    accurate value for the potential (u_array), h_array needs to be
    evaluated to a distance far enough that the force, f is close to zero
    - mmax (mesh density)
Output
    - u_array (reduced potential between the spheres as a function of h_array)
    - f_array (reduced force between the spheres as a function of h_array)

%}


f_array = zeros(size(h_array)); %array of forces with same size as h_array

for i = 1:length(h_array) %for
    f_array(i) = forces_two_spheres_charge(h_array(i),r_red,psi_0_red, charge_0_red, mmax);
    disp("Calculating (constant charge):" + i + " of " + length(h_array) + " steps done.")
end

% Calculate potential by integrating force from furthest to that point

h_array_flip = flipud(h_array);
f_array_flip = flipud(f_array);
u_array = cumtrapz(h_array_flip, -f_array_flip)./2; %integrate negative force and divide two (between spheres)
u_array = flipud(u_array); %flip U back to proper 

end

function [f2] = forces_two_spheres_charge(h,r,psi_0,charge_0_red, mmax)
    %intermediate parameters
    theta_0 = pi;
    cent_to_cent = h+2*r; %center to center distance
    u0 = psi_0; %initial guess

    % constants for cartesian to bipolar transformation
    eta_0 = acosh(1+(h/(2*r)));
    c_0 = (sqrt(cent_to_cent^2-4*r^2))/2; % or c_0 = r*sinh(eta_0)

    %parameters for pdenonlin solver: -(del)(c(del)u)+au = f
    a = 0;
    c = '-sin(y)./(cosh(x)-cos(y))'; %x is eta, y is theta
    f = strcat(num2str(c_0^2),'.*sin(y).*sinh(u)./(cosh(x)-cos(y)).^3'); 

    % Create model, modify and solve
    %create model instance
    model = createpde();

    %create square object
    R1 = [3,4,0,eta_0,eta_0,0,0,0,theta_0,theta_0]';
    gm = [R1];
    sf = 'R1';
    ns = char('R1');
    ns = ns';
    g = decsg(gm,sf,ns);

    %insert object into model
    geometryFromEdges(model,g);

    % figure
    % pdegplot(model,'EdgeLabels','on');
    % axis equal

    function bc = bcfunc(location, state)
        bc = charge_0_red*r*sinh(eta_0)./(cosh(eta_0) - cos(location.y));
    end
    applyBoundaryCondition(model,'neumann','Edge',2,'q',0, 'g',@bcfunc); % rest has zero e field
    applyBoundaryCondition(model,'neumann','Edge',[1,3,4],'q',0, 'g',0); % rest has zero e field

    %Generate mesh -- need linear mesh for pdenonlin?
    generateMesh(model,'GeometricOrder','linear','Hmax',mmax);

    % figure
    % pdemesh(model);     
    % axis equal

    %run the numerical solver pdenonlin to get the potential u
    u = pdenonlin(model,c,a,f,'U0',u0,'Jacobian','full'); 

    %extract the coordinates for the potential
    eta = model.Mesh.Nodes(1,:)';
    theta = model.Mesh.Nodes(2,:)';

    % figure
    % pdeplot(model,'XYData',u,'ZData',u)
    % title('Numerical Solution');
    % xlabel('eta')
    % ylabel('theta')

    % Extract potential at median plane (eta=0) and get force between spheres

    % Extract potential at median plane (eta=0) and arrange
    pos_etazero = find(eta==0); %positions in vector for eta=0
    theta_etazero = theta(pos_etazero); %extract positions which theta = 0
    u_etazero = u(pos_etazero);
    theta_u_etazero = [theta_etazero'; u_etazero']'; %combine u and theta
    theta_u_etazero = sortrows(theta_u_etazero, 1); %sort theta and u in increasing theta
    theta_etazero = theta_u_etazero(:,1); % split up again
    u_etazero = theta_u_etazero(:,2);

    %first derivative of u at midplane
    u_etazero_der = gradient(u_etazero, theta_etazero); %take first derivative
    theta_etazero = theta_etazero(2:end); %remove first value to remove asymptotes
    u_etazero = u_etazero(2:end);
    u_etazero_der = u_etazero_der(2:end);

    % figure
    % scatter(theta_etazero, u_etazero,5, 'filled')

    % integrand for force
    f_integrand = c_0^2.*sin(theta_etazero)./((1-cos(theta_etazero)).^2)...
        .*2.*(cosh(u_etazero)-1) + sin(theta_etazero).*(u_etazero_der).^2;

    % figure
    % plot(theta_etazero,u_etazero, theta_etazero, f_integrand)

    % integrate to get reduced force , f*=-kappa^2/(nkT) * dU/dR where
    %   U=potential energy
    f = 2*pi.*trapz(theta_etazero,f_integrand); %integrate (over 0 to pi) wrt theta to get reduced force f*
    f2 = f/(2*pi*r*(psi_0)^2); %double reduced force  
end

