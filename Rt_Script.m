%%                   DECLARE CONSTANTS               %%
% To Start our calcuation of the Total resistance, we first need to define
% our plasma temperature and our plasma density. As these figures are
% usually measured, it is best to obtain experimental findings that
% conclude what these values are for specific conditions. Although some
% methods to estimate both plasma temperature and plasma density do exist
% using the two line method
% (https://www.osapublishing.org/ao/abstract.cfm?uri=ao-59-10-3002), the
% stark broadening method for plasma temperature
% (https://www.scielo.br/j/bjp/a/zPKTQY83bg3yMVh3WcmkKsj/?lang=en), and
% using boltzmann plots
% (https://learning.edx.org/course/course-v1:EPFLx+PlasmaIntroductionX+1T_2018/home)
% they are not within the scope of this work.

% We start with a value plasma density.
% This is the plasma density obtained by a group from Taiwan who worked with Georg and 
% Chris montag. The citation: https://researchoutput.ncku.edu.tw/en/publications/plasma-behavior-in-a-solid-fed-pulsed-plasma-thruster
Ne = 8e20; 

% This is the temperature of the plasma and it is cited from a paper by
% Keidar as 2electronvolts (https://arc.aiaa.org/doi/10.2514/1.8985) but it
% was first cited by yung an in his review of thermal low power co-axial
% thrusters in https://www.researchgate.net/publication/279975945_Review_of_Thermal_Pulsed_Plasma_Thruster_-_Design_Characterization_and_Application
% in the second sentence of page 9 paragraph 1. 
Te = 2; 

% Next we define the inner and outer radii of our electrodes and we must
% also specify the thickness of our outer electrode:
Ri = 0.008;
Ro = 0.026;
thickness = 0.004; %This is an assumed thickness of 4mm for the outer electrode
opethick = Ro-thickness; %This is the radius of the inner hollow cylinder electrode

%Next we define the length of the co-axial geometry

l = 0.05;

Re = 0; %This is the resistance of the wires and leads. 
Rc = 0.03; % This is the equivalent series resistance of the capacitor, since we are using DC current and not AC, I am assuming this is 0. 

% Resistivity of various materials: %
Silver = 1.59e-08;
Copper = 1.68e-08;
Brass = 0.9e-09;
PTFE = 10e22;

mat = Copper; % Select the material you are using for your electrodes. 

% The final constant we must define is the characteristic pulse time, Tau.
% Tau is the estimated length of time it takes for a pulse to travel from
% the start of the ablation process near the teflon to the nozzle exhaust
% of the thruster head. It is determined by the square of the capacitance and inductance
% of the circuit. Where Tau = sqrt(LC). It has been shown for low power
% ablative PPT's that thrust efficiency is greatly affected by the
% characteristic pulse time wherein a pulse time of 2microseconds was found
% to be most promising as shown in yung-an chan's review:
% (https://www.researchgate.net/publication/279975945_Review_of_Thermal_Pulsed_Plasma_Thruster_-_Design_Characterization_and_Application).
% 

Tau = 4e-7;

%%                          Perform Calculation of R_Total              %%

[R_Total, Rp, Rpe] = Total_Resistance(Te,Ne,Ri,Ro,opethick,l,Re,Rc,mat,Tau);


