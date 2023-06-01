%% Co-Axial Pulsed Plasma Thruster %%
% Coded by: Dillon O'Reilly 
% With reference to: https://digital.wpi.edu/pdfviewer/m900nt53h,
% https://www.sciencedirect.com/science/article/abs/pii/S0094576519313931
% and Jahn. 
% This model incorporates coaxial thruster characteristic to calculate
% the performance characteristics of a coaxial pulsed plasma thruster. 


%% CLEAR ALL PARAMETERS EACH TIME CODE IS RUN %%
clear all;
clc;

%% DEFINE THE CALCULATION BOUNDARIES %%
Nsteps = 1000; %number of times the calculation repeats itself
t_0 = 0; % the time when the calculation starts 
t_end = 4e-5; %the time when the calculation ends
tspan = linspace(t_0,t_end,Nsteps); %Interval definition

% Specification of our initial values for the parameters we want to iterate
% our calculation over
Xs_0 = 0; % The position of the current sheet at t=0 
Xdots_0 = 0; % The speed of the current sheet at t=0 i.e. the derivative of the current sheet 
I_0 = 0; % The current at time t=0 
Idot_0 = 0; % The rate of change of current at time t=0

%% DECLARE CONSTANTS %%

% We first need to define our plasma temperature and our plasma density. As these figures are
% usually measured, it is best to obtain experimental findings that
% conclude what these values are for specific conditions. Although some
% methods to estimate both plasma temperature and plasma density do exist
% using the two line method
% (method), the
% stark broadening method for plasma temperature
% (https://www.scielo.br/j/bjp/a/zPKTQY83bg3yMVh3WcmkKsj/?lang=en), and
% using boltzmann plots
% (https://learning.edx.org/course/course-v1:EPFLx+PlasmaIntroductionX+1T_2018/home)
% they are not within the scope of this work. We start with a value for the plasma density.
% This is the plasma density obtained by a group from Taiwan who worked with Georg and 
% Chris montag. The citation: 8
Ne = 8e20; 

% Next, we add the temperature of the plasma in electron volts which is converted in our function to kelvin and it is cited from a paper by
% Keidar as 2electronvolts (https://arc.aiaa.org/doi/10.2514/1.8985) but it
% was first cited by yung an in his review of thermal low power co-axial
% thrusters in https://www.researchgate.net/publication/279975945_Review_of_Thermal_Pulsed_Plasma_Thruster_-_Design_Characterization_and_Application
% in the second sentence of page 9 paragraph 1. 
Te = 2; 

% Next we define the inner and outer radii of our electrodes and we must
% also specify the thickness of our outer electrode:
Ri = 0.003;
Ro = 0.012;
thickness = 0.00055; %This is an assumed thickness of 4mm for the outer electrode
opethick = Ro-thickness; %This is the radius of the inner hollow cylinder electrode

%Next we define the length of the co-axial geometry

l = 0.05;

Re = 0; %This is the resistance of the wires and leads. 
Rc = 0.03; % This is the equivalent series resistance of the capacitor. 

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
Tau = 2e-7;

% Finally we must specify the initial voltage on the capacitor banks 
V_0 = 1200; 
C = 8e-6; %capacitance

mu_0 =  1.2566370614e-06; % This is the permitivity of free space
m_0 = 2e-8;   %   Mass of current sheet at t = 0 picture at page 57    

gamma = 1.23; %ratio of specific heats which defines the impulse bit due to neutral gas expansion as described by Guman in 1968 - https://apps.dtic.mil/sti/pdfs/AD0845757.pdf page22
E_g = 5; %energy released per discharge 
M_g = 1e-8; %mass of the neutral gas
%% BRING IN OUR FUNCTIONS %% 
[LT,Lc,Le,Lce] = Total_Inductance(Ri, Ro); % My Function for calculating total inductance 

[R_Total, Rp, Rpe] = Total_Resistance(Te,Ne,Ri,Ro,opethick,l,Re,Rc,mat,Tau); % My function for calculating the total resistance

%% DEFINE OUR Co-Efficients %%

C1 = ((mu_0/(4*pi))*log(Ro/Ri))/m_0;  % eq. 2.70 page 52
C2 = -(1/C);         % eq. 2.71 page 52
C3 = ((mu_0/(4*pi))*log(Ro/Ri));       %Lce = mu_0*(2*pi)*log(Ro/Ri);
C4 = Rc + Re + Rpe +Rp; % eq. 2.71 page 52
C5= Lc+Le;
C6 = V_0;

%% PERFORM Ordinary Differential Equation CALCULATION %%
[t,x] = ode45(@(t,x) odefun3(t,x,C1,C2,C3,C4,C5,C6), tspan, [Xs_0; Xdots_0; I_0; Idot_0]);
%ode15s, ode23s

%% CALCULATE CIRCUIT PARAMETERS %%
V=V_0-x(:,2)/C;  % Voltage
E_c=0.5*C*V.^2;
E_m=0.5*(Lc+Le+x(:,1)*Lce).*x(:,4).^2;
E_ke=0.5*m_0*x(:,3).^2;
d =x(:,3).*t;

Lce_t = (mu_0/(4*pi))*log(Ro/Ri)*x(:,1);

%% CALCULATE PERFORMANCE CHARACTERISTICS OF THRUSTER %%
Imp_array = mu_0/(4*pi)*log(Ro/Ri)*x(:,2); %an array of all the impulses, to set up for trapezoidal summation

Imp_p = trapz(Imp_array); %The impulse bit of the plasma sheet

Imp_ng = (((8*(gamma-1))/((gamma^2)*gamma+1))*M_g*E_g)^(1/2);
Ibit = Imp_p + Imp_ng;

%% PLOTTING %%
figure(3)
yyaxis left
plot(t,V,'LineWidth',2);
ylabel('Voltage (V)');
xlabel('Time (s)');
yyaxis right
plot(t,x(:,4),'-r','LineWidth',2);
yyaxis right
ylabel('Current (A)');
grid on;
legend('Voltage','Current')
% ----- figures from page 54
% 
% figure(4)
% plot(t,E_c, t, E_m, t, E_ke,'LineWidth',2);
% xlabel('Time (s)');ylabel('Energy (J)');
% legend('Cpacitor Energy', 'Inductance field Energy','Kinetic Energy');
% grid on;
% 
% figure(5)
% yyaxis left
% plot(t, x(:,3),'LineWidth',2);
% ylabel('Velocity (m/s)');
% xlabel('Time (s)');
% yyaxis right
% plot(t,x(:,1),'-r','LineWidth',2);
% ylabel('Distance (m)');
% grid on;
% legend('Velocity','Distance')
% 
% figure(6)
% yyaxis left
% plot(t, x(:,1),'LineWidth',2);
% ylabel('Distance (m)');
% xlabel('Time (s)');
% yyaxis right
% plot(t,Lce_t,'-r','LineWidth',2);
% ylabel('Inductance (nH)');
% grid on;
% legend('Inductance','Distance')

% Save file data to Excel
%filename = 'PETRUSModel.xlsx';
%writecell(V,filename,'Sheet',1,'V',false);









