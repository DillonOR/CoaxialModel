%% Program for the Qausi Resonant Capacitor Charger %%
%
% By Dr. Dillon O'Reilly
% Date: May 2022
% 
%% Brief Description of this Algorithm %%
% This file is intended to be used for the design of a quasi resonant
% flyback converter for pulsed plasma thrusters. The goal of this design is
% to implement a design reference for the parameters of the circuit such as
% switching control and feedback of the power MOSFET device. The PPU itself
% contains a charging circuit coupled to a spark plug circuit that acts as
% a switch for the discharge of the capacitor bank by which this circuit
% charges. 

% This code is split into four different modes of operation pertinent to
% the modes of operation of the QR flyback converter itself. 

%% Detailed Description of this Algorithm %% 
% MODE 1 
% In mode 1 from time t0-t1 thw power MOSFET is switched ON. 
% With the Initial current equal to the ending current of the previous charging cycle. 
% The input voltage is applied across the primary winding of the flyback
% transformer and it charges the magnetizing inductance and the leakage
% inductance of the transformer winding. It continues to do this until the
% Inpit current is equal to the current limit set by Ilim. The time
% duration of mode one is calculated via an equation. So is the primary
% side voltage and the current on the primary side which is equal to across
% all inductance values (iLm = ilk = ip = iini)

% MODE 2 
% In mode two from time t1-t2 the MOSFET is switched off and the quasi
% resonance begins between the magnetisinf inductance current, the leakage
% inductance current and the capacitor which is the sum of capacitance
% values Cr=Cin+Coss+Cw. During this mode the secondary side diode is not
% conducting since the primary side voltage is less than the transformer
% ratio i.e. Vp < Vo(k-1)/mn. During this mode we calculate the time
% duration of this mode using an equation 
clear all;
clc;


%% One Cycle Through K with Vo(0) = 0, Iini(0) = 0 %%
% Constants Declaration 
Vin = 5; Lm = 56.1e-06; Llk = 0.3e-06; npri = 8; 
nsec = 40; Cw = 0.142e-09; Cin = 200e-09; Coss = 418e-12;
Co = 2e-6; Ilim = 3.8; m = 2; Vo_0 = 0; Iini = 0; 

n = nsec/npri; a = Lm/(Lm+Llk); Cr = Cw+Coss+Cin;
Z = sqrt((Lm+Llk)/Cr); z2 = sqrt(Llk/Cr); z3 = (1/n)*(sqrt(Lm/Co));
w = 1/(sqrt((Lm+Llk)*Cr)); w2 = 1/sqrt(Llk*Cr); w3 = 1/(n*sqrt(Lm*Co));

It2 = -Ilim*cos(w)+(Vin/Z)*sin(w);
Vp = a*Vin;
Vomax = m*n*a*sqrt((Vin^2)+(Z*Ilim)^2); 
Vozero = sqrt(((Cr*Vin^2)+(Lm+Llk)*(Ilim^2))/(Cr/((m*n*a)^2)));

%% Arrays Definition 
t_ay=[]; iLm_ay = [];  t4_ay = []; Vds_ay = []; Vgate_ay = [];

Is = 0;
t=0;
Nsteps = 500 ;
k_0 = 0; % the time when the calculation starts 
k_end = 300; %the time when the calculation ends
kspan = linspace(k_0,k_end,Nsteps);

%% Perform Calculations %% 
[k,Vo] = ode45(@(k,Vo) odefun3(k,Vo,Cr,Vin,m,n,a,Lm,Llk,Ilim,Co), kspan, [1]);
[Vf,Vfind] = max(Vo);
VoZVS = find(Vo>=(m*n*a*Vin));
  
for kt = 2:1:4
    %% TIME DURATIONS %%
    T_1 = ((Lm+Llk)*(Ilim - Iini))/Vin;
    T_1 = t+T_1;
    T_2 = (1/w)*(pi-acos(Vo(kt-1)/(m*n*a*sqrt(Vin^2 + (Z*Ilim)^2)))-atan((Z*Ilim)/Vin));
    T_12 = T_1 + T_2;
    
    %% MODE 1
    t1 = t:1e-9:T_1;
    iLm1 = Iini+((Vin/(Lm+Llk))*(t1-t));
    Vds1 = 0*t1;
    Vgate1 = 5+(t1*0);
    %% MODE 2
    t2 = T_1:1e-6:T_12;
    Ip_2 = (Ilim*cos(w*(t2-T_1))+(Vin/Z)*sin(w*(t2-T_1)));
    Vp2 = Vp*cos(w*(t2-T_1))-a*Z*Ilim*sin(w*(t2-T_1));
    iLm23 = (Ilim*cos(w*(t2-T_1))+(Vin/Z)*sin(w*(t2-T_1)));
    iLm2 = linspace(iLm23(1),3.3,111);
    Vds2 = 0*t2;
    Vgate2 = 0*t2;
    iLm12 = [iLm1,iLm23];
    t12 = [t1,t2];
    Vds12 = [Vds1,Vds2];
    Vgate12 = [Vgate1,Vgate2];
    %
    T_3 = 1/w3*((pi/2)-atan(Vo(kt)/(m*n*z3*Ip_2(end))));
    T_123 = T_12+T_3;
    T_4VVS = pi/w;
    T_4ZVS = (1/w)*(pi-(acos((m*n*a*Vin)/(Vo(kt))))); %
    T_4 = T_4VVS;
    if Vo(kt) >= m*n*a*Vin
        T_4 = T_4ZVS;
        Iini = (1/(m*n*w*Lm))*sqrt((Vo(kt)^2)-(m*n*a*Vin)^2);
    end
    T_1234 =T_123+T_4;
    %% MODE 3
    t3 = T_12:1e-9:T_123;
    iLm3 = (Ip_2(end)*cos(w3*(t3-T_12)))-((Vo(kt-1)/(n*z3))*sin(w3*(t3-T_12)));
    Ip3 = (iLm3.*cos(w2*(t3-T_12))+(1/z2)*(Vin+(Vo(kt-1)/(n)))*sin(w2*(t3-T_12)));
    iss = abs(iLm3);
    Is3 = iLm3-abs(Ip3);
    Is3(Is3<=0) = 0;
    
    Vds3 = (((Vin)+(Vo(kt)/(m*n))+((Llk/Lm)*cos(w2*(t3-T_12)))+(z2*iss.*sin(w2*(t3-T_12)))));
    Vgate3 = 0*t3;
    iLm123 = [iLm12,Is3];
    t123 = [t12,t3];
    Vds123 = [Vds12,Vds3];
    Vgate123 = [Vgate12,Vgate3];
    %% MODE 4
    t4 = T_123:1e-6:T_1234;
    %Vp = (Vo(kt)/(n))*cos(w*(t4-T_123));
    Vds4 = Vin+((Vo(kt)/(m*n*a))*cos(w*(t4-T_123)));
    iLm4 = (Vo(kt)/(n*a*Z))*sin(w*(t4-T_123));
    % VVS implemented when Vo(k) < n*a*Vin
    % The MOSFET is turned on at the Valley Voltage
    VallVol = Vin - ((Vo(kt)) /(n*a));
    % if we are on valley voltage switching the initial current of the next
    % charging cycle will be zero i.e. Iini(k+1) = 0
    % When Vo(k) >= n*a*Vin, we turn on MOSFET at ZVS at Vds(t) = 0
    %
    %the ending current will be equal to
    %
    iLm = [iLm123,iLm4];
    t13 = [t123,t4];
    Vgate4 = 0*t4;
    Vdsf = [Vds123,Vds4];
    Vgatef = [Vgate123,Vgate4];
    
    t4_ay = [t4_ay,t4(end)];
    iLm_ay = [iLm_ay,iLm];
    t_ay = [t_ay,t13];
    t = t4_ay(kt-1);
    
    Vds_ay = [Vds_ay,Vdsf];
    Vgate_ay = [Vgate_ay, Vgatef];
end
tf = Vfind*T_1234;
%% Add Sub internal MODE Results together 
% t = [t1,t2,t3,t4];
% iLm = [iLm1,iLm2,iLm3,iLm4];
% Vds = [Vds1,Vds2,Vds3,Vds4];
% Vgate = [Vgate1,Vgate2,Vgate3,Vgate4];


%% PLOT Figures %%
figure(1)
plot(k, Vo);
ylabel('Capacitor Voltage');
xlabel('Charging Cycle (k)');
grid on;
legend('Capacitor Voltage','Charging Cycle (k)')

figure(2) 
yyaxis left
plot(t_ay,iLm_ay);
yyaxis right
plot(t_ay,Vgate_ay);

figure(3) 
yyaxis left
plot(t_ay,Vds_ay);
yyaxis right
plot(t_ay,Vgate_ay);


%% Functions at Bottom

function dVodk = odefun3(k,Vo,Cr,Vin,m,n,a,Lm,Llk,Ilim,Co)

dVodk = [ 
 ((Cr*((Vin^2) - ((Vo*(k))/(m*n*a))^2)) + ((Lm+Llk)*(Ilim*(k))^2)) /( Co*(((Vo^2)*(k))-((Vo^2)*(k-1))));
   ]; 
end


