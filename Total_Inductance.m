function [LT,Lc,Le,Lce] = Total_Inductance(Ri, Ro)
%UNTITLED3 Summary of this function goes here
%   Total Inductance = Lc + Le + Lce

mu_0 =  1.2566370614e-06; % This is the permitivity of free space

Lce = mu_0*(2*pi)*log(Ro/Ri);
Lc = 3.4e-8; 
Le = 0;

LT = Lc+Le+Lce;

end
