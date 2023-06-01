function [R_Total, Rp, Rpe] = Total_Resistance(Te,Ne,Ri,Ro,opethick,l,Re,Rc,mat,Tau)
%This function will work so that I pass in the above inputs to get the
%required output i.e. The total resistance. 
%UNTITLED2 Summary of this function goes here  outputArg2 = inputArg2;
%   Detailed explanation goes here

mu_0 =  1.2566370614e-06; % This is the permitivity of free space
%Find the resistance of the plasma as a function of plasma temperature,
%electron density and the characteristic pulse time. 
Te=Te*11604.5250061657; % We multiply by this value to convert the electron temperature to Kelvin
Rp_part1 = (Ro-Ri)/((Te^(3/4))*(Ro+Ri));
Rp_part2 = sqrt((mu_0*log((1.24e07*((Te^3)/Ne)^(1/2))))/Tau);
Rp = 2.57*Rp_part1*Rp_part2;

%Rp=(8.08*0.03/(Te^(3/4)*0.01))*sqrt(mu_0*log(1.24*1e7*(Te^3/Ne)^(1/2))/Tau);

%Find the resistance of the inner and outer plate electrodes by first
%finding the area. 
AofS = (pi*(Ri*2))/4;
AofH = pi*((Ro^2)-(opethick^2)); %
%Now we find the resistance of the plate electrodes depending on material
%used where rho is the material. 
Rpe_outer = (mat*l)/AofH;
Rpe_inner = (mat*l)/AofS;
Rpe = Rpe_outer + Rpe_inner;

R_Total = Re+Rc+Rpe+Rp;
end

