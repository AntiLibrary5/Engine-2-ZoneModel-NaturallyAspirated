function [dQu_dt, dQb_dt] = ComputeHeatLosses(Tu, Tb, P, RunCond, EngGeom, theta, V, Vb)
%UNTITLED2 Summary of this function goes here
%   Tu [K] : Temperature of the Unburnt gases
%   Tb [K] : Temperature of the Burnt gases
%   P [Pa] : In-cylinder Pressure
%   RunCond : Running Conditions (struct)
%   EngGeom : Engine Geometry (struct)
%   V [m3] : Instantaneous chamber volume
%   Vb [m3] : Instantaneous volume occupied by the burnt gases
global CA
CA_SoCb = RunCond.SparkAngle;
T_bdc = RunCond.T_plenum;
P_bdc = RunCond.P_plenum;
V_bdc = Vcyl(-180,EngGeom,RunCond);

MPS = 2*EngGeom.S*(RunCond.N/60); % [m/s] Mean Piston Speed
hChamber = V/((pi*(EngGeom.B/2)^2));% [m] Instantaneous Height of the chamber
radiusBurntGases = sqrt(Vb/pi/hChamber); % [m] Radius of the virtual cylinder containing the burnt gases
A_piston = pi*(EngGeom.B/2)^2; % [m2] Area of the piston

A_piston_B = pi*radiusBurntGases^2; % [m2] Area Burnt gases/piston
A_head_B = A_piston_B; % [m2] Area Burnt gases/cylinder head

A_piston_U = A_piston-A_head_B; % [m2] Area Unburnt gases/piston
A_head_U = A_piston_U; % [m2] Area Unburnt gases/cylinder head

A_liner = hChamber*pi*EngGeom.B; % [m2] Area of the liner
% A_BurntUnburnt = hChamber*2*pi*radiusBurntGases; % [m2] Area of the interface burnt gases / unburnt gases
A_BurntUnburnt =0;

T_head = 400;% [K] Temperature of the cylinder head
T_liner = 475;% [K] Temperature of the liner
T_piston = 600;% [K] Temperature of the piston

if ~exist('P_motored','var')
    P_ref = P;
else
    i = find(CA==round(theta));
    P_ref = P_motored(i);
end

%% C2
if theta >= CA_SoCb && theta <= 180
    C2 = 3.24e-3;
else
    C2 = 0;
end

%% C1
if theta >= -180 && theta <= 180
    C1 = 2.28; % for compression and expansion
else
    C1 = 6.18;
end

alphaU = 110*(EngGeom.B^(-0.2)*(P/1e5)^0.8*( C1*MPS + C2*V*T_bdc/P_bdc/V_bdc*(P-P_ref) )^0.8*Tu^(-0.53)); % [W/(m².K)] Woschni heat transfer coef. for the Unburnt gases
alphaB = 110*(EngGeom.B^(-0.2)*(P/1e5)^0.8*( C1*MPS + C2*V*T_bdc/P_bdc/V_bdc*(P-P_ref) )^0.8*Tb^(-0.53)); % [W/(m².K)] Woschni heat transfer coef. for the Burnt gases

dQu_dt = (alphaU*A_piston_U*(T_piston-Tu)...
    +alphaU*A_head_U*(T_head-Tu)...
    +alphaU*A_liner*(T_liner-Tu)...
    -alphaB*A_BurntUnburnt*(Tu-Tb)); % [J/s] Heat exchanged by the Unburnt gases
dQb_dt = (alphaB*A_piston_B*(T_piston-Tb)...
    +alphaB*A_head_B*(T_head-Tb)...
    +alphaB*A_BurntUnburnt*(Tu-Tb)); % [J/s] Heat exchanged by the Burnt gases


end

