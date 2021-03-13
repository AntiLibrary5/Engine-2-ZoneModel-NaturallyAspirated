function [V, dVdt] = Vcyl(CA,EngGeom,RunCond)
% Computes instantaneous volume of combustion chamber from a given Crank Angle value.
% CA [°] is the Crank Angle. CA=0 is TDC.
% EngGeom is a struct containing Engine Geometry data.

S = EngGeom.S;      % Stroke [m]
B = EngGeom.B;      % Bore [m]
R = EngGeom.R;      % Drive shaft length [m]
L = EngGeom.L;      % Connecting rod length [m]
eta = EngGeom.eta;  % Compression ratio
% Derived data
alpha = L/R;
A = pi/4*B^2;       % Piston Surface [m^2]
Vs = S*A;           % Displaced Volume [m^3]
Vc = 1/(eta-1)*Vs;  % Compression volume [m^3]

CA = (2*pi/360).*CA;% CA has to be converted from [°] to [rad]

V = Vc + (Vs/2).*(alpha+1-cos(CA)-sqrt(alpha.^2-(sin(CA)).^2)); % [m3]
dVdt = (Vs/2).*sin(CA).*(1+cos(CA)./sqrt(alpha^2 - sin(CA).^2))*pi*RunCond.N/30; % [m3/s]
end