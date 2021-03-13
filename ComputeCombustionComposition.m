function [Mreac_ave,Mprod_ave,Yreac,Yprod,Mi] = ComputeCombustionComposition(ResMassFrac, RunCond)
%UNTITLED Summary of this function goes here
% Input arguments:
% - Fuel.x: number of C atoms in the fuel molecule
% - Fuel.y: number of H atoms in the fuel molecule
% - Fuel.phi: Fuel to Air equivalent ratio
% - RunCond.ResFrac = mResiduals/TotalTrappedMass
% Output arguments:
% - Mreac_ave [kg/mol] Average molar mass of the fresh gases

% ResMassFrac = RunCond.ResMassFrac;
x = RunCond.Fuel.x;   %No. Carbon atoms
y = RunCond.Fuel.y;   %No. of Hydrogen atoms
phi = RunCond.phi;    %Fuel to Air equivalence ratio

C=12;   % [g/mol] Molar mass of Carbon 
O=16;   % [g/mol] Molar mass of Oxygen 
H=1;    % [g/mol] Molar mass of Hydrogen
N=14;   % [g/mol] Molar mass of Nitrogen

%% Molar fractions
Xreac=[1 (x+(y/4))*(1/phi) (x+(y/4))*(1/phi)*(0.79/0.21) 0 0];
Xreac=Xreac./sum(Xreac); % Molar fraction in the fresh gases [fuel O2 N2 CO2 H2O]

if phi > 1
	Xprod=[phi-1, 0, (x+(y/4))*0.79/0.21, x, y/2];
else
	Xprod=[0, (x+(y/4))*((1/phi)-1), (x+(y/4))*(0.79/(0.21*phi)), x, 0.5*y];
end

Xprod=Xprod./sum(Xprod); % Molar fraction in the burnt gases [fuel O2 N2 CO2 H2O]

%% Mean molecular mass
Mi=[x*C+y*H 2*O 2*N C+(2*O) (2*H)+O].*1e-3; % [kg/mol] Molar masses of the species [fuel O2 N2 CO2 H2O]


Mreac_ave=sum(Xreac*Mi');%[kg/mol] Average molar mass of the fresh gases
Mprod_ave=sum(Xprod*Mi');%[kg/mol] Average molar mass of the burnt gases

%% Mass Fractions

Yreac = (Xreac.*Mi)./Mreac_ave;% Mass fraction in the fresh gases [fuel O2 N2 CO2 H2O]
Yprod = (Xprod.*Mi)./Mprod_ave;% Mass fraction in the burnt gases [fuel O2 N2 CO2 H2O]
Yreac = (1-ResMassFrac)*Yreac+ResMassFrac*Yprod;


% CxHy + (x+y/4)*(1/phi)*(O2+4*N2) -> z CO2 + N2 + (phi-1)*CxHy
% phi*CxHy + (x+y/4)*(O2+4*N2) -> (phi-1) CO2 + N2 + (phi-1)*CxHy