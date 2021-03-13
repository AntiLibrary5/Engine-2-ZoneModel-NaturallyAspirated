function [FP] = ComputeFuelProperties(FuelType, phi)
%Compute Lower and Higher heating values.
%   Yreac: Mass fractions of [fuel O2 N2 CO2 H2O] in the reactants
%   Yprod: Mass fractions of [fuel O2 N2 CO2 H2O] in the products
% Output :
%   FP : Struct containing the fuel properties

switch FuelType
    case 'Gasoline (High Octane)'
        FP.Type = 'Gasoline (High Octane)';
        FP.x = 7.77; % carbon fraction of fuel (Gasoline/Diesel/LPG)
        FP.y = 13.08; % hydrogen fraction of fuel (Gasoline/Diesel/LPG)
        FP.OctaneNumber = 93;
%         FuelProperties.phi=1; % fuel-air equiv. ratio (Gasoline/Diesel/LPG)
    case 'Diesel'
        FP.Type = 'Diesel';
        FP.x = 10.82; % carbon fraction of fuel (Gasoline/Diesel/LPG)
        FP.y = 18.73; % hydrogen fraction of fuel (Gasoline/Diesel/LPG)
%         FuelProperties.phi=0.8; % fuel-air equiv. ratio (Gasoline/Diesel/LPG)
    case 'LPG'
        FP.Type = 'LPG';
        FP.x = 3.53; % carbon fraction of fuel (Gasoline/Diesel/LPG)
        FP.y = 8.98; % hydrogen fraction of fuel (Gasoline/Diesel/LPG)
%         FuelProperties.phi=1; % fuel-air equiv. ratio (Gasoline/Diesel/LPG)
    otherwise
        disp('Unknown fuel')
end

%% Computation of Yreac & Yprod
RunCond.Fuel.x = FP.x;   %No. Carbon atoms
RunCond.Fuel.y = FP.y;   %No. of Hydrogen atoms
RunCond.phi = phi;    %Fuel to Air equivalence ratio
[~,~,Yreac,Yprod] = ComputeCombustionComposition(0, RunCond);

%% Computation of LHV
Tfuel = 298.15; % [K] Fuel Temperature
hi=[Compute_hi('Gasoline (High Octane)',Tfuel) Compute_hi('O2',Tfuel) Compute_hi('N2',Tfuel) Compute_hi('CO2',Tfuel) Compute_hi('H2O',Tfuel)]; % [J/kg] Specific enthalpy of the species [fuel O2 N2 CO2 H2O] @ T1
h0_reac = sum(hi.*Yreac); % [J/kg] Specific enthalpy of the unburnt gases @ T1
h0_prod = sum(hi.*Yprod); % [J/kg] Specific enthalpy of the burnt gases @ T1
FP.LHV = abs( 1/Yreac(1)*(h0_prod-h0_reac) ); % [J/kg] Lower Heating Value of the fuel

%% Computation of HHV
FP.HHV = FP.LHV+2.443e6*(Yprod(5)/Yreac(1)); % [J/kg] Higher Heating Value of the fuel

%% Computation of AFRst
FP.AFRst = (Yreac(2)+Yreac(3))/Yreac(1); % [-] Stoechiometric Air-to-Fuel Ratio
end