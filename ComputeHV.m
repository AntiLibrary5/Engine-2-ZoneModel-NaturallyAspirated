function [FuelProperties] = ComputeHV(Yreac,Yprod)
%Compute Lower and Higher heating values.
%   Yreac: Mass fractions of [fuel O2 N2 CO2 H2O] in the reactants
%   Yprod: Mass fractions of [fuel O2 N2 CO2 H2O] in the products

T=298.15; % [K] Temperature

%% Computation of LHV
hi=[Compute_hi('Gasoline (High Octane)',T) Compute_hi('O2',T) Compute_hi('N2',T) Compute_hi('CO2',T) Compute_hi('H2O',T)]; % [J/kg] Specific enthalpy of the species [fuel O2 N2 CO2 H2O] @ T1
h0_reac = sum(hi.*Yreac); % [J/kg] Specific enthalpy of the unburnt gases @ T1
h0_prod = sum(hi.*Yprod); % [J/kg] Specific enthalpy of the burnt gases @ T1
FuelProperties.LHV = abs( 1/Yreac(1)*(h0_prod-h0_reac) ); % [J/kg] Lower Heating Value of the fuel

%% Computation of HHV
FuelProperties.HHV = FuelProperties.LHV+2.443e6*(Yprod(5)/Yreac(1)); % [J/kg] Higher Heating Value of the fuel

FuelProperties.AFRst = (Yreac(2)+Yreac(3))/Yreac(1);
end

