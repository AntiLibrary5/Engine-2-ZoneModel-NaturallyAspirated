function [hi] = Compute_hi(species,T)
global NASA
%   species = Name of the species (example: 'H2')
%   T [Kelvin] = Temperature
%   hi [J/kg] = NASA.Specific total enthalpy 

switch species
    case 'LPG'
        i = 54;
    case 'Gasoline (High Octane)'
        i = 55;
    case 'Diesel'
        i = 56;
    case 'O2'
        i = 4;
    case 'N2'
        i = 48;
    case 'CO2'
        i = 16;
    case 'H2O'
        i = 6;
    otherwise
        i = find(strcmp({NASA.Sp.Name},species));
end
% SpeciesNames = {NASA.Sp.Name};
% i = find(strcmp(SpeciesNames,species)); % Index of the considered species
R = 8.3145; % [J/mol.K] Universal ideal gas constant
M = NASA.Sp(i).Mass; % [kg/mol]

%% Temperature threshold
if T <= NASA.Sp(i).Ts
    j=1; % T is lower than Ts -> use 1st line
else
    j=2; % T is higher than Ts -> use 2nd line
end
hi=(R/M)*(NASA.Sp(i).Pol(j,6) + NASA.Sp(i).Pol(j,1)*T + NASA.Sp(i).Pol(j,2)/2*T^2 + NASA.Sp(i).Pol(j,3)/3*T^3 + NASA.Sp(i).Pol(j,4)/4*T^4+...
    NASA.Sp(i).Pol(j,5)/5*T^5); % [J/kg] NASA.Specific total enthalpy
end

