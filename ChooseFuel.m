function [FuelProperties] = ChooseFuel(FuelName)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
switch FuelName
    case 'Gasoline (High Octane)'
        FuelProperties.Type = 'Gasoline (High Octane)';
        FuelProperties.x = 7.77; % carbon fraction of fuel (Gasoline/Diesel/LPG)
        FuelProperties.y = 13.08; % hydrogen fraction of fuel (Gasoline/Diesel/LPG)
%         FuelProperties.phi=1; % fuel-air equiv. ratio (Gasoline/Diesel/LPG)
    case 'Diesel'
        FuelProperties.Type = 'Diesel';
        FuelProperties.x=10.82; % carbon fraction of fuel (Gasoline/Diesel/LPG)
        FuelProperties.y=18.73; % hydrogen fraction of fuel (Gasoline/Diesel/LPG)
%         FuelProperties.phi=0.8; % fuel-air equiv. ratio (Gasoline/Diesel/LPG)
    case 'LPG'
        FuelProperties.Type = 'LPG';
        FuelProperties.x=3.53; % carbon fraction of fuel (Gasoline/Diesel/LPG)
        FuelProperties.y=8.98; % hydrogen fraction of fuel (Gasoline/Diesel/LPG)
%         FuelProperties.phi=1; % fuel-air equiv. ratio (Gasoline/Diesel/LPG)
    otherwise
        disp('Unknown fuel')
end
end

