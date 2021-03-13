function [cp,cv,gamma] = cp_cv_mixture(T,Yi)
global NASA
%Computes the cp and cv average values of a given gases mixture.
%   T [K] : Temperature of the gases
%   Yi [-] : Mass fractions of [fuel O2 N2 CO2 H2O] in the mixture

%% Compute cv
cv = [Compute_cv('Gasoline (High Octane)',T) Compute_cv('O2',T) Compute_cv('N2',T) Compute_cv('CO2',T) Compute_cv('H2O',T)];
cv = sum(Yi.*cv); % [J/kgK] Specific heat at constant volume

%% Compute cp
cp = [Compute_cp('Gasoline (High Octane)',T) Compute_cp('O2',T) Compute_cp('N2',T) Compute_cp('CO2',T) Compute_cp('H2O',T)];
cp = sum(Yi.*cp); % [J/kgK] Specific heat at constant pressure

gamma=cp/cv;

    function [cp] = Compute_cp(species,T)
        %   species = name of the species (example: 'H2')
        %   T = temperature [Kelvin]
        %   cp = specific heat at constant pressure [J/kgK]
        switch species
            case 'LPG'
                index = 54;
            case 'Gasoline (High Octane)'
                index = 55;
            case 'Diesel'
                index = 56;
            case 'O2'
                index = 4;
            case 'N2'
                index = 48;
            case 'CO2'
                index = 16;
            case 'H2O'
                index = 6;
            otherwise
                index = find(strcmp({NASA.Sp.Name},species));
        end
        % SpeciesNames = {NASA.Sp.Name};
        % index = find(strcmp(SpeciesNames,species));
        R = 8.3145;
        M = NASA.Sp(index).Mass; %kg/mol
        
        if T<=NASA.Sp(index).Ts
            LineNo=1;
        else
            LineNo=2;
        end
        
        cp = (R/M)*((NASA.Sp(index).Pol(LineNo,1))+NASA.Sp(index).Pol(LineNo,2)*T+NASA.Sp(index).Pol(LineNo,3)*T^2+...
            NASA.Sp(index).Pol(LineNo,4)*T^3+NASA.Sp(index).Pol(LineNo,5)*T^4);
    end

    function [cv] = Compute_cv(species,T)
        %UNTITLED Summary of this function goes here
        %   species = name of the species (example: 'H2')
        %   T [K] = temperature
        %   cv [J/kgK] = specific heat at constant volume
        switch species
            case 'LPG'
                index = 54;
            case 'Gasoline (High Octane)'
                index = 55;
            case 'Diesel'
                index = 56;
            case 'O2'
                index = 4;
            case 'N2'
                index = 48;
            case 'CO2'
                index = 16;
            case 'H2O'
                index = 6;
            otherwise
                index = find(strcmp({NASA.Sp.Name},species));
        end
        % SpeciesNames = {NASA.Sp.Name};
        % index = find(strcmp(SpeciesNames,species));
        
        R = 8.3145;
        M = NASA.Sp(index).Mass; %kg/mol
        
        if T<=NASA.Sp(index).Ts
            LineNo=1;
        else
            LineNo=2;
        end
        
        cv=(R/M)*((NASA.Sp(index).Pol(LineNo,1)-1)+NASA.Sp(index).Pol(LineNo,2)*T+NASA.Sp(index).Pol(LineNo,3)*T^2+NASA.Sp(index).Pol(LineNo,4)*T^3+...
            NASA.Sp(index).Pol(LineNo,5)*T^4);
    end
end
