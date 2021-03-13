function [dmIntake_dt, dmExhaust_dt, A_Int, A_Exh] = ComputeMassFlows(theta, P, Tu, Tb, ResMassFrac, EngGeom, RunCond, ValveLift)
%Computes the massflows going through the intake and exhaust valves.
% Input arguments:
%   theta [CAD]: Considered crank angle
%   p [Pa] : In-cylinder pressure
%   Tu [K]: Temperature of the in-cylinder Uburnt gases
%   Tb [K]: Temperature of the in-cylinder Burnt gases
%   EngGeom: Struct containing the Engine specifications
%   RunCond: Struc containig the running conditions of the engine
%
% Output arguments:
%   dmIntake_dt [kg/s]
%   dmExhaust_dt [kg/s]

%% Load / define required data

R = 8.3145; % [J/mol.K] Universal ideal gas constant

pIntake  = RunCond.P_plenum;
pExhaust = RunCond.P_exhaust;

NoIntValves = EngGeom.NoIntValves; % Number of intake valves per cylinder
NoExhValves = EngGeom.NoExhValves; % Number of exhaust valves per cylinder
IntValveDiam = EngGeom.IntValveDiam; % [m] Diameter of the intake valve(s)
ExhValveDiam = EngGeom.ExhValveDiam; % [m] Diameter of the exhaust valve(s)

[IntakeLift,ExhaustLift] = ComputeValveLifts(theta,ValveLift); % [m] Current lift of the intake and exhaust valves

if isempty(IntakeLift)
    IntakeLift = 0;
end

if isempty(ExhaustLift)
    ExhaustLift = 0;
end


%% Intake valves
if P <= pIntake % "Normal" case (Gases are going into the cylinder)
%     disp("Intake normal Flow")
    coef = 1; % [-] In-cylinder mass should increase
    p0 = pIntake;% [Pa] Upstream pressure
    p1 = P;% [Pa] Downstream pressure
    T0 = RunCond.T_plenum;
    [M0,~,Yi] = ComputeCombustionComposition(ResMassFrac, RunCond);
else            % "Reverse flow" case (Gases are going out from the cylinder)
%     disp("Intake reverse Flow")
    coef = -1; % [-] In-cylinder mass should decrease
    p0 = P;% [Pa] Upstream pressure
    p1 = pIntake;% [Pa] Downstream pressure
    if theta < -270 % Gases are considered as burnt gases
        T0 = Tb;
        [~,M0,~,Yi] = ComputeCombustionComposition(ResMassFrac, RunCond);
    else             % Gases are considered as unburnt gases
        T0 = Tu;
        [M0,~,Yi] = ComputeCombustionComposition(ResMassFrac, RunCond);
    end
end

Cd_Int = 0.7; % [-] Discharge coefficient at the intake valve

A_IntValve = pi*IntValveDiam*IntakeLift; % [m2] Flow cross section Area at the level of the intake valve
A_IntPort = pi*(0.9*IntValveDiam/2)^2;% [m2] Flow cross section Area at the level of the intake port
A_Int = NoIntValves*min(A_IntValve,A_IntPort);% [m2] Area of the smallest section crossed by the intake flow

rho0 = p0*M0/R/T0; % [kg/m3] Density of the gases
[~,~,g] = cp_cv_mixture(T0,Yi); % [-] Gamma (cp/cv)

if p1/p0 >= (2/(g+1))^(g/(g-1))
    dmIntake_dt = coef*Cd_Int*A_Int * rho0*(p1/p0)^(1/g)*sqrt((2*g/(g-1))...
        *(p0/rho0)*(1-(p1/p0)^((g-1)/g))); % [kg/s] Mass flow rate going through the intake valve(s)
else
    a0 = sqrt(g*(R/M0)*T0); % [m/s] Speed of sound
    dmIntake_dt = coef*Cd_Int*A_Int * rho0*a0*(2/(g+1))^((g+1)/(2*g-2)); % [kg/s] Mass flow rate going through the intake valve(s)
end

%% Exhaust valves
if P >= pExhaust%   "Normal" case (Gases are going out from the cylinder)
%     disp("Exhaust normal Flow")
    coef = -1; % [-] In-cylinder mass should decrease
    p0 = P;% [Pa] Upstream pressure
    p1 = pExhaust;% [Pa] Downstream pressure
else%               "Reverse flow" case (Gases are going into the cylinder)
%     disp("Exhaust reverse Flow")
    coef = 1;%      [-] In-cylinder mass should increase
    p0 = pExhaust;% [Pa] Upstream pressure
    p1 = P;%        [Pa] Downstream pressure
end

% For now, we consider that gases going through the exhaust valve(s) are
% only burnt gases.
T0 = Tb;
[~,M0,~,Yi] = ComputeCombustionComposition(ResMassFrac, RunCond); % Mass fractions in burnt gases
[~,~,g] = cp_cv_mixture(T0,Yi); % [-] Gamma (cp/cv)

Cd_Exh = 0.6; % [-] Discharge coefficient at the exhaust valve
A_ExhValve = pi*ExhValveDiam*ExhaustLift; % [m2] Flow cross section Area at the level of the exhaust valve
A_ExhPort = pi*(0.9*ExhValveDiam/2)^2;% [m2] Flow cross section Area at the level of the exhaust port
A_Exh = NoExhValves*min(A_ExhValve,A_ExhPort);% [m2] Area of the smallest section crossed by the exhaust flow

rho0 = p0*M0/R/T0; % [kg/m3] Density of the gases

if p1/p0 >= (2/(g+1))^(g/(g-1))
    dmExhaust_dt = coef*Cd_Exh*A_Exh * rho0*(p1/p0)^(1/g)*sqrt((2*g/(g-1))...
        *(p0/rho0)*(1-(p1/p0)^((g-1)/g))); % [kg/s] Mass flow rate going through the exhaust valve(s)
else
    a0 = sqrt(g*(R/M0)*T0); % [m/s] Speed of sound
    dmExhaust_dt = coef*Cd_Exh*A_Exh * rho0*a0*(2/(g+1))^((g+1)/(2*g-2)); % [kg/s] Mass flow rate going through the exhaust valve(s)
end

% figure(21)
% hold on
% plot(theta,dmIntake_dt,'o')
% title('dmIntake')
% 
% figure(22)
% hold on
% plot(theta,dmExhaust_dt,'o')
% title('dmExhaust')

end



