% clc
% clear all
% close all

% function [KnockOrNot, KI, Efficiency, IMEP] = MainScript(ON,N)
global NASA CA

dbstop if error

NASA = load('NasaThermDatFull.mat');

CA = -360:360; % [°] Crank Angles
R = 8.3145; % [J/kgK] Universal gas constant

%% Engine geometry
EngGeom.S = 80.5e-3;     % [m] Stroke
EngGeom.B = 79.5e-3;     % [m] Bore
EngGeom.L = 128e-3;      % [m] Connecting rod length
EngGeom.R = EngGeom.S/2; % [m] Crankshaft Radius
EngGeom.eta = 9.8;       % [-] Compression Ratio
EngGeom.NoIntValves = 2; % [-] Number of intake valve(s) per cylinder
EngGeom.NoExhValves = 2; % [-] Number of exhaust valve(s) per cylinder
EngGeom.IntValveDiam = 28.5e-3; % [m] Diameter intake valve(s)
EngGeom.ExhValveDiam = 23.5e-3; % [m] Diameter exhaust valve(s)

%% Running Conditions (RunCond structure array)
RunCond.N = 4000; % [RPM] Engine speed
% RunCond.N = N; % [RPM] Engine speed
RunCond.phi = 1; % fuel-air equiv. ratio
RunCond.T_plenum = 400; % [K] Temperature of the intake gases in the plenum
RunCond.P_plenum = 1.1e5; % [Pa] Pressure of the intake gases in the plenum
RunCond.T_exhaust = 1000;% [K] Temperature of the exhaust gases in the manifold
RunCond.P_exhaust = 1.75e5; % [Pa] Pressure in the exhaust manifold
RunCond.IntakeShift  = 0; % [CAD] Shift of the intake camshaft
RunCond.ExhaustShift = 0;% [CAD] Shift of the exhaust camshaft
RunCond.IntLiftMultiplier = 1;
RunCond.ExhLiftMultiplier = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RunCond.SparkAngle = -10; % [CAD] SparkAngle>0 means after TDC
% RunCond.SparkAngle = SA; % [CAD] SparkAngle>0 means after TDC
RunCond.CombDuration = 30; % [CAD] Combustion angular duration
% RunCond.CombDuration = CD; % [CAD] Combustion angular duration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Valve Lifts
ValveData = load('ValveData.mat');
[ValveLift,ValveOverlap] = AdaptValveData(ValveData,RunCond);
%% Fuel
FuelType = 'Gasoline (High Octane)';%'LPG','Gasoline (High Octane)';
[RunCond.Fuel] = ComputeFuelProperties(FuelType, RunCond.phi);
% RunCond.Fuel.OctaneNumber = ON;
%% Definition of the limits of the different strokes
i_SoCp = find(CA==max( round(max(ValveLift.IntakeCA)), -180));  % Index of Start of Compression (either Intake Valve Closure or BDC)
% i.EoCp = find(CA==0)-1;   % Index of End of Compression
i_SoCb = find(CA==RunCond.SparkAngle);% Index of Start of Combustion
i_EoCb = find(CA==RunCond.SparkAngle+RunCond.CombDuration);     % Index of End of Combustion
% i.SoExp = find(CA==0)+1;  % Index of Start of Expansion
% i.EoExp = find(CA==180)-1;% Index of End of Expansion
% i.SoExh = find(CA==180);  % Index of Start ofCQ Exhaust

%% Preallocation of the memory


V = Vcyl(CA,EngGeom,RunCond);       % [m3] Volume of the chamber

tFinal = 2*60/RunCond.N; % [s] Total duration of the complete cycle (2 engine rotations for a 4-stroke engine)
tStep = tFinal/(numel(CA)-1); % [s] Time step
tSpan = 0 : tStep : tFinal; % [s] Time span

P_init = RunCond.P_plenum;
Tu_init = RunCond.T_plenum;
Tb_init = RunCond.T_exhaust;
Vb_init = 0.01*min(V);% [m3] Initial Volume of Burnt gases

y0 = [P_init ; Tu_init ; Tb_init ; Vb_init]; % Initial values
options = odeset('RelTol',1e-5,'Mass',@CreateMassMatrix2zones);
[t,y] = ode45(@CreateForceMatrix2zones,tSpan,y0,options,RunCond,EngGeom,ValveLift);

P  = y(:,1);
Tu = y(:,2);
Tb = y(:,3);
Vb = y(:,4);

[Mu,Mb] = ComputeCombustionComposition(0, RunCond);
Vu = V'-Vb;
mu = P.*Vu./((R./Mu).*Tu);% [kg] Mass of unburnt gases
mb = P.*Vb./((R./Mb).*Tb);% [kg] Mass of burnt gases
iter = 1;
ResMassFrac(iter) = mb(i_SoCb)/(mu(i_SoCb)+mb(i_SoCb));

Variation_P(iter) = 100*abs( P(1)/P(end)-1 );% [%]
Variation_Tu(iter) = 100*abs( Tu(1)/Tu(end)-1 );% [%]
Variation_Tb(iter) = 100*abs( Tb(1)/Tb(end)-1 );% [%]
Variation_Vb(iter) = 100*abs( Vb(1)/Vb(end)-1 );% [%]
Variation_mu(iter) = 100*abs( mu(1)/mu(end)-1 );% [%]
Variation_mb(iter) = 100*abs( mb(1)/mb(end)-1 );% [%]

%% Iterative loop
while Variation_P(iter) > 5 || Variation_Tu(iter) > 5 || Variation_Tb(iter) > 5 || Variation_Vb(iter) > 5 || Variation_mu(iter) > 5 || Variation_mb(iter) > 5
    iter = iter+1;
    P_init = P(end);
    Tu_init = Tu(end);
    Tb_init = Tb(end);
    Vb_init = Vb(end);% [m3] Initial Volume of Burnt gases
    y0 = [P_init ; Tu_init ; Tb_init ; Vb_init]; % Initial values
    
    options = odeset('RelTol',1e-5,'Mass',@CreateMassMatrix2zones);
    [t,y] = ode45(@CreateForceMatrix2zones,tSpan,y0,options,RunCond,EngGeom,ValveLift);
    
    P  = y(:,1);
    Tu = y(:,2);
    Tb = y(:,3);
    Vb = y(:,4);
    
    Variation_P(iter) = 100*abs( P(1)/P(end)-1 );% [%]
    Variation_Tu(iter) = 100*abs( Tu(1)/Tu(end)-1 );% [%]
    Variation_Tb(iter) = 100*abs( Tb(1)/Tb(end)-1 );% [%]
    Variation_Vb(iter) = 100*abs( Vb(1)/Vb(end)-1 );% [%]
    Variation_mu(iter) = 100*abs( mu(1)/mu(end)-1 );% [%]
    Variation_mb(iter) = 100*abs( mb(1)/mb(end)-1 );% [%]

    Vu = V'-Vb;
    
    mu = P.*Vu./((R./Mu).*Tu);% [kg] Mass of unburnt gases
    mb = P.*Vb./((R./Mb).*Tb);% [kg] Mass of burnt gases
    ResMassFrac(iter) = mb(i_SoCb)/(mu(i_SoCb)+mb(i_SoCb));
end

[KnockOrNot, KI] = CheckKnock(t(i_SoCp:i_EoCb), P(i_SoCp:i_EoCb), Tu(i_SoCp:i_EoCb), RunCond.Fuel.OctaneNumber);
%% Performance computation
MaxP = max(P); % [Pa]
MaxTu = max(Tu); % [K]
MaxTb = max(Tb); % [K]

Wcycle = trapz(V,P); % [J]
IMEP = Wcycle/(max(V)-min(V)); % [Pa]
PwrOut = 1e-3*Wcycle*RunCond.N/120; % [kW] Power Output
SpecPwrOut = PwrOut/( 1e3*(max(V)-min(V)) );% [kW/L] Specific Power Output

mFresh = mu(i_SoCb);% [kg] Mass of fresh gases trapped into the cylinder
mRes = mb(i_SoCb);% [kg] Mass of residuals gases
mTrapped = mu(i_SoCb)+mb(i_SoCb);% [kg] Total Trapped Mass
ResMassFrac = mRes/mTrapped;% [-] Mass fraction of residual gases (burnt gases remaining from the previous cycle)
mFuel = mFresh/(1+RunCond.Fuel.AFRst/RunCond.phi); % [kg] Mass of fuel consumed in each cycle
Efficiency = Wcycle/(mFuel*RunCond.Fuel.LHV); % [-] Efficiency of the cycle
mFreshTheor = RunCond.P_plenum*(max(V)-min(V))*Mu/(R*RunCond.T_plenum);% [kg] Theoretical Mass of intake fresh gases
VolEff = mFresh/mFreshTheor; % [-] Volumetric Efficiency

dQu_dt      = NaN(numel(t),1);
dQb_dt      = NaN(numel(t),1);
dmIntake_dt = NaN(numel(t),1);
dmExhaust_dt= NaN(numel(t),1);
A_Int       = NaN(numel(t),1);
A_Exh       = NaN(numel(t),1);
for i = 1:numel(t)
    [~, dQu_dt(i), dQb_dt(i), dmIntake_dt(i), dmExhaust_dt(i), A_Int(i), A_Exh(i)] = CreateForceMatrix2zones(tSpan(i),y(i,:),RunCond,EngGeom,ValveLift);
end
