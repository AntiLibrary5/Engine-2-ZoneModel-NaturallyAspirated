function [ForceMatrix, dQu_dt, dQb_dt, dmIntake_dt, dmExhaust_dt, A_Int, A_Exh] = CreateForceMatrix2zones(t,y,RunCond,EngGeom,ValveLift)

persistent mResiduals mU_left ResMassFrac
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

P  = y(1);% [Pa]
Tu = y(2);% [K]
Tb = y(3);% [K]
Vb = y(4);% [m3]

theta = t*6*RunCond.N-360; % [CAD]
[V,dV_dt] = Vcyl(theta,EngGeom,RunCond);

R = 8.3145; % [J/mol.K] Universal ideal gas constant

[Mu,Mb] = ComputeCombustionComposition(0, RunCond);
Vu = V-Vb;% [m3] Volume of Unburnt gases
mu = P*Vu*Mu/(R*Tu);% [kg] Mass of Unburnt gases
mb = P*Vb*Mb/(R*Tb);% [kg] Mass of Burnt gases

if theta < RunCond.SparkAngle
    mResiduals = mb;% [kg] Mass of residual gases (burnt gases remaining from the previous cycle)
    ResMassFrac = mResiduals/(mu+mb);% [-]
    mU_left = 0.01*mu;% [kg] Mass of fresh gases that will not burn
end

[~,~,Yi_u,Yi_b] = ComputeCombustionComposition(ResMassFrac, RunCond);

%% Unburnt (u) gases properties
hi_u = [Compute_hi(RunCond.Fuel.Type,Tu) Compute_hi('O2',Tu) Compute_hi('N2',Tu) Compute_hi('CO2',Tu)...
    Compute_hi('H2O',Tu)]; % [J/kg] Specific enthalpy of [fuel O2 N2 CO2 H2O] @ Tu
hu = sum(Yi_u.*hi_u); % Average Specific Enthalpy of the Unburnt gases


eu = hu - P*Vu/mu;% [J/kg] Specific internal energy of Unburnt gases

%% Burnt (b) gases properties
hi_b = [Compute_hi(RunCond.Fuel.Type,Tb) Compute_hi('O2',Tb) Compute_hi('N2',Tb) Compute_hi('CO2',Tb)...
    Compute_hi('H2O',Tb)]; % [J/kg] Specific enthalpy of the species [fuel O2 N2 CO2 H2O] @ T3
hb = sum(Yi_b.*hi_b);% [J/kg] Average Specific Enthalpy of the Burnt gases        

eb = hb - P*Vb/mb;% [J/kg] Specific internal energy of Burnt gases

%% Mass flows going through the valves
% hExh = hb; % Specific Enthalpy of the (burnt or unburnt) gases going through the Exhaust valve
% hInt = hu; % Specific Enthalpy of the (burnt or unburnt) gases going through the Intake valve

[dmIntake_dt, dmExhaust_dt, A_Int, A_Exh] = ComputeMassFlows(theta, P, Tu, Tb, ResMassFrac, EngGeom, RunCond, ValveLift);

if dmIntake_dt >= 0 % normal flow through intake valve(s)
    dmuInt_dt = dmIntake_dt;
    T_Int = RunCond.T_plenum;% [K] Temperature of the gases in the plenum
    hu_Int_i = [Compute_hi(RunCond.Fuel.Type,T_Int) Compute_hi('O2',T_Int) Compute_hi('N2',T_Int)...
        Compute_hi('CO2',T_Int) Compute_hi('H2O',T_Int)]; % [J/kg] Specific enthalpy of the species [fuel O2 N2 CO2 H2O] @ T_int
    hu_Int = sum(Yi_u.*hu_Int_i);% [J/kg] Specific Enthalpy of the Unburnt gases going through the Intake valve(s) computed with intake port temperature
    dmbInt_dt = 0;
    hb_Int = 0;
else% reverse flow through intake valve(s)
    dmuInt_dt = dmIntake_dt*mu/(mu+mb);
    hu_Int = hu;% [J/kg] Specific Enthalpy of the Unburnt gases going through the Intake valve(s)
    dmbInt_dt = dmIntake_dt*mb/(mu+mb);
    hb_Int = hb;% [J/kg] Specific Enthalpy of the Burnt gases going through the Intake valve(s)
end

if dmExhaust_dt <= 0 % normal flow through exhaust valve(s)
    dmuExh_dt = dmExhaust_dt*mu/(mu+mb);
    hu_Exh = hu;% [J/kg] Specific Enthalpy of the Unburnt gases going through the Exhaust valve(s)
    dmbExh_dt = dmExhaust_dt*mb/(mu+mb);
    hb_Exh = hb;% [J/kg] Specific Enthalpy of the Burnt gases going through the Exhaust valve(s)
else % reverse flow through exhaust valve(s)
    dmuExh_dt = 0;
    hu_Exh = 0;
    dmbExh_dt = dmExhaust_dt;
    T_Exh = RunCond.T_exhaust;
    hu_Exh_i = [Compute_hi(RunCond.Fuel.Type,T_Exh) Compute_hi('O2',T_Exh) Compute_hi('N2',T_Exh)...
        Compute_hi('CO2',T_Exh) Compute_hi('H2O',T_Exh)]; % [J/kg] Specific enthalpy of the species [fuel O2 N2 CO2 H2O] @ T_Exh
    hb_Exh = sum(Yi_b.*hu_Exh_i);% [J/kg] Specific Enthalpy of the Burnt gases going through the Exhaust valve(s), computed with exhaust port temperature
end


%% Combustion
[~, dYb_dTheta] = CombustionModel(theta, RunCond);
Omega = RunCond.N*6; % [CAD/s] Engine speed
mTotal = mu + mb;% [kg] Total mass into the cylinder
dmbReac_dt = (mTotal-mResiduals-mU_left)*dYb_dTheta*Omega;% [kg/s] 



dmb_dt = (dmbExh_dt+dmbInt_dt) + dmbReac_dt; % [kg/s] (Gas Exch.) + (Reac.)
dmuReac_dt = -dmbReac_dt;
dmu_dt = (dmuExh_dt+dmuInt_dt) + dmuReac_dt; % [kg/s] (Gas Exch.) + (Reac.)

%% Heat losses
[dQu_dt, dQb_dt] = ComputeHeatLosses(Tu, Tb, P, RunCond, EngGeom, theta, V, Vb);

% figure(17)
% hold on
% plot(theta,dQu_dt,'o')
% title('dQu_{dt}')
% 
% figure(18)
% hold on
% plot(theta,dQb_dt,'o')
% title('dQb_{dt}')

%% Force Matrix
RH1 = (-eu*dmu_dt + dQu_dt -P*dV_dt + hu_Exh*dmuExh_dt + hu_Int*dmuInt_dt -hu*dmbReac_dt);
RH2 = (-eb*dmb_dt + dQb_dt + hb_Exh*dmbExh_dt + hb_Int*dmbInt_dt +hu*dmbReac_dt);
RH3 = (-P*dV_dt + Tu*R/Mu*dmu_dt);
RH4 = (Tb*R/Mb*dmb_dt);

ForceMatrix = [RH1; RH2; RH3; RH4];

% figure(17)
% hold on
% plot(theta,mb,'o')
% ylim([0 0.0006])
% title('m_b')
% 
% figure(18)
% hold on
% plot(theta,mu,'o')
% ylim([0 0.0006])
% title('m_u')
end

