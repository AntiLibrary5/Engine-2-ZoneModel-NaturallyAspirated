function MassMatrix = CreateMassMatrix2zones(t, y, RunCond, EngGeom,~)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

persistent ResMassFrac

P  = y(1);% [Pa]
Tu = y(2);% [K]
Tb = y(3);% [K]
Vb = y(4);% [m3]
R = 8.3145; % [J/mol.K] Universal ideal gas constant

theta = t*6*RunCond.N-360; % [CAD]
V = Vcyl(theta,EngGeom,RunCond); % [m3]

[Mu,Mb] = ComputeCombustionComposition(0, RunCond);
Vu = V-Vb;% [m3] Volume of Unburnt gases
mu = P*Vu*Mu/(R*Tu);% [kg] Mass of Unburnt gases
mb = P*Vb*Mb/(R*Tb);% [kg] Mass of Burnt gases

if theta < RunCond.SparkAngle
    ResMassFrac = mb/(mu+mb);% [-] Mass fraction of residual gases (burnt gases remaining from the previous cycle)
end

[~,~,Yi_u,Yi_b] = ComputeCombustionComposition(ResMassFrac, RunCond);
[~,cv_u] = cp_cv_mixture(Tu,Yi_u);
[~,cv_b] = cp_cv_mixture(Tb,Yi_b); % [J/kg.K] Specific heat of

MassMatrix =   [0, mu*cv_u, 0, -P;...
                0, 0, mb*cv_b, P;...
                Vu, -mu*R/Mu, 0, -P;...
                Vb, 0, -mb*R/Mb, P];
end

