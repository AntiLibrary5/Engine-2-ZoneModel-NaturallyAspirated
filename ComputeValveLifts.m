function [IntakeLift,ExhaustLift] = ComputeValveLifts(theta,VL)
%UNTITLED Summary of this function goes here
%   theta [CAD]: crank angle value (between -360 and +360°)
%   VL (ValveLift): struct containing the lift laws

%% Intake
IVO = min(VL.IntakeCA); % [CAD] Crank Angle value of IVO (Intake Valve Opening)
if IVO<-360 && theta>IVO+720
    thetaInt = theta-720;
else
    thetaInt = theta;
end
IntakeLift = interp1(VL.IntakeCA, VL.IntakeLifts, thetaInt,'linear',0);

%% Exhaust
EVC = max(VL.ExhaustCA); % [CAD] Crank Angle value of EVC (Exhaust Valve Closure)
if EVC>360 && theta<EVC-720
    thetaExh = theta+720;
else
    thetaExh = theta;
end
ExhaustLift = interp1(VL.ExhaustCA,VL.ExhaustLifts,thetaExh,'linear',0);
end

% figure
% plot(thetas,IntakeLift,thetas,ExhaustLift)
% hold on
% plot(IntakeCA,IntakeLifts,ExhaustCA,ExhaustLifts)