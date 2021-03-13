function [ValveLift,ValveOverlap] = AdaptValveData(VD,RunCond)
%UNTITLED Summary of this function goes here
%   VD: struct containing original Valve Data
%   RunCond: struct containing Running Conditions

IntDephasing = RunCond.IntakeShift; % [CAD] Shift of the intake camshaft
ExhDephasing = RunCond.ExhaustShift;% [CAD] Shift of the exhaust camshaft
IntLiftMultiplier = RunCond.IntLiftMultiplier;
ExhLiftMultiplier = RunCond.ExhLiftMultiplier;


i_IVO = find(VD.LiftValveAdmission,1,'first'); % Find the index of 1st lift value >0
i_IVC = find(VD.LiftValveAdmission,1,'last'); % Find the index of last lift value >0
ValveLift.IntakeLifts = IntLiftMultiplier.*VD.LiftValveAdmission(i_IVO:i_IVC);% [m] Intake valve lift
ValveLift.IntakeCA = VD.AngleValveAdmission(i_IVO:i_IVC) + IntDephasing;
IVO = min(ValveLift.IntakeCA); % [CAD] Crank Angle value of IVO (Intake Valve Opening)
% IVC = max(ValveLift.IntakeCA); % [CAD] Crank Angle value of IVC (Intake Valve Closure)

i_EVO = find(VD.LiftValveExhaust,1,'first'); % Find the index of 1st lift value >0
i_EVC = find(VD.LiftValveExhaust,1,'last'); % Find the index of last lift value >0
ValveLift.ExhaustLifts = ExhLiftMultiplier.*VD.LiftValveExhaust(i_EVO:i_EVC);% [m] Intake valve lift
ValveLift.ExhaustCA = VD.AngleValveExhaust(i_EVO:i_EVC) + ExhDephasing;
% EVO = min(ValveLift.ExhCA);
EVC = max(ValveLift.ExhaustCA); % [CAD] Crank Angle value of EVC (Exhaust Valve Closure)
ValveOverlap = EVC-(720+IVO);% [CAD] Angular duration of valve overlap

% figure
% plot(IntakeCA,IntakeLift, ExhhaustCA,ExhaustLift)
% 
% figure
% plot(IntakeCA,IntakeLift)
% figure
% plot(ExhhaustCA,ExhaustLift)

end

