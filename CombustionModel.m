function [Yb, dYb_dTheta] = CombustionModel(theta, RunCond)
%UNTITLED4 Summary of this function goes here
%   Yb [-] : Mass fraction of burnt gases
%   dYb_dTheta [CAD^-1] : Heat release rate
deltaTheta = RunCond.CombDuration; % [CAD] Combustion angular duration
thetaIgnition = RunCond.SparkAngle; % [CAD] SparkAngle>0 means after TDC
a1 = 5;
a2 = 2.5;
if theta > thetaIgnition
    Yb = 1-exp( -a1*((theta-thetaIgnition)/deltaTheta)^(a2+1) ); % Mass fraction of burnt gases
    dYb_dTheta = a1*(a2+1)/deltaTheta * ( (theta-thetaIgnition)/deltaTheta )^a2...
        *exp( -a1*((theta-thetaIgnition)/deltaTheta)^(a2+1) ); % [CAD^-1] Heat release rate
else
    Yb = 0;
    dYb_dTheta = 0;
end

% figure(13)
% title('dYb_dTheta = f( theta )')
% hold on
% plot(theta,dYb_dTheta,'o')

end

