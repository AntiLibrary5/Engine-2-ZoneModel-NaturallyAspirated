% t=0:100Tu = 500:600;
% P(1:101) = 6e6;
% ON = 93;

function [KnockOrNot, KI] = CheckKnock(t,P,Tu,ON)
%UNTITLED4 Summary of this function goes here
%   t [s] : time
%   P [Pa] : in-cylinder Pressure during combustion
%   Tu [K] : Temperature of the Unburnt gases during combustion
%   ON : Octane Number of the fuel

tau = 17.68e-3 * (ON/100)^3.402 * (P/1.01325e5).^-1.7 .*exp(3800./Tu);
KI = trapz(t,1./tau);
if KI < 1
    KnockOrNot = false;
else
    KnockOrNot = true;
end

