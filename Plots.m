% Plots
figure
plot(CA,P*1e-5)
hold on

figure
plot(CA,Tu, CA,Tb)
legend('T_u','T_b')

figure
plot(V,P)
yyaxis right
plot(V,T)
legend('Pressure','Temperature')

figure
subplot(2,1,1)
plot(CA,P)
yyaxis right
plot(CA,Tu, CA,Tb)
legend('Pressure','T_u','T_b')
subplot(2,1,2)
plot(CA,mu, CA,mb)
legend('m_{unburnt}','m_{burnt}')

set(0,'DefaultAxesFontSize',18.0)
figure
subplot(2,1,1)
plot(CA,dmIntake_dt, CA,dmExhaust_dt,'LineWidth',2)
yyaxis right
plot(CA,A_Int, CA,A_Exh,'LineWidth',2)
grid on; grid minor;
ylabel('Mass flow rate through the valves [kg/s]')
xlabel('Crank angle [CAD]')
legend('Intake Mass Flow','Exhaust Mass Flow','Intake Area','Exhaust Area')
subplot(2,1,2)
plot(CA,mu, CA,mb,'LineWidth',2)
grid on; grid minor;

%% Variations
figure
plot(iter,Variation_P, iter,Variation_Tu, iter,Variation_Tb, iter,Variation_Vb, iter,Variation_mu, iter,Variation_mb)
legend('Variation_P','Variation_Tu','Variation_Tb','Variation_Vb','Variation_mu','Variation_mb')






