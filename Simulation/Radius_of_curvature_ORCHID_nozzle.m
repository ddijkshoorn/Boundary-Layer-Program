%% Radius of Curvature - ORCHID nozzle
% radius of curvature: take analytical relation from Wikipedia

% Last checked on 9-03-2022

close all
clear all
clc

%% Settings
Color = [0 0 0];
LineStyle = {'-','--','-.',':','-',':'};
Marker = {'none','none','none','none','none','none'}; % Marker = {'+','x','s','d','o'};
LineWidth = 0.5;
FontSize = 11;

%% Format figure
% figure 11: Wall Turbulent Prandtl-number PrTw
h1 = figure;
ax1 = gca;

ax1.FontSizeMode = 'manual';
ax1.FontSize = FontSize; %%
ax1.FontWeight = 'normal';
ax1.TickLabelInterpreter = 'latex';
ax1.XScale = 'linear';
ax1.YScale = 'linear';

ax1.XLabel.Interpreter = 'latex';
ax1.XLabel.FontSize = FontSize; %%
ax1.XLabel.FontWeight = 'normal';
ax1.XLabel.Color = Color;
ax1.XLabel.String = '$X \, [\mathrm{m}]$';
ax1.XLim = [0 0.1];
ax1.XLimMode = 'manual';
ax1.XTickLabelMode = 'manual';
ax1.XTickMode = 'manual';
ax1.XTick = [0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.10];
ax1.XTickLabel = {'0'; '0.01'; '0.02'; '0.03'; '0.04'; '0.05'; '0.06'; '0.07'; '0.08'; '0.09';'0.10'};

ax1.YLabel.Interpreter = 'latex';
ax1.YLabel.FontSize = FontSize; %%
ax1.YLabel.FontWeight = 'normal';
ax1.YLabel.Color = Color;
ax1.YLabel.String = '$\frac{\delta}{R_{\mathrm{c}}} \, [-]$';
ax1.YLim = [0 0.014];
ax1.YLimMode = 'manual';
ax1.YTickLabelMode = 'manual';
ax1.YTickMode = 'manual';
ax1.YTick = [0 0.002 0.004 0.006 0.008 0.010 0.012 0.014];
ax1.YTickLabel = {'0'; '0.002'; '0.004'; '0.006'; '0.008'; '0.010'; '0.012'; '0.014'};
ax1.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',FontSize)
box on
hold on
%%

cd ..
% load geometry with RefProp laminar simulation results
load('Stored_Sim_Data/ORCHID_Nozzle_RefProp_PsE_Interpolated_lam')

% calculation of radius of curvature
x = INP.x;
y = INP.y;
yp = GRADNT(y,x);
p1 = pchip([x(1:32) x(39:44) x(64:70)],[yp(1:32) yp(39:44) yp(64:70)],x(1:70));
figure
hold on
plot(yp)
plot(p1);
yp(1:70) = p1;
legend('yp','p1')
% plot(yp)
p2 = pchip([x(1:31) x(44:53) x(71:74)],[yp(1:31) yp(44:53) yp(71:74)],x(1:74));
plot(p2);
yp(1:74) = p2;
ypp = GRADNT(yp,x);
% figure;plot(x,ypp)
R_c = abs(((1 + yp.^2).^(3/2))./ypp);

[A, B] = min(R_c);

fprintf('Minimum radius of curvature is %2.4f meter at station %g.\n',A,B);

% Laminar
delta_Rc_lam = BLC.delta./R_c;
delta_lam_scaled = BLC.delta/BLC.delta(end)*0.013;

% Turbulent
% load('Stored_Sim_Data/ORCHID_Nozzle_RefProp_PsE_Smoothed_turb')
load('Stored_Sim_Data/ORCHID_Nozzle_RefProp_PsE_Interpolated_turb')
% why Smoothed and not Interpolated? Mistake?

delta_Rc_turb = BLC.delta./R_c;

% plot results
figure
hold on
plot(x,y)
plot(x,yp)
plot(x,ypp)
legend('y','yp','ypp')
title('geometry and its derivatives')
xlabel('X [m]')
ylabel('y, y'', y''''')

figure(h1)
hold on
plot(x,delta_Rc_lam,'DisplayName','$\frac{\delta}{R_{\mathrm{c}}}$ Lam')
plot(x,delta_Rc_turb,'DisplayName','$\frac{\delta}{R_{\mathrm{c}}}$ Turb')
% add
plot(x,y/y(1)*0.013,'DisplayName','Nozzle Geom scaled') % scaled nozzle geometry
plot(x,delta_lam_scaled,'DisplayName','$\delta_{\mathrm{Lam}}$ scaled')
plot(x,BLC.delta/BLC.delta(end)*0.013,'DisplayName','$\delta_{\mathrm{Turb}}$ scaled')
% legend('\delta/R_c Lam','\delta/R_c Turb','scaled Nozzle Geom','scaled \delta Lam','scaled \delta Turb','Location','north')

% centrifugal forces
figure
hold on
plot(x,EDG.rhoE.*EDG.UE.^2./R_c)
plot(x,EDG.PsE)
title('Centrifugal forces')
xlabel('x [m]')
ylabel('F [N/m3]')
legend('F_c','PsE')

%% Add Legend

figure(h1) % Turbulent Prandtl-number at the wall %%% no legend to show
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'north');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

%% Calculate centrifugal force profile at station 124:
Fc = EDG.rhoE(124).*EDG.UE(124).^2./SOL{124}.c.*SOL{124}.u.^2./(R_c(124) + BLC.y(1:length(SOL{124}.u),124));

figure
hold on
plot(BLC.y(1:length(SOL{124}.u),124),Fc)
xlabel('Y [m]')
ylabel('Fc [N/m3]')
title('centrifugal forces at station 124')

Int = cumtrapz(BLC.y(1:length(SOL{124}.u),124),Fc);

figure
hold on
plot(BLC.y(1:length(SOL{124}.u),124),Int)
xlabel('Y [m]')
ylabel('Int(Fc) [N/m2]')
title('Centrifugal Pressure correction inside BL at station 124')

pY = EDG.PsE(124) - Int; %%% ??? %%%%%

figure
hold on
plot(BLC.y(1:length(SOL{124}.u),124),pY)
xlabel('Y [m]')
ylabel('P [N/m2]')
title('Pressure inside BL at station 124')

%% Save figure
% keyboard
saveas(h1,'./Thesis_Report_Figures/Radius_curv','epsc')

fprintf('Figure generated and saved successfully\n')
