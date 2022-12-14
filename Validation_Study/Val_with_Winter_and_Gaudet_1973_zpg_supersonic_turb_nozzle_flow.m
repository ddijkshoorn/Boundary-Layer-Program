%% MAIN-file
% plotting case data
% last run 04-05-2020
% Adapted slightly for upload and checked: 03-03-2022

clear all
close all
clc

pause(0.1)

%% Input
% load('../Data/NozzleGeometryFlexibleWalledNozzleCaseWinter_Gaudet_1973')
load('./Data/NozzleGeometryFlexibleWalledNozzleCaseWinter_Gaudet_1973')
Geom = 0.0254.*Geom; % convert from Inches to meter

%% Settings
Color = [0 0 0];
LineStyle = {'-','--',':','-.'};
Marker = {'o','^','s','d','p','v','h','>','+','x','*'}; %%% open symbols (first four (fill), second four (fill), last three)
% Marker = {'+','x','s','d','o'};
% Marker = {'o','^','s','d'};
LineWidth = 0.5;
FontSize = 11;
LgndPos1 = [0.20 0.55 0.20 0.20]; % normalized position
LgndPos2 = [0.55 0.3 0.3 0.2];
LgndPos3 = [0.6 0.2 0.2 0.2];
LegendName = {'Ma 0.2','Ma 1.4','Ma 2.2','Ma 2.8'};
LegendName2 = {'','','','$\delta^{\ast}$','','','','$\theta$'};
LegendName3 = {'','','','$c_{\mathrm{w}}$','','','','$C_{\mathrm{w}}$'};

% figure 1 (two y-axes) - Mach-density plot (pressure history)
% n = 3;
% h1 = figure(n);
h0 = figure;

ax0 = gca;
% ax.FontSizeMode = 'manual';
% ax.FontSize = FontSize; % werkt niet

yyaxis left
ax0.YAxis(1).Color = Color;
ax0.YAxis(1).FontSize = FontSize;
ax0.YAxis(1).Limits = [0 12];
ax0.YAxis(1).TickLabelsMode = 'manual';
ax0.YAxis(1).TickValues = [0 2 4 6 8 10 12];
ax0.YAxis(1).TickLabels = {'0'; '2.0'; '4.0'; '6.0'; '8.0'; '10.0'; '12.0'};
ax0.YAxis(1).TickLabelInterpreter = 'latex';
ax0.YAxis(1).Label.FontSize = FontSize; % geen effect
ax0.YAxis(1).Label.FontWeight = 'normal';
ax0.YAxis(1).Label.Color = Color;
ax0.YAxis(1).Label.LineWidth = LineWidth;
ax0.YAxis(1).Label.Interpreter = 'latex';
ax0.YAxis(1).Label.String = '$\rho_{\mathrm{e}} / \rho_{\mathrm{e,outlet}} \, [-]$';

yyaxis right % or: yyaxis right; set(ax.YAxis(2).Label, 'Interpreter', 'latex', 'FontSize', 16)
ax0.YAxis(2).Color = Color;
ax0.YAxis(2).FontSize = FontSize;
ax0.YAxis(2).Limits = [0 3];
ax0.YAxis(2).TickLabelsMode = 'manual';
ax0.YAxis(2).TickValues = [0 0.5 1 1.5 2 2.5 3];
ax0.YAxis(2).TickLabels = {'0.0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'; '3.0'};
ax0.YAxis(2).TickLabelInterpreter = 'latex';
ax0.YAxis(2).Label.FontSize = FontSize; % geen effect
ax0.YAxis(2).Label.FontWeight = 'normal';
ax0.YAxis(2).Label.Color = Color;
ax0.YAxis(2).Label.LineWidth = LineWidth;
ax0.YAxis(2).Label.Interpreter = 'latex';
ax0.YAxis(2).Label.String = '$\mathrm{Ma} \, [-]$';

ax0.XLabel.Color = Color;
ax0.XLabel.FontSize = FontSize;
ax0.XLimMode = 'manual';
ax0.XLim = [0 35];
ax0.XTickLabelMode = 'auto';
ax0.XTick = [0 5 10 15 20 25 30 35];
ax0.XTickLabel = {'0.0'; '5.0'; '10.0'; '15.0'; '20.0'; '25.0'; '30.0'; '35.0'};
ax0.TickLabelInterpreter = 'latex';
ax0.XLabel.FontSize = FontSize; % geen effect
ax0.XLabel.FontWeight = 'normal';
ax0.XLabel.Color = Color;
ax0.XLabel.LineWidth = LineWidth;
ax0.XLabel.Interpreter = 'latex';
ax0.XLabel.String = '$X \, [\mathrm{m}]$';
% ax0.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11) % because fontsize doesn't have any effect on results
box on
hold on

% annotation(h0,'arrow',[0.5 0.8],[0.2 0.8]);
% annotation(h0,'textbox',[0.35 0.30 0.2 0.1],'String',{'Compressibility'},'Interpreter','latex','FontSize',FontSize,'LineStyle','none','FitBoxToText','off','VerticalAlignment','bottom')

% figure 0 - Mach-number profiles
h1 = figure;
ax1 = gca;

ax1.FontSizeMode = 'manual';
ax1.FontSize = FontSize; %%
ax1.FontWeight = 'normal';
ax1.TickLabelInterpreter = 'latex';

ax1.XLabel.Interpreter = 'latex';
ax1.XLabel.FontSize = FontSize; %%
ax1.XLabel.FontWeight = 'normal';
ax1.XLabel.Color = Color;
ax1.XLabel.String = '$Y \, [\mathrm{m}]$';
ax1.XLim = [-0.005 0.18]; % change to 0.15?
ax1.XLimMode = 'manual';
ax1.XTickLabelMode = 'manual';
ax1.XTickMode = 'manual';
ax1.XTick = [0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18];
ax1.XTickLabel = {'0.00'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'; '0.12'; '0.14'; '0.16'; '0.18'};

ax1.YLabel.Interpreter = 'latex';
ax1.YLabel.FontSize = FontSize; %%
ax1.YLabel.FontWeight = 'normal';
ax1.YLabel.Color = Color;
ax1.YLabel.String = '$\mathrm{Ma} \, [-]$';
ax1.YLim = [0 2.9];
ax1.YLimMode = 'manual';
ax1.YTickLabelMode = 'manual';
ax1.YTickMode = 'manual';
ax1.YTick = [0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8];
ax1.YTickLabel = {0; 0.2; 0.4; 0.6; 0.8; 1.0; 1.2; 1.4; 1.6; 1.8; 2.0; 2.2; 2.4; 2.6; 2.8};
ax1.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on
hold on

% figure 2 - Velocity profiles
% h2 = figure(n+1);
h2 = figure;
ax2 = gca;

ax2.FontSizeMode = 'manual';
ax2.FontSize = FontSize; %%
ax2.FontWeight = 'normal';
ax2.TickLabelInterpreter = 'latex';

ax2.XLabel.Interpreter = 'latex';
ax2.XLabel.FontSize = FontSize; %%
ax2.XLabel.FontWeight = 'normal';
ax2.XLabel.Color = Color;
ax2.XLabel.String = '$Y \, [\mathrm{m}]$';
ax2.XLim = [-0.005 0.18]; % change to 0.15?
ax2.XLimMode = 'manual';
ax2.XTickLabelMode = 'manual';
ax2.XTickMode = 'manual';
ax2.XTick = [0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18];
ax2.XTickLabel = {'0.00'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'; '0.12'; '0.14'; '0.16'; '0.18'};

ax2.YLabel.Interpreter = 'latex';
ax2.YLabel.FontSize = FontSize; %%
ax2.YLabel.FontWeight = 'normal';
ax2.YLabel.Color = Color;
ax2.YLabel.String = '$u/u_{\mathrm{e}} \, [-]$';
ax2.YLim = [0.4 1.35];
ax2.YLimMode = 'manual';
ax2.YTickLabelMode = 'manual';
ax2.YTickMode = 'manual';
ax2.YTick = [0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3];
ax2.YTickLabel = {'0.4'; '0.4'; '0.4'; '0.4'; '0.5'; '0.6'; '0.7'; '0.8'; '0.9'; '1.0'};
ax2.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on
hold on

% figure 3 - Density profiles
h3 = figure;
ax3 = gca;

ax3.FontSizeMode = 'manual';
ax3.FontSize = FontSize; %%
ax3.FontWeight = 'normal';
ax3.TickLabelInterpreter = 'latex';

ax3.XLabel.Interpreter = 'latex';
ax3.XLabel.FontSize = FontSize; %%
ax3.XLabel.FontWeight = 'normal';
ax3.XLabel.Color = Color;
ax3.XLabel.String = '$Y \, [\mathrm{m}]$';
ax3.XLim = [-0.005 0.18]; % change to 0.15?
ax3.XLimMode = 'manual';
ax3.XTickLabelMode = 'manual';
ax3.XTickMode = 'manual';
ax3.XTick = [0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18];
ax3.XTickLabel = {'0.00'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'; '0.12'; '0.14'; '0.16'; '0.18'};

ax3.YLabel.Interpreter = 'latex';
ax3.YLabel.FontSize = FontSize; %%
ax3.YLabel.FontWeight = 'normal';
ax3.YLabel.Color = Color;
ax3.YLabel.String = '$c = \rho_{\mathrm{e}} / \rho = T / T_{\mathrm{e}} \, [-]$';
ax3.YLim = [0.95 2.8];
ax3.YLimMode = 'manual';
ax3.YTickLabelMode = 'manual';
ax3.YTickMode = 'manual';
ax3.YTick = [1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7];
ax3.YTickLabel = {'1.0'; '1.0'; '1.0'; '1.0'; '1.1'; '1.2'; '1.3'; '1.4'; '1.5'; '1.6'; '1.7'; '1.8'; '1.9'; '2.0'; '2.1'; '2.2'; '2.3'; '2.4'};
ax3.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on
hold on

% figure 3b - Density profiles
h3b = figure;
ax3b = gca;

ax3b.FontSizeMode = 'manual';
ax3b.FontSize = FontSize; %%
ax3b.FontWeight = 'normal';
ax3b.TickLabelInterpreter = 'latex';

ax3b.XLabel.Interpreter = 'latex';
ax3b.XLabel.FontSize = FontSize; %%
ax3b.XLabel.FontWeight = 'normal';
ax3b.XLabel.Color = Color;
ax3b.XLabel.String = '$Y \, [\mathrm{m}]$';
ax3b.XLim = [0 0.18]; % change to 0.15?
ax3b.XLimMode = 'manual';
ax3b.XTickLabelMode = 'manual';
ax3b.XTickMode = 'manual';
ax3b.XTick = [0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18];
ax3b.XTickLabel = {'0.00'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'; '0.12'; '0.14'; '0.16'; '0.18'};

ax3b.YLabel.Interpreter = 'latex';
ax3b.YLabel.FontSize = FontSize; %%
ax3b.YLabel.FontWeight = 'normal';
ax3b.YLabel.Color = Color;
ax3b.YLabel.String = '$\rho / \rho_{\mathrm{e}} = T_{\mathrm{e}} / T \, [-]$';
ax3b.YLim = [0.25 1.02]; % ax3b.YLim = [0.4 1.02];
ax3b.YLimMode = 'manual';
ax3b.YTickLabelMode = 'manual';
ax3b.YTickMode = 'manual';
ax3b.YTick = [0.25 0.35 0.45 0.55 0.65 0.75 0.85 0.9 0.95 1.0];
ax3b.YTickLabel = {'0.4'; '0.5'; '0.6'; '0.7'; '0.8'; '0.9'; '1.0'; '1.0'; '1.0'; '1.0'};
ax3b.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on
hold on

% figure 3c - Density profiles - Log Plot
h3c = figure;
ax3c = gca;

ax3c.FontSizeMode = 'manual';
ax3c.FontSize = FontSize; %%
ax3c.FontWeight = 'normal';
ax3c.TickLabelInterpreter = 'latex';
ax3c.XScale = 'linear';
ax3c.YScale = 'log';

ax3c.XLabel.Interpreter = 'latex';
ax3c.XLabel.FontSize = FontSize; %%
ax3c.XLabel.FontWeight = 'normal';
ax3c.XLabel.Color = Color;
ax3c.XLabel.String = '$Y \, [\mathrm{m}]$';
ax3c.XLim = [0 0.18]; % change to 0.15?
ax3c.XLimMode = 'manual';
ax3c.XTickLabelMode = 'manual';
ax3c.XTickMode = 'manual';
ax3c.XTick = [0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18];
ax3c.XTickLabel = {'0.00'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'; '0.12'; '0.14'; '0.16'; '0.18'};

ax3c.YLabel.Interpreter = 'latex';
ax3c.YLabel.FontSize = FontSize; %%
ax3c.YLabel.FontWeight = 'normal';
ax3c.YLabel.Color = Color;
ax3c.YLabel.String = '$c = \rho_{\mathrm{e}} / \rho = T / T_{\mathrm{e}} \, [-]$';
% ax3c.YLim = [0.90 2.8];
% ax3c.YLimMode = 'manual';
% ax3c.YTickLabelMode = 'manual';
% ax3c.YTickMode = 'manual';
% ax3c.YTick = [1.0 1.1 1.2 1.3 1.5 1.7 1.9 2.1 2.3 2.5 2.7 2.9];
% ax3c.YTickLabel = {'1.0'; '1.0'; '1.0'; '1.0'; '1.2'; '1.4'; '1.6'; '1.8'; '2.0'; '2.2'; '2.4'; '2.6'};
ax3c.YLimMode = 'manual';
ax3c.YTickLabelMode = 'manual';
ax3c.YTickMode = 'manual';
ax3c.YLim = [10^-5 10^0.2];
ax3c.YTick = [10^-5 10^-4 10^-3 10^-2 10^-1.5 10^-1 10^-0.5 10^0 10^0.1 10^0.2];
ax3c.YTickLabel = {'$10^{-5}$'; '$10^{-4}$'; '$10^{-3}$'; '$10^{-2}$'; '$10^{-1.5}$'; '$10^{-1}$'; '$10^{-0.5}$'; '$10^{0}$'; '$10^{0.1}$'; '$10^{0.2}$'};
ax3c.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on
hold on

% figure 4a - log-law velocity profiles (normal, with wall density estimate
% of Fernholz & Finley (1970)
h4a = figure;
ax4a = gca;

ax4a.FontSizeMode = 'manual';
ax4a.FontSize = FontSize; %%
ax4a.FontWeight = 'normal';
ax4a.TickLabelInterpreter = 'latex';
ax4a.XScale = 'log';
ax4a.YScale = 'linear';

ax4a.XLabel.Interpreter = 'latex';
ax4a.XLabel.FontSize = FontSize; %%
ax4a.XLabel.FontWeight = 'normal';
ax4a.XLabel.Color = Color;
ax4a.XLabel.String = '$y^+_{\mathrm{n}} \, [-]$';
ax4a.XLabel.String = '$y^{\ast} \, [-]$';
ax4a.XLim = [10^-1 10^5];
ax4a.XLimMode = 'manual';
ax4a.XTickLabelMode = 'manual';
ax4a.XTickMode = 'manual';
ax4a.XTick = [10^-1 10^0 10^1 10^2 10^3 10^4 10^5];
ax4a.XTickLabel = {'$10^{-1}$'; '$10^0$'; '$10^1$'; '$10^2$'; '$10^3$'; '$10^4$'; '$10^5$'};

ax4a.YLabel.Interpreter = 'latex';
ax4a.YLabel.FontSize = FontSize; %%
ax4a.YLabel.FontWeight = 'normal';
ax4a.YLabel.Color = Color;
ax4a.YLabel.String = '$u^+_{\mathrm{n}} \, [-]$';
ax4a.YLabel.String = '$u^{\ast} \, [-]$';
% ax4a.YLim = [0 35];
% ax4a.YLimMode = 'manual';
% ax4a.YTickLabelMode = 'manual';
% ax4a.YTickMode = 'manual';
% ax4a.YTick = [0 5 10 15 20 25 30 35];
% ax4a.YTickLabel = {'0'; '5'; '10'; '15'; '20'; '25'; '30'; '35'};
ax4a.YLim = [0 58];
ax4a.YLimMode = 'manual';
ax4a.YTickLabelMode = 'manual';
ax4a.YTickMode = 'manual';
ax4a.YTick = [0 5 10 15 20 25 30 35 40 45 50 55];
ax4a.YTickLabel = {'0'; '0'; '0'; '0'; '5'; '10'; '15'; '20'; '25'; '30'; '35'; '40'};
ax4a.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on
hold on

% figure 4b - log-law velocity profiles (not using any estimate of wall
% density, only measured data; there is a possibility to estimate 
% compressibility factor Fc to correct for compressibility)
h4b = figure;
ax4b = gca;

ax4b.FontSizeMode = 'manual';
ax4b.FontSize = FontSize; %%
ax4b.FontWeight = 'normal';
ax4b.TickLabelInterpreter = 'latex';
ax4b.XScale = 'log';
ax4b.YScale = 'linear';

ax4b.XLabel.Interpreter = 'latex';
ax4b.XLabel.FontSize = FontSize; %%
ax4b.XLabel.FontWeight = 'normal';
ax4b.XLabel.Color = Color;
ax4b.XLabel.String = '$y^+ \, [-]$';
ax4b.XLim = [10^-1 10^5]; % change to 0.15?
ax4b.XLimMode = 'manual';
ax4b.XTickLabelMode = 'manual';
ax4b.XTickMode = 'manual';
ax4b.XTick = [10^-1 10^0 10^1 10^2 10^3 10^4 10^5];
ax4b.XTickLabel = {'$10^{-1}$'; '$10^0$'; '$10^1$'; '$10^2$'; '$10^3$'; '$10^4$'; '$10^5$'};

ax4b.YLabel.Interpreter = 'latex';
ax4b.YLabel.FontSize = FontSize; %%
ax4b.YLabel.FontWeight = 'normal';
ax4b.YLabel.Color = Color;
ax4b.YLabel.String = '$u^+ \, [-]$';
% ax4b.YLim = [0 35];
% ax4b.YLimMode = 'manual';
% ax4b.YTickLabelMode = 'manual';
% ax4b.YTickMode = 'manual';
% ax4b.YTick = [0 5 10 15 20 25 30 35];
% ax4b.YTickLabel = {'0'; '5'; '10'; '15'; '20'; '25'; '30'; '35'};
ax4b.YLim = [0 43];
ax4b.YLimMode = 'manual';
ax4b.YTickLabelMode = 'manual';
ax4b.YTickMode = 'manual';
ax4b.YTick = [0 5 10 15 20 25 30 35 40];
ax4b.YTickLabel = {'0'; '0'; '0'; '0'; '5'; '10'; '15'; '20'; '25'};
ax4b.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on
hold on

% figure 5 - log-law density (inverse temperature) profiles (normal, with wall density estimate
% of Fernholz & Finley (1970)
h5 = figure;
ax5 = gca;

ax5.FontSizeMode = 'manual';
ax5.FontSize = FontSize; %%
ax5.FontWeight = 'normal';
ax5.TickLabelInterpreter = 'latex';
ax5.XScale = 'log';
ax5.YScale = 'linear';

ax5.XLabel.Interpreter = 'latex';
ax5.XLabel.FontSize = FontSize; %%
ax5.XLabel.FontWeight = 'normal';
ax5.XLabel.Color = Color;
ax5.XLabel.String = '$y^+ \, [-]$';
ax5.XLim = [10^-1 10^5]; % change to 0.15?
ax5.XLimMode = 'manual';
ax5.XTickLabelMode = 'manual';
ax5.XTickMode = 'manual';
ax5.XTick = [10^-1 10^0 10^1 10^2 10^3 10^4 10^5];
ax5.XTickLabel = {'$10^{-1}$'; '$10^0$'; '$10^1$'; '$10^2$'; '$10^3$'; '$10^4$'; '$10^5$'};

ax5.YLabel.Interpreter = 'latex';
ax5.YLabel.FontSize = FontSize; %%
ax5.YLabel.FontWeight = 'normal';
ax5.YLabel.Color = Color;
ax5.YLabel.String = '$\rho^* = \frac{\rho - \rho_{\mathrm{w}}}{\rho_{\mathrm{e}} - \rho_{\mathrm{w}}} \, [-]$';
ax5.YLim = [0 1.65];
ax5.YLimMode = 'manual';
ax5.YTickLabelMode = 'manual';
ax5.YTickMode = 'manual';
% ax5.YTick = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3];
% ax5.YTickLabel = {'0'; '0'; '0'; '0'; '0.1'; '0.2'; '0.3'; '0.4'; '0.5'; '0.6'; '0.7'; '0.8'; '0.9'; '1.0'};
ax5.YTick = [0 0.2 0.4 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6];
ax5.YTickLabel = {'0'; '0'; '0'; '0'; '0.1'; '0.2'; '0.3'; '0.4'; '0.5'; '0.6'; '0.7'; '0.8'; '0.9'; '1.0'};
ax5.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on
hold on

% figure 6 - Nozzle Geometry
% h6 = figure;
% L=length(Geom(:,1));
% plot3(Geom(:,1),Geom(:,2),Geom(:,3),'Color',Color,'Linestyle','-')
% hold on
% axis 'equal'
% plot3(Geom(:,1),-Geom(:,2),Geom(:,3),'Color',Color,'Linestyle','-')
% plot3(Geom(:,1),-Geom(:,2),2*Geom(end,2)*ones(L,1),'Color',Color,'Linestyle','-')
% plot3(Geom(:,1),Geom(:,2),2*Geom(end,2)*ones(L,1),'Color',Color,'Linestyle','-')
% plot3([Geom(1,1) Geom(1,1)],[-Geom(1,2) Geom(1,2)],[0 0],'Color',Color,'Linestyle','-')
% plot3([Geom(end,1) Geom(end,1)],[-Geom(end,2) Geom(end,2)],[0 0],'Color',Color,'Linestyle','-')
% plot3([Geom(1,1) Geom(1,1)],[-Geom(1,2) Geom(1,2)],[2*Geom(end,2)*ones(L,1) 2*Geom(end,2)*ones(L,1)],'Color',Color,'Linestyle','-')
% plot3([Geom(end,1) Geom(end,1)],[-Geom(end,2) Geom(end,2)],[2*Geom(end,2)*ones(L,1) 2*Geom(end,2)*ones(L,1)],'Color',Color,'Linestyle','-')
% plot3([Geom(1,1) Geom(1,1)],[-Geom(1,2) -Geom(1,2)],[0 2*Geom(end,2)],'Color',Color,'Linestyle','-')
% plot3([Geom(1,1) Geom(1,1)],[Geom(1,2) Geom(1,2)],[0 2*Geom(end,2)],'Color',Color,'Linestyle','-')
% plot3([Geom(end,1) Geom(end,1)],[-Geom(end,2) -Geom(end,2)],[0 2*Geom(end,2)],'Color',Color,'Linestyle','-')
% plot3([Geom(end,1) Geom(end,1)],[Geom(end,2) Geom(end,2)],[0 2*Geom(end,2)],'Color',Color,'Linestyle','-')
% plot3([Geom(1,1) Geom(end,1)],[0 0],[2*Geom(end,2) 2*Geom(end,2)],'Color',Color,'LineStyle','-.')
% plot3([Geom(end,1) Geom(end,1)],[0 0],[2*Geom(end,2) 0.9*Geom(end,2)],'Color',Color,'Linestyle','-')
% plot3(Geom(end,1),0,0.9*Geom(end,2),'Color',Color,'Marker','v')

% % % % added (18-03-2020)
% % % Geom(:,1) = Geom(:,1) - Geom(1,1);
% % % % added (18-03-2020)

% figure 6 - Nozzle Geometry
h6 = figure;
ax6 = gca;
L1 = Geom(1,1); %%% added (19-03-2020)
Geom(:,1) = Geom(:,1) - L1; %%% added (19-03-2020)
L=length(Geom(:,1));
plot3(Geom(:,1),Geom(:,3),Geom(:,2),'Color',Color,'Linestyle','-')
hold on
axis 'equal'
plot3(Geom(:,1),Geom(:,3),-Geom(:,2),'Color',Color,'Linestyle','-')
plot3(Geom(:,1),2*Geom(end,2)*ones(L,1),-Geom(:,2),'Color',Color,'Linestyle','-')
plot3(Geom(:,1),2*Geom(end,2)*ones(L,1),Geom(:,2),'Color',Color,'Linestyle','-')
plot3([Geom(1,1) Geom(1,1)],[0 0],[-Geom(1,2) Geom(1,2)],'Color',Color,'Linestyle','-')
plot3([Geom(end,1) Geom(end,1)],[0 0],[-Geom(end,2) Geom(end,2)],'Color',Color,'Linestyle','-')
plot3([Geom(1,1) Geom(1,1)],[2*Geom(end,2)*ones(L,1) 2*Geom(end,2)*ones(L,1)],[-Geom(1,2) Geom(1,2)],'Color',Color,'Linestyle','-')
plot3([Geom(end,1) Geom(end,1)],[2*Geom(end,2)*ones(L,1) 2*Geom(end,2)*ones(L,1)],[-Geom(end,2) Geom(end,2)],'Color',Color,'Linestyle','-')
plot3([Geom(1,1) Geom(1,1)],[0 2*Geom(end,2)],[-Geom(1,2) -Geom(1,2)],'Color',Color,'Linestyle','-')
plot3([Geom(1,1) Geom(1,1)],[0 2*Geom(end,2)],[Geom(1,2) Geom(1,2)],'Color',Color,'Linestyle','-')
plot3([Geom(end,1) Geom(end,1)],[0 2*Geom(end,2)],[-Geom(end,2) -Geom(end,2)],'Color',Color,'Linestyle','-')
plot3([Geom(end,1) Geom(end,1)],[0 2*Geom(end,2)],[Geom(end,2) Geom(end,2)],'Color',Color,'Linestyle','-')
plot3([Geom(1,1) Geom(end,1)],[2*Geom(end,2) 2*Geom(end,2)],[0 0],'Color',Color,'LineStyle','-.')
plot3([Geom(end,1) Geom(end,1)],[2*Geom(end,2) 0.9*Geom(end,2)],[0 0],'Color',Color,'Linestyle','-')
plot3(Geom(end,1),2*Geom(end,2),0,'Color',Color,'Marker','x')
plot3(Geom(end,1),0.9*Geom(end,2),0,'Color',Color,'Marker','>')
% COLOR SURFACES
C(:,:,1) = 0.9400.*ones(L); % red
C(:,:,2) = 0.9400.*ones(L); % green
C(:,:,3) = 0.9400.*ones(L); % blue
% Bottom surface
Xb = Geom(:,1)'.*ones(L,L);
Yb = linspace(Geom(end,3),2*Geom(end,2),L)'.*ones(L,L);
Zb = -Geom(:,2)'.*ones(L,L);
sb = surf(Xb,Yb,Zb,C);
sb.EdgeColor = 'none';
% Top surface
Xt = Geom(:,1)'.*ones(L,L);
Yt = linspace(Geom(end,3),2*Geom(end,2),L)'.*ones(L,L);
Zt = Geom(:,2)'.*ones(L,L);
st = surf(Xt,Yt,Zt,C);
st.EdgeColor = 'none';

% Camera view 1 %%%%%%%%%%%%%%%%%%%%%%%%%%
% ax6.CameraPositionMode = 'manual';
% ax6.CameraTargetMode = 'manual';
% ax6.CameraUpVectorMode = 'manual';
% ax6.CameraViewAngleMode = 'manual';
% ax6.CameraPosition = [-21.1352  -93.2874   -3.6962];
% ax6.CameraTarget = [26.8195 1.0714 0];
% ax6.CameraUpVector = [0 0 1];
% ax6.CameraViewAngle = 10.3018;
% Camera view 2 %%%%%%%%%%%%%%%%%%%%%%%%%%
% ax6.CameraPositionMode = 'manual';
% ax6.CameraTargetMode = 'manual';
% ax6.CameraUpVectorMode = 'manual';
% ax6.CameraViewAngleMode = 'manual';
% ax6.CameraPosition = [-50.4727 - L1  -69.2996   17.0509]; % only this one changes!
% ax6.CameraTarget = [26.8195 - L1 1.0714 0];
% % Xnew = 26.8195 - Geom(1,1);
% ax6.CameraUpVector = [0 0 1];
% ax6.CameraViewAngle = 10.3018;

ax6.TickLabelInterpreter = 'latex';

% added (18-03-2020)
% ax6.XLim = [0 20];
ax6.XLimMode = 'manual';
ax6.XTickLabelMode = 'manual';
ax6.XTickMode = 'manual';
ax6.XLim = [0 20.33];
ax6.XTick = [0 5 10 15 20];
ax6.XTickLabel = {'0'; '5'; '10'; '15'; '20'};

ax6.ZLim = [-2.8 2.8];
ax6.ZLimMode = 'manual';
ax6.ZTickLabelMode = 'manual';
ax6.ZTickMode = 'manual';
ax6.ZTick = [-2.8 0 2.8];
ax6.ZTickLabel = {'-2.8'; '0'; '2.8'};

% ax6.ZLim = [0 0.1];
% ax6.ZLimMode = 'manual';
ax6.YTickLabelMode = 'manual';
ax6.YTickMode = 'manual';
ax6.YTick = [];
% ax6.ZTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};
% added (18-03-2020)
% ax6.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

% already the case?
% ax6.DataAspectRatioMode = 'manual';
% ax6.DataAspectRatio = [1 1 1]; % axis = 'equal';

% ax6.XLabel.Interpreter = 'latex';
% ax6.XLabel.FontSize = FontSize; %%
% ax6.XLabel.FontWeight = 'normal';
% ax6.XLabel.Color = Color;
% ax6.XLabel.String = 'X [m]';
% 
% ax6.YLabel.Interpreter = 'latex';
% ax6.YLabel.FontSize = FontSize; %%
% ax6.YLabel.FontWeight = 'normal';
% ax6.YLabel.Color = Color;
% ax6.YLabel.String = 'Y [m]';
% 
% ax6.ZLabel.Interpreter = 'latex';
% ax6.ZLabel.FontSize = FontSize; %%
% ax6.ZLabel.FontWeight = 'normal';
% ax6.ZLabel.Color = Color;
% ax6.ZLabel.String = 'Z [m]';

ax6.CameraPositionMode = 'manual';
ax6.CameraTargetMode = 'manual';
ax6.CameraUpVectorMode = 'manual';
ax6.CameraViewAngleMode = 'manual';
ax6.CameraUpVector = [0 0 1];
ax6.CameraViewAngle = 10.3018;
ax6.CameraPosition = [-50.4727 - L1  -69.2996   17.0509]; % only this one changes!
ax6.CameraTarget = [26.8195 - L1 1.0714 0];

% saveas(h6,'Thesis_Report_Figures/NozzleGeometry','epsc')
% keyboard

% % figure 6 - Dimensionless Velocity profiles
% % h2 = figure(n+1);
% h6 = figure;
% ax6 = gca;
% 
% ax6.FontSizeMode = 'manual';
% ax6.FontSize = FontSize; %%
% ax6.FontWeight = 'normal';
% ax6.TickLabelInterpreter = 'latex';
% 
% ax6.XLabel.Interpreter = 'latex';
% ax6.XLabel.FontSize = FontSize; %%
% ax6.XLabel.FontWeight = 'normal';
% ax6.XLabel.Color = Color;
% ax6.XLabel.String = '$Y \, [\mathrm{m}]$';%'$y/\delta \, [-]$';
% ax6.XLim = [0 0.18]; % change to 0.15?
% ax6.XLimMode = 'manual';
% ax6.XTickLabelMode = 'manual';
% ax6.XTickMode = 'manual';
% ax6.XTick = [0 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18];
% ax6.XTickLabel = {'0.00'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'; '0.12'; '0.14'; '0.16'; '0.18'};
% 
% ax6.YLabel.Interpreter = 'latex';
% ax6.YLabel.FontSize = FontSize; %%
% ax6.YLabel.FontWeight = 'normal';
% ax6.YLabel.Color = Color;
% ax6.YLabel.String = '$u/u_{\mathrm{e}} \, [-]$';
% ax6.YLim = [0.4 1.2];
% ax6.YLimMode = 'manual';
% ax6.YTickLabelMode = 'manual';
% ax6.YTickMode = 'manual';
% ax6.YTick = [0.4 0.6 0.8 1.0 1.2];
% ax6.YTickLabel = {'0.4'; '0.6'; '0.8'; '1.0'; '1.2'};
% ax6.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
% 
% set(gca,'FontSize',11)
% box on
% hold on

% % figure 7 - Boundary Layer Thicknesses
% h7 = figure;
% ax7 = gca;
% 
% ax7.FontSizeMode = 'manual';
% ax7.FontSize = FontSize; %%
% ax7.FontWeight = 'normal';
% ax7.TickLabelInterpreter = 'latex';
% 
% ax7.XLabel.Interpreter = 'latex';
% ax7.XLabel.FontSize = FontSize; %%
% ax7.XLabel.FontWeight = 'normal';
% ax7.XLabel.Color = Color;
% ax7.XLabel.String = '$\mathrm{Ma} \, [-]$';
% ax7.XLim = [0 3];
% ax7.XLimMode = 'manual';
% ax7.XTickLabelMode = 'manual';
% ax7.XTickMode = 'manual';
% ax7.XTick = [0 0.5 1.0 1.5 2.0 2.5 3.0];
% ax7.XTickLabel = {'0.0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'; '3.0'};
% 
% ax7.YLabel.Interpreter = 'latex';
% ax7.YLabel.FontSize = FontSize; %%
% ax7.YLabel.FontWeight = 'normal';
% ax7.YLabel.Color = Color;
% ax7.YLabel.String = '$\delta^{\ast}, \, \theta \, \times \, 10^{-3} \, [\mathrm{m}]$';
% ax7.YLim = [0 30];
% ax7.YLimMode = 'manual';
% ax7.YTickLabelMode = 'manual';
% ax7.YTickMode = 'manual';
% ax7.YTick = [0 5 10 15 20 25 30];
% ax7.YTickLabel = {'0'; '5'; '10'; '15'; '20'; '25'; '30'};
% ax7.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
% 
% set(gca,'FontSize',11)
% box on
% hold on

% figure 7 - Boundary Layer Thicknesses
h7 = figure;
ax7 = gca;

ax7.FontSizeMode = 'manual';
ax7.FontSize = FontSize; %%
ax7.FontWeight = 'normal';
ax7.TickLabelInterpreter = 'latex';

ax7.XLabel.Interpreter = 'latex';
ax7.XLabel.FontSize = FontSize; %%
ax7.XLabel.FontWeight = 'normal';
ax7.XLabel.Color = Color;
ax7.XLabel.String = '$\mathrm{Ma} \, [-]$';
ax7.XLim = [0 3];
ax7.XLimMode = 'manual';
ax7.XTickLabelMode = 'manual';
ax7.XTickMode = 'manual';
ax7.XTick = [0 0.5 1.0 1.5 2.0 2.5 3.0];
ax7.XTickLabel = {'0.0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'; '3.0'};

yyaxis left
ax7.YAxis(1).Color = Color;
ax7.YAxis(1).FontSize = FontSize;
ax7.YAxis(1).Limits = [0 30];
ax7.YAxis(1).TickLabelsMode = 'manual';
ax7.YAxis(1).TickValues = [0 5 10 15 20 25 30];
ax7.YAxis(1).TickLabels = {'0'; '5'; '10'; '15'; '20'; '25'; '30'};
ax7.YAxis(1).TickLabelInterpreter = 'latex';
ax7.YAxis(1).Label.FontSize = FontSize; % geen effect
ax7.YAxis(1).Label.FontWeight = 'normal';
ax7.YAxis(1).Label.Color = Color;
ax7.YAxis(1).Label.LineWidth = LineWidth;
ax7.YAxis(1).Label.Interpreter = 'latex';
ax7.YAxis(1).Label.String = '$\delta^{\ast}, \, \theta \, \times \, 10^{-3} \, [\mathrm{m}]$';

yyaxis right
ax7.YAxis(2).Color = Color;
ax7.YAxis(2).FontSize = FontSize;
ax7.YAxis(2).Limits = [0 30];
ax7.YAxis(2).TickLabelsMode = 'manual';
ax7.YAxis(2).TickValues = [0 5 10 15 20 25 30];
ax7.YAxis(2).TickLabels = {};
ax7.YAxis(2).TickLabelInterpreter = 'latex';
ax7.YAxis(2).Label.FontSize = FontSize; % geen effect
ax7.YAxis(2).Label.FontWeight = 'normal';
ax7.YAxis(2).Label.Color = Color;
ax7.YAxis(2).Label.LineWidth = LineWidth;
ax7.YAxis(2).Label.Interpreter = 'latex';
ax7.YAxis(2).Label.String = '';
ax7.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on
hold on

% % figure X (two y-axes) - cf-Cd plot (two y-axes)
% h8 = figure;
% 
% ax8 = gca;
% 
% yyaxis left
% ax8.YAxis(1).Color = Color;
% ax8.YAxis(1).FontSize = FontSize;
% ax8.YAxis(1).Limits = [1 2];
% ax8.YAxis(1).TickLabelsMode = 'manual';
% ax8.YAxis(1).TickValues = [1 1.2 1.4 1.6 1.8 2];
% ax8.YAxis(1).TickLabels = {'1.0'; '1.2'; '1.4'; '1.6'; '1.8'; '2.0'};
% ax8.YAxis(1).TickLabelInterpreter = 'latex';
% ax8.YAxis(1).Label.FontSize = FontSize; % geen effect
% ax8.YAxis(1).Label.FontWeight = 'normal';
% ax8.YAxis(1).Label.Color = Color;
% ax8.YAxis(1).Label.LineWidth = LineWidth;
% ax8.YAxis(1).Label.Interpreter = 'latex';
% ax8.YAxis(1).Label.String = '$c_{\mathrm{f}} \, \times \, 10^{-3} \, [-]$';
% 
% yyaxis right
% ax8.YAxis(2).Color = Color;
% ax8.YAxis(2).FontSize = FontSize;
% ax8.YAxis(2).Limits = [0 1];
% ax8.YAxis(2).TickLabelsMode = 'manual';
% ax8.YAxis(2).TickValues = [0 0.2 0.4 0.6 0.8 1.0];
% ax8.YAxis(2).TickLabels = {'0.0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'};
% ax8.YAxis(2).TickLabelInterpreter = 'latex';
% ax8.YAxis(2).Label.FontSize = FontSize; % geen effect
% ax8.YAxis(2).Label.FontWeight = 'normal';
% ax8.YAxis(2).Label.Color = Color;
% ax8.YAxis(2).Label.LineWidth = LineWidth;
% ax8.YAxis(2).Label.Interpreter = 'latex';
% ax8.YAxis(2).Label.String = '$C_{\mathrm{d}} \, \times \, 10^{-3} \, [-]$';
% 
% ax8.XLabel.Color = Color;
% ax8.XLabel.FontSize = FontSize;
% ax8.XLimMode = 'manual';
% ax8.XLim = [0 3];
% ax8.XTickLabelMode = 'auto';
% ax8.XTick = [0 0.5 1.0 1.5 2.0 2.5 3.0];
% ax8.XTickLabel = {'0.0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'; '3.0'};
% ax8.TickLabelInterpreter = 'latex';
% ax8.XLabel.FontSize = FontSize; % geen effect
% ax8.XLabel.FontWeight = 'normal';
% ax8.XLabel.Color = Color;
% ax8.XLabel.LineWidth = LineWidth;
% ax8.XLabel.Interpreter = 'latex';
% ax8.XLabel.String = '$\mathrm{Ma} \, [-]$';
% ax8.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
% 
% set(gca,'FontSize',11) % because fontsize doesn't have any effect on results
% box on
% hold on

% figure X (single y-axis) - cf-Cd plot
h8 = figure;
ax8 = gca;

ax8.FontSizeMode = 'manual';
ax8.FontSize = FontSize; %%
ax8.FontWeight = 'normal';
ax8.TickLabelInterpreter = 'latex';

ax8.XLabel.Interpreter = 'latex';
ax8.XLabel.FontSize = FontSize; %%
ax8.XLabel.FontWeight = 'normal';
ax8.XLabel.Color = Color;
ax8.XLabel.String = '$\mathrm{Ma} \, [-]$';
ax8.XLim = [0 3];
ax8.XLimMode = 'manual';
ax8.XTickLabelMode = 'manual';
ax8.XTickMode = 'manual';
ax8.XTick = [0 0.5 1.0 1.5 2.0 2.5 3.0];
ax8.XTickLabel = {'0.0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'; '3.0'};

ax8.YLabel.Interpreter = 'latex';
ax8.YLabel.FontSize = FontSize; %%
ax8.YLabel.FontWeight = 'normal';
ax8.YLabel.Color = Color;
% ax8.YLabel.String = '$c_{\mathrm{f}}, \, C_{\mathrm{d}} \, \times \, 10^{-3} \, [-]$';
ax8.YLabel.String = '$c_{\mathrm{f}} \, \times \, 10^{-3} \, [-]$';
ax8.YLim = [0 2];
ax8.YLimMode = 'manual';
ax8.YTickLabelMode = 'manual';
ax8.YTickMode = 'manual';
ax8.YTick = [0 0.5 1.0 1.5 2.0];
ax8.YTickLabel = {'0.0'; '0.5'; '1.0'; '1.5'; '2.0'};
ax8.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on
hold on

% % figure 9 - Wall properties (single axes)
% h9 = figure;
% ax9 = gca;
% 
% ax9.FontSizeMode = 'manual';
% ax9.FontSize = FontSize; %%
% ax9.FontWeight = 'normal';
% ax9.TickLabelInterpreter = 'latex';
% 
% ax9.XLabel.Interpreter = 'latex';
% ax9.XLabel.FontSize = FontSize; %%
% ax9.XLabel.FontWeight = 'normal';
% ax9.XLabel.Color = Color;
% ax9.XLabel.String = '$\mathrm{Ma} \, [-]$';
% ax9.XLim = [0 3]; % change to 0.15?
% ax9.XLimMode = 'manual';
% ax9.XTickLabelMode = 'manual';
% ax9.XTickMode = 'manual';
% ax9.XTick = [0 0.5 1.0 1.5 2.0 2.5 3.0];
% ax9.XTickLabel = {'0.0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'; '3.0'};
% 
% ax9.YLabel.Interpreter = 'latex';
% ax9.YLabel.FontSize = FontSize; %%
% ax9.YLabel.FontWeight = 'normal';
% ax9.YLabel.Color = Color;
% ax9.YLabel.String = '$c_{\mathrm{w}}, \, C_{\mathrm{w}} \, [-]$';
% ax9.YLim = [0 2.5];
% ax9.YLimMode = 'manual';
% ax9.YTickLabelMode = 'manual';
% ax9.YTickMode = 'manual';
% ax9.YTick = [0 0.5 1.0 1.5 2.0 2.5];
% ax9.YTickLabel = {'0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'};
% ax9.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
% 
% set(gca,'FontSize',11)
% box on
% hold on

% figure 9 - Wall properties (double axes)
h9 = figure;

ax9 = gca;

yyaxis left
ax9.YAxis(1).Color = Color;
ax9.YAxis(1).FontSize = FontSize;
ax9.YAxis(1).Limits = [0.8 2.6];
ax9.YAxis(1).TickLabelsMode = 'manual';
ax9.YAxis(1).TickValues = [0.9 1.0 1.5 2.0 2.5];
ax9.YAxis(1).TickLabels = {'0.9'; '1.0'; '1.5'; '2.0'; '2.5'};
ax9.YAxis(1).TickLabelInterpreter = 'latex';
ax9.YAxis(1).Label.FontSize = FontSize; % geen effect
ax9.YAxis(1).Label.FontWeight = 'normal';
ax9.YAxis(1).Label.Color = Color;
ax9.YAxis(1).Label.LineWidth = LineWidth;
ax9.YAxis(1).Label.Interpreter = 'latex';
% ax9.YAxis(1).Label.String = '$c_{\mathrm{w}} \, [-]$';
ax9.YAxis(1).Label.String = '$c_{\mathrm{w}}, \, C_{\mathrm{w}} \, [-]$';

yyaxis right
ax9.YAxis(2).Color = Color;
ax9.YAxis(2).FontSize = FontSize;
% ax9.YAxis(2).Limits = [0.88 1.02];
% ax9.YAxis(2).TickLabelsMode = 'manual';
% ax9.YAxis(2).TickValues = [0.90 0.92 0.94 0.96 0.98 1.00];
% ax9.YAxis(2).TickLabels = {'0.90'; '0.92'; '0.94'; '0.96'; '0.98'; '1.00'};
ax9.YAxis(2).Limits = [0.8 2.6];
ax9.YAxis(2).TickLabelsMode = 'manual';
ax9.YAxis(2).TickValues = [0.9 1.0 1.5 2.0 2.5];
ax9.YAxis(2).TickLabels = {};
ax9.YAxis(2).TickLabelInterpreter = 'latex';
ax9.YAxis(2).Label.FontSize = FontSize; % geen effect
ax9.YAxis(2).Label.FontWeight = 'normal';
ax9.YAxis(2).Label.Color = Color;
ax9.YAxis(2).Label.LineWidth = LineWidth;
ax9.YAxis(2).Label.Interpreter = 'latex';
% ax9.YAxis(2).Label.String = '$C_{\mathrm{w}} \, [-]$';
ax9.YAxis(2).Label.String = '';

ax9.XLabel.Color = Color;
ax9.XLabel.FontSize = FontSize;
ax9.XLimMode = 'manual';
ax9.XLim = [0 3];
ax9.XTickLabelMode = 'auto';
ax9.XTick = [0 0.5 1.0 1.5 2.0 2.5 3.0];
ax9.XTickLabel = {'0.0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'; '3.0'};
ax9.TickLabelInterpreter = 'latex';
ax9.XLabel.FontSize = FontSize; % geen effect
ax9.XLabel.FontWeight = 'normal';
ax9.XLabel.Color = Color;
ax9.XLabel.LineWidth = LineWidth;
ax9.XLabel.Interpreter = 'latex';
ax9.XLabel.String = '$\mathrm{Ma} \, [-]$';
ax9.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11) % because fontsize doesn't have any effect on results
box on
hold on

% figure 10 - Shape factor
h10 = figure;
ax10 = gca;

ax10.FontSizeMode = 'manual';
ax10.FontSize = FontSize; %%
ax10.FontWeight = 'normal';
ax10.TickLabelInterpreter = 'latex';

ax10.XLabel.Interpreter = 'latex';
ax10.XLabel.FontSize = FontSize; %%
ax10.XLabel.FontWeight = 'normal';
ax10.XLabel.Color = Color;
ax10.XLabel.String = '$\mathrm{Ma} \, [-]$';
ax10.XLim = [0 3];
ax10.XLimMode = 'manual';
ax10.XTickLabelMode = 'manual';
ax10.XTickMode = 'manual';
ax10.XTick = [0 0.5 1.0 1.5 2.0 2.5 3.0];
ax10.XTickLabel = {'0.0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'; '3.0'};

ax10.YLabel.Interpreter = 'latex';
ax10.YLabel.FontSize = FontSize; %%
ax10.YLabel.FontWeight = 'normal';
ax10.YLabel.Color = Color;
ax10.YLabel.String = '$H \, [-]$';
ax10.YLim = [0 5];
ax10.YLimMode = 'manual';
ax10.YTickLabelMode = 'manual';
ax10.YTickMode = 'manual';
ax10.YTick = [0 1.0 2.0 3.0 4.0 5.0];
ax10.YTickLabel = {'0'; '1.0'; '2.0'; '3.0'; '4.0'; '5.0'};
ax10.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on
hold on

% Figure 11 - Temperature recovery factor
h11 = figure;
ax11 = gca;

ax11.FontSizeMode = 'manual';
ax11.FontSize = FontSize; %%
ax11.FontWeight = 'normal';
ax11.TickLabelInterpreter = 'latex';

ax11.XLabel.Interpreter = 'latex';
ax11.XLabel.FontSize = FontSize; %%
ax11.XLabel.FontWeight = 'normal';
ax11.XLabel.Color = Color;
ax11.XLabel.String = '$\mathrm{Ma} \, [-]$';
ax11.XLim = [0 3];
ax11.XLimMode = 'manual';
ax11.XTickLabelMode = 'manual';
ax11.XTickMode = 'manual';
ax11.XTick = [0 0.5 1.0 1.5 2.0 2.5 3.0];
ax11.XTickLabel = {'0.0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'; '3.0'};

ax11.YLabel.Interpreter = 'latex';
ax11.YLabel.FontSize = FontSize; %%
ax11.YLabel.FontWeight = 'normal';
ax11.YLabel.Color = Color;
ax11.YLabel.String = '$r \, [-]$';
ax11.YLim = [0.8 1.0];
ax11.YLimMode = 'manual';
ax11.YTickLabelMode = 'manual';
ax11.YTickMode = 'manual';
ax11.YTick = [0.80 0.85 0.90 0.95 1.0];
ax11.YTickLabel = {'0.80'; '0.85'; '0.90'; '0.95'; '1.0'};
ax11.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on
hold on

% Figure 12 - cf and Cd
h12 = figure();
ax12 = copyobj(ax8,h12); % the fast route
ax12.YLabel.String = '$c_{\mathrm{f}}, \, C_{\mathrm{d}} \, \times \, 10^{-3} \, [-]$';

% Figure 13 - f''(0)
h13 = figure();
ax13 = copyobj(ax8,h13); % the fast route
ax13.YLabel.String = '$f''''(0) \, [-]$';
ax13.YLim = [9 14];
ax13.YLimMode = 'manual';
ax13.YTickLabelMode = 'manual';
ax13.YTickMode = 'manual';
% ax13.YTick = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35];
% ax13.YTickLabel = {'0.0'; '0.05'; '0.10'; '0.15'; '0.20'; '0.25'; '0.30'; '0.35'};
ax13.YTick = [9 10 11 12 13 14];
ax13.YTickLabel = {'9'; '10'; '11'; '12'; '13'; '14'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure hp - Mach-number profiles for Presentation
hp = figure;
axp = gca;
axp.TickLabelInterpreter = 'latex';

axp.XLabel.Interpreter = 'latex';
axp.XLabel.Color = Color;
axp.XLabel.String = '$X \, [\mathrm{m}]$';
axp.XLim = [0 35];
axp.XTick = [0 5 10 15 20 25 30 35];
axp.XTickLabel = {'0.0'; '5.0'; '10.0'; '15.0'; '20.0'; '25.0'; '30.0'; '35.0'};

axp.YLabel.Interpreter = 'latex';
axp.YLabel.Color = Color;
axp.YLabel.String = '$\mathrm{Ma} \, [-]$';
axp.YLim = [0 3];
axp.YTick = [0 0.5 1 1.5 2 2.5 3];
axp.YTickLabel = {'0.0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'; '3.0'};

axp.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
set(gca,'FontSize',11)
box on
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% cf and cd cannot be distinguished
% % Figure 12 (two y-axes) - cf-Cd plot (two y-axes)
% h12 = figure;
% 
% ax12 = gca;
% 
% yyaxis left
% ax12.YAxis(1).Color = Color;
% ax12.YAxis(1).FontSize = FontSize;
% ax12.YAxis(1).Limits = [1 2];
% ax12.YAxis(1).TickLabelsMode = 'manual';
% ax12.YAxis(1).TickValues = [1 1.2 1.4 1.6 1.8 2];
% ax12.YAxis(1).TickLabels = {'1.0'; '1.2'; '1.4'; '1.6'; '1.8'; '2.0'};
% ax12.YAxis(1).TickLabelInterpreter = 'latex';
% ax12.YAxis(1).Label.FontSize = FontSize; % geen effect
% ax12.YAxis(1).Label.FontWeight = 'normal';
% ax12.YAxis(1).Label.Color = Color;
% ax12.YAxis(1).Label.LineWidth = LineWidth;
% ax12.YAxis(1).Label.Interpreter = 'latex';
% ax12.YAxis(1).Label.String = '$c_{\mathrm{f}} \, \times \, 10^{-3} \, [-]$';
% 
% yyaxis right
% ax12.YAxis(2).Color = Color;
% ax12.YAxis(2).FontSize = FontSize;
% ax12.YAxis(2).Limits = [0 1];
% ax12.YAxis(2).TickLabelsMode = 'manual';
% ax12.YAxis(2).TickValues = [0 0.2 0.4 0.6 0.8 1.0];
% ax12.YAxis(2).TickLabels = {'0.0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'};
% ax12.YAxis(2).TickLabelInterpreter = 'latex';
% ax12.YAxis(2).Label.FontSize = FontSize; % geen effect
% ax12.YAxis(2).Label.FontWeight = 'normal';
% ax12.YAxis(2).Label.Color = Color;
% ax12.YAxis(2).Label.LineWidth = LineWidth;
% ax12.YAxis(2).Label.Interpreter = 'latex';
% ax12.YAxis(2).Label.String = '$C_{\mathrm{d}} \, \times \, 10^{-3} \, [-]$';
% 
% ax12.XLabel.Color = Color;
% ax12.XLabel.FontSize = FontSize;
% ax12.XLimMode = 'manual';
% ax12.XLim = [0 3];
% ax12.XTickLabelMode = 'auto';
% ax12.XTick = [0 0.5 1.0 1.5 2.0 2.5 3.0];
% ax12.XTickLabel = {'0.0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'; '3.0'};
% ax12.TickLabelInterpreter = 'latex';
% ax12.XLabel.FontSize = FontSize; % geen effect
% ax12.XLabel.FontWeight = 'normal';
% ax12.XLabel.Color = Color;
% ax12.XLabel.LineWidth = LineWidth;
% ax12.XLabel.Interpreter = 'latex';
% ax12.XLabel.String = '$\mathrm{Ma} \, [-]$';
% ax12.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
% 
% set(gca,'FontSize',11) % because fontsize doesn't have any effect on results
% box on
% hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% labelfontsizemultiplier ???
% TitleFontSizeMultiplier ???
% TitleFontWeight
% XColor: [1×3 double]
% XColorMode: 'auto'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run from 1 to 4, always! otherwise you use the same data in the next
% calculation!!!
n = 4; % number of cases (and thus number of tables)
% variables for table data storage parameters
delta = NaN*ones(1,n);
% delta_ast = NaN*ones(n,1);
% theta = NaN*ones(n,1);
% Re = NaN*ones(n,1);
% Re_x = NaN*ones(n,1);
% Re_theta = NaN*ones(n,1);
% Cf = NaN*ones(n,1);
% Cd = NaN*ones(n,1);
% r = NaN*ones(n,1);
cw = NaN*ones(1,n);
Cw = NaN*ones(1,n);
shear_par = NaN*ones(1,n);
%
%% Tables

NewFolder = cd;
cd('..\Thesis_Report_Tables')
% fid = fopen('tables/tabMa02-14-22-28.tex','w');
fid = fopen('tabMa02-14-22-28.tex','w');
% % % fid = fopen('tabMa02-14-22-28_Test.tex','w');
% fid = fopen('Test_Validation_Case.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Winter & Gaudet Validation case tables (4 in 1 file) for presenting simulation results with experimental data \n');

%% Calculation
fprintf('Starting %1.0f calculations\n',n)

counter = 0;
SET.ITMAX0 = 10;
OPT.RLAM = 0;

for s = n:-1:1%1:4 %! since P24 and P26 have the same Mach-density input, currently P26 is used, which has nicer/smoother experimental data
    
    counter = counter + 1;
    
    cd('..\Validation_Study')
%     OldFolder = cd(NewFolder);
    
    if s == 1
        run('INPUT/INPUT_Profile02_Winter1973');
        run('Data/Exp_Data_Profile02_Winter1973');
%         table_caption = '    \\caption{Mach $0.2$ Comparison of simulation results with experimental data of profile 2, with recovery factor $r=0.89$.}\n';
        table_label = '    \\label{tab:Ma02}\n';
    elseif s == 2
        run('INPUT/INPUT_Profile12_Winter1973');
        run('Data/Exp_Data_Profile12_Winter1973');
%         table_caption = '    \\caption{Mach $1.4$ Comparison of simulation results with experimental data of profile 12.}\n';
        table_label = '    \\label{tab:Ma14}\n';
    elseif s == 3 % density of this one is incorrect
        run('INPUT/INPUT_Profile19_Winter1973');
        run('Data/Exp_Data_Profile19_Winter1973');
%         table_caption = '    \\caption{Mach $2.2$ Comparison of simulation results with experimental data of profile 19.}\n';
        table_label = '    \\label{tab:Ma22}\n';
    elseif s == 4
%         run('INPUT/INPUT_Profile24_Winter1973');
%         run('Data/Exp_Data_Profile24_Winter1973');
%     elseif s == 5
        run('INPUT/INPUT_Profile26_Winter1973');
        run('Data/Exp_Data_Profile26_Winter1973');
%         table_caption = '    \\caption{Mach $2.8$ Comparison of simulation results with experimental data of profile 26.}\n';
        table_label = '    \\label{tab:Ma26}\n';
    end
    
    %%%%%%%%%%%%%%%
    cd ..\.
    
    % PRECAL
    [X,NST,HVR,EDG,FRS,FLD,TCC,SOL,BLC,MON] = PRECAL(INP,OPT,SET,FLD,TCC);
    
    % Grid generation
    [NP,GRD] = GRID(GRD,SET);
    
    NS = 1;
    
    % initial profiles
    [sol,solprev] = IVPL(NS,GRD,HVR,FLD,SET); % needs to be here because of use of HVR
    
    % Run code (note: sol and solprev remain the same outside CSM-file and are re-used)
    [BLC,FLP,SOL,TCC,MON] = CSM(X,NP,NS,NST,GRD,HVR,EDG,FRS,FLD,TCC,sol,solprev,SOL,BLC,INP,OPT,SET,MON);
    
%     cd Validation_Study
%     %%%%%%%%%%%%%%%
    
    % Mach- density ratio plot
    figure(h0) % Mach-density ratio (pressure history)
    yyaxis left
    plot(X,EDG.rhoE./EDG.rhoE(end),'Color',Color,'LineStyle',LineStyle{s},'Marker','none','DisplayName',LegendName{s},'HandleVisibility','on')
    yyaxis right
    plot(X,EDG.MaE,'Color',Color,'LineStyle',LineStyle{s},'Marker','none','HandleVisibility','off')
    plot(X(SET.NTR),EDG.MaE(SET.NTR),'Color',Color,'Marker','*','HandleVisibility','off')
    
    figure(hp)
    plot(X,EDG.MaE,'Color',Color,'LineStyle',LineStyle{s},'Marker','none','DisplayName',LegendName{s},'HandleVisibility','on')
    plot(X(SET.NTR),EDG.MaE(SET.NTR),'Color',Color,'Marker','*','HandleVisibility','off')
%     plot(X,EDG.MaE,'Color',Color,'LineStyle',LineStyle{s},'Marker','none','HandleVisibility','on')
%     plot(X(SET.NTR),EDG.MaE(SET.NTR),'Color',Color,'Marker','*','HandleVisibility','off')
    
    figure(h1) % Mach-number profiles (obtained from total pressure measurements)
    plot(BLC.y(1:length(SOL{end}.u),end),EDG.UE(end).*SOL{end}.u./(sqrt(FLD.gamma.*FLD.Rsg.*FLP{end}.T)),'Color',Color,'LineStyle',LineStyle{s},'Marker','none','HandleVisibility','off')
    plot(EXP.y,EXP.Ma,'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName',LegendName{s})
    
    figure(h2) % Velocity profiles (obtained from total temperature measurements)
    plot(BLC.y(1:length(SOL{end}.u)),SOL{end}.u + (s - 1)*0.1,'Color',Color,'LineStyle','-','Marker','none','HandleVisibility','off')
%     if s == 1
%         plot(y(2:end),EDG.UE(end).*u_ue089(2:end),'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName','Ma 0.2, $r=0.89$','HandleVisibility','on')
%         plot(y(2:end),EDG.UE(end).*u_ue100(2:end),'Color',Color,'LineStyle','none','Marker',Marker{s+4},'DisplayName','Ma 0.2, $r=1.00$','HandleVisibility','on')
%     else
    plot(EXP.y,EXP.u_ue + (s - 1)*0.1,'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName',LegendName{s})
%     end
    
    figure(h3) % Density profiles (obtained from total temperature measurements)
    plot(BLC.y(1:length(SOL{end}.c)),SOL{end}.c + (s - 1)*0.1,'Color',Color,'LineStyle','-','Marker','none','HandleVisibility','off')
%     if s == 1
%         plot(y,rhoE.*rho_rhoe089,'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName','Ma 0.2, $r=0.89$','HandleVisibility','on')
%         plot(y,rhoE.*rho_rhoe100,'Color',Color,'LineStyle','none','Marker',Marker{s+4},'DisplayName','Ma 0.2, $r=1.00$','HandleVisibility','on')
%     else
        plot(EXP.y,1./EXP.rho_rhoe + (s - 1)*0.1,'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName',LegendName{s},'HandleVisibility','on')
%     end

figure(h3b) % Density profiles (obtained from total temperature measurements)
    plot(BLC.y(1:length(SOL{end}.c)),1./SOL{end}.c - (s - 1)*0.05,'Color',Color,'LineStyle','-','Marker','none','HandleVisibility','off')
%     if s == 1 %%% - (s - 1)*0.05
%         plot(y,rhoE.*rho_rhoe089,'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName','Ma 0.2, $r=0.89$','HandleVisibility','on')
%         plot(y,rhoE.*rho_rhoe100,'Color',Color,'LineStyle','none','Marker',Marker{s+4},'DisplayName','Ma 0.2, $r=1.00$','HandleVisibility','on')
%     else %%% - (s - 1)*0.05
        plot(EXP.y,EXP.rho_rhoe - (s - 1)*0.05,'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName',LegendName{s},'HandleVisibility','on')
%     end
    
    figure(h3c) % Density profiles - Log Plot (obtained from total temperature measurements)
    semilogy(BLC.y(1:length(SOL{end}.c)),SOL{end}.c - 1,'Color',Color,'LineStyle','-','Marker','none','HandleVisibility','off')
    semilogy(EXP.y,1./EXP.rho_rhoe - 1,'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName',LegendName{s},'HandleVisibility','on')
    %  + (s - 1)*0.05
    
        figure(h4a) % Log-law velocity profiles with rhoE instead of rhoW
    % Y+ and U+:
    %     semilogx(BLC.y_plus(1:length(SOL{end}.c),end),BLC.u_plus(1:length(SOL{end}.u),end) + (s - 1)*5,'Color',Color,'LineStyle',LineStyle{s},'Marker','none','HandleVisibility','off')
    % with rhoE:    
    semilogx(BLC.y_plus(1:length(SOL{end}.c),end).*sqrt(SOL{end}.c)./sqrt(SOL{end}.c(1)),BLC.u_plus(1:length(SOL{end}.c),end)./sqrt(SOL{end}.c).*sqrt(SOL{end}.c(1)) + (s - 1)*5,'Color',Color,'LineStyle',LineStyle{s},'Marker','none','HandleVisibility','off')
    
% % %     if s == 1
% % %         semilogx(y.*sqrt((Cf*rhoE*EDG.UE(end)^2/2)./(rhoE*rho_rhoe089(1)))./((1.45e-6*((1./rho_rhoe089*TE).^1.5)./((1./rho_rhoe089*TE) + FLD.SLV))./(rhoE*rho_rhoe089)),u_ue089*EDG.UE(end)/sqrt((Cf*rhoE*EDG.UE(end)^2/2)/(rhoE*rho_rhoe089(1))),...
% % %             'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName','Ma 0.2, $r=0.89$','HandleVisibility','on')
% % %         semilogx(y.*sqrt((Cf*rhoE*EDG.UE(end)^2/2)./(rhoE*rho_rhoe100(1)))./((1.45e-6*((1./rho_rhoe100*TE).^1.5)./((1./rho_rhoe100*TE) + FLD.SLV))./(rhoE*rho_rhoe100)),u_ue100*EDG.UE(end)/sqrt((Cf*rhoE*EDG.UE(end)^2/2)/(rhoE*rho_rhoe100(1))),...
% % %             'Color',Color,'LineStyle','none','Marker',Marker{s+4},'DisplayName','Ma 0.2, $r=1.00$','HandleVisibility','on')
% % %     else
% % %         semilogx(y.*sqrt((Cf*rhoE*EDG.UE(end)^2/2)./(rhoE*rho_rhoe(1)))./((1.45e-6*((1./rho_rhoe*TE).^1.5)./((1./rho_rhoe*TE) + FLD.SLV))./(rhoE*rho_rhoe)),u_ue*EDG.UE(end)/sqrt((Cf*rhoE*EDG.UE(end)^2/2)/(rhoE*rho_rhoe(1))) + (s - 1)*5,...
% % %             'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName',LegendName{s})
% % %     end
    % U+ correction:
%     Fu = 1;%sqrt(EXP.TW/EXP.TE);
%     Fu = 1;
    semilogx(sqrt(1/2)*EXP.rho_rhoe.^2.*EXP.ReE.*sqrt(EXP.cf).*EXP.y.*((1./EXP.rho_rhoe + EXP.S*(1 + (EXP.gamma - 1)/2*EXP.MaE^2)/EXP.T0)/(1 + EXP.S*(1 + (EXP.gamma - 1)/2*EXP.MaE^2)/EXP.T0)),sqrt(2).*sqrt(EXP.rho_rhoe).*EXP.u_ue./sqrt(EXP.cf)  + (s - 1)*5,...
        'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName',LegendName{s})

    figure(h4b) % Log-law velocity profiles (compressible! with estimated TW/TE from Fernholz & Finley (1970))
    % Y+ and U+:
    %     semilogx(BLC.y_plus(1:length(SOL{end}.c),end),BLC.u_plus(1:length(SOL{end}.u),end) + (s - 1)*5,'Color',Color,'LineStyle',LineStyle{s},'Marker','none','HandleVisibility','off')
    % with rhoE:    
    semilogx(BLC.y_plus(1:length(SOL{end}.c),end),BLC.u_plus(1:length(SOL{end}.u),end) + (s - 1)*5,'Color',Color,'LineStyle',LineStyle{s},'Marker','none','HandleVisibility','off')
    
% % %     if s == 1
% % %         semilogx(y.*sqrt((Cf*rhoE*EDG.UE(end)^2/2)./(rhoE*rho_rhoe089(1)))./((1.45e-6*((1./rho_rhoe089*TE).^1.5)./((1./rho_rhoe089*TE) + FLD.SLV))./(rhoE*rho_rhoe089)),u_ue089*EDG.UE(end)/sqrt((Cf*rhoE*EDG.UE(end)^2/2)/(rhoE*rho_rhoe089(1))),...
% % %             'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName','Ma 0.2, $r=0.89$','HandleVisibility','on')
% % %         semilogx(y.*sqrt((Cf*rhoE*EDG.UE(end)^2/2)./(rhoE*rho_rhoe100(1)))./((1.45e-6*((1./rho_rhoe100*TE).^1.5)./((1./rho_rhoe100*TE) + FLD.SLV))./(rhoE*rho_rhoe100)),u_ue100*EDG.UE(end)/sqrt((Cf*rhoE*EDG.UE(end)^2/2)/(rhoE*rho_rhoe100(1))),...
% % %             'Color',Color,'LineStyle','none','Marker',Marker{s+4},'DisplayName','Ma 0.2, $r=1.00$','HandleVisibility','on')
% % %     else
% % %         semilogx(y.*sqrt((Cf*rhoE*EDG.UE(end)^2/2)./(rhoE*rho_rhoe(1)))./((1.45e-6*((1./rho_rhoe*TE).^1.5)./((1./rho_rhoe*TE) + FLD.SLV))./(rhoE*rho_rhoe)),u_ue*EDG.UE(end)/sqrt((Cf*rhoE*EDG.UE(end)^2/2)/(rhoE*rho_rhoe(1))) + (s - 1)*5,...
% % %             'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName',LegendName{s})
% % %     end
    % U+ correction:
    Fu = sqrt(EXP.TE/EXP.TW)./sqrt(EXP.rho_rhoe');
%     Fu = 1;
    semilogx(1./Fu.*sqrt(1/2).*EXP.rho_rhoe'.^2.*EXP.ReE.*sqrt(EXP.cf).*EXP.y'.*((1./EXP.rho_rhoe' + EXP.S*(1 + (EXP.gamma - 1)/2*EXP.MaE^2)/EXP.T0)/(1 + EXP.S*(1 + (EXP.gamma - 1)/2*EXP.MaE^2)/EXP.T0)),Fu.*sqrt(2).*sqrt(EXP.rho_rhoe').*EXP.u_ue'./sqrt(EXP.cf)  + (s - 1)*5,...
        'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName',LegendName{s})
        
    figure(h5) % Log-law density profiles (compressible! with estimated TW/TE from Fernholz & Finley (1970))
    % Y+ and U+:
    %     semilogx(BLC.y_plus(1:length(SOL{end}.c),end),BLC.u_plus(1:length(SOL{end}.u),end) + (s - 1)*5,'Color',Color,'LineStyle',LineStyle{s},'Marker','none','HandleVisibility','off')
    % with rhoE:    
    semilogx(BLC.y_plus(1:length(SOL{end}.c),end),(SOL{end}.c(1) - SOL{end}.c)./(SOL{end}.c.*(SOL{end}.c(1) - 1)) + (s - 1)*0.2,'Color',Color,'LineStyle',LineStyle{s},'Marker','none','HandleVisibility','off')
    
% % %     if s == 1
% % %         semilogx(y.*sqrt((Cf*rhoE*EDG.UE(end)^2/2)./(rhoE*rho_rhoe089(1)))./((1.45e-6*((1./rho_rhoe089*TE).^1.5)./((1./rho_rhoe089*TE) + FLD.SLV))./(rhoE*rho_rhoe089)),u_ue089*EDG.UE(end)/sqrt((Cf*rhoE*EDG.UE(end)^2/2)/(rhoE*rho_rhoe089(1))),...
% % %             'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName','Ma 0.2, $r=0.89$','HandleVisibility','on')
% % %         semilogx(y.*sqrt((Cf*rhoE*EDG.UE(end)^2/2)./(rhoE*rho_rhoe100(1)))./((1.45e-6*((1./rho_rhoe100*TE).^1.5)./((1./rho_rhoe100*TE) + FLD.SLV))./(rhoE*rho_rhoe100)),u_ue100*EDG.UE(end)/sqrt((Cf*rhoE*EDG.UE(end)^2/2)/(rhoE*rho_rhoe100(1))),...
% % %             'Color',Color,'LineStyle','none','Marker',Marker{s+4},'DisplayName','Ma 0.2, $r=1.00$','HandleVisibility','on')
% % %     else
% % %         semilogx(y.*sqrt((Cf*rhoE*EDG.UE(end)^2/2)./(rhoE*rho_rhoe(1)))./((1.45e-6*((1./rho_rhoe*TE).^1.5)./((1./rho_rhoe*TE) + FLD.SLV))./(rhoE*rho_rhoe)),u_ue*EDG.UE(end)/sqrt((Cf*rhoE*EDG.UE(end)^2/2)/(rhoE*rho_rhoe(1))) + (s - 1)*5,...
% % %             'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName',LegendName{s})
% % %     end
    % U+ correction:
    Fu = sqrt(EXP.TE/EXP.TW)./sqrt(EXP.rho_rhoe');
%     Fu = 1;
    semilogx(1./Fu.*sqrt(1/2).*EXP.rho_rhoe'.^2.*EXP.ReE.*sqrt(EXP.cf).*EXP.y'.*((1./EXP.rho_rhoe' + EXP.S*(1 + (EXP.gamma - 1)/2*EXP.MaE^2)/EXP.T0)/(1 + EXP.S*(1 + (EXP.gamma - 1)/2*EXP.MaE^2)/EXP.T0)),(EXP.TW/EXP.TE - 1./EXP.rho_rhoe')./(1./EXP.rho_rhoe'.*(EXP.TW/EXP.TE - 1))  + (s - 1)*0.2,...
        'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName',LegendName{s})
    
% % %     figure(h6) % Dimensionless Velocity profiles
% % %     plot(BLC.y(1:length(SOL{end}.u),SOL{end}.u,'Color',Color,'LineStyle','-','Marker','none','HandleVisibility','off')
% % % %     if s == 1 % /BLC.delta(end)
% % % %         plot(y(2:end),u_ue089(2:end),'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName','Ma 0.2, $r=0.89$','HandleVisibility','on')
% % % %         plot(y(2:end),u_ue100(2:end),'Color',Color,'LineStyle','none','Marker',Marker{s+4},'DisplayName','Ma 0.2, $r=1.00$','HandleVisibility','on')
% % % %     else
% % %         plot(y(2:end),u_ue(2:end) + (s - 1)*5,'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName',LegendName{s})
% % % %     end
    
% % %     figure(h7) % Dimensionless Density profiles
% % %     plot(BLC.y(1:length(SOL{end}.c),1./SOL{end}.c,'Color',Color,'LineStyle','-','Marker','none','HandleVisibility','off')
% % % %     if s == 1
% % % %         plot(y,rho_rhoe089,'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName','Ma 0.2, $r=0.89$','HandleVisibility','on')
% % % %         plot(y,rho_rhoe100,'Color',Color,'LineStyle','none','Marker',Marker{s+4},'DisplayName','Ma 0.2, $r=1.00$','HandleVisibility','on')
% % % %     else
% % %         plot(y,rho_rhoe,'Color',Color,'LineStyle','none','Marker',Marker{s},'DisplayName',LegendName{s},'HandleVisibility','on')
% % % %     end
    
    % Calculate and save table data (for plots) that is calculated from experimental data
    
    % Estimate BL velocity thickness with linear interpolation:
    p = find(EXP.u_ue>0.99);
%     p = P(1);
    delta(1,s) = (0.99 - EXP.u_ue(p(1) - 1)) / ((EXP.u_ue(p(1)) - EXP.u_ue(p(1) - 1))/(EXP.y(p(1)) - EXP.y(p(1) - 1))) + EXP.y(p(1) - 1);
    
    % Estimate cw and Cw:
    cw(1,s) = EXP.TW/EXP.TE; % [-], data from Fernholz & Finley (1977), cw = rhoE/rhoW = TW/TE
    Cw(1,s) = cw(1,s)^(3/2)*(EXP.TE + EXP.S)/(EXP.TW + EXP.S)/cw(1,s); % [-], [-], data from Fernholz & Finley (1977), Cw = (rhoW muW)/(rhoE muE) = Suth(TW)/Suth(TE) / cw
    
    % store shear parameter at wall
    shear_par(1,s) = SOL{end}.v(1);
    
%     figure(h8) % cf (left) - Cd (right) as function of Mach
%     yyaxis left
%     plot(EDG.MaE(end),BLC.Cf_L(end)*10^3,'Color',Color,'Marker',Marker{s},'MarkerFaceColor',Color,'HandleVisibility','off')
%     plot(EXP.MaE,EXP.cf,'Color',Color,'Marker',Marker{s},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName{s},'HandleVisibility','on')
%     yyaxis right
%     plot(EDG.MaE(end),BLC.Cd(end)*10^3,'Color',Color,'Marker',Marker{s},'MarkerFaceColor',Color,'HandleVisibility','off')
% %     plot(EXP.MaE,,'Color',Color,'Marker',Marker{s},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName{s},'HandleVisibility','on')

%     figure(h7) % delta_ast and theta (delta skipped due to large uncertainty)
%     % delta_ast
%     plot(EDG.MaE(end),BLC.delta_ast(end)*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor',Color,'HandleVisibility','off')
%     plot(EXP.MaE,EXP.delta_ast*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName2{s})
%     % theta
%     plot(EDG.MaE(end),BLC.theta(end)*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s + 4},'MarkerFaceColor',Color,'HandleVisibility','off')
%     plot(EXP.MaE,EXP.theta*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s + 4},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName2{s + 4})

    figure(h7) % delta_ast and theta (delta skipped due to large uncertainty)
    % delta_ast
    yyaxis left
    plot(EDG.MaE(end),BLC.delta_ast(end)*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor',Color,'HandleVisibility','off')
    plot(EXP.MaE,EXP.delta_ast*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName2{s})
    % theta
    yyaxis right
    plot(EDG.MaE(end),BLC.theta(end)*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s + 4},'MarkerFaceColor',Color,'HandleVisibility','off')
    plot(EXP.MaE,EXP.theta*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s + 4},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName2{s + 4})

    figure(h8) % cf (left) - Cd (right) as function of Mach
    plot(EDG.MaE(end),BLC.Cf_L(end)*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor',Color,'HandleVisibility','off')
    plot(EXP.MaE,EXP.cf*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName{s})
%     plot(EDG.MaE(end),BLC.Cd(end)*10^3,'Color',Color,'Marker',Marker{s + 4},'MarkerFaceColor',Color,'HandleVisibility','on','DisplayName','Cd')
    %     plot(EDG.MaE(end),BLC.Cd(end)*10^3,'Color',Color,'Marker',Marker{s},'MarkerFaceColor',Color,'HandleVisibility','off')
%     plot(EXP.MaE,,'Color',Color,'Marker',Marker{s},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName{s},'HandleVisibility','on')
    
    figure(h9)
    % c_w
    yyaxis left
    plot(EDG.MaE(end),SOL{end}.c(1),'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor',Color,'HandleVisibility','off')
    plot(EXP.MaE,cw(1,s),'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName3{s})
    % C_w
    yyaxis right
    plot(EDG.MaE(end),FLP{end}.C(1),'Color',Color,'LineStyle','none','Marker',Marker{s + 4},'MarkerFaceColor',Color,'HandleVisibility','off')
    plot(EXP.MaE,Cw(1,s),'Color',Color,'LineStyle','none','Marker',Marker{s + 4},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName3{s + 4})
    
    figure(h10) % shape factor H
    plot(EDG.MaE(end),BLC.H(end),'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor',Color,'HandleVisibility','off')
    plot(EXP.MaE,EXP.H12,'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName{s})
    
    figure(h11) % recovery factor
    plot(EDG.MaE(end),BLC.Trecovery(end),'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor',Color,'HandleVisibility','off')
    plot(EXP.MaE,0.896,'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName{s})
    
    figure(h12) % cf and Cd as function of Mach
    plot(EDG.MaE(end),BLC.Cf_L(end)*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor',Color,'HandleVisibility','off')
    plot(EXP.MaE,EXP.cf*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName{s})
    plot(EDG.MaE(end),BLC.Cd(end)*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s + 4},'MarkerFaceColor',Color,'HandleVisibility','on','DisplayName','Cd')
    
    figure(h13) % shear parameter at wall
    plot(EDG.MaE(end),SOL{end}.v(1),'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor',Color,'HandleVisibility','off')
    
    % cf and cd cannot be distinguished
%     figure(h12) % cf (left) - Cd (right) as function of Mach
%     % cf
%     yyaxis left
%     plot(EDG.MaE(end),BLC.Cf_L(end)*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor',Color,'HandleVisibility','off')
%     plot(EXP.MaE,EXP.cf*10^3,'Color',Color,'LineStyle','none','Marker',Marker{s},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName{s})
%     % Cd
%     yyaxis right
%     plot(EDG.MaE(end),BLC.Cd(end)*10^3,'Color',Color,'Marker',Marker{s + 4},'MarkerFaceColor',Color,'HandleVisibility','off')
% %     plot(EXP.MaE,,'Color',Color,'Marker',Marker{s},'MarkerFaceColor','none','HandleVisibility','on','DisplayName',LegendName{s},'HandleVisibility','on')
    
    % PLOTFILE
%     cd ..
% keyboard % check data name string!!! (make sure newly stored data is exactly the same as old data!
    STORDATA
%     cd ./Validation_Study
        
%     cd Validation_Study
    %%%%%%%%%%%%%%%
    
%     % Save table data (for plots) that is calculated from experimental data
%     EXP.delta
%     delta_ast
%     theta
%     Re_x
%     Re_theta
%     cw
%     Cw

    % Print Table
    cd('Thesis_Report_Tables')
    
    if counter > 1
        fprintf(fid, '\n');
        fprintf(fid, '\n');
        fprintf(fid, '\\vspace{0.15cm}\n'); % Problem: vertical spacing doesn't work when using landscape
    end
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    fprintf(fid, '\\begin{adjustbox}{width=1.1\\textwidth,center=\\textwidth}\n');
    fprintf(fid, '\\begin{threeparttable}\n');
    fprintf(fid, '    \\centering\n');
    fprintf(fid, '    \\caption{Mach $%1.1f$: Comparison of simulated results with experimental data of profile %2.0f with Mach-number $%1.1f$ and initial total conditions $P_0^{\\ast} = %1.5f \\times \\, 10^{5} \\, \\si{\\pascal}$ and $T_0 = \\SI{%3.2f}{\\kelvin}$ from Winter \\& Gaudet\\cite{winter1970turbulent} and complemented by Fernholz \\& Finley\\cite{fernholz1977critical}.}\n',EXP.MaE,EXP.Pn,EXP.MaE,EXP.P0/10^5,EXP.T0);
%     fprintf(fid, table_caption);
    fprintf(fid, table_label);
%     fprintf(fid, '    \\begin{tabular}{l S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.3] S[table-format=-2.2] S[table-format=-2.2]}\n'); %%% including Cd
    fprintf(fid, '    \\begin{tabular}{l S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.3] S[table-format=-2.2] S[table-format=-2.2]}\n'); % excluding Cd
    fprintf(fid, '        \\toprule\n\n');
%     fprintf(fid, '        \\multicolumn{12}{c}{\textnormal{Simulation results (Mach $0.2035$) compared with (derived) experimental data profile 2 (Mach $0.2007$, $r=1.00$, $\rho_e = 1.7667 kg/m^3$}} \\\\\n');
%     fprintf(fid, '        \\cmidrule(lr){1-12}\n\n');
    fprintf(fid, '        \\multicolumn{1}{c}{}                                             &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{\\bm{$\\delta$}}                              &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{\\bm{$\\delta^*$}}                            &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{\\bm{$\\theta$}}                              &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{\\bm{$H$}}                                    &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{\\bm{$\\mathrm{Re}$}}                         &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{\\bm{$\\mathrm{Re}_x$}}                       &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{\\bm{$\\mathrm{Re}_{\\theta}$}}               &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{\\bm{$C_{\\mathrm{f}}$}}                      &\n');
%     fprintf(fid, '        \\multicolumn{1}{c}{\\bm{$C_{\\mathrm{d}}$}}                      &\n'); %%% Cd
    fprintf(fid, '        \\multicolumn{1}{c}{\\bm{$r$}}                                    &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{\\bm{$c_{\\mathrm{w}}$}}                      &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{\\bm{$C_{\\mathrm{w}}$}}                   \\\\\n');
    fprintf(fid, '        \\multicolumn{1}{c}{}                                             &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\, [\\si{\\metre}]$}    &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\, [\\si{\\metre}]$}    &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\, [\\si{\\metre}]$}    &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{$[-]$}                                        &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{7} \\, [\\si{\\per\\metre}]$}   &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{7} \\, [-]$}                 &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{3} \\, [-]$}                 &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\, [-]$}                &\n');
%     fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\, [-]$}                &\n'); % Cd
    fprintf(fid, '        \\multicolumn{1}{c}{$[-]$}                                        &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{$[-]$}                                        &\n');
    fprintf(fid, '        \\multicolumn{1}{c}{$[-]$}                                     \\\\\n');
    fprintf(fid, '        \\midrule\n\n');
    fprintf(fid, '        %% Data\n');

    % including Cd
%     fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{Experimental data}}    &   %2.2f\\tnote{$\\diamond$} &   %2.2f   &   %2.2f   & %1.2f   &   %2.2f   &   %2.2f   &  %2.2f   & %1.3f   &   NA      &   0.896\\tnote{$\\ast$}	&   %2.2f\\tnote{$\\ast$}   &   %2.2f\\tnote{$\\ast$}	\\\\\n',delta(1,n)*10^3,EXP.delta_ast*10^3,EXP.theta*10^3,EXP.H12,EXP.ReE*10^-8,EXP.ReE*X(end)*10^-7,EXP.ReE*EXP.theta*10^-3,EXP.cf*10^3,cw(1,n),Cw(1,n));
%     fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{Simulation}}           &   %2.2f                     &   %2.2f   &   %2.2f   & %1.2f   &   %2.2f   &   %2.2f   &  %2.2f   & %1.3f   &   %1.3f   &   %1.3f                    &   %2.2f                   &   %2.2f                   \\\\\n',BLC.delta(end)*10^3,BLC.delta_ast(end)*10^3,BLC.theta(end)*10^3,BLC.H(end),BLC.Re_x(end)/X(end)*10^-8,BLC.Re_x(end)*10^-7,BLC.Re_theta(end)*10^-3,BLC.Cf_L(end)*10^3,BLC.Cd(end)*10^3,BLC.Trecovery(end),SOL{end}.c(1),FLP{end}.C(1));
%     fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{Deviation [\\%%]}}       &   %2.2f                     &   %2.2f   &   %2.2f   & %2.2f   &   %2.2f   &   %2.2f   &  %2.2f   & %2.2f   &   NA      &   NA                       &   %2.2f                   &   %2.2f                   \\\\\n',(BLC.delta(end) - delta(1,n))/delta(1,n)*100,(BLC.delta_ast(end) - EXP.delta_ast)/EXP.delta_ast*100,(BLC.theta(end) - EXP.theta)/EXP.theta*100,(BLC.H(end) - EXP.H12)/EXP.H12*100,(BLC.Re_x(end)/X(end) - EXP.ReE)/EXP.ReE*100,(BLC.Re_x(end) - EXP.ReE*X(end))/(EXP.ReE*X(end))*100,(BLC.Re_theta(end) - EXP.ReE*EXP.theta)/(EXP.ReE*EXP.theta)*100,(BLC.Cf_L(end) - EXP.cf)/EXP.cf*100,(SOL{end}.c(1) - cw(1,n))/cw(1,n)*100,(FLP{end}.C(1) - Cw(1,n))/Cw(1,n)*100);

    % without Cd (since not measured)
    fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{Experimental data}}    &   %2.2f\\tnote{$\\diamond$} &   %2.2f   &   %2.2f   & %1.2f   &   %2.2f   &   %2.2f   &  %2.2f   & %1.3f   &   0.896\\tnote{$\\ast$}	&   %2.2f\\tnote{$\\ast$}   &   %2.2f\\tnote{$\\ast$}	\\\\\n',delta(1,s)*10^3,EXP.delta_ast*10^3,EXP.theta*10^3,EXP.H12,EXP.ReE*10^-7,EXP.ReE*X(end)*10^-7,EXP.ReE*EXP.theta*10^-3,EXP.cf*10^3,cw(1,s),Cw(1,s));
    fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{Simulation}}           &   %2.2f                     &   %2.2f   &   %2.2f   & %1.2f   &   %2.2f   &   %2.2f   &  %2.2f   & %1.3f   &   %1.3f                    &   %2.2f                   &   %2.2f                   \\\\\n',BLC.delta(end)*10^3,BLC.delta_ast(end)*10^3,BLC.theta(end)*10^3,BLC.H(end),BLC.Re_x(end)/X(end)*10^-7,BLC.Re_x(end)*10^-7,BLC.Re_theta(end)*10^-3,BLC.Cf_L(end)*10^3,BLC.Trecovery(end),SOL{end}.c(1),FLP{end}.C(1));
    fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{Deviation [\\%%]}}       &   %2.2f                     &   %2.2f   &   %2.2f   & %2.2f   &   %2.2f   &   %2.2f   &  %2.2f   & %2.2f   &   NA                       &   %2.2f                   &   %2.2f                   \\\\\n',(BLC.delta(end) - delta(1,s))/delta(1,s)*100,(BLC.delta_ast(end) - EXP.delta_ast)/EXP.delta_ast*100,(BLC.theta(end) - EXP.theta)/EXP.theta*100,(BLC.H(end) - EXP.H12)/EXP.H12*100,(BLC.Re_x(end)/X(end) - EXP.ReE)/EXP.ReE*100,(BLC.Re_x(end) - EXP.ReE*X(end))/(EXP.ReE*X(end))*100,(BLC.Re_theta(end) - EXP.ReE*EXP.theta)/(EXP.ReE*EXP.theta)*100,(BLC.Cf_L(end) - EXP.cf)/EXP.cf*100,(SOL{end}.c(1) - cw(1,s))/cw(1,s)*100,(FLP{end}.C(1) - Cw(1,s))/Cw(1,s)*100);
    
    fprintf(fid, '    \\bottomrule\n\n');
    fprintf(fid, '    \\end{tabular}\n');
%     if counter == n
    fprintf(fid, '    \\begin{tablenotes}\n');
    fprintf(fid, '        \\item[$\\diamond$] Estimated by linear interpolation of the experimental data.\\\\\n');
    fprintf(fid, '        \\item[$\\ast$] Estimated by Fernholz \\& Finley\\cite{fernholz1977critical}.\n');
    fprintf(fid, '    \\end{tablenotes}\n');
%     end
    fprintf(fid, '\\end{threeparttable}\n');
    fprintf(fid, '\\end{adjustbox}\n');
    
end

% fprintf(fid, '%%\\end{sidewaystable}\n');
% fprintf(fid, '%%\\end{table\n');
% fprintf(fid, '\\end{landscape}\n');
fclose(fid);
OldFolder = cd(NewFolder);

fprintf('Tables successfully generated\n')

%% PLOTS

% Legend
figure(h0) % Mach-density plot
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Position', LgndPos1); % 'Location', 'best');
% or: lgd = legend; lgd.Interpreter = 'latex'; lgd.FontSize = FontSize;
% etc. OR: set(lgd, 'Interpreter', 'latex', 'FontSize', FontSize);
% show only first four legend entries:
% % % legend('Ma 0.2','Ma 1.4','Ma 2.2','Ma 2.8')
legend('show')

% Legend
figure(hp) % Mach-number plot
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest');
legend('show')

figure(h1) % Mach-number profiles
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Position', LgndPos3); % 'Location', 'best');
% or: lgd = legend; lgd.Interpreter = 'latex'; lgd.FontSize = FontSize;
% etc. OR: set(lgd, 'Interpreter', 'latex', 'FontSize', FontSize);
% show only first four legend entries:
% % % legend('Ma 0.2','Ma 1.4','Ma 2.2','Ma 2.8')
legend('show')

figure(h2)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Position', LgndPos3); % 'Location', 'best');
legend('show')
% legend()
% https://stackoverflow.com/questions/39103748/matlab-change-order-of-entries-in-figure-legend/39104494#39104494
% plot(...,'HandleVisibility','off')
% plot(...,'DisplayName','Pressure side')
% legend('show')

figure(h3)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northeast'); %, 'Position', LgndPos2); % 'Location', 'best');
legend('show')

figure(h3b)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'southeast'); %, 'Position', LgndPos2); % 'Location', 'best');
legend('show')

figure(h3c)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'southeast'); %, 'Position', LgndPos2); % 'Location', 'best');
legend('show')

figure(h4a)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'southeast'); % 'Location', 'best');
legend('show')

figure(h4b)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'southeast'); % 'Location', 'best');
legend('show')

figure(h5)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'southeast'); % 'Location', 'best');
legend('show')

figure(h6)
legend('off')

figure(h7) % delta_ast and theta
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'east'); % 'Location', 'best');
legend('show')

figure(h8)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'southwest'); % 'Location', 'best');
legend('show')

figure(h9)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest'); % 'Location', 'best');
legend('show')

figure(h10)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest'); % 'Location', 'best');
legend('show')

figure(h11)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest'); % 'Location', 'best');
legend('show')

figure(h12)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest'); % 'Location', 'best');
legend('show')

% figure(h13) % legend not needed/useful
% set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest'); % 'Location', 'best');
% legend('show')

%% Store plots
% keyboard

saveas(h0,'../Thesis_Report_Figures/Mach-density','epsc')
saveas(hp,'../Thesis_Report_Figures/Mach-number','epsc')
saveas(h1,'../Thesis_Report_Figures/Mach','epsc')
saveas(h2,'../Thesis_Report_Figures/Velocity','epsc')
saveas(h3,'../Thesis_Report_Figures/Density','epsc')
% saveas(h3b,'../Thesis_Report_Figures/Density','epsc')
% saveas(h3c,'../Thesis_Report_Figures/Density','epsc')
saveas(h4a,'../Thesis_Report_Figures/LogLawVelocity_Pure','epsc')
saveas(h4b,'../Thesis_Report_Figures/LogLawVelocity','epsc')
saveas(h5,'../Thesis_Report_Figures/LogLawDensity','epsc')
saveas(h6,'../Thesis_Report_Figures/NozzleGeometry','epsc')
saveas(h7,'../Thesis_Report_Figures/BLThickness','epsc')
saveas(h8,'../Thesis_Report_Figures/SkinFriction','epsc')
saveas(h9,'../Thesis_Report_Figures/WallProp','epsc')
saveas(h10,'../Thesis_Report_Figures/ShapeFactor','epsc')
% saveas(h11,'../Thesis_Report_Figures/Trec','epsc')
% saveas(h12,'../Thesis_Report_Figures/cf_cd','epsc')
% saveas(h13,'../Thesis_Report_Figures/shear_par_wall','epsc')
OldFolder = cd(NewFolder);

fprintf('Figures successfully generated\n')

%% Clean-up FluidProp (if used)
% PLOTCHART % Pv-, and Ts-diagram, calculate critical properties inside
% if OPT.CHRT > 0
%     FPCHARTS
% end
if OPT.GASM == 3 % Real Gas
    Cleanup_FluidProp
end

fprintf('End of run\n')

%% Plot Results and generate graphs and tables
% STORDATA
% PLOTFILE
% TABLEGEN
% clear NP tr % not really needed, overwritten? // and other free variables? in case of multiple calculations
