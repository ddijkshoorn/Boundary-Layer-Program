%% Plot ORCHID Data Results Comparing RefProp and StanMix
% Last changes made on 01-04-2020 by DDD
% Last run on 18-05-2020
% StanMix is iPRSV
% RefProp is Helmholtz multiparameter from NIST
% RefProp is taken as reference ('truth') here
% StanMix results in r > 1 and Pr > 1 !!! (which is thought to be bullshit)
% Adapted slightly for upload and checked: 21-01-2022

close all
clear all
clc

%% Settings
% https://www.rapidtables.com/web/color/gray-color.html
Color = [0 0 0];
Color2 = [0.94 0.94 0.94]; % ligth grey (too light for thin lines, used for surfaces)
Color3 = [112 128 144]/255; % slate (darker) grey
Color3 = [192 192 192]/255; % silver
% Color3 = [0.94 0.94 0.94]; % ligth grey TEST
LineStyle = {'-','--','-.',':','-',':'};
Marker = {'none','none','none','none','none','none'}; % Marker = {'+','x','s','d','o'};
LineWidth = 0.5;
LineWidth2 = 1.2;1.0;
LineWidth3 = 0.8; % Linewidth corresponding to Color3
FontSize1 = 9;
FontSize = 11;
FontSize2 = 12;
FontSize3 = 14;
% LgndPos1 = [0.20 0.55 0.20 0.20]; % normalized position
% LgndPos2 = [0.55 0.3 0.3 0.2];
% LegendName = {'Ma 0.2','Ma 1.4','Ma 2.2','Ma 2.8'};

%% Initialize Plots

% figure 1: ZPR-graph
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
ax1.XLim = [0 0.1]; % change to 0.15?
ax1.XLimMode = 'manual';
ax1.XTickLabelMode = 'manual';
ax1.XTickMode = 'manual';
ax1.XTick = [0 0.02 0.04 0.06 0.08 0.10];
ax1.XTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};

ax1.YLabel.Interpreter = 'latex';
ax1.YLabel.FontSize = FontSize; %%
ax1.YLabel.FontWeight = 'normal';
ax1.YLabel.Color = Color;
ax1.YLabel.String = '$Z_{\mathrm{e}}, \, \mathrm{Pr}_{\mathrm{e}}, \, r \, [-]$';
ax1.YLim = [0.5 1.3];
ax1.YLimMode = 'manual';
ax1.YTickLabelMode = 'manual';
ax1.YTickMode = 'manual';
ax1.YTick = [0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3];
ax1.YTickLabel = {'0.5'; '0.6'; '0.7'; '0.8'; '0.9'; '1.0'; '1.1'; '1.2'; '1.3'};
ax1.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
% %
% ax1.YAxis.Exponent = 6;

set(gca,'FontSize',FontSize1)
box on
hold on

% figure 2: (two y-axes) - Static Pressure (input) and Mach-number along BL edge
h2 = figure;
ax2 = gca;

% ax1.YLimMode = 'manual';
% ax1.YTickLabelMode = 'manual';
% ax1.YTickMode = 'manual';

%%% plot first, then adapt to style!?
yyaxis left
ax2.YAxis(1).FontSize = FontSize3;%FontSize2;
ax2.YAxis(1).FontWeight = 'normal';
ax2.YAxis(1).Color = Color;
ax2.YAxis(1).LimitsMode = 'manual';
ax2.YAxis(1).TickLabelInterpreter = 'latex';
ax2.YAxis(1).TickLabelsMode = 'manual';
ax2.YAxis(1).TickValuesMode = 'manual';
% ax2.YAxis(1).ExponentMode = 'manual';
% ax2.YAxis(1).Exponent = 6;
ax2.YAxis(1).Limits = [0 2e6];
% ax2.YAxis(1).TickValues = [0 0.2e6 0.4e6 0.6e6 0.8e6 1.0e6 1.2e6 1.4e6 1.6e6 1.8e6 2.0e6];
ax2.YAxis(1).TickValues = [0 2e5 4e5 6e5 8e5 10e5 12e5 14e5 16e5 18e5 20e5];
% ax2.YAxis(1).YRuler.TickLabelFormat = '%.3f';
ax2.YAxis(1).TickLabels = {'0'; '2'; '4'; '6'; '8'; '10'; '12'; '14'; '16'; '18'; '20'};
ax2.YAxis(1).ExponentMode = 'manual';
ax2.YAxis(1).Exponent = 6;
ax2.YAxis(1).Label.Interpreter = 'latex';
% ax2.YAxis(1).Label.FontSize = FontSize; % geen effect
ax2.YAxis(1).Label.FontWeight = 'normal';
ax2.YAxis(1).Label.Color = Color;
ax2.YAxis(1).Label.LineWidth = LineWidth3;
ax2.YAxis(1).Label.String = '$P_{\mathrm{s,e}} \, \times 10^{5} \, [\mathrm{Pa}]$';

yyaxis right % or: yyaxis right; set(ax.YAxis(2).Label, 'Interpreter', 'latex', 'FontSize', 16)
ax2.YAxis(2).FontSize = FontSize3;%FontSize2;
ax2.YAxis(2).FontWeight = 'normal';
ax2.YAxis(2).Color = Color;
ax2.YAxis(2).LimitsMode = 'manual';
ax2.YAxis(2).TickLabelInterpreter = 'latex';
ax2.YAxis(2).TickLabelsMode = 'manual';
ax2.YAxis(2).TickValuesMode = 'manual';
ax2.YAxis(2).Limits = [0 2.5];
ax2.YAxis(2).TickLabelsMode = 'manual';
ax2.YAxis(2).TickValues = [0 0.5 1 1.5 2 2.5];
ax2.YAxis(2).TickLabels = {'0.0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'};
ax2.YAxis(2).Label.Interpreter = 'latex';
% ax2.YAxis(2).Label.FontSize = FontSize; % geen effect
ax2.YAxis(2).Label.FontWeight = 'normal';
ax2.YAxis(2).Label.Color = Color;
ax2.YAxis(2).Label.LineWidth = LineWidth3;
ax2.YAxis(2).Label.String = 'Mach $[-]$';

ax2.XLabel.Interpreter = 'latex';
ax2.FontSize = FontSize2;
% ax2.XLabel.FontSize = FontSize2;
ax2.XLabel.FontWeight = 'normal';
ax2.XLabel.Color = Color;
% ax2.XLabel.String = '$X \, [\mathrm{m}]$';
ax2.XLabel.String = '$x \, [\mathrm{m}]$';
ax2.XLimMode = 'manual';
ax2.XLim = [0 0.1];
% ax2.XLim = [0 0.09];
ax2.TickLabelInterpreter = 'latex';
ax2.XTickLabelMode = 'manual';
ax2.XTickMode = 'manual';
ax2.XTick = [0 0.02 0.04 0.06 0.08 0.10];
ax2.XTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};
% ax2.XTick = [0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09];
% ax2.XTickLabel = {'0'; '0.01'; '0.02'; '0.03'; '0.04'; '0.05'; '0.06'; '0.07'; '0.08';'0.09'};
% % % ax2.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
set(gcf,'Position',[488 342 1120 420]);
% ax2.XLabel.LineWidth = LineWidth; % what is this?
% set(gca,'FontSize',11) % because fontsize doesn't have any effect on results
set(gca,'FontSize',FontSize3) % alternative method
box on
hold on

% figure 3: Zoom inlet pressure
h3 = figure;
ax3 = gca;

ax3.FontSizeMode = 'manual';
% ax3.FontSize = FontSize2; % like this, or alterantive method below
ax3.FontWeight = 'normal';
ax3.TickLabelInterpreter = 'latex';
ax3.XScale = 'linear';
ax3.YScale = 'linear';

ax3.XLabel.Interpreter = 'latex';
% ax3.XLabel.FontSize = FontSize;
ax3.XLabel.FontWeight = 'normal';
ax3.XLabel.Color = Color;
% ax3.XLabel.String = '$X \, [\mathrm{m}]$';
ax3.XLabel.String = '$x \, [\mathrm{m}]$';
ax3.XLim = [0 0.025];
ax3.XLimMode = 'manual';
ax3.XTickLabelMode = 'manual';
ax3.XTickMode = 'manual';
ax3.XTick = [0 0.005 0.01 0.015 0.020 0.025];
ax3.XTickLabel = {'0'; '0.005'; '0.010'; '0.015'; '0.02'; '0.025'};

ax3.YLabel.Interpreter = 'latex';
% ax3.YLabel.FontSize = FontSize; %%
ax3.YLabel.FontWeight = 'normal';
ax3.YLabel.Color = Color;
ax3.YLabel.String = '$P_{\mathrm{s,e}} \, \times 10^{5} \, [\mathrm{Pa}]$';
% ax3.YLim = [1.805e6 1.825e6];
ax3.YLim = [1.810e6 1.825e6];
ax3.YLimMode = 'manual';
ax3.YTickLabelMode = 'manual';
ax3.YTickMode = 'manual';
% ax3.YTick = [1.805e6 1.810e6 1.815e6 1.820e6 1.825e6];
% ax3.YTickLabel = {'1.805'; '1.810'; '1.815'; '1.820'; '1.825'};
% ax3.YTick = [1.810e6 1.815e6 1.820e6 1.825e6];
% ax3.YTickLabel = {'1.810'; '1.815'; '1.820'; '1.825'};
ax3.YTick = [18.10e5 18.15e5 18.20e5 18.25e5];
ax3.YTickLabel = {'18.10'; '18.15'; '18.20'; '18.25'};
ax3.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',FontSize2) % alternative method
box on
hold on

% figure 4: Zoom outlet pressure
h4 = figure;
ax4 = gca;

ax4.FontSizeMode = 'manual';
% ax4.FontSize = FontSize2;
ax4.FontWeight = 'normal';
ax4.TickLabelInterpreter = 'latex';
ax4.XScale = 'linear';
ax4.YScale = 'linear';

ax4.XLabel.Interpreter = 'latex';
% ax4.XLabel.FontSize = FontSize;
ax4.XLabel.FontWeight = 'normal';
ax4.XLabel.Color = Color;
% ax4.XLabel.String = '$X \, [\mathrm{m}]$';
ax4.XLabel.String = '$x \, [\mathrm{m}]$';
% ax4.XLim = [0.085 0.100];
% ax4.XLim = [0.080 0.100];
ax4.XLim = [0.082 0.098];
ax4.XLimMode = 'manual';
ax4.XTickLabelMode = 'manual';
ax4.XTickMode = 'manual';
% ax4.XTick = [0.085 0.090 0.095 0.100];
% ax4.XTickLabel = {'0.085'; '0.090'; '0.095'; '0.100'};
% ax4.XTick = [0.080 0.085 0.090 0.095 0.100];
% ax4.XTickLabel = {'0.080'; '0.085'; '0.090'; '0.095'; '0.100'};
ax4.XTick = [0.085 0.090 0.095];
ax4.XTickLabel = {'0.085'; '0.090'; '0.095'};

ax4.YLabel.Interpreter = 'latex';
% ax4.YLabel.FontSize = FontSize; %%
ax4.YLabel.FontWeight = 'normal';
ax4.YLabel.Color = Color;
ax4.YLabel.String = '$P_{\mathrm{s,e}} \, \times 10^{5} \, [\mathrm{Pa}]$';
ax4.YLim = [0.188e6 0.198e6];
ax4.YLimMode = 'manual';
ax4.YTickLabelMode = 'manual';
ax4.YTickMode = 'manual';
% ax4.YTick = [0.188e6 0.190e6 0.192e6 0.194e6 0.196e6 0.198e6];
% ax4.YTickLabel = {'0.188'; '0.190'; '0.192'; '0.194'; '0.196'; '0.198'};
ax4.YTick = [1.88e5 1.90e5 1.92e5 1.94e5 1.96e5 1.98e5];
ax4.YTickLabel = {'1.88'; '1.90'; '1.92'; '1.94'; '1.96'; '1.98'};
ax4.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',FontSize2)
box on
hold on

% figure 5: Form factor
h5 = figure;
ax5 = gca;

ax5.FontSizeMode = 'manual';
ax5.FontSize = FontSize; %%
ax5.FontWeight = 'normal';
ax5.TickLabelInterpreter = 'latex';
ax5.XScale = 'linear';
ax5.YScale = 'linear';

ax5.XLabel.Interpreter = 'latex';
ax5.XLabel.FontSize = FontSize; %%
ax5.XLabel.FontWeight = 'normal';
ax5.XLabel.Color = Color;
ax5.XLabel.String = '$X \, [\mathrm{m}]$';
ax5.XLim = [0 0.1]; % change to 0.15?
ax5.XLimMode = 'manual';
ax5.XTickLabelMode = 'manual';
ax5.XTickMode = 'manual';
ax5.XTick = [0 0.02 0.04 0.06 0.08 0.10];
ax5.XTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};

ax5.YLabel.Interpreter = 'latex';
ax5.YLabel.FontSize = FontSize; %%
ax5.YLabel.FontWeight = 'normal';
ax5.YLabel.Color = Color;
ax5.YLabel.String = '$H \, [-]$';
ax5.YLim = [1 3];
ax5.YLimMode = 'manual';
ax5.YTickLabelMode = 'manual';
ax5.YTickMode = 'manual';
ax5.YTick = [1 1.5 2 2.5 3];
ax5.YTickLabel = {'1.0'; '1.5'; '2.0'; '2.5'; '3.0'};
ax5.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',FontSize)
box on
hold on

% figure 6: Displacement thickness (separate graph for velocity thickness?)
h6 = figure;
ax6 = gca;

ax6.FontSizeMode = 'manual';
ax6.FontSize = FontSize; %%
ax6.FontWeight = 'normal';
ax6.TickLabelInterpreter = 'latex';
ax6.XScale = 'linear';
ax6.YScale = 'linear';

ax6.XLabel.Interpreter = 'latex';
ax6.XLabel.FontSize = FontSize; %%
ax6.XLabel.FontWeight = 'normal';
ax6.XLabel.Color = Color;
ax6.XLabel.String = '$X \, [\mathrm{m}]$';
ax6.XLim = [0 0.1]; % change to 0.15?
ax6.XLimMode = 'manual';
ax6.XTickLabelMode = 'manual';
ax6.XTickMode = 'manual';
ax6.XTick = [0 0.02 0.04 0.06 0.08 0.10];
ax6.XTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};

ax6.YLabel.Interpreter = 'latex';
ax6.YLabel.FontSize = FontSize; %%
ax6.YLabel.FontWeight = 'normal';
ax6.YLabel.Color = Color;
ax6.YLabel.String = '$\delta^* \, \times 10^{-3} \, [\mathrm{m}]$';
ax6.YLim = [0 1.05e-4];
ax6.YLimMode = 'manual';
ax6.YTickLabelMode = 'manual';
ax6.YTickMode = 'manual';
ax6.YTick = [0 0.02e-3 0.04e-3 0.06e-3 0.08e-3 0.1e-3];
ax6.YTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};
ax6.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',FontSize)
box on
hold on

% figure 7: Velocity thickness
h7 = figure;
ax7 = gca;

ax7.FontSizeMode = 'manual';
ax7.FontSize = FontSize; %%
ax7.FontWeight = 'normal';
ax7.TickLabelInterpreter = 'latex';
ax7.XScale = 'linear';
ax7.YScale = 'linear';

ax7.XLabel.Interpreter = 'latex';
ax7.XLabel.FontSize = FontSize; %%
ax7.XLabel.FontWeight = 'normal';
ax7.XLabel.Color = Color;
ax7.XLabel.String = '$X \, [\mathrm{m}]$';
ax7.XLim = [0 0.1]; % change to 0.15?
ax7.XLimMode = 'manual';
ax7.XTickLabelMode = 'manual';
ax7.XTickMode = 'manual';
ax7.XTick = [0 0.02 0.04 0.06 0.08 0.10];
ax7.XTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};

ax7.YLabel.Interpreter = 'latex';
ax7.YLabel.FontSize = FontSize; %%
ax7.YLabel.FontWeight = 'normal';
ax7.YLabel.Color = Color;
ax7.YLabel.String = '$\delta \, \times 10^{-3} \, [\mathrm{m}]$';
ax7.YLim = [0 0.8e-3];
ax7.YLimMode = 'manual';
ax7.YTickLabelMode = 'manual';
ax7.YTickMode = 'manual';
% ax7.YTick = [0 5e-4 1e-3];
% ax7.YTickLabel = {'0'; '0.5'; '1'};
ax7.YTick = [0 0.2e-3 0.4e-3 0.6e-3 0.8e-3];
ax7.YTickLabel = {'0'; '0.20'; '0.40'; '0.60'; '0.80'};
ax7.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',FontSize)
box on
hold on

% figure 8: Loss coefficient
h8 = figure;
ax8 = gca;

ax8.FontSizeMode = 'manual';
ax8.FontSize = FontSize; %%
ax8.FontWeight = 'normal';
ax8.TickLabelInterpreter = 'latex';
ax8.XScale = 'linear';
ax8.YScale = 'linear';

ax8.XLabel.Interpreter = 'latex';
ax8.XLabel.FontSize = FontSize; %%
ax8.XLabel.FontWeight = 'normal';
ax8.XLabel.Color = Color;
ax8.XLabel.String = '$X \, [\mathrm{m}]$';
ax8.XLim = [0 0.1];
ax8.XLimMode = 'manual';
ax8.XTickLabelMode = 'manual';
ax8.XTickMode = 'manual';
ax8.XTick = [0 0.02 0.04 0.06 0.08 0.10];
ax8.XTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};

ax8.YLabel.Interpreter = 'latex';
ax8.YLabel.FontSize = FontSize; %%
ax8.YLabel.FontWeight = 'normal';
ax8.YLabel.Color = Color;
ax8.YLabel.String = '$C_{\mathrm{d}} \, \times 10^{-3} \, [-]$';
ax8.YLim = [0 3e-3];
ax8.YLimMode = 'manual';
ax8.YTickLabelMode = 'manual';
ax8.YTickMode = 'manual';
ax8.YTick = [0 0.5e-3 1e-3 1.5e-3 2e-3 2.5e-3 3e-3];
ax8.YTickLabel = {'0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'; '3.0'};
ax8.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',FontSize)
box on
hold on

% figure 9: Entropy generation rate
h9 = figure;
ax9 = gca;

ax9.FontSizeMode = 'manual';
ax9.FontSize = FontSize; %%
ax9.FontWeight = 'normal';
ax9.TickLabelInterpreter = 'latex';
ax9.XScale = 'linear';
ax9.YScale = 'linear';
  
ax9.XLabel.Interpreter = 'latex';
ax9.XLabel.FontSize = FontSize; %%
ax9.XLabel.FontWeight = 'normal';
ax9.XLabel.Color = Color;
ax9.XLabel.String = '$X \, [\mathrm{m}]$';
ax9.XLim = [0 0.1];
ax9.XLimMode = 'manual';
ax9.XTickLabelMode = 'manual';
ax9.XTickMode = 'manual';
ax9.XTick = [0 0.02 0.04 0.06 0.08 0.10];
ax9.XTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};
  
ax9.YLabel.Interpreter = 'latex';
ax9.YLabel.FontSize = FontSize; %%
ax9.YLabel.FontWeight = 'normal';
ax9.YLabel.Color = Color;
% ax9.YLabel.String = '$\dot{S}_{\mathrm{A}} \, [\mathrm{J}\mathrm{s}^{-1}\mathrm{K}^{-1}\mathrm{m}^{-2}]$';
ax9.YLabel.String = '$\dot{S}_{\mathrm{A}} \, [\mathrm{J}/\mathrm{s}\mathrm{K}\mathrm{m}^{2}]$';
ax9.YLim = [0 700];
ax9.YLimMode = 'manual';
ax9.YTickLabelMode = 'manual';
ax9.YTickMode = 'manual';
ax9.YTick = [0 100 200 300 400 500 600 700];
ax9.YTickLabel = {'0'; '100'; '200'; '300'; '400'; '500'; '600'; '700'};
ax9.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',FontSize)
box on
hold on

% figure 10: Cumulative entropy generated inside BL
h10 = figure;
ax10 = gca;

ax10.FontSizeMode = 'manual';
ax10.FontSize = FontSize; %%
ax10.FontWeight = 'normal';
ax10.TickLabelInterpreter = 'latex';
ax10.XScale = 'linear';
ax10.YScale = 'linear';

ax10.XLabel.Interpreter = 'latex';
ax10.XLabel.FontSize = FontSize; %%
ax10.XLabel.FontWeight = 'normal';
ax10.XLabel.Color = Color;
ax10.XLabel.String = '$X \, [\mathrm{m}]$';
ax10.XLim = [0 0.1];
ax10.XLimMode = 'manual';
ax10.XTickLabelMode = 'manual';
ax10.XTickMode = 'manual';
ax10.XTick = [0 0.02 0.04 0.06 0.08 0.10];
ax10.XTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};

ax10.YLabel.Interpreter = 'latex';
ax10.YLabel.FontSize = FontSize; %%
ax10.YLabel.FontWeight = 'normal';
ax10.YLabel.Color = Color;
% ax10.YLabel.String = '$\dot{S} \, [\mathrm{J}\mathrm{s}^{-1}\mathrm{K}^{-1}\mathrm{m}^{-1}]$';
ax10.YLabel.String = '$\dot{S} \, [\mathrm{J}/\mathrm{s}\mathrm{K}\mathrm{m}]$';
ax10.YLim = [0 35];
ax10.YLimMode = 'manual';
ax10.YTickLabelMode = 'manual';
ax10.YTickMode = 'manual';
ax10.YTick = [0 5 10 15 20 25 30 35];
ax10.YTickLabel = {'0'; '5'; '10'; '15'; '20'; '25'; '30'; '35'};
ax10.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',FontSize)
box on
hold on

% figure 11: Wall Turbulent Prandtl-number PrTw
h11 = figure;
ax11 = gca;

ax11.FontSizeMode = 'manual';
ax11.FontSize = FontSize; %%
ax11.FontWeight = 'normal';
ax11.TickLabelInterpreter = 'latex';
ax11.XScale = 'linear';
ax11.YScale = 'linear';
  
ax11.XLabel.Interpreter = 'latex';
ax11.XLabel.FontSize = FontSize; %%
ax11.XLabel.FontWeight = 'normal';
ax11.XLabel.Color = Color;
ax11.XLabel.String = '$X \, [\mathrm{m}]$';
ax11.XLim = [0 0.1];
ax11.XLimMode = 'manual';
ax11.XTickLabelMode = 'manual';
ax11.XTickMode = 'manual';
ax11.XTick = [0 0.02 0.04 0.06 0.08 0.10];
ax11.XTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};
  
ax11.YLabel.Interpreter = 'latex';
ax11.YLabel.FontSize = FontSize; %%
ax11.YLabel.FontWeight = 'normal';
ax11.YLabel.Color = Color;
ax11.YLabel.String = '$\mathrm{Pr}_{\mathrm{T}} \, [-], \, \, \mathrm{(RefProp)}$';
ax11.YLim = [1.24 1.29];
ax11.YLimMode = 'manual';
ax11.YTickLabelMode = 'manual';
ax11.YTickMode = 'manual';
ax11.YTick = [1.24 1.25 1.26 1.27 1.28 1.29];
ax11.YTickLabel = {'1.24'; '1.25'; '1.26'; '1.27'; '1.28'; '1.29'};
ax11.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',FontSize)
box on
hold on

%% Load Data and Simulation Results for RefProp laminar
% load('../Stored_Sim_Data/ORCHID_Nozzle_RefProp_PsE_Smoothed')
% laminar or turbulent does not matter for the Edge properties:
% load('../Stored_Sim_Data/ORCHID_Nozzle_RefProp_PsE_Smoothed_lam')
load('../Stored_Sim_Data/ORCHID_Nozzle_RefProp_PsE_Interpolated_lam') %%% NEW (31-03-2020)
% load('../Stored_Sim_Data/ORCHID_Nozzle_RefProp_PsE_Smoothed_turb')
load('../Data/ORCHID_Nozzle_Data_RefProp')
% load('../Data/ORCHID_PsE_adjusted_RefProp')
load('../Data/ORCHID_PsE_Interpolated_RefProp') %%% NEW (31-03-2020)
PsE_RP = transpose(table2array(ORCHID_Nozzle_Data_RefProp(1:249,4)));
MaE_RP = transpose(table2array(ORCHID_Nozzle_Data_RefProp(1:249,6)));
ZE_RP = EDG.ZE;
CpE_RP = EDG.CpE;
muE_RP = EDG.muE;
kE_RP = EDG.kE;
% Re_RP = BLC.Re_x./X; % unit Reynolds-number
UE_RP = EDG.UE;
rhoE_RP = EDG.rhoE;
% second table
delta_RP_lam = BLC.delta;
delta_ast_RP_lam = BLC.delta_ast;
delta_theta_RP_lam = BLC.theta;
H_RP_lam = BLC.H;
Re_x_RP = EDG.Re_x;
Re_theta_RP_lam = BLC.Re_theta;
Cf_RP_lam = BLC.Cf;
Cd_RP_lam = BLC.Cd;
r_RP_lam = BLC.Hrecovery;
for i = 1:length(SOL)
Cw_RP_lam(i) = FLP{i}.C(1);
cw_RP_lam(i) = SOL{i}.c(1);
end
% This, or later just load BLC.S_cum and plot it (included in OUTPUT-file
% on 16-04-2020)
S_cum_lam(1,1) = 0;
NS = 2;
S_cum_lam(1,NS) = S_cum_lam(NS - 1) + (EDG.rhoE(NS)*EDG.UE(NS)^3/EDG.TsE(NS)*BLC.Cd(NS))/2*(X(NS) - X(NS - 1));
for NS = 3:length(X)
    S_cum_lam(NS) = S_cum_lam(NS - 1) + (EDG.rhoE(NS)*EDG.UE(NS)^3/EDG.TsE(NS)*BLC.Cd(NS) + EDG.rhoE(NS - 1)*EDG.UE(NS - 1)^3/EDG.TsE(NS - 1)*BLC.Cd(NS - 1))/2*(X(NS) - X(NS - 1));
end  % units of entropy?

%% Plot RefProp laminar
figure(h1) % ZPR
PrE_RP = EDG.muE.*EDG.CpE./EDG.kE;
plot(X,EDG.ZE,'Color',Color3,'LineStyle',LineStyle{1},'LineWidth',LineWidth2,'Marker',Marker{1},'DisplayName','$Z_{\mathrm{e}}$ RefProp');
plot(X,PrE_RP,'Color',Color3,'LineStyle',LineStyle{2},'LineWidth',LineWidth2,'Marker',Marker{2},'DisplayName','$\mathrm{Pr}_{\mathrm{e}}$ RefProp');
plot(X,BLC.Hrecovery,'Color',Color3,'LineStyle',LineStyle{3},'LineWidth',LineWidth2,'Marker',Marker{3},'DisplayName','$r_{\mathrm{Lam}}$ RefProp');

figure(h2) % Pressure - Mach
yyaxis left
% plot(X,PsE_RP,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','$P_{\mathrm{s,e}}$ SU2 RefProp');
plot(INP.x,PsE_RP,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','$P_{\mathrm{s,e}}$ SU2 RefProp');

yyaxis right
% plot(X,MaE_RP,'Color',Color,'LineStyle',LineStyle{2},'LineWidth',LineWidth3,'Marker',Marker{2},'DisplayName','Mach SU2 RefProp');
plot(INP.x,MaE_RP,'Color',Color,'LineStyle',LineStyle{2},'LineWidth',LineWidth3,'Marker',Marker{2},'DisplayName','Mach SU2 RefProp');

figure(h3) % Pressure inlet
% % plot(X,ORCHID_PsE_Interpolated_RefProp,'Color',Color3,'LineStyle',LineStyle{2},'LineWidth',LineWidth2,'Marker',Marker{1},'DisplayName','interpolated pressure RefProp');
% % plot(X,PsE_RP,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','$P_{\mathrm{s,e}}$ SU2 RefProp');
plot(INP.x,ORCHID_PsE_Interpolated_RefProp,'Color',Color3,'LineStyle',LineStyle{2},'LineWidth',LineWidth2,'Marker',Marker{1},'DisplayName','interpolated pressure RefProp');
plot(INP.x,PsE_RP,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','$P_{\mathrm{s,e}}$ SU2 RefProp');

% plot(X,ORCHID_PsE_adjusted_RefProp,'Color','b','LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName','interpolated pressure RefProp');
% plot(X,ORCHID_PsE_Interpolated_RefProp,'Color','b','LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','interpolated pressure RefProp');

figure(h4) % Pressure outlet
% % plot(X,ORCHID_PsE_Interpolated_RefProp,'Color',Color3,'LineStyle',LineStyle{2},'LineWidth',LineWidth2,'Marker',Marker{1},'DisplayName','interpolated pressure RefProp');
% % plot(X,PsE_RP,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','SU2 RefProp');
plot(INP.x,ORCHID_PsE_Interpolated_RefProp,'Color',Color3,'LineStyle',LineStyle{2},'LineWidth',LineWidth2,'Marker',Marker{1},'DisplayName','interpolated pressure RefProp');
plot(INP.x,PsE_RP,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','SU2 RefProp');

% plot(X,ORCHID_PsE_adjusted_RefProp,'Color','b','LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName','interpolated pressure RefProp');
% plot(X,ORCHID_PsE_Interpolated_RefProp,'Color','b','LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','interpolated pressure RefProp');

figure(h5) % Shape factor
plot(X,BLC.H,'Color',Color3,'LineStyle',LineStyle{2},'LineWidth',LineWidth2,'Marker',Marker{2},'DisplayName','laminar RefProp');

figure(h6) % displacement thickness
% % plot(X,BLC.delta_ast,'Color',Color3,'LineStyle',LineStyle{2},'LineWidth',LineWidth2,'Marker',Marker{2},'DisplayName','laminar RefProp');
plot(X,BLC.delta_ast,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth2,'Marker',Marker{2},'DisplayName','laminar RefProp');

figure(h7) % velocity thickness
plot(X,BLC.delta,'Color',Color3,'LineStyle',LineStyle{2},'LineWidth',LineWidth2,'Marker',Marker{2},'DisplayName','laminar RefProp');

% Loss Coefficient Cd (Denton, 1993)
% figure(h8) % see below
% plot(X,BLC.Cd,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','laminar RefProp')

%% Load Data and Simulation Results for StanMix laminar
% load('../Stored_Sim_Data/ORCHID_Nozzle_StanMix_PsE_Smoothed')
% load('../Stored_Sim_Data/ORCHID_Nozzle_StanMix_PsE_Smoothed_lam')
load('../Stored_Sim_Data/ORCHID_Nozzle_StanMix_PsE_Interpolated_lam')
load('../Data/ORCHID_Nozzle_Data_StanMix')
% load('../Data/ORCHID_PsE_adjusted_StanMix')
% load('../Data/ORCHID_PsE_Interpolated_StanMix') % does not exist, same
% RefProp interpolated value is used.
PsE_SM = transpose(table2array(ORCHID_Nozzle_Data_StanMix(1:249,4)));
MaE_SM = transpose(table2array(ORCHID_Nozzle_Data_StanMix(1:249,6)));
ZE_SM = EDG.ZE;
CpE_SM = EDG.CpE;
muE_SM = EDG.muE;
kE_SM = EDG.kE;
Re_SM = BLC.Re_x./X; % unit Reynolds-number
% second table
% delta_SM = BLC.delta;
% delta_ast_SM = BLC.delta_ast;
% delta_theta_SM = BLC.theta;
UE_SM = EDG.UE;
rhoE_SM = EDG.rhoE;
delta_SM_lam = BLC.delta;
delta_ast_SM_lam = BLC.delta_ast;
delta_theta_SM_lam = BLC.theta;
H_SM_lam = BLC.H;
% Re_x_SM_lam = BLC.Re_x;
Re_x_SM = EDG.Re_x;
Re_theta_SM_lam = BLC.Re_theta;
Cf_SM_lam = BLC.Cf;
Cd_SM_lam = BLC.Cd;
r_SM_lam = BLC.Hrecovery;
for i = 1:length(SOL)
Cw_SM_lam(i) = FLP{i}.C(1);
cw_SM_lam(i) = SOL{i}.c(1);
end
% laminar or turbulent does not matter for the Edge properties:
% load('../Stored_Sim_Data/ORCHID_Nozzle_StanMix_PsE_Smoothed_lam')
% load('../Stored_Sim_Data/ORCHID_Nozzle_StanMix_PsE_Smoothed_turb')

%% Plot StanMix laminar

figure(h1)
PrE_SM = EDG.muE.*EDG.CpE./EDG.kE;
plot(X,EDG.ZE,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth,'Marker',Marker{4},'DisplayName','$Z_{\mathrm{e}}$ StanMix');
plot(X,PrE_SM,'Color',Color,'LineStyle',LineStyle{2},'LineWidth',LineWidth,'Marker',Marker{5},'DisplayName','$\mathrm{Pr}_{\mathrm{e}}$ StanMix');
plot(X,BLC.Hrecovery,'Color',Color,'LineStyle',LineStyle{3},'LineWidth',LineWidth,'Marker',Marker{6},'DisplayName','$r_{\mathrm{Lam}}$ StanMix');

% figure(h2)
% yyaxis left
% plot(X,PsE_SM,'Color',Color,'LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName','$P_{\mathrm{s,e}}$ SU2 StanMix');
% 
% yyaxis right
% plot(X,MaE_SM,'Color',Color,'LineStyle',LineStyle{2},'Marker',Marker{2},'DisplayName','Mach SU2 StanMix');

% Barely visible:
% figure(h3) % Pressure inlet
% plot(X,PsE_SM,'Color',Color,'LineStyle',LineStyle{3},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','SU2 StanMix');
% % plot(X,ORCHID_PsE_adjusted_StanMix,'Color',Color,'LineStyle',LineStyle{4},'Marker',Marker{4},'DisplayName','interpolated StanMix');

figure(h4) % Pressure outlet
% plot(X,PsE_SM,'Color',Color,'LineStyle',LineStyle{3},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','SU2 StanMix');
plot(INP.x,PsE_SM,'Color',Color,'LineStyle',LineStyle{3},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','SU2 StanMix');

% plot(X,ORCHID_PsE_adjusted_StanMix,'Color',Color,'LineStyle',LineStyle{4},'Marker',Marker{4},'DisplayName','interpolated StanMix');

figure(h5)
plot(X,BLC.H,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','laminar StanMix');

figure(h6)
% % plot(X,BLC.delta_ast,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','laminar StanMix');

figure(h7)
plot(X,BLC.delta,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','laminar StanMix');

%% Load RefProp turbulent
% % % load('../Stored_Sim_Data/ORCHID_Nozzle_RefProp_PsE_Smoothed_turb') %% Why this one?
load('../Stored_Sim_Data/ORCHID_Nozzle_RefProp_PsE_Interpolated_turb') % NEW (18-05-2020)

% second table
delta_RP_turb = BLC.delta;
delta_ast_RP_turb = BLC.delta_ast;
delta_theta_RP_turb = BLC.theta;
H_RP_turb = BLC.H;
% Re_x_SM_lam = BLC.Re_x;
% Re_x_SM = EDG.Re_x;
Re_theta_RP_turb = BLC.Re_theta;
Cf_RP_turb = BLC.Cf;
Cd_RP_turb = BLC.Cd;
r_RP_turb = BLC.Hrecovery;
for i = 1:length(SOL)
Cw_RP_turb(i) = FLP{i}.C(1);
cw_RP_turb(i) = SOL{i}.c(1);
PrTw_RP(i) = FLP{i}.PrT(1);
end
% This, or later just load BLC.S_cum and plot it (included in OUTPUT-file
% on 16-04-2020)
S_cum_turb(1,1) = 0;
NS = 2;
S_cum_turb(1,NS) = S_cum_turb(NS - 1) + (EDG.rhoE(NS)*EDG.UE(NS)^3/EDG.TsE(NS)*BLC.Cd(NS))/2*(X(NS) - X(NS - 1));
for NS = 3:length(X)
    S_cum_turb(NS) = S_cum_turb(NS - 1) + (EDG.rhoE(NS)*EDG.UE(NS)^3/EDG.TsE(NS)*BLC.Cd(NS) + EDG.rhoE(NS - 1)*EDG.UE(NS - 1)^3/EDG.TsE(NS - 1)*BLC.Cd(NS - 1))/2*(X(NS) - X(NS - 1));
end

%% Plot RefProp turbulent

figure(h1)
plot(X,BLC.Hrecovery,'Color',Color3,'LineStyle',LineStyle{4},'LineWidth',LineWidth2,'Marker',Marker{4},'DisplayName','$r_{\mathrm{Turb}}$ RefProp');

figure(h5)
plot(X,BLC.H,'Color',Color3,'LineStyle',LineStyle{2},'LineWidth',LineWidth2,'Marker',Marker{2},'DisplayName','turbulent RefProp');

figure(h6)
% % plot(X,BLC.delta_ast,'Color',Color3,'LineStyle',LineStyle{2},'LineWidth',LineWidth2,'Marker',Marker{2},'DisplayName','turbulent RefProp');
plot(X,BLC.delta_ast,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth2,'Marker',Marker{2},'DisplayName','turbulent RefProp');

figure(h7)
plot(X,BLC.delta,'Color',Color3,'LineStyle',LineStyle{2},'LineWidth',LineWidth2,'Marker',Marker{2},'DisplayName','turbulent RefProp');

% Loss Coefficient Cd (Denton, 1993)
% figure(h8) % see below
% plot(X,BLC.Cd,'DisplayName','RefProp Turb') % OLD (13-05-2020)
% plot(X,BLC.Cd,'Color',Color3,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','turbulent RefProp')

%% Load StanMix turbulent
% load('../Stored_Sim_Data/ORCHID_Nozzle_StanMix_PsE_Smoothed_turb') %% Why this one?
load('../Stored_Sim_Data/ORCHID_Nozzle_StanMix_PsE_Interpolated_turb') % NEW (18-05-2020)

% second table
delta_SM_turb = BLC.delta;
delta_ast_SM_turb = BLC.delta_ast;
delta_theta_SM_turb = BLC.theta;
H_SM_turb = BLC.H;
% Re_x_SM_lam = BLC.Re_x;
% Re_x_SM = EDG.Re_x;
Re_theta_SM_turb = BLC.Re_theta;
Cf_SM_turb = BLC.Cf;
Cd_SM_turb = BLC.Cd;
r_SM_turb = BLC.Hrecovery;
for i = 1:length(SOL)
Cw_SM_turb(i) = FLP{i}.C(1);
cw_SM_turb(i) = SOL{i}.c(1);
end

%% Plot StanMix turbulent

figure(h1)
plot(X,BLC.Hrecovery,'Color',Color,'LineStyle',LineStyle{4},'LineWidth',LineWidth,'Marker',Marker{2},'DisplayName','$r_{\mathrm{Turb}}$ StanMix');

figure(h5)
plot(X,BLC.H,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','turbulent StanMix');

figure(h6)
% % plot(X,BLC.delta_ast,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','turbulent StanMix');

figure(h7)
plot(X,BLC.delta,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','turbulent StanMix');

%% Plot remaining stored properties (S_a, S, Pr_T,w)
cd .. % for use of GRADNT-file to calculate the derivatives

figure(h8)  % Loss Coefficient Cd (Denton, 1993)
plot(X,Cd_RP_turb,'Color',Color3,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','turbulent RefProp')
plot(X,Cd_RP_lam,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','laminar RefProp')

figure(h9)  % Entropy generation rate per surface area
plot(X(1:length(S_cum_turb)),GRADNT(S_cum_turb,X(1:length(S_cum_turb))),'Color',Color3,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','turb. RefProp')
plot(X(1:length(S_cum_lam)),GRADNT(S_cum_lam,X(1:length(S_cum_lam))),'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','laminar RefProp')

figure(h10) % Cumulative entropy generated inside BL (per meter depth (z-direction))
plot(X(1:length(S_cum_turb)),S_cum_turb,'Color',Color3,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','turbulent RefProp')
plot(X(1:length(S_cum_lam)),S_cum_lam,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'DisplayName','laminar RefProp')

figure(h11) % Turbulent Prandtl-number at the wall
plot(X,PrTw_RP,'Color',Color3,'LineStyle',LineStyle{1},'LineWidth',LineWidth3,'Marker',Marker{1},'HandleVisibility','off')%DisplayName','turbulent RefProp')

%% Add legends
figure(h1)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize1, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northeast');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')
% legend('Pr RefProp','Pr StanMix','Z StanMix','Z RefProp','Enthalpy Recovery factor StanMix')
% legend('Pr RefProp','Pr StanMix','Z StanMix','Z RefProp','Enthalpy Recovery factor StanMix','Enthalpy Recovery factor RefProp')

figure(h2)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize2, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northeast');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

figure(h3)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize2, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'southwest');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

figure(h4)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize2, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northeast');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

figure(h5)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'north');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

figure(h6)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

figure(h7)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

figure(h8)  % Loss Coefficient Cd (Denton, 1993)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

figure(h9)  % Entropy generation rate per surface area
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

figure(h10) % Cumulative entropy generated inside BL (per meter depth (z-direction))
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

% figure(h11) % Turbulent Prandtl-number at the wall %%% no legend to show
% set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest');%'Location', 'best'); %%%  'Position', LgndPos2); %
% legend('show')

%% Store plots %%% moved downstream
% keyboard %%% uncomment save-commands
% NewFolder = cd;
% saveas(h1,'Thesis_Report_Figures/ORCHID_ZPR','epsc')
% saveas(h2,'Thesis_Report_Figures/ORCHID_PsE_Ma','epsc')
% saveas(h3,'Thesis_Report_Figures/ORCHID_PsE_inlet','epsc')
% saveas(h4,'Thesis_Report_Figures/ORCHID_PsE_outlet','epsc')
% saveas(h5,'Thesis_Report_Figures/ORCHID_H','epsc')
% saveas(h6,'Thesis_Report_Figures/ORCHID_delta_ast','epsc')
% saveas(h7,'Thesis_Report_Figures/ORCHID_delta','epsc')
% saveas(h8,'Thesis_Report_Figures/ORCHID_Cd','epsc')
% OldFolder = cd(NewFolder);

% saveas(h1,'Thesis_Report_Figures/ORCHID_ZPR','epsc')

% fprintf('Figure(s) generated and saved successfully\n')

% Fluid property plots
h12 = figure;
hold on
plot(X,CpE_RP)
plot(X,CpE_SM)
title('Constant pressure specific heat (Cp) along BL edge')
xlabel('X [m]')
ylabel('CpE [J/kg/K]')
legend('RefProp','StanMix')

h13 = figure;
hold on
plot(X,muE_RP)
plot(X,muE_SM)
title('Dynamic viscosity')
xlabel('X [m]')
ylabel('muE [Pa s]')
legend('RefProp','StanMix')

h14 = figure;
hold on
plot(X,kE_RP)
plot(X,kE_SM)
title('Thermal conductivity')
xlabel('X [m]')
ylabel('kE [W/m/K]')
legend('RefProp','StanMix')

h15 = figure;
hold on
plot(X,PrE_RP)
plot(X,PrE_SM)
title('Prandtl-number')
xlabel('X [m]')
ylabel('PrE [-]')
legend('RefProp','StanMix')

% bonus-plot 1: compressibility factor (Zucrow)
Fc = (INP.PtI-EDG.PsE)./(0.5*EDG.rhoE.*EDG.UE.^2);
figure
plot(X,Fc)
title('Compressibility factor (Zucrow): F_c = p_0 - p / 1/2 \rho v^2')
xlabel('X [m]')
ylabel('F_c [-]')

% bonus-plot 2: efficiency? of pressure converted to dynamic pressure?
% (inverse compressibility factor (Zucrow))
figure
plot(X,1./Fc)
title('1/2 \rho v^2/(p_0 - p)')
xlabel('X [m]')
ylabel('F_c [-]')

figure
hold on
plot(INP.x,INP.y)
title('ORCHID Nozzle geometry')
axb = gca;
axb.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

figure % Turbulent Prandtl-number at wall
hold on
plot(X,PrTw_RP)
title('PrT wall')
xlabel('X [m]')
ylabel('PrT_{w} [-]')

figure % units of entropy?
hold on
plot(X(1:length(S_cum_lam)),S_cum_lam,'DisplayName','RefProp Lam')
plot(X(1:length(S_cum_turb)),S_cum_turb,'DisplayName','RefProp Turb')
title('Cumulative entropy generated inside the BL')
xlabel('X [m]')
ylabel('S [J/s/K/m]')
legend('show')

% cd ..

figure % units of entropy?
hold on
plot(X(1:length(S_cum_lam)),GRADNT(S_cum_lam,X(1:length(S_cum_lam))),'DisplayName','RefProp Lam')
plot(X(1:length(S_cum_turb)),GRADNT(S_cum_turb,X(1:length(S_cum_turb))),'DisplayName','RefProp Turb')
title('Entropy generation rate along the BL')
xlabel('X [m]')
ylabel('S [J/s/K/m^2]')
legend('show')

%% Store plots %%% moved here from above
% keyboard %%% uncomment save-commands
% NewFolder = cd;

saveas(h1,'Thesis_Report_Figures/ORCHID_ZPR','epsc')
saveas(h2,'Thesis_Report_Figures/ORCHID_PsE_Ma','epsc')
saveas(h3,'Thesis_Report_Figures/ORCHID_PsE_inlet','epsc')
saveas(h4,'Thesis_Report_Figures/ORCHID_PsE_outlet','epsc')
saveas(h5,'Thesis_Report_Figures/ORCHID_H','epsc')
saveas(h6,'Thesis_Report_Figures/ORCHID_delta_ast','epsc')
saveas(h7,'Thesis_Report_Figures/ORCHID_delta','epsc')
saveas(h8,'Thesis_Report_Figures/ORCHID_Cd','epsc')
saveas(h9,'Thesis_Report_Figures/ORCHID_SA','epsc')
saveas(h10,'Thesis_Report_Figures/ORCHID_S','epsc')
saveas(h11,'Thesis_Report_Figures/ORCHID_PrT','epsc')
% and not formatted fluid properties for Appendix
saveas(h12,'Thesis_Report_Figures/ORCHID_Cp','epsc')
saveas(h13,'Thesis_Report_Figures/ORCHID_mu','epsc')
saveas(h14,'Thesis_Report_Figures/ORCHID_k','epsc')
saveas(h15,'Thesis_Report_Figures/ORCHID_Pr','epsc')
% OldFolder = cd(NewFolder);

fprintf('Figure(s) generated and saved successfully\n')

%% Tables - Write data to Tex-file (in Table format)

%% Boundary Layer Edge Comparison
% Comparison of simulation results of the ORCHID nozzle with two different thermophysical models: RefProp and StanMix
% keyboard

NewFolder = cd;
cd('Thesis_Report_Tables') %%% new, added (23-01-2021)fid = fopen('Comparison_ORCHID_Sim_Th_Models.tex','w');
fid = fopen('Comparison_ORCHID_Sim_Th_Models.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Comparison of simulation results of the ORCHID nozzle with two different thermophysical models: RefProp and StanMix. Where StanMix is considered as a reference.\n');
% fprintf(fid, '\\begin{landscape}\n');
fprintf(fid, '%%\\begin{table}\n');
fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
% fprintf(fid, '    \\caption{Comparison with tabulated data of laminar compressible similar flows with constant nonzero pressure gradients and heat transfer for calorically perfect ideal gas with $C = 1$ (constant) and $\\mathrm{Pr} = 1$ taken from Rogers\\cite{rogers1992laminar} table C-25. Transformed uniform (vertical) grid spacing (in Illingworth-Levy transformed y-coordinate) of $\\mathrm{d} \\eta = \\sqrt{\\frac{m_2 + 1}{2C}} %1.4f$ and height of $\\eta_{\\mathrm{e}} = \\sqrt{\\frac{m_2 + 1}{2C}} %1.1f$. Note that separation occurred when the table entry shows ''sep''.}\n',GRD.Deta(1),GRD.etaE);
fprintf(fid, '    \\caption{Comparison of ORCHID simulations with two different thermodynamic models: boundary layer edge properties. Results calculated using StanMix and RefProp. StanMix utilizes the iPRSV, whereas RefProp utilizes unknow industrial standards? StanMix is used here as starting point for convenience, since the fluid property retrieval is much faster, whereas RefProp is thought to be more accurate. Note also that edge properties are the same for laminar or turbulent flow. Static pressure obtained from SU2 CFD simulation second-order convergence is used as input.}\n');
fprintf(fid, '    \\label{tab:cosep}\n'); % Comparison of ORCHID Simulations edge properties
fprintf(fid, '    \\begin{tabular}{l S[table-format=-2.3] S[table-format=-2.3] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=-2.2] S[table-format=3.3] S[table-format=3.3] S[table-format=-2.3] S[table-format=-2.3] S[table-format=-2.3] S[table-format=-2.3] S[table-format=-2.3] S[table-format=-2.3] S[table-format=-2.3] S[table-format=-2.3]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{Boundary}                  &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$P_{\\mathrm{s,e}}$}                  &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$\\mathrm{Ma}_{\\mathrm{e}}$}         &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$Z_{\\mathrm{e}}$}                    &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$\\mathrm{Pr}_{\\mathrm{e}}$}         &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$C_{\\mathrm{p,e}}$}                  &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$\\mu_{\\mathrm{e}}$}                 &\n');        
fprintf(fid, '        \\multicolumn{2}{c}{$k_{\\mathrm{e}}$}                    &\n');        
fprintf(fid, '        \\multicolumn{2}{c}{$\\mathrm{Re}_{\\mathrm{e}}$}      \\\\\n');
fprintf(fid, '        \\multicolumn{1}{c}{layer edge}                  &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$\\times \\, 10^5 \\,$ [Pa]}                     &\n');
fprintf(fid, '        \\multicolumn{2}{c}{[-]}                                &\n');
fprintf(fid, '        \\multicolumn{2}{c}{[-]}                    &\n');
fprintf(fid, '        \\multicolumn{2}{c}{[-]}         &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$\\times \\, 10^{3} \\,$ [J/kg/K]}         &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$\\times \\, 10^{-5} \\,$ [kg/m/s]}                  &\n');       
fprintf(fid, '        \\multicolumn{2}{c}{$\\times \\, 10^{-2} \\,$ [W/m/K]}                 &\n');        
fprintf(fid, '        \\multicolumn{2}{c}{$\\times \\, 10^7 \\,$ [1/m]}                 \\\\\n'); % 10^8
fprintf(fid, '        \\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9} \\cmidrule(lr){10-11} \\cmidrule(lr){12-13} \\cmidrule(lr){14-15} \\cmidrule(lr){16-17} \\cmidrule(lr){18-19}\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{properties}                  &\n');
fprintf(fid, '        \\multicolumn{1}{c}{in}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{out}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{in}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{out}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{in}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{out}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{in}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{out}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{in}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{out}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{in}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{out}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{in}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{out}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{in}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{out}                                      \\\\\n');
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');

fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{RefProp}}            &    %2.3f       &  %2.3f       &   %1.1f       & %1.1f       &   %1.2f          &  %1.2f   & %1.2f          &  %1.2f   &  %3.3f &  %3.3f  &  %1.3f  &  %1.3f  &  %1.3f  &  %1.3f &  %2.0f  &  %2.0f  \\\\\n',PsE_RP(1)/10^5,PsE_RP(end)/10^5,MaE_RP(1),MaE_RP(end),ZE_RP(1),ZE_RP(end),PrE_RP(1),PrE_RP(end),CpE_RP(1)*10^-3,CpE_RP(end)*10^-3,muE_RP(1)*10^5,muE_RP(end)*10^5,kE_RP(1)*10^2,kE_RP(end)*10^2,UE_RP(1)*rhoE_RP(1)/muE_RP(1)/10^7,UE_RP(end)*rhoE_RP(end)/muE_RP(end)/10^7);
fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{StanMix}}            &    %2.3f       &  %2.3f       &   %1.1f       & %1.1f       &   %1.2f          &  %1.2f   & %1.2f          &  %1.2f   &  %3.3f &  %3.3f  &  %1.3f  &  %1.3f  &  %1.3f  &  %1.3f &  %2.0f  &  %2.0f  \\\\\n',PsE_SM(1)/10^5,PsE_SM(end)/10^5,MaE_SM(1),MaE_SM(end),ZE_SM(1),ZE_SM(end),PrE_SM(1),PrE_SM(end),CpE_SM(1)*10^-3,CpE_SM(end)*10^-3,muE_SM(1)*10^5,muE_SM(end)*10^5,kE_SM(1)*10^2,kE_SM(end)*10^2,UE_SM(1)*rhoE_SM(1)/muE_SM(1)/10^7,UE_SM(end)*rhoE_SM(end)/muE_SM(end)/10^7);
fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{Deviation [\\%%]}}   &    %2.3f       &  %2.2f       &   %2.2f       & %2.2f       &   %2.2f          &  %2.2f   & %2.2f          &  %2.2f   &  %2.2f &  %2.2f  &  %2.2f  &  %2.2f  &  %2.2f  &  %2.2f &  %2.1f  &  %2.1f  \\\\\n',(PsE_SM(1)-PsE_RP(1))/PsE_RP(1)*100,(PsE_SM(end)-PsE_RP(end))/PsE_RP(end)*100,(MaE_SM(1)-MaE_RP(1))/MaE_RP(1)*100,(MaE_SM(end)-MaE_RP(end))/MaE_RP(end)*100,(ZE_SM(1)-ZE_RP(1))/ZE_RP(1)*100,(ZE_SM(end)-ZE_RP(end))/ZE_RP(end)*100,(PrE_SM(1)-PrE_RP(1))/PrE_RP(1)*100,(PrE_SM(end)-PrE_RP(end))/PrE_RP(end)*100,(CpE_SM(1)-CpE_RP(1))/CpE_RP(1)*100,(CpE_SM(end)-CpE_RP(end))/CpE_RP(end)*100,(muE_SM(1)-muE_RP(1))/muE_RP(1)*100,(muE_SM(end)-muE_RP(end))/muE_RP(end)*100,(kE_SM(1)-kE_RP(1))/kE_RP(1)*100,(kE_SM(end)-kE_RP(end))/kE_RP(end)*100,(UE_SM(1)*rhoE_SM(1)/muE_SM(1)-UE_RP(1)*rhoE_RP(1)/muE_RP(1))/(UE_RP(1)*rhoE_RP(1)/muE_RP(1))*100,(UE_SM(end)*rhoE_SM(end)/muE_SM(end)-UE_RP(end)*rhoE_RP(end)/muE_RP(end))/(UE_RP(end)*rhoE_RP(end)/muE_RP(end))*100);

fprintf(fid, '        \\bottomrule\n\n');
fprintf(fid, '    \\end{tabular}\n');
% fprintf(fid, '    \\begin{tablenotes}\n');
% fprintf(fid, '        \\item[*] The solution does not converge for $m_2 = \\infty$, and therefore $m_2 = 1000$ is taken as approach.\n');
% fprintf(fid, '    \\end{tablenotes}\n');
fprintf(fid, '\\end{threeparttable}\n');
% fprintf(fid, '%%\\end{sidewaystable}\n');
% fprintf(fid, '%%\\end{table\n');
% fprintf(fid, '\\end{landscape}\n');
fclose(fid);

%% Boundary Layer Characteristics Laminar Flow
% keyboard

fid = fopen('Comparison_ORCHID_Sim_Results_lam.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Comparison of laminar simulation results of the ORCHID nozzle with two different thermophysical models: RefProp and StanMix. Where StanMix is considered as a reference.\n');
% fprintf(fid, '\\begin{landscape}\n');
fprintf(fid, '%%\\begin{table}\n');
fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Comparison of ORCHID laminar simulations with two different thermodynamic models: boundary layer edge properties. Results calculated using StanMix and RefProp. StanMix utilizes the iPRSV, whereas RefProp utilizes unknow industrial standards? StanMix is used here as starting point for convenience, since the fluid property retrieval is much faster, whereas RefProp is thought to be more accurate. Note also that edge properties are the same for laminar or turbulent flow. Static pressure obtained from SU2 CFD simulation second-order convergence is used as input.}\n');
fprintf(fid, '    \\label{tab:coslam}\n'); % Comparison of ORCHID Simulations BL Results
fprintf(fid, '    \\begin{tabular}{l S[table-format=1.3] S[table-format=1.3] S[table-format=1.3] S[table-format=1.2] S[table-format=1.1] S[table-format=1.1] S[table-format=1.2] S[table-format=1.2] S[table-format=1.2] S[table-format=1.2] S[table-format=1.2] S[table-format=1.2]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{Laminar flow}                  &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\delta$}                  &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\delta^*$}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\theta$}                    &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$H$}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{Re}_{\\mathrm{e}}$}                  &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{Re}_x$}                  &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{Re}_{\\theta}$}                &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$C_{\\mathrm{f}}$}                    &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$C_{\\mathrm{d}}$}                    &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$r$}      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$c_{\\mathrm{w}}$}                    &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$C_{\\mathrm{w}}$}                    \\\\\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                  &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\,$ [m]}                     &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\,$ [m]}                                &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\,$ [m]}                    &\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{7} \\,$ [1/m]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{7} \\,$ [-]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{3} \\,$ [-]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\,$ [-]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\,$ [-]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}                      \\\\\n');
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
% S[table-format=1.3] S[table-format=1.3] S[table-format=1.3]
% S[table-format=1.2] S[table-format=1.1] S[table-format=1.2] S[table-format=1.2] S[table-format=1.2] S[table-format=1.2] S[table-format=1.2] delta_RP_lam  delta_ast_RP_lam  delta_theta_RP_lam  H_RP_lam  Re_x_RP  Re_theta_RP_lam  Cf_RP_lam  r_RP_lam  cw_RP_lam  Cw_RP_lam
fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{RefProp}}            &    %1.3f       &  %1.3f       &   %1.3f       & %1.2f       &   %2.0f   &   %1.1f          &  %1.2f   & %1.3f     & %1.3f     &  %1.3f   &  %1.3f &  %1.3f  \\\\\n',delta_RP_lam(end)*10^3,delta_ast_RP_lam(end)*10^3,delta_theta_RP_lam(end)*10^3,H_RP_lam(end),Re_x_RP(end)/X(end)*10^-7,Re_x_RP(end)*10^-7,Re_theta_RP_lam(end)*10^-3,Cf_RP_lam(end)*10^3,Cd_RP_lam(end)*10^3,r_RP_lam(end),cw_RP_lam(end),Cw_RP_lam(end));
fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{StanMix}}            &    %1.3f       &  %1.3f       &   %1.3f       & %1.2f       &   %2.0f   &   %1.1f          &  %1.2f   & %1.3f     & %1.3f     &  %1.3f   &  %1.3f &  %1.3f  \\\\\n',delta_SM_lam(end)*10^3,delta_ast_SM_lam(end)*10^3,delta_theta_SM_lam(end)*10^3,H_SM_lam(end),Re_x_SM(end)/X(end)*10^-7,Re_x_SM(end)*10^-7,Re_theta_SM_lam(end)*10^-3,Cf_SM_lam(end)*10^3,Cd_SM_lam(end)*10^3,r_SM_lam(end),cw_SM_lam(end),Cw_SM_lam(end));
fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{Deviation [\\%%]}}   &    %2.2f       &  %2.2f       &   %2.2f       & %2.2f       &   %2.1f   &   %2.1f          &  %2.1f   & %2.2f     & %2.2f     &  %2.2f   &  %2.2f &  %2.2f  \\\\\n',(delta_SM_lam(end) - delta_RP_lam(end))/delta_RP_lam(end)*100,(delta_ast_SM_lam(end) - delta_ast_RP_lam(end))/delta_ast_RP_lam(end)*100,(delta_theta_SM_lam(end) - delta_theta_RP_lam(end))/delta_theta_RP_lam(end)*100,(H_SM_lam(end) - H_RP_lam(end))/H_RP_lam(end)*100,(Re_x_SM(end) - Re_x_RP(end))/Re_x_RP(end)*100,(Re_x_SM(end) - Re_x_RP(end))/Re_x_RP(end)*100,(Re_theta_SM_lam(end) - Re_theta_RP_lam(end))/Re_theta_RP_lam(end)*100,(Cf_SM_lam(end) - Cf_RP_lam(end))/Cf_RP_lam(end)*100,(Cd_SM_lam(end) - Cd_RP_lam(end))/Cd_RP_lam(end)*100,(r_SM_lam(end) - r_RP_lam(end))/r_RP_lam(end)*100,(cw_SM_lam(end) - cw_RP_lam(end))/cw_RP_lam(end)*100,(Cw_SM_lam(end) - Cw_RP_lam(end))/Cw_RP_lam(end)*100);

fprintf(fid, '    \\bottomrule\n\n');
fprintf(fid, '    \\end{tabular}\n');
% fprintf(fid, '    \\begin{tablenotes}\n');
% fprintf(fid, '        \\item[*] The solution does not converge for $m_2 = \\infty$, and therefore $m_2 = 1000$ is taken as approach.\n');
% fprintf(fid, '    \\end{tablenotes}\n');
fprintf(fid, '\\end{threeparttable}\n');
% fprintf(fid, '%%\\end{sidewaystable}\n');
% fprintf(fid, '%%\\end{table\n');
% fprintf(fid, '\\end{landscape}\n');
fclose(fid);

%% Boundary Layer Characteristics Turbulent Flow
% keyboard

fid = fopen('Comparison_ORCHID_Sim_Results_turb.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Comparison of simulation results of the ORCHID nozzle with two different thermophysical models: RefProp and StanMix. Where StanMix is considered as a reference.\n');
% fprintf(fid, '\\begin{landscape}\n');
fprintf(fid, '%%\\begin{table}\n');
fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Comparison of ORCHID turbulent simulations with two different thermodynamic models: boundary layer edge properties. Results calculated using StanMix and RefProp. StanMix utilizes the iPRSV, whereas RefProp utilizes unknow industrial standards? StanMix is used here as starting point for convenience, since the fluid property retrieval is much faster, whereas RefProp is thought to be more accurate. Note also that edge properties are the same for laminar or turbulent flow. Static pressure obtained from SU2 CFD simulation second-order convergence is used as input.}\n');
fprintf(fid, '    \\label{tab:costurb}\n'); % Comparison of ORCHID Simulations BL Results
fprintf(fid, '    \\begin{tabular}{l S[table-format=1.3] S[table-format=1.3] S[table-format=1.3] S[table-format=1.2] S[table-format=1.1] S[table-format=1.1] S[table-format=1.2] S[table-format=1.2] S[table-format=1.2] S[table-format=1.2] S[table-format=1.2] S[table-format=1.2]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{Turbulent flow}                  &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\delta$}                  &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\delta^*$}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\theta$}                    &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$H$}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{Re}_{\\mathrm{e}}$}                  &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{Re}_x$}                  &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{Re}_{\\theta}$}                &\n');        
fprintf(fid, '        \\multicolumn{1}{c}{$C_{\\mathrm{f}}$}                    &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$C_{\\mathrm{d}}$}                    &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$r$}      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$c_{\\mathrm{w}}$}                    &\n');        
fprintf(fid, '        \\multicolumn{1}{c}{$C_{\\mathrm{w}}$}                    \\\\\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                  &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\,$ [m]}                     &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\,$ [m]}                                &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\,$ [m]}                    &\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{7} \\,$ [1/m]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{7} \\,$ [-]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{3} \\,$ [-]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\,$ [-]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\times \\, 10^{-3} \\,$ [-]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}                      \\\\\n');
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
% S[table-format=1.3] S[table-format=1.3] S[table-format=1.3]
% S[table-format=1.2] S[table-format=1.1] S[table-format=1.2]
% S[table-format=1.2] S[table-format=1.2] S[table-format=1.2]
% S[table-format=1.2] delta_RP_lam  delta_ast_RP_lam  delta_theta_RP_lam  H_RP_lam  Re_x_RP  Re_theta_RP_lam  Cf_RP_lam  r_RP_lam  cw_RP_lam  Cw_RP_lam
fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{RefProp}}            &    %1.3f       &  %1.3f       &   %1.3f       & %1.2f       &   %2.0f    &   %1.1f          &  %1.2f   & %1.3f     & %1.3f     &  %1.3f   &  %1.3f &  %1.3f  \\\\\n',delta_RP_turb(end)*10^3,delta_ast_RP_turb(end)*10^3,delta_theta_RP_turb(end)*10^3,H_RP_turb(end),Re_x_RP(end)/X(end)*10^-7,Re_x_RP(end)*10^-7,Re_theta_RP_turb(end)*10^-3,Cf_RP_turb(end)*10^3,Cd_RP_turb(end)*10^3,r_RP_turb(end),cw_RP_turb(end),Cw_RP_turb(end));
fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{StanMix}}            &    %1.3f       &  %1.3f       &   %1.3f       & %1.2f       &   %2.0f    &   %1.1f          &  %1.2f   & %1.3f     & %1.3f     &  %1.3f   &  %1.3f &  %1.3f  \\\\\n',delta_SM_turb(end)*10^3,delta_ast_SM_turb(end)*10^3,delta_theta_SM_turb(end)*10^3,H_SM_turb(end),Re_x_SM(end)/X(end)*10^-7,Re_x_SM(end)*10^-7,Re_theta_SM_turb(end)*10^-3,Cf_SM_turb(end)*10^3,Cd_SM_turb(end)*10^3,r_SM_turb(end),cw_SM_turb(end),Cw_SM_turb(end));
fprintf(fid, '\\multicolumn{1}{l}{\\textnormal{Deviation [\\%%]}}   &    %2.2f       &  %2.2f       &   %2.2f       & %2.2f       &   %2.1f    &   %2.1f          &  %2.1f   & %2.2f     & %2.2f     &  %2.2f   &  %2.2f &  %2.2f  \\\\\n',(delta_SM_turb(end) - delta_RP_turb(end))/delta_RP_turb(end)*100,(delta_ast_SM_turb(end) - delta_ast_RP_turb(end))/delta_ast_RP_turb(end)*100,(delta_theta_SM_turb(end) - delta_theta_RP_turb(end))/delta_theta_RP_turb(end)*100,(H_SM_turb(end) - H_RP_turb(end))/H_RP_turb(end)*100,(Re_x_SM(end) - Re_x_RP(end))/Re_x_RP(end)*100,(Re_x_SM(end) - Re_x_RP(end))/Re_x_RP(end)*100,(Re_theta_SM_turb(end) - Re_theta_RP_turb(end))/Re_theta_RP_turb(end)*100,(Cf_SM_turb(end) - Cf_RP_turb(end))/Cf_RP_turb(end)*100,(Cd_SM_turb(end) - Cd_RP_turb(end))/Cd_RP_turb(end)*100,(r_SM_turb(end) - r_RP_turb(end))/r_RP_turb(end)*100,(cw_SM_turb(end) - cw_RP_turb(end))/cw_RP_turb(end)*100,(Cw_SM_turb(end) - Cw_RP_turb(end))/Cw_RP_turb(end)*100);

fprintf(fid, '    \\bottomrule\n\n');
fprintf(fid, '    \\end{tabular}\n');
% fprintf(fid, '    \\begin{tablenotes}\n');
% fprintf(fid, '        \\item[*] The solution does not converge for $m_2 = \\infty$, and therefore $m_2 = 1000$ is taken as approach.\n');
% fprintf(fid, '    \\end{tablenotes}\n');
fprintf(fid, '\\end{threeparttable}\n');
% fprintf(fid, '%%\\end{sidewaystable}\n');
% fprintf(fid, '%%\\end{table\n');
% fprintf(fid, '\\end{landscape}\n');
fclose(fid);
OldFolder = cd(NewFolder);

%%
fprintf('Table(s) generated successfully\n')
