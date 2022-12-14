%% Generate graphs for chapter ORCHID Nozzle simulations
% Created (18-06-2020) (2nd time completely revisited)
% Last run ()

% Description: load simulation data and plot to compare Ce, m1, m2, m3, cw,
% Cw; to characterize MM BL in current (inviscid) ORCHID Nozzle design.

close all
clear all
clc

cd ..

%% ORCHID Nozzle - Load and Store data

% ORCHID Nozzle - Fully Laminar - Load and Store data
load('./Stored_Sim_Data/ORCHID_Nozzle_RefProp_PsE_Interpolated_lam')

for j = 1:length(SOL)
    cw_ORCHID_Lam(j) = SOL{j}.c(1);     % cw
    Cw_ORCHID_Lam(j) = FLP{j}.C(1);     % Cw
    Tw_ORCHID_Lam(j) = FLP{j}.T(1);     % Tw
end

% MaE, Ce, m1, m2, m3 are exactly the same for both laminar and turbulent flow: 
X_ORCHID = X;                            % X [m]
MaE_ORCHID = EDG.MaE;                    % Ma [-]
aE_ORCHID = EDG.aE;                      % a [m/s]
Ce_ORCHID = EDG.rhoE.*EDG.muE./FRS.rhotI./FRS.mutI; % Ce [-]
m1_ORCHID = HVR.P1;                      % m1 [-]
m2_ORCHID = HVR.P2;                      % m2 [-]
m3_ORCHID = 2*HVR.P1 - 1 - HVR.P2;       % m3 [-]
delta_ast_Lam = BLC.delta_ast;           % delta_ast [m]
delta_Lam = BLC.delta;                   % delta [m]
H_Lam = BLC.H;                           % H [-]

% ORCHID Nozzle - Fully Turbulent - Load and Store data
load('./Stored_Sim_Data/ORCHID_Nozzle_RefProp_PsE_Interpolated_turb')

for j = 1:length(SOL)
    cw_ORCHID_Turb(j) = SOL{j}.c(1);     % cw
    Cw_ORCHID_Turb(j) = FLP{j}.C(1);     % Cw
    Tw_ORCHID_Turb(j) = FLP{j}.T(1);     % Tw
end

delta_ast_Turb = BLC.delta_ast;          % delta_ast [m]
delta_Turb = BLC.delta;                  % delta [m]
H_Turb = BLC.H;                          % H [-]

% ORCHID Nozzle - Fully Laminar Calorically Perfect Ideal Gas (CPIG) - Load and Store data
load('./Stored_Sim_Data/ORCHID_Nozzle_PsE_Interpolated_lam_IG_CalPerf')

for j = 1:length(SOL)
    cw_ORCHID_Lam_CPIG(j) = SOL{j}.c(1); % cw
    Cw_ORCHID_Lam_CPIG(j) = FLP{j}.C(1); % Cw
    Tw_ORCHID_Lam_CPIG(j) = FLP{j}.T(1); % Tw
end

% Ce, m1, m2, m3 are exactly the same for laminar and turbulent flow (not influenced by turbulent properties)
MaE_ORCHID_CPIG = EDG.MaE;               % Ma
aE_ORCHID_CPIG = EDG.aE;                 % a
Ce_ORCHID_CPIG = EDG.rhoE.*EDG.muE./FRS.rhotI./FRS.mutI; % Ce [-]
delta_ast_Lam_CPIG = BLC.delta_ast;      % delta_ast [m]
delta_Lam_CPIG = BLC.delta;              % delta [m]
H_Lam_CPIG = BLC.H;                      % H [-]

% ORCHID Nozzle - Fully Laminar Thermally Perfect Ideal Gas (TPIG) - Load and Store data
load('./Stored_Sim_Data/ORCHID_Nozzle_PsE_Interpolated_lam_IG_Thermally_Perf')

for j = 1:length(SOL)
    cw_ORCHID_Lam_TPIG(j) = SOL{j}.c(1); % cw
    Cw_ORCHID_Lam_TPIG(j) = FLP{j}.C(1); % Cw
    Tw_ORCHID_Lam_TPIG(j) = FLP{j}.T(1); % Tw
end

% Ce, m1, m2, m3 are exactly the same for laminar and turbulent flow (not influenced by turbulent properties)
MaE_ORCHID_TPIG = EDG.MaE;               % Ma
aE_ORCHID_TPIG = EDG.aE;                 % a
Ce_ORCHID_TPIG = EDG.rhoE.*EDG.muE./FRS.rhotI./FRS.mutI; % Ce [-]
delta_ast_Lam_TPIG = BLC.delta_ast;      % delta_ast [m]
delta_Lam_TPIG = BLC.delta;              % delta [m]
H_Lam_TPIG = BLC.H;                      % H [-]

% ORCHID Nozzle - Fully Laminar Thermally Perfect Ideal Gas (TPIG) - Load and Store data
% NB Cp polynomial fitted to TsE (static temperature at BL edge along ORCHID nozzle expansion)
load('./Stored_Sim_Data/ORCHID_Nozzle_PsE_Interpolated_lam_IG_Thermally_Perf_Cp_poly_edge')

for j = 1:length(SOL)
    cw_ORCHID_Lam_TPIG_CpE_poly(j) = SOL{j}.c(1); % cw
    Cw_ORCHID_Lam_TPIG_CpE_poly(j) = FLP{j}.C(1); % Cw
    Tw_ORCHID_Lam_TPIG_CpE_poly(j) = FLP{j}.T(1); % Tw
end

MaE_ORCHID_TPIG_CpE_poly = EDG.MaE;      % Ma
aE_ORCHID_TPIG_CpE_poly = EDG.aE;        % a

% NB 2 This case does not add higher accurcy to cw or Cw and therefore it
% is not considered here, only in bonus plot

% ORCHID Nozzle StanMix Simulation - Fully Laminar - Load and Store data
load('Stored_Sim_Data/ORCHID_Nozzle_StanMix_PsE_Interpolated_lam')

for j = 1:length(SOL)
    cw_ORCHID_Lam_SM(j) = SOL{j}.c(1);   % cw
    Cw_ORCHID_Lam_SM(j) = FLP{j}.C(1);   % Cw
    Tw_ORCHID_Lam_SM(j) = FLP{j}.T(1);   % Tw
end

% MaE, Ce, m1, m2, m3 are exactly the same for both laminar and turbulent flow: 
% X_ORCHID_SM = X;                        % X [m] %%% equal
MaE_ORCHID_SM = EDG.MaE;                 % Ma
aE_ORCHID_SM = EDG.aE;                   % a
Ce_ORCHID_SM = EDG.rhoE.*EDG.muE./FRS.rhotI./FRS.mutI; % Ce [-]
m1_ORCHID_SM = HVR.P1;                   % m1
m2_ORCHID_SM = HVR.P2;                   % m2
m3_ORCHID_SM = 2*HVR.P1 - 1 - HVR.P2;    % m3
delta_ast_Lam_SM = BLC.delta_ast;        % delta_ast [m]
delta_Lam_SM = BLC.delta;                % delta [m]
H_Lam_SM = BLC.H;                        % H [-]

%% Winter profile 19 - Load and Store data
load('Stored_Sim_Data/ADA045367_7302_Winter1973_P19')

for j = 1:length(SOL)
    cw_CPIG_AIR_Winter_P19(j) = SOL{j}.c(1); % cw
    Cw_CPIG_AIR_Winter_P19(j) = FLP{j}.C(1); % Cw
    Tw_CPIG_AIR_Winter_P19(j) = FLP{j}.T(1); % Tw
end

% MaE, Ce, m1, m2, m3 are exactly the same for both laminar and turbulent flow: 
X_Winter = X;                       % X [m]
MaE_Winter = EDG.MaE;               % Ma
Ce_Winter = EDG.rhoE.*EDG.muE./FRS.rhotI./FRS.mutI; % Ce [-]
m1_Winter = HVR.P1;                 % m1
m2_Winter = HVR.P2;                 % m2
m3_Winter = 2*HVR.P1 - 1 - HVR.P2;  % m3

% Winter profile 19 - Fully Laminar - Load and Store data
load('Stored_Sim_Data/ADA045367_7302_Winter1973_P19_Lam')

for j = 1:length(SOL)
    cw_CPIG_AIR_Winter_P19_Lam(j) = SOL{j}.c(1); % cw
    Cw_CPIG_AIR_Winter_P19_Lam(j) = FLP{j}.C(1); % Cw
    Tw_CPIG_AIR_Winter_P19_Lam(j) = FLP{j}.T(1); % Tw
end

% Winter profile 19 - Fully Turbulent - Load and Store data
load('Stored_Sim_Data/ADA045367_7302_Winter1973_P19_Turb')

for j = 1:length(SOL)
    cw_CPIG_AIR_Winter_P19_Turb(j) = SOL{j}.c(1); % cw
    Cw_CPIG_AIR_Winter_P19_Turb(j) = FLP{j}.C(1); % Cw
    Tw_CPIG_AIR_Winter_P19_Turb(j) = FLP{j}.T(1); % Tw
end

%% Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
Color = [0 0 0]; % black
Color1 = [0.94 0.94 0.94]; % ligth grey (too light for thin lines, used for surfaces)
Color2 = [112 128 144]/255; % slate (darker) grey
Color3 = [192 192 192]/255; % silver
Color4 = 'r';
LineStyle = {'-','--','-.',':','-',':'};
Marker = {'none','none','none','none','none','none'}; % Marker = {'+','x','s','d','o'};
LineWidth = 0.5;
LineWidth1 = 1;
LineWidth2 = 1.2;1.0;
LineWidth3 = 0.8; % Linewidth corresponding to Color3
FontSize = 11;
FontSize1 = 9;
FontSize2 = 12;
FontSize3 = 14;
LegendPosition4b = [0.421190476000666,0.493969006360160,0.390261695365674,0.180714278811500];
% https://www.rapidtables.com/web/color/gray-color.html

%% Format Figures

% Figure 1 - Ce
h1 = figure;
ax1 = gca;

ax1.XLabel.Interpreter = 'latex';
ax1.XLabel.Color = Color;
ax1.XLabel.String = '$\mathrm{Ma} \, [-]$';
ax1.XLim = [-0.03 2.03];
ax1.XTick = [0 0.5 1.0 1.5 2.0];
ax1.XTickLabel = {'0'; '0.5'; '1.0'; '1.5'; '2.0'};

ax1.YLabel.Interpreter = 'latex';
ax1.YLabel.Color = Color;
ax1.YLabel.String = '$C_{\mathrm{e}} \, [-]$';
ax1.YLim = [-0.01 1.01];
ax1.YTick = [0 0.2 0.4 0.6 0.8 1.0];
ax1.YTickLabel = {'0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'};

ax1.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
ax1.TickLabelInterpreter = 'latex';
set(gca,'FontSize',FontSize)
box on
hold on

% Figure 2 - m1, m2, m3
h2 = figure;
ax2 = gca;

ax2.XLabel.Interpreter = 'latex';
ax2.XLabel.Color = Color;
ax2.XLabel.String = '$X \, [\mathrm{m}]$';
ax2.XLim = [0 0.1];
ax2.XTick = [0 0.02 0.04 0.06 0.08 0.10];
ax2.XTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};

ax2.YLabel.Interpreter = 'latex';
ax2.YLabel.Color = Color;
ax2.YLabel.String = '${\mathrm{PG}} \, [-]$';
ax2.YLim = [-7 4];
ax2.YTick = [-7 -6 -4 -2 0 2 4];
ax2.YTickLabel = {'-7'; '-6'; '-4'; '-2'; '0'; '2'; '4'};

ax2.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
ax2.TickLabelInterpreter = 'latex';
set(gca,'FontSize',FontSize)
box on
hold on

% Figure 3 - cw
h3 = figure;
ax3 = gca;

ax3.XLabel.Interpreter = 'latex';
ax3.XLabel.Color = Color;
ax3.XLabel.String = '$\mathrm{Ma} \, [-]$';
ax3.XLim = [-0.03 2.03];
ax3.XTick = [0 0.5 1.0 1.5 2.0];
ax3.XTickLabel = {'0'; '0.5'; '1.0'; '1.5'; '2.0'};

ax3.YLabel.Interpreter = 'latex';
ax3.YLabel.Color = Color;
ax3.YLabel.String = '$c_{\mathrm{w}} \, [-]$';
ax3.YLim = [0.99 1.81];
ax3.YTick = [1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8];
ax3.YTickLabel = {'1.0'; '1.1'; '1.2'; '1.3'; '1.4'; '1.5'; '1.6'; '1.7'; '1.8'};

ax3.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
ax3.TickLabelInterpreter = 'latex';
set(gca,'FontSize',FontSize)
box on
hold on

% Figure 3 - cw Zoom - ORCHID Nozzle only
h3b = figure;
ax3b = gca;

ax3b.XLabel.Interpreter = 'latex';
ax3b.XLabel.Color = Color;
ax3b.XLabel.String = '$\mathrm{Ma} \, [-]$';
% ax3b.XLim = [-0.03 2.03];
% ax3b.XTick = [0 0.5 1.0 1.5 2.0];
% ax3b.XTickLabel = {'0'; '0.5'; '1.0'; '1.5'; '2.0'};
ax3b.XLim = [0.11 2.03];
ax3b.XTick = [0.15 0.3 0.5 1.0 1.5 2.0];
ax3b.XTickLabel = {'0.15'; '0.3'; '0.5'; '1.0'; '1.5'; '2.0'};

ax3b.YLabel.Interpreter = 'latex';
ax3b.YLabel.Color = Color;
ax3b.YLabel.String = '$c_{\mathrm{w}} \, [-]$';
ax3b.YLim = [0.998 1.0602];
ax3b.YTick = [1.0 1.01 1.02 1.03 1.04 1.05 1.06];
ax3b.YTickLabel = {'1.00'; '1.01'; '1.02'; '1.03'; '1.04'; '1.05'; '1.06'};

ax3b.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
ax3b.TickLabelInterpreter = 'latex';
set(gca,'FontSize',FontSize)
box on
hold on

% Figure 4 - Cw
h4 = figure;
ax4 = gca;

ax4.XLabel.Interpreter = 'latex';
ax4.XLabel.Color = Color;
ax4.XLabel.String = '$\mathrm{Ma} \, [-]$';
ax4.XLim = [0 2.03];
ax4.XTick = [0 0.5 1.0 1.5 2.0];
ax4.XTickLabel = {'0'; '0.5'; '1.0'; '1.5'; '2.0'};

ax4.YLabel.Interpreter = 'latex';
ax4.YLabel.Color = Color;
ax4.YLabel.String = '$C_{\mathrm{w}} \, [-]$';
ax4.YLim = [0.915 1.005];
ax4.YTick = [0.86 0.88 0.90 0.92 0.94 0.96 0.98 1.00 1.02];
ax4.YTickLabel = {'0.86'; '0.88'; '0.90'; '0.92'; '0.94'; '0.96'; '0.98'; '1.00'; '1.02'};

ax4.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
ax4.TickLabelInterpreter = 'latex';
set(gca,'FontSize',FontSize)
box on
hold on

% Figure 4 - Cw
h4b = figure;
ax4b = gca;

ax4b.XLabel.Interpreter = 'latex';
ax4b.XLabel.Color = Color;
ax4b.XLabel.String = '$\mathrm{Ma} \, [-]$';
ax4b.XLim = [0.11 2.03];
ax4b.XTick = [0.15 0.3 0.5 1.0 1.5 2.0];
ax4b.XTickLabel = {'0.15'; '0.3'; '0.5'; '1.0'; '1.5'; '2.0'};

ax4b.YLabel.Interpreter = 'latex';
ax4b.YLabel.Color = Color;
ax4b.YLabel.String = '$C_{\mathrm{w}} \, [-]$';
ax4b.YLim = [0.985 1.0005];
ax4b.YTick = [0.986 0.988 0.990 0.992 0.994 0.996 0.998 1.00];
ax4b.YTickLabel = {'0.986'; '0.988'; '0.990'; '0.992'; '0.994'; '0.996'; '0.998'; '1.00'};

ax4b.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
ax4b.TickLabelInterpreter = 'latex';
set(gca,'FontSize',FontSize)
box on
hold on

% figure 5: Form factor
h5 = figure;
ax5 = gca;

ax5.XLabel.Interpreter = 'latex';
ax5.XLabel.Color = Color;
ax5.XLabel.String = '$X \, [\mathrm{m}]$';
ax5.XLim = [0 0.1];
ax5.XTick = [0 0.02 0.04 0.06 0.08 0.10];
ax5.XTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};

ax5.YLabel.Interpreter = 'latex';
ax5.YLabel.Color = Color;
ax5.YLabel.String = '$H \, [-]$';
% ax5.YLim = [2 2.8];
% ax5.YTick = [2.0 2.2 2.4 2.6 2.8];
% ax5.YTickLabel = {'2.0'; '2.2'; '2.4'; '2.6'; '2.8'};
ax5.YLim = [2 2.7];
ax5.YTick = [2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7];
ax5.YTickLabel = {'2.0'; '2.1'; '2.2'; '2.3'; '2.4'; '2.5'; '2.6'; '2.7'};

ax5.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
ax5.TickLabelInterpreter = 'latex';
set(gca,'FontSize',FontSize)
box on
hold on

% figure 6: Displacement thickness
h6 = figure;
ax6 = gca;

ax6.XLabel.Interpreter = 'latex';
ax6.XLabel.Color = Color;
ax6.XLabel.String = '$X \, [\mathrm{m}]$';
ax6.XLim = [0 0.1];
ax6.XTick = [0 0.02 0.04 0.06 0.08 0.10];
ax6.XTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};

ax6.YLabel.Interpreter = 'latex';
ax6.YLabel.Color = Color;
ax6.YLabel.String = '$\delta^* \, \times 10^{-3} \, [\mathrm{m}]$';
ax6.YLim = [0 2.65e-5];
ax6.YTick = [0 0.5e-5 1.0e-5 1.5e-5 2.0e-5 2.5e-5];
ax6.YTickLabel = {'0'; '0.05'; '0.10'; '0.15'; '0.20'; '0.25'};

ax6.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
ax6.TickLabelInterpreter = 'latex';
set(gca,'FontSize',FontSize)
box on
hold on

% figure 7: Velocity thickness
h7 = figure;
ax7 = gca;

ax7.XLabel.Interpreter = 'latex';
ax7.XLabel.Color = Color;
ax7.XLabel.String = '$X \, [\mathrm{m}]$';
ax7.XLim = [0 0.1];
ax7.XTick = [0 0.02 0.04 0.06 0.08 0.10];
ax7.XTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};

ax7.YLabel.Interpreter = 'latex';
ax7.YLabel.Color = Color;
ax7.YLabel.String = '$\delta \, \times 10^{-3} \, [\mathrm{m}]$';
ax7.YLim = [0 8e-5];
ax7.YTick = [0 2.0e-5 4.0e-5 6.0e-5 8.0e-5];
ax7.YTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'};

ax7.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';
ax7.TickLabelInterpreter = 'latex';
set(gca,'FontSize',FontSize)
box on
hold on

%% Plot data
% in right order for Legends
% Rules/guidelines regarding Colours:
% dotted line for Calorically Perfect Ideal Gas (CPIG)
% dashed line for Thermally Perfect Ideal Gas (TPIG)
% Solid line for Nonideal Gas (NIG)
% Rules/guidelines regarding Colours:
% ORCHID: black
% Air: grey

% Ce
figure(h1)
% NB CPIG and TPIG practically overlap/superpose, leave one out
% plot(MaE_ORCHID_CPIG,Ce_ORCHID_CPIG,'Color',Color,'LineStyle',LineStyle{4},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM CPIG')
plot(MaE_ORCHID_TPIG,Ce_ORCHID_TPIG,'Color',Color,'LineStyle',LineStyle{2},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM TPIG')
plot(MaE_Winter,Ce_Winter,'Color','r','LineStyle',LineStyle{3},'LineWidth',LineWidth1,'DisplayName','AIR CPIG') % Winter p19
plot(MaE_ORCHID,Ce_ORCHID,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM RefProp')

% m1, m2, m3
figure(h2)

plot(X_ORCHID,m2_ORCHID,'Color',Color,'LineStyle',LineStyle{2},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','$m_2$')
plot(X_ORCHID,m1_ORCHID,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','$m_1$')
plot(X_ORCHID,m3_ORCHID,'Color',Color,'LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','$m_3$')

% cw
figure(h3)

plot(MaE_Winter,cw_CPIG_AIR_Winter_P19_Turb,'Color',Color2,'LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','AIR CPIG Turb')
plot(MaE_Winter,cw_CPIG_AIR_Winter_P19_Lam,'Color',Color2,'LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','AIR CPIG Winter Lam')
plot(MaE_Winter,cw_CPIG_AIR_Winter_P19,'Color','r','LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','AIR CPIG tr')
% plot(MaE_ORCHID,cw_ORCHID_MM,'Color',Color2,'LineStyle',LineStyle{1},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','ORCHID MM Lam')
plot(MaE_ORCHID,cw_ORCHID_Turb,'Color',Color2,'LineStyle',LineStyle{1},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM RefProp Turb')
plot(MaE_ORCHID,cw_ORCHID_Lam,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM RefProp Lam')
plot(MaE_ORCHID_TPIG,cw_ORCHID_Lam_TPIG,'Color',Color,'LineStyle',LineStyle{2},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM TPIG Lam')
plot(MaE_ORCHID_CPIG,cw_ORCHID_Lam_CPIG,'Color',Color,'LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM CPIG Lam')

% cw - Zoom - ORCHID Nozzle only
figure(h3b)

% plot(MaE_Winter,cw_CPIG_AIR_Winter_P19,'Color',Color2,'LineStyle',LineStyle{4},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','CPIG AIR Winter')
plot(MaE_ORCHID,cw_ORCHID_Turb,'Color',Color2,'LineStyle',LineStyle{1},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM RefProp Turb')
plot(MaE_ORCHID,cw_ORCHID_Lam,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM RefProp Lam')
plot(MaE_ORCHID_TPIG,cw_ORCHID_Lam_TPIG,'Color',Color,'LineStyle',LineStyle{2},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM TPIG Lam')
plot(MaE_ORCHID_CPIG,cw_ORCHID_Lam_CPIG,'Color',Color,'LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM CPIG Lam')

% Cw
figure(h4)

plot(MaE_ORCHID_CPIG,Cw_ORCHID_Lam_CPIG,'Color',Color,'LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM CPIG Lam')
plot(MaE_ORCHID_TPIG,Cw_ORCHID_Lam_TPIG,'Color',Color,'LineStyle',LineStyle{2},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM TPIG Lam')
plot(MaE_ORCHID,Cw_ORCHID_Lam,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM RefProp Lam')
plot(MaE_ORCHID,Cw_ORCHID_Turb,'Color',Color2,'LineStyle',LineStyle{1},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM RefProp Turb')
plot(MaE_Winter,Cw_CPIG_AIR_Winter_P19_Lam,'Color',Color2,'LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','AIR CPIG Lam')
plot(MaE_Winter,Cw_CPIG_AIR_Winter_P19_Turb,'Color',Color2,'LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','AIR CPIG Turb')
plot(MaE_Winter,Cw_CPIG_AIR_Winter_P19,'Color','r','LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','AIR CPIG tr')

% Cw - Zoom - ORCHID Nozzle only
figure(h4b)

plot(MaE_ORCHID_CPIG,Cw_ORCHID_Lam_CPIG,'Color',Color,'LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM CPIG Lam')
plot(MaE_ORCHID_TPIG,Cw_ORCHID_Lam_TPIG,'Color',Color,'LineStyle',LineStyle{2},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM TPIG Lam')
plot(MaE_ORCHID,Cw_ORCHID_Lam,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM RefProp Lam')
plot(MaE_ORCHID,Cw_ORCHID_Turb,'Color',Color2,'LineStyle',LineStyle{1},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','MM RefProp Turb')
% plot(MaE_Winter,Cw_CPIG_AIR_Winter_P19_Lam,'Color',Color2,'LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','CPIG AIR Winter Lam')
% plot(MaE_Winter,Cw_CPIG_AIR_Winter_P19,'Color','r','LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','CPIG AIR Winter')
% plot(MaE_Winter,Cw_CPIG_AIR_Winter_P19_Turb,'Color',Color2,'LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','CPIG AIR Winter Turb')

% H - Form Factor ORCHID Nozzle Laminar
figure(h5)

plot(X_ORCHID,H_Lam_TPIG,'Color',Color,'LineStyle',LineStyle{2},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','TPIG')
plot(X_ORCHID,H_Lam,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','RefProp')
plot(X_ORCHID,H_Lam_SM,'Color',Color,'LineStyle',LineStyle{4},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','StanMix')
plot(X_ORCHID,H_Lam_CPIG,'Color',Color,'LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','CPIG')
% plot(X_ORCHID,H_Turb)

% delta_ast - Displacement Thickness ORCHID Nozzle Laminar
figure(h6)

% plot(X_ORCHID,delta_ast_Turb)
plot(X_ORCHID,delta_ast_Lam_SM,'Color',Color,'LineStyle',LineStyle{4},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','StanMix')
plot(X_ORCHID,delta_ast_Lam_TPIG,'Color',Color,'LineStyle',LineStyle{2},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','TPIG')
plot(X_ORCHID,delta_ast_Lam_CPIG,'Color',Color,'LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','CPIG')
plot(X_ORCHID,delta_ast_Lam,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','RefProp')

% delta - Velocity Thickness ORCHID Nozzle Laminar
figure(h7)

% plot(X_ORCHID,delta_Turb)
plot(X_ORCHID,delta_Lam_SM,'Color',Color,'LineStyle',LineStyle{4},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','StanMix')
plot(X_ORCHID,delta_Lam_TPIG,'Color',Color,'LineStyle',LineStyle{2},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','TPIG')
plot(X_ORCHID,delta_Lam_CPIG,'Color',Color,'LineStyle',LineStyle{3},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','CPIG')
plot(X_ORCHID,delta_Lam,'Color',Color,'LineStyle',LineStyle{1},'LineWidth',LineWidth1,'Marker',Marker{1},'MarkerFaceColor','none','HandleVisibility','on','DisplayName','RefProp')

%% Add Legends

% Ce
figure(h1)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northeast');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

% m1, m2, m3
figure(h2)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northeast');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

% cw
figure(h3)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

% cw - Zoom
figure(h3b)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

% Cw
figure(h4)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'southwest');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

% Cw - Zoom
figure(h4b)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Position', LegendPosition4b);%'Location', 'southwest');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

% H
figure(h5)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'north');%'Location', 'southwest');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

% delta_ast
figure(h6)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest');%'Location', 'southwest');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

% delta
figure(h7)
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest');%'Location', 'southwest');%'Location', 'best'); %%%  'Position', LgndPos2); %
legend('show')

%% Bonus Plots

% Mach-number - ORCHID IG vs NIG
hx = figure;
hold on
plot(X_ORCHID,MaE_ORCHID)
plot(X_ORCHID,MaE_ORCHID_CPIG)
plot(X_ORCHID,MaE_ORCHID_TPIG)
plot(X_ORCHID,MaE_ORCHID_TPIG_CpE_poly)
plot(X_ORCHID,MaE_ORCHID_SM)
xlabel('X [m]')
ylabel('Ma [-]')
legend('RefProp','CPIG','TPIG','TPIG Cp(TsE)','StanMix', 'Location', 'northwest');

% Speed of Sound - ORCHID IG vs NIG
hy = figure;
hold on
plot(X_ORCHID,aE_ORCHID)
plot(X_ORCHID,aE_ORCHID_CPIG)
plot(X_ORCHID,aE_ORCHID_TPIG)
plot(X_ORCHID,aE_ORCHID_TPIG_CpE_poly)
plot(X_ORCHID,aE_ORCHID_SM)
xlabel('X [m]')
ylabel('SoS [m/s]')
legend('RefProp','CPIG','TPIG','TPIG Cp(TsE)','StanMix', 'Location', 'southeast');

%% Save (Thesis Report) Figures
% keyboard

saveas(h1,'./Thesis_Report_Figures/ORCHID_Ce','epsc')
saveas(h2,'./Thesis_Report_Figures/ORCHID_PGmx','epsc')
saveas(h3,'./Thesis_Report_Figures/ORCHID_c_w','epsc')
saveas(h3b,'./Thesis_Report_Figures/ORCHID_c_w_Zoom','epsc')
saveas(h4,'./Thesis_Report_Figures/ORCHID_Cw','epsc')
saveas(h4b,'./Thesis_Report_Figures/ORCHID_Cw_Zoom','epsc')
saveas(h5,'./Thesis_Report_Figures/ORCHID_H_Lam','epsc')
saveas(h6,'./Thesis_Report_Figures/ORCHID_delta_ast_Lam','epsc')
saveas(h7,'./Thesis_Report_Figures/ORCHID_delta_Lam','epsc')

fprintf('Figures created and saved succesfully\n')
%