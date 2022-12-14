%% Validation with NACA0012 McNally 1970

clear all
close all
clc

pause(0.1)

run('INPUT/INPUT_Nozzle_FlatWall_McNally1970')
cd ..

% PRECAL
[X,NST,HVR,EDG,FRS,FLD,TCC,SOL,BLC,MON] = PRECAL(INP,OPT,SET,FLD,TCC);

% Grid generation
[NP,GRD] = GRID(GRD,SET);

NS = 1;

% IVPL
[sol,solprev] = IVPL(NS,GRD,HVR,FLD,SET); %%% solprev = []; (empty, maybe needed for some case? move to PRECAL?)
flp = [];

tic

[BLC,FLP,SOL,TCC,MON] = CSM(X,NP,NS,NST,GRD,HVR,EDG,FRS,FLD,TCC,sol,solprev,SOL,BLC,INP,OPT,SET,MON);

toc

%% Plot Results and generate graphs and tables
fprintf('Simulation finished with success \n')

STORDATA

%% Experimental Data
load('Validation_Study/Data/delta_ast_FlatWalledNozzle_McNally1970')
load('Validation_Study/Data/theta_FlatWalledNozzle_McNally1970')

% %% PLOT
% figure
% hold on
% % plot(X(1:FST),BLC.edge)                     % grid
% % plot(X(1:FST),BLC.delta)                    % velocity
% plot(X(1:length(SOL)),BLC.delta_ast,'-k','DisplayName','\delta^* McNally')        % displacement
% plot(X(1:length(SOL)),BLC.theta,'-k','DisplayName','\theta McNally')          % momentum
% plot(0.3048*delta_ast_FlatWalledNozzle_McNally1970(:,1) + INP.L_FP,0.3048*delta_ast_FlatWalledNozzle_McNally1970(:,2),'ok','DisplayName','\delta^* McNally')        % displacement
% plot(0.3048*theta_FlatWalledNozzle_McNally1970(:,1) + INP.L_FP,0.3048*theta_FlatWalledNozzle_McNally1970(:,2),'dk','DisplayName','\theta McNally')          % momentum

%% PLOT
h1 = figure;
hold on

% plot(X(1:FST),BLC.edge)                     % grid
% plot(X(1:FST),BLC.delta)                    % velocity

% Experimental data Lewis Research Center
plot(0.3048*delta_ast_FlatWalledNozzle_McNally1970(:,1) + INP.L_FP,0.3048*delta_ast_FlatWalledNozzle_McNally1970(:,2),'ok','DisplayName','$\delta^* \, \mathrm{McNally}$') % displacement
plot(0.3048*theta_FlatWalledNozzle_McNally1970(:,1) + INP.L_FP,0.3048*theta_FlatWalledNozzle_McNally1970(:,2),'dk','DisplayName','$\theta \, \mathrm{McNally}$')           % momentum

% CSM
plot(X(1:length(SOL)),BLC.delta_ast,'-k','DisplayName','CSM')        % displacement
plot(X(1:length(SOL)),BLC.theta,'-k','HandleVisibility','off')       % momentum

%% Settings
Color = [0 0 0];
LineStyle = {'-','--',':','-.'};
Marker = {'o','^','s','d','p','v','h','>','+','x','*'}; %%% open symbols (first four (fill), second four (fill), last three)
% Marker = {'+','x','s','d','o'};
% Marker = {'o','^','s','d'};
LineWidth = 0.5;
FontSize = 11;

%%
ax1 = gca;

ax1.FontSizeMode = 'manual';
ax1.FontSize = FontSize; %%
ax1.FontWeight = 'normal';
ax1.TickLabelInterpreter = 'latex';

ax1.XLabel.Interpreter = 'latex';
ax1.XLabel.FontSize = FontSize; %%
ax1.XLabel.FontWeight = 'normal';
ax1.XLabel.Color = Color;
ax1.XLabel.String = '$X \, [\mathrm{m}]$';
ax1.XLim = [0 1.2];
ax1.XLimMode = 'manual';
ax1.XTickLabelMode = 'manual';
ax1.XTickMode = 'manual';
ax1.XTick = [0 0.2 0.4 0.6 0.8 1.0 1.2];
ax1.XTickLabel = {'0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'; '1.2'};

ax1.YLabel.Interpreter = 'latex';
ax1.YLabel.FontSize = FontSize; %%
ax1.YLabel.FontWeight = 'normal';
ax1.YLabel.Color = Color;
ax1.YLabel.String = '$\delta^{\ast}, \, \theta \, \times \, 10^{-3} [\mathrm{m}]$';
ax1.YLim = [0 1.5e-3];
ax1.YLimMode = 'manual';
ax1.YTickLabelMode = 'manual';
ax1.YTickMode = 'manual';
ax1.YTick = [0 0.5e-3 1.0e-3 1.5e-3];
ax1.YTickLabel = {'0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'; '3.0'; '3.5'};
% ax1.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on

set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest'); % 'Location', 'best');
legend('show')

%% Store plots
% keyboard
NewFolder = cd;
% Figure (complete) for Appendix

saveas(h1,'./Thesis_Report_Figures/NozzleFlatWall_McNally1970_Complete','epsc')

% Figure for Chapter 4
ax1.XLim = [0.68 1.06];
ax1.XLimMode = 'manual';
ax1.XTickLabelMode = 'manual';
ax1.XTickMode = 'manual';
ax1.XTick = [0.0 + 0.6851 0.05 + 0.6851 0.10 + 0.6851 0.15 + 0.6851 0.20 + 0.6851 0.25 + 0.6851 0.30 + 0.6851 0.35 + 0.6851];
ax1.XTickLabel = {'0'; '0.05'; '0.10'; '0.15'; '0.20'; '0.25'; '0.30'; '0.35'};
ax1.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'southwest'); % 'Location', 'best');
legend('show')

saveas(h1,'./Thesis_Report_Figures/NozzleFlatWall_McNally1970','epsc')
OldFolder = cd(NewFolder);

fprintf('Figures successfully generated and saved\n')

%%
% keyboard
%%% uncomment if needed
% PLOTFILE

% TABLEGEN
% PLOTCHART % Pv-, and Ts-diagram, calculate critical properties inside
% if OPT.CHRT > 0
%     FPCHARTS
% end

%% Clear variables or objects
% clear NP, sol, solprev % not really needed, overwritten? // and other free variables? in case of multiple calculations
% Clean-up FluidProp (if used)
% if OPT.GASM == 3 % Nonideal Gas
%     Cleanup_FluidProp
% end