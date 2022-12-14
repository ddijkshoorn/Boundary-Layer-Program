%% Validation with NACA0012 McNally 1970

clear all
close all
clc

pause(0.1)

run('INPUT/INPUT_NACA0012_McNally1970')
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
load('Validation_Study/Data/delta_ast_NACA0012_McNally1970')
load('Validation_Study/Data/theta_NACA0012_McNally1970')

%% PLOT
h1 = figure;
hold on
% plot(X(1:FST),BLC.edge)                     % grid
% plot(X(1:FST),BLC.delta)                    % velocity

plot(0.3048*delta_ast_NACA0012_McNally1970(:,1),0.3048*delta_ast_NACA0012_McNally1970(:,2),'ok','DisplayName','$\delta^* \, \mathrm{McNally}$')        % displacement
plot(0.3048*theta_NACA0012_McNally1970(:,1),0.3048*theta_NACA0012_McNally1970(:,2),'dk','DisplayName','$\theta \, \mathrm{McNally}$')         % momentum
plot(X(1:length(SOL)),BLC.delta_ast,'-k','DisplayName','CSM')        % displacement
plot(X(1:length(SOL)),BLC.theta,'-k','HandleVisibility','off')          % momentum

%% Settings
Color = [0 0 0];
Color3 = [112 128 144]/255;
LineStyle = {'-','--',':','-.'};
Marker = {'o','^','s','d','p','v','h','>','+','x','*'}; %%% open symbols (first four (fill), second four (fill), last three)
% Marker = {'+','x','s','d','o'};
% Marker = {'o','^','s','d'};
LineWidth = 0.5;
FontSize = 11;

plot([X(MON.NTR) X(MON.NTR)],[0 1.0e-3],'Color',Color3,'DisplayName','Transition forced')

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
ax1.XLim = [0 1.6];
ax1.XLimMode = 'manual';
ax1.XTickLabelMode = 'manual';
ax1.XTickMode = 'manual';
ax1.XTick = [0 0.2 0.4 0.6 0.8 1.0 1.2 1.4 1.6];
ax1.XTickLabel = {'0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'; '1.2'; '1.4'; '1.6'};

ax1.YLabel.Interpreter = 'latex';
ax1.YLabel.FontSize = FontSize; %%
ax1.YLabel.FontWeight = 'normal';
ax1.YLabel.Color = Color;
ax1.YLabel.String = '$\delta^{\ast}, \, \theta \, \times \, 10^{-3} [\mathrm{m}]$';
ax1.YLim = [0 4e-3];
ax1.YLimMode = 'manual';
ax1.YTickLabelMode = 'manual';
ax1.YTickMode = 'manual';
ax1.YTick = [0 0.5e-3 1.0e-3 1.5e-3 2.0e-3 2.5e-3 3.0e-3 3.5e-3 4.0e-3];
ax1.YTickLabel = {'0'; '0.5'; '1.0'; '1.5'; '2.0'; '2.5'; '3.0'; '3.5'; '4.0'};
ax1.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on

set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest'); % 'Location', 'best');
legend('show')

%% Store plots
% keyboard
NewFolder = cd;
saveas(h1,'./Thesis_Report_Figures/NACA0012_McNally1970','epsc')
OldFolder = cd(NewFolder);

fprintf('Figures successfully generated and saved\n')

%% 
% PLOTFILE % uncomment this one to see all characteristics

%%
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