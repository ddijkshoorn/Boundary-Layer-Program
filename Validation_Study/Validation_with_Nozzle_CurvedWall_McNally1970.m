%% Validation with NACA0012 McNally 1970

clear all
close all
clc

pause(0.1)

run('INPUT/INPUT_Nozzle_CurvedWall_McNally1970')
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
load('Validation_Study/Data/delta_ast_CurvedWalledNozzle_McNally1970')
load('Validation_Study/Data/theta_CurvedWalledNozzle_McNally1970')

% Predicted by FORTRAN code McNally
load('Validation_Study/Data/delta_ast_CurvedNozzle_FORTRAN_CODE_McNally1970')
load('Validation_Study/Data/theta_CurvedNozzle_FORTRAN_CODE_McNally1970')

%% PLOT
h1 = figure;
hold on

% plot(X(1:FST),BLC.edge)                     % grid
% plot(X(1:FST),BLC.delta)                    % velocity

% Experimental data Lewis Research Center
plot(0.3048*delta_ast_CurvedWalledNozzle_McNally1970(:,1) + INP.L_FP,0.3048*delta_ast_CurvedWalledNozzle_McNally1970(:,2),'ok','DisplayName','$\delta^* \, \mathrm{McNally}$') % displacement
plot(0.3048*theta_CurvedWalledNozzle_McNally1970(:,1) + INP.L_FP,0.3048*theta_CurvedWalledNozzle_McNally1970(:,2),'dk','DisplayName','$\theta \, \mathrm{McNally}$')           % momentum

% Predicted by FORTRAN code McNally
plot(0.3048*delta_ast_CurvedNozzle_FORTRAN_CODE_McNally1970(:,1) + INP.L_FP,0.3048*delta_ast_CurvedNozzle_FORTRAN_CODE_McNally1970(:,2),'r-','DisplayName','McNally')
plot(0.3048*theta_CurvedNozzle_FORTRAN_CODE_McNally1970(:,1) + INP.L_FP,0.3048*theta_CurvedNozzle_FORTRAN_CODE_McNally1970(:,2),'r-','HandleVisibility','off')

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
ax1.YLim = [0 3.5e-3];
ax1.YLimMode = 'manual';
ax1.YTickLabelMode = 'manual';
ax1.YTickMode = 'manual';
ax1.YTick = [0 0.5e-3 1.0e-3 1.5e-3 2.0e-3 2.5e-3 3.0e-3 3.5e-3];
ax1.YTickLabel = {0; 0.5; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5};
% ax1.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on

set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest'); % 'Location', 'best');
legend('show')

%% Store plots
% keyboard
NewFolder = cd;
% Figure (complete) for Appendix
saveas(h1,'./Thesis_Report_Figures/NozzleCurvedWall_McNally1970_Complete','epsc')

% Figure for Chapter 4
ax1.XLim = [0.68 1.06];
ax1.XLimMode = 'manual';
ax1.XTickLabelMode = 'manual';
ax1.XTickMode = 'manual';
ax1.XTick = [0.0 + 0.6851 0.05 + 0.6851 0.10 + 0.6851 0.15 + 0.6851 0.20 + 0.6851 0.25 + 0.6851 0.30 + 0.6851 0.35 + 0.6851];
ax1.XTickLabel = {'0'; '0.05'; '0.10'; '0.15'; '0.20'; '0.25'; '0.30'; '0.35'};
ax1.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northeast'); % 'Location', 'best');
legend('show')

saveas(h1,'./Thesis_Report_Figures/NozzleCurvedWall_McNally1970','epsc')

OldFolder = cd(NewFolder);

% fprintf('Figures successfully generated and saved\n')

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


%% Interpolation
% added for Appendix: closer to prediction by McNally, but wavy results
% because of interpolation.

% Interpolation of input
xq = [0:0.001:INP.x(end) INP.x(end)];
p = pchip(INP.x,INP.y,xq);
% q = pchip(INP.x,INP.PsE,xq);
q = pchip(INP.x,INP.psE,xq); % changed to pse on 5/03/2022
INP.x = xq;
INP.y = p;
% INP.PsE = q;
INP.psE = q; % changed to pse on 05/03/2022
INP.BCW = ones(1,length(INP.x));   % [-], Tw/TtI = 1
INP.AWD = [];

% Calculation
[X,NST,HVR,EDG,FRS,FLD,TCC,SOL,BLC,MON] = PRECAL(INP,OPT,SET,FLD,TCC);
[NP,GRD] = GRID(GRD,SET);
NS = 1;
[sol,solprev] = IVPL(NS,GRD,HVR,FLD,SET); %%% solprev = []; (empty, maybe needed for some case? move to PRECAL?)
flp = [];
[BLC,FLP,SOL,TCC,MON] = CSM(X,NP,NS,NST,GRD,HVR,EDG,FRS,FLD,TCC,sol,solprev,SOL,BLC,INP,OPT,SET,MON);

% Plot
h2 = figure;
hold on

% Experimental data Lewis Research Center
plot(0.3048*delta_ast_CurvedWalledNozzle_McNally1970(:,1) + INP.L_FP,0.3048*delta_ast_CurvedWalledNozzle_McNally1970(:,2),'ok','DisplayName','$\delta^* \, \mathrm{McNally}$') % displacement
plot(0.3048*theta_CurvedWalledNozzle_McNally1970(:,1) + INP.L_FP,0.3048*theta_CurvedWalledNozzle_McNally1970(:,2),'dk','DisplayName','$\theta \, \mathrm{McNally}$')           % momentum

% Predicted by FORTRAN code McNally
plot(0.3048*delta_ast_CurvedNozzle_FORTRAN_CODE_McNally1970(:,1) + INP.L_FP,0.3048*delta_ast_CurvedNozzle_FORTRAN_CODE_McNally1970(:,2),'r-','DisplayName','McNally')
plot(0.3048*theta_CurvedNozzle_FORTRAN_CODE_McNally1970(:,1) + INP.L_FP,0.3048*theta_CurvedNozzle_FORTRAN_CODE_McNally1970(:,2),'r-','HandleVisibility','off')

% CSM
plot(X(1:length(SOL)),BLC.delta_ast,'-k','DisplayName','CSM')        % displacement
plot(X(1:length(SOL)),BLC.theta,'-k','HandleVisibility','off')       % momentum

% Figure
ax2 = gca;

ax2.FontSizeMode = 'manual';
ax2.FontSize = FontSize; %%
ax2.FontWeight = 'normal';
ax2.TickLabelInterpreter = 'latex';

ax2.XLabel.Interpreter = 'latex';
ax2.XLabel.FontSize = FontSize; %%
ax2.XLabel.FontWeight = 'normal';
ax2.XLabel.Color = Color;
ax2.XLabel.String = '$X \, [\mathrm{m}]$';
ax2.XLim = [0 1.2];
ax2.XLimMode = 'manual';
ax2.XTickLabelMode = 'manual';
ax2.XTickMode = 'manual';
ax2.XTick = [0 0.2 0.4 0.6 0.8 1.0 1.2];
ax2.XTickLabel = {'0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'; '1.2'};

ax2.YLabel.Interpreter = 'latex';
ax2.YLabel.FontSize = FontSize; %%
ax2.YLabel.FontWeight = 'normal';
ax2.YLabel.Color = Color;
ax2.YLabel.String = '$\delta^{\ast}, \, \theta \, \times \, 10^{-3} [\mathrm{m}]$';
ax2.YLim = [0 3.5e-3];
ax2.YLimMode = 'manual';
ax2.YTickLabelMode = 'manual';
ax2.YTickMode = 'manual';
ax2.YTick = [0 0.5e-3 1.0e-3 1.5e-3 2.0e-3 2.5e-3 3.0e-3 3.5e-3];
ax2.YTickLabel = {0; 0.5; 1.0; 1.5; 2.0; 2.5; 3.0; 3.5};
% ax2.PlotBoxAspectRatio = [1 1 1]; % axis = 'square';

set(gca,'FontSize',11)
box on

set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', 'northwest'); % 'Location', 'best');
legend('show')

% Store plot
% keyboard
NewFolder = cd;
% Figure Interpolated (complete) for Appendix
saveas(h2,'./Thesis_Report_Figures/NozzleCurvedWall_McNally1970_Complete_Interpolated','epsc')
OldFolder = cd(NewFolder);
fprintf('Figures successfully generated and saved\n')
%%