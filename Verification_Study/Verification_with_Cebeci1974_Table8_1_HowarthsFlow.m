%% Verification with Howarths Flow - Table 6 from Smith and Clutter (1961) p.59
% Howarths flow: Nonsimilar decelarating flow
% adiabatic incompressible adverse pressure gradient flow with separation
% Nonsimilar solutions
% 19-11-2019 D.D.D.
% last changes on 19-11-2019
% Last run on 04-05-2020
% Simulation takes about ? second
% Adapted slightly for upload and checked: 23-01-2022

%% Note: only for single station calculations

clear all
close all
clc

pause(0.1)

%% INPUT
load('./Data/VerData_Cebeci1974_Table8_1_HowarthsFlow')
run('./INPUT/INPUT_Verification_with_Cebeci1974_Howarth')
cd ..

SET.ITMAX0 = 20;    %%% added (23-01-2022)
OPT.RLAM = 0;       %%% added (23-01-2022)

%% Calculation
% Boundary layer edge calculations
[X,NST,HVR,EDG,FRS,FLD,TCC,SOL,BLC,MON] = PRECAL(INP,OPT,SET,FLD,TCC);

% Grid generation
[NP,GRD] = GRID(GRD,SET);

NS = 1;

% IVPL
[sol,solprev] = IVPL(NS,GRD,HVR,FLD,SET); %%% solprev = []; (empty, maybe needed for some case? move to PRECAL?)
flp = [];

tic % comment in future?

[BLC,FLP,SOL,TCC,MON] = CSM(X,NP,NS,NST,GRD,HVR,EDG,FRS,FLD,TCC,sol,solprev,SOL,BLC,INP,OPT,SET,MON);

toc % comment in future?

fprintf('Simulation finished with success \n')
% PLOTFILE

% Cf_2 = BLC.Cf_L2*sqrt(EDG.UE(1)*1/EDG.muE(1)*EDG.rhoE(1));
Cf_2 = BLC.Cf_L2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI);

% keyboard

%% TABLE 6 - Write data to Tex-file (in Table format)

NewFolder = cd; % currently in 'boundary-layer-code'-folder
cd('Thesis_Report_Tables') % (23-01-2021)
fid = fopen('Verification_Cebeci1974_Table8_1_Howarth.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Verification with Howarths Flow case with Table 8-1 from Cebeci\\cite{cebeci1974analysis}: incompressible adiabatic decelarating flow)\n');
fprintf(fid, '%%\\begin{table}\n');
fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Comparison with tabulated data of nonsimilar Howarths Flow: incompressible adiabatic decelarating flow (adverse pressure gradient) taken from Cebeci\\cite{cebeci1974analysis} table 8-1. All values in compressible Falkner-Skan transformed y-coordinate. Grid used by the CS-method (CSM): uniform (vertical) grid spacing of $\\mathrm{d} \\eta = %1.4f$ and height of $\\eta_{\\mathrm{e}} = %1.1f$. Note that separation occurred when the table entry shows ''sep''.}\n',GRD.Deta(1),GRD.etaE,GRD.Deta(1),GRD.etaE);
fprintf(fid, '    \\label{tab:CT8_1}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=1.4] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\bar{X}$}                               &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$f''''(0)$}                                &\n');
fprintf(fid, '        \\multicolumn{5}{c}{$\\frac{\\tau(0)}{\\rho u_{\\infty}^2}$}     \\\\\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}                                     &\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}                                     &\n');
fprintf(fid, '        \\multicolumn{5}{c}{[-]}                                     \\\\\n');
fprintf(fid, '        \\cmidrule(lr){3-7}\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                        &\n');
fprintf(fid, '        \\multicolumn{2}{c}{CSM}                                     &\n');
fprintf(fid, '        \\multicolumn{1}{c}{Cebeci}                                  &\n');
fprintf(fid, '        \\multicolumn{1}{c}{Howarth}                                 &\n');
fprintf(fid, '        \\multicolumn{1}{c}{Smith-Clutter}                           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{G\\"{o}rtler}                             \\\\\n');
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:length(X)
    kk = ll - 1;
    if ll == 1
        fprintf(fid, '        %s   &   %1.6f   &   %1.6f   &      &        &       &   \\\\\n',num2str(X(1)),SOL{ll}.v(1),Cf_2(ll));
    elseif ll >= MON.STR % separation has occurred, print 'sep'
        fprintf(fid, '        %s   &   \\multicolumn{1}{c}{sep}   &   \\multicolumn{1}{c}{sep}   &   %s   &   %s   &   %s   &   %s   \\\\\n',VerData_Cebeci1974_Table8_1_HowarthsFlow(kk,1),VerData_Cebeci1974_Table8_1_HowarthsFlow(kk,2),VerData_Cebeci1974_Table8_1_HowarthsFlow(kk,3),VerData_Cebeci1974_Table8_1_HowarthsFlow(kk,4),VerData_Cebeci1974_Table8_1_HowarthsFlow(kk,5));
    else
        fprintf(fid, '        %s   &   %1.6f   &   %1.6f   &   %s   &   %s   &  %s   &  %s   \\\\\n',VerData_Cebeci1974_Table8_1_HowarthsFlow(kk,1),SOL{ll}.v(1),Cf_2(ll),VerData_Cebeci1974_Table8_1_HowarthsFlow(kk,2),VerData_Cebeci1974_Table8_1_HowarthsFlow(kk,3),VerData_Cebeci1974_Table8_1_HowarthsFlow(kk,4),VerData_Cebeci1974_Table8_1_HowarthsFlow(kk,5));
    end
end
% fprintf(fid, '        %% Data\n');
fprintf(fid, '        \\bottomrule\n\n');
fprintf(fid, '    \\end{tabular}\n');
% fprintf(fid, '    \\begin{tablenotes}\n');
% fprintf(fid, '        \\item[*] The solution does not converge for $m_2 = \\infty$, and therefore $m_2 = 1000$ is taken as approach.\n');
% fprintf(fid, '    \\end{tablenotes}\n');
fprintf(fid, '\\end{threeparttable}\n');
% fprintf(fid, '%%\\end{sidewaystable}\n');
% fprintf(fid, '%%\\end{table\n');
fclose(fid);
OldFolder = cd(NewFolder);

%% PLOT

i = 0;
for j = 1:length(VerData_Cebeci1974_Table8_1_HowarthsFlow)
    if ~isnan(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(j,4)))
        i = i + 1;
    Hartree(i,1) = VerData_Cebeci1974_Table8_1_HowarthsFlow(j,1);
    Hartree(i,2) = VerData_Cebeci1974_Table8_1_HowarthsFlow(j,2);
    end
end

for i = 1:length(SOL)
    fpp(i) = SOL{i}.v(1);
end

%% Thesis figure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Thesis figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To do:
% add year to legend (of source). also for 1961; and try to put all in one
% graph again (should work, since same CSM calculation)
% add grey lines, plotted first, to indicate zoom? with circle instead of
% rectangle?

% Settings
Color = [0 0 0];
LineStyle = {'--','none','none','-.','-'};
Marker = {'none','*','+','none','none'};
LineWidth = 0.5;
FontSize = 11;
FontSizeZoom = 8;
LegendPos = 'southwest';%[0.20 0.55 0.20 0.20]; % normalized position
LegendName = {'Cebeci','Howarth','Clutter','G\"{o}rtler','CSM'}; % Clutter (1963) is first author of this report!

% Figure
h = figure;
ax = gca;

ax.FontSizeMode = 'manual';
ax.FontSize = FontSize; %%
ax.FontWeight = 'normal';
ax.TickLabelInterpreter = 'latex';

ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = FontSize; %%
ax.XLabel.FontWeight = 'normal';
ax.XLabel.Color = Color;
ax.XLabel.String = '$X$-coordinate [-]';
ax.XLim = [0 1.0];
ax.XLimMode = 'manual';
ax.XTickLabelMode = 'manual';
ax.XTickMode = 'manual';
ax.XTick = [0 0.2 0.4 0.6 0.8 1.0];
ax.XTickLabel = {'0.0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'};

ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = FontSize; %%
ax.YLabel.FontWeight = 'normal';
ax.YLabel.Color = Color;
ax.YLabel.String = '$C_{\mathrm{f}} \; [-]$';
ax.YLim = [0 1.0];
ax.YLimMode = 'manual';
ax.YTickLabelMode = 'manual';
ax.YTickMode = 'manual';
ax.YTick = [0 0.2 0.4 0.6 0.8 1.0];
ax.YTickLabel = {'0.0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'};
set(gca,'FontSize',FontSize)
box on
hold on

% Data
X1 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1));
X2 = X(1:length(SOL));
Y1 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,2));
Y2 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,3));
Y3 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,4));
Y4 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,5));
Y5 = BLC.Cf2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI);

% Plotting
plot(X1,Y1,'Color',Color,'LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName',LegendName{1},'HandleVisibility','on')
plot(X1,Y2,'Color',Color,'LineStyle',LineStyle{2},'Marker',Marker{2},'DisplayName',LegendName{2},'HandleVisibility','on')
plot(X1,Y3,'Color',Color,'LineStyle',LineStyle{3},'Marker',Marker{3},'DisplayName',LegendName{3},'HandleVisibility','on')
plot(X1,Y4,'Color',Color,'LineStyle',LineStyle{4},'Marker',Marker{4},'DisplayName',LegendName{4},'HandleVisibility','on')
plot(X2,Y5,'Color',Color,'LineStyle',LineStyle{5},'Marker',Marker{5},'DisplayName',LegendName{5},'HandleVisibility','on')

% Legend
% set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Position', LegendPos); % 'Location', 'best');
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', LegendPos); % 'Location', 'best');
legend('show')

% Zoom inside figure
% x = [0.88 0.96 0.96 0.88 0.88]; y = [0.0 0.0 0.10 0.10 0.0]; plot(x,y,'r', 'LineWidth',LineWidth)
% % a = annotation('line',[0.1 0.2],[0.1 0.2]);
% % a.Color = 'red';
% annotation('line',[0.5 0.8],[0.5 0.15],'Color',[0.50 0.50 0.50])
% annotation('line',[0.85 0.9],[0.5 0.15],'Color',[0.50 0.50 0.50])
axes('position',[0.5 0.5 0.35 .35])
% axes('Location', ZoomPlotPos)
axz = gca;

axz.FontSizeMode = 'manual';
axz.FontSize = FontSizeZoom; %%
axz.FontWeight = 'normal';
axz.TickLabelInterpreter = 'latex';

axz.XLim = [0.88 0.96];
axz.XLimMode = 'manual';
axz.XTickLabelMode = 'manual';
axz.XTickMode = 'manual';
axz.XTick = [0.88 0.90 0.92 0.94 0.96];
axz.XTickLabel = {'0.88'; '0.90'; '0.92'; '0.94'; '0.96'};

axz.YLim = [0 0.10];
axz.YLimMode = 'manual';
axz.YTickLabelMode = 'manual';
axz.YTickMode = 'manual';
axz.YTick = [0 0.02 0.04 0.06 0.08 0.10];
axz.YTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};
set(gca,'FontSize',FontSizeZoom)
box on
hold on

% X1 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1));
% X2 = X(1:length(SOL));
% Y1 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2));
% Y2 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3));
% Y3 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4));
% Y4 = fpp;

IndexOfInterest1 = (X1 < 0.97) & (X1 > 0.87);
IndexOfInterest2 = (X2 < 0.97) & (X2 > 0.87);

plot(X1(IndexOfInterest1),Y1(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName',LegendName{1},'HandleVisibility','off')
plot(X1(IndexOfInterest1),Y2(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{2},'Marker',Marker{2},'DisplayName',LegendName{2},'HandleVisibility','off')
plot(X1(IndexOfInterest1),Y3(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{3},'Marker',Marker{3},'DisplayName',LegendName{3},'HandleVisibility','off')
plot(X1(IndexOfInterest1),Y4(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{4},'Marker',Marker{4},'DisplayName',LegendName{4},'HandleVisibility','off')
plot(X2(IndexOfInterest2),Y5(IndexOfInterest2),'Color',Color,'LineStyle',LineStyle{5},'Marker',Marker{5},'DisplayName',LegendName{5},'HandleVisibility','off')
% % % axis tight % use this command to find a nice coordinate range

% For reference:
% https://nl.mathworks.com/matlabcentral/answers/33779-zooming-a-portion-of-figure-in-a-figure
% indexOfInterest = (t < 11*pi/8) & (t > 9*pi/8); % range of t near perturbation
% plot(t(indexOfInterest),signal(indexOfInterest)) % plot on new axes
% axis tight
%%% Thesis figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Thesis figure 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Thesis figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Thesis figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To do:
% add year to legend (of source). also for 1961; and try to put all in one
% graph again (should work, since same CSM calculation)
% add grey lines, plotted first, to indicate zoom? with circle instead of
% rectangle?

% Settings
Color = [0 0 0];
ColorZoom = [0.65 0.65 0.65];
LineStyle = {'--','none','none','-.','-'};
Marker = {'none','*','+','none','none'};
LineWidth = 0.5;
FontSize = 11;
FontSizeZoom = 8;
LegendPos = 'southwest';%[0.20 0.55 0.20 0.20]; % normalized position
LegendName = {'Cebeci','Howarth','Clutter','G\"{o}rtler','CSM'}; % Clutter (1963) is first author of this report!

% Figure
h = figure;
axis square
ax = gca;

ax.FontSizeMode = 'manual';
ax.FontSize = FontSize; %%
ax.FontWeight = 'normal';
ax.TickLabelInterpreter = 'latex';

ax.XLabel.Interpreter = 'latex';
ax.XLabel.FontSize = FontSize; %%
ax.XLabel.FontWeight = 'normal';
ax.XLabel.Color = Color;
ax.XLabel.String = '$X$-coordinate [-]';
ax.XLim = [0 1.0];
ax.XLimMode = 'manual';
ax.XTickLabelMode = 'manual';
ax.XTickMode = 'manual';
ax.XTick = [0 0.2 0.4 0.6 0.8 1.0];
ax.XTickLabel = {'0.0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'};

ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = FontSize; %%
ax.YLabel.FontWeight = 'normal';
ax.YLabel.Color = Color;
ax.YLabel.String = '$C_{\mathrm{f}} \; [-]$';
ax.YLim = [0 1.0];
ax.YLimMode = 'manual';
ax.YTickLabelMode = 'manual';
ax.YTickMode = 'manual';
ax.YTick = [0 0.2 0.4 0.6 0.8 1.0];
ax.YTickLabel = {'0.0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'};
set(gca,'FontSize',FontSize)
box on
hold on

% Data
X1 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1));
X2 = X(1:length(SOL));
Y1 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,2));
Y2 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,3));
Y3 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,4));
Y4 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,5));
Y5 = BLC.Cf2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI);

% Zoom box and lines
% XZ1 = 0.92;
% XZ2 = 0.96;
% YZ1 = 0;
% YZ2 = 0.05;
plot([0.88 0.98 0.98 0.88 0.88],[0 0 0.10 0.10 0],'Color',ColorZoom,'HandleVisibility','off')
plot([0.89 0.45],[0.10 0.50],'Color',ColorZoom,'HandleVisibility','off') % left line
plot([0.97 0.96],[0.10 0.50],'Color',ColorZoom,'HandleVisibility','off') % right line

% Plotting
plot(X1,Y1,'Color',Color,'LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName',LegendName{1},'HandleVisibility','on')
plot(X1,Y2,'Color',Color,'LineStyle',LineStyle{2},'Marker',Marker{2},'DisplayName',LegendName{2},'HandleVisibility','on')
plot(X1,Y3,'Color',Color,'LineStyle',LineStyle{3},'Marker',Marker{3},'DisplayName',LegendName{3},'HandleVisibility','on')
plot(X1,Y4,'Color',Color,'LineStyle',LineStyle{4},'Marker',Marker{4},'DisplayName',LegendName{4},'HandleVisibility','on')
plot(X2,Y5,'Color',Color,'LineStyle',LineStyle{5},'Marker',Marker{5},'DisplayName',LegendName{5},'HandleVisibility','on')

% Legend
% set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Position', LegendPos); % 'Location', 'best');
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', LegendPos); % 'Location', 'best');
legend('show')

% Zoom inside figure
dSX = 1.0 - 0.0;
dSY = 1.0 - 0.0;
SY = 0.35;
SX = SY*dSY/dSX;
axes('position',[0.45 0.55 SX SY])
% axes('Location', ZoomPlotPos)
axz = gca;

axz.FontSizeMode = 'manual';
axz.FontSize = FontSizeZoom; %%
axz.FontWeight = 'normal';
axz.TickLabelInterpreter = 'latex';

axz.XLim = [0.88 0.98];
axz.XLimMode = 'manual';
axz.XTickLabelMode = 'manual';
axz.XTickMode = 'manual';
axz.XTick = [0.88 0.90 0.92 0.94 0.96 0.98];
axz.XTickLabel = {'0.88'; '0.90'; '0.92'; '0.94'; '0.96'; '0.98'};

axz.YLim = [0 0.10];
axz.YLimMode = 'manual';
axz.YTickLabelMode = 'manual';
axz.YTickMode = 'manual';
axz.YTick = [0 0.02 0.04 0.06 0.08 0.10];
axz.YTickLabel = {'0'; '0.02'; '0.04'; '0.06'; '0.08'; '0.10'};
set(gca,'FontSize',FontSizeZoom)
box on
hold on

% X1 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1));
% X2 = X(1:length(SOL));
% Y1 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2));
% Y2 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3));
% Y3 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4));
% Y4 = fpp;

IndexOfInterest1 = (X1 < 0.97) & (X1 > 0.87);
IndexOfInterest2 = (X2 < 0.97) & (X2 > 0.87);

plot(X1(IndexOfInterest1),Y1(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName',LegendName{1},'HandleVisibility','off')
plot(X1(IndexOfInterest1),Y2(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{2},'Marker',Marker{2},'DisplayName',LegendName{2},'HandleVisibility','off')
plot(X1(IndexOfInterest1),Y3(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{3},'Marker',Marker{3},'DisplayName',LegendName{3},'HandleVisibility','off')
plot(X1(IndexOfInterest1),Y4(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{4},'Marker',Marker{4},'DisplayName',LegendName{4},'HandleVisibility','off')
plot(X2(IndexOfInterest2),Y5(IndexOfInterest2),'Color',Color,'LineStyle',LineStyle{5},'Marker',Marker{5},'DisplayName',LegendName{5},'HandleVisibility','off')
% % % axis tight % use this command to find a nice coordinate range

% Store plot
% NewFolder = cd;

saveas(h,'Thesis_Report_Figures/Cebeci1974_Howarth','epsc')

% OldFolder = cd(NewFolder);

% For reference:
% https://nl.mathworks.com/matlabcentral/answers/33779-zooming-a-portion-of-figure-in-a-figure
% indexOfInterest = (t < 11*pi/8) & (t > 9*pi/8); % range of t near perturbation
% plot(t(indexOfInterest),signal(indexOfInterest)) % plot on new axes
% axis tight
%%% Thesis figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Thesis figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
title('Howarth''s Flow')
xlabel('X [-]')
ylabel('shear parameter []')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,2)),'k')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,3)),'g')
% plot(str2double(Hartree(:,1)),str2double(Hartree(:,2)))
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(9,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(9,4)),'m*','HandleVisibility','off')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,4)),'m')
plot(X(1:length(SOL)),fpp,'b') % separation might have occurred
% legend('Smith 1','Smith 2','Hartree','Hartree','CSM')
legend('Smith 1','Smith 2','Hartree','CSM')

figure
title('Howarth''s Flow')
xlabel('X [-]')
ylabel('Cf [-]')
hold on
% plot(X(1:length(BLC.Cf)),Cf_2)
plot(X(1:length(BLC.Cf)),BLC.Cf_L)
plot(X(1:length(BLC.Cf)),BLC.Cf_L2)
plot(X(1:length(BLC.Cf)),BLC.Cf)
plot(X(1:length(BLC.Cf)),BLC.Cf2)
% legend('Cf_2','Cf_L','Cf_L2','Cf','Cf2')
legend('Cf_L','Cf_L2','Cf','Cf2')

%%% Option 1: a few points not plotted due to empty spaces %%%%%%%%%%%%%%%%
figure
hold on
title('Howarth''s Flow')
xlabel('X [-]')
ylabel('Cf [-]')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,2)),'k')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,3)),'g')
% plot(str2double(Hartree(:,1)),str2double(Hartree(:,2)))
% plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(9,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(9,4)),'m*','HandleVisibility','off')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,4)),'m')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,5)),'c')
plot(X(1:length(SOL)),BLC.Cf_L2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'b') % separation might have occurred
plot(X(1:length(SOL)),BLC.Cf2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'r') % separation might have occurred
% legend('Smith 1','Smith 2','Hartree','Hartree','CSM')
legend('Cebeci','Howarth','Smith-Clutter','Gortler','CSM Cf_L/2','CSM Cf/2')

%%% Option 2: plot points instead of lines %%%%%%%%%%%%%%%%
figure
hold on
title('Howarth''s Flow')
xlabel('X [-]')
ylabel('Cf [-]')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,2)),'k')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,3)),'g*')
% plot(str2double(Hartree(:,1)),str2double(Hartree(:,2)))
% plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(9,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(9,4)),'m*','HandleVisibility','off')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,4)),'m+')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,5)),'c')
plot(X(1:length(SOL)),BLC.Cf_L2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'b') % separation might have occurred
plot(X(1:length(SOL)),BLC.Cf2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'r') % separation might have occurred
% legend('Smith 1','Smith 2','Hartree','Hartree','CSM')
legend('Cebeci','Howarth','Smith-Clutter','Gortler','CSM Cf_L/2','CSM Cf/2')

%%% Option 3: combine plotted points with table from Smith (1961) %%%%%%%%%%%%%%%%
load('Verification_Study/Data/VerData_Smith1961_Table6_HowarthsFlow')
figure
hold on
title('Howarth''s Flow')
xlabel('X [-]')
ylabel('Cf [-]')
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2)),'k')
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3)),'g')
% plot(str2double(Hartree(:,1)),str2double(Hartree(:,2)))
% plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(9,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(9,4)),'m*','HandleVisibility','off')
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4)),'mo')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,2))/sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'k')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,3))/sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'g*')
% plot(str2double(Hartree(:,1)),str2double(Hartree(:,2)))
% plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(9,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(9,4)),'m*','HandleVisibility','off')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,4))/sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'m+')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,5))/sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'c')
% plot(X(1:length(SOL)),BLC.Cf_L2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'b') % separation might have occurred
plot(X(1:length(SOL)),BLC.Cf2,'r') % separation might have occurred
% legend('Smith 1','Smith 2','Hartree','Hartree','CSM')
legend('Smith 1','Smith 2','Hartree','Cebeci','Howarth','Smith-Clutter','Gortler','CSM Cf/2')

%% Test
for j = 1:length(SOL)
    Cw(j) = FLP{j}.C(1);
end
%%% Option 3: combine plotted points with table from Smith (1961) %%%%%%%%%%%%%%%%
% load('VerData_Smith1961_Table6_HowarthsFlow')
figure
hold on
title('Howarth''s Flow')
xlabel('X [-]')
ylabel('Cf [-]')
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2)),'k')
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3)),'g')
% plot(str2double(Hartree(:,1)),str2double(Hartree(:,2)))
% plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(9,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(9,4)),'m*','HandleVisibility','off')
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4)),'mo')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(1:14,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(1:14,2))'./sqrt(FRS.UI*1/FRS.muI*FRS.rhoI).*sqrt(EDG.rhoE(2:15).*EDG.UE(2:15).*str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(1:14,1))'./EDG.muE(2:15))./Cw./(EDG.rhoE(2:15)/FRS.rhoI).*(EDG.UE(2:15)/FRS.UI).^2,'k')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,3))/sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'g*')
% plot(str2double(Hartree(:,1)),str2double(Hartree(:,2)))
% plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(9,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(9,4)),'m*','HandleVisibility','off')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,4))/sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'m+')
plot(str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1)),str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,5))/sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'c')
% plot(X(1:length(SOL)),BLC.Cf_L2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'b') % separation might have occurred
plot(X(1:length(SOL)),BLC.Cf2,'r') % separation might have occurred
% legend('Smith 1','Smith 2','Hartree','Hartree','CSM')
legend('Smith 1','Smith 2','Hartree','Cebeci','Howarth','Smith-Clutter','Gortler','CSM Cf/2')

% %%%
% load('VerData_Smith1961_Table6_HowarthsFlow')
% 
% figure
% hold on
% title('Howarth''s Flow')
% xlabel('Dimensionless X-coordinate [-]')
% ylabel('shear parameter []')
% plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2))'./sqrt(EDG.Re_x(2:end)).*EDG.rhoE(2:end)/FRS.rhoI.*(EDG.UE(2:end)/FRS.UI).^2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'k')
% plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3))'./sqrt(EDG.Re_x(2:end)).*EDG.rhoE(2:end)/FRS.rhoI.*(EDG.UE(2:end)/FRS.UI).^2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'g')
% % plot(str2double(Hartree(:,1)),str2double(Hartree(:,2)))
% plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4))'./sqrt(EDG.Re_x(2:end)).*EDG.rhoE(2:end)/FRS.rhoI.*(EDG.UE(2:end)/FRS.UI).^2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'m*')
% % plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4)),'m')
% plot(X(1:length(SOL)),Cf_2,'b') % separation might have occurred
% % legend('Smith 1','Smith 2','Hartree','Hartree','CSM')
% legend('Smith 1','Smith 2','Hartree','CSM')

PLOTFILE
