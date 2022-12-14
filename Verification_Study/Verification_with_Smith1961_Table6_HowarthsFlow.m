%% Verification with Howarths Flow - Table 6 from Smith and Clutter (1961) p.59
% Howarths flow: Nonsimilar decelarating flow
% adiabatic incompressible adverse pressure gradient flow with separation
% Nonsimilar solutions
% 18-11-2019 D.D.D.
% last changes on 18-11-2019
% last run on 04-05-2020
% Simulation takes about a second
% Adapted slightly for upload and checked: 23-01-2022

%%% It takes 346.410033 seconds for etaE=9.0; h1=0.001;X=0:0.001:1.0
%%% 'Separation occurred at station 960/1001'
%%% and changing X=0:0.001:1.0 to X=0:0.0001:1.0 gives:
%%% Separation occurred at station 9584/10001 
%%% Elapsed time is 4733.569593 seconds.

%% Note: only for single station calculations

clear all
close all
clc

pause(0.1)

%% INPUT
load('./Data/VerData_Smith1961_Table6_HowarthsFlow')
run('./INPUT/INPUT_Verification_with_Smith1961_Howarth') % changed name, added Howarth (12-02-2020)
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
% keyboard
%% TABLE 6 - Write data to Tex-file (in Table format)

NewFolder = cd; % currently in 'boundary-layer-code'-folder
cd('Thesis_Report_Tables') %%% new, added (23-01-2021)
fid = fopen('Verification_Smith1961_Table6_Howarth.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Verification with Howarths Flow case with Table 6 from Smith\\cite{smith1961solution}: incompressible adiabatic decelarating flow)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Comparison with tabulated data of nonsimilar Howarths Flow: incompressible adiabatic decelarating flow (adverse pressure gradient) taken from Smith\\cite{smith1961solution} table 6. All values in compressible Falkner-Skan transformed $y$-coordinate. Grid used by the CS-method (CSM): uniform (vertical) grid spacing of $\\mathrm{d} \\eta = %1.4f$ and height of $\\eta_{\\mathrm{e}} = %1.1f$. Note that separation occurred when the table entry shows ''sep''.}\n',GRD.Deta(1),GRD.etaE);
fprintf(fid, '    \\label{tab:CT6}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=1.4] S[table-format=1.6] S[table-format=1.5] S[table-format=1.5] S[table-format=1.6]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$X$}                       &\n');
fprintf(fid, '        \\multicolumn{4}{c}{$f''''(0)$}                  \\\\\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}                       &\n');
fprintf(fid, '        \\multicolumn{4}{c}{[-]}                       \\\\\n');
fprintf(fid, '        \\cmidrule(lr){2-5}\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                          &\n');
fprintf(fid, '        \\multicolumn{2}{c}{Douglas}                   &\n');
fprintf(fid, '        \\multicolumn{1}{c}{Hartree}                   &\n');
fprintf(fid, '        \\multicolumn{1}{c}{CSM}                       \\\\\n\n');
fprintf(fid, '        \\cmidrule(lr){2-3}\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                          &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\epsilon = 0.000001$}   &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\epsilon = 0.00001$}    &\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                          &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\epsilon = 0.00001$}                          \\\\\n');
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:length(X)
    kk = ll - 1;
    if ll == 1
        fprintf(fid, '        %s   &        &        &       &  %1.6f   \\\\\n',num2str(X(1)),SOL{ll}.v(1));
    elseif ll >= MON.STR % separation has occurred, print 'sep'
        if kk < 9 || kk == 10 || kk == 18
            fprintf(fid, '        %s   &   %s   &   %s   &       &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3));
        else
            fprintf(fid, '        %s   &   %s   &   %s   &  %s   &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),VerData_Smith1961_Table6_HowarthsFlow(kk,4));
        end
        elseif kk < 9 || kk == 10 || kk == 18
        if ll >= MON.STR % separation has occurred, print 'sep'
            fprintf(fid, '        %s   &   %s   &   %s   &       &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3));
        else
            fprintf(fid, '        %s   &   %s   &   %s   &       &   %1.6f  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),SOL{ll}.v(1));
        end
    else
        fprintf(fid, '        %s   &   %s   &   %s   &  %s   &  %1.6f   \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),VerData_Smith1961_Table6_HowarthsFlow(kk,4),SOL{ll}.v(1));
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

% keyboard

%% PLOT

i = 0;
for j = 1:length(VerData_Smith1961_Table6_HowarthsFlow)
    if ~isnan(str2double(VerData_Smith1961_Table6_HowarthsFlow(j,4)))
        i = i + 1;
    Hartree(i,1) = VerData_Smith1961_Table6_HowarthsFlow(j,1);
    Hartree(i,2) = VerData_Smith1961_Table6_HowarthsFlow(j,2);
    end
end

for i = 1:length(SOL)
    fpp(i) = SOL{i}.v(1);
end

for i = 1:length(SOL)
    Cw(i) = FLP{i}.C(1);
end

Cf_2 = BLC.Cf_L2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI);
% flp.C(1)*sol.v(1)/sqrt(EDG.Re_x(NS))
% EDG.Re_x(i) = EDG.rhoE(i)*EDG.UE(i)*X(i)/EDG.muE(i)
% % % Cw(2:end).*fpp(2:end)./sqrt(EDG.Re_x(2:end-1))*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI);
% Cw=1; everywhere

%%% Thesis figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
Color = [0 0 0];
ColorZoom = [0.65 0.65 0.65];
LineStyle = {'--','-.','none','-'};
Marker = {'none','none','*','none'};
LineWidth = 0.5;
FontSize = 11;
FontSizeZoom = 8;
LegendPos = 'northeast';%[0.20 0.55 0.20 0.20]; % normalized position
% LgndPos2 = [0.55 0.3 0.3 0.2];
LegendName = {'Smith A','Smith B','Hartree','CSM'};

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
ax.XLim = [0 1.0]; % change to 0.15?
ax.XLimMode = 'manual';
ax.XTickLabelMode = 'manual';
ax.XTickMode = 'manual';
ax.XTick = [0 0.2 0.4 0.6 0.8 1.0];
ax.XTickLabel = {'0.0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'};

ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = FontSize; %%
ax.YLabel.FontWeight = 'normal';
ax.YLabel.Color = Color;
ax.YLabel.String = 'Shear parameter $f''''(0) \; [-]$';
ax.YLim = [0 0.35];
ax.YLimMode = 'manual';
ax.YTickLabelMode = 'manual';
ax.YTickMode = 'manual';
ax.YTick = [0 0.05 0.10 0.15 0.20 0.25 0.30 0.35];
ax.YTickLabel = {'0'; '0.05'; '0.10'; '0.15'; '0.20'; '0.25'; '0.30'; '0.35'};
set(gca,'FontSize',FontSize)
box on
hold on

% Data
X1 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1));
X2 = X(1:length(SOL));
Y1 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2));
Y2 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3));
Y3 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4));
Y4 = fpp;

% Plotting
% 'Color',Color,'LineStyle',LineStyle,'Marker',Marker,'DisplayName',LegendName,'HandleVisibility','on'
% % plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2)),'Color',Color,'LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName',LegendName{1},'HandleVisibility','on')
% % plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3)),'Color',Color,'LineStyle',LineStyle{2},'Marker',Marker{2},'DisplayName',LegendName{2},'HandleVisibility','on')
% % plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4)),'Color',Color,'LineStyle',LineStyle{3},'Marker',Marker{3},'DisplayName',LegendName{3},'HandleVisibility','on')
% % plot(X(1:length(SOL)),fpp,'Color',Color,'LineStyle',LineStyle{4},'Marker',Marker{4},'DisplayName',LegendName{4},'HandleVisibility','on') % separation might have occurred

% Zoom box and lines
plot([0.92 0.98 0.98 0.92 0.92],[0 0 0.04 0.04 0],'Color',ColorZoom)
plot([0.54 0.92],[0.188 0.04],'Color',ColorZoom)
plot([0.54 0.92],[0.04 0],'Color',ColorZoom)
% % Zoom box and lines
% plot([0.92 0.96 0.96 0.92 0.92],[0 0 0.05 0.05 0],'Color',ColorZoom)
% plot([0.54 0.92],[0.188 0.05],'Color',ColorZoom)
% plot([0.54 0.92],[0.04 0],'Color',ColorZoom)
% % a = annotation('line',[0.1 0.2],[0.1 0.2]);
% % a.Color = 'red';

plot(X1,Y1,'Color',Color,'LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName',LegendName{1},'HandleVisibility','on')
plot(X1,Y2,'Color',Color,'LineStyle',LineStyle{2},'Marker',Marker{2},'DisplayName',LegendName{2},'HandleVisibility','on')
plot(X1,Y3,'Color',Color,'LineStyle',LineStyle{3},'Marker',Marker{3},'DisplayName',LegendName{3},'HandleVisibility','on')
plot(X2,Y4,'Color',Color,'LineStyle',LineStyle{4},'Marker',Marker{4},'DisplayName',LegendName{4},'HandleVisibility','on')

% Legend
% set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Position', LegendPos); % 'Location', 'best');
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', LegendPos); % 'Location', 'best');
legend('show')

% Zoom inside figure
axes('position',[0.2 0.2 0.35 .35])
% axes('Location', ZoomPlotPos)
axz = gca;

axz.FontSizeMode = 'manual';
axz.FontSize = FontSizeZoom; %%
axz.FontWeight = 'normal';
axz.TickLabelInterpreter = 'latex';

axz.XLim = [0.92 0.96];
axz.XLimMode = 'manual';
axz.XTickLabelMode = 'manual';
axz.XTickMode = 'manual';
axz.XTick = [0.92 0.93 0.94 0.95 0.96];
axz.XTickLabel = {'0.92'; '0.93'; '0.94'; '0.95'; '0.96'};

axz.YLim = [0 0.05];
axz.YLimMode = 'manual';
axz.YTickLabelMode = 'manual';
axz.YTickMode = 'manual';
axz.YTick = [0 0.01 0.02 0.03 0.04 0.05];
axz.YTickLabel = {'0'; '0.01'; '0.02'; '0.03'; '0.04'; '0.05'};
set(gca,'FontSize',FontSizeZoom)
box on
hold on

% X1 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1));
% X2 = X(1:length(SOL));
% Y1 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2));
% Y2 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3));
% Y3 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4));
% Y4 = fpp;

IndexOfInterest1 = (X1 < 0.96) & (X1 > 0.90);
IndexOfInterest2 = (X2 < 0.96) & (X2 > 0.90);

plot(X1(IndexOfInterest1),Y1(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName',LegendName{1},'HandleVisibility','off')
plot(X1(IndexOfInterest1),Y2(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{2},'Marker',Marker{2},'DisplayName',LegendName{2},'HandleVisibility','off')
plot(X1(IndexOfInterest1),Y3(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{3},'Marker',Marker{3},'DisplayName',LegendName{3},'HandleVisibility','off')
plot(X2(IndexOfInterest2),Y4(IndexOfInterest2),'Color',Color,'LineStyle',LineStyle{4},'Marker',Marker{4},'DisplayName',LegendName{4},'HandleVisibility','off')
% % % axis tight % use this command to find a nice coordinate range

% For reference:
% https://nl.mathworks.com/matlabcentral/answers/33779-zooming-a-portion-of-figure-in-a-figure
% indexOfInterest = (t < 11*pi/8) & (t > 9*pi/8); % range of t near perturbation
% plot(t(indexOfInterest),signal(indexOfInterest)) % plot on new axes
% axis tight
%%% Thesis figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Thesis figure2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
Color = [0 0 0];
ColorZoom = [0.65 0.65 0.65];
LineStyle = {'--','-.','none','-'};
Marker = {'none','none','*','none'};
LineWidth = 0.5;
FontSize = 11;
FontSizeZoom = 8;
LegendPos = 'northeast';%[0.20 0.55 0.20 0.20]; % normalized position
% LgndPos2 = [0.55 0.3 0.3 0.2];
LegendName = {'Smith A','Smith B','Hartree','CSM'};

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
ax.XLim = [0 1.0]; % change to 0.15?
ax.XLimMode = 'manual';
ax.XTickLabelMode = 'manual';
ax.XTickMode = 'manual';
ax.XTick = [0 0.2 0.4 0.6 0.8 1.0];
ax.XTickLabel = {'0.0'; '0.2'; '0.4'; '0.6'; '0.8'; '1.0'};

ax.YLabel.Interpreter = 'latex';
ax.YLabel.FontSize = FontSize; %%
ax.YLabel.FontWeight = 'normal';
ax.YLabel.Color = Color;
ax.YLabel.String = 'Shear parameter $f''''(0) \; [-]$';
ax.YLim = [0 0.35];
ax.YLimMode = 'manual';
ax.YTickLabelMode = 'manual';
ax.YTickMode = 'manual';
ax.YTick = [0 0.05 0.10 0.15 0.20 0.25 0.30 0.35];
ax.YTickLabel = {'0'; '0.05'; '0.10'; '0.15'; '0.20'; '0.25'; '0.30'; '0.35'};
set(gca,'FontSize',FontSize)
box on
hold on

% Data
X1 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1));
X2 = X(1:length(SOL));
Y1 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2));
Y2 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3));
Y3 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4));
Y4 = fpp;

% Plotting
% 'Color',Color,'LineStyle',LineStyle,'Marker',Marker,'DisplayName',LegendName,'HandleVisibility','on'
% % plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2)),'Color',Color,'LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName',LegendName{1},'HandleVisibility','on')
% % plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3)),'Color',Color,'LineStyle',LineStyle{2},'Marker',Marker{2},'DisplayName',LegendName{2},'HandleVisibility','on')
% % plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4)),'Color',Color,'LineStyle',LineStyle{3},'Marker',Marker{3},'DisplayName',LegendName{3},'HandleVisibility','on')
% % plot(X(1:length(SOL)),fpp,'Color',Color,'LineStyle',LineStyle{4},'Marker',Marker{4},'DisplayName',LegendName{4},'HandleVisibility','on') % separation might have occurred

% Zoom box and lines
% XZ1 = 0.92;
% XZ2 = 0.96;
% YZ1 = 0;
% YZ2 = 0.05;
plot([0.92 0.98 0.98 0.92 0.92],[0 0 0.04 0.04 0],'Color',ColorZoom,'HandleVisibility','off')
plot([0.35 0.92],[0.23 0.035],'Color',ColorZoom,'HandleVisibility','off') % top line
plot([0.35 0.92],[0.03 0.005],'Color',ColorZoom,'HandleVisibility','off') % bottom line

% % plot([0.92 0.98 0.98 0.92 0.92],[0 0 0.04 0.04 0],'Color',ColorZoom)
% % plot([0.54 0.92],[0.188 0.04],'Color',ColorZoom)
% % plot([0.54 0.92],[0.04 0],'Color',ColorZoom)
% % Zoom box and lines
% % plot([0.92 0.96 0.96 0.92 0.92],[0 0 0.05 0.05 0],'Color',ColorZoom)
% % plot([0.54 0.92],[0.188 0.05],'Color',ColorZoom)
% % plot([0.54 0.92],[0.04 0],'Color',ColorZoom)
% % a = annotation('line',[0.1 0.2],[0.1 0.2]);
% % a.Color = 'red';

plot(X1,Y1,'Color',Color,'LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName',LegendName{1},'HandleVisibility','on')
plot(X1,Y2,'Color',Color,'LineStyle',LineStyle{2},'Marker',Marker{2},'DisplayName',LegendName{2},'HandleVisibility','on')
plot(X1,Y3,'Color',Color,'LineStyle',LineStyle{3},'Marker',Marker{3},'DisplayName',LegendName{3},'HandleVisibility','on')
plot(X2,Y4,'Color',Color,'LineStyle',LineStyle{4},'Marker',Marker{4},'DisplayName',LegendName{4},'HandleVisibility','on')

% Legend
% set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Position', LegendPos); % 'Location', 'best');
set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', LegendPos); % 'Location', 'best');
legend('show')

% Zoom inside figure
% axes('position',[0.2 0.2 0.35 .35])
% dSX = 0.96 - 0.92;
% dSY = 0.05 - 0;
% SY = 0.4;
% SX = SY*dSX/dSY;
dSX = 1.0 - 0.0;
dSY = 0.35 - 0.0;
SY = 0.5;
SX = SY*dSY/dSX;
axes('position',[0.26 0.16 SX SY])
% axes('Location', ZoomPlotPos)
axz = gca;

axz.FontSizeMode = 'manual';
axz.FontSize = FontSizeZoom; %%
axz.FontWeight = 'normal';
axz.TickLabelInterpreter = 'latex';

axz.XLim = [0.92 0.98];
axz.XLimMode = 'manual';
axz.XTickLabelMode = 'manual';
axz.XTickMode = 'manual';
axz.XTick = [0.92 0.94 0.96 0.98];
axz.XTickLabel = {'0.92'; '0.94'; '0.96'; '0.98'};

axz.YLim = [0 0.04];
axz.YLimMode = 'manual';
axz.YTickLabelMode = 'manual';
axz.YTickMode = 'manual';
axz.YTick = [0 0.01 0.02 0.03 0.04 0.05];
axz.YTickLabel = {'0'; '0.01'; '0.02'; '0.03'; '0.04'; '0.05'};
set(gca,'FontSize',FontSizeZoom)
box on
hold on

% X1 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1));
% X2 = X(1:length(SOL));
% Y1 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2));
% Y2 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3));
% Y3 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4));
% Y4 = fpp;

IndexOfInterest1 = (X1 < 0.96) & (X1 > 0.90);
IndexOfInterest2 = (X2 < 0.96) & (X2 > 0.90);

plot(X1(IndexOfInterest1),Y1(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName',LegendName{1},'HandleVisibility','off')
plot(X1(IndexOfInterest1),Y2(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{2},'Marker',Marker{2},'DisplayName',LegendName{2},'HandleVisibility','off')
plot(X1(IndexOfInterest1),Y3(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{3},'Marker',Marker{3},'DisplayName',LegendName{3},'HandleVisibility','off')
plot(X2(IndexOfInterest2),Y4(IndexOfInterest2),'Color',Color,'LineStyle',LineStyle{4},'Marker',Marker{4},'DisplayName',LegendName{4},'HandleVisibility','off')
% % % axis tight % use this command to find a nice coordinate range

% Store plot
% NewFolder = cd;
% cd('Thesis_Report_Tables')

saveas(h,'Thesis_Report_Figures/Clutter1961_Howarth','epsc')

% OldFolder = cd(NewFolder);

% For reference:
% https://nl.mathworks.com/matlabcentral/answers/33779-zooming-a-portion-of-figure-in-a-figure
% indexOfInterest = (t < 11*pi/8) & (t > 9*pi/8); % range of t near perturbation
% plot(t(indexOfInterest),signal(indexOfInterest)) % plot on new axes
% axis tight
%%% Thesis figure2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keyboard

figure
hold on
title('Howarth''s Flow')
xlabel('Dimensionless $X$-coordinate [-]')
ylabel('shear parameter []')
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2)),'k')
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3)),'g')
% plot(str2double(Hartree(:,1)),str2double(Hartree(:,2)))
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(9,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(9,4)),'m*','HandleVisibility','off')
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4)),'m')
plot(X(1:length(SOL)),fpp,'b') % separation might have occurred
% legend('Smith 1','Smith 2','Hartree','Hartree','CSM')
legend('Smith 1','Smith 2','Hartree','CSM')

% figure
% hold on
% title('Howarths Flow')
% xlabel('Dimensionless X-coordinate [-]')
% ylabel('shear parameter []')
% plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2))/FRS.rhoI/FRS.UI^2,'k')
% plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3))/FRS.rhoI/FRS.UI^2,'g')
% % plot(str2double(Hartree(:,1)),str2double(Hartree(:,2)))
% plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(9,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(9,4))/FRS.rhoI/FRS.UI^2,'m*','HandleVisibility','off')
% plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4))/FRS.rhoI/FRS.UI^2,'m')
% plot(X(1:length(SOL)),fpp/FRS.rhoI/FRS.UI^2,'b') % separation might have occurred
% % legend('Smith 1','Smith 2','Hartree','Hartree','CSM')
% legend('Smith 1','Smith 2','Hartree','CSM')

%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
title('Howarth''s Flow')
xlabel('Dimensionless $X$-coordinate [-]')
ylabel('shear parameter []')
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2))'./sqrt(EDG.Re_x(2:end))*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'k')
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3))'./sqrt(EDG.Re_x(2:end))*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'g')
% plot(str2double(Hartree(:,1)),str2double(Hartree(:,2)))
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4))'./sqrt(EDG.Re_x(2:end))*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'m*')
% plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4)),'m')
plot(X(1:length(SOL)),Cf_2,'b') % separation might have occurred
% legend('Smith 1','Smith 2','Hartree','Hartree','CSM')
legend('Smith 1','Smith 2','Hartree','CSM')

%%% TEST2
% flp.C(1)*sol.v(1)/sqrt(EDG.Re_x(NS)) %%% NB Cw=1
% *EDG.rhoE(NS)/FRS.rhoI*(EDG.UE(NS)/FRS.UI)^2
% *sqrt(FRS.UI*1/FRS.muI*FRS.rhoI)
% ./sqrt(EDG.Re_x(2:end)).*EDG.rhoE(2:end)/FRS.rhoI.*(EDG.UE(2:end)/FRS.UI).^2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI)
figure
hold on
title('Howarth''s Flow')
xlabel('Dimensionless $X$-coordinate [-]')
ylabel('shear parameter []')
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,2))'./sqrt(EDG.Re_x(2:end)).*EDG.rhoE(2:end)/FRS.rhoI.*(EDG.UE(2:end)/FRS.UI).^2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'k')
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,3))'./sqrt(EDG.Re_x(2:end)).*EDG.rhoE(2:end)/FRS.rhoI.*(EDG.UE(2:end)/FRS.UI).^2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'g')
% plot(str2double(Hartree(:,1)),str2double(Hartree(:,2)))
plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4))'./sqrt(EDG.Re_x(2:end)).*EDG.rhoE(2:end)/FRS.rhoI.*(EDG.UE(2:end)/FRS.UI).^2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI),'m*')
% plot(str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1)),str2double(VerData_Smith1961_Table6_HowarthsFlow(:,4)),'m')
plot(X(1:length(SOL)),Cf_2,'b') % separation might have occurred
% legend('Smith 1','Smith 2','Hartree','Hartree','CSM')
legend('Smith 1','Smith 2','Hartree','CSM')

PLOTFILE
