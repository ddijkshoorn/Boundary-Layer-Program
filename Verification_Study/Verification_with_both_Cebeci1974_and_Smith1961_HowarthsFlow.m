%% Verification with Howarths Flow - Table 6 from Smith and Clutter (1961) p.59
% Howarths flow: Nonsimilar decelarating flow
% adiabatic incompressible adverse pressure gradient flow with separation
% Nonsimilar solutions
% 28-01-2020 D.D.D.
% last changes on 30-01-2020
% last run on 04-05-2020
% Simulation takes about 1 second

% put all in one graph again

%% Note: only for single station calculations

clear all
close all
clc

pause(0.1)

%% First calculation: Cebeci 1974 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT
load('./Data/VerData_Cebeci1974_Table8_1_HowarthsFlow')
run('./INPUT/INPUT_Verification_with_Cebeci1974_Howarth')
cd ..

% Calculation
% Boundary layer edge calculations
[X,NST,HVR,EDG,FRS,FLD,TCC,SOL,BLC,MON] = PRECAL(INP,OPT,SET,FLD,TCC);

% Grid generation
[NP,GRD] = GRID(GRD,SET);

NS = 1;

% IVPL
[sol,solprev] = IVPL(NS,GRD,HVR,FLD); %%% solprev = []; (empty, maybe needed for some case? move to PRECAL?)
flp = [];

tic % comment in future?

[BLC,FLP,SOL,TCC,MON] = CSM(X,NP,NS,NST,GRD,HVR,EDG,FRS,FLD,TCC,sol,solprev,SOL,BLC,INP,OPT,SET,MON);

toc % comment in future?

fprintf('Simulation finished with success \n')

% PLOTFILE

% Cf_2 = BLC.Cf_L2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI); % for reference

% Data
X1 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,1));
X2 = X(1:length(SOL));
Y1 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,2));
Y2 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,3));
Y3 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,4));
Y4 = str2double(VerData_Cebeci1974_Table8_1_HowarthsFlow(:,5));
Y5 = BLC.Cf2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI);

%% Second simulation: Clutter 1961 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all
clear BLC FRS FLD MON EDG GRD HVR INP NP NST NS OPT SET SOL TCC X
% close all
clc
pause(0.1)

% INPUT
load('./Verification_Study/Data/VerData_Smith1961_Table6_HowarthsFlow')
run('./Verification_Study/INPUT/INPUT_Verification_with_Smith1961_Howarth') % changed name, added Howarth (12-02-2020)
cd ..

% Calculation
% Boundary layer edge calculations
[X,NST,HVR,EDG,FRS,FLD,TCC,SOL,BLC,MON] = PRECAL(INP,OPT,SET,FLD,TCC);

% Grid generation
[NP,GRD] = GRID(GRD,SET);

NS = 1;

% IVPL
[sol,solprev] = IVPL(NS,GRD,HVR,FLD); %%% solprev = []; (empty, maybe needed for some case? move to PRECAL?)
flp = [];

tic % comment in future?

[BLC,FLP,SOL,TCC,MON] = CSM(X,NP,NS,NST,GRD,HVR,EDG,FRS,FLD,TCC,sol,solprev,SOL,BLC,INP,OPT,SET,MON);

toc % comment in future?

fprintf('Simulation finished with success \n')

% Cf_2 = BLC.Cf_L2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI); % for reference

% Data
for i = 1:length(FLP)
    Cw(i) = FLP{i}.C(1);
    solv(i) = SOL{i}.v(1);
end

X3 = str2double(VerData_Smith1961_Table6_HowarthsFlow(:,1));
X4 = X(1:length(SOL));
Y6 = str2double(VerData_Smith1961_Table6_HowarthsFlow(1:end-1,2)).*Cw(2:end)'./sqrt(EDG.Re_x(2:length(Cw))').*EDG.rhoE(2:length(Cw))'/FRS.rhoI.*(EDG.UE(2:length(Cw))'/FRS.UI).^2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI);
Y7 = str2double(VerData_Smith1961_Table6_HowarthsFlow(1:end-1,3)).*Cw(2:end)'./sqrt(EDG.Re_x(2:length(Cw))').*EDG.rhoE(2:length(Cw))'/FRS.rhoI.*(EDG.UE(2:length(Cw))'/FRS.UI).^2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI);
Y8 = str2double(VerData_Smith1961_Table6_HowarthsFlow(1:end-1,4)).*Cw(2:end)'./sqrt(EDG.Re_x(2:length(Cw))').*EDG.rhoE(2:length(Cw))'/FRS.rhoI.*(EDG.UE(2:length(Cw))'/FRS.UI).^2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI);
Y9 = solv(2:end)'.*Cw(2:end)'./sqrt(EDG.Re_x(2:length(Cw))').*EDG.rhoE(2:length(Cw))'/FRS.rhoI.*(EDG.UE(2:length(Cw))'/FRS.UI).^2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI);
Y10 = BLC.Cf2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI);
Y11 = solv(2:end)'.*Cw(2:end)'./sqrt(EDG.Re_x(2:length(Cw))').*EDG.rhoE(2:length(Cw))'/FRS.rhoI.*(EDG.UE(2:length(Cw))'/FRS.UI).^2*sqrt(FRS.UI*1/FRS.muI*FRS.rhoI); % check with Y10 for correctness, should be the same
% BLC.Cf2 = sol.v(1)*flp.C(1)/sqrt(EDG.Re_x(NS))*EDG.rhoE(NS)/FRS.rhoI*(EDG.UE(NS)/FRS.UI)^2;

%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
Color = [0 0 0];
ColorZoom = [0.65 0.65 0.65];
LineStyle = {'--','none','none','-.','-','--','--','none',':'};
Marker = {'none','*','+','none','none','none','none','x','none'};
LineWidth = 0.5;
FontSize = 11;
FontSizeZoom = 8;
LegendPos = 'southwest';%[0.20 0.55 0.20 0.20]; % normalized position
LegendName = {'Cebeci','Howarth','Clutter','G\"{o}rtler','CSM A','Smith A','Smith B','Hartree','CSM B'}; % Clutter (1963) is first author of this report!

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

% Zoom box and lines
plot([0.94 0.96 0.96 0.94 0.94],[0 0 0.04 0.04 0],'Color',ColorZoom,'HandleVisibility','off')
plot([0.94 0.45],[0.04 0.50],'Color',ColorZoom,'HandleVisibility','off') % left line
plot([0.96 0.96],[0.04 0.50],'Color',ColorZoom,'HandleVisibility','off') % right line

plot(X1,Y1,'Color',Color,'LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName',LegendName{1},'HandleVisibility','on')
plot(X1,Y2,'Color',Color,'LineStyle',LineStyle{2},'Marker',Marker{2},'DisplayName',LegendName{2},'HandleVisibility','on')
plot(X1,Y3,'Color',Color,'LineStyle',LineStyle{3},'Marker',Marker{3},'DisplayName',LegendName{3},'HandleVisibility','on')
plot(X1,Y4,'Color',Color,'LineStyle',LineStyle{4},'Marker',Marker{4},'DisplayName',LegendName{4},'HandleVisibility','on')
plot(X2,Y5,'Color',Color,'LineStyle',LineStyle{5},'Marker',Marker{5},'DisplayName',LegendName{5},'HandleVisibility','on')
plot(X3(1:end-1),Y6,'Color',Color,'LineStyle',LineStyle{6},'Marker',Marker{6},'DisplayName',LegendName{6},'HandleVisibility','on')
plot(X3(1:end-1),Y7,'Color',Color,'LineStyle',LineStyle{7},'Marker',Marker{7},'DisplayName',LegendName{7},'HandleVisibility','on')
plot(X3(1:end-1),Y8,'Color',Color,'LineStyle',LineStyle{8},'Marker',Marker{8},'DisplayName',LegendName{8},'HandleVisibility','on')
plot(X4(2:end),Y9,'Color',Color,'LineStyle',LineStyle{9},'Marker',Marker{9},'DisplayName',LegendName{9},'HandleVisibility','on')

set(legend, 'Interpreter', 'latex', 'FontSize', FontSize, 'LineWidth', LineWidth, 'Units', 'normalized', 'Location', LegendPos); % 'Location', 'best');
legend('show')

% keyboard

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

axz.XLim = [0.94 0.96];
axz.XLimMode = 'manual';
axz.XTickLabelMode = 'manual';
axz.XTickMode = 'manual';
axz.XTick = [0.94 0.95 0.96];
axz.XTickLabel = {'0.94'; '0.95'; '0.96'};

axz.YLim = [0 0.04];
axz.YLimMode = 'manual';
axz.YTickLabelMode = 'manual';
axz.YTickMode = 'manual';
axz.YTick = [0 0.01 0.02 0.03 0.04 ];
axz.YTickLabel = {'0'; '0.01'; '0.02'; '0.03'; '0.04'};

set(gca,'FontSize',FontSizeZoom)
box on
hold on

IndexOfInterest1 = (X1 < 0.97) & (X1 > 0.87);
IndexOfInterest2 = (X2 < 0.97) & (X2 > 0.87);
IndexOfInterest3 = (X3 < 0.97) & (X3 > 0.87);
IndexOfInterest4 = (X4 < 0.97) & (X4 > 0.87);

plot(X1(IndexOfInterest1),Y1(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{1},'Marker',Marker{1},'DisplayName',LegendName{1},'HandleVisibility','off')
plot(X1(IndexOfInterest1),Y2(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{2},'Marker',Marker{2},'DisplayName',LegendName{2},'HandleVisibility','off')
plot(X1(IndexOfInterest1),Y3(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{3},'Marker',Marker{3},'DisplayName',LegendName{3},'HandleVisibility','off')
plot(X1(IndexOfInterest1),Y4(IndexOfInterest1),'Color',Color,'LineStyle',LineStyle{4},'Marker',Marker{4},'DisplayName',LegendName{4},'HandleVisibility','off')
plot(X2(IndexOfInterest2),Y5(IndexOfInterest2),'Color',Color,'LineStyle',LineStyle{5},'Marker',Marker{5},'DisplayName',LegendName{5},'HandleVisibility','off')
plot(X3(IndexOfInterest3(1:end-1)),Y6(IndexOfInterest3(1:end-1)),'Color',Color,'LineStyle',LineStyle{6},'Marker',Marker{6},'DisplayName',LegendName{6},'HandleVisibility','off')
plot(X3(IndexOfInterest3(1:end-1)),Y7(IndexOfInterest3(1:end-1)),'Color',Color,'LineStyle',LineStyle{7},'Marker',Marker{7},'DisplayName',LegendName{7},'HandleVisibility','off')
plot(X3(IndexOfInterest3(1:end-1)),Y8(IndexOfInterest3(1:end-1)),'Color',Color,'LineStyle',LineStyle{8},'Marker',Marker{8},'DisplayName',LegendName{8},'HandleVisibility','off')
plot(X4(IndexOfInterest4),Y9(IndexOfInterest4(2:end)),'Color',Color,'LineStyle',LineStyle{9},'Marker',Marker{9},'DisplayName',LegendName{9},'HandleVisibility','off')
% % axis tight % use this command to find a nice coordinate range

% keyboard

% Store plot
NewFolder = cd; % currently in 'boundary-layer-code'-folder
cd('Thesis_Report_Figures') % (25-11-2022)
saveas(h,'Thesis_Report_Figures/Howarth_both','epsc')

OldFolder = cd(NewFolder);

% For reference:
% https://nl.mathworks.com/matlabcentral/answers/33779-zooming-a-portion-of-figure-in-a-figure
% indexOfInterest = (t < 11*pi/8) & (t > 9*pi/8); % range of t near perturbation
% plot(t(indexOfInterest),signal(indexOfInterest)) % plot on new axes
% axis tight
%%% Thesis figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%