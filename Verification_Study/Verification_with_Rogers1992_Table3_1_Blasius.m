%% Verification with Rogers (1992) table C-3.1: Blasius profile
% Similar solutions
% 24-09-2019 D.D.D.
% Last changes on 25-10-2019
% last run on 04-05-2020
% Simulation takes about three seconds (including graphs)
% Adapted slightly for upload and checked: 23-01-2022

%% Note: only for single station calculations
% test if data was correctly transformed and applied in string format
% \num{} command centers numbers and thus points are not aligned!
% how to place table in middle of page (centering sidewaystable?):
% https://tex.stackexchange.com/questions/102363/centring-a-sidewaystable-both-horizontally-and-vertically
% answer: add '\centering' after sidewaystable

clear all
close all
clc

% %% for testing (still to be done)
% ROG = [...
% 0    0.4696005     0            0
% 0.2  0.46930657    9.3905401e-2 9.391422e-3
% 0.4  0.4672547     0.18760534   0.03754924
% 0.6  0.46173493    0.28057576   8.4385663e-2
% 0.8  0.45119049    0.37196365   0.14967468
% 1.   0.43437958    0.46063307   0.23299035
% 1.2  0.41056575    0.54524709   0.33365774
% 1.4  0.37969252    0.62438697   0.4507241
% 1.6  0.34248737    0.69670022   0.58295692
% 1.8  0.3004455     0.76105808   0.728873
% 2.   0.25566929    0.8166954    0.88679774
% 2.2  0.21057998    0.86330501   1.0549482
% 2.4  0.16756036    0.90106625   1.2315289
% 2.6  0.12861283    0.93060206   1.4148256
% 2.8  9.5113386e-2  0.95287624   1.6032851
% 3.   6.7710286e-2  0.96905538   1.7955696
% 3.2  4.6370361e-2  0.98036575   1.9905828
% 3.4  3.0535214e-2  0.98797122   2.1874693
% 3.6  1.9328694e-2  0.99288865   2.3855926
% 3.8  1.1758678e-2  0.99594502   2.5845011
% 4.   6.8740853e-3  0.99777083   2.7838889
% 4.2  3.861352e-3   0.99881904   2.9835579
% 4.4  2.0840747e-3  0.99939734   3.1833855
% 4.6  1.0807525e-3  0.99970394   3.3832989
% 4.8  5.3848399e-4  0.99986013   3.5832571
% 5.   2.5778052e-4  0.99993659   3.7832377
% 5.2  1.1856508e-4  0.99997256   3.9832291
% 5.4  5.2395285e-5  0.99998882   4.1832254
% 5.6  2.2246211e-5  0.99999588   4.383224
% 5.8  9.0750329e-6  0.99999883   4.5832235
% 6.   3.556875e-6   1.           4.7832234...
% ];

%% INPUT
load('./Data/VerData_Rogers1992_Table3_1_Blasius')     % format string
run('./INPUT/INPUT_Verification_with_Rogers1992')   % case specific INPUT (options and settings)

% Changes in input for Blasius profile:
OPT.BCEE = 0;                                   % adiabatic wall
OPT.COMP = 1;                                   % (in)compressible

% Pressure gradient (parameter) input
PGm2 = 0;                                       % Pressure gradient parameter

% Set size of vectors for initilization (vertical/eta-direction)
SET.NPT = 6001;                                 % preallocation
SET.ITMAX = 6;                                  % maximum number of iterations
SET.ITMAX0 = 20;    %%% added (23-01-2022)
OPT.RLAM = 0;       %%% added (23-01-2022)

% Grid Parameters - adapted for laminar flows to constant VGP
GRD.etaE = 6.0*sqrt((2*FLD.C)/(PGm2 + 1));      % non-dimensional BL edge location
GRD.VGP = 1.0;                                  % variable grid parameter (spacing ratio)
GRD.Deta = 0.001*sqrt((2*FLD.C)/(PGm2 + 1));    % first grid spacing: h1

%% Initialization
Size = length(VerData_Rogers1992_Table3_1_Blasius);
% Beta_dact = zeros(Size,1);
% PGm2 = zeros(Size,1);
% gw = zeros(Size,1);
% CSM_DataTableC25 = NaN*ones(Size,2);
% CSM_DataTableC25 = NaN*ones(Size,5);
% sep = NaN*ones(Size,1);
% kk = NaN*ones(Size,1);
NS = 1;
NST = length(INP.x);
TCC = []; % only laminar flow
MON.SEP = 0; %% added on 03-09-2019
MON.tr = 0;
SOL = [];
BLC.tau = NaN*ones(SET.NPT,NST);
BLC.Sv = NaN*ones(SET.NPT,NST);

%% Precalculations (boundary conditions (edge and wall) and free stream conditions)
X = INP.x;

% Free stream conditions
FRS.PtI = INP.PtI;                                  % [Pa]
FRS.TtI = INP.TtI;                                  % [K]
FRS.MaI = INP.MaI;                                  % [-]
FRS.gammaI = FLD.gamma;                             % [-]
FRS.CpI = FLD.Cp;                                   % [J/kg/K]
FRS.HtI = FRS.CpI*INP.TtI;                          % [J/kg/K]
FRS.TsI = INP.TtI/(1 + (FRS.gammaI - 1)/2*INP.MaI^2); % [K]
FRS.HsI = FRS.CpI*FRS.TsI;                          % [J/kg]
FRS.aI = sqrt(FRS.gammaI*FLD.Rsg*FRS.TsI);          % [m/s]
FRS.UI = INP.MaI*FRS.aI;                            % [m/s]
FRS.PsI = INP.PtI*(FRS.TsI/INP.TtI)^(FRS.gammaI/(FRS.gammaI - 1)); % [Pa]
FRS.rhoI = FRS.PsI/FLD.Rsg/FRS.TsI;                 % [kg/m3], static density
FRS.muI = 1.45e-6*(FRS.TsI^1.5)/(FRS.TsI + FLD.SLV);% [Pas], static viscosity
FRS.kI = FRS.muI*FRS.CpI/FLD.Pr;                    % [W/m/K], thermal conductivity

% Edge conditions
EDG.MaE = INP.MaE;                        % [-]
EDG.CpE = FLD.Cp;                            % [J/kg/K]
EDG.gammaE = FLD.gamma;                      % [-]
EDG.HtE = FRS.HtI;                           % [J/kg]
EDG.TtE = INP.TtI;                           % [K]
EDG.PtE = FRS.PtI;                           % [Pa]
EDG.TsE = EDG.TtE/(1 + (EDG.gammaE - 1)/2*INP.MaE^2);% [K]
%             EDG.PsE(i) = FRS.PsI*(EDG.TsE(i)/FRS.TsI)^(EDG.gammaE(i)/(EDG.gammaE(i) - 1));% [Pa], gives WRONG RESULTS!!!
EDG.PsE = INP.PtI*(EDG.TsE/INP.TtI)^(EDG.gammaE/(EDG.gammaE - 1));% [Pa], works also for SP-flow?
EDG.aE = sqrt(EDG.gammaE*FLD.Rsg*EDG.TsE);% [m/s]
EDG.UE = EDG.MaE*EDG.aE;               % [m/s]
EDG.HsE = EDG.CpE*EDG.TsE;             % [J/kg]
EDG.rhoE = EDG.PsE/FLD.Rsg/EDG.TsE;    % [kg/m3], static density
EDG.muE = 1.45e-6*(EDG.TsE^1.5)/(EDG.TsE + FLD.SLV); % [Pas], static viscosity
EDG.kE = EDG.muE*EDG.CpE/FLD.Pr;       % [W/m/K], thermal conductivity
EDG.Re_x = EDG.rhoE*EDG.UE*X/EDG.muE;

% Boundary conditions
if OPT.BCEE > 0
    INP.BCW = gw(jj)*ones(1,NST);   % [-], scaled h ratio or it's derivative
    INP.AWD = []; % only g_aw or T_aw ((static) ratio? -> see OUTPUT) should be given here, calculated after adiabatic simulation
else
    INP.BCW = zeros(1,NST);   % [-], zero wall heat flux: derivative at wall of enthalpy ratio equal to zero
end

% Pressure gradient parameters
HVR.P2 = PGm2;
HVR.P1 = 0.5*(HVR.P2(1) + 1);
HVR.P1P = HVR.P1;
HVR.P2P = HVR.P2;
HVR.CEL = 0; % defined above
if  OPT.BCEE == 0 % adiabatic: derivative of enthalpy ratio = 0 [-]
    HVR.alpha0 = 0;
    HVR.alpha1 = 1;
    HVR.WW = INP.BCW; % zeros
elseif OPT.BCEE == 1 % enthalpy ratio [-]
    HVR.alpha0 = 1;
    HVR.alpha1 = 0;
    HVR.WW = INP.BCW;
end

%% Grid
[NP,GRD] = GRID(GRD,SET);

%% Initial profiles
[sol,solprev] = IVPL(NS,GRD,HVR,FLD,SET);

%% Run solver (note: sol and solprev remain the same outside CSM-file and are re-used)
[BLC,FLP,SOL,TCC,MON] = CSM(X,NP,NS,NST,GRD,HVR,EDG,FRS,FLD,TCC,sol,solprev,SOL,BLC,INP,OPT,SET,MON);

%% Store data - Transformed
CSM_TransfDataTable3_1(:,1) = GRD.eta(1:200:200*Size)/sqrt((2*FLD.C)/(PGm2 + 1)); % \eta
CSM_TransfDataTable3_1(:,2) = SOL{end}.v(1:200:200*Size)*sqrt(2*FLD.C/(PGm2 + 1)); % f''
CSM_TransfDataTable3_1(:,3) = SOL{end}.u(1:200:200*Size);                          % f'
CSM_TransfDataTable3_1(:,4) = SOL{end}.f(1:200:200*Size)/sqrt(2*FLD.C/(PGm2 + 1)); % f
CSM_TransfDataTable3_1(:,5) = SOL{end}.p(1:200:200*Size)*sqrt(2*FLD.C/(PGm2 + 1)); % p
CSM_TransfDataTable3_1(:,6) = SOL{end}.g(1:200:200*Size);                          % g

%% Differences

% % absolute
% for j = 1:Size
%     eta(j,1) = GRD.eta(1+(j-1)*200,1);
%     C(j,1) = SOL{end}.v(1+(j-1)*200,2) - ROG(j,2);
%     C(j,2) = SOL{end}.u(1+(j-1)*200,3) - ROG(j,3);
%     C(j,3) = SOL{end}.f(1+(j-1)*200,4) - ROG(j,4);
% end
% 
% % relative
% for j = 1:Size
% %     eta(j,1) = CSM_TransfDataTable3_1(1+(j-1)*200,1);
%     D(j,1) = (SOL{end}.v(1+(j-1)*200,2) - ROG(j,2))/ROG(j,2)*100;
%     D(j,2) = (SOL{end}.u(1+(j-1)*200,3) - ROG(j,3))/ROG(j,3)*100;
%     D(j,3) = (SOL{end}.f(1+(j-1)*200,4) - ROG(j,4))/ROG(j,4)*100;
% end

% % absolute
% for j = 1:Size
%     eta(j,1) = CSM_TransfDataTable3_1(1+(j-1)*200,1);
%     A(j,1) = CSM_TransfDataTable3_1(1+(j-1)*200,2) - str2num(VerData_Rogers1992_Table3_1_Blasius(j,2));
%     A(j,2) = CSM_TransfDataTable3_1(1+(j-1)*200,3) - str2num(VerData_Rogers1992_Table3_1_Blasius(j,3));
%     A(j,3) = CSM_TransfDataTable3_1(1+(j-1)*200,4) - str2num(VerData_Rogers1992_Table3_1_Blasius(j,4));
% end

% % relative
% for j = 1:Size
%     eta(j,1) = CSM_TransfDataTable3_1(1+(j-1)*200,1);
%     B(j,1) = (CSM_TransfDataTable3_1(1+(j-1)*200,2) - str2num(VerData_Rogers1992_Table3_1_Blasius(j,2)))/str2num(VerData_Rogers1992_Table3_1_Blasius(j,2))*100;
%     B(j,2) = (CSM_TransfDataTable3_1(1+(j-1)*200,3) - str2num(VerData_Rogers1992_Table3_1_Blasius(j,3)))/str2num(VerData_Rogers1992_Table3_1_Blasius(j,3))*100;
%     B(j,3) = (CSM_TransfDataTable3_1(1+(j-1)*200,4) - str2num(VerData_Rogers1992_Table3_1_Blasius(j,4)))/str2num(VerData_Rogers1992_Table3_1_Blasius(j,4))*100;
% end

eta = CSM_TransfDataTable3_1(:,1);
% absolute
A = CSM_TransfDataTable3_1(:,2:4) - str2double(VerData_Rogers1992_Table3_1_Blasius(:,2:4));
% relative
B = (CSM_TransfDataTable3_1(:,2:4) - str2double(VerData_Rogers1992_Table3_1_Blasius(:,2:4)))./str2double(VerData_Rogers1992_Table3_1_Blasius(:,2:4))*100;
% A = CSM_TransfDataTable3_1(1+(j-1)*200,3) - str2num(VerData_Rogers1992_Table3_1_Blasius(j,3));
% A = CSM_TransfDataTable3_1(1+(j-1)*200,4) - str2num(VerData_Rogers1992_Table3_1_Blasius(j,4));



%% Print results
fprintf('Input Pressure gradient m2 = %7.9f.\n', PGm2)
fprintf('Input Pressure gradient Beta = %7.9f.\n', 2*PGm2/(PGm2 + 1)*(1+(FLD.gamma - 1)/2*INP.MaI^2))
if HVR.alpha0 > 0
    fprintf('Input gw = %7.9f.\n', INP.BCW)
else
    fprintf('Input g''w = %7.9f.\n', INP.BCW)
end
fprintf('CSM transformed fw'''' = %7.9f.\n', SOL{end}.v(1)*sqrt(2*FLD.C/(PGm2 + 1)))
fprintf('CSM transformed gw = %7.9f.\n', SOL{end}.g(1))

%% Plot results

% absolute differences
h1 = figure;
hold on
plot(eta,A(:,1))
plot(eta,A(:,2))
plot(eta,A(:,3))
title('absolute differences')
xlabel('\eta [-]')
ylabel('d [-]')
legend('dv','du','df')
title('absolute differences')

% relative differences
h2 = figure;
hold on
plot(eta,B(:,1))
plot(eta,B(:,2))
plot(eta,B(:,3))
title('relative differences: (CSM - Rogers)/Rogers*100%')
xlabel('\eta [-]')
ylabel('d/d [-]')
legend('dv','du','df')

% %% TEST
% % absolute differences
% figure
% hold on
% plot(eta,C(:,1))
% plot(eta,C(:,2))
% plot(eta,C(:,3))
% title('absolute differences')
% xlabel('\eta [-]')
% ylabel('d [-]')
% legend('dv','du','df')
% title('absolute differences')
% 
% % relative differences
% figure
% hold on
% plot(eta,D(:,1))
% plot(eta,D(:,2))
% plot(eta,D(:,3))
% title('relative differences: (CSM - Rogers)/Rogers*100%')
% xlabel('\eta [-]')
% ylabel('d/d [-]')
% legend('dv','du','df')
% %% TEST

% %% Format numbers in scientific notation
% % formatting scientific notation before printing numbers
% pattern = 'e';
% Table3_1 = VerData_Rogers1992_Table3_1_Blasius;
% % check for scientific notation (needs to be formatted before printing to table)
% for i = 1:Size
%     for j = 1:size(Table3_1,2)
%         if contains(Table3_1(i,j),pattern)
% %             Table3_1(i,j) = fprintf('\\num{%s}',Table3_1(i,j));
% %             Table3_1(i,j) = fprintf('\\num{%s}',Table3_1(i,j));
%             Table3_1(i,j) = strcat('\\num{',Table3_1(i,j),'}');
%         end
%     end
% end

%% Format numbers in scientific notation
% formatting scientific notation before printing numbers
pattern = 'e';
% check for scientific notation (needs to be formatted before printing to table)
for i = 1:Size
    for j = 1:size(VerData_Rogers1992_Table3_1_Blasius,2)
        if contains(VerData_Rogers1992_Table3_1_Blasius(i,j),pattern)
%             Table3_1(i,j) = fprintf('\\num{%s}',Table3_1(i,j));
%             Table3_1(i,j) = fprintf('\\num{%s}',Table3_1(i,j));
%             VerData_Rogers1992_Table3_1_Blasius(i,j) = strcat('\num{',VerData_Rogers1992_Table3_1_Blasius(i,j),'}');
            VerData_Rogers1992_Table3_1_Blasius(i,j) = strcat('\multicolumn{1}{l}{\num{',VerData_Rogers1992_Table3_1_Blasius(i,j),'}}');
        end
    end
end

% keyboard

%% Print Table 3-1 - Blasius' solution

NewFolder = cd; % currently in 'boundary-layer-code'-folder
cd('Thesis_Report_Tables') %%% new, added (23-01-2021)
fid = fopen('Verification_Rogers1992_Table3_1_Blasius.tex','w');
fprintf(fid, '%%%% Verification with Table 3-1 from Rogers\\cite{rogers1992laminar}: rho=constant, Pr=1, C=1 (constant)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
% fprintf(fid, '    \\centering\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Comparison with tabulated data of the Blasius'' solution taken from Rogers\\cite{rogers1992laminar} table 3-1. With uniform (vertical) grid spacing (in compressible Falkner-Skan transformed y-coordinate) of $\\mathrm{d} \\eta = \\sqrt{\\frac{2C}{m_2 + 1}} %1.3f$ and heigth of $\\eta_{\\mathrm{e}} = \\sqrt{\\frac{2C}{m_2 + 1}} %1.1f$.}\n',GRD.Deta(1)/sqrt(2*FLD.C/(PGm2 + 1)),GRD.etaE/sqrt(2*FLD.C/(PGm2 + 1)));
fprintf(fid, '    \\label{tab:31B}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=1.2] S[table-format=1.10] S[table-format=1.10] S[table-format=1.10] S[table-format=1.10] S[table-format=1.10] S[table-format=1.10]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\eta$}                                &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$f''''$}                                 &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$f''$}                                  &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$f$}                                   \\\\\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          \\\\\n');
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for i = 1:Size
    %         fprintf(fid, '   %s   &   %s   &   %1.6f   &   %s   &   %1.6f   &   %s   &   %1.6f   \\\\\n',VerData_Rogers1992_Table3_1_Blasius(i,1),VerData_Rogers1992_Table3_1_Blasius(i,2),CSM_TransfDataTable3_1(1+(i-1)*200,2),VerData_Rogers1992_Table3_1_Blasius(i,3),CSM_TransfDataTable3_1(1+(i-1)*200,3),VerData_Rogers1992_Table3_1_Blasius(i,4),CSM_TransfDataTable3_1(1+(i-1)*200,4));
    % working original %	fprintf(fid, '        %1.1f   &   %s   &   %1.6f   &   %s   &   %1.6f   &   %s   &   %1.6f   \\\\\n',VerData_Rogers1992_Table3_1_Blasius(i,1),VerData_Rogers1992_Table3_1_Blasius(i,2),CSM_TransfDataTable3_1(i,2),VerData_Rogers1992_Table3_1_Blasius(i,3),CSM_TransfDataTable3_1(i,3),VerData_Rogers1992_Table3_1_Blasius(i,4),CSM_TransfDataTable3_1(i,4));
%     for j = 1:size(VerData_Rogers1992_Table3_1_Blasius,2)
%         if contains(VerData_Rogers1992_Table3_1_Blasius(i,j),pattern)
%             %             Table3_1(i,j) = fprintf('\\num{%s}',Table3_1(i,j));
%             Table3_1(i,j) = fprintf('\\num{%s}',Table3_1(i,j));
%         end
%     end
    fprintf(fid, '        %1.1f   &   %s   &   %1.8f   &   %s   &   %1.8f   &   %s   &   %1.8f   \\\\\n',VerData_Rogers1992_Table3_1_Blasius(i,1),VerData_Rogers1992_Table3_1_Blasius(i,2),CSM_TransfDataTable3_1(i,2),VerData_Rogers1992_Table3_1_Blasius(i,3),CSM_TransfDataTable3_1(i,3),VerData_Rogers1992_Table3_1_Blasius(i,4),CSM_TransfDataTable3_1(i,4));
end
% fprintf(fid, '        %% Data\n');
fprintf(fid, '        \\bottomrule\n\n');
fprintf(fid, '    \\end{tabular}\n');
fprintf(fid, '\\end{threeparttable}\n');
% fprintf(fid, '%%\\end{sidewaystable}\n');
% fprintf(fid, '%%\\end{table\n');
fclose(fid);
OldFolder = cd(NewFolder);
