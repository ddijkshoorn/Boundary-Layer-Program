%% Verification with Rogers (1992) table C-1: Falkner-Skan wedge flows
% Verification of similar solutions with tabulated data from Rogers (1992) table C1
% adiabatic nonzero pressure gradient flow with constant density and
% constant fluid properties
% last changes on 25-10-2019
% last run on 04-05-2020
% Simulation takes about 4 seconds
% Adapted slightly for upload and checked: 23-01-2022

%% Note: only for single station calculations
% 1) if grid is adapted, step-size in caption also needs to be adapted:
% Be aware of the step-size and etaE given below (before calculation loop)
% and in the captions. In Illingworth-Levy (IL) coordinates they don't
% change, in compressible Falkner-Skan (FS) coordinates they do with every
% change in Pressure Gradient (PG).
% 2) Rogers changes eta_max (etaE) for Beta = -0.05 from 6.0 to 9.0. This
% code then predicts separation (etaE value in compr FS coordinates seems
% too high for this program)

clear all
close all
clc

%% INPUT
load('./Data/VerData_Rogers1992_TableC1_FalknerSkan') %%% added on 24-10-2019
run('./INPUT/INPUT_Verification_with_Rogers1992')         % case specific INPUT (options and settings; INPUT-file is standardized for case C-25)

% Change of input (compared with Table C-25 and C-26, where INPUT-file is standardized for case C-25):
OPT.COMP = 0;
OPT.BCEE = 0;

%% Initialization
Size = length(VerData_Rogers1992_TableC1_FalknerSkan); % for entry 21 separation occurs?
Beta = zeros(Size,1);
PGm2 = zeros(Size,1);
Deta = zeros(Size,1);
etaE = zeros(Size,1);
etaR = zeros(Size,1);
CSM_DataTableC1 = NaN*ones(Size,1); %%% f''(0), I1, I2
sep = NaN*ones(Size,1);
kk = NaN*ones(Size,1);
NS = 1;
NST = length(INP.x);
TCC = []; % only laminar flow
MON.SEP = 0;
MON.tr = 0;
SET.ITMAX0 = 20;    %%% added (23-01-2022)
OPT.RLAM = 0;       %%% added (23-01-2022)

%% Precalculations (boundary layer edge and free stream conditions)
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
EDG.MaE = INP.MaE;                                  % [-]
EDG.CpE = FLD.Cp;                                   % [J/kg/K]
EDG.gammaE = FLD.gamma;                             % [-]
EDG.HtE = FRS.HtI;                                  % [J/kg]
EDG.TtE = INP.TtI;                                  % [K]
EDG.PtE = FRS.PtI;                                  % [Pa]
EDG.TsE = EDG.TtE/(1 + (EDG.gammaE - 1)/2*INP.MaE^2);% [K]
%             EDG.PsE(i) = FRS.PsI*(EDG.TsE(i)/FRS.TsI)^(EDG.gammaE(i)/(EDG.gammaE(i) - 1));% [Pa], gives WRONG RESULTS!!!
EDG.PsE = INP.PtI*(EDG.TsE/INP.TtI)^(EDG.gammaE/(EDG.gammaE - 1));% [Pa], works also for SP-flow?
EDG.aE = sqrt(EDG.gammaE*FLD.Rsg*EDG.TsE);          % [m/s]
EDG.UE = EDG.MaE*EDG.aE;                            % [m/s]
EDG.HsE = EDG.CpE*EDG.TsE;                          % [J/kg]
EDG.rhoE = EDG.PsE/FLD.Rsg/EDG.TsE;                 % [kg/m3], static density
EDG.muE = 1.45e-6*(EDG.TsE^1.5)/(EDG.TsE + FLD.SLV);% [Pas], static viscosity
EDG.kE = EDG.muE*EDG.CpE/FLD.Pr;                    % [W/m/K], thermal conductivity
EDG.Re_x = EDG.rhoE*EDG.UE*X/EDG.muE;

%% see inside loop
% Grid
% GRD.etaE1 = 6.0;
GRD.h1 = 0.001; % #2
% [NP,GRD] = GRID(GRD,SET);	%%% OLD method (commenting all '#2')
% grid spacing/step-size used by Rogers: h=0.01 in Illingworth-Levy coordinates (Rogers (1992) Appendix D: 'D-13 Main Program for the Falkner-Skan Equation' p.369)

%% Initial profiles
% calculated inside loop for every new Beta again, since dependent on Beta

%% Main Calculation
ii = 0; % table line (row) index (skip separated flows in table!)

for jj = 1:Size
    
    Beta(jj) = str2double(VerData_Rogers1992_TableC1_FalknerSkan(jj,1));
    PGm2(jj) = Beta(jj)*(1+(FLD.gamma - 1)/2*INP.MaI^2)^-1/(2 - Beta(jj)*(1+(FLD.gamma - 1)/2*INP.MaI^2)^-1);
    
    if jj > 1 && (Beta(jj) > Beta(jj - 1) && sign(Beta(jj)) == sign(Beta(jj - 1)))
        % skip line
        fprintf('Input Beta = %7.9f is skipped since it is part of the second solution for decelarating flow which cannot be calculated with this program.\n', Beta(jj))
    else % run calculation
        ii = ii + 1; % next line in table
        kk(ii) = jj; % store location entry
        
        % Extraordinary cases (Beta_dact(jj) = 2.0)
            if Beta(jj) == 2.0
                PGm2(jj,1) = 1000;
            end
% % %             SET.NPT = 100000;
% % %             GRD.Deta(1) = 0.0001;
% % %         else
% % %             % normal:
% % %             SET.NPT = 1000;
% % %             GRD.Deta = 0.01;
% % %         end
        
        % Initialize (because of new grid parameters only, fix this in another way?)
        SOL = [];
        BLC.tau = NaN*ones(SET.NPT,NST);
        BLC.Sv = NaN*ones(SET.NPT,NST);
        
% #2 %  % Grid (useful: reducing the amount of grid points for large PG flow)
        if jj == 1
            GRD.etaE1 = str2double(VerData_Rogers1992_TableC1_FalknerSkan(jj,5));      % non-dimensional BL edge location in Illingworth-Levy coordinates
%             GRD.etaE1 = 8; %%% (30-10-2019) gives exactly same results
%             (amount of digits) for Beta=-0.16, but separation predicted
%             for other three cases.
        % eta_max (eta_E) is changed (line 14):
%         elseif abs(str2double(VerData_Rogers1992_TableC1_FalknerSkan(jj,5)) - str2double(VerData_Rogers1992_TableC1_FalknerSkan(jj - 1,5))) > 0
%             GRD.etaE1 = str2double(VerData_Rogers1992_TableC1_FalknerSkan(jj,5));      % non-dimensional BL edge location in Illingworth-Levy coordinates
% % % %             GRD.etaE1 = 8;
% % % %         elseif jj == 16
% % % %             GRD.etaE1 = 8;
% % % %         elseif jj == 17
% % % %             GRD.etaE1 = 9;
% % % %         elseif jj == 18
% % % %             GRD.etaE1 = 9;
% % % %             GRD.etaE1 = 11;
% % % %         elseif jj == 19
% % % %             GRD.etaE1 = 12;
% % % %         elseif jj == 20
% % % %             GRD.etaE1 = 9;
        end
        GRD.Deta = GRD.h1*sqrt((2*FLD.C)/(PGm2(jj) + 1));    % first grid spacing: h1 in compressible Falkner-Skan coordinates
        GRD.etaE = GRD.etaE1*sqrt((2*FLD.C)/(PGm2(jj) + 1));
        [NP,GRD] = GRID(GRD,SET);
        Deta(jj) = GRD.Deta(1);
        etaE(jj) = GRD.etaE;
%         etaR(jj) = GRD.Deta(1)/GRD.etaE; %%% added (28-10-2019), same ratio of course, since linear transformation
% #2 %  %
        if OPT.BCEE > 0
            INP.BCW = gw(jj)*ones(1,length(INP.x));   % [-], scaled h ratio or it's derivative
            INP.AWD = []; % only g_aw or T_aw ((static) ratio? -> see OUTPUT) should be given here, calculated after adiabatic simulation
        else
            INP.BCW = zeros(1,length(INP.x));   % [-], zero wall heat flux: derivative at wall of enthalpy ratio equal to zero
        end
        
        HVR.P2 = PGm2(jj);
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
        
        % initial profiles %%% (28-10-2019) note that solprev is initialized here and is empty []
%         if ~exist('sol','var')
        [sol,solprev] = IVPL(NS,GRD,HVR,FLD,SET); % needs to be here because of use of HVR
%         end
        
        % Run code (note: sol and solprev remain the same outside CSM-file and are re-used)
        [BLC,FLP,SOL,TCC,MON] = CSM(X,NP,NS,NST,GRD,HVR,EDG,FRS,FLD,TCC,sol,solprev,SOL,BLC,INP,OPT,SET,MON);
        
        % Store results (in table format/size, which is smaller than Table from Rogers)
        if MON.SEP == 1
            CSM_DataTableC1(ii,1) = NaN;
            CSM_DataTableC1(ii,2) = NaN;
            CSM_DataTableC1(ii,3) = NaN;
            CSM_DataTableC1(ii,4) = NaN;
            CSM_DataTableC1(ii,5) = NaN;
            sep(ii,1) = 1;
        else
            CSM_DataTableC1(ii,1) = SOL{end}.v(1); % storing shear parameter; or in transformed form: SOL{end}.v(1)*sqrt(2/(m+1));
            %%% Added for Extended table (J1, J2, J3) #1
            % Caclulate J1, J2, J3, J4 and J5? #1 NB J5 defined by myself (25-10-2019)
            TERMP1 = 0;
            TERMP2 = 0;     % momentum thickness
            TERMP3 = 0;
            TERMP4 = 0;
            TERMP5 = 0;     % added (25-10-2019)
            SUM1 = 0;
            SUM2 = 0;       % momentum thickness
            SUM3 = 0;
            SUM4 = 0;
            SUM5 = 0;       % added (25-10-2019)
            for mm = 2:NP
                TERM1 = SOL{end}.g(mm) - SOL{end}.u(mm)^2;
                TERM2 = SOL{end}.u(mm)*(1 - SOL{end}.u(mm));    % momentum thickness
                TERM3 = 1 - SOL{end}.g(mm);
                TERM4 = 1 - SOL{end}.u(mm)^2;
                TERM5 = (1 - SOL{end}.u(mm));                   % added (25-10-2019)
                SUM1 = SUM1 + GRD.A(mm)*(TERM1 + TERMP1);
                SUM2 = SUM2 + GRD.A(mm)*(TERM2 + TERMP2);       % momentum thickness
                SUM3 = SUM3 + GRD.A(mm)*(TERM3 + TERMP3);
                SUM4 = SUM4 + GRD.A(mm)*(TERM4 + TERMP4);
                SUM5 = SUM5 + GRD.A(mm)*(TERM5 + TERMP5);       % added (25-10-2019)
            end
            J1 = SUM1; % NOTE! NOT: in Illingworth-Levy transformed plane (/sqrt(2*FLD.C/(PGm2(jj)+1))); BUT: J1_tr = J1*2/sqrt(2*FLD.C/(PGm2(jj)+1))
            J2 = SUM2; % NOTE! NOT: in Illingworth-Levy transformed plane (/sqrt(2*FLD.C/(PGm2(jj)+1))); BUT: J2_tr = J2*2/sqrt(2*FLD.C/(PGm2(jj)+1))
            J3 = SUM3; % NOTE! NOT: in Illingworth-Levy transformed plane (/sqrt(2*FLD.C/(PGm2(jj)+1))); BUT: J3_tr = J4*2/sqrt(2*FLD.C/(PGm2(jj)+1))
            J4 = SUM4; % NOTE! NOT: in Illingworth-Levy transformed plane (/sqrt(2*FLD.C/(PGm2(jj)+1))); BUT: J4_tr = J?*2/sqrt(2*FLD.C/(PGm2(jj)+1))
            J5 = SUM5; % added (25-10-2019)
            CSM_DataTableC1(ii,2) = J1;
            CSM_DataTableC1(ii,3) = J2;
            CSM_DataTableC1(ii,4) = J3;
            CSM_DataTableC1(ii,5) = J4;
            CSM_DataTableC1(ii,6) = J5;
            sep(ii,1) = 0;
        end
        
        % Print results
        if MON.SEP == 1
            if  OPT.BCEE == 0 % adiabatic: derivative of enthalpy ratio = 0 [-] (HVR.alpha0 = 0;)
                fprintf('Separation occurred. No results for input Beta(jj) = %2.6f \n', Beta(jj), INP.BCW)
            else
                fprintf('Separation occurred. No results for input Beta(jj) = %2.6f \n', Beta(jj), gw(jj))
            end
        else
            fprintf('Input Pressure gradient m2 = %7.9f.\n', PGm2(jj))
            fprintf('Input Pressure gradient Beta = %7.9f.\n', 2*PGm2(jj)/(PGm2(jj) + 1)*(1+(FLD.gamma - 1)/2*INP.MaI^2))

            if  OPT.BCEE == 0 % adiabatic: derivative of enthalpy ratio = 0 [-] (HVR.alpha0 = 0;)
                fprintf('Input g''w = %7.9f.\n', INP.BCW)
            else
                fprintf('Input gw = %7.9f.\n', gw(jj))
            end
            fprintf('CSM transformed fw'''' = %7.9f.\n', SOL{end}.v(1)*sqrt(2*FLD.C/(PGm2(jj) + 1)))
            if  OPT.BCEE == 0 % adiabatic: derivative of enthalpy ratio = 0 [-] (HVR.alpha0 = 0;)
                fprintf('CSM gw = %7.9f.\n', SOL{end}.g(1))
            else
                fprintf('CSM transformed gw'' = %7.9f.\n', SOL{end}.p(1)*sqrt(2*FLD.C/(PGm2(jj) + 1)))
            end
            % use result as initial profile for next calculation:
%             sol = SOL{end};
%             solprev = [];
        end
        
% % %         % Settings back to normal
% % %         if Beta_dact(jj) == 2 || gw(jj) == 0
% % %             SET.NPT = 1000;
% % %             GRD.Deta = 0.01;
% % %         end
        
        % clear output
        clear BLC FLP SOL MON
        
        % and redefine for next run
        MON.SEP = 0;
        MON.tr = 0;
        
    end
    
end

% keyboard

%% 2 tables without etaE
%% Table C-1 - Write data to Tex-file (in Table format)

NewFolder = cd; % currently in 'boundary-layer-code'-folder
cd('Thesis_Report_Tables') %%% new, added (23-01-2022)
fid = fopen('Verification_Rogers1992_TableC1_FalknerSkan.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Verification with Table C-1 from Rogers\\cite{rogers1992laminar}: Falkner-Skan wedge flows (adiabatic incompressible nonzero pressure gradient flows with constant fluid properties)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
% fprintf(fid, '    \\caption{Verification of the Falkner-Skan wedge flows with tabulated data obtained from Rogers\\cite{rogers1992laminar} table C-1. Similar (nonzero) pressure gradient flows in adiabatic constant density constant fluid property flows. The values obtained with the CS-method (CSM) are transformed from the compressible Falkner-Skan transformed y-coordinate with uniform (vertical) grid spacing of $\\mathrm{d} \\eta = \\sqrt{\\frac{m_2 + 1}{2C}} %1.4f$ and height of $\\eta_{\\mathrm{e}} = \\sqrt{\\frac{m_2 + 1}{2C}} %1.1f$ to the Illingworth-Levy coordinates ($\\mathrm{d} \\eta = %1.4f$ and $\\eta_{\\mathrm{e}} = %1.1f$).}\n',GRD.Deta(1),GRD.etaE,GRD.Deta(1),GRD.etaE);
fprintf(fid, '    \\caption{Comparison with tabulated data of the Falkner-Skan wedge flows (adiabatic incompressible similar flows with nonzero constant pressure gradients for calorically perfect ideal gas with constant fluid properties) taken from Rogers\\cite{rogers1992laminar} table C-1. The values obtained with the CS-method (CSM) are transformed from the compressible Falkner-Skan transformed y-coordinate with uniform (vertical) grid spacing of $\\mathrm{d} \\eta = \\sqrt{\\frac{2C}{m_2 + 1}} %1.4f$ and height of $\\eta_{\\mathrm{e}} = \\sqrt{\\frac{2C}{m_2 + 1}} %1.1f$ to the Illingworth-Levy coordinates ($\\mathrm{d} \\eta = %1.4f$ and $\\eta_{\\mathrm{e}} = %1.1f$). Note that separation occurred when the table entry shows ''sep''.}\n',GRD.h1,GRD.etaE1,GRD.h1,GRD.etaE1);
fprintf(fid, '    \\label{tab:FSF}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=-1.7] S[table-format=4.6] S[table-format=1.8] S[table-format=1.7]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\beta$}                               &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$m_2$}                                 &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$f''''(0)$}                              \\\\\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          \\\\\n');
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:ii % ii is counter inside loop above (counting valid table rows), ll is index of resulting CSM_data, mm is index of table data from Rogers 1992
    mm = kk(ll);
    if Beta(mm) == 2.0
        fprintf(fid, '         %s       &   %4.0f\\tnote{*}  &  %s   &  %1.7f \\\\\n',VerData_Rogers1992_TableC1_FalknerSkan(mm,1),PGm2(mm),VerData_Rogers1992_TableC1_FalknerSkan(mm,2),CSM_DataTableC1(mm,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)));
    else
        fprintf(fid, '         %s       &   %1.6f  &  %s   &  %1.7f \\\\\n',VerData_Rogers1992_TableC1_FalknerSkan(mm,1),PGm2(mm),VerData_Rogers1992_TableC1_FalknerSkan(mm,2),CSM_DataTableC1(mm,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)));
    end
end
% fprintf(fid, '        %% Data\n');
fprintf(fid, '        \\bottomrule\n\n');
fprintf(fid, '    \\end{tabular}\n');
fprintf(fid, '    \\begin{tablenotes}\n');
fprintf(fid, '        \\item[*] The solution does not converge for $m_2 = \\infty$, and therefore $m_2 = 1000$ is taken as approach.\n');
fprintf(fid, '    \\end{tablenotes}\n');
fprintf(fid, '\\end{threeparttable}\n');
% fprintf(fid, '%%\\end{sidewaystable}\n');
% fprintf(fid, '%%\\end{table\n');
fclose(fid);
% OldFolder = cd(NewFolder);

% keyboard

%% Table C-1 extended - Write data to Tex-file (in Table format)

fid = fopen('Verification_Rogers1992_TableC1_FalknerSkan_extended.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Verification with Table C-1 from Rogers\\cite{rogers1992laminar}: Falkner-Skan wedge flows (adiabatic incompressible nonzero pressure gradient flows with constant fluid properties)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
% fprintf(fid, '    \\centering\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
% fprintf(fid, '    \\caption{Verification of the Falkner-Skan wedge flows with tabulated data obtained from Rogers\\cite{rogers1992laminar} table C-1. Similar (nonzero) pressure gradient flows in adiabatic constant density constant fluid property flows. The values obtained with the CS-method (CSM) are transformed from the compressible Falkner-Skan transformed y-coordinate with uniform (vertical) grid spacing of $\\mathrm{d} \\eta = \\sqrt{\\frac{m_2 + 1}{2C}} %1.4f$ and height of $\\eta_{\\mathrm{e}} = \\sqrt{\\frac{m_2 + 1}{2C}} %1.1f$ to the Illingworth-Levy coordinates ($\\mathrm{d} \\eta = %1.4f$ and $\\eta_{\\mathrm{e}} = %1.1f$).}\n',GRD.Deta(1),GRD.etaE,GRD.Deta(1),GRD.etaE);
fprintf(fid, '    \\caption{Comparison with tabulated data of the Falkner-Skan wedge flows (adiabatic incompressible similar flows with nonzero constant pressure gradients for calorically perfect ideal gas with constant fluid properties) taken from Rogers\\cite{rogers1992laminar} table C-1. The values obtained with the CS-method (CSM) are transformed from the compressible Falkner-Skan transformed y-coordinate with uniform (vertical) grid spacing of $\\mathrm{d} \\eta = \\sqrt{\\frac{2C}{m_2 + 1}} %1.4f$ and height of $\\eta_{\\mathrm{e}} = \\sqrt{\\frac{2C}{m_2 + 1}} %1.1f$ to the Illingworth-Levy coordinates ($\\mathrm{d} \\eta = %1.4f$ and $\\eta_{\\mathrm{e}} = %1.1f$). Note that separation occurred when the table entry shows ''sep''.}\n',GRD.h1,GRD.etaE1,GRD.h1,GRD.etaE1);
fprintf(fid, '    \\label{tab:FSFE}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=-1.7] S[table-format=4.6] S[table-format=1.8] S[table-format=1.7] S[table-format=1.7] S[table-format=1.7] S[table-format=1.7] S[table-format=1.7]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\beta$}                               &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$m_2$}                                 &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$f''''(0)$}                              &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          \\\\\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$I_1$}                                 &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$2 J_{\\mathrm{x}}$}                    &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$I_2$}                                 &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$2 J_2$}                               \\\\\n');
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:ii % ii is counter inside loop above (counting valid table rows), ll is index of resulting CSM_data, mm is index of table data from Rogers 1992
    mm = kk(ll);
    if Beta(mm) == 2.0
        fprintf(fid, '         %s       &   %4.0f\\tnote{*}  &  %s   &  %1.6f  &  %s   &  %1.7f  &  %s   &  %1.7f \\\\\n',VerData_Rogers1992_TableC1_FalknerSkan(mm,1),PGm2(mm),VerData_Rogers1992_TableC1_FalknerSkan(mm,2),CSM_DataTableC1(mm,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC1_FalknerSkan(mm,3),CSM_DataTableC1(mm,6)*2/sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC1_FalknerSkan(mm,4),CSM_DataTableC1(mm,3)*2/sqrt(2*FLD.C/(PGm2(mm) + 1)));
    else
        fprintf(fid, '         %s       &   %1.6f  &  %s   &  %1.6f  &  %s   &  %1.7f  &  %s   &  %1.7f \\\\\n',VerData_Rogers1992_TableC1_FalknerSkan(mm,1),PGm2(mm),VerData_Rogers1992_TableC1_FalknerSkan(mm,2),CSM_DataTableC1(mm,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC1_FalknerSkan(mm,3),CSM_DataTableC1(mm,6)*2/sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC1_FalknerSkan(mm,4),CSM_DataTableC1(mm,3)*2/sqrt(2*FLD.C/(PGm2(mm) + 1)));
    end
end
% fprintf(fid, '        %% Data\n');
fprintf(fid, '        \\bottomrule\n\n');
fprintf(fid, '    \\end{tabular}\n');
fprintf(fid, '    \\begin{tablenotes}\n');
fprintf(fid, '        \\item[*] The solution does not converge for $m_2 = \\infty$, and therefore $m_2 = 1000$ is taken as approach.\n');
fprintf(fid, '    \\end{tablenotes}\n');
fprintf(fid, '\\end{threeparttable}\n');
% fprintf(fid, '%%\\end{sidewaystable}\n');
% fprintf(fid, '%%\\end{table\n');
fclose(fid);
OldFolder = cd(NewFolder);
