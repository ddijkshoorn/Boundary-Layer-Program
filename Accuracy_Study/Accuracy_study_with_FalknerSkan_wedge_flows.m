%% Verification: accuracy study of program refining eta-grid comparing solutions of Falkner-Skan wedge flows
% Accuracy Study
% adiabatic nonzero pressure gradient flow with constant density and
% constant fluid properties
% last changes on 28-11-2019
% Last run on 05-05-2020
% Simulation takes about 10 minutes
% Adapted slightly for upload and checked: 21-01-2022

%% Note: only for single station calculations
% calculated on compr FS grid (since we want to study the accuracy of the
% current code)
% transformed to Illingworth-Levy coordinates (because of Cebeci (1974) and
% Rogers (1992)
% one etaE = 8.0?

clear all
close all
clc

c1 = clock;

%% INPUT
run('../Verification_Study/INPUT/INPUT_Verification_with_Rogers1992')     % case specific INPUT (options and settings; INPUT-file is standardized for case C-25)
% is suitable

% Change of input (compared with Table C-25 and C-26, where INPUT-file is standardized for case C-25):
OPT.COMP = 0;
OPT.BCEE = 0;

%% Initialization

% input
% Beta = [1.0; 0.5; 0.3; 0.1; 0; -0.1; -0.1988376];
Beta = ["1.0"; "0.5"; "0.3"; "0.1"; "0"; "-0.1"; "-0.1988376"];
% h1 = [1; 0.5; 0.2; 0.1; 0.01; 0.001; 0.0001; 0.00001]; % if this vector is changed, change table header accordingly manually
h1 = [1; 0.5; 0.2; 0.1; 0.01; 0.001; 0.0001; 0.00001];%; 0.0001; 0.00001]; % if this vector is changed, change table header accordingly manually

% initialization
Size = length(Beta);
PGm2 = zeros(Size,1);
% Deta = zeros(Size,1);
% etaE = zeros(Size,1);
% etaR = zeros(Size,1);
CSM_Data = NaN*ones(Size,6); %%% f''(0), J1, J2, J3, J4, J5
sep = NaN*ones(length(h1),1);
npoints = zeros(length(h1),1);
% kk = NaN*ones(Size,1);
NS = 1;
NST = length(INP.x);
TCC = []; % only laminar flow
MON.SEP = 0;
MON.tr = 0;
SET.ITMAX0 = 20;    %%% added (21-01-2022)
OPT.RLAM = 0;       %%% added (21-01-2022)

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

%% Initial profiles
% calculated inside loop for every new Beta again, since dependent on Beta

%% Main Calculation

for ii = 1:length(h1)

    % Grid
    GRD.Deta(1) = h1(ii);
    if ii == 8
        SET.NPT = 1000000;
    end
    [NP,GRD] = GRID(GRD,SET);
    npoints(ii) = NP;
    
    for jj = 1:Size % size of Beta
        
        PGm2(jj) = str2double(Beta(jj))*(1+(FLD.gamma - 1)/2*INP.MaI^2)^-1/(2 - str2double(Beta(jj))*(1+(FLD.gamma - 1)/2*INP.MaI^2)^-1); % Mach = 0;
        
        % Extraordinary cases (Beta_dact(jj) = 2.0)
        if str2double(Beta(jj)) == 2.0
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
        
%         % initial profiles %%% (28-10-2019) note that solprev is initialized here and is empty []
%         %         if ~exist('sol','var')
%         [sol,solprev] = IVPL(NS,GRD,HVR,FLD); % needs to be here because of use of HVR
%         %         end
        
%%%%%%%%% Circumvene stop due to SET.ITMAX0 in IVPL; copy IVPL here: %%%%%%%%%%%%%%%%%%%%%%%
        NP = GRD.NP; % since NP is a local variable in MAIN-file
        etaNPQ = 0.25*GRD.eta(NP);
        etaU15 = 1.5/GRD.eta(NP);
        
        for j = 1:NP
            etaR = GRD.eta(j)/GRD.eta(NP);
            etaR2 = etaR^2;
            sol.f(j,1) = etaNPQ*etaR2*(3 - 0.5*etaR2);
            sol.u(j,1) = 0.5*etaR*(3 - etaR2);
            sol.v(j,1) = etaU15*(1 - etaR2);
            sol.g(j,1) = 1;
            sol.p(j,1) = 0;
            sol.b(j,1) = 1;
            sol.c(j,1) = 1;
            sol.d(j,1) = 0;
            sol.e(j,1) = 1/FLD.Pr; % ideal constant Prandtl-number
        end
        
        solprev = []; % not used for first station (remains empty, but general initialization is needed for COEF-file)

%         figure(jj)  % #1 checking convergence of initial profile calculation
%         hold on  % #1 checking convergence of initial profile calculation
        IT = 0;
        DvW = 1;
        while abs(DvW) > 1e-5 && IT < SET.ITMAX0 %% objective: 20 iterations when convergence not reached (imitate Return function here, since I don't want to comment the keyboard function in IVPL)
            IT = IT + 1;
            if IT > SET.ITMAX0 - 1 % set in case specific INPUT-file
                fprintf('Calculation of initial profiles did not converge.\n')
%                 keyboard
%                 return
            end
            [S,B,R,sol] = COEF(IT,NS,NP,GRD,HVR,sol,solprev);
            [sol, DvW] = SOLV5(NP,S,B,R,GRD,HVR,sol);
%             plot(sol.v,'DisplayName',sprintf('%d',IT)) % #1 checking convergence of initial profile calculation
        end
%         legend('show') % #1 checking convergence of initial profile calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Run code (note: sol and solprev remain the same outside CSM-file and are re-used)
        [BLC,FLP,SOL,TCC,MON] = CSM(X,NP,NS,NST,GRD,HVR,EDG,FRS,FLD,TCC,sol,solprev,SOL,BLC,INP,OPT,SET,MON);
        
        % Store results (in table format/size, which is smaller than Table from Rogers)
        if MON.SEP == 1
%             CSM_Data(jj + (ii - 1)*Size,1) = NaN;
            CSM_Data(jj + (ii - 1)*Size,1) = NaN;
            CSM_Data(jj + (ii - 1)*Size,2) = NaN;
            CSM_Data(jj + (ii - 1)*Size,3) = NaN;
            CSM_Data(jj + (ii - 1)*Size,4) = NaN;
            CSM_Data(jj + (ii - 1)*Size,5) = NaN;
            CSM_Data(jj + (ii - 1)*Size,6) = NaN;
            sep(ii,1) = 1;
        else
%             TEST(jj + (ii - 1)*Size,1) = SOL{end}.v(1);
            CSM_Data(jj + (ii - 1)*Size,1) = SOL{end}.v(1); % storing shear parameter; or in transformed form: SOL{end}.v(1)*sqrt(2/(m+1));
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
            end % NB H = J5/J2 (Form Factor)
            J1 = SUM1; % NOTE! NOT: in Illingworth-Levy transformed plane (/sqrt(2*FLD.C/(PGm2(jj)+1))); BUT: J1_tr = J1*2/sqrt(2*FLD.C/(PGm2(jj)+1))
            J2 = SUM2; % NOTE! NOT: in Illingworth-Levy transformed plane (/sqrt(2*FLD.C/(PGm2(jj)+1))); BUT: J2_tr = J2*2/sqrt(2*FLD.C/(PGm2(jj)+1))
            J3 = SUM3; % NOTE! NOT: in Illingworth-Levy transformed plane (/sqrt(2*FLD.C/(PGm2(jj)+1))); BUT: J3_tr = J4*2/sqrt(2*FLD.C/(PGm2(jj)+1))
            J4 = SUM4; % NOTE! NOT: in Illingworth-Levy transformed plane (/sqrt(2*FLD.C/(PGm2(jj)+1))); BUT: J4_tr = J?*2/sqrt(2*FLD.C/(PGm2(jj)+1))
            J5 = SUM5; % added (25-10-2019)
            CSM_Data(jj + (ii - 1)*Size,2) = J1;
            CSM_Data(jj + (ii - 1)*Size,3) = J2;
            CSM_Data(jj + (ii - 1)*Size,4) = J3;
            CSM_Data(jj + (ii - 1)*Size,5) = J4;
            CSM_Data(jj + (ii - 1)*Size,6) = J5;
            sep(ii,1) = 0;
        end
        
        % Print results
        if MON.SEP == 1
            if  OPT.BCEE == 0 % adiabatic: derivative of enthalpy ratio = 0 [-] (HVR.alpha0 = 0;)
                fprintf('Separation occurred. No results for input Beta(jj) = %2.6f \n', str2double(Beta(jj)), INP.BCW)
            else
                fprintf('Separation occurred. No results for input Beta(jj) = %2.6f \n', str2double(Beta(jj)), gw(jj))
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
    
    % instead, add ii*length(Beta) to jj in CSM_Data(jj,2)
%     CSM_Results{ii} = CSM_Data;
    
%     clear CSM_Data
    
end

c2 = clock;

% keyboard

%% Table - Write data to Tex-file (in Table format)

NewFolder = cd; % currently in 'boundary-layer-code'-folder
cd('Thesis_Report_Tables') %%% new, added (06-01-2022)
fid = fopen('Accuracy_Study_FalknerSkan.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Accuracy Study with Falkner-Skan wedge flows (adiabatic incompressible nonzero pressure gradient flows with constant fluid properties)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Accuracy Study with Falkner-Skan wedge flows (adiabatic incompressible similar flows with nonzero constant pressure gradients for calorically perfect ideal gas with constant fluid properties) following the example of Cebeci\\cite{cebeci1974analysis} table 8-1. The convergence criterion used was $\\left| \\mathrm{d} f''''(0) \\right| < 1e-5$. The values obtained with the CS-method (CSM) are transformed from the compressible Falkner-Skan transformed y-coordinate with uniform (vertical) grid spacing and height of $\\eta_{\\mathrm{e}} = %1.1f$ to the Illingworth-Levy coordinates.}\n',GRD.etaE);
fprintf(fid, '    \\label{tab:FSA}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=-1.7] S[table-format=-1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\beta$}                                           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$m_2$}                                             &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$f''''_{\\mathrm{w}}(\\mathrm{d} \\eta = %1.1f)$}         &\n',h1(1));
fprintf(fid, '        \\multicolumn{1}{c}{$f''''_{\\mathrm{w}}(\\mathrm{d} \\eta = %1.1f)$}         &\n',h1(2));
fprintf(fid, '        \\multicolumn{1}{c}{$f''''_{\\mathrm{w}}(\\mathrm{d} \\eta = %1.1f)$}         &\n',h1(3));
fprintf(fid, '        \\multicolumn{1}{c}{$f''''_{\\mathrm{w}}(\\mathrm{d} \\eta = %1.1f)$}         &\n',h1(4));
fprintf(fid, '        \\multicolumn{1}{c}{$f''''_{\\mathrm{w}}(\\mathrm{d} \\eta = %1.2f)$}        &\n',h1(5));
fprintf(fid, '        \\multicolumn{1}{c}{$f''''_{\\mathrm{w}}(\\mathrm{d} \\eta = %1.3f)$}        &\n',h1(6));
fprintf(fid, '        \\multicolumn{1}{c}{$f''''_{\\mathrm{w}}(\\mathrm{d} \\eta = %1.4f)$}        &\n',h1(7));
fprintf(fid, '        \\multicolumn{1}{c}{$f''''_{\\mathrm{w}}(\\mathrm{d} \\eta = %1.5f)$}       \\\\\n',h1(8));
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:Size
    if str2double(Beta(ll)) == 2.0
        fprintf(fid, '         %s       &   %4.0f\\tnote{*}   &   %1.6f   &   %1.6f   &   %1.6f   &   %1.6f   &   %1.6f   &   %1.6f  &   %1.6f   &   %1.6f  \\\\\n',Beta(ll),PGm2(ll),CSM_Data(ll,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)),CSM_Data(ll + Size,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)),CSM_Data(ll + 2*Size,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)),CSM_Data(ll + 3*Size,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)),CSM_Data(ll + 4*Size,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)),CSM_Data(ll + 5*Size,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)),CSM_Data(ll + 6*Size,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)),CSM_Data(ll + 7*Size,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)));
    else
        fprintf(fid, '         %s       &   %1.6f             &   %1.6f   &   %1.6f   &   %1.6f   &   %1.6f   &   %1.6f   &   %1.6f  &   %1.6f   &   %1.6f  \\\\\n',Beta(ll),PGm2(ll),CSM_Data(ll,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)),CSM_Data(ll + Size,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)),CSM_Data(ll + 2*Size,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)),CSM_Data(ll + 3*Size,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)),CSM_Data(ll + 4*Size,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)),CSM_Data(ll + 5*Size,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)),CSM_Data(ll + 6*Size,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)),CSM_Data(ll + 7*Size,1)*sqrt(2*FLD.C/(PGm2(ll) + 1)));
    end
end
fprintf(fid, '        \\addlinespace\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{n}_{\\mathrm{points}}$}   &      &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   \\\\\n',num2str(npoints(1)),num2str(npoints(2)),num2str(npoints(3)),num2str(npoints(4)),num2str(npoints(5)),num2str(npoints(6)),num2str(npoints(7)),num2str(npoints(8)));
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

%% Table - form factor H

fid = fopen('Accuracy_Study_FalknerSkan2.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Accuracy Study with Falkner-Skan wedge flows (adiabatic incompressible nonzero pressure gradient flows with constant fluid properties)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Accuracy Study of form factor $H$ with Falkner-Skan wedge flows (adiabatic incompressible similar flows with nonzero constant pressure gradients for calorically perfect ideal gas with constant fluid properties) following the example of Cebeci\\cite{cebeci1974analysis} table 8-1. The convergence criterion used was $\\left| \\mathrm{d} f''''(0) \\right| < 1e-5$. The values obtained with the CS-method (CSM) are transformed from the compressible Falkner-Skan transformed y-coordinate with uniform (vertical) grid spacing and height of $\\eta_{\\mathrm{e}} = %1.1f$ to the Illingworth-Levy coordinates.}\n',GRD.etaE);
fprintf(fid, '    \\label{tab:FSA2}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=-1.7] S[table-format=-1.6] S[table-format=1.4] S[table-format=1.4] S[table-format=1.4] S[table-format=1.4] S[table-format=1.4] S[table-format=1.4] S[table-format=1.4] S[table-format=1.4]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\beta$}                            &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$m_2$}                              &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$H(\\mathrm{d} \\eta = %1.1f)$}         &\n',h1(1));
fprintf(fid, '        \\multicolumn{1}{c}{$H(\\mathrm{d} \\eta = %1.1f)$}         &\n',h1(2));
fprintf(fid, '        \\multicolumn{1}{c}{$H(\\mathrm{d} \\eta = %1.1f)$}         &\n',h1(3));
fprintf(fid, '        \\multicolumn{1}{c}{$H(\\mathrm{d} \\eta = %1.1f)$}         &\n',h1(4));
fprintf(fid, '        \\multicolumn{1}{c}{$H(\\mathrm{d} \\eta = %1.2f)$}        &\n',h1(5));
fprintf(fid, '        \\multicolumn{1}{c}{$H(\\mathrm{d} \\eta = %1.3f)$}        &\n',h1(6));
fprintf(fid, '        \\multicolumn{1}{c}{$H(\\mathrm{d} \\eta = %1.4f)$}        &\n',h1(7));
fprintf(fid, '        \\multicolumn{1}{c}{$H(\\mathrm{d} \\eta = %1.5f)$}       \\\\\n',h1(8));
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:Size
    if str2double(Beta(ll)) == 2.0
        fprintf(fid, '         %s       &   %4.0f\\tnote{*}   &   %1.4f   &   %1.4f   &   %1.4f   &   %1.4f   &   %1.4f   &   %1.4f  &   %1.4f   &   %1.4f  \\\\\n',Beta(ll),PGm2(ll),CSM_Data(ll,6)/CSM_Data(ll,3),CSM_Data(ll + Size,6)/CSM_Data(ll + Size,3),CSM_Data(ll + 2*Size,6)/CSM_Data(ll + 2*Size,3),CSM_Data(ll + 3*Size,6)/CSM_Data(ll + 3*Size,3),CSM_Data(ll + 4*Size,6)/CSM_Data(ll + 4*Size,3),CSM_Data(ll + 5*Size,6)/CSM_Data(ll + 5*Size,3),CSM_Data(ll + 6*Size,6)/CSM_Data(ll + 6*Size,3),CSM_Data(ll + 7*Size,6)/CSM_Data(ll + 7*Size,3));
    else
        fprintf(fid, '         %s       &   %1.6f             &   %1.4f   &   %1.4f   &   %1.4f   &   %1.4f   &   %1.4f   &   %1.4f  &   %1.4f   &   %1.4f  \\\\\n',Beta(ll),PGm2(ll),CSM_Data(ll,6)/CSM_Data(ll,3),CSM_Data(ll + Size,6)/CSM_Data(ll + Size,3),CSM_Data(ll + 2*Size,6)/CSM_Data(ll + 2*Size,3),CSM_Data(ll + 3*Size,6)/CSM_Data(ll + 3*Size,3),CSM_Data(ll + 4*Size,6)/CSM_Data(ll + 4*Size,3),CSM_Data(ll + 5*Size,6)/CSM_Data(ll + 5*Size,3),CSM_Data(ll + 6*Size,6)/CSM_Data(ll + 6*Size,3),CSM_Data(ll + 7*Size,6)/CSM_Data(ll + 7*Size,3));
    end
end
fprintf(fid, '        \\addlinespace\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{n}_{\\mathrm{points}}$}   &      &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   \\\\\n',num2str(npoints(1)),num2str(npoints(2)),num2str(npoints(3)),num2str(npoints(4)),num2str(npoints(5)),num2str(npoints(6)),num2str(npoints(7)),num2str(npoints(8)));
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
