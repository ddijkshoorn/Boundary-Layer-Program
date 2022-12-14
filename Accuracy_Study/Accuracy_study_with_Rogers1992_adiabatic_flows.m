%% Accuracy Study with Rogers (1992) table C-28
% Nonzero pressure gradient in adiabatic flows for Pr=0.723, C=1 (C needs to be a constant), Ma=inf (sigma=2), calorically perfect ideal gas.
% Similar solutions
% 13-12-2019 D.D.D.
% last changes on 16-12-2019
% Last run on 05-05-2020
% Simulation takes about 3 minutes
% Adapted slightly for upload and checked: 23-01-2022

%% Note: only for single station calculations
%%% string-replace is used to print scientific notation (for small and large m2);
%%% for the first m2 (Beta_dact=1.0; sep) strrep is changed to substitute for a - (instead of + for Beta_dact=2.0)

clear all
close all
clc

c1 = clock;

%% INPUT
load('../Verification_Study/Data/VerData_Rogers1992_TableC28')     % format string
run('../Verification_Study/INPUT/INPUT_Verification_with_Rogers1992')      % case specific INPUT (options and settings)
deta = ["1.0"; "0.5"; "0.2"; "0.1"; "0.01"; "0.001"; "0.0001"]; % h1 = [1; 0.5; 0.2; 0.1; 0.01; 0.001; 0.0001; 0.00001];

% Change of input (compare with Table C-25, C-26 and C-27, where INPUT-file is standardized for case C-25):
FLD.Pr = 0.723;
INP.MaI = 1000; % since not possible to simulate Mach = inf
INP.MaE = INP.MaI;

% Change of options
OPT.BCEE = 0; % adiabatic (wall) flow

% Settings
SET.ITMAX = 40;
SET.ITMAX0 = 20;    %%% added (21-01-2022)
OPT.RLAM = 0;       %%% added (21-01-2022)

%% Initialization
Size = length(VerData_Rogers1992_TableC28);
Beta_dact = zeros(Size,1);
PGm2 = zeros(Size,1);
% gw = zeros(Size,1);
fpp0 = NaN*ones(Size,length(deta));
gp0 = NaN*ones(Size,length(deta));
MON_SEP = zeros(Size,length(deta));
CSM_DataTableC28 = NaN*ones(Size,6); %%% f''(0), g'(0) or g(0), J1, J2, J3, J4
sep = NaN*ones(Size,1);
kk = NaN*ones(Size,1);
NS = 1;
NST = length(INP.x);
TCC = []; % only laminar flow
MON.SEP = 0; %% added on 03-09-2019
MON.tr = 0;

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

% % Grid
% [NP,GRD] = GRID(GRD,SET);

%% Initial profiles
% [sol,solprev] = IVPL(NS,GRD,HVR,FLD);
% not applied here since initial profile depends on gw and Beta

%% Main Calculation
ii = 0; % table line (row) index

for hh = 1:length(deta)
    
    GRD.Deta(1) = str2double(deta(hh));
    if hh == 7
        SET.NPT = 1e5;%1e6;
    elseif hh == 6
        SET.NPT = 1e5;
%     else
%         SET.NPT = 1000;
    end
    
    % Grid
    [NP,GRD] = GRID(GRD,SET);
    
    for jj = 1:Size %% pressure gradients, from first to last line of table C-28 Rogers (1992)
        
        Beta_dact(jj) = str2double(VerData_Rogers1992_TableC28(jj,1)); %%% changed for Table C-28!
        %     gw(jj) = str2double(VerData_Rogers1992_TableC28(jj,1));
        PGm2(jj) = Beta_dact(jj)*(1+(FLD.gamma - 1)/2*INP.MaI^2)^-1/(2 - Beta_dact(jj)*(1+(FLD.gamma - 1)/2*INP.MaI^2)^-1);
        
        if jj > 1 && (Beta_dact(jj) > Beta_dact(jj - 1) && sign(Beta_dact(jj)) == sign(Beta_dact(jj - 1))) %|| gw(jj) == 1.0 % 1) 2nd solution, %%% not applicable here, since Pr not equal to 1 and Ma not equal to zero: or 2) since already taken into account by other table (Falkner-Skan)
            % skip line
            %         if gw(jj) == 1.0
            %             ii = ii + 1; % next line in table
            %             kk(ii) = jj; % store location entry
            %             fprintf('Input gw = %7.9f is skipped since it is already included in Table C-1.\n', gw(jj))
            %         else
            fprintf('Input Beta_dact = %7.9f is skipped since it is part of the second solution for decelarating flow which cannot be calculated with this program.\n', Beta_dact(jj))
            %         end
        else % run calculation
            ii = ii + 1; % next line in table
            kk(ii) = jj; % store location entry
            
            % Extraordinary cases (gw = 0.0 or Beta_dact(jj) = 2.0)
            % % %         if Beta_dact(jj) == 2 || gw(jj) == 0
            if Beta_dact(jj) == 2
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
            
            % Grid (preferrably only one grid calculation)
            %         [NP,GRD] = GRID(GRD,SET);
            
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
%             [sol,solprev] = IVPL(NS,GRD,HVR,FLD); % needs to be here because of use of HVR
            %         end

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
                CSM_DataTableC28(ii,1) = NaN;
                CSM_DataTableC28(ii,2) = NaN;
                CSM_DataTableC28(ii,3) = NaN; % #1
                CSM_DataTableC28(ii,4) = NaN; % #1
                CSM_DataTableC28(ii,5) = NaN; % #1
                CSM_DataTableC28(ii,6) = NaN; % #1
                sep(ii,1) = 1;
            else
                CSM_DataTableC28(ii,1) = SOL{end}.v(1); % storing shear parameter; or in transformed form: SOL{end}.v(1)*sqrt(2/(m+1));
                if  OPT.BCEE == 0 % adiabatic: derivative of enthalpy ratio = 0 [-]
                    CSM_DataTableC28(ii,2) = SOL{end}.g(1); % storing adiabatic wall enthalpy ratio
                else
                    CSM_DataTableC28(ii,2) = SOL{end}.p(1); % storing HT parameter; or in transformed form: SOL{end}.p(1)*sqrt(2/(m+1));
                end
                %%% Added for Extended table (J1, J2, J3) #1
                % Caclulate J1, J2, J3, J4 and J5? #1
                TERMP1 = 0;
                TERMP2 = 0;     % momentum thickness
                TERMP3 = 0;
                TERMP4 = 0;
                SUM1 = 0;
                SUM2 = 0;       % momentum thickness
                SUM3 = 0;
                SUM4 = 0;
                for j = 2:NP
                    TERM1 = SOL{end}.g(j) - SOL{end}.u(j)^2;
                    TERM2 = SOL{end}.u(j)*(1 - SOL{end}.u(j));     % momentum thickness
                    TERM3 = 1 - SOL{end}.g(j);
                    TERM4 = 1 - SOL{end}.u(j)^2;
                    SUM1 = SUM1 + GRD.A(j)*(TERM1 + TERMP1);
                    SUM2 = SUM2 + GRD.A(j)*(TERM2 + TERMP2);     % momentum thickness
                    SUM3 = SUM3 + GRD.A(j)*(TERM3 + TERMP3);
                    SUM4 = SUM4 + GRD.A(j)*(TERM4 + TERMP4);
                end
                J1 = SUM1; % NOTE! NOT: in Illingworth-Levy transformed plane (/sqrt(2*FLD.C/(PGm2(jj)+1))); BUT: J1_tr = J1*2/sqrt(2*FLD.C/(PGm2(jj)+1))
                J2 = SUM2; % NOTE! NOT: in Illingworth-Levy transformed plane (/sqrt(2*FLD.C/(PGm2(jj)+1))); BUT: J2_tr = J2*2/sqrt(2*FLD.C/(PGm2(jj)+1))
                J3 = SUM3; % NOTE! NOT: in Illingworth-Levy transformed plane (/sqrt(2*FLD.C/(PGm2(jj)+1))); BUT: J3_tr = J4*2/sqrt(2*FLD.C/(PGm2(jj)+1))
                J4 = SUM4; % NOTE! NOT: in Illingworth-Levy transformed plane (/sqrt(2*FLD.C/(PGm2(jj)+1))); BUT: J4_tr = J?*2/sqrt(2*FLD.C/(PGm2(jj)+1))
                CSM_DataTableC28(ii,3) = J1;
                CSM_DataTableC28(ii,4) = J2;
                CSM_DataTableC28(ii,5) = J3;
                CSM_DataTableC28(ii,6) = J4;
                % % %             *sqrt(2*FLD.C/(PGm2(mm) + 1))
                %%% Added for Extended table (J1, J2, J3) #1
                %             CSM_DataTableC28(ii,1) = SOL{end}.v(1); % storing shear parameter; or in transformed form: SOL{end}.v(1)*sqrt(2/(m+1));
                %             CSM_DataTableC28(ii,2) = SOL{end}.p(1); % storing HT parameter; or in transformed form: SOL{end}.p(1)*sqrt(2/(m+1));
                sep(ii,1) = 0;
            end
            
            % Print results
            if MON.SEP == 1
                %             if HVR.alpha0 > 0
                %                 fprintf('Separation occurred. No results for input Beta_dact(jj) = %2.6f and gw = %1.2f \n', Beta_dact(jj), gw(jj))
                %             else
                %                 fprintf('Separation occurred. No results for input Beta_dact(jj) = %2.6f and g''w = %1.2f \n', Beta_dact(jj), INP.BCW)
                %             end
                if  OPT.BCEE == 0 % adiabatic: derivative of enthalpy ratio = 0 [-] (HVR.alpha0 = 0;)
                    fprintf('Separation occurred. No results for input Beta_dact(jj) = %2.6f and g''w = %1.2f \n', Beta_dact(jj), INP.BCW)
                else
                    fprintf('Separation occurred. No results for input Beta_dact(jj) = %2.6f and gw = %1.2f \n', Beta_dact(jj), gw(jj))
                end
            else
                fprintf('Input Pressure gradient m2 = %7.9f.\n', PGm2(jj))
                fprintf('Input Pressure gradient Beta_dact = %7.9f.\n', 2*PGm2(jj)/(PGm2(jj) + 1)*(1+(FLD.gamma - 1)/2*INP.MaI^2))
                %             if HVR.alpha0 > 0
                %                 fprintf('Input gw = %7.9f.\n', gw(jj))
                %             else
                %                 fprintf('Input g''w = %7.9f.\n', INP.BCW)
                %             end
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
            
            if MON.SEP < 1
                fpp0(jj,hh) = SOL{end}.v(1);
                g0(jj,hh) = SOL{end}.g(1);
                MON_SEP(jj,hh) = MON.SEP;
            else
                fpp0(jj,hh) = NaN;
                g0(jj,hh) = NaN;
            end
            
            % clear output
            clear BLC FLP SOL MON
            
            % and redefine for next run
            MON.SEP = 0; %% added on 12-09-2019
            MON.tr = 0;
            
        end
        
    end
    
end

c2 = clock;

% keyboard

ii = ii/length(deta); % because of multiple simulations per table line

%% Table C-28 - Write data to Tex-file (in Table format)

NewFolder = cd;
cd('Thesis_Report_Tables')
fid = fopen('Accuracy_Study_Rogers1992_TableC28_gw.tex','w');
% cd 'C:\Users\Dominic Dijkshoorn\Documents\boundary-layer-code\tables_temporary' % temporary
% fid = fopen('Accuracy_Study_Rogers1992_TableC28_gw_test.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Accuracy study with Table C-28 from Rogers\\cite{rogers1992laminar}: nonzero pressure gradient adiaabtic flows with Pr=0.723, C=1 (constant), sigma=2.0 (Ma=inf)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Accuracy study with adiabatic laminar compressible similar flows with constant nonzero pressure gradients for calorically perfect ideal gas with $C = 1$ (constant) and $\\mathrm{Pr} = 0.723$ and $\\bar{\\sigma} = 2.0 \\quad (\\mathrm{Ma} = \\infty)$ based on tabulated data obtained from Rogers\\cite{rogers1992laminar} table C-28 and the method following the example of Cebeci\\cite{cebeci1974analysis} table 8-3. The convergence criterion used was $\\left| \\mathrm{d} f''''(0) \\right| < 1e-5$. The values obtained with the CS-method (CSM) are in transformed compressible Falkner-Skan coordinates with grid height of $\\eta_{\\mathrm{e}} = %1.1f$. Note that separation occurred when the table entry shows ''sep''.}\n',GRD.etaE);
fprintf(fid, '    \\label{tab:C28Acc}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=-1.6] S[table-format=-1.6e-1] S[table-format=1.6] S[table-format=1.6] S[table-format=-1.6] S[table-format=-1.6] S[table-format=1.6] S[table-format=-1.6] S[table-format=-1.6]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\hat{\\beta}$}                	&\n');
fprintf(fid, '        \\multicolumn{1}{c}{$m_2$}                          	&\n');
fprintf(fid, '        \\multicolumn{7}{c}{$g_{\\mathrm{aw}}$}             	\\\\\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                              	&\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                             	&\n');
fprintf(fid, '        \\multicolumn{7}{c}{[-]}                              \\\\\n');
fprintf(fid, '        \\cmidrule(lr){3-9}\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                 &\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                 &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}      	&\n',deta(1));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}      	&\n',deta(2));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(3));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(4));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(5));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(6));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}        \\\\\n',deta(7));
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:ii % ii is counter inside loop above (counting valid table rows), ll is index of resulting CSM_data, mm is index of table data from Rogers 1992
    mm = kk(ll);
%     if gw(mm) == 1
%         fprintf(fid, '         \\multicolumn{1}{l}{%s} & \\multicolumn{6}{l}{Solution to the Falkner-Skan equations, see Table \\ref{tab:FSF}}   \\\\\n',VerData_Rogers1992_TableC28(mm,1));
%     elseif Beta_dact(mm) == 2.0
%     if Beta_dact(mm) == 2.0
%         fprintf(fid, '         %s       &   %s       &   %4.0f\\tnote{1}          &  %s   &  %1.6f   &  %s   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC28(mm,1),VerData_Rogers1992_TableC28(mm,2),PGm2(mm),VerData_Rogers1992_TableC28(mm,3),CSM_DataTableC28(ll,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC28(mm,4),CSM_DataTableC28(ll,2)*sqrt(2*FLD.C/(PGm2(mm) + 1)));
%     elseif sep(ll) == 1 % separation has occurred, print 'sep'
%     if sep(ll) == 1 % separation has occurred, print 'sep'
%         fprintf(fid, '         %s       &      %s               &  %s   &  \\multicolumn{1}{c}{sep}   &  %s   &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Rogers1992_TableC28(mm,1),strrep(num2str(PGm2(mm),'%1.6e'), 'e-0','e-'),VerData_Rogers1992_TableC28(mm,2),VerData_Rogers1992_TableC28(mm,3));
%     else
%         fprintf(fid, '         %s       &      %s               &  %s   &  %1.6f   &  %s   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC28(mm,1),strrep(num2str(PGm2(mm),'%1.6e'), 'e+0','e+'),VerData_Rogers1992_TableC28(mm,2),CSM_DataTableC28(ll,2),VerData_Rogers1992_TableC28(mm,3),CSM_DataTableC28(ll,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)));
%     end
    if sep(ll) == 1 % separation has occurred, print 'sep'
        fprintf(fid, '         %s       &      %s               &  \\multicolumn{1}{c}{sep}   &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Rogers1992_TableC28(mm,1),num2str(PGm2(mm),'%1.6e'));
    else
        fprintf(fid, '         %s       &      %s               &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC28(mm,1),num2str(PGm2(mm),'%1.6e'),g0(ll,1),g0(ll,2),g0(ll,3),g0(ll,4),g0(ll,5),g0(ll,6),g0(ll,7));
    end
end
% fprintf(fid, '        %% Data\n');
fprintf(fid, '        \\bottomrule\n\n');
fprintf(fid, '    \\end{tabular}\n');
% fprintf(fid, '    \\begin{tablenotes}\n');
% fprintf(fid, '        \\item[1] The solution does not converge for $m_2 = \\infty$, and therefore $m_2 = 1000$ is taken as approach.\n');
% fprintf(fid, '    \\end{tablenotes}\n');
fprintf(fid, '\\end{threeparttable}\n');
% fprintf(fid, '%%\\end{sidewaystable}\n');
% fprintf(fid, '%%\\end{table\n');
fclose(fid);
% OldFolder = cd(NewFolder);

%% Table C-28 - Appendix - Write data to Tex-file (in Table format)

% NewFolder = cd;
% cd('Thesis_Report_Tables')
fid = fopen('Accuracy_Study_Rogers1992_TableC28_gw_A.tex','w');
% cd 'C:\Users\Dominic Dijkshoorn\Documents\boundary-layer-code\tables_temporary' % temporary
% fid = fopen('Accuracy_Study_Rogers1992_TableC28_gw_test.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Accuracy study with Table C-28 from Rogers\\cite{rogers1992laminar}: nonzero pressure gradient adiaabtic flows with Pr=0.723, C=1 (constant), sigma=2.0 (Ma=inf)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Accuracy study with adiabatic laminar compressible similar flows with constant nonzero pressure gradients for calorically perfect ideal gas with $C = 1$ (constant) and $\\mathrm{Pr} = 0.723$ and $\\bar{\\sigma} = 2.0 \\quad (\\mathrm{Ma} = \\infty)$ based on tabulated data obtained from Rogers\\cite{rogers1992laminar} table C-28 and the method following the example of Cebeci\\cite{cebeci1974analysis} table 8-3. The convergence criterion used was $\\left| \\mathrm{d} f''''(0) \\right| < 1e-5$. The values obtained with the CS-method (CSM) are in transformed compressible Falkner-Skan coordinates with grid height of $\\eta_{\\mathrm{e}} = %1.1f$. Note that separation occurred when the table entry shows ''sep''.}\n',GRD.etaE);
fprintf(fid, '    \\label{tab:C28Acc_A}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=-1.6] S[table-format=-1.6e-1] S[table-format=1.6] S[table-format=1.6] S[table-format=-1.6] S[table-format=-1.6] S[table-format=1.6] S[table-format=-1.6] S[table-format=-1.6]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\hat{\\beta}$}                	&\n');
fprintf(fid, '        \\multicolumn{1}{c}{$m_2$}                          	&\n');
fprintf(fid, '        \\multicolumn{7}{c}{$g_{\\mathrm{aw}}$}             	\\\\\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                              	&\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                             	&\n');
fprintf(fid, '        \\multicolumn{7}{c}{[-]}                              \\\\\n');
fprintf(fid, '        \\cmidrule(lr){3-9}\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                 &\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                 &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}      	&\n',deta(1));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}      	&\n',deta(2));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(3));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(4));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(5));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(6));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}        \\\\\n',deta(7));
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:ii % ii is counter inside loop above (counting valid table rows), ll is index of resulting CSM_data, mm is index of table data from Rogers 1992
    mm = kk(ll);
%     if gw(mm) == 1
%         fprintf(fid, '         \\multicolumn{1}{l}{%s} & \\multicolumn{6}{l}{Solution to the Falkner-Skan equations, see Table \\ref{tab:FSF}}   \\\\\n',VerData_Rogers1992_TableC28(mm,1));
%     elseif Beta_dact(mm) == 2.0
%     if Beta_dact(mm) == 2.0
%         fprintf(fid, '         %s       &   %s       &   %4.0f\\tnote{1}          &  %s   &  %1.6f   &  %s   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC28(mm,1),VerData_Rogers1992_TableC28(mm,2),PGm2(mm),VerData_Rogers1992_TableC28(mm,3),CSM_DataTableC28(ll,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC28(mm,4),CSM_DataTableC28(ll,2)*sqrt(2*FLD.C/(PGm2(mm) + 1)));
%     elseif sep(ll) == 1 % separation has occurred, print 'sep'
%     if sep(ll) == 1 % separation has occurred, print 'sep'
%         fprintf(fid, '         %s       &      %s               &  %s   &  \\multicolumn{1}{c}{sep}   &  %s   &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Rogers1992_TableC28(mm,1),strrep(num2str(PGm2(mm),'%1.6e'), 'e-0','e-'),VerData_Rogers1992_TableC28(mm,2),VerData_Rogers1992_TableC28(mm,3));
%     else
%         fprintf(fid, '         %s       &      %s               &  %s   &  %1.6f   &  %s   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC28(mm,1),strrep(num2str(PGm2(mm),'%1.6e'), 'e+0','e+'),VerData_Rogers1992_TableC28(mm,2),CSM_DataTableC28(ll,2),VerData_Rogers1992_TableC28(mm,3),CSM_DataTableC28(ll,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)));
%     end
    if sep(ll) == 1 % separation has occurred, print 'sep'
        fprintf(fid, '         %s       &      %s               &  \\multicolumn{1}{c}{sep}   &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Rogers1992_TableC28(mm,1),num2str(PGm2(mm),'%1.6e'));
    else
        fprintf(fid, '         %s       &      %s               &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC28(mm,1),num2str(PGm2(mm),'%1.6e'),g0(ll,1),g0(ll,2),g0(ll,3),g0(ll,4),g0(ll,5),g0(ll,6),g0(ll,7));
    end
end
% fprintf(fid, '        %% Data\n');
fprintf(fid, '        \\bottomrule\n\n');
fprintf(fid, '    \\end{tabular}\n');
% fprintf(fid, '    \\begin{tablenotes}\n');
% fprintf(fid, '        \\item[1] The solution does not converge for $m_2 = \\infty$, and therefore $m_2 = 1000$ is taken as approach.\n');
% fprintf(fid, '    \\end{tablenotes}\n');
fprintf(fid, '\\end{threeparttable}\n');
% fprintf(fid, '%%\\end{sidewaystable}\n');
% fprintf(fid, '%%\\end{table\n');
fclose(fid);
% OldFolder = cd(NewFolder);

%% Table C-28 - Write data to Tex-file (in Table format)

% NewFolder = cd;
% cd('Thesis_Report_Tables')
fid = fopen('Accuracy_Study_Rogers1992_TableC28_fppw.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Accuracy study with Table C-28 from Rogers\\cite{rogers1992laminar}: nonzero pressure gradient adiaabtic flows with Pr=0.723, C=1 (constant), sigma=2.0 (Ma=inf)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Accuracy study with adiabatic laminar compressible similar flows with constant nonzero pressure gradients for calorically perfect ideal gas with $C = 1$ (constant) and $\\mathrm{Pr} = 0.723$ and $\\bar{\\sigma} = 2.0 \\quad (\\mathrm{Ma} = \\infty)$ based on tabulated data obtained from Rogers\\cite{rogers1992laminar} table C-28 and the method following the example of Cebeci\\cite{cebeci1974analysis} table 8-3. The convergence criterion used was $\\left| \\mathrm{d} f''''(0) \\right| < 1e-5$. The values obtained with the CS-method (CSM) are in transformed compressible Falkner-Skan coordinates with grid height of $\\eta_{\\mathrm{e}} = %1.1f$. Note that separation occurred when the table entry shows ''sep''.}\n',GRD.etaE);
fprintf(fid, '    \\label{tab:C28Bcc}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=-1.6] S[table-format=-1.6e-1] S[table-format=1.6] S[table-format=1.6] S[table-format=-1.6] S[table-format=-1.6] S[table-format=1.6] S[table-format=-1.6] S[table-format=-1.6]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\hat{\\beta}$}                	&\n');
fprintf(fid, '        \\multicolumn{1}{c}{$m_2$}                          	&\n');
fprintf(fid, '        \\multicolumn{7}{c}{$f''''_{0}$}             	\\\\\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                              	&\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                             	&\n');
fprintf(fid, '        \\multicolumn{7}{c}{[-]}                              \\\\\n');
fprintf(fid, '        \\cmidrule(lr){3-9}\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                 &\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                 &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}      	&\n',deta(1));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}      	&\n',deta(2));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(3));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(4));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(5));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(6));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}        \\\\\n',deta(7));
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:ii % ii is counter inside loop above (counting valid table rows), ll is index of resulting CSM_data, mm is index of table data from Rogers 1992
    mm = kk(ll);
%     if gw(mm) == 1
%         fprintf(fid, '         \\multicolumn{1}{l}{%s} & \\multicolumn{6}{l}{Solution to the Falkner-Skan equations, see Table \\ref{tab:FSF}}   \\\\\n',VerData_Rogers1992_TableC28(mm,1));
%     elseif Beta_dact(mm) == 2.0
%     if Beta_dact(mm) == 2.0
%         fprintf(fid, '         %s       &   %s       &   %4.0f\\tnote{1}          &  %s   &  %1.6f   &  %s   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC28(mm,1),VerData_Rogers1992_TableC28(mm,2),PGm2(mm),VerData_Rogers1992_TableC28(mm,3),CSM_DataTableC28(ll,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC28(mm,4),CSM_DataTableC28(ll,2)*sqrt(2*FLD.C/(PGm2(mm) + 1)));
%     elseif sep(ll) == 1 % separation has occurred, print 'sep'
%     if sep(ll) == 1 % separation has occurred, print 'sep'
%         fprintf(fid, '         %s       &      %s               &  %s   &  \\multicolumn{1}{c}{sep}   &  %s   &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Rogers1992_TableC28(mm,1),strrep(num2str(PGm2(mm),'%1.6e'), 'e-0','e-'),VerData_Rogers1992_TableC28(mm,2),VerData_Rogers1992_TableC28(mm,3));
%     else
%         fprintf(fid, '         %s       &      %s               &  %s   &  %1.6f   &  %s   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC28(mm,1),strrep(num2str(PGm2(mm),'%1.6e'), 'e+0','e+'),VerData_Rogers1992_TableC28(mm,2),CSM_DataTableC28(ll,2),VerData_Rogers1992_TableC28(mm,3),CSM_DataTableC28(ll,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)));
%     end
    if sep(ll) == 1 % separation has occurred, print 'sep'
        fprintf(fid, '         %s       &      %s               &  \\multicolumn{1}{c}{sep}   &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Rogers1992_TableC28(mm,1),num2str(PGm2(mm),'%1.6e'));
    else
        fprintf(fid, '         %s       &      %s               &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC28(mm,1),num2str(PGm2(mm),'%1.6e'),fpp0(ll,1),fpp0(ll,2),fpp0(ll,3),fpp0(ll,4),fpp0(ll,5),fpp0(ll,6),fpp0(ll,7));
    end
end
% fprintf(fid, '        %% Data\n');
fprintf(fid, '        \\bottomrule\n\n');
fprintf(fid, '    \\end{tabular}\n');
% fprintf(fid, '    \\begin{tablenotes}\n');
% fprintf(fid, '        \\item[1] The solution does not converge for $m_2 = \\infty$, and therefore $m_2 = 1000$ is taken as approach.\n');
% fprintf(fid, '    \\end{tablenotes}\n');
fprintf(fid, '\\end{threeparttable}\n');
% fprintf(fid, '%%\\end{sidewaystable}\n');
% fprintf(fid, '%%\\end{table\n');
fclose(fid);
% OldFolder = cd(NewFolder);

%% Table C-28 - Appendix - Write data to Tex-file (in Table format)

% NewFolder = cd;
% cd('Thesis_Report_Tables')
fid = fopen('Accuracy_Study_Rogers1992_TableC28_fppw_A.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Accuracy study with Table C-28 from Rogers\\cite{rogers1992laminar}: nonzero pressure gradient adiaabtic flows with Pr=0.723, C=1 (constant), sigma=2.0 (Ma=inf)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Accuracy study with adiabatic laminar compressible similar flows with constant nonzero pressure gradients for calorically perfect ideal gas with $C = 1$ (constant) and $\\mathrm{Pr} = 0.723$ and $\\bar{\\sigma} = 2.0 \\quad (\\mathrm{Ma} = \\infty)$ based on tabulated data obtained from Rogers\\cite{rogers1992laminar} table C-28 and the method following the example of Cebeci\\cite{cebeci1974analysis} table 8-3. The convergence criterion used was $\\left| \\mathrm{d} f''''(0) \\right| < 1e-5$. The values obtained with the CS-method (CSM) are in transformed compressible Falkner-Skan coordinates with grid height of $\\eta_{\\mathrm{e}} = %1.1f$. Note that separation occurred when the table entry shows ''sep''.}\n',GRD.etaE);
fprintf(fid, '    \\label{tab:C28Bcc_A}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=-1.6] S[table-format=-1.6e-1] S[table-format=1.6] S[table-format=1.6] S[table-format=-1.6] S[table-format=-1.6] S[table-format=1.6] S[table-format=-1.6] S[table-format=-1.6]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\hat{\\beta}$}                	&\n');
fprintf(fid, '        \\multicolumn{1}{c}{$m_2$}                          	&\n');
fprintf(fid, '        \\multicolumn{7}{c}{$f''''_{0}$}             	\\\\\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                              	&\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                             	&\n');
fprintf(fid, '        \\multicolumn{7}{c}{[-]}                              \\\\\n');
fprintf(fid, '        \\cmidrule(lr){3-9}\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                 &\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                 &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}      	&\n',deta(1));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}      	&\n',deta(2));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(3));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(4));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(5));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}       	&\n',deta(6));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}        \\\\\n',deta(7));
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:ii % ii is counter inside loop above (counting valid table rows), ll is index of resulting CSM_data, mm is index of table data from Rogers 1992
    mm = kk(ll);
%     if gw(mm) == 1
%         fprintf(fid, '         \\multicolumn{1}{l}{%s} & \\multicolumn{6}{l}{Solution to the Falkner-Skan equations, see Table \\ref{tab:FSF}}   \\\\\n',VerData_Rogers1992_TableC28(mm,1));
%     elseif Beta_dact(mm) == 2.0
%     if Beta_dact(mm) == 2.0
%         fprintf(fid, '         %s       &   %s       &   %4.0f\\tnote{1}          &  %s   &  %1.6f   &  %s   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC28(mm,1),VerData_Rogers1992_TableC28(mm,2),PGm2(mm),VerData_Rogers1992_TableC28(mm,3),CSM_DataTableC28(ll,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC28(mm,4),CSM_DataTableC28(ll,2)*sqrt(2*FLD.C/(PGm2(mm) + 1)));
%     elseif sep(ll) == 1 % separation has occurred, print 'sep'
%     if sep(ll) == 1 % separation has occurred, print 'sep'
%         fprintf(fid, '         %s       &      %s               &  %s   &  \\multicolumn{1}{c}{sep}   &  %s   &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Rogers1992_TableC28(mm,1),strrep(num2str(PGm2(mm),'%1.6e'), 'e-0','e-'),VerData_Rogers1992_TableC28(mm,2),VerData_Rogers1992_TableC28(mm,3));
%     else
%         fprintf(fid, '         %s       &      %s               &  %s   &  %1.6f   &  %s   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC28(mm,1),strrep(num2str(PGm2(mm),'%1.6e'), 'e+0','e+'),VerData_Rogers1992_TableC28(mm,2),CSM_DataTableC28(ll,2),VerData_Rogers1992_TableC28(mm,3),CSM_DataTableC28(ll,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)));
%     end
    if sep(ll) == 1 % separation has occurred, print 'sep'
        fprintf(fid, '         %s       &      %s               &  \\multicolumn{1}{c}{sep}   &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Rogers1992_TableC28(mm,1),num2str(PGm2(mm),'%1.6e'));
    else
        fprintf(fid, '         %s       &      %s               &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC28(mm,1),num2str(PGm2(mm),'%1.6e'),fpp0(ll,1),fpp0(ll,2),fpp0(ll,3),fpp0(ll,4),fpp0(ll,5),fpp0(ll,6),fpp0(ll,7));
    end
end
% fprintf(fid, '        %% Data\n');
fprintf(fid, '        \\bottomrule\n\n');
fprintf(fid, '    \\end{tabular}\n');
% fprintf(fid, '    \\begin{tablenotes}\n');
% fprintf(fid, '        \\item[1] The solution does not converge for $m_2 = \\infty$, and therefore $m_2 = 1000$ is taken as approach.\n');
% fprintf(fid, '    \\end{tablenotes}\n');
fprintf(fid, '\\end{threeparttable}\n');
% fprintf(fid, '%%\\end{sidewaystable}\n');
% fprintf(fid, '%%\\end{table\n');
fclose(fid);
OldFolder = cd(NewFolder);

%% PLOT Results
