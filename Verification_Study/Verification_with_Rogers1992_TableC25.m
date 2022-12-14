%% Verification with Rogers (1992) table C-25
% Nonzero pressure gradient including heat transfer for Pr=1, C=1 (C needs to be a constant), calorically perfect ideal gas.
% Similar solutions
% 03-09-2019 D.D.D.
% last changes: 24-10-2019
% last run: 04-05-2020
% Simulation takes about 7 minutes
% Adapted slightly for upload and checked: 6-01-2022

%% Note: only for single station calculations

clear all
close all
clc

%% INPUT
% load('Rogers1992_TableC25')                       % format double, not needed anymore
load('./Data/VerData_Rogers1992_TableC25')          % format string
run('./INPUT/INPUT_Verification_with_Rogers1992')   % case specific INPUT (options and settings)

%% Initialization
Size = length(VerData_Rogers1992_TableC25);
Beta_dact = zeros(Size,1);
PGm2 = zeros(Size,1);
gw = zeros(Size,1);
% CSM_DataTableC25 = NaN*ones(Size,2);
CSM_DataTableC25 = NaN*ones(Size,6); %%% f''(0), g'(0) or g(0), J1, J2, J3, J4
sep = NaN*ones(Size,1);
kk = NaN*ones(Size,1);
NS = 1;
NST = length(INP.x);
TCC = []; % only laminar flow
MON.SEP = 0; %% added on 03-09-2019
MON.tr = 0;
SET.ITMAX0 = 20;    %%% added (31-12-2021)
OPT.RLAM = 0;       %%% added (31-12-2021)

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

% Grid
[NP,GRD] = GRID(GRD,SET);

%% Initial profiles
% [sol,solprev] = IVPL(NS,GRD,HVR,FLD);
% not applied here since initial profile depends on gw and Beta

%% Main Calculation
ii = 0; % table line (row) index

for jj = 1:Size

    Beta_dact(jj) = str2double(VerData_Rogers1992_TableC25(jj,2));
    gw(jj) = str2double(VerData_Rogers1992_TableC25(jj,1));
    PGm2(jj) = Beta_dact(jj)*(1+(FLD.gamma - 1)/2*INP.MaI^2)^-1/(2 - Beta_dact(jj)*(1+(FLD.gamma - 1)/2*INP.MaI^2)^-1);
    
    if jj > 1 && (Beta_dact(jj) > Beta_dact(jj - 1) && sign(Beta_dact(jj)) == sign(Beta_dact(jj - 1))) || gw(jj) == 1.0 % 1) 2nd solution, or 2) since already taken into account by other table (Falkner-Skan)
        % skip line
        if gw(jj) == 1.0
            ii = ii + 1; % next line in table
            kk(ii) = jj; % store location entry
            fprintf('Input gw = %7.9f is skipped since it is already included in Table C-1.\n', gw(jj))
        else
            fprintf('Input Beta_dact = %7.9f is skipped since it is part of the second solution for decelarating flow which cannot be calculated with this program.\n', Beta_dact(jj))
        end
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
        
        % initial profiles
%         [sol,solprev] = IVPL(NS,GRD,HVR,FLD,SET); % needs to be here because of use of HVR

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
            CSM_DataTableC25(ii,1) = NaN;
            CSM_DataTableC25(ii,2) = NaN;
            CSM_DataTableC25(ii,3) = NaN; % #1
            CSM_DataTableC25(ii,4) = NaN; % #1
            CSM_DataTableC25(ii,5) = NaN; % #1
            CSM_DataTableC25(ii,6) = NaN; % #1
            sep(ii,1) = 1;
        else
            CSM_DataTableC25(ii,1) = SOL{end}.v(1); % storing shear parameter; or in transformed form: SOL{end}.v(1)*sqrt(2/(m+1));
            CSM_DataTableC25(ii,2) = SOL{end}.p(1); % storing HT parameter; or in transformed form: SOL{end}.p(1)*sqrt(2/(m+1));
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
            CSM_DataTableC25(ii,3) = J1;
            CSM_DataTableC25(ii,4) = J2;
            CSM_DataTableC25(ii,5) = J3;
            CSM_DataTableC25(ii,6) = J4;
            %% Test %%%%%%%%%%%%%%%%%%%%
%             if Beta_dact(jj) == -0.20 && gw(jj) == 0.6
%                 keyboard
%             end
%             %%% correct 'transformation for comparison' %%%%%%%%%%%%%%%%%%%
%             J_tr(ii,1) = J1*2/sqrt(2*FLD.C/(PGm2(jj)+1));%*sqrt(2*FLD.C/(PGm2(mm) + 1))
%             J_tr(ii,2) = J2*2/sqrt(2*FLD.C/(PGm2(jj)+1));
%             J_tr(ii,3) = J3*2/sqrt(2*FLD.C/(PGm2(jj)+1));
%             J_tr(ii,4) = J4*2/sqrt(2*FLD.C/(PGm2(jj)+1));
%             %%% correct 'transformation for comparison' %%%%%%%%%%%%%%%%%%%
%             J_tr(ii,1) = J1/sqrt(2*FLD.C/(PGm2(jj)+1));%*sqrt(2*FLD.C/(PGm2(mm) + 1))
%             J_tr(ii,2) = J2/sqrt(2*FLD.C/(PGm2(jj)+1));
%             J_tr(ii,3) = J3/sqrt(2*FLD.C/(PGm2(jj)+1));
%             J_tr(ii,4) = J4/sqrt(2*FLD.C/(PGm2(jj)+1));
            %% Test %%%%%%%%%%%%%%%%%%%%
% % %             *sqrt(2*FLD.C/(PGm2(mm) + 1))
            %%% Added for Extended table (J1, J2, J3) #1
%             CSM_DataTableC25(ii,1) = SOL{end}.v(1); % storing shear parameter; or in transformed form: SOL{end}.v(1)*sqrt(2/(m+1));
%             CSM_DataTableC25(ii,2) = SOL{end}.p(1); % storing HT parameter; or in transformed form: SOL{end}.p(1)*sqrt(2/(m+1));
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
        end
        
% % %         % Settings back to normal
% % %         if Beta_dact(jj) == 2 || gw(jj) == 0
% % %             SET.NPT = 1000;
% % %             GRD.Deta = 0.01;
% % %         end

        % clear output
        clear BLC FLP SOL MON
        
        % and redefine for next run
        MON.SEP = 0; %% added on 12-09-2019
        MON.tr = 0;
        
    end
    
end

% keyboard

%% Table C-25 - Write data to Tex-file (in Table format)

NewFolder = cd; % currently in 'boundary-layer-code'-folder
cd('Thesis_Report_Tables') %%% new, added (06-01-2021)
fid = fopen('Verification_Rogers1992_TableC25.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Verification with Table C-25 from Rogers\\cite{rogers1992laminar}: Pr=1, C=1 (constant)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
% fprintf(fid, '    \\caption{Comparison with tabulated data of laminar compressible similar flows with constant nonzero pressure gradients and heat transfer for calorically perfect ideal gas with $C = 1$ (constant) and $\\mathrm{Pr} = 1$ taken from Rogers\\cite{rogers1992laminar} table C-25. Transformed uniform (vertical) grid spacing (in Illingworth-Levy transformed y-coordinate) of $\\mathrm{d} \\eta = \\sqrt{\\frac{m_2 + 1}{2C}} %1.4f$ and height of $\\eta_{\\mathrm{e}} = \\sqrt{\\frac{m_2 + 1}{2C}} %1.1f$. Note that separation occurred when the table entry shows ''sep''.}\n',GRD.Deta(1),GRD.etaE);
fprintf(fid, '    \\caption{Comparison with tabulated data of laminar compressible similar flows with constant nonzero pressure gradients and heat transfer for calorically perfect ideal gas with $C = 1$ (constant) and $\\mathrm{Pr} = 1$ taken from Rogers\\cite{rogers1992laminar} table C-25. The values obtained with the CS-method (CSM) are transformed from the compressible Falkner-Skan transformed y-coordinate with uniform (vertical) grid spacing of $\\mathrm{d} \\eta = \\sqrt{\\frac{2C}{m_2 + 1}} %1.4f$ and height of $\\eta_{\\mathrm{e}} = \\sqrt{\\frac{2C}{m_2 + 1}} %1.1f$ to the Illingworth-Levy coordinates ($\\mathrm{d} \\eta = %1.4f$ and $\\eta_{\\mathrm{e}} = %1.1f$). Note that separation occurred when the table entry shows ''sep''.}\n',GRD.Deta(1),GRD.etaE,GRD.Deta(1),GRD.etaE);
fprintf(fid, '    \\label{tab:C25}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=1.2] S[table-format=-1.6] S[table-format=4.6] S[table-format=1.6] S[table-format=1.6] S[table-format=-1.6] S[table-format=-1.6]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$g_{\\mathrm{w}}$}                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\hat{\\beta}$}                         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$m_2$}                                 &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$f''''(0)$}                              &\n');        
fprintf(fid, '        \\multicolumn{2}{c}{$g''(0)$}                               \\\\\n');
fprintf(fid, '        \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          \\\\\n');
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:ii % ii is counter inside loop above (counting valid table rows), ll is index of resulting CSM_data, mm is index of table data from Rogers 1992
    mm = kk(ll);
    if gw(mm) == 1
        fprintf(fid, '         \\multicolumn{1}{l}{%s} & \\multicolumn{6}{l}{Solution to the Falkner-Skan equations, see Table \\ref{tab:FSF}}   \\\\\n',VerData_Rogers1992_TableC25(mm,1));
    elseif Beta_dact(mm) == 2.0
        fprintf(fid, '         %s       &   %s       &   %4.0f\\tnote{*}          &  %s   &  %1.6f   &  %s   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC25(mm,1),VerData_Rogers1992_TableC25(mm,2),PGm2(mm),VerData_Rogers1992_TableC25(mm,3),CSM_DataTableC25(ll,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC25(mm,4),CSM_DataTableC25(ll,2)*sqrt(2*FLD.C/(PGm2(mm) + 1)));
    elseif sep(ll) == 1 % separation has occurred, print 'sep'
        fprintf(fid, '                  &   %s       &   %1.6f            &  %s   &  \\multicolumn{1}{c}{sep}   &  %s   &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Rogers1992_TableC25(mm,2),PGm2(mm),VerData_Rogers1992_TableC25(mm,3),VerData_Rogers1992_TableC25(mm,4));
    else
        fprintf(fid, '                  &   %s       &   %1.6f            &  %s   &  %1.6f   &  %s   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC25(mm,2),PGm2(mm),VerData_Rogers1992_TableC25(mm,3),CSM_DataTableC25(ll,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC25(mm,4),CSM_DataTableC25(ll,2)*sqrt(2*FLD.C/(PGm2(mm) + 1)));
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
% % % OldFolder = cd(NewFolder);

% keyboard

%% Table C-25 extended (J1, J2, J3) - Write data to Tex-file (in Table format) #1

% % % NewFolder = cd;
% cd('Thesis_Report_Tables') %%% new, added (06-01-2021)
fid = fopen('Verification_Rogers1992_TableC25_extended.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Verification with Table C-25 from Rogers\\cite{rogers1992laminar}: Pr=1, C=1 (constant)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '\\%%begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Comparison with tabulated data of laminar compressible similar flows with constant nonzero pressure gradients and heat transfer for calorically perfect ideal gas with $C = 1$ (constant) and $\\mathrm{Pr} = 1$ taken from Rogers\\cite{rogers1992laminar} table C-25. The values obtained with the CS-method (CSM) are transformed from the compressible Falkner-Skan transformed y-coordinate with uniform (vertical) grid spacing of $\\mathrm{d} \\eta = \\sqrt{\\frac{2C}{m_2 + 1}} %1.4f$ and height of $\\eta_{\\mathrm{e}} = \\sqrt{\\frac{2C}{m_2 + 1}} %1.1f$ to the Illingworth-Levy coordinates ($\\mathrm{d} \\eta = %1.4f$ and $\\eta_{\\mathrm{e}} = %1.1f$). Note that separation occurred when the table entry shows ''sep''.}\n',GRD.Deta(1),GRD.etaE,GRD.Deta(1),GRD.etaE);
fprintf(fid, '    \\label{tab:C25E}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=1.2] S[table-format=-1.6] S[table-format=4.6] S[table-format=1.6] S[table-format=1.6] S[table-format=-1.6] S[table-format=-1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$g_{\\mathrm{w}}$}                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\hat{\\beta}$}                         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$m_2$}                                 &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$f''''(0)$}                              &\n');
fprintf(fid, '        \\multicolumn{2}{c}{$g''(0)$}                               &\n');
% fprintf(fid, '        \\multicolumn{2}{c}{$J_1$}                                 &\n');
% fprintf(fid, '        \\multicolumn{2}{c}{$J_2$}                                 &\n');
% fprintf(fid, '        \\multicolumn{2}{c}{$J_3$}                                 \\\\\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          \\\\\n');
fprintf(fid, '        \\cmidrule(lr){4-5} \\cmidrule(lr){6-7}\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                                      &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          &\n');
% fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
% fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          &\n');
% fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
% fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          &\n');
% fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{Rogers}}                       &\n');
% fprintf(fid, '        \\multicolumn{1}{c}{\\textbf{CSM}}                          \\\\\n');
fprintf(fid, '        \\multicolumn{1}{c}{$J_1$}                                 &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$2 J_1$}                               &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$J_2$}                                 &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$2 J_2$}                               &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$J_3$}                                 &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$2 J_4$}                               \\\\\n\n');
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:ii % ii is counter inside loop above (counting valid table rows), ll is index of resulting CSM_data, mm is index of table data from Rogers 1992
    mm = kk(ll);
    if gw(mm) == 1
        fprintf(fid, '         \\multicolumn{1}{l}{%s} & \\multicolumn{12}{l}{Solution to the Falkner-Skan equations, see Table \\ref{tab:FSF}}   \\\\\n',VerData_Rogers1992_TableC25(mm,1));
    elseif Beta_dact(mm) == 2.0
        fprintf(fid, '         %s       &   %s       &   %4.0f\\tnote{*}          &  %s   &  %1.6f   &  %s   &  %1.6f   &  %s   &  %1.6f   &  %s   &  %1.6f   &  %s   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC25(mm,1),VerData_Rogers1992_TableC25(mm,2),PGm2(mm),VerData_Rogers1992_TableC25(mm,3),CSM_DataTableC25(ll,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC25(mm,4),CSM_DataTableC25(ll,2)*sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC25(mm,5),CSM_DataTableC25(ll,3)*2/sqrt(2*FLD.C/(PGm2(mm)+1)),VerData_Rogers1992_TableC25(mm,6),CSM_DataTableC25(ll,4)*2/sqrt(2*FLD.C/(PGm2(mm)+1)),VerData_Rogers1992_TableC25(mm,7),CSM_DataTableC25(ll,6)*2/sqrt(2*FLD.C/(PGm2(mm)+1)));
    elseif Beta_dact(mm) == -0.129507
        fprintf(fid, '                  &   %s       &   %1.6f            &  %s   &  %1.6f   &  %s   &  %1.6f   &  %s\\tnote{$\\diamond$}   &  %1.6f   &  %s\\tnote{$\\diamond$}   &  %1.6f   &  %s   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC25(mm,2),PGm2(mm),VerData_Rogers1992_TableC25(mm,3),CSM_DataTableC25(ll,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC25(mm,4),CSM_DataTableC25(ll,2)*sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC25(mm,5),CSM_DataTableC25(ll,3)*2/sqrt(2*FLD.C/(PGm2(mm)+1)),VerData_Rogers1992_TableC25(mm,6),CSM_DataTableC25(ll,4)*2/sqrt(2*FLD.C/(PGm2(mm)+1)),VerData_Rogers1992_TableC25(mm,7),CSM_DataTableC25(ll,6)*2/sqrt(2*FLD.C/(PGm2(mm)+1)));
    elseif sep(ll) == 1 % separation has occurred, print 'sep'
        fprintf(fid, '                  &   %s       &   %1.6f            &  %s   &  \\multicolumn{1}{c}{sep}   &  %s   &  \\multicolumn{1}{c}{sep}   &  %s   &  \\multicolumn{1}{c}{sep}   &  %s   &  \\multicolumn{1}{c}{sep}   &  %s   &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Rogers1992_TableC25(mm,2),PGm2(mm),VerData_Rogers1992_TableC25(mm,3),VerData_Rogers1992_TableC25(mm,4),VerData_Rogers1992_TableC25(mm,5),VerData_Rogers1992_TableC25(mm,6),VerData_Rogers1992_TableC25(mm,7));
    else
        fprintf(fid, '                  &   %s       &   %1.6f            &  %s   &  %1.6f   &  %s   &  %1.6f   &  %s   &  %1.6f   &  %s   &  %1.6f   &  %s   &  %1.6f  \\\\\n',VerData_Rogers1992_TableC25(mm,2),PGm2(mm),VerData_Rogers1992_TableC25(mm,3),CSM_DataTableC25(ll,1)*sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC25(mm,4),CSM_DataTableC25(ll,2)*sqrt(2*FLD.C/(PGm2(mm) + 1)),VerData_Rogers1992_TableC25(mm,5),CSM_DataTableC25(ll,3)*2/sqrt(2*FLD.C/(PGm2(mm)+1)),VerData_Rogers1992_TableC25(mm,6),CSM_DataTableC25(ll,4)*2/sqrt(2*FLD.C/(PGm2(mm)+1)),VerData_Rogers1992_TableC25(mm,7),CSM_DataTableC25(ll,6)*2/sqrt(2*FLD.C/(PGm2(mm)+1)));
    end
end
% fprintf(fid, '        %% Data\n');
fprintf(fid, '        \\bottomrule\n\n');
fprintf(fid, '    \\end{tabular}\n');
fprintf(fid, '    \\begin{tablenotes}\n');
fprintf(fid, '        \\item[*] The solution does not converge for $m_2 = \\infty$, and therefore $m_2 = 1000$ is taken as approach.\\\\\n');
fprintf(fid, '        \\item[$\\diamond$] Same value as in row above! Probably a wrong value is printed here considering the CSM result.\n');
fprintf(fid, '    \\end{tablenotes}\n');
fprintf(fid, '\\end{threeparttable}\n');
% fprintf(fid, '%%\\end{sidewaystable}\n');
% fprintf(fid, '%%\\end{table\n');
fclose(fid);
OldFolder = cd(NewFolder);

%% PLOT Results
