%% Verification with Howarths Flow - Table 6 from Smith and Clutter (1961) p.59
% Howarths flow: Nonsimilar decelarating flow
% adiabatic incompressible adverse pressure gradient flow with separation
% Nonsimilar solutions
% 29-11-2019 D.D.D.
% last changes on 12-12-2019
% Last run on 05-05-2020
% Simulation takes about 76 minutes
% Adapted slightly for upload and checked: 22-01-2022

%% Note: only for single station calculations

clear all
close all
clc
pause(0.1)

c1 = clock;

%% INPUT
load('../Verification_Study/Data/VerData_Smith1961_Table6_HowarthsFlow')
run('../Verification_Study/INPUT/INPUT_Verification_with_Smith1961_Howarth') % changed name, added Howarth (12-02-2020)
cd ..

SET.ITMAX0 = 20;    %%% added (21-01-2022)
OPT.RLAM = 0;       %%% added (21-01-2022)

% for the table (same as INP.x in INPUT-file):
x_print = ["0" "0.0125" "0.025" "0.050" "0.075" "0.100" "0.150" "0.200" "0.300" "0.400" "0.600" "0.800" "0.840" "0.880" "0.920" "0.948" "0.956" "0.958" "0.9589"];                    % [-]
deta = ["1.0"; "0.5"; "0.2"; "0.1"; "0.01"; "0.001"; "0.0001"]; % h1 = [1; 0.5; 0.2; 0.1; 0.01; 0.001; 0.0001; 0.00001];
dx = "0.0001"; % only used for length of fpp0, makes no sense actually % INP.x = [0 0.0125 0.025 0.050 0.075 0.100 0.150 0.200 0.300 0.400 0.600 0.800 0.840 0.880 0.920 0.948 0.956 0.958 0.9589]; % [-]
Size = length(deta); % is equal to length(dx)
SizeX = length(INP.x);
Xoriginal = INP.x;
fpp0 = NaN*ones((INP.x(end) - INP.x(1))/str2double(dx) + 1,Size);
MON_SEP = NaN*ones(3,Size); % store separation location for each succesive calculation
% POS1 = NaN*ones(SizeX,Size); % store location index of original x-coordinate for table 1
% POS2 = NaN*ones(SizeX,Size); % store location index of original x-coordinate for table 2
% POS3 = NaN*ones(SizeX,Size); % store location index of original x-coordinate for table 3
c = NaN*ones(SizeX,Size); % store location index of original x-coordinate for tables 1 to 3
c_1 = c;
c_2 = c;
c_3 = c;
npointsX = NaN*ones(3,Size);
npointsEta = NaN*ones(3,Size);

for ii = 1:3 %1:3%1:2 % three cases: 1) refine eta; 2) refine x; 3) refine both
    
    if ii == 2
        SET.NPT = 1000;
        GRD.Deta = 0.01;
        detaX = GRD.Deta;
    end
    
    for jj = 1:Size % length of refinements
% % %   % #1  if ii == 3 && jj == Size  % added to exclude ii=3,jj=7 (because of exceeding RAM)
% % %             % skip (out of RAM because of 1e6 X 1e6 matrix
% % %         else % added to exclude ii=3,jj=7 (because of exceeding RAM)
            % Grid refinement eta-coordinate
            if ii ~= 2
                GRD.Deta(1) = str2double(deta(jj));
                if jj == 7
                    SET.NPT = 1e5;%1e6;
                elseif jj == 6
                    SET.NPT = 1e5;
                else
                    SET.NPT = 1000;
                end
            end
            
            if ii > 1 % add this statement to selection of right indices, where needed?
                % Grid refinement x-coordinate
                if jj == 1
                    INP.x = Xoriginal;
                    c(1:SizeX,jj) = 1:SizeX;
                    % #1
% % %                 elseif jj == Size
% % %                     INP.x = 0:str2double(dx):Xoriginal(end);
% % %                     for mm = 1:SizeX
% % %                         c(mm,jj) = round(Xoriginal(mm)/str2double(dx)) + 1;
% % %                     end
                    %             c(1:SizeX,jj) = (Xoriginal(end) - 0)/str2double(dx) + 1;
                else % refine x-grid
                    for kk = 1:SizeX
                        if kk == 1
                            b = Xoriginal(1);
                            c(kk,jj) = 1; % coordinate (for table)
                            %                     else % OLD
                            %                         a = linspace(Xoriginal(kk-1),Xoriginal(kk),jj+1);
                            %                         b = [b a(2:end)];
                            %                         c(kk,jj) = c(kk-1,jj) + jj; % coordinate (for table)
                            %                     else % NEW
                            %                         a = linspace(Xoriginal(kk-1),Xoriginal(kk),2*nprev - 1);
                            %                         b = [b a(2:end)];
                            %                         c(kk,jj) = 2*c(kk-1,jj) - 1; % coordinate (for table)
                        else % NEW2
                            a = linspace(Xoriginal(kk-1),Xoriginal(kk),2^(jj - 1) + 1);
                            b = [b a(2:end)];
                            c(kk,jj) = c(kk-1,jj) + 2^(jj - 1); % coordinate (for table)
                        end
                    end
                    INP.x = b;
                end
                INP.y = zeros(1,length(INP.x)); % [-]
                INP.uE = (1 - 1/8*INP.x)./INP.UI;  % [m/s]
                INP.BCW = zeros(1,length(INP.x));   % [-], zero wall heat flux: derivative at wall of enthalpy ratio equal to zero
            else
                c(1:SizeX,jj) = 1:SizeX;
            end
            npointsX(ii,jj) = length(INP.x);
            
            % find right indices for plotting
            %         for kk = 1:SizeX
            %             if ii == 5
            %                 POS1(kk,jj) = find(INP.x==Xoriginal(kk));
            %             elseif ii == 2
            %                 POS2(kk,jj) = find(INP.x==Xoriginal(kk));
            %             elseif ii == 3
            %                 POS3(kk,jj) = find(INP.x==Xoriginal(kk));
            %             end
            %         end
            
            %% Calculation
            % Boundary layer edge calculations
            [X,NST,HVR,EDG,FRS,FLD,TCC,SOL,BLC,MON] = PRECAL(INP,OPT,SET,FLD,TCC);
            
            % Grid generation
            [NP,GRD] = GRID(GRD,SET);
            npointsEta(ii,jj) = NP;
            NS = 1;
            
            % IVPL
%             [sol,solprev] = IVPL(NS,GRD,HVR,FLD); %%% solprev = []; (empty, maybe needed for some case? move to PRECAL?)
            
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

            flp = [];
            
            tic % comment in future?
            
            [BLC,FLP,SOL,TCC,MON] = CSM(X,NP,NS,NST,GRD,HVR,EDG,FRS,FLD,TCC,sol,solprev,SOL,BLC,INP,OPT,SET,MON);
            
            for j = 1:length(INP.x)
                if isfield(MON,'STR') && j >= MON.STR % separation has occurred, skip station
                    fpp0(j,jj) = NaN;
                else
                    fpp0(j,jj) = SOL{j}.v(1);
                end
            end
            if isfield(MON,'STR')
                MON_SEP(ii,jj) = MON.STR;
            else
                MON_SEP(ii,jj) = NaN;
            end
            
            toc % comment in future?
            
            fprintf('Simulation finished with success\n')
            
            clearvars BLC FLP SOL MON sol solprev NP X NST HVR EDG FRS
            
%  % #1   end % added to exclude ii=3,jj=7 (because of exceeding RAM)
    end
    if ii == 1 % table 1
        HowarthData1 = fpp0;
        c_1 = c;
    elseif ii == 2 % table 2
        HowarthData2 = fpp0;
        c_2 = c;
    elseif ii == 3 % table 3
        HowarthData3 = fpp0;
        c_3 = c;
    end
    
    clearvars fpp0 c
    
    %     end
end

c2 = clock;

% keyboard

% % check where separation occurs
% for kk = 1:SizeX
%     for jj = 1:Size
%         %     if POS1(kk,jj) >= MON_SEP(1,jj)
%         if MON_SEP(1,jj) <= POS1(end-1,jj)
%             % warning
%             fprintf('separation occurs unexpected at an earlier station then last station of original x-coordinate, change table output accordingly for Table 1 refinement number %1.0f',jj)
%         end
%         if MON_SEP(2,jj) <= POS2(end-1,jj)
%             % warning
%             fprintf('separation occurs unexpected at an earlier station then last station of original x-coordinate, change table output accordingly for Table 2 refinement number %1.0f',jj)
%         end
%          if MON_SEP(3,jj) <= POS3(end-1,jj)
%             % warning
%             fprintf('separation occurs unexpected at an earlier station then last station of original x-coordinate, change table output accordingly for Table 3 refinement number %1.0f',jj)
%         end
%     end
% end

% this routine only checks if separation has occured AT the SECOND-LAST
% station
% check where separation occurs: check MON_SEP afterwards to be sure
for kk = 1:SizeX
    for jj = 1:Size
        %     if POS1(kk,jj) >= MON_SEP(1,jj)
        if MON_SEP(1,jj) <= c_1(end-1,jj)
            % warning
            fprintf('separation occurs unexpected at an earlier station then last station of original x-coordinate, change table output accordingly for Table 1 refinement number %1.0f\n',jj)
        end
        if MON_SEP(2,jj) <= c_2(end-1,jj)
            % warning
            fprintf('separation occurs unexpected at an earlier station then last station of original x-coordinate, change table output accordingly for Table 2 refinement number %1.0f\n',jj)
        end
        % uncommented #1 % commented because of omitting last jj
        if MON_SEP(3,jj) <= c_3(end-1,jj)
            % warning
            fprintf('separation occurs unexpected at an earlier station then last station of original x-coordinate, change table output accordingly for Table 3 refinement number %1.0f\n',jj)
        end
    end
end

% keyboard

%% TABLE - Write data to Tex-file (in Table format)

NewFolder = cd; % currently in 'boundary-layer-code'-folder
cd('Thesis_Report_Tables') %%% new, added (21-01-2022)
fid = fopen('Accuracy_study1_Howarth.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Accuracy Study with Howarths Flow case with Table 6 from Smith\\cite{Smith1961laminar}: incompressible adiabatic decelarating flow)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Accuracy study with nonsimilar Howarth''s Flow: incompressible adiabatic decelarating flow (adverse pressure gradient) following the example of Cebeci\\cite{cebeci1974analysis} table 8-2. The convergence criterion used was $\\left| \\mathrm{d} f''''(0) \\right| < 1e-5$. The values obtained with the CS-method (CSM) are in transformed compressible Falkner-Skan coordinates with grid height of $\\eta_{\\mathrm{e}} = %1.1f$.}\n',GRD.etaE);
fprintf(fid, '    \\label{tab:AccH1}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=1.4] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$X$}                       &\n');
fprintf(fid, '        \\multicolumn{7}{c}{$f''''(0)$}                  \\\\\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}                       &\n');
fprintf(fid, '        \\multicolumn{7}{c}{[-]}                       \\\\\n');
fprintf(fid, '        \\cmidrule(lr){2-8}\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(1));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(2));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(3));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(4));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(5));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(6));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}    \\\\\n',deta(7));
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:SizeX
    %     for mm = 1:Size % refinements
    if c_1(ll,7) >= MON_SEP(1,end)  % separation has occurred, print 'sep' (if separation occurs at last station of most refined grid, then probably it will occur at all other last stations; educated guess
        fprintf(fid, '        %s   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}   \\\\\n',num2str(x_print(ll)),HowarthData1(c_1(ll,1),1),HowarthData1(c_1(ll,2),2),HowarthData1(c_1(ll,3),3),HowarthData1(c_1(ll,4),4));
    else
        fprintf(fid, '        %s   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   \\\\\n',num2str(x_print(ll)),HowarthData1(c_1(ll,1),1),HowarthData1(c_1(ll,2),2),HowarthData1(c_1(ll,3),3),HowarthData1(c_1(ll,4),4),HowarthData1(c_1(ll,5),5),HowarthData1(c_1(ll,6),6),HowarthData1(c_1(ll,7),7));
    end
    %     end
end
%             if kk < 9 || kk == 10 || kk == 18
%             fprintf(fid, '        %s   &   %s   &   %s   &       &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3));
%         else
%             fprintf(fid, '        %s   &   %s   &   %s   &  %s   &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),VerData_Smith1961_Table6_HowarthsFlow(kk,4));
%         end
%         elseif kk < 9 || kk == 10 || kk == 18
%         if ll >= MON.STR % separation has occurred, print 'sep'
%             fprintf(fid, '        %s   &   %s   &   %s   &       &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3));
%         else
%             fprintf(fid, '        %s   &   %s   &   %s   &       &   %1.6f  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),SOL{ll}.v(1));
%         end
%     else
%         fprintf(fid, '        %s   &   %s   &   %s   &  %s   &  %1.6f   \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),VerData_Smith1961_Table6_HowarthsFlow(kk,4),SOL{ll}.v(1));
%     end
% end
fprintf(fid, '        \\addlinespace\n');
fprintf(fid, '        \\multicolumn{1}{c}{$n_{\\mathrm{\\eta}_{\\mathrm{points}}}$}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   \\\\\n',num2str(npointsEta(1,1)),num2str(npointsEta(1,2)),num2str(npointsEta(1,3)),num2str(npointsEta(1,4)),num2str(npointsEta(1,5)),num2str(npointsEta(1,6)),num2str(npointsEta(1,7)));
fprintf(fid, '        \\multicolumn{1}{c}{$n_{\\mathrm{X}_{\\mathrm{stations}}}$}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   \\\\\n',num2str(npointsX(1,1)),num2str(npointsX(1,2)),num2str(npointsX(1,3)),num2str(npointsX(1,4)),num2str(npointsX(1,5)),num2str(npointsX(1,6)),num2str(npointsX(1,7)));
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

% keyboard

%%

fid = fopen('Accuracy_study2_Howarth.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Accuracy Study with Howarths Flow case with Table 6 from Smith\\cite{Smith1961laminar}: incompressible adiabatic decelarating flow)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Accuracy study with nonsimilar Howarth''s Flow: incompressible adiabatic decelarating flow (adverse pressure gradient) following the example of Cebeci\\cite{cebeci1974analysis} table 8-3. The convergence criterion used was $\\left| \\mathrm{d} f''''(0) \\right| < 1e-5$. The values obtained with the CS-method (CSM) are in transformed compressible Falkner-Skan coordinates with grid height of $\\eta_{\\mathrm{e}} = %1.1f$.}\n',GRD.etaE);
fprintf(fid, '    \\label{tab:AccH2}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=1.4] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$X$}                       &\n');
fprintf(fid, '        \\multicolumn{7}{c}{$f''''(0)$}                  \\\\\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}                       &\n');
fprintf(fid, '        \\multicolumn{7}{c}{[-]}                       \\\\\n');
fprintf(fid, '        \\cmidrule(lr){2-8}\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',num2str(detaX));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',num2str(detaX));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',num2str(detaX));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',num2str(detaX));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',num2str(detaX));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',num2str(detaX));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}    \\\\\n\n',num2str(detaX));
fprintf(fid, '        \\multicolumn{1}{c}{}                         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X$}                &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /2 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /4 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /8 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /16 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /32 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /64 $}           \\\\\n');%X = %s $}    \\\\\n',dx);
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:SizeX
    %     for mm = 1:Size % refinements
    if c_2(ll,7) >= MON_SEP(2,end)  % separation has occurred, print 'sep' (if separation occurs at last station of most refined grid, then probably it will occur at all other last stations; educated guess
        fprintf(fid, '        %s   &    \\multicolumn{1}{c}{sep}    &    \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}   \\\\\n',num2str(x_print(ll)));
    else
        fprintf(fid, '        %s   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   \\\\\n',num2str(x_print(ll)),HowarthData2(c_2(ll,1),1),HowarthData2(c_2(ll,2),2),HowarthData2(c_2(ll,3),3),HowarthData2(c_2(ll,4),4),HowarthData2(c_2(ll,5),5),HowarthData2(c_2(ll,6),6),HowarthData2(c_2(ll,7),7));
    end
    %     end
end
%             if kk < 9 || kk == 10 || kk == 18
%             fprintf(fid, '        %s   &   %s   &   %s   &       &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3));
%         else
%             fprintf(fid, '        %s   &   %s   &   %s   &  %s   &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),VerData_Smith1961_Table6_HowarthsFlow(kk,4));
%         end
%         elseif kk < 9 || kk == 10 || kk == 18
%         if ll >= MON.STR % separation has occurred, print 'sep'
%             fprintf(fid, '        %s   &   %s   &   %s   &       &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3));
%         else
%             fprintf(fid, '        %s   &   %s   &   %s   &       &   %1.6f  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),SOL{ll}.v(1));
%         end
%     else
%         fprintf(fid, '        %s   &   %s   &   %s   &  %s   &  %1.6f   \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),VerData_Smith1961_Table6_HowarthsFlow(kk,4),SOL{ll}.v(1));
%     end
% end
fprintf(fid, '        \\addlinespace\n');
fprintf(fid, '        \\multicolumn{1}{c}{$n_{\\mathrm{\\eta}_{\\mathrm{points}}}$}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   \\\\\n',num2str(npointsEta(2,1)),num2str(npointsEta(2,2)),num2str(npointsEta(2,3)),num2str(npointsEta(2,4)),num2str(npointsEta(2,5)),num2str(npointsEta(2,6)),num2str(npointsEta(2,7)));
fprintf(fid, '        \\multicolumn{1}{c}{$n_{\\mathrm{X}_{\\mathrm{stations}}}$}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   \\\\\n',num2str(npointsX(2,1)),num2str(npointsX(2,2)),num2str(npointsX(2,3)),num2str(npointsX(2,4)),num2str(npointsX(2,5)),num2str(npointsX(2,6)),num2str(npointsX(2,7)));
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

% keyboard

%%

fid = fopen('Accuracy_study3_Howarth.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Accuracy Study with Howarths Flow case with Table 6 from Smith\\cite{Smith1961laminar}: incompressible adiabatic decelarating flow)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Accuracy study with nonsimilar Howarth''s Flow: incompressible adiabatic decelarating flow (adverse pressure gradient) following the example of Cebeci\\cite{cebeci1974analysis} table 8-4. The convergence criterion used was $\\left| \\mathrm{d} f''''(0) \\right| < 1e-5$. The values obtained with the CS-method (CSM) are in transformed compressible Falkner-Skan coordinates with grid height of $\\eta_{\\mathrm{e}} = %1.1f$.}\n',GRD.etaE);
fprintf(fid, '    \\label{tab:AccH3}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=1.4] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$X$}                       &\n');
fprintf(fid, '        \\multicolumn{7}{c}{$f''''(0)$}                  \\\\\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}                       &\n');
fprintf(fid, '        \\multicolumn{7}{c}{[-]}                       \\\\\n');
fprintf(fid, '        \\cmidrule(lr){2-8}\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(1));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(2));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(3));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(4));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(5));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(6));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}    \\\\\n\n',deta(7));
fprintf(fid, '        \\multicolumn{1}{c}{}                         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X$}                &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /2 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /4 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /8 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /16 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /32 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /64 $}           \\\\\n');%X = %s $}    \\\\\n',dx);
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:SizeX
    %     for mm = 1:Size % refinements
%     if c_3(ll,7) >= MON_SEP(3,end)  % separation has occurred, print 'sep' (if separation occurs at last station of most refined grid, then probably it will occur at all other last stations; educated guess
%         fprintf(fid, '        %s   &    \\multicolumn{1}{c}{sep}    &    \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}   \\\\\n',num2str(x_print(ll)));
    if ll == SizeX  % separation has occurred, print 'sep' (if separation occurs at last station of most refined grid, then probably it will occur at all other last stations; educated guess
        fprintf(fid, '        %s   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    \\\\\n',num2str(x_print(ll)),HowarthData3(c_3(ll,1),1),HowarthData3(c_3(ll,2),2),HowarthData3(c_3(ll,3),3),HowarthData3(c_3(ll,4)));
    else
        fprintf(fid, '        %s   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   \\\\\n',num2str(x_print(ll)),HowarthData3(c_3(ll,1),1),HowarthData3(c_3(ll,2),2),HowarthData3(c_3(ll,3),3),HowarthData3(c_3(ll,4),4),HowarthData3(c_3(ll,5),5),HowarthData3(c_3(ll,6),6),HowarthData3(c_3(ll,7),7));
    end
    %     end
end
%             if kk < 9 || kk == 10 || kk == 18
%             fprintf(fid, '        %s   &   %s   &   %s   &       &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3));
%         else
%             fprintf(fid, '        %s   &   %s   &   %s   &  %s   &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),VerData_Smith1961_Table6_HowarthsFlow(kk,4));
%         end
%         elseif kk < 9 || kk == 10 || kk == 18
%         if ll >= MON.STR % separation has occurred, print 'sep'
%             fprintf(fid, '        %s   &   %s   &   %s   &       &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3));
%         else
%             fprintf(fid, '        %s   &   %s   &   %s   &       &   %1.6f  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),SOL{ll}.v(1));
%         end
%     else
%         fprintf(fid, '        %s   &   %s   &   %s   &  %s   &  %1.6f   \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),VerData_Smith1961_Table6_HowarthsFlow(kk,4),SOL{ll}.v(1));
%     end
% end
fprintf(fid, '        \\addlinespace\n');
fprintf(fid, '        \\multicolumn{1}{c}{$n_{\\mathrm{\\eta}_{\\mathrm{points}}}$}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   \\\\\n',num2str(npointsEta(3,1)),num2str(npointsEta(3,2)),num2str(npointsEta(3,3)),num2str(npointsEta(3,4)),num2str(npointsEta(3,5)),num2str(npointsEta(3,6)),num2str(npointsEta(3,7)));
fprintf(fid, '        \\multicolumn{1}{c}{$n_{\\mathrm{X}_{\\mathrm{stations}}}$}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   \\\\\n',num2str(npointsX(3,1)),num2str(npointsX(3,2)),num2str(npointsX(3,3)),num2str(npointsX(3,4)),num2str(npointsX(3,5)),num2str(npointsX(3,6)),num2str(npointsX(3,7)));
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
% OldFolder = cd(NewFolder);

% keyboard

%%

fid = fopen('Accuracy_study3_Howarth_A.tex','w');
fprintf(fid, '%%%%%%%%%% TABLE PRINTED BY MATLAB PROGRAM %%%%%%%%%%\n');
fprintf(fid, '%%%% Accuracy Study with Howarths Flow case with Table 6 from Smith\\cite{Smith1961laminar}: incompressible adiabatic decelarating flow)\n');
% fprintf(fid, '%%\\begin{table}\n');
% fprintf(fid, '%%\\begin{sidewaystable}\n');
fprintf(fid, '\\begin{threeparttable}\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption{Accuracy study with nonsimilar Howarth''s Flow: incompressible adiabatic decelarating flow (adverse pressure gradient) following the example of Cebeci\\cite{cebeci1974analysis} table 8-4. The convergence criterion used was $\\left| \\mathrm{d} f''''(0) \\right| < 1e-5$. The values obtained with the CS-method (CSM) are in transformed compressible Falkner-Skan coordinates with grid height of $\\eta_{\\mathrm{e}} = %1.1f$.}\n',GRD.etaE);
fprintf(fid, '    \\label{tab:AccH3_A}\n');
fprintf(fid, '    \\begin{tabular}{S[table-format=1.4] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6] S[table-format=1.6]}\n');
fprintf(fid, '        \\toprule\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{$X$}                       &\n');
fprintf(fid, '        \\multicolumn{7}{c}{$f''''(0)$}                  \\\\\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{[-]}                       &\n');
fprintf(fid, '        \\multicolumn{7}{c}{[-]}                       \\\\\n');
fprintf(fid, '        \\cmidrule(lr){2-8}\n\n');
fprintf(fid, '        \\multicolumn{1}{c}{}                         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(1));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(2));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(3));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(4));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(5));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}                         &\n',deta(6));
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} \\eta = %s $}    \\\\\n\n',deta(7));
fprintf(fid, '        \\multicolumn{1}{c}{}                         &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X$}                &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /2 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /4 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /8 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /16 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /32 $}           &\n');
fprintf(fid, '        \\multicolumn{1}{c}{$\\mathrm{d} X /64 $}           \\\\\n');%X = %s $}    \\\\\n',dx);
fprintf(fid, '        \\midrule\n\n');
fprintf(fid, '        %% Data\n');
for ll = 1:SizeX
    %     for mm = 1:Size % refinements
%     if c_3(ll,7) >= MON_SEP(3,end)  % separation has occurred, print 'sep' (if separation occurs at last station of most refined grid, then probably it will occur at all other last stations; educated guess
%         fprintf(fid, '        %s   &    \\multicolumn{1}{c}{sep}    &    \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}   \\\\\n',num2str(x_print(ll)));
    if ll == SizeX  % separation has occurred, print 'sep' (if separation occurs at last station of most refined grid, then probably it will occur at all other last stations; educated guess
        fprintf(fid, '        %s   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    &   \\multicolumn{1}{c}{sep}    \\\\\n',num2str(x_print(ll)),HowarthData3(c_3(ll,1),1),HowarthData3(c_3(ll,2),2),HowarthData3(c_3(ll,3),3),HowarthData3(c_3(ll,4)));
    else
        fprintf(fid, '        %s   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   &  %1.6f   \\\\\n',num2str(x_print(ll)),HowarthData3(c_3(ll,1),1),HowarthData3(c_3(ll,2),2),HowarthData3(c_3(ll,3),3),HowarthData3(c_3(ll,4),4),HowarthData3(c_3(ll,5),5),HowarthData3(c_3(ll,6),6),HowarthData3(c_3(ll,7),7));
    end
    %     end
end
%             if kk < 9 || kk == 10 || kk == 18
%             fprintf(fid, '        %s   &   %s   &   %s   &       &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3));
%         else
%             fprintf(fid, '        %s   &   %s   &   %s   &  %s   &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),VerData_Smith1961_Table6_HowarthsFlow(kk,4));
%         end
%         elseif kk < 9 || kk == 10 || kk == 18
%         if ll >= MON.STR % separation has occurred, print 'sep'
%             fprintf(fid, '        %s   &   %s   &   %s   &       &  \\multicolumn{1}{c}{sep}  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3));
%         else
%             fprintf(fid, '        %s   &   %s   &   %s   &       &   %1.6f  \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),SOL{ll}.v(1));
%         end
%     else
%         fprintf(fid, '        %s   &   %s   &   %s   &  %s   &  %1.6f   \\\\\n',VerData_Smith1961_Table6_HowarthsFlow(kk,1),VerData_Smith1961_Table6_HowarthsFlow(kk,2),VerData_Smith1961_Table6_HowarthsFlow(kk,3),VerData_Smith1961_Table6_HowarthsFlow(kk,4),SOL{ll}.v(1));
%     end
% end
fprintf(fid, '        \\addlinespace\n');
fprintf(fid, '        \\multicolumn{1}{c}{$n_{\\mathrm{\\eta}_{\\mathrm{points}}}$}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   \\\\\n',num2str(npointsEta(3,1)),num2str(npointsEta(3,2)),num2str(npointsEta(3,3)),num2str(npointsEta(3,4)),num2str(npointsEta(3,5)),num2str(npointsEta(3,6)),num2str(npointsEta(3,7)));
fprintf(fid, '        \\multicolumn{1}{c}{$n_{\\mathrm{X}_{\\mathrm{stations}}}$}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   &   \\multicolumn{1}{c}{%s}   \\\\\n',num2str(npointsX(3,1)),num2str(npointsX(3,2)),num2str(npointsX(3,3)),num2str(npointsX(3,4)),num2str(npointsX(3,5)),num2str(npointsX(3,6)),num2str(npointsX(3,7)));
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
% PLOTFILE
