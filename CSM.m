function [BLC,FLP,SOL,TCC,MON] = CSM(X,NP,NS,NST,GRD,HVR,EDG,FRS,FLD,TCC,sol,solprev,SOL,BLC,INP,OPT,SET,MON)
% CSM = Cebeci-Smith-method
% The CS-method is an approach to predict compressible boundary layer
% flows. The method is described in detail first in the book
% 'Analysis of Turbulent Boundary Layers' (1974) by Tuncer Cebeci and
% A.M.O. Smith.

% Changes to original (FORTRAN code paragraph 11.2 in Convective Heat Transfer by Cebeci, 2002):
% - added a missing coefficient to the D-coefficient (was missing in FORTRAN-code source); and,
% - updated D-coefficient after calling EDDY-file (was missing in FORTRAN-code source)
% - update Fluid Properties (call FLDPRS-function) after last iteration (right after last solver call)
% - added 2 transition checks, engineering estimates (in separate TRANS-file)
% - added relaminarization check; forced, and 1 engineering estimate (in separate RELAM-file)
% - added turbulent Prandtl-number model (in EDDY-file)

% Initialization
flp = [];
FLP = [];
MGE = 0; % count number of mesh grid extensions at every station
ITT = 0; % keep track of total number of iterations when mesh is extended

% Start of BL calculations
while NS <= NST

    DvW = 1;
    IT = 0;
    while MON.tr < 1 && abs(DvW) > 1e-5 || MON.tr > 0 && abs(DvW/(sol.v(1) + 0.5*DvW)) > 0.02
        IT = IT + 1;
        if IT > SET.ITMAX % 10 %%% Nonideal Gas Simulations can require more iterations (for example 8 for station 1 of acc. case D6)
            fprintf('Calculations stopped: number of iterations IT exceeded ITMAX = %g at station %g \n',SET.ITMAX, NS)
            return
        end

        % Fluid properties
        [sol,flp] = FLDPRS(NS,NP,EDG,sol,flp,FLD,OPT);

        if MON.tr > 0
            [flp] = EDDY(NS,NP,sol,GRD,EDG,HVR,TCC,flp,OPT);
            sol.b(1,1) = flp.C(1);
            % leave 'sol.d(1,1)' as is
            sol.e(1,1) = flp.C(1)/flp.Pr(1);
            for j = 2:NP
                sol.b(j,1) = (1.0 + flp.EV(j))*flp.C(j);
                sol.d(j,1) = flp.C(j)*(EDG.UE(NS)^2)*(1 - 1/flp.Pr(j) + flp.EV(j)*(1 - 1/flp.PrT(j)))/EDG.HtE(NS); % added, missing in FORTRAN-code from Cebeci (2002)
                sol.e(j,1) = flp.C(j)*((1.0 + flp.EV(j)*flp.Pr(j)/flp.PrT(j))/flp.Pr(j));
            end
        else
            for j = 1:NP
                sol.b(j,1) = flp.C(j);
                % leave 'sol.d(j,1)' as is
                sol.e(j,1) = flp.C(j)/flp.Pr(j);
            end
        end

        [S,B,R,sol] = COEF(IT,NS,NP,GRD,HVR,sol,solprev);
        [sol, DvW] = SOLV5(NP,S,B,R,GRD,HVR,sol);

        if sol.v(1) < 0
            fprintf('Separation occurred at station %g/%g \n',NS, NST)
            MON.SEP = 1;
            MON.STR = NS;
            % NB results of this station are not stored
            return
        end
    end
    if abs(sol.v(NP)) <= 1e-3 || NS == 1 % if grid too small at NX=1, then stop and change ETAE!
        [sol,flp] = FLDPRS(NS,NP,EDG,sol,flp,FLD,OPT); % update fluid properties; added to prevent spikes in wall temperature (temperature calculation in FLDPRS is based on previous velocity and enthalpy profiles)
%         [flp] = EDDY(NS,NP,sol,GRD,EDG,HVR,TCC,flp,OPT); % see NB line below
        % NB EV and PrT not updated here, since result would diverge from original FORTRAN code (at station 68 for Cebeci2002)
        
        % store variables
        [BLC,MON,FLP,SOL,MGE,ITT] = OUTPUT(X,NS,NP,MGE,IT,ITT,OPT,GRD,EDG,FRS,FLD,FLP,INP,HVR,SOL,BLC,MON,sol,flp);
        solprev = sol;
        
        if MON.tr < 1 && OPT.TRME > 0
            %         if (tr < 1 && OPT.TRME > 1 && NS~=1) % check if possible: NTR = 1,2,3,4,5,etc.
            [MON,TCC] = TRANS(NS,X,TCC,EDG,BLC,OPT,SET,MON);

            if MON.tr > 0 % happens only once
                if OPT.TRME == 1
                    fprintf('Transition occurred at station %g/%g \n',NS + 1, NST)
                    NS = NS + 1; % go to next station (turbulent flow)
                else
                    fprintf('Transition occurred at station %g/%g \n',NS, NST)
                    % Wazzan and Michel's methods, restart calculation at this
                    % station, but now for turbulent flow
                end
            else
                NS = NS + 1;
            end

        elseif OPT.RLAM > 0 && MON.tr > 0 && MON.rl < 1
            % Relaminarization check, for air only! (no experimental data available for other fluids)
            [MON,TCC] = RELAM(NS,X,HVR,TCC,SOL,EDG,FLP,BLC,OPT,SET,MON);
            if MON.rl > 0 % happens only once
                if OPT.RLAM == 1
                    fprintf('Relaminarization occurred at station %g/%g \n',NS + 1, NST)
                    NS = NS + 1;
                else
                    fprintf('Relaminarization occurred at station %g/%g \n',NS, NST)
                end
            else
                NS = NS + 1;
            end
        else
            NS = NS + 1;
        end
    else
        % Mesh growth
        ITT = IT; % restart IT after mesh growth
        MGE = MGE + 1;
        NPO = NP;
        NP1 = NP + 1;
        NP = NP + 1; % grid is extended with only one single point
        if NP > SET.NPT
            fprintf('Warning: NP exceeded NPT = %g at station = %g.\nGrid is extended and calculation continued.\nMaximum vector size was reached.\nPreallocate larger vector size for this case in INPUT-file.\n',SET.NPT,NS)
            NPTO = SET.NPT;
            SET.NPT = SET.NPT + 20;
            for j = NPTO + 1:SET.NPT
                GRD.Deta(j,1) = GRD.VGP*GRD.Deta(j-1);        % preallocate
                GRD.A(j,1) = 0.5*GRD.Deta(j-1);               % preallocate
                GRD.eta(j,1) = GRD.eta(j-1) + GRD.Deta(j-1);  % preallocate
            end
        end

        % Define new profiles
        for j = NP1:NP
            % current station
            sol.f(j) = sol.f(NPO) + GRD.eta(j) - GRD.eta(NPO);
            sol.u(j) = 1;
            sol.v(j) = 0;
            sol.g(j) = sol.g(NPO);
            sol.p(j) = 0;
            sol.b(j) = sol.b(NPO);
            sol.c(j) = sol.c(NPO);
            sol.d(j) = sol.d(NPO);
            sol.e(j) = sol.e(NPO);
            
            % previous station
            solprev.f(j) = solprev.f(NPO) + GRD.eta(j) - GRD.eta(NPO);
            solprev.u(j) = 1;
            solprev.v(j) = 0;
            solprev.g(j) = solprev.g(NPO);
            solprev.p(j) = 0;
            solprev.b(j) = solprev.b(NPO);
            solprev.c(j) = solprev.c(NPO);
            solprev.d(j) = solprev.d(NPO);
            solprev.e(j) = solprev.e(NPO);
        end
    end
end
