function [MON,TCC] = RELAM(NS,X,HVR,TCC,SOL,EDG,FLP,BLC,OPT,SET,MON)
% Relaminarization check, analogous to transition check for lam to turb
% flow
% NB turbulence model is still being used, since program uses 1 - gamma_tr 
% from Cebeci 2002 for smooth transitional flow from turbulent to laminar

% OPT.RLAM = 0,1,2;     % As set in INPUT-file
% Relaminarization method, 0=no relam.; 1=prescribed location NRL; 2=simple engineering estimate based on exp. data for air only!

% OPT.RLAM == 0 % no Relaminarization check
if OPT.RLAM == 1 % 1) Prescribed point of transition; looking one station ahead
    if NS >= SET.NRL - 1
        MON.rl = 1;
        MON.NRL = NS + 1;
        MON.Re_rl = EDG.Re_x(NS + 1);
%     else % superfluous, but otherwise warning: 'output arguments not assigned during call'
%         MON.rl = 0;
    end
elseif OPT.RLAM == 2 % 2) check if relaminarization has occured with method of Nash-Webber & Oates 1972

    Re_d2 = 1/FLP{NS}.C(1)/SOL{NS}.c(1)*BLC.Re_theta(NS);

    if Re_d2 > 4e3
        fprintf('Value of Reynolds momentum thickness out of range for relaminarization check at station %g: \n',NS);
    end

    % K = X(NS)^2/EDG.Re_x(NS)*FLP{NS}.C(1)*SOL{NS}.c(1)^2*HVR.P2(NS);
    K = FLP{NS}.C(1)*SOL{NS}.c(1)^2/EDG.Re_x(NS)*HVR.P2(NS);
    K_relam = 1.2e-6 + 1.1e-10*Re_d2 + 1e-13*Re_d2^2;

    if K > K_relam
        %     RELAM = 1;
        MON.rl = 1; % (uncommented on 11-05-2020)
        MON.NRL = NS;
        MON.Re_rl = EDG.Re_x(NS);
        %     MON.tr = 0; % ? wrong! since we first need to calculate EV up to gamma_tr(NS) ~ 0!
        %     fprintf('Relaminarization might occur at station %g: \n',NS);
        % else
        %     MON.rl = 0;
        %     % RELAM = 0;
    end

end

%% Reverse gamma-transition
% NB the following does not work if turbulence model is not used anymore

if MON.rl > 0 % if MON.rl > 0 was this on (11-05-2020) %%% MON.tr < 1 %%% (23-02-2020) what if we change this to MON.rl > 0? No, and what to change inside CSm-file?
    if OPT.RLAM == 1 % uses prescribed point of transition in a smart way
        NS = NS + 1; %%% Note! Is not an output of this function, therefore value outside remains equal to current NS.
    end %%% This is done (09-09-2019) to take into account prescribed transition for next station
    %     G_tr = 8.33e-4*UE(NTR-1)^3/(CNUE(NTR-1)^2*(RXNTR^1.34)); %%% added (NTR-1)  to CNUE (since it used to be: CNUE which is CNUE(NTR), why?)
%     G_tr = 8.33e-4*EDG.UE(NS - 1)/X(NS - 1)^2*EDG.Re_x(NS - 1)^0.66; %%% CHANGED WRT ORIGINAL!!! NEW: substituted fluid properties
    G_tr = 8.33e-4*EDG.UE(NS - 1)^3*EDG.rhoE(NS - 1)^2/EDG.muE(NS - 1)^2*EDG.Re_x(NS - 1)^-1.34;
%     fprintf('G_tr0 - G_tr = %g\n', G_tr0 - G_tr)
    for i = NS:length(X) % NS = point of transition station
        UE_G = 0;
        U1 = 1/EDG.UE(NS - 1);
        for j = NS:i
            U2 = 1/EDG.UE(j);
            UE_G = UE_G + 0.5*(U1 + U2)*(X(j) - (X(j - 1)));
            U1 = U2;
        end
        dX_G = X(i) - X(NS - 1);
        EXPT = G_tr*dX_G*UE_G;
        if EXPT <= 10
            gamma_tr = 1 - exp(-EXPT);
        else
            gamma_tr = 1;
        end
        TCC.gamma_tr(i) = 1 - gamma_tr;
%         if i == length(X)
%             keyboard
%         end
    end
%     MON.rl = 1;
end

%% Output of function needed when transition has not occurred
MON.rl_last_check = NS; % last station checked is NS
% TCC.rl = MON.rl;
