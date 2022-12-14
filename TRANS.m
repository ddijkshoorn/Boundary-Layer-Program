function [MON,TCC] = TRANS(NS,X,TCC,EDG,BLC,OPT,SET,MON)
% Prediction of transition from laminar to turbulent flow

% OPT.TRME = 0,1,2,3,4;     % As set in INPUT-file
% Transition method: 0=no transition; 1=prescribed transition location NTR; 2=Wazzan correlation; 3=Michel's Method (Cebeci (2002)); 4=fully turbulent (for air only!)

% The following three methods are used:

% 1) Prescribed point of transition (forced transition)

% 2) Method of Wazzan (1981): 
%    a method for a rough estimate of transition location derived from the e^9 method 

% 3) Michel's Method (see Cebeci (2002)):
%    See Cebeci (2002) p175/176

if OPT.TRME == 1 % 1) Prescribed point of transition; looking one station ahead
    if NS >= SET.NTR - 1
        MON.tr = 1;
        MON.NTR = NS + 1;
        MON.Re_tr = EDG.Re_x(NS + 1); %%% this is taken at the first turbulent station, Re_x for G_tr for previous station
%     else % superfluous, but otherwise warning: 'output arguments not assigned during call'
%         MON.tr = 0;
    end
    
elseif OPT.TRME == 2 % 2) Method of Wazzan
    if BLC.H(NS) >= 2.1 && BLC.H(NS) <= 2.9
        logRe_tr = -40.4557 + 64.8066*BLC.H(NS) - 26.7538*BLC.H(NS)^2 + 3.3819*BLC.H(NS)^3;
    elseif BLC.H(NS) > 2.9
        logRe_tr = 5;
        fprintf('Value of form factor H is outside validity range of Wazzan''s method (2.1-2.8) at station %g: \n H = %g \n ', NS, BLC.H(NS));
    else % < 2.1
        logRe_tr = 9.5;
        fprintf('Value of form factor H is outside validity range of Wazzan''s method (2.1-2.8) at station %g: \n H = %g \n ', NS, BLC.H(NS));
    end

    Re_tr = 10^logRe_tr;

    if Re_tr < BLC.Re_x(NS)
        MON.tr = 1;
        MON.NTR = NS;
        MON.Re_tr = Re_tr; %%% this is taken at the first turbulent station, Re_x for G_tr for previous station
        fprintf('Transition to turbulent flow occurred at station %g, with Re_x = %g while Re_x_tr = %g \n', NS, BLC.Re_x(NS), Re_tr);
%     else % superfluous, but otherwise warning: 'output arguments not assigned during call'
%         MON.tr = 0;
    end
    
elseif OPT.TRME == 3 % Michel's Method (Cebeci (2002) p.175)
    if BLC.Re_theta(NS) > 1.174*(1 + 22400/BLC.Re_x(NS))*BLC.Re_x(NS)^0.46
        MON.tr = 1;
        MON.NTR = NS;
        MON.Re_tr = BLC.Re_x(NS); %%% this is taken at the first turbulent station, Re_x for G_tr for previous station
        fprintf('Transition to turbulent flow occurred at station %g, with Re_theta = %g and Re_x = %g \n', NS, BLC.Re_theta(NS), BLC.Re_x(NS));
%     else % superfluous, but otherwise warning: 'output arguments not assigned during call'
%         MON.tr = 0;
    end
end

%% Calculation of gamma_transition for (numerical) smooth transition
% the following relation for modeling the transition flow region is obtained from the CS-method from Cebeci (2002).
% Note that it is only valid for adiabatic flows with Ma < 5 (see Cebeci (2004) p274)

if MON.tr > 0 % only once
    if OPT.TRME == 1 % uses prescribed point of transition in a smart way
        NS = NS + 1; %%% Note! Is not an output of this function, therefore value outside remains equal to current NS.
    end %%% This is done (09-09-2019) to take into account prescribed transition for next station
    %     G_tr = 8.33e-4*UE(NTR-1)^3/(CNUE(NTR-1)^2*(RXNTR^1.34)); %%% added (NTR-1)
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
            TCC.gamma_tr(i) = 1 - exp(-EXPT);
        else
            TCC.gamma_tr(i) = 1;
        end
    end
end

%% Output needed when transition has not occurred
MON.tr_last_check = NS; % last station checked is NS
% TCC.tr = MON.tr;
