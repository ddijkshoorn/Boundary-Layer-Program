function [S,B,R,sol] = COEF(IT,NS,NP,GRD,HVR,sol,solprev)
% Code (FORTRAN) obtained from DVD enclosed with 'Convective Heat Transfer' by Cebeci (2002)
% Few changes made (see comments)
% NB consider preallocating of variables to increase speed

%% Initialization
% Initialization
Deta = GRD.Deta;
alpha0 = HVR.alpha0;
P1 = HVR.P1;
P2 = HVR.P2;
WW = HVR.WW;
CEL = HVR.CEL(NS);
P1P = HVR.P1P(NS);
P2P = HVR.P2P(NS);
S = zeros(NP,8); %%% S = zeros(NPT,8);
B = zeros(NP,10);%%% B = zeros(NPT,8);
R = zeros(NP,5); %%% R = zeros(NPT,8);

% % Calculation of coefficients %%% Move to PRECAL -> change into HVAR
if IT == 1
    if alpha0 > 0.1 % if OPT.BCEE == 0;
        sol.g(1) = WW(NS); % Input wall values EE
    else % ALPHA1 > 0 % else % OPT.BCEE == 1;
        sol.p(1) = WW(NS);
    end
end

for j = 2:NP
    % Present station
    USB = 0.5*(sol.u(j)^2 + sol.u(j-1)^2);
    FVB = 0.5*(sol.f(j)*sol.v(j) + sol.f(j-1)*sol.v(j-1));
    FPB = 0.5*(sol.f(j)*sol.p(j) + sol.f(j-1)*sol.p(j-1));
    UGB = 0.5*(sol.u(j)*sol.g(j) + sol.u(j-1)*sol.g(j-1));
    UB = 0.5*(sol.u(j) + sol.u(j-1));
    VB = 0.5*(sol.v(j) + sol.v(j-1));
    FB = 0.5*(sol.f(j) + sol.f(j-1));
    GB = 0.5*(sol.g(j) + sol.g(j-1));
    PB = 0.5*(sol.p(j) + sol.p(j-1));
    CB = 0.5*(sol.c(j) + sol.c(j-1));
    DERBV = (sol.b(j)*sol.v(j) - sol.b(j-1)*sol.v(j-1))/Deta(j-1);
    DEREP = (sol.e(j)*sol.p(j) - sol.e(j-1)*sol.p(j-1))/Deta(j-1);
    DRDUV = (sol.d(j)*sol.u(j)*sol.v(j) - sol.d(j-1)*sol.u(j-1)*sol.v(j-1))/Deta(j-1);
    
    % Previous station
    if NS == 1 % if first station, then...
        CFB    =   0;
        CVB    =   0;
        CPB    =   0;
        CUB    =   0;
        CGB    =   0;
%         CUGB   =   0;     % not used
%         CFPB   =   0;     % not used
%         CFVB   =   0;     % not used
%         CUSB   =   0;     % not used
%         CDERBV =   0;     % not used
%         CDEREP =   0;     % not used
        CRB    =  -P2(NS)*CB;
        CTB    =   0;
    else
        CFB    = 0.5*(solprev.f(j) + solprev.f(j-1));
        CVB    = 0.5*(solprev.v(j) + solprev.v(j-1));
        CPB    = 0.5*(solprev.p(j) + solprev.p(j-1));
        CUB    = 0.5*(solprev.u(j) + solprev.u(j-1));
        CGB    = 0.5*(solprev.g(j) + solprev.g(j-1));
        CUGB   = 0.5*(solprev.u(j)*solprev.g(j) + solprev.u(j-1)*solprev.g(j-1));
        CFPB   = 0.5*(solprev.f(j)*solprev.p(j) + solprev.f(j-1)*solprev.p(j-1));
        CFVB   = 0.5*(solprev.f(j)*solprev.v(j) + solprev.f(j-1)*solprev.v(j-1));
        CUSB   = 0.5*(solprev.u(j)^2 + solprev.u(j-1)^2);
        CCB    = 0.5*(solprev.c(j) + solprev.c(j-1));
        CDERBV = (solprev.b(j)*solprev.v(j) - solprev.b(j-1)*solprev.v(j-1))/Deta(j-1);
        CDEREP = (solprev.e(j)*solprev.p(j) - solprev.e(j-1)*solprev.p(j-1))/Deta(j-1);
        CDRDUV = (solprev.d(j)*solprev.u(j)*solprev.v(j) - solprev.d(j-1)*solprev.u(j-1)*solprev.v(j-1))/Deta(j-1);
        
        CLB    = CDERBV + P1(NS-1)*CFVB + P2(NS-1)*(CCB - CUSB);
        CRB    = -CLB - P2(NS)*CB - CEL*CUSB + CEL*CFVB;
        CMB    = CDEREP + CDRDUV + P1(NS-1)*CFPB;
        CTB    = -CMB + CEL*(CFPB - CUGB);
    end
    
    % Coefficients of differenced Momentum Equation (ME)
    % (sj(1) is not used)
    S(j,1) = sol.b(j)/Deta(j-1) + 0.5*P1P*sol.f(j) - 0.5*CEL*CFB;
    S(j,2) = -sol.b(j-1)/Deta(j-1) + 0.5*P1P*sol.f(j-1) - 0.5*CEL*CFB;
    S(j,3) = 0.5*(P1P*sol.v(j) + CEL*CVB);
    S(j,4) = 0.5*(P1P*sol.v(j-1) + CEL*CVB);
    S(j,5) = -P2P*sol.u(j);
    S(j,6) = -P2P*sol.u(j-1);
    S(j,7) = 0;
    S(j,8) = 0;
    
    R(j,2) = CRB - (DERBV + P1P*FVB - P2P*USB + CEL*(FB*CVB - VB*CFB));
    
    % Coefficients of differenced Energy Equation (EE) %%% last solution
    % (bj(1) is not used)
    B(j,1) = sol.e(j)/Deta(j-1) + 0.5*P1P*sol.f(j) - 0.5*CEL*CFB;
    B(j,2) = -sol.e(j-1)/Deta(j-1) + 0.5*P1P*sol.f(j-1) - 0.5*CEL*CFB;
    B(j,3) = 0.5*(P1P*sol.p(j) + CEL*CPB);
    B(j,4) = 0.5*(P1P*sol.p(j-1) + CEL*CPB);
    B(j,5) = sol.d(j)*sol.v(j)/Deta(j-1) - 0.5*CEL*(sol.g(j) - CGB);
    B(j,6) = -sol.d(j-1)*sol.v(j-1)/Deta(j-1) - 0.5*CEL*(sol.g(j-1) - CGB);
    B(j,7) = -0.5*CEL*(sol.u(j) + CUB);
    B(j,8) = -0.5*CEL*(sol.u(j-1) + CUB);
    B(j,9) = sol.d(j)*sol.u(j)/Deta(j-1);
    B(j,10) = -sol.d(j-1)*sol.u(j-1)/Deta(j-1);
    
    R(j,3) = CTB - (DEREP + DRDUV + P1P*FPB - CEL*(UGB - CGB*UB + CUB*GB) + CEL*(CPB*FB - CFB*PB));
    
    % Definitions of Rj %%% last solution
    R(j,1) = sol.f(j-1) - sol.f(j) + Deta(j-1)*UB;
    R(j-1,4) = sol.u(j-1) - sol.u(j) + Deta(j-1)*VB;
    R(j-1,5) = sol.g(j-1) - sol.g(j) + Deta(j-1)*PB;
end

R(1,1) = 0;
R(1,2) = 0;
R(1,3) = 0;
R(NP,4) = 0;
R(NP,5) = 0;
