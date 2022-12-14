function [flp] = EDDY(NS,NP,sol,GRD,EDG,HVR,TCC,flp,OPT) % OPT only for constant PrT
% EV-model based on Cebeci (2002), PrT model based on Cebeci (1974)

%% Notes
% gamma_int currently not stored

IED = 0; % start with inner layer

% Initial values
PPLUS = HVR.P2(NS)/(EDG.Re_x(NS)^0.25)*(flp.C(1)*sol.c(1)*sol.v(1))^-1.5;
CN = sqrt(1 - TCC.ints*PPLUS*flp.C(1)*sol.c(1)^3);
CYA = CN*(EDG.Re_x(NS)^0.25)*sqrt(flp.C(1)*sol.v(1))/TCC.A_plus;  
BPLUS2 = 0;
for i = 1:length(TCC.Bcoeff)
    BPLUS2 = BPLUS2 + TCC.Bcoeff(i)*log10(flp.Pr(1))^(i-1);
end
BPLUS = BPLUS2/flp.Pr(1)^0.5;
if OPT.CPRT > 0
    flp.PrT(1,1) = 0.4/TCC.kappa_h*BPLUS/TCC.A_plus;
else
    flp.PrT(1,1) = TCC.PrT; %%% constant Turbulent Pr-number for original case
end

TERMP = sol.c(1);
SUM = 0;
CII = 0;
CI = 0;

for j = 2:NP
    CII = CII + GRD.A(j)*(sol.c(j-1) + sol.c(j));
    TERM = sol.c(j)*(1 - sol.u(j));
    SUM = SUM + GRD.A(j)*(TERM + TERMP);
    TERMP = TERM;
end

for j = 2:NP
    CI = CI + GRD.A(j)*(sol.c(j-1) + sol.c(j));
    EDVO = TCC.alpha*SUM*(EDG.Re_x(NS)^0.5)/flp.C(j)/sol.c(j)^2;
    gamma_int = 1.0/(1.0 + 5.5*(CI/CII)^6); %%% Note! Written in this form it depends on the eta grid!
    EDVO = EDVO*TCC.gamma_tr(NS)*gamma_int;
    %% Added - moved from below (for calculation of PrT)
    YOA = CYA*CI/flp.C(j)/sol.c(j)^1.5; %% for use in PrT
    EL = 1;
    if YOA <= 10
        EL = 1 - exp(-YOA);
    end
    %% Added
    BPLUS2 = 0;
    for i = 1:length(TCC.Bcoeff)
        BPLUS2 = BPLUS2 + TCC.Bcoeff(i)*log10(flp.Pr(j))^(i-1);
    end
    BPLUS = BPLUS2/flp.Pr(j)^0.5;
    CYB = CN*(EDG.Re_x(NS)^0.25)*sqrt(flp.C(1)*sol.v(1))/BPLUS; 
    YOB = CYB*CI/flp.C(j)/sol.c(j)^1.5;
    ELh = 1;
    if YOB <= 10
        ELh = 1 - exp(-YOB);
    end
    if OPT.CPRT > 0
        flp.PrT(j,1) = 0.4*EL/(TCC.kappa_h*ELh);
    else
        flp.PrT(j,1) = TCC.PrT; %%% for ogiringal case
    end
    if IED == 1
        flp.EV(j,1) = EDVO;
    else
        EDVI = TCC.kappa^2*sqrt(EDG.Re_x(NS))*TCC.gamma_tr(NS)*sol.v(j)*(CI*EL)^2/flp.C(j)/sol.c(j)^3;
        if EDVI < EDVO
            flp.EV(j,1) = EDVI;
        else
            IED = 1; % MON.IO = j; % inner-outer boundary inside BL
            flp.EV(j,1) = EDVO;
        end
    end
%     TCC.gamma_int(NP,NS) = gamma_int;

end