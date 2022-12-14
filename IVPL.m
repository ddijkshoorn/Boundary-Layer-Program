function [sol,solprev] = IVPL(NS,GRD,HVR,FLD,SET)
% IVPL calculates the initial profiles using Calorically Perfect Ideal Gas,
% and constant fluid properties for laminar zero PG incompressible adiabatic flow.
% Code (FORTRAN) obtained from DVD enclosed with 'Convective Heat Transfer' by Cebeci (2002)
% No changes made (see comments)

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
IT = 0;
DvW = 1;
while abs(DvW) > 1e-5
    IT = IT + 1;
    if IT > SET.ITMAX0 % set in case specific INPUT-file
        fprintf('Calculation of initial profiles did not converge.\n')
        keyboard
        return % return and try with current (not converged) solution
    end
    [S,B,R,sol] = COEF(IT,NS,NP,GRD,HVR,sol,solprev);
    [sol, DvW] = SOLV5(NP,S,B,R,GRD,HVR,sol);
end
