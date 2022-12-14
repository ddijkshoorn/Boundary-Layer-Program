function [BLC,MON,FLP,SOL,MGE,ITT] = OUTPUT(X,NS,NP,MGE,IT,ITT,OPT,GRD,EDG,FRS,FLD,FLP,INP,HVR,SOL,BLC,MON,sol,flp)
% OUTPUT calculates Boundary Layer Characteristics (BLC),
% and dimensionless numbers
% Calculation of integral boundary layer properties based on Cebeci (2002)
% for comparison

%% Heat transfer

if OPT.BCEE == 0 % || abs(sol.p(1)) < 1e-6 % adiabatic wall; in case OPT.BCEE = 4 and INP.AWD = 0; or close to zero 
    BLC.g_aw(NS) = sol.g(1);
    BLC.H_aw(NS) = EDG.HtE(NS)*sol.g(1);
    BLC.T_aw(NS) = flp.T(1);
    % only for adiabatic case, and only stationary frame (stator, not valid for rotor: since work is extracted), in other words: total temperature needs to remain constant
    BLC.Trecovery(NS) = (flp.T(1) - EDG.TsE(NS))/(EDG.TtE(NS) - EDG.TsE(NS));
    BLC.Hrecovery(NS) = (EDG.HtE(NS)*sol.g(1) - EDG.HsE(NS))/(EDG.HtE(NS) - EDG.HsE(NS));
elseif ~isempty(INP.AWD)
    % adiabatic wall enthalpy ratio
    if OPT.BCEE == 1 || OPT.BCEE == 2 % g_aw as input
        g_aw = INP.AWD(NS);
        H_aw = EDG.HtE(NS)*INP.AWD(NS);
        T_aw = H_aw/flp.Cp(1);
    else % OPT.BCEE == 3 || OPT.BCEE == 4 % T_aw [K] given as input (adiabatic wall temperature)
        T_aw = INP.AWD(NS);
        if OPT.GASM == 3
            % approximation: Cp = Cp(T)
%             H_aw = flp.Cp(1)*T_aw;
            H_aw = 1e3*FLD.FP.Enthalpy('PT',EDG.PsE(NS)/1e5,T_aw - 273.15); % [J/kg]
            g_aw = H_aw/EDG.HtE(NS);
        elseif OPT.GASM == 2 % OPT2
            H_aw = FLD.Cpcoeff(1)*T_aw + FLD.Cpcoeff(2)/2*T_aw^2 + FLD.Cpcoeff(3)/3*T_aw^3 + FLD.Cpcoeff(4)/4*T_aw^4 + FLD.Cpcoeff(5)/5*T_aw^5; % [J/kg], Andrews (1981) IG relation; range: 100 - 590 K
            g_aw = H_aw/EDG.HtE(NS);
        else % OPT1
            H_aw = FLD.Cp*T_aw; % exact in this case, since Cp = Cp(T) or Cp = constant
            g_aw = H_aw/EDG.HtE(NS);
        end
    end
else % estimate adiabatic wall enthalpy/temperature in case INP.AWD has been left empty
    if OPT.COMP > 0
        if MON.tr > 0 % NS >= NTR (transition location dependency prevents this to be moved to PRECAL-file)
            r = flp.Pr(NP)^(1/3); % turbulent - for AIR! (Pr < 1)
        else
            r = flp.Pr(NP)^(1/2); % laminar - for AIR! (Pr < 1)
        end
    else % incompressible approximation
        r = 1;
    end
    % estimate recovery temperature/enthalpy (adiabatic wall)
%     Tr = EDG.TsE(NS)*(1 + r*(FLD.gamma - 1)/2*EDG.MaE(NS)^2); % compressible
    Tr = EDG.TsE(NS)*(1 + r*(EDG.gammaE(NS) - 1)/2*EDG.MaE(NS)^2); % compressible
    Hr = EDG.HsE(NS) + r*EDG.UE(NS)^2/2; % Hr = EDG.HsE(NS)*(1 + r*(flp.gamma(NP) - 1)/2*EDG.MaE(NS)^2);
    T_aw = Tr;
    H_aw = Hr;
    g_aw = H_aw/EDG.HtE(NS);
end

% Enthalpy driving force for heat transfer
if OPT.BCEE == 0 %|| abs(sol.p(1)) < 1e-6 % adiabatic wall / 1e-6 small enough?
    dHF = 0;
else
    dHF = (g_aw - sol.g(1)); % static! scaled with HtE
end

% Heat flux
if OPT.BCEE == 4
    qw = INP.BCW(NS);
else
    qw = -flp.C(1)*sol.p(1)/flp.Pr(1)*EDG.HtE(NS)*EDG.muE(NS)*sqrt(EDG.Re_x(NS))/X(NS);
end

%% Boundary layer thicknesses and Dimensionless Numbers
% find BL thicknesses and dimensionless numbers, characterizing the BL

% Boundary layer velocity thickness is defined here as: 0.99*U_free-stream
% NB velocity overshoot, relative to free-stream velocity, is NOT taken into account
etaI = GRD.eta(NP);                     % eta_infinity [-]
etaref = linspace(0,etaI,1001);         % interpolated grid, spanning entire grid
soluref = pchip(GRD.eta(1:NP),sol.u(1:NP),etaref);
rhoref = EDG.rhoE(NS)./pchip(GRD.eta(1:NP),sol.c(1:NP),etaref);

% find boundary layer edge and y-coordinate [m]
[~, index] = min(abs(soluref - 0.99));
delta = sqrt(EDG.rhoE(NS)*EDG.muE(NS)*X(NS)/EDG.UE(NS))*trapz(etaref(1:index),1./rhoref(1:index));
edge = sqrt(EDG.rhoE(NS)*EDG.muE(NS)*X(NS)/EDG.UE(NS))*trapz(etaref(1:end),1./rhoref(1:end));
y = sqrt(EDG.rhoE(NS)*EDG.muE(NS)*X(NS)/EDG.UE(NS))*cumtrapz(GRD.eta(1:NP),sol.c./EDG.rhoE(NS));

% Find other BL thicknesses: BL integral properties
TERMP2 = 0;     % momentum thickness
TERMP3 = 0;     % kinetic energy thickness
TERMP4 = 0;     % enthalpy thickness
TERMP4a = 0;    % total enthalpy thickness (Schlichting)
TERMP4b = 0;    % total enthalpy loss thickness
SUM1 = 0;       % displacement thickness
SUM2 = 0;       % momentum thickness
SUM3 = 0;       % kinetic energy thickness
SUM4 = 0;       % enthalpy thickness
SUM4a = 0;      % total enthalpy thickness (Schlichting)
SUM4b = 0;      % total enthalpy loss thickness

for j = 2:NP
    TERM2 = sol.u(j)*(1 - sol.u(j));
    TERM3 = sol.u(j)*(1 - sol.u(j)^2);
    TERM4 = sol.u(j)*((EDG.HtE(NS)*sol.g(j) - (EDG.UE(NS)*sol.u(j))^2/2) - EDG.HsE(NS))/(EDG.HtE(NS)*sol.g(1) - EDG.HsE(NS)); % enthalpy thickness
    TERM4a = sol.u(j)*(sol.g(j) - 1); % enthalpy thickness (Schlichting) (total quantities, see 'Delta' in Winter and Gaudet (1973)
    TERM4b = sol.u(j)*(sol.g(j) - 1)/(sol.g(1) - 1); % total enthalpy loss thickness
    SUM1 = SUM1 + GRD.A(j)*(sol.c(j) + sol.c(j-1));%%%trapz(BLC.y(:,end),1./SOL{end}.c.*SOL{end}.u.*(SOL{end}.g-1)./(SOL{end}.g(1)-1))
    SUM2 = SUM2 + GRD.A(j)*(TERM2 + TERMP2);
    SUM3 = SUM3 + GRD.A(j)*(TERM3 + TERMP3);
    SUM4 = SUM4 + GRD.A(j)*(TERM4 + TERMP4);
    SUM4a = SUM4a + GRD.A(j)*(TERM4a + TERMP4a);
    SUM4b = SUM4b + GRD.A(j)*(TERM4b + TERMP4b);
    TERMP2 = TERM2;
    TERMP3 = TERM3;
    TERMP4 = TERM4;
    TERMP4a = TERM4a;
    TERMP4b = TERM4b;
end

BLC.y(1:NP,NS) = y;     % y-coordinate [m], obtained from transformed eta-coordinate [-]
BLC.edge(NS) = edge;    % edge of grid (for comparison with velocity thickness)
BLC.delta(NS) = delta;                                                  % velocity thickness
BLC.delta_ast(NS) = X(NS)/sqrt(EDG.Re_x(NS))*(SUM1 - sol.f(NP));        % displacement thickness
BLC.theta(NS) = X(NS)/sqrt(EDG.Re_x(NS))*SUM2;                          % momentum thickness
BLC.delta3(NS) = X(NS)/sqrt(EDG.Re_x(NS))*SUM3;                         % kinetic energy thickness
BLC.delta4(NS) = X(NS)/sqrt(EDG.Re_x(NS))*SUM4;                         % enthalpy thickness % NaN when BCEE=1/ones;
BLC.delta4a(NS) = X(NS)/sqrt(EDG.Re_x(NS))*SUM4a;                       % enthalpy thickness (Schlichting)
BLC.delta4b(NS) = X(NS)/sqrt(EDG.Re_x(NS))*SUM4b;                       % total enthalpy loss thickness % NaN when BCEE=1/ones;
BLC.delta5(1) = NaN;                                                    % Entropy thickness
BLC.H(NS) = BLC.delta_ast(NS)/BLC.theta(NS);                            % form factor
BLC.tau_wall(1) = NaN;                                                  % wall shear stress
BLC.Cf_L(NS) = 2*flp.C(1)*sol.v(1)/sqrt(EDG.Re_x(NS));                  % Local Skin friction coefficient Cf (scaled with UE(NS)); preferrably change name to Cf
BLC.Cf(NS) = BLC.Cf_L(NS)*EDG.rhoE(NS)/FRS.rhoI*(EDG.UE(NS)/FRS.UI)^2;  % formal Cf (using UI) (aerodynamics definition)
BLC.Cf_L2(NS) = BLC.Cf_L(NS)/2;
BLC.Cf2(NS) = BLC.Cf(NS)/2;
BLC.q_wall(NS) = qw;
BLC.dHF(NS) = dHF;

% Local dimensionless numbers
BLC.Re_x(NS) = EDG.Re_x(NS);                                                % reference
BLC.Re_delta_ast(NS) = sqrt(EDG.Re_x(NS))*(SUM1 - sol.f(NP));               % Reynolds displacement thickness
BLC.Re_theta(NS) = sqrt(EDG.Re_x(NS))*SUM2;                                 % Reynolds momentum thickness
BLC.Pr_x(NS) = flp.Pr(NP);                                                  % local reference Pr-number chosen at BL edge
BLC.Pe_x(NS) = BLC.Re_x(NS)*BLC.Pr_x(NS);                                   % Peclet-number along BL edge
BLC.St_x(NS) = flp.C(1)*sol.p(1)/flp.Pr(1)/sqrt(EDG.Re_x(NS))/dHF;          % Stanton-number; (INP.g_aw(NS) - sol.g(1)); % or change with dHF
BLC.Nu_x(NS) = BLC.St_x(NS)*BLC.Re_x(NS)*BLC.Pr_x(NS);                      % Nusselt-number
BLC.Ec_x(NS) = EDG.UE(NS)^2/EDG.HtE(NS)/-dHF;                               % Eckert-number, local; (sol.g(1) - INP.g_aw(NS)); % Incompressible: Ec0=EDG.MaE(NS)^2*(flp.gamma - 1); % IG; check also Kluwick
BLC.Ec_e(NS) = EDG.UE(NS)^2/EDG.TsE(NS)/EDG.CpE(NS);                        % Eckert-number; adopted from Kluwick (2017)
% or, in case of adiabatic flow: Ec = U^2/(H0 - Hw) OR Ec = U^2/(Hw - HsE)
% and in case of heat transfer: Ec = U^2/(Hw - Haw)
BLC.Br_x(NS) = BLC.Ec_x(NS)*BLC.Pr_x(NS);                                   % Brinkman-number

%% Shear stress

for j = 1:NP
    if MON.tr > 0
        BLC.tau(j,NS) = EDG.rhoE(NS)*EDG.UE(NS)^2*flp.C(j)*(1 + flp.EV(j))*sol.v(j)/sqrt(EDG.Re_x(NS)); % [N/m2]
    else
        BLC.tau(j,NS) = EDG.rhoE(NS)*EDG.UE(NS)^2*flp.C(j)*sol.v(j)/sqrt(EDG.Re_x(NS)); % [N/m2]
    end
end
BLC.tau_wall(NS) = BLC.tau(1,NS);           % [N/m2]

%% Dimensionless wall coordinates (law of the wall)

BLC.u_shear(NS) = sqrt(BLC.tau_wall(NS)*sol.c(1)/EDG.rhoE(NS)); % shear velocity, [m/s]
BLC.y_plus(1:NP,NS) = BLC.y(1:NP,NS).*BLC.u_shear(NS)./(flp.mu./EDG.rhoE(NS).*sol.c); % [-]
BLC.u_plus(1:NP,NS) = sol.u*EDG.UE(NS)./BLC.u_shear(NS); % [-]

%% Loss Coefficient - Denton (1993) - compressible flows including non-zero PG and heat transfer

% Local entropy generation inside BL per unit volume
if MON.tr > 0
    % turbulent
    BLC.Sv(1:NP,NS) = flp.C./sol.c.*(1 + flp.EV).*sol.v.^2*EDG.UE(NS)^4*EDG.rhoE(NS)^2/EDG.muE(NS)./EDG.Re_x(NS)./flp.T; % [J/K/s/m3]
else
    % laminar
    BLC.Sv(1:NP,NS) = flp.C./sol.c.*sol.v.^2*EDG.UE(NS)^4*EDG.rhoE(NS)^2/EDG.muE(NS)./EDG.Re_x(NS)./flp.T; % [J/K/s/m3]
end

% Entropy generation of BL + entropy accompanied with wall heat transfer per unit surface area
BLC.Sa(NS) = X(NS)/sqrt(EDG.Re_x(NS))*trapz(GRD.eta(1:NP),sol.c.*BLC.Sv(1:NP,NS)) + qw/flp.T(1); % [J/K/s/m2], qw = BLC.q_wall(NS); opposed to Denton, where heat transfer to blade is positive
BLC.Cd(NS) = (EDG.TsE(NS)*BLC.Sa(NS))/EDG.rhoE(NS)/EDG.UE(NS)^3; % [-]
BLC.Cd_Twall(NS) = (flp.T(1)*BLC.Sa(NS))/EDG.rhoE(NS)/EDG.UE(NS)^3;  % [-], Test for comparison (in case of heat transfer)

if NS > 1 % in case of heat transfer replace TsE with Twall! (Denton (1993))
    BLC.delta5(NS) = EDG.TsE(NS)/EDG.rhoE(NS)/EDG.UE(NS)^3*trapz(X(1:NS),BLC.Sa(1:NS)); % [m]
    % Total entropy produced (cumulative term) inside the boundary layer:
    if NS == 2 % effectively setting BLC.Cd(1) to zero:
        BLC.S_cum(NS) = BLC.S_cum(NS - 1) + (EDG.rhoE(NS)*EDG.UE(NS)^3/EDG.TsE(NS)*BLC.Cd(NS))/2*(X(NS) - X(NS - 1)); % [J/K/s/m]
    else
        BLC.S_cum(NS) = BLC.S_cum(NS - 1) + (EDG.rhoE(NS)*EDG.UE(NS)^3/EDG.TsE(NS)*BLC.Cd(NS) + EDG.rhoE(NS - 1)*EDG.UE(NS - 1)^3/EDG.TsE(NS - 1)*BLC.Cd(NS - 1))/2*(X(NS) - X(NS - 1)); % [J/K/s/m]
    end
else
    BLC.Sa(1) = 0; % change from NaN to 0; otherwise delta5 is not working
    BLC.delta5(1) = 0;
    BLC.S_cum(NS) = 0;
end

%% Pressure Gradient parameters
% For reference

BLC.Beta(NS) = 2*HVR.P2(NS)/(HVR.P2(NS) + 1); % used in Cohen and Reshotko (1956)
% Beta = +/- 2 -> infinity:
if abs(BLC.Beta(NS)) > 2
    BLC.Beta(NS) = NaN;
end

% Pohlhausen pressure gradient parameter
if MON.tr > 0 % use velocity (original) or momentum thickness (modified); see Brown and Martin (1976)
    BLC.Lambda(NS) = BLC.delta(NS)*X(NS)*HVR.P2(NS)*EDG.UE(NS)*EDG.rhoE(NS)/(EDG.muE(NS)*(1 + flp.EV(NP))); % Original Pohlhausen PG parameter
    BLC.Lambda2(NS) = BLC.theta(NS)*X(NS)*HVR.P2(NS)*EDG.UE(NS)*EDG.rhoE(NS)/(EDG.muE(NS)*(1 + flp.EV(NP))); % Modified Pohlhausen PG parameter
else
    BLC.Lambda(NS) = BLC.delta(NS)*X(NS)*HVR.P2(NS)*EDG.UE(NS)*EDG.rhoE(NS)/EDG.muE(NS); % Original Pohlhausen PG parameter
    BLC.Lambda2(NS) = BLC.theta(NS)*X(NS)*HVR.P2(NS)*EDG.UE(NS)*EDG.rhoE(NS)/EDG.muE(NS); % Modified Pohlhausen PG parameter
end

%% Store variables

% tracking variables
MON.ITE(NS) = IT + ITT; % total number of iterations at current station NS
MON.NP(NS) = NP;        % number of grid-points at current station NS
MON.MGE(NS) = MGE;      % number of Mesh Grid Extensions at current station NS

% reset parameters for next station
MGE = 0;
ITT = 0;

% solver output and fluid properties stored (calculated variables)
SOL{NS} = sol;
FLP{NS} = flp;
