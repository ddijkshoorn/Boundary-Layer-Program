function [X,NST,HVR,EDG,FRS,FLD,TCC,SOL,BLC,MON] = PRECAL(INP,OPT,SET,FLD,TCC)
% Name of this function is derived from McNally (1970)
% Precalculation of parameters (Initial Conditions and Boundary Conditions)
% for solving the system of PDE's, and Initialization of variables

% Structure:
% Initialization of variables
% Geometry
% Calculation of BL Edge BCs (for all gas models (1 to 3)  and options (1 to 6) separately)
% Calculation of BL Edge gradients 
% Calculation of BL Wall BCs, Energy Equation (EE) only (for all gas models (1 to 3)  and options (1 to 6) separately)

%% Initialization of variables

NST = length(INP.x);
SOL = [];
BLC = [];
BLC.tau = NaN(SET.NPT,NST);
BLC.DIV = NaN(SET.NPT,NST);
BLC.Sv = NaN(SET.NPT,NST);
BLC.y = NaN(SET.NPT,NST);
BLC.y_plus = NaN(SET.NPT,NST);
BLC.u_plus = NaN(SET.NPT,NST);
BLC.DIV_Loss = NaN(SET.NPT,NST);
MON = [];
MON.SEP = 0;
MON.ITE = NaN(1,NST);       % keep track of number of iterations when mesh is extended
MON.MGE = NaN(1,NST);       % keep track of number of mesh grid extensions
MON.NST = NST;              % store for plotting

%% Fully laminar, transitional or fully turbulent flow

if OPT.TRME == 4
    MON.tr = 1;
    TCC.gamma_tr = ones(1,NST);
else
    MON.tr = 0;
    TCC.gamma_tr = zeros(1,NST); % to be changed inside Transition-file
end
MON.rl = 0; % initialization for relaminarization check
MON.tr_last_check = NaN;    % last transition check at station NS
MON.rl_last_check = NaN;    % last relaminarization check at station NS
TCC.gamma_int = NaN(SET.NPT,NST);
TCC.tr = MON.tr; % superfluous but otherwise error: 'output arguments not assigned during call' (to TCC) in Transition.m

% Check input
if OPT.TRME == 1 && isempty(SET.NTR)
    fprintf('Required forced point of transition NTR for option OPT.TRME=1 is not given! Calculation stopped.')
    keyboard
    return
% elseif % add for relam check
end

if OPT.RLAM == 1 && isempty(SET.NRL)
    fprintf('Required forced point of Relaminarization NRL for option OPT.RLAM=1 is not given! Calculation stopped.')
    keyboard
    return
end

%% Calculation of input parameters - geometry

% Surface length - Surface parameter X [-]
X = NaN(1,NST);     % [-]
X(1) = 0;           % [-]
for i = 2:NST
    X(i) = X(i-1) + sqrt((INP.x(i) - INP.x(i-1))^2 + (INP.y(i) - INP.y(i-1))^2); % [-]
end

% Radius of curvature
%%% to be implemented

%% Fluid properties along Boundary Layer Edge (E)
if OPT.GASM == 1 % Calorically perfect IG
    if OPT.BCIE == 1 % given: PsI, TsI, UI, uE [-]!
        FRS.PsI = INP.PsI;                                  % [Pa]
        FRS.TsI = INP.TsI;                                  % [K]
        FRS.UI = INP.UI;                                    % [m/s]
        FRS.gammaI = FLD.gamma;                             % [-]
        FRS.CpI = FLD.Cp;                                   % [J/kg/K]
        FRS.HsI = FRS.CpI*FRS.TsI;                          % [J/kg]
        FRS.aI = sqrt(FRS.gammaI*FLD.Rsg*FRS.TsI);          % [m/s], speed of sound (SoS)
        FRS.MaI = FRS.UI/FRS.aI;                            % [-]
        FRS.TtI = FRS.TsI + FRS.UI^2/2/FLD.Cp;              % [K]
        FRS.HtI = FLD.Cp*FRS.TtI;                           % [J/kg]
        FRS.PtI = FRS.PsI*(FRS.TtI/FRS.TsI)^(FRS.gammaI/(FRS.gammaI - 1)); % [Pa]
        FRS.rhoI = FRS.PsI/FLD.Rsg/FRS.TsI;                 % [kg/m3]
        FRS.muI = FLD.mu_ref*(FRS.TsI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(FRS.TsI + FLD.SLV); % [Pas], static viscosity calc. with Sutherland's law - viscosity
        if OPT.CPRN < 1 % constant Pr-number
            FRS.kI = FRS.muI*FRS.CpI/FLD.Pr;                % [W/m/K], thermal conductivity
        else
            FRS.kI = FLD.k_ref*(INP.TsI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(INP.TsI + FLD.SLC); % [W/m/K], Sutherland's law - thermal conductivity
            FRS.PrI = FRS.muI*FRS.CpI/FRS.kI;
        end
        FRS.rhotI = FRS.PtI/FRS.TtI/FLD.Rsg;
        FRS.mutI = FLD.mu_ref*(FRS.TtI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(FRS.TtI + FLD.SLV); % [Pas], dynamic viscosity calc. with Sutherland's law - viscosity
        
        % Edge conditions
        for i = 1:NST % Note! Input uE is velocity UE(i) [m/s] scaled with UI [m/s]
            EDG.HtE(i) = FRS.HtI;                           % [J/kg]
            EDG.TtE(i) = FRS.TtI;                           % [K]
            EDG.PtE(i) = FRS.PtI;                           % [Pa]
            EDG.TsE(i) = INP.TsI*(1 - (FLD.gamma - 1)/2*(FRS.MaI^2)*(INP.uE(i)^2 - 1)); % [K], static temperature, but at SP should be equal to TtI
            EDG.UE(i) = INP.UI*INP.uE(i);                   % [m/s], non-scaled velocity
            EDG.CpE(i) = FLD.Cp;                            % [J/kg/K]
            EDG.HsE(i) = EDG.CpE(i)*EDG.TsE(i);             % [J/kg]
            EDG.muE(i) = FLD.mu_ref*(EDG.TsE(i)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(EDG.TsE(i) + FLD.SLV); % Sutherland's law - viscosity
            if OPT.CPRN < 1 % constant Pr-number
                EDG.kE(i) = EDG.muE(i)*EDG.CpE(i)/FLD.Pr;   % [W/m/K], thermal conductivity
            else
                % option for k=Suth(T) and Pr=f(mu,k,cp)
                EDG.kE(i) = FLD.k_ref*(EDG.TsE(i)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(EDG.TsE(i) + FLD.SLC); % Sutherland's law - thermal conductivity
                EDG.PrE(i) = EDG.muE(i)*EDG.CpE(i)/EDG.kE(i); % [-]
            end
            EDG.gammaE(i) = FLD.gamma;                      % [-]
            EDG.aE(i) = sqrt(EDG.gammaE(i)*FLD.Rsg*EDG.TsE(i)); % [m/s], speed of sound (SoS)
            EDG.MaE(i) = EDG.UE(i)/EDG.aE(i);               % [-], Mach-number
            EDG.PsE(i) = INP.PsI*(EDG.TsE(i)/INP.TsI)^(FLD.gamma/(FLD.gamma - 1)); % [Pa], static pressure
            EDG.rhoE(i) = EDG.PsE(i)/FLD.Rsg/EDG.TsE(i);    % [kg/m3], static density
        end
        
    elseif OPT.BCIE == 2 % given: PtI, TtI, MaI, MaE
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
        FRS.muI = FLD.mu_ref*(FRS.TsI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(FRS.TsI + FLD.SLV); % [Pas], static viscosity calc. with Sutherland's law - viscosity
        if OPT.CPRN < 1 % constant Pr-number
            FRS.kI = FRS.muI*FRS.CpI/FLD.Pr;                % [W/m/K], thermal conductivity
        else
            FRS.kI = FLD.k_ref*(FRS.TsI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(FRS.TsI + FLD.SLC); % [W/m/K], Sutherland's law - thermal conductivity
            FRS.PrI = FRS.muI*FRS.CpI/FRS.kI;
        end
        FRS.rhotI = INP.PtI/INP.TtI/FLD.Rsg;
        FRS.mutI = FLD.mu_ref*(INP.TtI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(INP.TtI + FLD.SLV); % [Pas], dynamic viscosity calc. with Sutherland's law - viscosity
        
        % Edge conditions
        for i=1:NST
            EDG.MaE(i) = INP.MaE(i);                        % [-]
            EDG.CpE(i) = FLD.Cp;                            % [J/kg/K]
            EDG.gammaE(i) = FLD.gamma;                      % [-]
            EDG.HtE(i) = FRS.HtI;                           % [J/kg]
            EDG.TtE(i) = INP.TtI;                           % [K]
            EDG.PtE(i) = FRS.PtI;                           % [Pa]
            EDG.TsE(i) = EDG.TtE(i)/(1 + (EDG.gammaE(i) - 1)/2*INP.MaE(i)^2);% [K]
            EDG.PsE(i) = INP.PtI*(EDG.TsE(i)/INP.TtI)^(EDG.gammaE(i)/(EDG.gammaE(i) - 1));% [Pa], works also for SP-flow?
            EDG.aE(i) = sqrt(EDG.gammaE(i)*FLD.Rsg*EDG.TsE(i));% [m/s]
            EDG.UE(i) = EDG.MaE(i)*EDG.aE(i);               % [m/s]
            EDG.HsE(i) = EDG.CpE(i)*EDG.TsE(i);             % [J/kg]
            EDG.rhoE(i) = EDG.PsE(i)/FLD.Rsg/EDG.TsE(i);    % [kg/m3], static density
            EDG.muE(i) = FLD.mu_ref*(EDG.TsE(i)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(EDG.TsE(i) + FLD.SLV); % Sutherland's law - viscosity
            if OPT.CPRN < 1 % constant Pr-number
                EDG.kE(i) = EDG.muE(i)*EDG.CpE(i)/FLD.Pr;   % [W/m/K], thermal conductivity
            else
                % option for k=Suth(T) and Pr=f(mu,k,cp)
                EDG.kE(i) = FLD.k_ref*(EDG.TsE(i)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(EDG.TsE(i) + FLD.SLC); % Sutherland's law - thermal conductivity
                EDG.PrE(i) = EDG.muE(i)*EDG.CpE(i)/EDG.kE(i); % [-]
            end
        end
        
    elseif OPT.BCIE == 3 % given: PtI, TtI, MaI, psE [-] (static (surface) pressure scaled with PtI)
        % MaI = f(UI,TsI), with UI and TsI unknown -> iterational procedure
        FRS.PtI = INP.PtI;                                  % [Pa]
        FRS.TtI = INP.TtI;                                  % [K]
        FRS.MaI = INP.MaI;                                  % [-]
        FRS.gammaI = FLD.gamma;                             % [-]
        FRS.CpI = FLD.Cp;                                   % [J/kg/K]
        FRS.HtI = FRS.CpI*INP.TtI;                          % [J/kg/K]
        FRS.TsI = INP.TtI/(1 + (FRS.gammaI - 1)/2*INP.MaI^2); % [K]
        FRS.aI = sqrt(FRS.gammaI*FLD.Rsg*FRS.TsI);          % [m/s]
        FRS.UI = FRS.MaI*FRS.aI;                            % [m/s]
        FRS.HsI = FRS.CpI*FRS.TsI;                          % [J/kg]
        FRS.PsI = INP.PtI*(FRS.TsI/INP.TtI)^(FRS.gammaI/(FRS.gammaI - 1));
        FRS.rhoI = FRS.PsI/FLD.Rsg/FRS.TsI;                 % [kg/m3], static density
        FRS.muI = FLD.mu_ref*(FRS.TsI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(FRS.TsI + FLD.SLV); % [Pas], static viscosity calc. with Sutherland's law - viscosity
        if OPT.CPRN < 1 % constant Pr-number
            FRS.kI = FRS.muI*FRS.CpI/FLD.Pr;                % [W/m/K], thermal conductivity
        else
            FRS.kI = FLD.k_ref*(FRS.TsI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(FRS.TsI + FLD.SLC); % [W/m/K], Sutherland's law - thermal conductivity
            FRS.PrI = FRS.muI*FRS.CpI/FRS.kI;
        end
        FRS.rhotI = INP.PtI/INP.TtI/FLD.Rsg;
        FRS.mutI = FLD.mu_ref*(INP.TtI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(INP.TtI + FLD.SLV); % [Pas], dynamic viscosity calc. with Sutherland's law - viscosity
        
        % Edge conditions
        for i=1:NST
            EDG.PsE(i) = INP.PtI*INP.psE(i);                % [Pa], static pressure at BL edge (SU2 case: static surface pressure = PsE)
            EDG.HtE(i) = FRS.HtI;                           % [J/kg]
            EDG.TtE(i) = INP.TtI;                           % [K]
            EDG.PtE(i) = FRS.PtI;                           % [Pa]
            EDG.CpE(i) = FLD.Cp;                            % [J/kg/K]
            EDG.gammaE(i) = FLD.gamma;                       % [-]
            EDG.TsE(i) = INP.TtI*(EDG.PsE(i)/INP.PtI)^((EDG.gammaE(i) - 1)/EDG.gammaE(i));% [K]
            EDG.HsE(i) = EDG.CpE(i)*EDG.TsE(i);             % [J/kg]
            EDG.UE(i) = sqrt(2*(EDG.HtE(i) - EDG.HsE(i)));  % [m/s]
            % set UE(1) to zero in case of SP flow to account for
            % rounding errors (due to the extensive set of equations above)
            if i==1 && abs(1 - EDG.PsE(i)/INP.PtI) < 1e-5
                EDG.UE(i) = 0; % set UE(1) to zero at SP
            end
            EDG.rhoE(i) = EDG.PsE(i)/FLD.Rsg/EDG.TsE(i);    % [kg/m3], static density
            EDG.aE(i) = sqrt(EDG.gammaE(i)*FLD.Rsg*EDG.TsE(i));% [m/s]
            EDG.MaE(i) = EDG.UE(i)/EDG.aE(i);               % [-], Added for computation of Cd
            EDG.muE(i) = FLD.mu_ref*(EDG.TsE(i)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(EDG.TsE(i) + FLD.SLV); % Sutherland's law - viscosity
            if OPT.CPRN < 1 % constant Pr-number
                EDG.kE(i) = EDG.muE(i)*EDG.CpE(i)/FLD.Pr;   % [W/m/K], thermal conductivity
            else
                % option for k=Suth(T) and Pr=f(mu,k,cp)
                EDG.kE(i) = FLD.k_ref*(EDG.TsE(i)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(EDG.TsE(i) + FLD.SLC); % Sutherland's law - thermal conductivity
                EDG.PrE(i) = EDG.muE(i)*EDG.CpE(i)/EDG.kE(i); % [-]
            end
        end
    end
    
elseif OPT.GASM == 2 % Thermally perfect IG
    if OPT.BCIE == 1 % given: PsI, TsI, UI, uE (scaled: UE/UI)
        FRS.PsI = INP.PsI;                                  % [Pa], static Pressure
        FRS.TsI = INP.TsI;                                  % [K], static Temperature
        FRS.UI = INP.UI;                                    % [m/s]
        FRS.HsI = FLD.Cpcoeff(1)*INP.TsI + FLD.Cpcoeff(2)/2*INP.TsI^2 + FLD.Cpcoeff(3)/3*INP.TsI^3 + FLD.Cpcoeff(4)/4*INP.TsI^4 + FLD.Cpcoeff(5)/5*INP.TsI^5; % [J/kg], static enthalpy, reference: 0 J/kg at 0 K
        FRS.HtI = FRS.HsI + INP.UI^2/2;                     % [J/kg], total enthalpy
        FRS.CpI = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*INP.TsI + FLD.Cpcoeff(3)*INP.TsI^2 + FLD.Cpcoeff(4)*INP.TsI^3 + FLD.Cpcoeff(5)*INP.TsI^4; % [J/kg/K], Andrews (1981) IG relation; range: 100 - 590 K
        FRS.gammaI = FRS.CpI/(FRS.CpI - FLD.Rsg);           % [-]
        FRS.aI = sqrt(FRS.gammaI*FLD.Rsg*INP.TsI);          % [m/s], Speed of Sound
        FRS.MaI = INP.UI/FRS.aI;                            % [-]
        FRS.rhoI = INP.PsI/FLD.Rsg/INP.TsI;                 % [kg/m3], static density
        FRS.muI = FLD.mu_ref*(INP.TsI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(INP.TsI + FLD.SLV); % [Pas], static viscosity calc. with Sutherland's law - viscosity
        FRS.kI = FLD.k_ref*(INP.TsI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(INP.TsI + FLD.SLC); % [W/m/K], Sutherland's law - thermal conductivity
        FRS.PrI = FRS.muI*FRS.CpI/FRS.kI;

        % Iterate to find TtI (Newton-Raphson, solve H0(T0) - h(T) - u^2/2 = 0)
        U = INP.UI;
        T = INP.TsI;
        h = FRS.HsI;
        To = 0;
        Tn = T*(1 + (FRS.gammaI - 1)/2*FRS.MaI^2); % estimate T0
        IT2 = 0;
        while abs(Tn - To) > 1e-3%1e-3 % 1e-5
            IT2 = IT2 + 1;
            if IT2 > SET.ITMAX2
                fprintf('IT exceeded ITMAX = %g at station = %g \n',SET.ITMAX2, 1)
                keyboard
                return
            end
            To = Tn;
            H0 = FLD.Cpcoeff(1)*To + FLD.Cpcoeff(2)/2*To^2 + FLD.Cpcoeff(3)/3*To^3 + FLD.Cpcoeff(4)/4*To^4 + FLD.Cpcoeff(5)/5*To^5;
            F = H0 - h - U^2/2;
            dH0 = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*To + FLD.Cpcoeff(3)*To^2 + FLD.Cpcoeff(4)*To^3 + FLD.Cpcoeff(5)*To^4;
            dF = dH0;
            Tn = To - F/dF;
        end
        FRS.TtI = Tn;                                       % [K], total Temperature
        FRS.PtI = INP.PsI*exp(1/FLD.Rsg*(FLD.Cpcoeff(1)*log(FRS.TtI/INP.TsI) + FLD.Cpcoeff(2)*(FRS.TtI - INP.TsI) + FLD.Cpcoeff(3)/2*(FRS.TtI^2 - INP.TsI^2) + FLD.Cpcoeff(4)/3*(FRS.TtI^3 - INP.TsI^3) + FLD.Cpcoeff(5)/4*(FRS.TtI^4 - INP.TsI^4))); % [Pa]
        FRS.rhotI = FRS.PtI/FRS.TtI/FLD.Rsg;                % [kg/m3]
        FRS.mutI = FLD.mu_ref*(FRS.TtI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(FRS.TtI + FLD.SLV); % [Pas], dynamic viscosity calc. with Sutherland's law - viscosity
        FRS.ktI = FLD.k_ref*(FRS.TtI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(FRS.TtI + FLD.SLC); % [W/m/K], Sutherland's law - thermal conductivity
        FRS.CptI = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*FRS.TtI + FLD.Cpcoeff(3)*FRS.TtI^2 + FLD.Cpcoeff(4)*FRS.TtI^3 + FLD.Cpcoeff(5)*FRS.TtI^4; % [J/kg/K]
        FRS.gammatI = FRS.CptI/(FRS.CptI - FLD.Rsg);        % [-]
        
        % Edge conditions
        for i=1:NST
            EDG.UE(i) = INP.UI*INP.uE(i);                   % [m/s], non-scaled velocity
            EDG.PtE(i) = FRS.PtI;                           % [K], total pressure
            EDG.TtE(i) = FRS.TtI;                           % [K], total Temperature
            EDG.HtE(i) = FRS.HtI;                           % [J/kg]
            EDG.HsE(i) = EDG.HtE(i) - EDG.UE(i)^2/2;        % [J/kg]
            
            if i == 1 && EDG.UE(1) < 1e-5
                % Stagnation Point Flow
                EDG.TsE(i) = FRS.TtI;
                EDG.PsE(i) = FRS.PtI;
%                 EDG.UE(i) = 0; % recalibrate if needed (in this case not needed since UE is input)
                EDG.MaE(i) = 0;
                EDG.rhoE(i) = FRS.rhotI;
                EDG.CpE(i) = FRS.CptI;
                EDG.gammaE(i) = FRS.gammatI;
                EDG.aE(i) = sqrt(FRS.gammatI*FLD.Rsg*FRS.TtI);
                EDG.muE(i) = FRS.mutI;
                EDG.kE(i) = FRS.ktI;
                EDG.PrE(i) = FRS.mutI*FRS.CptI/FRS.ktI;
            else
                % Iterate to find TsE (Newton-Raphson, solve H0 - h(T) - u^2/2 = 0)
                if i == 1
                    Tn = INP.TsI; % Tn = T0*(1 + (FLD.gamma - 1)/2*Ma^2)^(-1); % estimate T
                else
                    Tn = EDG.TsE(i - 1);
                end
                U = INP.UI*INP.uE(i);
                H0 = FRS.HtI;
                To = 0;
                IT2 = 0;
                while abs(Tn - To) > 1e-3
                    IT2 = IT2 + 1;
                    if IT2 > SET.ITMAX2
                        fprintf('IT exceeded ITMAX = %g at station = %g \n',SET.ITMAX2, 1)
                        keyboard
                        return
                    end
                    To = Tn;
                    Cp = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*To + FLD.Cpcoeff(3)*To^2 + FLD.Cpcoeff(4)*To^3 + FLD.Cpcoeff(5)*To^4;
                    h = FLD.Cpcoeff(1)*To + FLD.Cpcoeff(2)/2*To^2 + FLD.Cpcoeff(3)/3*To^3 + FLD.Cpcoeff(4)/4*To^4 + FLD.Cpcoeff(5)/5*To^5;
                    F = H0 - h - U^2/2;
                    dF = -Cp;
                    Tn = To - F/dF;
                end
                
                EDG.TsE(i) = Tn; % [K], static Temperature
                % The following relation follows from: ds = 0 = Cp/T dT - R 1/P dP
                EDG.PsE(i) = FRS.PtI/(exp(1/FLD.Rsg*(FLD.Cpcoeff(1)*log(FRS.TtI/EDG.TsE(i)) + FLD.Cpcoeff(2)*(FRS.TtI - EDG.TsE(i)) + FLD.Cpcoeff(3)/2*(FRS.TtI^2 - EDG.TsE(i)^2) + FLD.Cpcoeff(4)/3*(FRS.TtI^3 - EDG.TsE(i)^3) + FLD.Cpcoeff(5)/4*(FRS.TtI^4 - EDG.TsE(i)^4))));
                EDG.rhoE(i) = EDG.PsE(i)/FLD.Rsg/EDG.TsE(i);    % [kg/m3], static density
                EDG.CpE(i) = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*EDG.TsE(i) + FLD.Cpcoeff(3)*EDG.TsE(i)^2 + FLD.Cpcoeff(4)*EDG.TsE(i)^3  + FLD.Cpcoeff(5)*EDG.TsE(i)^4; % [J/kg/K], Andrews (1981) IG relation; range: 100 - 590 K
                EDG.gammaE(i) = EDG.CpE(i)/(EDG.CpE(i) - FLD.Rsg); % [-]
                EDG.aE(i) = sqrt(EDG.gammaE(i)*FLD.Rsg*EDG.TsE(i));% [m/s]
                EDG.MaE(i) = EDG.UE(i)/EDG.aE(i);               % [-]
                EDG.muE(i) = FLD.mu_ref*(EDG.TsE(i)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(EDG.TsE(i) + FLD.SLV); % Sutherland's law - viscosity
                EDG.kE(i) = FLD.k_ref*(EDG.TsE(i)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(EDG.TsE(i) + FLD.SLC); % Sutherland's law - thermal conductivity
                EDG.PrE(i) = EDG.muE(i)*EDG.CpE(i)/EDG.kE(i); % [-]
            end
        end
        
    elseif OPT.BCIE == 2 % given: PtI, TtI, MaI, MaE
        % MaI = f(UI,TsI), with UI and TsI unknown -> iterational procedure
        FRS.PtI = INP.PtI;                                  % [Pa], total Pressure
        FRS.TtI = INP.TtI;                                  % [K], total Temperature
        FRS.MaI = INP.MaI;                                  % [-]
        FRS.HtI = FLD.Cpcoeff(1)*INP.TtI + FLD.Cpcoeff(2)/2*INP.TtI^2 + FLD.Cpcoeff(3)/3*INP.TtI^3 + FLD.Cpcoeff(4)/4*INP.TtI^4 + FLD.Cpcoeff(5)/5*INP.TtI^5; % [J/kg], Andrews (1981) IG relation; range: 100 - 590 K
        
        % Iterate to find TsI (Newton-Raphson, solve H0 - h(T) - u(T,Ma)^2/2 = 0)
        H0 = FRS.HtI;
        T0 = INP.TtI;
        Ma = INP.MaI;
        R = FLD.Rsg;
        To = 0;
        Tn = T0*(1 + (FLD.gamma - 1)/2*Ma^2)^(-1); % estimate T
        IT2 = 0;
        while abs(Tn - To) > 1e-5 % 1e-3
            IT2 = IT2 + 1;
            if IT2 > SET.ITMAX2
                fprintf('IT exceeded ITMAX = %g at station = %g \n',SET.ITMAX2, 1)
                keyboard
                return
            end
            To = Tn;
            h = FLD.Cpcoeff(1)*To + FLD.Cpcoeff(2)/2*To^2 + FLD.Cpcoeff(3)/3*To^3 + FLD.Cpcoeff(4)/4*To^4 + FLD.Cpcoeff(5)/5*To^5;
            Cp = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*To + FLD.Cpcoeff(3)*To^2 + FLD.Cpcoeff(4)*To^3 + FLD.Cpcoeff(5)*To^4;
            g = Cp - R;
            dg = FLD.Cpcoeff(2) + 2*FLD.Cpcoeff(3)*To + 3*FLD.Cpcoeff(4)*To^2 + 4*FLD.Cpcoeff(5)*To^3;
            f = Cp*To;
            df = FLD.Cpcoeff(1) + 2*FLD.Cpcoeff(2)*To + 3*FLD.Cpcoeff(3)*To^2 + 4*FLD.Cpcoeff(4)*To^3 + 5*FLD.Cpcoeff(5)*To^4;
            F = H0 - h - R*Ma^2/2*(Cp/(Cp - R)*To);
            dF = -Cp - R*Ma^2/2*((df*g - dg*f)/g^2);
            Tn = To - F/dF;
        end
        FRS.TsI = Tn;                                       % [K], total Temperature
        % The following relation follows from: ds = 0 = Cp/T dT - R 1/P dP
        FRS.PsI = INP.PtI/exp(1/FLD.Rsg*(FLD.Cpcoeff(1)*log(INP.TtI/FRS.TsI) + FLD.Cpcoeff(2)*(INP.TtI - FRS.TsI) + FLD.Cpcoeff(3)/2*(INP.TtI^2 - FRS.TsI^2) + FLD.Cpcoeff(4)/3*(INP.TtI^3 - FRS.TsI^3) + FLD.Cpcoeff(5)/4*(INP.TtI^4 - FRS.TsI^4))); % [Pa]
        FRS.CpI = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*FRS.TsI + FLD.Cpcoeff(3)*FRS.TsI^2 + FLD.Cpcoeff(4)*FRS.TsI^3 + FLD.Cpcoeff(5)*FRS.TsI^4; % [J/kg/K]
        FRS.gammaI = FRS.CpI/(FRS.CpI - FLD.Rsg);           % [-]
        FRS.aI = sqrt(FRS.gammaI*FLD.Rsg*FRS.TsI);          % [m/s], Speed of Sound
        FRS.UI = FRS.aI*INP.MaI;                            % [m/s]
        FRS.HsI = FRS.HtI - FRS.UI^2/2;                     % [J/kg], total enthalpy, reference is 0J/kg at 0K
        FRS.rhoI = FRS.PsI/FLD.Rsg/FRS.TsI;                 % [kg/m3], static density
        FRS.muI = FLD.mu_ref*(FRS.TsI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(FRS.TsI + FLD.SLV); % [Pas], static viscosity calc. with Sutherland's law - viscosity
        FRS.kI = FLD.k_ref*(FRS.TsI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(FRS.TsI + FLD.SLC); % [W/m/K], Sutherland's law - thermal conductivity
        FRS.PrI = FRS.muI*FRS.CpI/FRS.kI;
        FRS.rhotI = INP.PtI/INP.TtI/FLD.Rsg;
        FRS.mutI = FLD.mu_ref*(INP.TtI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(INP.TtI + FLD.SLV); % [Pas], dynamic viscosity calc. with Sutherland's law - viscosity
        FRS.ktI = FLD.k_ref*(INP.TtI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(INP.TtI + FLD.SLC); % [W/m/K], Sutherland's law - thermal conductivity
        FRS.CptI = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*INP.TtI + FLD.Cpcoeff(3)*INP.TtI^2 + FLD.Cpcoeff(4)*INP.TtI^3 + FLD.Cpcoeff(5)*INP.TtI^4; % [J/kg/K]
        FRS.gammatI = FRS.CptI/(FRS.CptI - FLD.Rsg);           % [-]
        
        % Edge conditions
        for i=1:NST
            EDG.PtE(i) = FRS.PtI;                           % [K], total pressure
            EDG.TtE(i) = FRS.TtI;                           % [K], total Temperature
            EDG.HtE(i) = FRS.HtI;                           % [J/kg]
            EDG.MaE(i) = INP.MaE(i);                        % [Pa]
            
            if i == 1 && INP.MaE(1) < 1e-5 %%% be careful with this condition
                % Stagnation Point Flow
                EDG.PsE(i) = INP.PtI;
                EDG.TsE(i) = INP.TtI;
                EDG.HsE(i) = FRS.HtI;
                EDG.UE(i) = 0;
                % EDG.MaE(i) = 0; % reset value if it was not equal to zero
                EDG.rhoE(i) = FRS.rhotI;
                EDG.CpE(i) = FRS.CptI;
                EDG.gammaE(i) = FRS.gammatI;
                EDG.aE(i) = sqrt(FRS.gammatI*FLD.Rsg*INP.TtI);
                EDG.muE(i) = FRS.mutI;
                EDG.kE(i) = FRS.ktI;
                EDG.PrE(i) = FRS.mutI*FRS.CptI/FRS.ktI;
                % Define for Newton-Raphson scheme applied at following stations
                H0 = FRS.HtI;
                R = FLD.Rsg;
                A0 = FLD.Cpcoeff(1)*log(INP.TtI) + FLD.Cpcoeff(2)*INP.TtI + FLD.Cpcoeff(3)/2*INP.TtI^2 + FLD.Cpcoeff(4)/3*INP.TtI^3 + FLD.Cpcoeff(5)/4*INP.TtI^4;
            else
                % Iterate to find TsE (Newton-Raphson, solve H0 - h(T) - u(T,Ma)^2/2 = 0)
                if i == 1
                    H0 = FRS.HtI;
                    R = FLD.Rsg;
                    A0 = FLD.Cpcoeff(1)*log(INP.TtI) + FLD.Cpcoeff(2)*INP.TtI + FLD.Cpcoeff(3)/2*INP.TtI^2 + FLD.Cpcoeff(4)/3*INP.TtI^3 + FLD.Cpcoeff(5)/4*INP.TtI^4;
                    Tn = FRS.TsI;
                else
                    Tn = EDG.TsE(i - 1);
                end
                Ma = EDG.MaE(i);
                To = 0;
                IT2 = 0;
                while abs(Tn - To) > 1e-5 % 1e-3
                    IT2 = IT2 + 1;
                    if IT2 > SET.ITMAX2
                        fprintf('IT exceeded ITMAX = %g at station = %g \n',SET.ITMAX2, 1)
                        keyboard
                        return
                    end
                    To = Tn;
                    h = FLD.Cpcoeff(1)*To + FLD.Cpcoeff(2)/2*To^2 + FLD.Cpcoeff(3)/3*To^3 + FLD.Cpcoeff(4)/4*To^4 + FLD.Cpcoeff(5)/5*To^5;
                    Cp = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*To + FLD.Cpcoeff(3)*To^2 + FLD.Cpcoeff(4)*To^3 + FLD.Cpcoeff(5)*To^4;
                    g = Cp - R;
                    dg = FLD.Cpcoeff(2) + 2*FLD.Cpcoeff(3)*To + 3*FLD.Cpcoeff(4)*To^2 + 4*FLD.Cpcoeff(5)*To^3;
                    f = Cp*To;
                    df = FLD.Cpcoeff(1) + 2*FLD.Cpcoeff(2)*To + 3*FLD.Cpcoeff(3)*To^2 + 4*FLD.Cpcoeff(4)*To^3 + 5*FLD.Cpcoeff(5)*To^4;
                    F = H0 - h - R*Ma^2/2*(Cp/(Cp - R)*To);
                    dF = -Cp - R*Ma^2/2*((df*g - dg*f)/g^2);
                    Tn = To - F/dF;
                end

                EDG.TsE(i) = Tn;                                % [K], total Temperature
                A = FLD.Cpcoeff(1)*log(Tn) + FLD.Cpcoeff(2)*Tn + FLD.Cpcoeff(3)/2*Tn^2 + FLD.Cpcoeff(4)/3*Tn^3 + FLD.Cpcoeff(5)/4*Tn^4;
                EDG.PsE(i) = INP.PtI/exp(1/FLD.Rsg*(A0 - A));   % [Pa]
                EDG.rhoE(i) = EDG.PsE(i)/FLD.Rsg/EDG.TsE(i);    % [kg/m3], static density
                EDG.CpE(i) = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*EDG.TsE(i) + FLD.Cpcoeff(3)*EDG.TsE(i)^2 + FLD.Cpcoeff(4)*EDG.TsE(i)^3 + FLD.Cpcoeff(5)*EDG.TsE(i)^4; % [J/kg/K], Andrews (1981) IG relation; range: 100 - 590 K
                EDG.gammaE(i) = EDG.CpE(i)/(EDG.CpE(i) - FLD.Rsg); % [-]
                EDG.aE(i) = sqrt(EDG.gammaE(i)*FLD.Rsg*EDG.TsE(i)); % [m/s]
                EDG.UE(i) = EDG.MaE(i)*EDG.aE(i);               % [-]
                EDG.HsE(i) = EDG.HtE(i) - EDG.UE(i)^2/2;        % [J/kg]
                EDG.muE(i) = FLD.mu_ref*(EDG.TsE(i)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(EDG.TsE(i) + FLD.SLV); % Sutherland's law - viscosity
                EDG.kE(i) = FLD.k_ref*(EDG.TsE(i)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(EDG.TsE(i) + FLD.SLC); % Sutherland's law - thermal conductivity
                EDG.PrE(i) = EDG.muE(i)*EDG.CpE(i)/EDG.kE(i);   % [-]
            end
        end
        
    elseif OPT.BCIE == 3 % given: PtI, TtI, MaI, pE [-] (Iterational procedure)
        % MaI = f(UI,TsI), with UI and TsI unknown -> iterational procedure
        FRS.PtI = INP.PtI;                                  % [Pa], total Pressure
        FRS.TtI = INP.TtI;                                  % [K], total Temperature
        FRS.MaI = INP.MaI;                                  % [-]
        FRS.HtI = FLD.Cpcoeff(1)*INP.TtI + FLD.Cpcoeff(2)/2*INP.TtI^2 + FLD.Cpcoeff(3)/3*INP.TtI^3 + FLD.Cpcoeff(4)/4*INP.TtI^4 + FLD.Cpcoeff(5)/5*INP.TtI^5; % [J/kg], Andrews (1981) IG relation; range: 100 - 590 K
        
        % Iterate to find TsI (Newton-Raphson, solve H0 - h(T) - u(T,Ma)^2/2 = 0)
        H0 = FRS.HtI;
        T0 = INP.TtI;
        Ma = INP.MaI;
        R = FLD.Rsg;
        To = 0;
        Tn = T0*(1 + (FLD.gamma - 1)/2*Ma^2)^(-1); % estimate T
        IT2 = 0;
        while abs(Tn - To) > 1e-4
            IT2 = IT2 + 1;
            if IT2 > SET.ITMAX2
                fprintf('IT exceeded ITMAX = %g at station = %g \n',SET.ITMAX2, 1)
                keyboard
                return
            end
            To = Tn;
            h = FLD.Cpcoeff(1)*To + FLD.Cpcoeff(2)/2*To^2 + FLD.Cpcoeff(3)/3*To^3 + FLD.Cpcoeff(4)/4*To^4 + FLD.Cpcoeff(5)/5*To^5;
            Cp = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*To + FLD.Cpcoeff(3)*To^2 + FLD.Cpcoeff(4)*To^3 + FLD.Cpcoeff(5)*To^4;
            g = Cp - R;
            dg = FLD.Cpcoeff(2) + 2*FLD.Cpcoeff(3)*To + 3*FLD.Cpcoeff(4)*To^2 + 4*FLD.Cpcoeff(5)*To^3;
            f = Cp*To;
            df = FLD.Cpcoeff(1) + 2*FLD.Cpcoeff(2)*To + 3*FLD.Cpcoeff(3)*To^2 + 4*FLD.Cpcoeff(4)*To^3 + 5*FLD.Cpcoeff(5)*To^4;
            F = H0 - h - R*Ma^2/2*(Cp/(Cp - R)*To);
            dF = -Cp - R*Ma^2/2*((df*g - dg*f)/g^2);
            Tn = To - F/dF;
        end

        FRS.TsI = Tn;                                       % [K], total Temperature
        FRS.PsI = INP.PtI/exp(1/FLD.Rsg*(FLD.Cpcoeff(1)*log(INP.TtI/FRS.TsI) + FLD.Cpcoeff(2)*(INP.TtI - FRS.TsI) + FLD.Cpcoeff(3)/2*(INP.TtI^2 - FRS.TsI^2) + FLD.Cpcoeff(4)/3*(INP.TtI^3 - FRS.TsI^3) + FLD.Cpcoeff(5)/4*(INP.TtI^4 - FRS.TsI^4))); % [Pa]
        FRS.CpI = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*FRS.TsI + FLD.Cpcoeff(3)*FRS.TsI^2 + FLD.Cpcoeff(4)*FRS.TsI^3 + FLD.Cpcoeff(5)*FRS.TsI^4; % [J/kg/K]
        FRS.gammaI = FRS.CpI/(FRS.CpI - FLD.Rsg);           % [-]
        FRS.aI = sqrt(FRS.gammaI*FLD.Rsg*FRS.TsI);          % [m/s], Speed of Sound
        FRS.UI = FRS.aI*INP.MaI;                            % [m/s]
        FRS.HsI = FRS.HtI - FRS.UI^2/2;                     % [J/kg], total enthalpy, reference is 0J/kg at 0K
        FRS.rhoI = FRS.PsI/FLD.Rsg/FRS.TsI;                 % [kg/m3], static density
        FRS.muI = FLD.mu_ref*(FRS.TsI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(FRS.TsI + FLD.SLV); % [Pas], static viscosity calc. with Sutherland's law - viscosity
        FRS.kI = FLD.k_ref*(FRS.TsI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(FRS.TsI + FLD.SLC); % [W/m/K], Sutherland's law - thermal conductivity
        FRS.PrI = FRS.muI*FRS.CpI/FRS.kI;
        FRS.rhotI = INP.PtI/INP.TtI/FLD.Rsg;
        FRS.mutI = FLD.mu_ref*(INP.TtI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(INP.TtI + FLD.SLV); % [Pas], dynamic viscosity calc. with Sutherland's law - viscosity
        FRS.ktI = FLD.k_ref*(INP.TtI/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(INP.TtI + FLD.SLC); % [W/m/K], Sutherland's law - thermal conductivity
        FRS.CptI = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*INP.TtI + FLD.Cpcoeff(3)*INP.TtI^2 + FLD.Cpcoeff(4)*INP.TtI^3 + FLD.Cpcoeff(5)*INP.TtI^4; % [J/kg/K]
        FRS.gammatI = FRS.CptI/(FRS.CptI - FLD.Rsg);           % [-]
        
        % Edge conditions
        for i=1:NST
            EDG.PtE(i) = FRS.PtI;                           % [K], total pressure
            EDG.TtE(i) = FRS.TtI;                           % [K], total Temperature
            EDG.HtE(i) = FRS.HtI;                           % [J/kg]
            EDG.PsE(i) = INP.PtI*INP.psE(i);                % [Pa]
            
            if i == 1 && abs(EDG.PsE(1) - INP.PtI)/INP.PtI < 1e-5  % difference must be smaller than...
                % Stagnation Point Flow
                EDG.TsE(i) = INP.TtI;
                EDG.HsE(i) = FRS.HtI;
                EDG.UE(i) = 0;   % set to zero
                EDG.MaE(i) = 0;  % set to zero
                EDG.rhoE(i) = FRS.rhotI;%INP.PsE(i)/EDG.TsE(i)/FLD.Rsg;
                EDG.CpE(i) = FRS.CptI;
                EDG.gammaE(i) = FRS.gammatI;
                EDG.aE(i) = sqrt(FRS.gammatI*FLD.Rsg*INP.TtI);
                EDG.muE(i) = FRS.mutI;
                EDG.kE(i) = FRS.ktI;
                EDG.PrE(i) = FRS.mutI*FRS.CptI/FRS.ktI;
                % Define for Newton-Raphson scheme applied at following stations
                T0 = INP.TtI;
                P0 = INP.PtI;
                R = FLD.Rsg;
                A0 = FLD.Cpcoeff(1)*log(T0) + FLD.Cpcoeff(2)*T0 + FLD.Cpcoeff(3)/2*T0^2 + FLD.Cpcoeff(4)/3*T0^3 + FLD.Cpcoeff(5)/4*T0^4;
            else % Newton-Raphson iteration scheme
                % Iterate to find TsE (Newton-Raphson, solve A0 - A(T) - R*ln(P/P0) = 0, where A'(T) = Cp/T)
                if i == 1 % in case of (Flat) Plate Flow (Leading Edge conditions)
                    T0 = INP.TtI;
                    P0 = INP.PtI;
                    R = FLD.Rsg;
                    A0 = FLD.Cpcoeff(1)*log(T0) + FLD.Cpcoeff(2)*T0 + FLD.Cpcoeff(3)/2*T0^2 + FLD.Cpcoeff(4)/3*T0^3 + FLD.Cpcoeff(5)/4*T0^4; % [J/kg]
                    Tn = FRS.TsI;
                else
                    Tn = EDG.TsE(i - 1);
                end
                P = EDG.PsE(i);
                To = 0;
                IT2 = 0;
                while abs(Tn - To) > 1e-3
                    IT2 = IT2 + 1;
                    if IT2 > SET.ITMAX2
                        fprintf('IT exceeded ITMAX = %g at station = %g \n',SET.ITMAX2, 1)
                        keyboard
                        return
                    end
                    To = Tn;
                    A = FLD.Cpcoeff(1)*log(To) + FLD.Cpcoeff(2)*To + FLD.Cpcoeff(3)/2*To^2 + FLD.Cpcoeff(4)/3*To^3 + FLD.Cpcoeff(5)/4*To^4;
                    F = A0 - A - R*log(P0/P);
                    dF = -1*(FLD.Cpcoeff(1)/To + FLD.Cpcoeff(2) + FLD.Cpcoeff(3)*To + FLD.Cpcoeff(4)*To^2 + FLD.Cpcoeff(5)*To^3);
                    Tn = To - F/dF;
                end
                EDG.TsE(i) = Tn;                                % [K], total Temperature
                EDG.rhoE(i) = EDG.PsE(i)/FLD.Rsg/EDG.TsE(i);    % [kg/m3], static density
                EDG.CpE(i) = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*EDG.TsE(i) + FLD.Cpcoeff(3)*EDG.TsE(i)^2 + FLD.Cpcoeff(4)*EDG.TsE(i)^3 + FLD.Cpcoeff(5)*EDG.TsE(i)^4; % [J/kg/K], Andrews (1981) IG relation; range: 100 - 590 K
                EDG.gammaE(i) = EDG.CpE(i)/(EDG.CpE(i) - FLD.Rsg); % [-]
                EDG.aE(i) = sqrt(EDG.gammaE(i)*FLD.Rsg*EDG.TsE(i)); % [m/s]
                EDG.HsE(i) = FLD.Cpcoeff(1)*EDG.TsE(i) + FLD.Cpcoeff(2)/2*EDG.TsE(i)^2 + FLD.Cpcoeff(3)/3*EDG.TsE(i)^3 + FLD.Cpcoeff(4)/4*EDG.TsE(i)^4 + FLD.Cpcoeff(5)/5*EDG.TsE(i)^5; % [J/kg], static enthalpy, reference: 0J/kg at 0K
                EDG.UE(i) = sqrt(2*(EDG.HtE(i) - EDG.HsE(i)));  % [m/s]
                EDG.MaE(i) = EDG.UE(i)/EDG.aE(i);               % [-]
                EDG.muE(i) = FLD.mu_ref*(EDG.TsE(i)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLV)/(EDG.TsE(i) + FLD.SLV); % Sutherland's law - viscosity
                EDG.kE(i) = FLD.k_ref*(EDG.TsE(i)/FLD.Tref_S)^(3/2)*(FLD.Tref_S + FLD.SLC)/(EDG.TsE(i) + FLD.SLC); % Sutherland's law - thermal conductivity
                EDG.PrE(i) = EDG.muE(i)*EDG.CpE(i)/EDG.kE(i);   % [-]
            end
        end
    end
    
elseif OPT.GASM == 3 % Non-Ideal Gas
    
    %% FluidProp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Init_FluidProp
    % Define the thermodynamic model in Fluidprop
    MON.FP_ErrorMsg = invoke(FLD.FP,'SetFluid_M',FLD.Model,FLD.nCmp,FLD.Cmp,FLD.Cnc);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if OPT.BCIE == 1 % given: PsI, TsI, UI, uE (ratio)
        FRS.PsI = INP.PsI;                                      % [Pa]
        FRS.TsI = INP.TsI;                                      % [K]
        FRS.UI = INP.UI;                                        % [m/s]
        FRS.HsI = 1e3*FLD.FP.Enthalpy('PT',INP.PsI/1e5,INP.TsI - 273.15); % [J/kg]
        FRS.aI = FLD.FP.SoundSpeed('PT',INP.PsI/1e5,INP.TsI - 273.15); % [m/s]
        FRS.MaI = INP.UI/FRS.aI;                                % [Ma]
        FRS.HtI = FRS.HsI + INP.UI^2/2;                         % [J/kg]
        FRS.s = FLD.FP.Entropy('PT',INP.PsI/1e5,INP.TsI - 273.15)*1e3;  % [J/kg/K]
        FRS.TtI = FLD.FP.Temperature('hs',FRS.HtI/1e3,FRS.s/1e3) + 273.15; % [K]
        FRS.PtI = FLD.FP.Pressure('hs',FRS.HtI/1e3,FRS.s/1e3)*1e5;      % [Pa] 
        FRS.rhotI = FLD.FP.Density('Ps',FRS.PtI/1e5,FRS.s/1e3);         % [kg/m3]
        FRS.mutI = FLD.FP.Viscosity('Ps',FRS.PtI/1e5,FRS.s/1e3);        % [Pas], static viscosity
        FRS.rhoI = FLD.FP.Density('PT',INP.PsI/1e5,INP.TsI - 273.15);  % [kg/m3]
        FRS.CpI = FLD.FP.HeatCapP('hs',FRS.HsI/1e3,FRS.s/1e3)*1e3;      % [J/kg/K]
        FRS.gammaI = FRS.CpI/(FRS.CpI - FLD.Rsg);               % [-]
        FRS.muI = FLD.FP.Viscosity('hs',FRS.HsI/1e3,FRS.s/1e3);         % [Pas], static viscosity
        FRS.kI = FLD.FP.ThermCond('hs',FRS.HsI/1e3,FRS.s/1e3);          % [W/m/K]
        FRS.GAMI = FLD.FP.Gamma('hs',FRS.HsI/1e3,FRS.s/1e3);            % [-]
        FRS.ZI = INP.PsI/FLD.Rsg/FRS.rhoI/INP.TsI;              % [-]
        
        for i=1:NST
            EDG.UE(i) = INP.UI*INP.uE(i);                           % [m/s], non-scaled velocity
            EDG.PtE(i) = FRS.PtI;                                   % [Pa], total Temperature
            EDG.TtE(i) = FRS.TtI;                                   % [K]
            EDG.HtE(i) = FRS.HtI;                                   % [J/kg]
            EDG.HsE(i) = EDG.HtE(i) - EDG.UE(i)^2/2;                % [J/kg]
            EDG.TsE(i) = FLD.FP.Temperature('hs',EDG.HsE(i)/1e3,FRS.s/1e3) + 273.15;% [K]
            EDG.PsE(i) = FLD.FP.Pressure('hs',EDG.HsE(i)/1e3,FRS.s/1e3)*1e5;% [Pa]
            EDG.CpE(i) = FLD.FP.HeatCapP('hs',EDG.HsE(i)/1e3,FRS.s/1e3)*1e3;% [J/kg/K]
            EDG.gammaE(i) = EDG.CpE(i)/(EDG.CpE(i) - FLD.Rsg);      % [-]
            EDG.rhoE(i) = FLD.FP.Density('hs',EDG.HsE(i)/1e3,FRS.s/1e3);    % [kg/m3]
            EDG.muE(i) = FLD.FP.Viscosity('hs',EDG.HsE(i)/1e3,FRS.s/1e3);   % [Pas], static viscosity
            EDG.kE(i) = FLD.FP.ThermCond('hs',EDG.HsE(i)/1e3,FRS.s/1e3);    % [W/m/K]
            EDG.PrE(i) = EDG.muE(i).*EDG.CpE(i)./EDG.kE(i);                 % [-]
            EDG.aE(i) = FLD.FP.SoundSpeed('hs',EDG.HsE(i)/1e3,FRS.s/1e3);   % [m/s]
            EDG.MaE(i) = EDG.UE(i)/EDG.aE(i);                       % [-]
            EDG.GAM(i) = FLD.FP.Gamma('hs',EDG.HsE(i)/1e3,FRS.s/1e3);       % [-], Fundamental Derivative of Gas Dynamics
            EDG.ZE(i) = EDG.PsE(i)/FLD.Rsg/EDG.rhoE(i)/EDG.TsE(i);  % [-]
        end
        
    elseif OPT.BCIE == 2 % given: PtI, TtI, MaI, MaE
        FRS.PtI = INP.PtI;                                      % [Pa]
        FRS.TtI = INP.TtI;                                      % [K]
        FRS.MaI = INP.MaI;                                      % [-]
        FRS.HtI = FLD.FP.Enthalpy('PT',INP.PtI/1e5,INP.TtI - 273.15)*1e3;% [J/kg]
        FRS.s = FLD.FP.Entropy('PT',INP.PtI/1e5,INP.TtI - 273.15)*1e3;  % [J/kg/K]
        FRS.rhotI = FLD.FP.Density('Ps',FRS.PtI/1e5,FRS.s/1e3);         % [kg/m3]
        FRS.mutI = FLD.FP.Viscosity('Ps',FRS.PtI/1e5,FRS.s/1e3);        % [Pas], static viscosity
        
        % determine: TsI, HsI and PsI
        PIn = INP.PtI;
        PIo = 0;
        IT2 = 0;
        while abs(PIn - PIo) > 1e-1
            IT2 = IT2 + 1;
            if IT2 > SET.ITMAX3
                fprintf('IT exceeded ITMAX = %g at station = %g \n',SET.ITMAX3, 1)
                keyboard
                return
            end
            PIo = PIn;
            aI = FLD.FP.SoundSpeed('Ps',PIo/10^5,FRS.s/1e3);
            HsI = FRS.HtI - (aI*INP.MaI)^2/2;
            PIn = FLD.FP.Pressure('hs',HsI/1e3,FRS.s/1e3)*10^5;
        end
        
        MON.IT_PRECAL_FRS = IT2;
        FRS.PsI = PIn;                                      % [Pa]
        FRS.aI = FLD.FP.SoundSpeed('Ps',FRS.PsI/1e5,FRS.s/1e3);     % [m/s]
        FRS.UI = FRS.MaI*FRS.aI;                            % [m/s]
        FRS.HsI = FRS.HtI - FRS.UI^2/2;                     % [J/kg]
        FRS.TsI = FLD.FP.Temperature('Ps',FRS.PsI/1e5,FRS.s/1e3) + 273.15;% [K]
        FRS.rhoI = FLD.FP.Density('Ps',FRS.PsI/1e5,FRS.s/1e3);      % [kg/m3]
        FRS.CpI = FLD.FP.HeatCapP('Ps',FRS.PsI/1e5,FRS.s/1e3)*1e3;  % [J/kg/K]
        FRS.gammaI = FRS.CpI/(FRS.CpI - FLD.Rsg);           % [-]
        FRS.muI = FLD.FP.Viscosity('Ps',FRS.PsI/1e5,FRS.s/1e3);     % [Pas], static viscosity
        FRS.kI = FLD.FP.ThermCond('Ps',FRS.PsI/1e5,FRS.s/1e3);      % [W/m/K]
        FRS.GAMI = FLD.FP.Gamma('Ps',FRS.PsI/1e5,FRS.s/1e3);        % [-]
        FRS.ZI = FRS.PsI/FLD.Rsg/FRS.rhoI/FRS.TsI;          % [-]
        
        for i = 1:NST
            EDG.MaE(i) = INP.MaE(i);                        % [-]
            EDG.PtE(i) = FRS.PtI;                           % [Pa], total Temperature
            EDG.TtE(i) = FRS.TtI;                           % [K]
            EDG.HtE(i) = FRS.HtI;                           % [J/kg]
            
            % Determine PsE by iterational procedure
            if i == 1
                PEn = INP.PtI;
            else
                % estimate (calorically perfect ideal gas)
                PEn = INP.PtI*(1 + (FRS.gammaI - 1)/2*INP.MaE(i)^2)^(-FRS.gammaI/(FRS.gammaI - 1));
            end
            PEo = 0;
            IT2 = 0;
            while abs(PEn - PEo) > 1e1;10;1e-1;
                IT2 = IT2 + 1;
                if IT2 > SET.ITMAX3
                    fprintf('IT exceeded ITMAX = %g at station = %g \n',SET.ITMAX3, i)
                    keyboard
                    return
                end
                PEo = PEn;
                aE = FLD.FP.SoundSpeed('Ps',PEo/10^5,FRS.s/1e3);
                HsE = FRS.HtI - (aE*INP.MaE(i))^2/2;
                PEn = FLD.FP.Pressure('hs',HsE/1e3,FRS.s/1e3)*10^5;
            end
            MON.IT_PRECAL_EDG(i) = IT2;
            EDG.PsE(i) = PEn;                                       % [Pa]
            EDG.aE(i) = FLD.FP.SoundSpeed('Ps',EDG.PsE(i)/10^5,FRS.s/1e3);  % [m/s]
            EDG.UE(i) = EDG.MaE(i)*EDG.aE(i);                       % [m/s]
            EDG.HsE(i) = EDG.HtE(i) - EDG.UE(i)^2/2;                % [J/kg]
            EDG.TsE(i) = FLD.FP.Temperature('Ph',EDG.PsE(i)/10^5,EDG.HsE(i)/1e3) + 273.15;% [K]
            EDG.CpE(i) = FLD.FP.HeatCapP('Ph',EDG.PsE(i)/10^5,EDG.HsE(i)/1e3)*1e3;% [J/kg/K]
            EDG.gammaE(i) = EDG.CpE(i)/(EDG.CpE(i) - FLD.Rsg);      % [-]
            EDG.rhoE(i) = FLD.FP.Density('Ph',EDG.PsE(i)/10^5,EDG.HsE(i)/1e3);% [kg/m3]
            EDG.muE(i) = FLD.FP.Viscosity('Ph',EDG.PsE(i)/10^5,EDG.HsE(i)/1e3);% [Pas], static viscosity
            EDG.kE(i) = FLD.FP.ThermCond('hs',EDG.HsE(i)/1e3,FRS.s/1e3);    % [W/m/K]
            EDG.PrE(i) = EDG.muE(i).*EDG.CpE(i)./EDG.kE(i);                 % [-]
            EDG.GAM(i) = FLD.FP.Gamma('Ph',EDG.PsE(i)/10^5,EDG.HsE(i)/1e3);% [-], Fundamental Derivative of Gas Dynamics
            EDG.ZE(i) = EDG.PsE(i)/FLD.Rsg/EDG.rhoE(i)/EDG.TsE(i);  % [-]
        end
        
    elseif OPT.BCIE == 3  % given: PtI, TtI, MaI, PsE or psE=PsE/PtI (stator blade simulation)
        FRS.PtI = INP.PtI;                                      % [Pa]
        FRS.TtI = INP.TtI;                                      % [K]
        FRS.MaI = INP.MaI;                                      % [-]
        FRS.HtI = FLD.FP.Enthalpy('PT',INP.PtI/1e5,INP.TtI - 273.15)*1e3;% [J/kg]
        FRS.s = FLD.FP.Entropy('PT',INP.PtI/1e5,INP.TtI - 273.15)*1e3;  % [J/kg/K]
        % Total Fluid properties added for C0 (26-03-2020)
        FRS.rhotI = FLD.FP.Density('Ps',FRS.PtI/1e5,FRS.s/1e3);      % [kg/m3]
        FRS.mutI = FLD.FP.Viscosity('Ps',FRS.PtI/1e5,FRS.s/1e3);     % [Pas], static viscosity
        PIn = INP.PtI;
        PIo = 0;
        IT2 = 0;
        while abs(PIn - PIo) > 1e-1
            IT2 = IT2 + 1;
                        if IT2 > SET.ITMAX3
                            fprintf('IT exceeded ITMAX = %g at station = %g \n',SET.ITMAX3, 1)
                            keyboard
                            return
                        end
            PIo = PIn;
            aI = FLD.FP.SoundSpeed('Ps',PIo/10^5,FRS.s/1e3);
            HsI = FRS.HtI - (aI*INP.MaI)^2/2;
            PIn = FLD.FP.Pressure('hs',HsI/1e3,FRS.s/1e3)*10^5;
        end
        
        MON.IT_PRECAL_FRS = IT2;
        FRS.PsI = PIn;                                      % [Pa]
        FRS.aI = FLD.FP.SoundSpeed('Ps',FRS.PsI/1e5,FRS.s/1e3);     % [m/s]
        FRS.UI = FRS.MaI*FRS.aI;                            % [m/s]
        FRS.HsI = FRS.HtI - FRS.UI^2/2;                     % [J/kg]
        FRS.TsI = FLD.FP.Temperature('Ps',FRS.PsI/1e5,FRS.s/1e3) + 273.15;% [K]
        FRS.rhoI = FLD.FP.Density('Ps',FRS.PsI/1e5,FRS.s/1e3);      % [kg/m3]
        FRS.CpI = FLD.FP.HeatCapP('Ps',FRS.PsI/1e5,FRS.s/1e3)*1e3;  % [J/kg/K]
        FRS.gammaI = FRS.CpI/(FRS.CpI - FLD.Rsg);           % [-]
        FRS.muI = FLD.FP.Viscosity('Ps',FRS.PsI/1e5,FRS.s/1e3);     % [Pas], static viscosity
        FRS.kI = FLD.FP.ThermCond('Ps',FRS.PsI/1e5,FRS.s/1e3);      % [W/m/K]
        FRS.GAMI = FLD.FP.Gamma('Ps',FRS.PsI/1e5,FRS.s/1e3);        % [-]
        FRS.ZI = FRS.PsI/FLD.Rsg/FRS.rhoI/FRS.TsI;          % [-]
        
        for i=1:NST
            EDG.PsE(i) = INP.PtI*INP.psE(i);                        % [Pa], static pressure
            EDG.PtE(i) = INP.PtI;                                   % [Pa]
            EDG.TtE(i) = INP.TtI;                                   % [K]
            EDG.HtE(i) = FRS.HtI;                                   % [J/kg]
            EDG.TsE(i) = FLD.FP.Temperature('Ps',EDG.PsE(i)/10^5,FRS.s/1e3) + 273.15;% [K]
            EDG.aE(i) = FLD.FP.SoundSpeed('Ps',EDG.PsE(i)/10^5,FRS.s/1e3);  % [m/s]
            EDG.HsE(i) = FLD.FP.Enthalpy('Ps',EDG.PsE(i)/10^5,FRS.s/1e3)*1e3;% [J/kg]
            if i == 1 && abs(EDG.PsE(i) - INP.PtI) < 1e-5 % then SP flow
                EDG.UE(i) = 0;
                EDG.MaE(i) = 0;
            else
                EDG.UE(i) = sqrt(2*(EDG.HtE(i) - EDG.HsE(i)));      % [m/s]
                EDG.MaE(i) = EDG.UE(i)/EDG.aE(i);                   % [-]
            end
            EDG.CpE(i) = FLD.FP.HeatCapP('Ps',EDG.PsE(i)/10^5,FRS.s/1e3)*1e3;% [J/kg/K]
            EDG.gammaE(i) = EDG.CpE(i)/(EDG.CpE(i) - FLD.Rsg);      % [-]
            EDG.rhoE(i) = FLD.FP.Density('Ps',EDG.PsE(i)/10^5,FRS.s/1e3);   % [kg/m3]
            EDG.muE(i) = FLD.FP.Viscosity('Ps',EDG.PsE(i)/10^5,FRS.s/1e3);  % [Pas], static viscosity
            EDG.kE(i) = FLD.FP.ThermCond('Ps',EDG.PsE(i)/10^5,FRS.s/1e3);   % [W/m/K]
            EDG.PrE(i) = EDG.muE(i).*EDG.CpE(i)./EDG.kE(i);                 % [-]
            EDG.GAM(i) = FLD.FP.Gamma('Ps',EDG.PsE(i)/10^5,FRS.s/1e3);      % [-], Fundamental Derivative of Gas Dynamics
            EDG.ZE(i) = EDG.PsE(i)/FLD.Rsg/EDG.rhoE(i)/EDG.TsE(i);  % [-], Compressibility factor Z=v/v0
        end
        
    elseif OPT.BCIE == 4 % given: TsI, s, UE (Howarth) (for Ts-diagram input)
        FRS.TsI = INP.TsI;                                          % [K]
        FRS.s = INP.s;                                              % [J/kg/K]
        FRS.MaI = INP.MaI;                                          % [Ma]
        FRS.HsI = FLD.FP.Enthalpy('Ts',INP.TsI - 273.15,INP.s)*1e3; % [J/kg]
        FRS.PsI = FLD.FP.Pressure('Ts',INP.TsI - 273.15,INP.s)*1e5; % [Pa]
        FRS.aI = FLD.FP.SoundSpeed('Ts',INP.TsI - 273.15,INP.s);    % [m/s]
        FRS.UI = INP.MaI*FRS.aI;                                    % [m/s]
        FRS.HtI = FRS.HsI + FRS.UI^2/2;                             % [J/kg]
        FRS.TtI = FLD.FP.Temperature('hs',FRS.HtI/1e3,INP.s) + 273.15;% [K]
        FRS.PtI = FLD.FP.Pressure('hs',FRS.HtI/1e3,INP.s)*1e5;      % [Pa]
        FRS.rhotI = FLD.FP.Density('hs',FRS.HtI/1e3,INP.s);         % [kg/m3]
        FRS.mutI = FLD.FP.Viscosity('hs',FRS.HtI/1e3,INP.s);        % [Pas], static viscosity
        FRS.rhoI = FLD.FP.Density('Ts',INP.TsI - 273.15,INP.s);     % [kg/m3]
        FRS.CpI = FLD.FP.HeatCapP('Ts',INP.TsI - 273.15,INP.s)*1e3; % [J/kg/K]
        FRS.gammaI = FRS.CpI/(FRS.CpI - FLD.Rsg);                   % [-]
        FRS.muI = FLD.FP.Viscosity('Ts',INP.TsI - 273.15,INP.s);    % [Pas], static viscosity
        FRS.kI = FLD.FP.ThermCond('Ts',INP.TsI - 273.15,INP.s);     % [W/m/K]
        FRS.GAMI = FLD.FP.Gamma('Ts',INP.TsI - 273.15,INP.s);       % [-]
        FRS.ZI = FRS.PsI/FLD.Rsg/FRS.rhoI/INP.TsI;                  % [-]
        
        for i = 1:NST
%             EDG.UE(i) = FRS.UI*(1 - X(i));                          % [m/s], Howarth's Flow velocity input
            EDG.UE(i) = INP.UI*INP.uE(i);                               % [m/s], non-scaled velocity
            EDG.PtE(i) = FRS.PtI;                                       % [Pa]
            EDG.TtE(i) = FRS.TtI;                                       % [K]
            EDG.HtE(i) = FRS.HtI;                                       % [J/kg]
            EDG.HsE(i) = EDG.HtE(i) - EDG.UE(i)^2/2;                    % [J/kg]
            EDG.TsE(i) = FLD.FP.Temperature('hs',EDG.HsE(i)/1e3,INP.s) + 273.15;% [K]
            EDG.PsE(i) = FLD.FP.Pressure('hs',EDG.HsE(i)/1e3,INP.s)*1e5;% [Pa], static pressure
            EDG.CpE(i) = FLD.FP.HeatCapP('hs',EDG.HsE(i)/1e3,INP.s)*1e3;% [J/kg/K]
            EDG.gammaE(i) = EDG.CpE(i)/(EDG.CpE(i) - FLD.Rsg);          % [-]
            EDG.rhoE(i) = FLD.FP.Density('hs',EDG.HsE(i)/1e3,INP.s);    % [kg/m3]
            EDG.muE(i) = FLD.FP.Viscosity('hs',EDG.HsE(i)/1e3,INP.s);   % [Pas], static viscosity
            EDG.kE(i) = FLD.FP.ThermCond('hs',EDG.HsE(i)/1e3,INP.s);    % [W/m/K]
            EDG.PrE(i) = EDG.muE(i).*EDG.CpE(i)./EDG.kE(i);             % [-]
            EDG.aE(i) = FLD.FP.SoundSpeed('hs',EDG.HsE(i)/1e3,INP.s);   % [m/s], SoS
            EDG.MaE(i) = EDG.UE(i)/EDG.aE(i);                           % [Ma]
            EDG.GAM(i) = FLD.FP.Gamma('hs',EDG.HsE(i)/1e3,INP.s);       % [-], Fundamental Derivative of Gas Dynamics
            EDG.ZE(i) = EDG.PsE(i)/FLD.Rsg/EDG.rhoE(i)/EDG.TsE(i);      % [-], Compressibility factor
        end
        
    elseif OPT.BCIE == 5 % given: PsI, TsI, MaI; UE (Howarth)
        FRS.PsI = INP.PsI;                                          % [Pa]
        FRS.TsI = INP.TsI;                                          % [K]
        FRS.MaI = INP.MaI;                                          % [-]
        FRS.HsI = 1e3*FLD.FP.Enthalpy('PT',INP.PsI/1e5,INP.TsI - 273.15);% [J/kg]
        FRS.aI = FLD.FP.SoundSpeed('PT',INP.PsI/1e5,INP.TsI - 273.15);% [m/s]
        FRS.UI = INP.MaI*FRS.aI;                                    % [m/s]
        FRS.HtI = FRS.HsI + FRS.UI^2/2;                             % [J/kg]
        FRS.s = FLD.FP.Entropy('PT',INP.PsI/1e5,INP.TsI - 273.15)*1e3;      % [J/kg/K]
        FRS.TtI = FLD.FP.Temperature('hs',FRS.HtI/1e3,FRS.s/1e3) + 273.15;  % [K], 'hs' is not working! (total enthalpy too high?)
        FRS.PtI = FLD.FP.Pressure('hs',FRS.HtI/1e3,FRS.s/1e3)*1e5;          % [Pa], 'hs' is not working! (total enthalpy too high?)
        FRS.rhotI = FLD.FP.Density('Ps',FRS.PtI/1e5,FRS.s/1e3);             % [kg/m3]
        FRS.mutI = FLD.FP.Viscosity('Ps',FRS.PtI/1e5,FRS.s/1e3);            % [Pas], static viscosity
        FRS.rhoI = FLD.FP.Density('PT',INP.PsI/1e5,INP.TsI - 273.15);% [kg/m3]
        FRS.CpI = FLD.FP.HeatCapP('PT',INP.PsI/1e5,INP.TsI - 273.15)*1e3;% [J/kg/K]
        FRS.gammaI = FRS.CpI/(FRS.CpI - FLD.Rsg);                   % [-]
        FRS.muI = FLD.FP.Viscosity('PT',INP.PsI/1e5,INP.TsI - 273.15);% [Pas], static viscosity
        FRS.kI = FLD.FP.ThermCond('PT',INP.PsI/1e5,INP.TsI - 273.15);% [W/m/K]
        FRS.GAMI = FLD.FP.Gamma('PT',INP.PsI/1e5,INP.TsI - 273.15); % [-]
        FRS.ZI = INP.PsI/FLD.Rsg/FRS.rhoI/INP.TsI;                  % [-], larger than 1.0 !!! ?
        
        for i=1:NST
%             EDG.UE(i) = FRS.UI*(1 - X(i));                          % [m/s], Howarth's Flow velocity input
            EDG.UE(i) = INP.UI*INP.uE(i);                           % [m/s], non-scaled velocity
            EDG.PtE(i) = FRS.PtI;                                   % [Pa]
            EDG.TtE(i) = FRS.TtI;                                   % [K]
            EDG.HtE(i) = FRS.HtI;                                   % [J/kg]
            EDG.HsE(i) = EDG.HtE(i) - EDG.UE(i)^2/2;                % [J/kg]
            EDG.TsE(i) = FLD.FP.Temperature('hs',EDG.HsE(i)/1e3,FRS.s/1e3) + 273.15; % [K]
            EDG.PsE(i) = FLD.FP.Pressure('hs',EDG.HsE(i)/1e3,FRS.s/1e3)*1e5;% [Pa], static pressure
            EDG.CpE(i) = FLD.FP.HeatCapP('hs',EDG.HsE(i)/1e3,FRS.s/1e3)*1e3;% [J/kg/K]
            EDG.gammaE(i) = EDG.CpE(i)/(EDG.CpE(i) - FLD.Rsg);      % [-]
            EDG.rhoE(i) = FLD.FP.Density('hs',EDG.HsE(i)/1e3,FRS.s/1e3);    % [kg/m3]
            EDG.muE(i) = FLD.FP.Viscosity('hs',EDG.HsE(i)/1e3,FRS.s/1e3);   % [Pas], static viscosity
            EDG.kE(i) = FLD.FP.ThermCond('hs',EDG.HsE(i)/1e3,FRS.s/1e3);    % [W/m/K]
            EDG.PrE(i) = EDG.muE(i).*EDG.CpE(i)./EDG.kE(i);                 % [-]
            EDG.aE(i) = FLD.FP.SoundSpeed('hs',EDG.HsE(i)/1e3,FRS.s/1e3);   % [m/s], SoS
            EDG.MaE(i) = EDG.UE(i)/EDG.aE(i);                       % [Ma]
            EDG.GAM(i) = FLD.FP.Gamma('hs',EDG.HsE(i)/1e3,FRS.s/1e3);       % [-], Fundamental Derivative of Gas Dynamics
            EDG.ZE(i) = EDG.PsE(i)/FLD.Rsg/EDG.rhoE(i)/EDG.TsE(i);  % [-], Compressibility factor
        end
        
    elseif OPT.BCIE == 6 % given: PsI, vI, MaI; UE (Howarth)
        FRS.PsI = INP.PsI;                                          % [Pa]
        FRS.vI = INP.vI;                                            % [kg/m3]
        FRS.MaI = INP.MaI;                                          % [Ma]
        FRS.HsI = 1e3*FLD.FP.Enthalpy('Pv',INP.PsI/1e5,INP.vI);     % [J/kg]
        FRS.TsI = FLD.FP.Temperature('Pv',INP.PsI/1e5,INP.vI) + 273.15; % [K]
        FRS.aI = FLD.FP.SoundSpeed('Pv',INP.PsI/1e5,INP.vI);        % [m/s]
        FRS.UI = INP.MaI*FRS.aI;                                    % [m/s]
        FRS.HtI = FRS.HsI + FRS.UI^2/2;                             % [J/kg]
        FRS.s = FLD.FP.Entropy('Pv',INP.PsI/1e5,INP.vI)*1e3;                % [J/kg/K]
        FRS.TtI = FLD.FP.Temperature('hs',FRS.HtI/1e3,FRS.s/1e3) + 273.15;  % [K]
        FRS.PtI = FLD.FP.Pressure('hs',FRS.HtI/1e3,FRS.s/1e3)*1e5;          % [Pa]
        FRS.rhotI = FLD.FP.Density('Ps',FRS.PtI/1e5,FRS.s/1e3);             % [kg/m3]
        FRS.mutI = FLD.FP.Viscosity('Ps',FRS.PtI/1e5,FRS.s/1e3);            % [Pas], static viscosity
        FRS.rhoI = FLD.FP.Density('Pv',INP.PsI/1e5,INP.vI);         % [kg/m3]
        FRS.CpI = FLD.FP.HeatCapP('Pv',INP.PsI/1e5,INP.vI)*1e3;     % [J/kg/K]
        FRS.gammaI = FRS.CpI/(FRS.CpI - FLD.Rsg);                   % [-]
        FRS.muI = FLD.FP.Viscosity('Pv',INP.PsI/1e5,INP.vI);        % [Pas], static viscosity
        FRS.kI = FLD.FP.ThermCond('Pv',INP.PsI/1e5,INP.vI);         % [W/m/K]
        FRS.GAMI = FLD.FP.Gamma('Pv',INP.PsI/1e5,INP.vI);           % [-]
        FRS.ZI = INP.PsI/FLD.Rsg/FRS.rhoI/FRS.TsI;                  % [-]
        
        for i=1:NST
%             EDG.UE(i) = FRS.UI*(1 - X(i));                          % [m/s], Howarth's Flow velocity input
            EDG.UE(i) = INP.UI*INP.uE(i);                           % [m/s], non-scaled velocity
            EDG.PtE(i) = FRS.PtI;                                   % [Pa]
            EDG.TtE(i) = FRS.TtI;                                   % [K]
            EDG.HtE(i) = FRS.HtI;                                   % [J/kg]
            EDG.HsE(i) = EDG.HtE(i) - EDG.UE(i)^2/2;                % [J/kg]
            EDG.TsE(i) = FLD.FP.Temperature('hs',EDG.HsE(i)/1e3,FRS.s/1e3) + 273.15; % [K]
            EDG.PsE(i) = FLD.FP.Pressure('hs',EDG.HsE(i)/1e3,FRS.s/1e3)*1e5;% [Pa], static pressure
            EDG.CpE(i) = FLD.FP.HeatCapP('hs',EDG.HsE(i)/1e3,FRS.s/1e3)*1e3;% [J/kg/K]
            EDG.gammaE(i) = EDG.CpE(i)/(EDG.CpE(i) - FLD.Rsg);      % [-]
            EDG.rhoE(i) = FLD.FP.Density('hs',EDG.HsE(i)/1e3,FRS.s/1e3);    % [kg/m3]
            EDG.muE(i) = FLD.FP.Viscosity('hs',EDG.HsE(i)/1e3,FRS.s/1e3);   % [Pas], static viscosity
            EDG.kE(i) = FLD.FP.ThermCond('hs',EDG.HsE(i)/1e3,FRS.s/1e3);    % [W/m/K]
            EDG.PrE(i) = EDG.muE(i).*EDG.CpE(i)./EDG.kE(i);                 % [-]
            EDG.aE(i) = FLD.FP.SoundSpeed('hs',EDG.HsE(i)/1e3,FRS.s/1e3);   % [m/s], SoS
            EDG.MaE(i) = EDG.UE(i)/EDG.aE(i);                       % [Ma]
            EDG.GAM(i) = FLD.FP.Gamma('hs',EDG.HsE(i)/1e3,FRS.s/1e3);       % [-], Fundamental Derivative of Gas Dynamics
            EDG.ZE(i) = EDG.PsE(i)/FLD.Rsg/EDG.rhoE(i)/EDG.TsE(i);  % [-], Compressibility factor
        end
    end
end

%% Calculation of Reynolds-number
for i = 1:NST
    EDG.Re_x(i) = EDG.rhoE(i)*EDG.UE(i)*X(i)/EDG.muE(i);
end

%% Calculation of dimensionless numbers in free-stream
FRS.Re_L = FRS.rhoI*FRS.UI*INP.L/FRS.muI;           % [-], characteristic length L
FRS.Pr_L = FRS.muI*FRS.CpI/FRS.kI;                  % [-]
FRS.Pe_L = FRS.Re_L*FRS.Pr_L;                       % [-]
FRS.Ec_L = FRS.UI^2/FRS.HsI;                        % [-]
FRS.Br_L = FRS.Ec_L*FRS.Pr_L;                       % [-]

%% Calculation of pressure-gradient parameter P=P2
HVR.P1 = zeros(1,NST);
HVR.P2 = zeros(1,NST);
% In case of Stagnation Point (SP) flow
if EDG.UE(1) == 0 % include a tolerance here for input of CFD data if needed
    HVR.P2(1) = 1;
end
HVR.P1(1) = 0.5*(HVR.P2(1) + 1); % correct for SP flow also

% Calculation of help variables
HVR.CEL = zeros(1,NST); % CEL(1) = 0; !!!
HVR.P1P = zeros(1,NST);
HVR.P2P = zeros(1,NST);
HVR.P1P(1) = HVR.P1(1) + HVR.CEL(1);
HVR.P2P(1) = HVR.P2(1) + HVR.CEL(1);

if OPT.GRAD == 1
    DUDX = LAGRANGE(EDG.UE,X); % seems best option, obtained from Cebeci (2002)
elseif OPT.GRAD == 2
    DUDX = GRADNT(EDG.UE,X);
elseif OPT.GRAD == 3
    % [DFDX,DF2DX2,yy,xx,yp] = SPLINEGRAD(EDG.UE,X,1);
    [DUDX,~,~,~,~] = SPLINEGRAD(EDG.UE,X,1); 
end

for i = 2:NST
    HVR.P2(i) = X(i)/EDG.UE(i)*DUDX(i);
    HVR.P1(i) = 0.5*(1 + HVR.P2(i) + X(i)*(EDG.rhoE(i)*EDG.muE(i) - EDG.rhoE(i-1)*EDG.muE(i-1))/(X(i) - X(i-1))/EDG.rhoE(i)/EDG.muE(i));
    % Calculation of help variables
    HVR.CEL(i) = 0.5*(X(i) + X(i-1))/(X(i) - X(i-1));
    HVR.P1P(i) = HVR.P1(i) + HVR.CEL(i);
    HVR.P2P(i) = HVR.P2(i) + HVR.CEL(i);
end

%% Wall Boundary Condition Energy Equation
if  OPT.BCEE == 0 % adiabatic: derivative of enthalpy ratio p = 0 [-]
    HVR.alpha0 = 0;
    HVR.alpha1 = 1;
    HVR.WW = INP.BCW; % zeros
elseif OPT.BCEE == 1 % enthalpy ratio g [-]
    HVR.alpha0 = 1;
    HVR.alpha1 = 0;
    HVR.WW = INP.BCW;
elseif OPT.BCEE == 2 % derivative of enthalpy ratio [-]
    HVR.alpha0 = 0;
    HVR.alpha1 = 1;
    HVR.WW = INP.BCW;
elseif OPT.BCEE == 3 % temperature ratio [-]
    HVR.alpha0 = 1;
    HVR.alpha1 = 0;
    if OPT.GASM == 1 % cp = constant, thus gwall = hwall/htE = Twall/TtI (since CpE = Cp_wall)
        % g(1) = T/T; assuming Cp = constant gives CpE = Cp_wall
        HVR.WW = INP.BCW; % original solution
    elseif OPT.GASM == 2 % Cp = Cp(T)
        Twall = FRS.TtI*INP.BCW;
        Cpwall = FLD.Cpcoeff(1) + FLD.Cpcoeff(2)*Twall + FLD.Cpcoeff(3)*Twall.^2 + FLD.Cpcoeff(4)*Twall.^3 + FLD.Cpcoeff(5)*Twall.^4; % [J/kg/K], Andrews (1981) IG relation
        HVR.WW = Cpwall.*Twall/FRS.HtI;
        % HVR.WW = % hwall/HtI = (Cpwall(Twall)*Twall)/(Cp(TtI)*TtI) -> hwall = HtI*h*(Cp*T)/(Cpwall*Twall)
    elseif OPT.GASM == 3
        % h = f(P,T); since dP/dy = 0 and T is given
        HVR.WW = NaN(1,NST);
        for i = 1:NST
            HVR.WW(i) = FLD.FP.Enthalpy('PT',EDG.PsE(i)/1e5,FRS.TtI*INP.BCW(i) - 273.15)*1e3/FRS.HtI; % [-], wall enthapy ratio
        end
    end
elseif OPT.BCEE == 4 % heat flux q [W/m2]; inaccurate for most cases, since wall temperature needs to be estimated
    HVR.alpha0 = 0;
    HVR.alpha1 = 1;
    if OPT.GASM == 1 % cp = constant
        if OPT.COMP > 0 % compressible flow, standard
            if OPT.CPRN < 1 % constant Pr-number
                if OPT.CCRP > 0 % General/variable Chapman-Rubesin parameter (depending on fluid properties)
                    HVR.WW = -X./sqrt(EDG.Re_x).*FLD.Pr./EDG.muE./FRS.HtI.*INP.BCW; % approximation, Cw not considered, assumed that Pr/Cw ~ Pr
                else % Predefined (and) constant Chapman-Rubesin parameter, exact!
                    HVR.WW = -X./sqrt(EDG.Re_x).*FLD.Pr./FLD.C./EDG.muE./FRS.HtI.*INP.BCW; % exact relation
                end
            else % variable Pr-number
                HVR.WW = -X./sqrt(EDG.Re_x).*EDG.PrE./EDG.muE./FRS.HtI.*INP.BCW; % surprisingly good approximation according to nonidael gas and air nozzle simulations: Pr/Cw ~ PrE
            end
        else % incompressible flow, Cw = 1 (given), exact!
            HVR.WW = -X./sqrt(EDG.Re_x).*FLD.Pr./EDG.muE./FRS.HtI.*INP.BCW; % exact relation
        end
    else % approximation only
%         HVR.WW = -FLD.Pr*INP.BCW./FRS.HtI./EDG.muE./sqrt(EDG.Re_x).*X; % OLD, missing Cw
        HVR.WW = -X./sqrt(EDG.Re_x).*EDG.PrE./EDG.muE./FRS.HtI.*INP.BCW; % surprisingly good approximation according to nonidael gas (MM) and air nozzle simulations: Pr/Cw ~ PrE
    end

    % set Stagnation Point
    if INP.BCW(1) == 0 % needed for SP
        HVR.WW(1) = 0;
    end

end
