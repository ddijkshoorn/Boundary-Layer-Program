%% 'FLUID' Fluid property data

% Molecular Fluid properties
FLD.R = 8.3145;                 % [J/mol/K], Universal gas constant
FLD.MW = ;                      % [g/mol = kg/kmol], Molecular Weight of
FLD.Rsg = FLD.R/FLD.MW*1e3;     % [J/kg/K], Specific Gas constant toluene

% Fluid Property data
FLD.Pcrit = ;                   % [Pa], [list model here]
FLD.Tcrit = ;                   % [K], [list model here]
FLD.vcrit =                     % [m3/kg], [list model here]
FLD.gamma = ;                   % [-], in case of calorically perfect ideal gas
FLD.SLV = ;                     % [K], Sutherland's constant for viscosity as function of temperature
FLD.SLC = ;                     % [K], Sutherland's constant for thermal conductivity as function of temperature
FLD.Pr = ;                      % [-], Molecular/Laminar Prandlt-number, or overrule this with Pr=f(mu,k,cp)
FLD.Cp = FLD.Rsg*FLD.gamma/(FLD.gamma - 1); % [J/kg/K], constant Cp for Calorically Perfect IG; or overrule this and define Cp directly:
% FLD.Cp = ;                    % [J/kg/K], constant
FLD.Cpcoeff = [];               % a0, a1, a2, a3, a4 -> Cp = a0 + a1*T + ... % Andrews (1981) IG relation; range: 100 - 590 K
FLD.lambda = ;					% [Pas], second coefficient of viscosity (dilatational viscosity)
FLD.mub = ;						% [Pas], bulk viscosity
FLD.mubcoeff = [];				% b0, b1, b2, b3, b4, etc. -> mub = b0 + b1*T + ... % Cramer (2012)?

% FluidProp
% FLD.FP = []; % Old? for initialization? or not needed anymore?
FLD.Model = ;                   % for example: 'StanMix' or 'REFPROP'
FLD.nCmp  = 1;
FLD.Cmp   = ;                   % for example: 'AIR';
FLD.Cnc   = [1 0];