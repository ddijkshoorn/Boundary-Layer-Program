%% AIR Fluid property data

% Molecular Fluid properties
FLD.R = 8.3145;                 % [J/mol/K], Universal gas constant
FLD.MW = 28.97;                 % [g/mol = kg/kmol], Molecular Weight of dry Air
% FLD.Rsg = 287.0;               % [J/kg/K], (Cebeci (2004)!) Specific Gas constant Air
% FLD.Rsg = FLD.R/FLD.MW;        % [kJ/kg/K], Specific Gas constant Air
FLD.Rsg = FLD.R/FLD.MW*1e3;     % [J/kg/K], Specific Gas constant Air

% Fluid Property data
FLD.Pcrit = 38.500352295368310e5;         % [Pa], Air 'REFPROP'
FLD.Tcrit = -140.3267060519201 + 273.15;  % [K], Air 'REFPROP'
FLD.vcrit = 0.002906056291283;  % [m3/kg], Air 'REFPROP' -> replace with StanMix?
FLD.gamma = 1.4;                % [-], Specific heat ratio (Cp/Cv)
FLD.SLV = 110.4;111;110.33;     % [K],  (Groot, 2018) Sutherland's constant (later ones: White (2006) and Cebeci (2002)!)
FLD.SLC = 194;                  % [K], (Groot, 2018; and White, 2006) Sutherland's constant for thermal conductivity
FLD.Tref_S = 273;               % [K], (Groot, 2018; and White, 2006) reference temperature Sutherland's Law
FLD.mu_ref = 1.716e-5;          % [Pas], (Groot, 2018; and White, 2006) reference viscosity Sutherland's Law
FLD.k_ref = 0.0241;             % [W/m/K], (Groot, 2018; and White, 2006) reference thermal conductivity Sutherland's Law
FLD.Pr = 0.72;                  % [-], (Cebeci, 2002) Prandlt-number
% FLD.Cp = 1004.3;                % [J/kg/K], (Cebeci, 2002) Constant Pressure Heat capacity
FLD.Cp = FLD.Rsg*FLD.gamma/(FLD.gamma - 1); % [J/kg/K]
% FLD.Cp = FLD.Rsg*FLD.gamma/(FLD.gamma - 1)*1e3; % [J/kg/K]
FLD.Cpcoeff = [1022.5294853 -0.1758625 4.020605e-4 -4.8640623e-8 0]; % a0, a1, a2, a3, a4, etc. -> Cp = a0 + a1*T + ... % Andrews (1981) IG relation; range: 100 - 590 K

% FluidProp
% FLD.FP = []; % not used anymore?
FLD.Model = 'REFPROP';%'StanMix'; %%% use RefProp in case of NACA0012
FLD.nCmp  = 1;
FLD.Cmp   = 'AIR';
FLD.Cnc   = [1 0];