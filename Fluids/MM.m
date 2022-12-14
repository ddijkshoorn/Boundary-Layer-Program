%% MM Fluid property data

% Molecular Fluid properties
FLD.R = 8.3145;                         % [J/mol/K], Universal gas constant
FLD.MW = 162.3768;                      % [g/mol = kg/kmol], Molecular Weight MM ('RefProp')
FLD.Rsg = FLD.R/FLD.MW*1e3;             % [J/kg/K], Specific Gas constant Air

% Fluid Property data
FLD.Pcrit = 19.3113e5;                  % [Pa], 'REFPROP' 1931.13 !Critical pressure [kPa]
FLD.Tcrit = 518.7;                      % [K], MM 'REFPROP'
FLD.vcrit = 0.003705922873401;          % [m3/kg], MM 'REFPROP': 0.003334265686713;
FLD.gamma = 1.0118;                     % [-], Low P=0.01 bar, high T=300 Celsius (IG-limit) Specific heat ratio (Cp/Cv); MM 'REFPROP': 1.011789571824677 with T=300 celsius, P=0.01 bar.
FLD.SLV = 393.2584;                     % [K], obtained from RefProp data with P=0.1 bar; T range: 311.15 - 541.15 K
FLD.SLC = 236322660091;                     	% [K]
FLD.Tref_S = 350;                       % [K], well chosen, near lower temperature limit
FLD.mu_ref = 7.6500e-06;                % [kg/m/s], obtained from RefProp data with P=0.1 bar
FLD.k_ref = 1.69421e-2;               	% [W/m/K] 0.016942052178071
% FLD.Tref_SLC = ;            	% [K], separate Tref needed for different Thermal conductivity relation?
FLD.Pr = 0.751553006708536;             % [-], Prandlt-number (same IG-limit), StanMix; REFPROP estiamte: 0.75;
% FLD.Cp = 2176;                      % [J/kg/K], Isobaric specific heat RefProp; Constant Pressure Heat capacity StanMix (same IG-limit): 2174.393487565100;
FLD.Cp = FLD.Rsg*FLD.gamma/(FLD.gamma - 1); % [J/kg/K]
% Based on curve-fit (MATLAB polyfit-function with n=4, CpE,TsE from ORCHID sim)
% FLD.Cpcoeff = [2.515385539556594e+07 -2.043093883596354e+05 6.223887307944882e+02 -0.842733774976201 4.279599164975811e-04]; % n=4 %%% 	% a0, a1, a2, a3, a4, etc. -> Cp = a0 + a1*T + ... 
% FLD.Cpcoeff = [-1.981176058677229e+06 1.211866172834274e+04 -24.691732596913965 0.016774988253622 0]; % n=3 %%% 
%
% Based on IG (P=0.1 bar) only (same as dyn visc and therm conductivity (file)):
FLD.Cpcoeff = [1.342980882726127e+03 -4.525540680281105 0.027684932803028 -4.413152424374620e-05 2.448539196059856e-08]; % n=4 %%% 	% a0, a1, a2, a3, a4, etc. -> Cp = a0 + a1*T + ...
%
% FLD.lambda = ;				% [Pas], second coefficient of viscosity (dilatational viscosity)
% FLD.mub = ;                   % [Pas], bulk viscosity
% FLD.mubcoeff = [];				% b0, b1, b2, b3, b4, etc. -> mub = b0 + b1*T + ... % Cramer (2012)?

% FluidProp
FLD.Model = 'REFPROP'; %'StanMix';%
FLD.nCmp  = 1;
FLD.Cmp   = 'MM';
FLD.Cnc   = [1 0];