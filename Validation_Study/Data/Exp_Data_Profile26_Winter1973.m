%% Experimental data Case_name,'ADA045367_7302_Winter1973_P26'
% Data from ADA045367 case 7302 Profile 26
% Original sources (used as much as possible): Winter & Gaudet (1973)
% Composed by D.Dijkshoorn (17-04-2020)

%% Description of case
% Objective was to find how to predict compressible
% boundary layers using incompressible prediction methods corrected with
% compressibility factors. This case is a nice comparison for MM, showing
% the relative compressibility effects inside the BL and in the core (edge)
% flow.
% Winter & Gaudet measured total pressure and total temperature seperately at same single
% station, just before nozzle outlet, at ZPG conditions. cf was also
% measured, T0 and Mach were varied to obtaind data for different unit
% Reynolds-numbers and Mach-numbers.

%% Data obtained from original source
% Mach-number distribution (pressure history) obtained from combining Fig. 6 with Fig. 1 in Winter & Gaudet (1973)
% profile data and integral properties obtained from original source:
% tables 1 (edge data and 'measured temperature' profile data), and 2c.
% The Mach-number obtained from the graph was scaled to make it reach Mach
% 2.8, etc.; The bypass air slots are said to be used above mach 2.2 (Ma > 2.2),
% and thus the simulation of profile 26 (Ma 2.8) might be the most accurate in terms
% of surface length simulated (the other cases might have a longer but unknown surface length, extending upstream).
% For Mach 0.2 they did not list the temperature measurements, since u/ue
% and rho/rhoe were indistuingishable from the tabulated values using
% r=0.89. Thus, this profile data is used here too.
% Furthermore, they say that the later simulations of mach 2.2 are the most
% reliable in their opinion. None of them was used here. I chose to use the
% data with increments of Mach-number for consistency.
% Since pressure and temperature were not measured at the same
% time/experiment, values of cf and P and T were interpolated on basis of
% unit Reynolds-number and location resepctively.

%% Data obtained from AGARD report
%
% Fernholz and Finley consider the temperature measurements as more
% consistent, and thus list the temperature measurements only (pressure
% measurements resulting in Mach-profiles were omitted by them). Measured
% temperatures resulted in listed density and velocity profiles, which were
% used to calculate the other properties. Their aim was to obtain a
% complete data set, I tried to use as much as possible the tabulated data
% from Winter & Gaudet, to avoid having their assumptions affect the
% comparison (validation).

% NB Fernholz & Finley assume a recovery factor. This is not useful for the
% current objective of validation, since it was not at all measured.

%% Doubts
% 

%% Data related to extraction of pressure history from drawings
EXP.SP = 480; % measured at station SP, also last station of Mach-number data

%% BL Edge and Wall data (table 1)
EXP.Pn = 26;                % Profile-number
EXP.MaE = 2.7996;           % [-], Mach-number
EXP.ReE = 1.997e6/0.3048;   % [1/m], Unit Reynolds-number
EXP.T0 = 18.7 + 273.15;     % [K]
EXP.cf = 121e-5;            % [-], local skin friction coefficient
EXP.S = 110.4;              % [K], Sutherland's constant for air
EXP.gamma = 1.4;            % [-], specific heat ratio
EXP.R = 287;                % [J/kg/K]

% Obtained from Fernholz & Finley (1977) (they assumed a recovery factor of
% 0.896):
% NB Estimated values!!! Not measured.
EXP.TE = 113.67;            % [K], BL Edge Temperature
EXP.TW = 272.25;            % [K], equals T_recovery (source: Agard)
EXP.P0 = 7.6454e4;          % [Pa], estimated by Fernholz & Finley

%% Profile data (table 1)
% Y-coordinate:
EXP.y = ... % [m]
     0.0254*[0.012,0.014,0.033,0.047,0.057,0.063,0.086,0.100,0.122,0.136,0.168,0.199,0.255,0.295,0.350,0.402,0.443,0.484,0.592,0.689,0.781,0.890,0.992,1.190,1.400,1.597,1.791,1.997,2.192,2.410,2.606,2.786,2.993,3.493,4.492,4.991,5.494];

% Profile data obtained from measured (static) temperature:
EXP.u_ue = ... % [-]
            [0.4864      0.4894	0.5770	0.5960	0.6111	0.6172	0.6452	0.6541	0.6713	0.6825	0.6955	0.7112	0.7320	0.7437	0.7545	0.767	0.7746	0.7815	0.7974	0.8095	0.8266	0.8345	0.8443	0.8651	0.8845	0.8996	0.9154	0.9301	0.9423	0.9556	0.9669	0.9757	0.9854	0.9855	0.9988	0.9998	1.0000];
EXP.rho_rhoe =  ... % [-]
            [0.4825	0.4835	0.5160	0.5240	0.5310	0.5332	0.5485	0.5538	0.5628	0.5700	0.5785	0.5885	0.6025	0.6124	0.6208	0.6312	0.6375	0.6432	0.6592	0.6707	0.6718	0.6965	0.7078	0.7330	0.7585	0.7816	0.8060	0.8303	0.8530	0.8770	0.9000	0.9175	0.9375	0.9376	1.0000	1.0000	1.0000];

% Profile data obtained from measured total temperature:
EXP.Ma = ... % [-]
            [0.9459,0.9526,1.1603,1.2078,1.2467,1.2617,1.3378,1.3628,1.4098,1.44250,1.4809,1.5274,1.5907,1.6292,1.6643,1.7061,1.7315,1.7547,1.8126,1.856,1.8968,1.9497,1.9885,2.0735,2.1566,2.2266,2.3008,2.3727,2.4364,2.5052,2.5681,2.6165,2.6712,2.6716,2.7963,2.7989,2.7996];

%% Integral properties (table 2c)
EXP.delta_ast = 0.0254*1.0340; % [m], displacement thickness (delta1)
EXP.theta = 0.0254*0.2271; % [m], displacement thickness (delta2)
EXP.delta3 = 0.0254*0.4178; % [m], displacement thickness
EXP.Delta = 0.0254*-0.0091; % [m], enthalpy thickness: integrate this term over BL: (rho u)/(rhoE uE)*(T0/T0E - 1)
EXP.H12 = 4.553; % [-]
EXP.H32 = 1.840; % [-]

% EXP.MaE = 2.7996; % [-], Mach-number, same as table 1
% EXP.ReE = 1.664e5/0.0254; % [1/m], Unit Reynolds-number, same as in table
% 1, but different units!
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% FOR REFERENCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Experimental data Case_name,'ADA045367_7302_Winter1973_P26'
% % Data from ADA045367 case 7302 Profile 26
% 
% % Velocity profile at last station (SP)
% EXP.SP = 480;               % (also) total amount of stations (obtained from figure 6 using grabit (MATLAB m-file))
% % (see Excel-file):
% EXP.rhoE = 0.086408;        % [kg/m3], obtained from IG: rho=p/RT, where p and T are measured at BL edge
% EXP.TE = 113.67;            % [K]
% EXP.TW = 272.25;            % [K], equals T_recovery (source: Agard)
% EXP.Cf = 121e-5;            % [-], skin friction coefficient
% % Load data case Winter 1973 profile 26, measured temperature:
% % (Note: zero wall velocity/wall coordinate is added. At the wall only the temperature
% % was measured, and added calculated from ADA data (see Excel-file)
% % note that all density values have a zero valued 5th digit! -> but temperature
% % is measured and 5 digit? U/UD is also 4 digit.
% EXP.y =         [0      0.0003048	0.0003556	0.0008382	0.0011938	0.0014478	0.0016002	0.0021844	0.00254	0.0030988	0.0034544	0.0042672	0.0050546	0.006477	0.007493	0.00889	0.0102108	0.0112522	0.0122936	0.0150368	0.0175006	0.0198374	0.022606	0.0251968	0.030226	0.03556	0.0405638	0.0454914	0.0507238	0.0556768	0.061214	0.0661924	0.0707644	0.0760222	0.0887222	0.1140968	0.1267714	0.1395476];
% EXP.u_ue =      [0      0.4864      0.4894	0.5770	0.5960	0.6111	0.6172	0.6452	0.6541	0.6713	0.6825	0.6955	0.7112	0.7320	0.7437	0.7545	0.767	0.7746	0.7815	0.7974	0.8095	0.8266	0.8345	0.8443	0.8651	0.8845	0.8996	0.9154	0.9301	0.9423	0.9556	0.9669	0.9757	0.9854	0.9855	0.9988	0.9998	1.0000];
% EXP.rho_rhoe =  [0.4175 0.4825	0.4835	0.5160	0.5240	0.5310	0.5332	0.5485	0.5538	0.5628	0.5700	0.5785	0.5885	0.6025	0.6124	0.6208	0.6312	0.6375	0.6432	0.6592	0.6707	0.6718	0.6965	0.7078	0.7330	0.7585	0.7816	0.8060	0.8303	0.8530	0.8770	0.9000	0.9175	0.9375	0.9376	1.0000	1.0000	1.0000];
% 
% %% Experimental data Case_name,'ADA045367_7302_Winter1973_P26'
% % Data from ADA045367 case 7302
% % velocity profile at last station (SP)
% SP = 480;               % (also) total amount of stations (obtained from figure 6 using grabit (MATLAB m-file))
% % (see Excel-file):
% rhoE = 0.086408;        % [kg/m3], obtained from IG: rho=p/RT, where p and T are measured at BL edge
% TE = 113.67;            % [K]
% TW = 272.25;            % [K], equals T_recovery (source: Agard)
% Cf = 121e-5;            % [-], skin friction coefficient
% % Load data case Winter 1973 profile 26, measured temperature:
% % (Note: zero wall velocity/wall coordinate is added. At the wall only the temperature
% % was measured, and added calculated from ADA data (see Excel-file)
% % note that all density values have a zero valued 5th digit! -> but temperature
% % is measured and 5 digit? U/UD is also 4 digit.
% y =         [0.0003048	0.0003556	0.0008382	0.0011938	0.0014478	0.0016002	0.0021844	0.00254	0.0030988	0.0034544	0.0042672	0.0050546	0.006477	0.007493	0.00889	0.0102108	0.0112522	0.0122936	0.0150368	0.0175006	0.0198374	0.022606	0.0251968	0.030226	0.03556	0.0405638	0.0454914	0.0507238	0.0556768	0.061214	0.0661924	0.0707644	0.0760222	0.0887222	0.1140968	0.1267714	0.1395476];
% u_ue =      [0.4864      0.4894	0.5770	0.5960	0.6111	0.6172	0.6452	0.6541	0.6713	0.6825	0.6955	0.7112	0.7320	0.7437	0.7545	0.767	0.7746	0.7815	0.7974	0.8095	0.8266	0.8345	0.8443	0.8651	0.8845	0.8996	0.9154	0.9301	0.9423	0.9556	0.9669	0.9757	0.9854	0.9855	0.9988	0.9998	1.0000];
% rho_rhoe =  [0.4825	0.4835	0.5160	0.5240	0.5310	0.5332	0.5485	0.5538	0.5628	0.5700	0.5785	0.5885	0.6025	0.6124	0.6208	0.6312	0.6375	0.6432	0.6592	0.6707	0.6718	0.6965	0.7078	0.7330	0.7585	0.7816	0.8060	0.8303	0.8530	0.8770	0.9000	0.9175	0.9375	0.9376	1.0000	1.0000	1.0000];