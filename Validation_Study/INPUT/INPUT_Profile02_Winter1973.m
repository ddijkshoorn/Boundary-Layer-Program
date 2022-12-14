%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   INPUT File                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Written by D.D.Dijkshoorn (DDD) on 17/04/2020                         %
%   Checked by ... on ././.                                               %
%   Last revision: 03/03/2022                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Case: Continuous running windtunnel with variable nozzle
%   Type: Flat Plate flow (zero-PG Plate flow)
%   Source: ADA045367 (1977) case entry 7302
%   Original source: Winter, K.G. and Gaudet, L.; 1973. Turbulent
%   boundary-layer studies at high Reynolds numbers at Mach numbers between
%   0.2 and 2.8. ARC (London) R&M 3712
%   Validation of: turbulence model by comparing velocity profile and ...

Case_name = 'ADA045367_7302_Winter1973_P02'; % tag: document name, case entry, first author and year

%% OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ideal Gas (IG); Nonideal Gas (NG), Real Gas (RG) (not applied here)

OPT.GASM = 1;                   % Gas model: 1=calorically perfect IG; 2=thermally perfect IG; 3=Nonideal Gas (NG);
%%% The following three options are only available with GASM=1, and in this
%%% order. Options are ignored when inappropriate, see manual. Note that
%%% constants for Sutherland's Law and the cp(T) polynomial are not
%%% available for all fluids
OPT.COMP = 1;                   % 0=incompressible flow; 1=compressible flow (standard); only available with GASM=1 and only valid for gases and incompressible input data; set FLD.C in fluid file
OPT.CPRN = 0;                   % 0=constant Pr-number (input) and k=f(Pr,mu,cp); 1=variable Pr-number: Pr=f(mu(Suth(T)),k(Suth(T)),cp=constant); only available with GASM=1
OPT.CCRP = 1;                   % 0=constant Chapman-Rubesin parameter, 1=general/variable Chapman-Rubesin parameter; only available with GASM=1
%%%
OPT.CPRT = 0;1;                   % 0=constant PrT; 1=variable turbulent Pr-number
OPT.BCIE = 2;                   % BC Input Edge (BL Edge input: 1=UE/uE; 2=MaE; 3=psE (ratio?))
OPT.BCEE = 0;                   % Wall BC of EE: 0=adiabatic wall; 1=enthalpy ratio [-]; 2=derivative of enthalpy ratio [-]; 3=Temperature ratio [-]; 4=heat flux [W/m2] (temperature or heat flux directly)
OPT.SMTH = 0;                   % Smoothing experimental data, entry is number of runs)
OPT.TRME = 1;                   % transition method: 0=no transition; 1=prescribed transition location NTR; 2=Wazzan correlation; 3=Michel's Method (Cebeci (2002)); 4=fully turbulent
OPT.RLAM = 0;                   % 1=Relaminarization check, 0=no check; simple engineering estimate based on exp. data for air
OPT.GRAD = 1;                   % 1=Lagrange; 2=Weighted-difference technique; 3=SPLINE
OPT.CHRT = 0;                   % Pv-diagram and Ts-diagram including range

%% PHYSICAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fluid possibilities: AIR, CO2, MDM2, D6, ...
% cd('../Fluids')
% AIR
% cd('../INPUT')

% Molecular Fluid properties
FLD.R = 8.3145;                 % [J/mol/K], Universal gas constant
FLD.MW = 28.97;                 % [g/mol = kg/kmol], Molecular Weight of dry Air
FLD.Rsg = FLD.R/FLD.MW*1e3;     % [J/kg/K], Specific Gas constant Air

% Fluid Property data
FLD.gamma = 1.4;                % [-], Specific heat ratio (Cp/Cv)
FLD.SLV = 110.4;                % [K],  (Groot, 2018) Sutherland's constant (later ones: White (2006) and Cebeci (2002)!)
FLD.SLC = 194;                  % [K], (Groot, 2018; and White, 2006) Sutherland's constant for thermal conductivity
FLD.Tref_S = 273;               % [K], (Groot, 2018; and White, 2006) reference temperature Sutherland's Law
FLD.mu_ref = 1.716e-5;          % [Pas], (Groot, 2018; and White, 2006) reference viscosity Sutherland's Law
FLD.k_ref = 0.0241;             % [W/m/K], (Groot, 2018; and White, 2006) reference thermal conductivity Sutherland's Law
FLD.Pr = 0.72;                  % [-], (Cebeci, 2002) Prandlt-number
FLD.Cp = FLD.Rsg*FLD.gamma/(FLD.gamma - 1); % [J/kg/K]
% FLD.Cpcoeff = [1022.5294853 -0.1758625 4.020605e-4 -4.8640623e-8]; % a0, a1, a2, a3, a4, etc. -> Cp = a0 + a1*T + ... % Andrews (1981) IG relation; range: 100 - 590 K

%% Case Specific Input Data

% Geometry (scaled)
INP.L = 1;                      % [m], scale factor: chord length/plate length/nozzle throat half-width; equal to one for non-scaled geometry
INP.x = ...                     % [-]
       [0	0.5	0.92046007	1.643064237	2.382873265	3.208706598	3.982925349	4.877578127	5.686206599	6.426015628	7.157663753	7.858022242	8.530366391	9.20271054	9.81902601	10.35129846	10.82754223	11.21974299	11.49988638	11.80804412	12.10297474	12.3629044	12.52536044	12.66021677	12.75512078	12.84791726	12.93438298	13.01031112	13.09676862	13.18113504	13.26126999	13.36250477	13.43209392	13.51856785	13.60714109	13.69361502	13.77587392	13.87078614	13.95302862	14.04583332	14.12174502	14.17869728	14.24194745	14.30099901	14.36638956	14.43387119	14.50980754	14.58995071	14.67640821	14.75867533	14.83670276	14.92739995	15.00965885	15.08347125	15.17415201	15.24374116	15.31966929	15.39557278	15.47783168	15.56009058	15.65288707	15.7414603	15.81948774	15.90173843	15.99453491	16.07679381	16.16746635	16.24973347	16.3277609	16.40789586	16.48593151	16.56187607	16.64202745	16.72218705	16.79182549	16.87410082	16.95636794	17.03020499	17.09559553	17.18419342	17.25591473	17.32129706	17.41200246	17.49003811	17.58075173	17.64614227	17.71997932	17.78328699	17.86133907	17.94783765	18.01112068	18.09340423	18.1567119	18.22844965	18.30862568	18.39931465	18.46891201	18.57222144	18.64604206	18.7303674	18.81685776	18.90122418	18.98560703	19.07208097	19.17121644	19.25138425	19.31255155	19.37792566	19.45173806	19.53819557	19.62045447	19.68161355	19.75966563	19.84192453	19.93052241	20.02751751	20.13085158	20.21944125	20.3080227	20.41979506	20.5125669	20.60113192	20.67495254	20.7698319	20.87106667	20.96173921	21.05240354	21.14940684	21.2421869	21.32652867	21.39821713	21.49309649	21.59220731	21.68709489	21.78406533	21.88947406	21.98013017	22.04971111	22.12564746	22.20576598	22.29640566	22.37231736	22.4629899	22.53676944	22.60632573	22.69697363	22.76865386	22.84455735	22.93730455	23.02372098	23.12492289	23.2197776	23.31463232	23.42002462	23.5191108	23.62238737	23.71515099	23.81422074	23.90272826	24.01658349	24.10509922	24.21469013	24.32009886	24.41072211	24.47606336	24.55826476	24.64257367	24.72475864	24.81960513	24.91445163	25.03036509	25.14628676	25.24956333	25.35705494	25.47508413	25.58256752	25.67951332	25.78488919	25.87340492	25.96401174	26.06517257	26.16633341	26.26117991	26.34125736	26.44031889	26.5246278	26.61945787	26.70799003	26.78597639	26.8576402	26.91453496	26.97566118	27.05151537	27.14422971	27.23274544	27.33179876	27.42243022	27.51514456	27.61631361	27.71115189	27.80177514	27.90293598	27.99567496	28.07786814	28.16426814	28.25277566	28.35604401	28.44877478	28.53097617	28.62792197	28.71432197	28.7901926	28.87870833	28.97147196	29.06209521	29.15693349	29.24968069	29.3234438	29.38667753	29.45833313	29.51734361	29.57846162	29.65853907	29.73019466	29.81451178	29.90724255	29.99152682	30.08004255	30.1791123	30.29713327	30.39197156	30.50576928	30.64697295	30.77342399	30.937835	31.05795528	31.17177765	31.275046	31.36567747	31.46892939	31.57432169	31.69023515	31.7935035	31.90311905	31.99162657	32.08436555	32.18973321	32.30564666	32.42577516	32.55221799	32.67233827	32.80931868	32.94000119	33.04537706	33.1254463	33.193];
INP.y = zeros(1,length(INP.x)); % [-]

%% Free-Stream and BL Edge conditions
INP.PtI = 1.4688e5;             % [Pa], only value obtained from Fernholz & Finley (1977)
INP.TtI = 10.8 + 273.15;        % [K]
% INP.MaI = 0.0247;               % [-]
INP.MaE = ...                   % [-]
       0.2/0.203502374*[0.0247	0.0247	0.024761441	0.024833416	0.024948575	0.025063734	0.025207683	0.025394817	0.025553161	0.025769084	0.026029457	0.026384291	0.026881059	0.027661695	0.028797166	0.029932636	0.031351974	0.032629378	0.03419065	0.035609988	0.03769505	0.040063161	0.042431271	0.043947875	0.045309176	0.046669818	0.047862026	0.049050943	0.050076692	0.0514347	0.052458474	0.053988209	0.054842232	0.056200899	0.057393766	0.058752433	0.060109783	0.061637543	0.062661976	0.064189077	0.065045075	0.066061606	0.066580736	0.067431466	0.068617091	0.069470455	0.07082583	0.072016063	0.073041812	0.074565622	0.075588737	0.07711518	0.078472531	0.079494329	0.080687854	0.081541877	0.082730793	0.083420332	0.084777682	0.086135033	0.087495675	0.088688542	0.089711657	0.090902548	0.092263191	0.093620541	0.094647607	0.096171417	0.097194532	0.098218306	0.09940788	0.100929715	0.102286407	0.103809558	0.105662336	0.107352604	0.108876414	0.11039759	0.111583214	0.113275458	0.114629517	0.115648682	0.117341584	0.118531159	0.12039052	0.121576144	0.12309732	0.124781663	0.126304156	0.128162201	0.129347167	0.131203895	0.132888238	0.134575215	0.136431284	0.137791268	0.13881175	0.139676307	0.140864565	0.141390278	0.143081863	0.144439872	0.146130799	0.147489467	0.149185003	0.150874613	0.151892461	0.152745167	0.153766966	0.154792715	0.156150065	0.157001454	0.158523947	0.159881298	0.161573541	0.162602582	0.163966517	0.165492302	0.166851627	0.168384654	0.169245919	0.170272327	0.171460584	0.172322508	0.173852243	0.174879309	0.175739915	0.176935415	0.17796314	0.178821771	0.179509993	0.180371916	0.181568075	0.182596457	0.183126121	0.183824877	0.184519024	0.185206588	0.186561963	0.187252819	0.187614048	0.188470046	0.189497112	0.189853074	0.19004126	0.190568948	0.191090711	0.19178025	0.192142137	0.19233559	0.193199489	0.193562035	0.193924581	0.194290419	0.1949872	0.19518592	0.195880726	0.196244589	0.196105782	0.19714009	0.197167743	0.197201979	0.197900735	0.197929046	0.198115915	0.198308052	0.198500847	0.198360065	0.198556152	0.198752239	0.19878845	0.198991121	0.199189841	0.199389878	0.199593208	0.199626786	0.199657072	0.199689991	0.199717644	0.199413036	0.199444639	0.199476242	0.199672329	0.199530889	0.199728292	0.199921087	0.199784256	0.200144826	0.200335646	0.200524491	0.200375808	0.200561361	0.200252144	0.199948195	0.199975848	0.200006792	0.200201562	0.199897613	0.200095675	0.200125303	0.200153614	0.200185217	0.200380645	0.200406322	0.200266857	0.200128051	0.200160312	0.200189281	0.200381418	0.200411704	0.200272239	0.200295941	0.200323593	0.201018399	0.20104671	0.201076338	0.201438225	0.201461269	0.20164748	0.201669865	0.2016883	0.201707394	0.201565953	0.201588339	0.201947593	0.201976562	0.201669979	0.201697632	0.202061495	0.202098365	0.202127992	0.201997086	0.202041199	0.202080702	0.202631434	0.202502504	0.202870975	0.202903236	0.203098006	0.202797349	0.203163187	0.203199399	0.20323166	0.203765274	0.203626467	0.203821896	0.203688356	0.203724568	0.203762096	0.20363514	0.20350621	0.203382546	0.203756285	0.203789204	0.203481305	0.203502374];
INP.MaI = INP.MaE(1);           % [-]

% Wall Boundary Condition (for specific wall enthalpy ratio, derivative of wall enthalpy ratio, wall
% temperature, or wall heat flux (other then adiabatic flow)
if OPT.BCEE > 0
    INP.BCW = 0.95*ones(1,length(INP.x));   % [-], scaled h ratio or it's derivative
    % optional %%% Check!!! %%%
    INP.AWD = []; % only g_aw or T_aw ((static) ratio? -> see OUTPUT) should be given here, calculated after adiabatic simulation
%     INP.AWD = 272.25*ones(1,length(INP.x));
else
    INP.BCW = zeros(1,length(INP.x));   % [-], zero wall heat flux: derivative at wall of enthalpy ratio equal to zero
end

%% Experimental data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% include in specific PLOT-file for generating thesis graphs!
% cd('./Data')
% Exp_Data_Profile02_Winter1973 % m-file, includes description of data
% cd('..')

%% TURBULENCE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Closure coeficients, for the CS-model: EV (Cebeci, 2002), and PrT (Cebeci, 1974)

% Momentum Equation; and,
TCC.kappa = 0.40;           % [-], Von Karman constant
TCC.A_plus = 26.0;          % [-], Van Driest damping factor constant
TCC.alpha = 0.0168;         % [-], Clauser(’s)/outer eddy viscosity constant
TCC.ints = 11.8;            % [-], Assumed intersection of viscous sublayer with (intermediate) log layer

% Energy Equation:
TCC.kappa_h = 0.44;         % [-], Heat transfer mixing-length constant
TCC.PrT = 0.9;              % [-], (Cebeci (2004)!) Constant Turbulent Prandtl-number (independent of fluid!) if PrT model is not used
TCC.Bcoeff = [34.96 28.79 33.95 6.33 -1.186]; % [-], (Cebeci (1974)!) Turbulent Prandtl-number coefficients

%% NUMERICAL PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set size of vectors for initilization (vertical/eta-direction)
SET.NPT = 201;  % 600;            % preallocation (size is 201 in FORTRAN code)
SET.ITMAX = 6;                   % maximum number of iterations (6 in MAIN code, 8 for IVPL in FORTRAN code)
SET.ITMAX0 = 10;
SET.NTR = 110;                    % Force transition at station NTR

% Grid Parameters
GRD.etaE = 8.0;                  % non-dimensional BL edge location
GRD.VGP = 1.14; % 1;1.02; for smooth BL velocity thickness      % variable grid parameter (spacing ratio)
GRD.Deta = 0.0001; % original: 0.01;    % first grid spacing: h1

%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the following profiles (velocity, dimensionless velocity, enthalpy)
PLT.NSP = [length(INP.x)];               % Plot profiles at stations NS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%