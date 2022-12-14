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

Case_name = 'ADA045367_7302_Winter1973_P12'; % tag: document name, case entry, first author and year

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
FLD.SLV = 110.4;                % [K],  Same as Winter & Gaudet (1974), (Groot, 2018) Sutherland's constant (later ones: White (2006) and Cebeci (2002)!)
FLD.SLC = 194;                  % [K], (Groot, 2018; and White, 2006) Sutherland's constant for thermal conductivity
FLD.Tref_S = 273;               % [K], (Groot, 2018; and White, 2006) reference temperature Sutherland's Law
FLD.mu_ref = 1.716e-5;          % [Pas], (Groot, 2018; and White, 2006) reference viscosity Sutherland's Law
FLD.k_ref = 0.0241;             % [W/m/K], (Groot, 2018; and White, 2006) reference thermal conductivity Sutherland's Law
FLD.Pr = 0.72;                  % [-], Winter & Gaudet?
FLD.Cp = FLD.Rsg*FLD.gamma/(FLD.gamma - 1); % [J/kg/K]
% FLD.Cpcoeff = [1022.5294853 -0.1758625 4.020605e-4 -4.8640623e-8]; % a0, a1, a2, a3, a4, etc. -> Cp = a0 + a1*T + ... % Andrews (1981) IG relation; range: 100 - 590 K

%% Case Specific Input Data

% Geometry (scaled)
INP.L = 1;                      % [m], scale factor: chord length/plate length/nozzle throat half-width; equal to one for non-scaled geometry
INP.x = ...                     % [-]
       [0.0000	0.5000	1.0000	1.7915	2.4732	3.1300	3.8034	4.4685	5.2001	5.8961	6.6279	7.3117	7.9132	8.4584	9.0003	9.4741	9.9495	10.3326	10.7358	11.0181	11.2600	11.5019	11.7439	11.9656	12.1777	12.3896	12.5582	12.7116	12.8063	12.8584	12.9267	12.9837	13.0357	13.0846	13.1431	13.1952	13.2505	13.3026	13.3612	13.4295	13.4751	13.5500	13.5972	13.6802	13.7518	13.8169	13.8788	13.9357	14.0106	14.0676	14.1148	14.1604	14.2109	14.2548	14.2971	14.3395	14.3785	14.4257	14.4713	14.5185	14.5722	14.6162	14.6732	14.7220	14.7709	14.8002	14.8376	14.8832	14.9369	14.9858	15.0411	15.0949	15.1535	15.2121	15.2626	15.3163	15.3700	15.4254	15.4791	15.5345	15.5947	15.6615	15.7071	15.7576	15.8097	15.8651	15.9188	15.9725	16.0312	16.0882	16.1338	16.1957	16.2462	16.3032	16.3651	16.4270	16.4889	16.5427	16.5932	16.6486	16.7007	16.7496	16.7887	16.8213	16.8539	16.9028	16.9517	16.9973	17.0478	17.0934	17.1407	17.1880	17.2303	17.2792	17.3330	17.3835	17.4373	17.4764	17.5285	17.5856	17.6426	17.6996	17.7567	17.8235	17.8838	17.9426	18.0029	18.0583	18.1187	18.1741	18.2279	18.2785	18.3127	18.3437	18.3763	18.4219	18.4758	18.5214	18.5687	18.6111	18.6601	18.7172	18.7710	18.8232	18.8722	18.9244	18.9701	19.0206	19.0631	19.1104	19.1659	19.2311	19.2850	19.3405	19.3976	19.4433	19.4890	19.5412	19.5885	19.6375	19.6800	19.7273	19.7697	19.8138	19.8710	19.9200	19.9722	20.0179	20.0636	20.1158	20.1697	20.2285	20.2889	20.3298	20.3869	20.4326	20.4735	20.4996	20.5420	20.5829	20.6253	20.6613	20.7103	20.7511	20.7952	20.8393	20.8802	20.9227	20.9668	21.0109	21.0484	21.0876	21.1302	21.1743	21.2282	21.2772	21.3213	21.3655	21.4178	21.4767	21.5404	21.5944	21.6418	21.6859	21.7317	21.7807	21.8379	21.8788	21.9262	21.9785	22.0358	22.0702	22.1078	22.1372	22.1634	22.2027	22.2354	22.2616	22.2927	22.3320	22.3696	22.4104	22.4448	22.4825	22.5283	22.5643	22.6035	22.6510	22.6935	22.7360	22.7753	22.8113	22.8587	22.9013	22.9471	22.9880	23.0338	23.0813	23.1336	23.1811	23.2383	23.2907	23.3348	23.3806	23.4248	23.4705	23.5147	23.5491	23.5916	23.6358	23.6799	23.7175	23.7584	23.8123	23.8499	23.8875	23.9235	23.9611	23.9939	24.0381	24.0952	24.1345	24.1835	24.2179	24.2539	24.2899	24.3357	24.3766	24.4175	24.4648	24.5040	24.5334	24.5760	24.6217	24.6610	24.7116	24.7606	24.8211	24.8734	24.9224	24.9861	25.0334	25.0906	25.1233	25.1478	25.1788	25.2294	25.2751	25.3208	25.3551	25.3895	25.4417	25.4792	25.5135	25.5559	25.6065	25.6457	25.6833	25.7241	25.7763	25.8138	25.8611	25.9101	25.9688	26.0226	26.0780	26.1204	26.1873	26.2460	26.3145	26.3830	26.4417	26.5020	26.5639	26.6291	26.6991	26.7562	26.8311	26.9061	26.9761	27.0543	27.1227	27.1927	27.2627	27.3327	27.4206	27.4890	27.5704	27.6533	27.7282	27.8144	27.9038	27.9722	28.0616	28.1559	28.2242	28.3006	28.3689	28.4437	28.5153	28.6047	28.6843	28.7510	28.8046	28.8680	28.9477	29.0306	29.0988	29.1898	29.2564	29.3312	29.4075	29.4774	29.5505	29.6236	29.7195	29.7990	29.8770	29.9566	30.0362	30.1256	30.2247	30.3270	30.4440	30.5659	30.6633	30.7624	30.8583	30.9460	31.0403	31.1264	31.2336	31.3148	31.3977	31.5017	31.5894	31.6918	31.7957	31.8769	31.9582	32.0427	32.1271	32.2019	32.2863	32.3773	32.4537	32.5284	32.6064	32.6925	32.7737	32.8485	32.9297	33.0109	33.0922	33.1425	33.193];
INP.y = zeros(1,length(INP.x)); % [-]

%% Free-Stream and BL Edge conditions
INP.PtI = 0.43688e5;            % [Pa], checked, correct, only value obtained from Fernholz & Finley (1977)
INP.TtI = 18.6 + 273.15;        % [K], 291.75;
% INP.MaI = 0.0702;               % [-]
INP.MaE = ...                   % [-]
       1.4/1.4017*[0.0702	0.0702	0.0703	0.0706	0.0707	0.0709	0.0712	0.0714	0.0716	0.0719	0.0723	0.0730	0.0738	0.0750	0.0765	0.0784	0.0810	0.0839	0.0870	0.0904	0.0933	0.0970	0.1008	0.1052	0.1095	0.1138	0.1180	0.1221	0.1250	0.1269	0.1292	0.1308	0.1326	0.1340	0.1356	0.1375	0.1397	0.1414	0.1436	0.1458	0.1475	0.1500	0.1518	0.1548	0.1573	0.1598	0.1625	0.1648	0.1679	0.1704	0.1727	0.1749	0.1766	0.1780	0.1798	0.1818	0.1831	0.1851	0.1870	0.1890	0.1908	0.1929	0.1957	0.1976	0.1996	0.2009	0.2025	0.2045	0.2070	0.2093	0.2116	0.2142	0.2168	0.2190	0.2216	0.2239	0.2266	0.2291	0.2314	0.2340	0.2365	0.2394	0.2419	0.2446	0.2474	0.2501	0.2530	0.2552	0.2588	0.2615	0.2641	0.2668	0.2698	0.2732	0.2767	0.2800	0.2831	0.2865	0.2894	0.2925	0.2953	0.2979	0.3002	0.3024	0.3046	0.3073	0.3108	0.3138	0.3171	0.3197	0.3228	0.3257	0.3283	0.3313	0.3344	0.3375	0.3411	0.3435	0.3464	0.3501	0.3540	0.3572	0.3617	0.3660	0.3698	0.3745	0.3791	0.3829	0.3876	0.3915	0.3956	0.4000	0.4019	0.4037	0.4067	0.4100	0.4140	0.4180	0.4219	0.4250	0.4291	0.4338	0.4384	0.4429	0.4474	0.4515	0.4556	0.4603	0.4638	0.4681	0.4733	0.4791	0.4847	0.4892	0.4939	0.4985	0.5024	0.5071	0.5114	0.5163	0.5209	0.5254	0.5295	0.5337	0.5393	0.5446	0.5496	0.5537	0.5584	0.5637	0.5692	0.5752	0.5812	0.5859	0.5917	0.5963	0.6006	0.6033	0.6078	0.6123	0.6168	0.6212	0.6263	0.6305	0.6351	0.6398	0.6450	0.6501	0.6548	0.6598	0.6638	0.6683	0.6741	0.6793	0.6854	0.6905	0.6961	0.7022	0.7085	0.7161	0.7241	0.7309	0.7381	0.7431	0.7489	0.7549	0.7619	0.7675	0.7735	0.7799	0.7875	0.7931	0.7982	0.8021	0.8063	0.8111	0.8153	0.8199	0.8252	0.8307	0.8359	0.8408	0.8458	0.8522	0.8581	0.8631	0.8682	0.8752	0.8818	0.8874	0.8928	0.8984	0.9048	0.9115	0.9175	0.9236	0.9307	0.9371	0.9442	0.9512	0.9584	0.9660	0.9722	0.9789	0.9844	0.9896	0.9952	1.0009	1.0067	1.0130	1.0189	1.0239	1.0294	1.0351	1.0400	1.0449	1.0500	1.0549	1.0611	1.0669	1.0734	1.0792	1.0846	1.0899	1.0948	1.1003	1.1062	1.1120	1.1175	1.1225	1.1263	1.1303	1.1358	1.1411	1.1466	1.1521	1.1581	1.1647	1.1708	1.1768	1.1833	1.1882	1.1940	1.1978	1.2006	1.2041	1.2088	1.2140	1.2180	1.2221	1.2260	1.2308	1.2341	1.2377	1.2416	1.2469	1.2508	1.2542	1.2578	1.2624	1.2662	1.2700	1.2740	1.2789	1.2828	1.2865	1.2901	1.2949	1.3002	1.3051	1.3100	1.3141	1.3181	1.3216	1.3258	1.3301	1.3337	1.3378	1.3423	1.3458	1.3498	1.3529	1.3565	1.3597	1.3622	1.3661	1.3691	1.3722	1.3745	1.3767	1.3788	1.3810	1.3828	1.3849	1.3870	1.3879	1.3894	1.3905	1.3915	1.3929	1.3939	1.3949	1.3960	1.3966	1.3975	1.3985	1.3990	1.3995	1.4000	1.3999	1.3999	1.4001	1.4001	1.4000	1.4000	1.4003	1.4001	1.4001	1.4001	1.4003	1.4003	1.4002	1.4001	1.4004	1.4003	1.4002	1.4000	1.4003	1.4006	1.4006	1.4006	1.4009	1.4010	1.4011	1.4014	1.4013	1.4013	1.4015	1.4011	1.4012	1.4013	1.4013	1.4013	1.4013	1.4014	1.4014	1.4014	1.4014	1.4014	1.4014	1.4016	1.4016	1.4016	1.4016	1.4017	1.4017];
INP.MaI = INP.MaE(1);           % [-]

% Wall Boundary Condition (for specific wall enthalpy ratio, derivative of wall enthalpy ratio, wall
% temperature, or wall heat flux (other then adiabatic flow)
if OPT.BCEE > 0
%     INP.BCW = 1.00001*ones(1,length(INP.x));   % [-], scaled h ratio or it's derivative
%     INP.BCW = 0.95*ones(1,length(INP.x));   % [-], scaled h ratio or it's derivative
    % optional %%% Check!!! %%%
    INP.AWD = []; % only g_aw or T_aw ((static) ratio? -> see OUTPUT) should be given here, calculated after adiabatic simulation
%     INP.AWD = 272.25*ones(1,length(INP.x));
else
    INP.BCW = zeros(1,length(INP.x));   % [-], zero wall heat flux: derivative at wall of enthalpy ratio equal to zero
end

%% Experimental data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% include in specific PLOT-file for generating thesis graphs!
% cd('./Data')
% Exp_Data_Profile12_Winter1973 % m-file, includes description of data
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
SET.NTR = 246;                    % Force transition at station NTR

% Grid Parameters
GRD.etaE = 8.0;                  % non-dimensional BL edge location
GRD.VGP = 1.14; % 1;1.02; for smooth BL velocity thickness      % variable grid parameter (spacing ratio)
GRD.Deta = 0.0001; % original: 0.01;    % first grid spacing: h1

%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the following profiles (velocity, dimensionless velocity, enthalpy)
PLT.NSP = [length(INP.x)];               % Plot profiles at stations NS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%