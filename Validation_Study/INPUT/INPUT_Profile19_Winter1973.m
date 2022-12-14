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
%   Validation of: turbulence model by turbulent BL characteristics

Case_name = 'ADA045367_7302_Winter1973_P19'; % tag: document name, case entry, first author and year

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
INP.x = ...                     % [-], rounded first (negative) X-coordinate to zero
       [0.0000	0.3000	0.6592	1.0313	1.4214	1.8415	2.3422	2.8051	3.2775	3.7757	4.2642	4.7135	5.1922	5.7784	6.3059	6.8139	7.3609	7.9178	8.4062	8.7677	9.0900	9.4222	9.7055	9.9790	10.2623	10.5366	10.8025	11.0568	11.3227	11.5770	11.8082	12.0162	12.2243	12.4439	12.5434	12.6003	12.6659	12.7556	12.8430	12.9065	12.9962	13.0815	13.1755	13.2893	13.3878	13.4731	13.5518	13.6306	13.7072	13.7969	13.8757	13.9523	14.0463	14.1382	14.2214	14.2958	14.3702	14.4424	14.5322	14.6197	14.6963	14.7663	14.8298	14.9020	14.9830	15.0705	15.1647	15.2544	15.3332	15.4295	15.5171	15.5959	15.6769	15.7447	15.8213	15.8980	15.9615	16.0228	16.0907	16.1695	16.2592	16.3337	16.4147	16.4870	16.5615	16.6359	16.7104	16.7761	16.8331	16.9053	16.9623	17.0411	17.1134	17.1747	17.2383	17.3062	17.3697	17.4377	17.4990	17.5647	17.6173	17.6568	17.7203	17.7861	17.8562	17.9285	17.9920	18.0578	18.1170	18.1762	18.2398	18.3099	18.3866	18.4502	18.5094	18.5840	18.6498	18.7112	18.7726	18.8340	18.8889	18.9525	19.0096	19.0688	19.1281	19.1786	19.2291	19.2752	19.3301	19.3674	19.4266	19.4903	19.5583	19.6066	19.6768	19.7295	19.7844	19.8436	19.9051	19.9424	19.9930	20.0325	20.0852	20.1270	20.1775	20.2302	20.2697	20.3159	20.3642	20.4016	20.4477	20.4894	20.5334	20.5905	20.6411	20.6916	20.7356	20.7774	20.8191	20.8631	20.9049	20.9488	20.9906	21.0368	21.0786	21.1182	21.1623	21.1997	21.2349	21.2767	21.3207	21.3603	21.3977	21.4439	21.4880	21.5276	21.5628	21.6025	21.6355	21.6751	21.7148	21.7478	21.7875	21.8227	21.8537	21.8801	21.9043	21.9220	21.9572	21.9859	22.0212	22.0477	22.0830	22.1183	22.1535	22.1910	22.2175	22.2593	22.2946	22.3300	22.3631	22.3962	22.4315	22.4690	22.5000	22.5353	22.5750	22.6038	22.6369	22.6612	22.6966	22.7340	22.7693	22.7981	22.8313	22.8600	22.8953	22.9241	22.9484	22.9793	23.0103	23.0390	23.0678	23.0944	23.1297	23.1650	23.1916	23.2313	23.2667	23.2911	23.3155	23.3486	23.3840	23.4106	23.4504	23.4683	23.5080	23.5390	23.5744	23.6011	23.6321	23.6632	23.6942	23.7273	23.7517	23.7718	23.8071	23.8272	23.8560	23.8841	23.9084	23.9457	23.9757	24.0000	24.0354	24.0634	24.0934	24.1178	24.1533	24.1888	24.2262	24.2580	24.2916	24.3216	24.3478	24.3851	24.3946	24.4170	24.4413	24.4675	24.4936	24.5217	24.5423	24.5741	24.5965	24.6319	24.6600	24.6843	24.7159	24.7347	24.7534	24.7815	24.8151	24.8449	24.8767	24.8954	24.9234	24.9532	24.9720	24.9944	25.0113	25.0373	25.0690	25.0914	25.1064	25.1250	25.1549	25.1810	25.2072	25.2389	25.2614	25.2838	25.3044	25.3399	25.3697	25.4033	25.4258	25.4557	25.4929	25.5283	25.5637	25.5954	25.6234	25.6551	25.6869	25.7167	25.7355	25.7598	25.7841	25.8009	25.8253	25.8421	25.8590	25.8832	25.9001	25.9170	25.9413	25.9563	25.9843	26.0012	26.0199	26.0462	26.0686	26.0874	26.1190	26.1544	26.1750	26.2011	26.2383	26.2607	26.2905	26.3148	26.3483	26.3652	26.3931	26.4285	26.4584	26.4975	26.5218	26.5571	26.5851	26.6075	26.6391	26.6708	26.7006	26.7397	26.7788	26.8178	26.8513	26.8886	26.9184	26.9685	27.0038	27.0372	27.0688	27.0967	27.1303	27.1711	27.2101	27.2491	27.2992	27.3308	27.3716	27.4143	27.4515	27.4941	27.5368	27.5758	27.6314	27.6760	27.7186	27.7557	27.8001	27.8335	27.8798	27.9262	27.9706	28.0188	28.0651	28.1077	28.1559	28.2077	28.2577	28.3225	28.3762	28.4354	28.4928	28.5612	28.6352	28.7147	28.7757	28.8645	28.9291	29.0049	29.0806	29.1730	29.2468	29.3225	29.3927	29.4684	29.5496	29.6309	29.7287	29.8080	29.8763	29.9686	30.0479	30.1364	30.2397	30.3522	30.4462	30.5477	30.6638	30.7579	30.8482	30.9294	31.0142	31.1027	31.2078	31.3148	31.4217	31.5286	31.6134	31.7185	31.8088	31.8972	32.0042	32.0742	32.1479	32.2438	32.3341	32.4152	32.5184	32.6142	32.7064	32.7911	32.8943	32.9699	33.0657	33.1302	33.193];

INP.y = zeros(1,length(INP.x)); % [-]

%% Free-Stream and BL Edge conditions
INP.PtI = 0.34522e5;            % [Pa], only value obtained from Fernholz & Finley (1977)
INP.TtI = 18.5 + 273.15;        % [K],  291.65
% INP.MaI = 0.0372;               % [-]
INP.MaE = ...                   % [-]
       2.2/2.1912*[0.0372	0.0372	0.0375	0.0376	0.0377	0.0379	0.0383	0.0387	0.0391	0.0398	0.0404	0.0410	0.0416	0.0424	0.0431	0.0439	0.0449	0.0457	0.0466	0.0475	0.0481	0.0490	0.0496	0.0505	0.0515	0.0526	0.0540	0.0550	0.0564	0.0581	0.0598	0.0611	0.0629	0.0649	0.0662	0.0669	0.0678	0.0689	0.0694	0.0703	0.0715	0.0724	0.0738	0.0756	0.0770	0.0782	0.0795	0.0809	0.0821	0.0840	0.0855	0.0870	0.0886	0.0904	0.0920	0.0934	0.0949	0.0967	0.0988	0.1009	0.1028	0.1042	0.1056	0.1072	0.1089	0.1109	0.1133	0.1154	0.1175	0.1200	0.1224	0.1245	0.1268	0.1289	0.1308	0.1334	0.1357	0.1372	0.1393	0.1416	0.1446	0.1473	0.1500	0.1531	0.1560	0.1585	0.1613	0.1639	0.1658	0.1688	0.1710	0.1740	0.1766	0.1792	0.1820	0.1851	0.1877	0.1910	0.1938	0.1961	0.1990	0.2011	0.2044	0.2075	0.2105	0.2140	0.2174	0.2209	0.2240	0.2279	0.2315	0.2353	0.2397	0.2431	0.2468	0.2516	0.2560	0.2603	0.2645	0.2683	0.2725	0.2771	0.2818	0.2865	0.2908	0.2953	0.3005	0.3052	0.3092	0.3127	0.3163	0.3208	0.3264	0.3305	0.3362	0.3402	0.3446	0.3491	0.3541	0.3583	0.3633	0.3669	0.3725	0.3768	0.3813	0.3858	0.3900	0.3950	0.3995	0.4036	0.4075	0.4120	0.4168	0.4222	0.4279	0.4339	0.4388	0.4435	0.4487	0.4532	0.4578	0.4629	0.4679	0.4736	0.4790	0.4845	0.4911	0.4964	0.5008	0.5056	0.5113	0.5163	0.5224	0.5285	0.5343	0.5406	0.5458	0.5515	0.5572	0.5620	0.5683	0.5736	0.5797	0.5852	0.5913	0.5963	0.5996	0.6037	0.6080	0.6136	0.6198	0.6255	0.6311	0.6378	0.6430	0.6496	0.6554	0.6611	0.6677	0.6745	0.6816	0.6878	0.6944	0.7016	0.7087	0.7158	0.7224	0.7296	0.7360	0.7429	0.7497	0.7563	0.7625	0.7692	0.7777	0.7841	0.7903	0.7971	0.8019	0.8080	0.8159	0.8230	0.8299	0.8372	0.8451	0.8517	0.8591	0.8661	0.8740	0.8809	0.8872	0.8951	0.9029	0.9096	0.9195	0.9273	0.9357	0.9439	0.9520	0.9603	0.9684	0.9783	0.9864	0.9936	1.0011	1.0085	1.0165	1.0251	1.0329	1.0412	1.0502	1.0590	1.0684	1.0764	1.0836	1.0917	1.1014	1.1113	1.1212	1.1317	1.1418	1.1520	1.1617	1.1716	1.1793	1.1885	1.1949	1.2007	1.2074	1.2161	1.2233	1.2314	1.2391	1.2477	1.2542	1.2632	1.2711	1.2788	1.2852	1.2920	1.2986	1.3066	1.3155	1.3236	1.3325	1.3391	1.3452	1.3525	1.3592	1.3655	1.3714	1.3770	1.3832	1.3895	1.3946	1.3998	1.4065	1.4137	1.4223	1.4294	1.4368	1.4437	1.4507	1.4596	1.4676	1.4763	1.4839	1.4915	1.4992	1.5076	1.5161	1.5231	1.5312	1.5393	1.5477	1.5551	1.5618	1.5693	1.5761	1.5815	1.5894	1.5952	1.6000	1.6054	1.6116	1.6186	1.6252	1.6313	1.6379	1.6444	1.6516	1.6601	1.6676	1.6738	1.6803	1.6877	1.6950	1.7009	1.7076	1.7141	1.7211	1.7273	1.7347	1.7417	1.7477	1.7549	1.7632	1.7716	1.7790	1.7854	1.7924	1.7980	1.8045	1.8123	1.8190	1.8255	1.8338	1.8421	1.8488	1.8563	1.8636	1.8716	1.8774	1.8824	1.8882	1.8939	1.9009	1.9076	1.9130	1.9204	1.9276	1.9335	1.9392	1.9453	1.9510	1.9573	1.9628	1.9704	1.9766	1.9834	1.9882	1.9933	1.9977	2.0025	2.0076	2.0142	2.0189	2.0240	2.0291	2.0342	2.0394	2.0445	2.0499	2.0559	2.0607	2.0667	2.0711	2.0764	2.0828	2.0888	2.0932	2.0995	2.1037	2.1077	2.1122	2.1169	2.1202	2.1233	2.1268	2.1305	2.1334	2.1375	2.1405	2.1431	2.1453	2.1483	2.1500	2.1517	2.1539	2.1561	2.1583	2.1605	2.1630	2.1652	2.1664	2.1678	2.1691	2.1715	2.1734	2.1749	2.1762	2.1777	2.1786	2.1797	2.1804	2.1812	2.1825	2.1830	2.1839	2.1846	2.1852	2.1858	2.1866	2.1874	2.1881	2.1886	2.1895	2.1899	2.1904	2.1912	2.1912];
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
% Exp_Data_Profile19_Winter1973 % m-file, includes description of data
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
SET.NPT = 201;  % 600;           % preallocation (size is 201 in FORTRAN code)
SET.ITMAX = 6;                   % maximum number of iterations (6 in MAIN code, 8 for IVPL in FORTRAN code)
SET.ITMAX0 = 10;
SET.NTR = 268;                   % Force transition at station NTR

% Grid Parameters
GRD.etaE = 8.0;                  % non-dimensional BL edge location
GRD.VGP = 1.14; % 1;1.02; for smooth BL velocity thickness      % variable grid parameter (spacing ratio)
GRD.Deta = 0.0001; % original: 0.01;    % first grid spacing: h1

%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the following profiles (velocity, dimensionless velocity, enthalpy)
PLT.NSP = [length(INP.x)];               % Plot profiles at stations NS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%