%% ORCHID Nozzle Simulations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% With this file all 7 INPUT files can be run again 

clear all
close all
clc

pause(0.1)

cd ./..

%% Nonideal gas: Interpolated pressure PsE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run('Simulation/INPUT/INPUT_ORCHID_Nozzle_laminar_RefProp_PsE_Interpolated');
% run('Simulation/INPUT/INPUT_ORCHID_Nozzle_turbulent_RefProp_PsE_Interpolated');
% run('Simulation/INPUT/INPUT_ORCHID_Nozzle_laminar_StanMix_PsE_Interpolated');
% run('Simulation/INPUT/INPUT_ORCHID_Nozzle_turbulent_StanMix_PsE_Interpolated');

%% IG simulation cases %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run('Simulation/INPUT/INPUT_ORCHID_Nozzle_laminar_PsE_Interpolated_IG_Cal_Perf');
% run('Simulation/INPUT/INPUT_ORCHID_Nozzle_laminar_PsE_Interpolated_IG_Thermally_Perf');
% run('Simulation/INPUT/INPUT_ORCHID_Nozzle_laminar_PsE_Interp_IG_Th_Perf_CpE_poly');
%%


% PRECAL (calculate Boundary Conditions)
[X,NST,HVR,EDG,FRS,FLD,TCC,SOL,BLC,MON] = PRECAL(INP,OPT,SET,FLD,TCC);

% Grid generation (generate grid)
[NP,GRD] = GRID(GRD,SET);

NS = 1;

% IVPL (calculate initial values)
[sol,solprev] = IVPL(NS,GRD,HVR,FLD,SET);

tic

% solve PDE's with the Cebeci-Smith-method (CS-method = CSM)
[BLC,FLP,SOL,TCC,MON] = CSM(X,NP,NS,NST,GRD,HVR,EDG,FRS,FLD,TCC,sol,solprev,SOL,BLC,INP,OPT,SET,MON);

toc

%% Plot Results and generate graphs and tables

fprintf('Simulation finished with success \n')

keyboard

STORDATA
PLOTFILE
% TABLEGEN
% PLOTCHART % Pv-, and Ts-diagram, calculate critical properties inside
if OPT.CHRT > 0
    FPCHARTS
end

%% Clear variables or objects
% Clean-up FluidProp (if used)
if OPT.GASM == 3 % Nonideal Gas
    Cleanup_FluidProp
end
