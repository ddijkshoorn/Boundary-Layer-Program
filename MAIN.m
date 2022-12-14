%% Boundary Layer Program
% For simulation of two-dimensional steady state boundary layers in
% nonideal gas flows

% Main program based on (for references, see thesis):
% Cebeci (2002) (method of solution)
% Cebeci (2002) EV-model
% Cebeci (1974) PrT-model
% McNally (1970) (structure, functionalities and partly nomenclature)

% In short, FORTRAN program by Cebeci (2002) was extended by D.D. Dijkshoorn with:
%  o Turbulent Prandtl-number model obtained from Cebeci (1974)
%  o Fluid property models (IG) and implementation of FluidProp (NG) /// Ideal Gas (IG); Nonideal Gas (NG)
% Improved by D.D. Dijkshoorn with:
%  o coefficient of PrT (omitted by Cebeci (2002))

%% START %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

pause(0.1)

%% INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replace with any other INPUT-file from 'INPUT'-folder or 'Validation_Study/INPUT'-folder

run('INPUT/INPUT_NACA0012_Cebeci2002') %%% original, reference and working case (14-11-2021)

%% Calculate boundary conditions and solve system of PDE's %%%%%%%%%%%%%%%%

% PRECAL (pre-calculation of Boundary Conditions)
[X,NST,HVR,EDG,FRS,FLD,TCC,SOL,BLC,MON] = PRECAL(INP,OPT,SET,FLD,TCC);

% Grid generation (generate grid)
[NP,GRD] = GRID(GRD,SET);

NS = 1;

% IVPL (calculate initial values)
[sol,solprev] = IVPL(NS,GRD,HVR,FLD,SET);

tic

% Solve PDE's
[BLC,FLP,SOL,TCC,MON] = CSM(X,NP,NS,NST,GRD,HVR,EDG,FRS,FLD,TCC,sol,solprev,SOL,BLC,INP,OPT,SET,MON);

toc

%% Plot Results and generate graphs and tables %%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Simulation finished with success \n')

STORDATA
PLOTFILE
% TABLEGEN
if OPT.CHRT > 0
    FPCHARTS        % Pv-, and Ts-diagram (not available)
end

%% Clear variables or objects %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean-up FluidProp (if used)
if OPT.GASM == 3 % Nonideal Gas
    Cleanup_FluidProp
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%