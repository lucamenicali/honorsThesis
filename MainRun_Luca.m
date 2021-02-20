% This script sets parameters and more

close all;
clear;
clc;

addpath C:\dynare\4.6.3\matlab

% setup parameters
ParameterSetup_Luca;

% solve for the steady state
SteadyStateSolver_Luca;

% DYANRE STEP
% uncomment the scenario you want to run

% KEY
% C - consumption tax
% G - direct government spending
% K - capital tax
% N - income tax 

dynare DynareStepBeforeC.mod noclearall
% dynare DynareStepBeforeG.mod noclearall
% dynare DynareStepBeforeK.mod noclearall
% dynare DynareStepBeforeN.mod noclearall
% dynare DynareStepAfterC.mod noclearall
% dynare DynareStepAfterG.mod noclearall
% dynare DynareStepAfterK.mod noclearall
% dynare DynareStepAfterN.mod noclearall
fprintf('\n')



#PlotIRFs;