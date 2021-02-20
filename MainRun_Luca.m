% This script sets parameters and more

close all;
clear;
clc;

addpath C:\dynare\4.6.3\matlab

% setup parameters
ParameterSetup_Luca;

% Dynare step
dynare DynareStepBeforeC.mod noclearall
fprintf('\n')
%VerifySteadyState;



#PlotIRFs;