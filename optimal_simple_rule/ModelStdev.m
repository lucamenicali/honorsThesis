% This function computes the sum of theoretical standard deviations of
% select variables

function [sum_stdev] = ModelStdev(invars,oo_,M_)
% Output
    % sum_stdev: the sum of standard deviations of variables in listofvars
    
% Input
    % listofvars: the list of variables of interest
    % oo_: Dynare output variable
    % M_: Dynare variable with the list of endogenous variables

% Initializing sum_stdev
sum_stdev = 0;
numvars = size(invars,2);
endonames = cellstr(M_.endo_names_long);
varind = nan(numvars);
for ii = 1:numvars
    vartemp = find(ismember(endonames,invars(ii)));
    if ~isempty(vartemp)
        varind(ii) = vartemp;
    end
    sum_stdev = sum_stdev + sqrt(oo_.var(varind(ii),varind(ii)));
end