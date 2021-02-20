% This script outputs a correlation matrix of select variables

function [mcorr] = ModelCorr(invars,oo_,M_)

numvars = size(invars,1);

endonames = cellstr(M_.endo_names_long);
varind = nan(numvars,1);

for ii = 1:numvars
    vartemp = find(ismember(endonames,invars(ii)));
    if ~isempty(vartemp)
        varind(ii) = vartemp;
    end
end

mcorr = oo_.var(varind,varind);

% Octave does not support covcorr, so I had to use cov2corr
% mcorr = covcorr(mcorr,1);
mcorr = cov2corr(mcorr);