NumExp = 30;
mpexpresults = nan(NumExp,NumExp,5);

gl_grid = linspace(1,3,NumExp);
parastr_gl = 'gamm_kl';
gy_grid = linspace(-6,6,NumExp);
parastr_gy = 'gamm_ky';

ystd = ModelStdev({'y_ht'},oo_,M_);
nstd = ModelStdev({'wt'},oo_,M_);
cstd = ModelStdev({'ct'},oo_,M_);
kdystd = ModelStdev({'kdy'},oo_,M_);
baselineresults = [ystd nstd cstd kdystd];

for ii = 1:NumExp
    for jj = 1:NumExp
        set_param_value(parastr_gl,gl_grid(ii));
        set_param_value(parastr_gy,gy_grid(jj));
        [dr,~,M,~,oo] = resol(0,M_,options_,oo_);
        info = stoch_simul(var_list_);
        
        % extracting model correlation
        modelcorr = ModelCorr({'y_ht';'nt';'ct';'kdy'},oo_,M_);
        % output volatility
        mpexpresults(ii,jj,1) = ModelStdev({'y_ht'},oo_,M_);
        % wage volatility
        mpexpresults(ii,jj,2) = ModelStdev({'wt'},oo_,M_);
        % consumption volatility
        mpexpresults(ii,jj,3) = ModelStdev({'ct'},oo_,M_);
        % capital-to-output ratio volatility
        mpexpresults(ii,jj,4) = ModelStdev({'kdy'},oo_,M_);
        % sum of volatilities
        mpexpresults(ii,jj,5) = ModelStdev({'y_ht';'wt';'ct';'kdy'},oo_,M_);
    end;
end;
