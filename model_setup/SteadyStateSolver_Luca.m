
%% steady state values 
a_ss = 1;
pi_ss = 1;
pi_h_ss = 1;
theta_ss = 1;
delta_ss = 1;
pi_star_ss = 1;
q_ss = 1/bet;
r_k_ss = (1/(1-tau_k))*((1/bet) - 1 + delt);
y_hdk_ss = (fi/(fi-1))*(r_k_ss/alfa);
ndk_ss = y_hdk_ss^(1/(1-alfa));
y_hdn_ss = y_hdk_ss/ndk_ss;
c_hdy_h_ss = (nuu/tau_c)*(s_b*(r_cap - 1) + s_f*(q_ss - 1) + s_g - tau_k*(1 - (1-alfa)*(fi - 1)/fi) - tau_n*(1-alfa)*(fi - 1)/fi - s_l);
c_f_stardy_h_ss = 1 - c_hdy_h_ss - delt*y_hdk_ss^(-1) - s_g;
c_fdy_h_ss = c_f_stardy_h_ss/0.9;
c_hdc_f_ss = c_hdy_h_ss/c_fdy_h_ss;
tt_ss = c_hdc_f_ss*(1-nuu)/nuu;
mc_ss = (fi-1)/fi*tt_ss^(nuu-1);
w_ss = mc_ss*(1-alfa)*y_hdn_ss;
cdy_h_ss = (c_hdy_h_ss^(nuu)*(c_fdy_h_ss)^(1 - nuu))/(nuu^(nuu)*(1 - nuu)^(1 - nuu));
cdn_ss = cdy_h_ss*y_hdn_ss;
n_ss = ((1/chi_n)*(1-tau_n)*w_ss*cdn_ss^(-sigma)*(1/(1+tau_c)))^(1/(eta+sigma));
y_h_ss = y_hdn_ss*n_ss;
k_ss = y_hdk_ss^(-1)*y_h_ss;
x_ss = delt*k_ss;
c_h_ss = c_hdy_h_ss*y_h_ss;
c_f_ss = c_fdy_h_ss*y_h_ss;
c_ss = (c_h_ss^(nuu) * c_f_ss^(1-nuu))/(nuu^(nuu)*(1-nuu)^(1-nuu));
f_tt_ss = s_f*y_h_ss;
z1_ss = 1/(1-bet*theta)*y_h_ss*tt_ss^(nuu-1);
z2_ss = 1/(1-bet*theta)*y_h_ss*mc_ss;
d_ss = tt_ss^(nuu-1)*y_h_ss-tt_ss^(nuu-1)*r_k_ss*k_ss-w_ss*n_ss;
m_ss = (c_ss^(-sigma)/((1+tau_c)*chi_m)*(1-bet))^(-1/muu);
b_ss = s_b*y_h_ss*(tt_ss^(nuu-1));
c_f_star = c_f_ss*0.9;
f_h_ss = (-c_f_star*tt_ss^(nuu-1) + c_f_ss*tt_ss^(nuu) + q_ss*f_tt_ss*(tt_ss^(nuu-1)) - f_tt_ss*tt_ss^(nuu-1) ) / ( q_ss*(tt_ss^(nuu_star+nuu-1)) - tt_ss^(nuu_star+nuu-1) );

%f_h = (-(0.86^(nuu-1)) * (0.9*0.22) + (0.86^(nuu)) * 0.22 + 1.0413 * (0.86^(nuu-1)) * s_f/0.59 - (s_f/0.59) * 0.86^(nuu-1)) / (1.0413 * (0.86^(nuu_star+nuu-1)) - (0.86^(nuu_star+nuu-1)) );

psi = (q_ss - q_star) / ( exp(s_f + s_b) - 1);

l_ss = ( r_cap * b_ss + q_ss * (tt_ss^(nuu-1)) * f_tt_ss ) / ( (tt_ss^(nuu-1)) * y_h_ss );

log_v_ss = ( 1 / (1-bet) ) * ( log(c_ss) - chi_n * log(n_ss) + chi_m * log(m_ss) + chi_g * s_g * log(y_h_ss) );
v_ss = exp(log_v_ss);

Test;
%% save parameters 
save parameterfile_Luca.mat;