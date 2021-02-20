
test = ones(24, 1);

test(1) = bet*r_cap - 1;
test(2) = bet*q_ss - 1;
test(3) = bet*(1 - delt + (1 - tau_k)*r_k_ss) - 1;
test(4) = c_ss^(-sigma)*(1/(1 + tau_c))*(1 - bet) - chi_m*(m_ss^(-muu));
test(5) = ((1-tau_n)/(1+tau_c))*w_ss*c_ss^(-sigma) - chi_n*(n_ss^(eta));
test(6) = c_h_ss/c_f_ss - (nuu/(1-nuu))*tt_ss;
test(7) = x_ss - delt*k_ss;  
test(8) = (c_h_ss^(nuu) * c_f_ss^(1-nuu))/(nuu^(nuu)*(1-nuu)^(1-nuu)) - c_ss;
test(9) = mc_ss*(1-alfa)*k_ss^(alfa)*n_ss^(-alfa) - w_ss;
test(10) = mc_ss*alfa*k_ss^(alfa-1)*n_ss^(1-alfa) - tt_ss^(nuu-1)*r_k_ss;
test(11) = tt_ss^(nuu-1)*y_h_ss - tt_ss^(nuu-1)*r_k_ss*k_ss - w_ss*n_ss - d_ss;
test(12) = z2_ss*fi/(fi-1) - z1_ss;
test(13) = k_ss^(alfa)*n_ss^(1-alfa) - y_h_ss;
test(14) = c_h_ss + x_ss + s_g*y_h_ss + c_f_star - y_h_ss;
test(15) = r_cap*s_b*y_h_ss*tt_ss^(nuu-1) + q_ss*tt_ss^(nuu-1)*f_tt_ss + tt_ss^(nuu-1)*s_g*y_h_ss - tau_c*(tt_ss^(nuu-1)*c_h_ss + tt_ss^(nuu)*c_f_ss) - tau_k*(r_k_ss*tt_ss^(nuu-1)*k_ss + d_ss) - tau_n*n_ss*w_ss - s_l*y_h_ss*tt_ss^(nuu-1) - s_b*y_h_ss*tt_ss^(nuu-1) - tt_ss^(nuu-1)*f_tt_ss;
test(16) = -c_f_star*tt_ss^(nuu-1) + c_f_ss*tt_ss^(nuu) + q_ss*tt_ss^(nuu_star+nuu-1)*((f_tt_ss/tt_ss^(nuu_star)) - f_h_ss) - (tt_ss^(nuu_star+nuu-1))*(f_tt_ss/tt_ss^(nuu_star) - f_h_ss);
test(17) = theta + (1-theta)*(theta_ss*pi_h_ss)^(1-theta) - (theta_ss)^(1-theta);
test(18) = pi_ss / pi_h_ss - (tt_ss/tt_ss)^(1-nuu);
test(19) = tt_ss / tt_ss -  epst*pi_h_st/pi_h_ss;
test(20) = pi_star_ss / pi_h_st - (tt_ss/tt_ss)^(1-nuu_star);
test(21) = -delta_ss + theta*delta_ss*(pi_h_ss)^(fi) + (1-theta)*(theta_ss)^(-fi);
test(22) = -q_ss + q_star + psi*( exp( (f_tt_ss/y_h_ss) + s_b - d_bar ) - 1 );
test(23) = z1_ss - y_h_ss*tt_ss^(nuu-1) - bet*theta*z1_ss;
test(24) = z2_ss - y_h_ss*mc_ss - bet*theta*z2_ss;
