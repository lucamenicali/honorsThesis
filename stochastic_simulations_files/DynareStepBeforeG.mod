// ----------------
// Author: Luca Menicali
// Date: February, 2019
// ----------------

// ----------------
// Endogenous variables
// ----------------
var
    at, y_ht, ct, c_ht, c_ft, nt, xt, kt, f_ht, mt, ttt, pit, pi_ht, thetat, deltat, wt, mct, dt, r_kt, qt, 
    
    f_ttt, pi_star, z1t, z2t, rt, c_f_start, lt, sgt;

// ----------------
// Exogeneous variables
// ----------------
varexo
    etaat;

// ----------------
// Parameters
// ----------------
parameters
    a_ss, alfa, bet, nuu, muu, delt, fi, eta, sigma, nuu_star, theta, psi, chi_m, chi_n,
    
    chi_g, d_bar, rho_a, sigma_a, gamm, xi, phi_g, phi_h, r_cap, tau_c, tau_k, 
    
    tau_n, s_g, s_b, s_f, s_l, q_star, pi_h_st, epst, tt_ss, f_h_ss, f_tt_ss, cdy_h_ss, y_hdk_ss,
    
    pi_ss, phi_pi, phi_y, sigma_iota, gamm_gl, gamm_gy, l_ss, y_h_ss, c_f_star, phi_s,
    
    gamm_cl, gamm_cy, gamm_kl, gamm_ky, gamm_nl, gamm_ny, lambd_h, lambd_i, lambd_y, pi_h_ss, phi_R;

// ----------------
// Calibrated parameters
// ----------------
// parameterfile_Luca.mat is a MAT file where we had stored the parameter values
// and the steady-state values.
load parameterfile_Luca.mat;

set_param_value('alfa', alfa);
set_param_value('bet', bet);
set_param_value('nuu', nuu);
set_param_value('muu', muu);
set_param_value('delt', delt);
set_param_value('fi', fi);
set_param_value('eta', eta);
set_param_value('sigma', sigma);
set_param_value('nuu_star', nuu_star);
set_param_value('theta', theta);
set_param_value('psi', psi);
set_param_value('chi_m', chi_m);
set_param_value('chi_n', chi_n);
set_param_value('chi_g', chi_g);
set_param_value('d_bar', d_bar);
set_param_value('rho_a', rho_a);
set_param_value('sigma_a', sigma_a);
set_param_value('gamm', gamm);
set_param_value('xi', xi);
set_param_value('phi_g', phi_g);
set_param_value('phi_h', phi_h);
set_param_value('r_cap', r_cap);
set_param_value('tau_c', tau_c);
set_param_value('tau_k', tau_k);
set_param_value('tau_n', tau_n);
set_param_value('s_b', s_b);
set_param_value('s_g', s_g);
set_param_value('s_f', s_f);
set_param_value('s_l', s_l);
set_param_value('q_star', q_star);
set_param_value('epst', epst);
set_param_value('pi_h_st', pi_h_st);
set_param_value('pi_h_ss', pi_h_ss);
set_param_value('tt_ss', tt_ss);
set_param_value('f_tt_ss', f_tt_ss);
set_param_value('f_h_ss', f_h_ss);
set_param_value('a_ss', a_ss);
set_param_value('c_f_star',c_f_star);
set_param_value('y_h_ss', y_h_ss);
set_param_value('l_ss', l_ss);
set_param_value('gamm_gl', gamm_gl);
set_param_value('gamm_gy', gamm_gy);
set_param_value('gamm_kl', gamm_kl);
set_param_value('gamm_ky', gamm_ky);
set_param_value('gamm_cl', gamm_cl);
set_param_value('gamm_cy', gamm_cy);
set_param_value('gamm_nl', gamm_nl);
set_param_value('gamm_ny', gamm_ny);
set_param_value('phi_pi', phi_pi);
set_param_value('phi_y', phi_y);
set_param_value('sigma_iota', sigma_iota);
set_param_value('pi_ss', pi_ss);
set_param_value('cdy_h_ss', cdy_h_ss);
set_param_value('y_hdk_ss', y_hdk_ss);
set_param_value('lambd_i', lambd_i);
set_param_value('lambd_h', lambd_h);
set_param_value('lambd_y', lambd_y);
set_param_value('phi_R', phi_R);
set_param_value('phi_s', phi_s);


// ----------------
// Model step
// ----------------
model;
// 1. Evolution of technology
    log(exp(at)) = (1 - rho_a) * log(a_ss) + rho_a * log(exp(at(-1))) + sigma_a * etaat;

// 2. Capital law of motion (426)
    exp(kt) = exp(xt) + (1 - delt) * exp(kt(-1));

// 3. Intratemporal Euler equation for wages (424)
    chi_n * (exp(nt)^(eta)) = ( 1 - tau_n ) * exp(wt) * (exp(ct)^(-sigma)) * ( 1 / ( 1 + tau_c ) );

// 4. Intertemporal Euler equation for money holdings (423)
    chi_m * (exp(mt)^(-muu)) = (exp(ct)^(-sigma)) * (1 / ( 1 + tau_c )) - bet * (exp(ct(+1))^(-sigma)) * ( 1 / (1 + tau_c )) * (1 / exp(pit(+1)));

// 5. Allocation between home and foreign consumption goods (425)
    ( exp(c_ht) / exp(c_ft) ) = exp(ttt) * (nuu / (1 - nuu));

// 6. aggregation of consumption (427)
    exp(ct) = ( (exp(c_ht)^(nuu)) * (exp(c_ft)^(1 - nuu)) ) / ( (nuu^(nuu)) * ((1 - nuu)^(1 - nuu)) );

// 7. Marginal product of labor (428)
    exp(wt) = exp(mct) * (1 - alfa) * exp(at) * ( exp(kt(-1))^(alfa) ) * ( exp(nt)^(-alfa) );

// 8. Marginal product of capital (429)
    ( 1 / (exp(ttt)^(1 - nuu)) ) * exp(r_kt) = exp(mct) * alfa * exp(at) * ( exp(kt(-1))^(alfa - 1) ) * ( exp(nt)^(1 - alfa) );

// 9. z1 and z2 (431)
    exp(z1t) = exp(z2t) * (fi / (fi - 1));

// 10. Production function (432)
    exp(y_ht) = ( 1 / exp(deltat) ) * exp(at) * ( exp(kt(-1))^(alfa) ) * ( exp(nt)^(1 - alfa) );

// 11. Production constraint (434);
    exp(y_ht) = exp(c_ht) + exp(xt) + ( exp(sgt) * exp(y_ht) ) + exp(c_f_start);

// 12. Evolution of inflation (436);
    exp(pi_ht)^(1-fi) = theta + (1 - theta) * ( ( exp(thetat) * exp(pi_ht) )^(1-fi) );

// 13. Home inflation (437)
    exp(pit) / exp(pi_ht) = ( exp(ttt) / exp(ttt(-1)) )^(1-nuu);

// 14. Terms of trade (438)
    exp(ttt) / exp(ttt(-1)) = pi_h_st / exp(pi_ht);

// 15. Foreign inflation (439)
    exp(pi_star) / pi_h_st = ( exp(ttt(-1)) / exp(ttt) )^(1 - nuu_star);

// 16. Delta (440)
    exp(deltat) = theta * ( exp(deltat(-1)) ) * ( exp(pi_ht)^(fi) ) + (1-theta) * ( (exp(thetat))^(-fi) );

// 17. Interest rate on foreign assets (441)
    exp(qt) = q_star + psi * ( exp( exp(f_ttt) / exp(y_ht) + s_b ) - 1 );

// 18. z1 (442)
    exp(z1t) = ( exp(thetat)^(1 - fi) ) * exp(y_ht) * ( exp(ttt)^(nuu-1) ) + bet * theta * ( (exp(ct(+1))^(-sigma)) / (exp(ct)^(-sigma)) ) * ( (1 + tau_c) / (1 + tau_c )) * ( ( exp(thetat) / exp(thetat(+1)) )^(1 - fi) ) * ( ( 1 / exp(pi_ht(+1)) )^(1-fi) ) * exp(z1t(+1));

// 19. z2 (443)
    exp(z2t) = ( exp(thetat)^(-fi) ) * exp(y_ht) * exp(mct) + bet * theta * ( (exp(ct(+1))^(-sigma)) / (exp(ct)^(-sigma)) ) * ( (1 + tau_c) / (1 + tau_c) ) * ( ( exp(thetat) / exp(thetat(+1)) )^(-fi) ) * ( ( 1 / exp(pi_ht(+1)) )^(-fi) ) * exp(z2t(+1));

// 20. Interest rate (monetary policy one) (420)
    bet * ( exp(ct(+1))^(-sigma) ) * exp(rt) * ( 1 / (1 + tau_c) ) * ( 1 / exp(pit(+1)) ) = ( 1 / (1 + tau_c) ) * ( exp(ct)^(-sigma) );

// 21. Intertemporal equation on int. rate on foreign assets (421)
    bet * ( exp(ct(+1))^(-sigma) ) * exp(qt) * ( 1 / (1 + tau_c) ) * ( 1 / exp(pi_star(+1)) ) * ( exp(ttt(+1))^(nuu_star + nuu - 1) ) = ( 1 / (1 + tau_c) ) * ( exp(ct)^(-sigma) ) * ( exp(ttt)^(nuu_star + nuu - 1) ) * (1 + phi_h * ( exp(f_ht) * ( exp(ttt)^(nuu_star+nuu-1) ) - f_h_ss * ( tt_ss^(nuu_star+nuu-1) ) ) );

// 22. Intertemporal E.E. with respect to capital (422)
    bet * ( exp(ct(+1))^(-sigma) ) * (exp(ttt(+1))^(nuu - 1)) * ( 1 / (1 + tau_c) ) * (1 - delt + ( 1 - tau_k ) * exp(r_kt(+1))) = (exp(ct)^(-sigma)) * (exp(ttt)^(nuu-1)) * (1 / (1 + tau_c));

// 23. Gov't budget constraint (433)
    s_b * exp(y_ht) * (exp(ttt)^(nuu-1)) + exp(mt) + (exp(ttt)^(nuu_star+nuu-1)) * exp(f_ttt) / exp(ttt)^(nuu_star) = (phi_g/2) * (( ( exp(f_ttt) / exp(ttt)^(nuu_star) ) * (exp(ttt)^(nuu_star + nuu - 1)) - ( f_tt_ss / tt_ss^(nuu_star) ) * (tt_ss^(nuu_star + nuu - 1)) )^(2)) + exp(rt(-1)) * (1 / exp(pit)) * s_b * exp(y_ht(-1)) * (exp(ttt(-1))^(nuu-1)) + (1 / exp(pit)) * exp(mt(-1)) + exp(qt(-1)) * (exp(ttt)^(nuu_star + nuu - 1)) * (1 / exp(pi_star)) * (exp(f_ttt(-1)) / (exp(ttt(-1))^(nuu_star))) + exp(ttt)^(nuu - 1) * exp(sgt) * exp(y_ht) - tau_c * ((exp(ttt)^(nuu - 1)) * exp(c_ht) + (exp(ttt)^(nuu)) * exp(c_ft)) - tau_k * ( exp(r_kt) * ( exp(ttt)^(nuu - 1) ) * exp(kt(-1)) + exp(dt) ) - tau_n * exp(nt) * exp(wt) - s_l * exp(y_ht) * (exp(ttt)^(nuu - 1));

// 24. Trade & gov't deficits (435)
   ( exp(ttt)^(nuu_star + nuu - 1) ) * ( exp(f_ttt) / (exp(ttt)^(nuu_star)) - exp(f_ht) ) = - ( exp(ttt)^(nuu - 1) ) * exp(c_f_start) + ( exp(ttt)^(nuu) ) * exp(c_ft) + exp(qt(-1)) * (exp(ttt)^(nuu_star + nuu - 1)) * ( 1 / exp(pi_star) ) * ( exp(f_ttt(-1)) / ( exp(ttt(-1))^(nuu_star) ) - exp(f_ht(-1)) ) + (phi_h/2) * ( (exp(f_ht) * (exp(ttt)^(nuu_star + nuu - 1)) - f_h_ss * (tt_ss^(nuu_star + nuu - 1)))^(2)) + (phi_g/2) * ((exp(f_ttt) * (exp(ttt)^(nuu - 1)) - f_tt_ss * (tt_ss^(nuu - 1)))^(2));
            
// 25. Profit function (430)
    exp(dt) = ( exp(ttt)^(nuu-1) ) * exp(y_ht) - ( exp(ttt)^(nuu - 1) ) * exp(r_kt) * exp(kt(-1)) - exp(wt) * exp(nt);

// 26. Monetary union interest rate
    exp(rt) = exp(qt);
    
// 27. Public liabilities
    exp(lt) = ( exp(rt) * s_b * exp(y_ht) * (exp(ttt)^(nuu-1)) + exp(qt) * (exp(ttt(+1))^(nuu_star+nuu-1)) * (exp(pit(+1)) / exp(pi_star(+1))) * (exp(f_ttt) / exp(ttt)^(nuu_star)) ) / ( (exp(ttt)^(nuu-1)) * exp(y_ht) );    
  
// 28. Endogenous government spending
    exp(sgt) / s_g = ( exp(lt(-1)) / l_ss )^(- gamm_gl) * ( exp(y_ht(-1)) / y_h_ss )^(gamm_gy);

// 29. Endogenous consumption tax
//    exp(tau_ct) - tau_c = gamm_cl * ( exp(lt(-1)) - l_ss ) + gamm_cy * ( exp(y_ht) - y_h_ss );

// 30. Endogenous capital tax
//    exp(tau_kt) - tau_k = gamm_kl * ( exp(lt(-1)) - l_ss ) + gamm_ky * ( exp(y_ht) - y_h_ss );

// 31. Endogenous income tax
//    exp(tau_nt) - tau_n = gamm_nl * ( exp(lt(-1)) - l_ss ) + gamm_ny * ( exp(y_ht) - y_h_ss );
 
// 33. Capital to output ratio
//      exp(kdy) = exp(kt) / exp(y_ht);

// 34. Taylor Rule
//      ( exp(rt) / r_cap ) = ((exp(rt(-1)) / r_cap)^(phi_R)) * ( ( exp(pit(-1)) / pi_ss )^(phi_pi) * ( exp(y_ht(-1)) / y_h_ss )^(phi_y) )^(1-phi_R);    
//      exp(rt) = (exp(rt(-1))^(phi_R)) * (pi_ss / bet)^(1-phi_R) * ( ((exp(pit)*exp(pit(-1))*exp(pit(-2))*exp(pit(-3)))^(1/4)) / pi_ss )^(phi_pi) * ( exp(y_ht) / y_h_ss )^(phi_y) * (exp(st) / epst)^(phi_s) ;  
//       exp(rt) / r_cap = (exp(rt(-1)) / r_cap)^(phi_R) * ( ( exp(pit) / pi_ss )^(phi_pi) * ( exp(y_ht) / y_h_ss )^(phi_y) * (exp(st) / epst)^(phi_s) )^(1 - phi_R) ;  

// 35. Endogenous exchange rate
//      exp(st) = ( exp(st(-1))^(lambd_i) ) * ( (exp(pi_ht) / bet )^(1-lambd_i) ) * (( exp(pi_ht) / pi_h_ss )^(1-lambd_h) ) * ( ( exp(y_ht) / y_h_ss )^(1-lambd_y) );
//      exp(st) = (exp(st(-1))^(lambd_i)) * (exp(rt) / r_cap)^(1-lambd_i) *  ( exp(pit) / pi_ss )^(lambd_h) *  ( exp(y_ht) / y_h_ss )^(lambd_y);
//        exp(rt) - exp(pit) = (exp(qt) - exp(pit)) * (exp(st(+1)) / exp(st));
end;

// ----------------
// Steady state step
// ----------------
initval;

    at = log(a_ss);
    y_ht = log(y_h_ss);
    ct = log(c_ss);
    c_ht = log(c_h_ss);
    c_ft = log(c_f_ss);
    nt = log(n_ss);
    xt = log(x_ss);
    kt = log(k_ss);
    f_ht = log(f_h_ss);
    mt = log(m_ss);
    ttt = log(tt_ss);
    pit = log(pi_ss); 
    pi_ht = log(pi_h_ss); 
    thetat = log(theta_ss);
    deltat = log(delta_ss);
    wt = log(w_ss);
    mct = log(mc_ss);
    dt = log(d_ss);
    r_kt = log(r_k_ss);
    qt = log(q_ss);
    f_ttt = log(f_tt_ss);
    pi_star = log(pi_star_ss); 
    z1t = log(z1_ss);
    z2t = log(z2_ss);
    rt = log(r_cap);
    c_f_start = log(c_f_star);
    lt = log(l_ss);
    sgt = log(s_g);
//    tau_ct = log(tau_c);
//    tau_nt = log(tau_n);
//    tau_kt = log(tau_k);
//    cdy =  log(cdy_h_ss);
//    st = log(epst);

end;

steady;
check;

// ----------------
// Shock step
// ----------------
shocks;

    var etaat; stderr sigma_a;

end;

// ----------------
// Stochastic simulation step
// ----------------
//stoch_simul(order=1, relative_irf, irf=0, pruning, noprint);

// Optimal Simple Rule

optim_weights;

y_ht 1;
pit 1;
rt 1;
//st 1;

end;


//phi_pi = 1.1;
//phi_y = 0.2;
//phi_R = 0.9;
//phi_s = 0.06;
gamm_gl = 0.2;
gamm_gy = -0.1;
// gamm_cl = 0.08;
// gamm_cy = 0.23;
// gamm_kl = 0.1;
// gamm_ky = 0.15;
// gamm_nl = 1;
// gamm_ny = 0.15;

osr_params gamm_gl gamm_gy;


osr_params_bounds;
gamm_gl, 0.05, 2.65;
//gamm_cl, 0.08, 11.02;
// gamm_kl, 0.1, 43.9;
//gamm_nl, 0.12, 4.9; 
//phi_R, 0, 0.999;
//phi_pi, 1.001, 100;
end;

osr(order=1, irf=0, opt_algo=9) y_ht pit rt;


