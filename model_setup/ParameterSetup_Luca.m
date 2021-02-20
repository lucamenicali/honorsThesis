% % This script sets up annual frequency parameters and data

% % Model parameters
alfa = 0.42;                % capital return
bet = 0.9603;               % household discount factor
nuu = 0.5;                  % home goods bias parameter 
muu = 3.42;                 % parameter related to money demand elasticity
delt = 0.1;                 % capital depreciation 
fi = 6;                     % price elasticity of demand
eta = 1;                    % Frish labor elasticity
sigma = 1;                  % elasticity of intertemporal substitution
nuu_star = 0.5;             % foreign goods bias parameter   
theta = 0.5;                % price rigidity parameter
% psi = 0.05;               % risk premia parameter
chi_m = 0.001;              % preference parameter related to real money holdings
chi_n = 7;                  % preference parameter related to hours worked
chi_g = 0.1;                % preference parameter related to public spending
d_bar = 0.9;                % threshold value for share of public debt
rho_a = 0.92;               % persistence of TFP
sigma_a = 0.1;              % st dev of TFP 
gamm = 0.9;                 % foreign imports parameter
xi = 2;                     % adjustment cost parameter on physical capital 
phi_g = 2;                  % adjustment cost parameter on foreign public debt
phi_h = 2;                  % adjustment cost parameter on foreign assets
r_cap = 1/bet;              % gross nominal interest rate
tau_c = 0.17;               % consumption tax
tau_k = 0.32;               % capital tax
tau_n = 0.42;               % income tax
s_g = 0.22;                 % government spending as a share of output 
s_b = 0.6;                 % domestic public debt as a share of output 
s_f = 0.5;                 % foreign public debt as a share of output 
s_l = -0.21;                % residually determined from s_g, s_b, and s_f              
epst = 1;                   % exogenous exchange rate 


% % rest of the world variables, given
q_star = 1.0303;
pi_h_st = 1;

% fiscal policy parameters
gamm_gl = 0.05;
gamm_cl = 0.2;
gamm_kl = 0.1;
gamm_nl = 0.12;
gamm_gy = 0.07;
gamm_cy = 0.4;
gamm_ky = 0.15;
gamm_ny = 0.15;

%Taylor Rule
phi_pi = 3;
phi_y = gamm_gy;
iota = 0.85;
sigma_iota = 0;
phi_R = 0.9;

% exchange rate 
lambd_i = 0.92;
lambd_h = 1.09;
lambd_y = 0.14;
phi_s = 0.01;

