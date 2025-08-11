clear

% material properties, BCC W
mu_SI = 161e9;   % [Pa]
b_SI = 0.27e-9; % [m]
KG_SI = 3.1e6;  % [Pa`m^0.5]
nu = 0.28;
tauP_SI = 2030e6;  % [Pa]
E_SI = 2*(1+nu)*mu_SI;  % [Pa]
rho_SI = 19250; % [kg/m^3]
cs_SI = sqrt(mu_SI/rho_SI); % [m/s], speed of shear wave

% unified
b = 1;
mu = 1;
tauP = tauP_SI / mu_SI;
KG = KG_SI / (mu_SI*sqrt(b_SI));
Elas = E_SI / mu_SI;

% mobility law parameters
dH0 = 1.63;     % [eV]
Tc = 800;       % [K]
Be_SI = 3.3e-7;    % [Pa`s/K]
Bk_SI = 8.3e-5;    % [Pa`s]
p = 0.999;
q = 1;

kB_eV = 8.617e-5;

DBTT = [600, 773];

% Derive parameters
unitSIF = mu_SI * sqrt(b_SI);
unitSIFrate = mu_SI * cs_SI / sqrt(b_SI);
unitTime = b_SI / cs_SI;

save('matPara_W.mat')
