
function [velocity, athermal] = mobilityLaw_W(tau_in, Temp)
% mobilityLaw_W calculates the velocity of dislocations in tungsten
% based on the applied stress (tau_in) and temperature (Temp).
%% material parameters - tungsten
    kB = 8.6173303e-5;  % Boltzmann constant, [eV/K]
    Tm = 3695; % Melting temperature of tungsten, [K]
    a_SI = 3.16e-10; % Lattice parameter of tungsten, [m]
    b_SI = a_SI * sqrt(3)/2;
    mu_SI = 161e9; % Shear modulus of tungsten, [Pa]
    h_SI = a_SI * sqrt(2/3);
    w_SI = a_SI * 25;
    rho_SI = 19250; % Density of tungsten, [kg/m^3]
    cs_SI = sqrt(mu_SI/rho_SI); % Speed of shear wave in tungsten, [m/s]

%% dislocation mobility parameters
    dH0 = 1.63; % Activation enthalpy, [eV]
    p = 0.86; % Fitting parameter
    q = 1.69; % Fitting parameter
    T0 = 0.8*Tm;
    tauP_SI = 2.03e9; % Peierls stress, [Pa]
    a0 = 1.5;
    Bk_SI = 8.3e-5; % Drag coefficient, [Pa/s]
    B0_SI = 9.8e-4; % Drag coefficient in athermal regime, [Pa/s]
    L_SI = 1e-7;   % length of the dislocation segment, [m]

%% Normalize units
    a = a_SI / b_SI;
    b = b_SI / b_SI;
    h = h_SI / b_SI;
    w = w_SI / b_SI;
    mu = mu_SI / mu_SI;
    tauP = tauP_SI / mu_SI;
    L = L_SI / b_SI;
    Bk = Bk_SI*cs_SI/(mu_SI*b_SI);
    B0 = B0_SI*cs_SI/(mu_SI*b_SI);

%% calculate dislocation velocity
    % dGkp = @(tau, T) dH0*( (1-(tau/tauP/a0).^p) .^q - T/T0);
    % B_screw = @(tau, T) a*( 2*a * exp(dGkp(tau, T)./(2*kB*T)) + L)*Bk/(2*h*L);
    % vs = @(tau, T) tau*b./B_screw(tau, T) .* exp(-dGkp(tau, T)./(2*kB*T));
    Theta = tau_in/tauP/a0;
    astress = Theta < 1; 
    dGkp = astress .* dH0.*( (1-Theta.^p) .^q - Temp/T0 );
    athermal = dGkp <= 0; % Check if in athermal regime
    B_screw = a*( 2*a * exp(dGkp./(2*kB*Temp)) + L)*Bk/(2*h*L);
    B_eff = athermal .* B0 + ~athermal .* B_screw; % Effective drag coefficient
    velocity = tau_in*b ./ B_eff .* ...
            (athermal + ~athermal .* exp(-dGkp./(2*kB*Temp)) ); % Dislocation velocity

end