clear
%% 2D dislocation dynamics
load("matPara_W.mat");

%% Define applied stress field
KappDot = 0.02e6;
Kapp0 = 0.0;

dt = 1e-6;
Nsteps = 1000;
time_curr = dt*linspace(0, Nsteps, Nsteps+1);
Kapp = Kapp0 + KappDot * dt;

sigma_app = @(r) Kapp / sqrt(2*pi*r);

%% Initiate dislocation configuration
Nd = 0;
x_dis = zeros(1, Nd);
v_dis = zeros(1, Nd);
% Initialize dislocation group as a structure array
disArr = struct('id', {}, 'position', {}, 'velocity', {});

for i = 1: Nd
    disArr(i).id = i; % Unique identifier
    disArr(i).position = x_dis(i); % Current position
    disArr(i).velocity = v_dis(i); % Current velocity
end

%% Dislocation nucleation
% Ke = 2e6; % Pa`m^1/2
tau_nuc = 100e6; % Pa
r_source = 30*b_SI; % m
rss_source = sigma_app(r_source); % Resolved shear stress at source
if rss_source >= tau_nuc
    Nd = Nd + 1; % Increment dislocation count
    disArr(Nd).id = Nd; % Assign new ID
    disArr(Nd).position = r_source; % Set nucleation position
    disArr(Nd).velocity = 0; % Initial velocity is zero
end


%% Resolved shear stress
for i = 1: Nd
    for j = 1: Nd
        if i ~= j
            r_ij = (disArr(i).position - disArr(j).position);
            rss_ij = sigma_app(r_ij);
            % Update resolved shear stress for dislocation i
            disArr(i).resolved_shear_stress(j) = rss_ij;
        end
    end
end


%% Mobility law

%% 