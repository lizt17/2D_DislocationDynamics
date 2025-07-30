clear
%% 2D dislocation dynamics
%% Define parameters
load("matPara_W.mat");
unitSIF = mu_SI*sqrt(b_SI);
unitSIFrate = mu_SI*cs_SI/sqrt(b_SI);
unitTime = b_SI / cs_SI;

% Ke = 2e6; % Pa`m^1/2
% Assming crack-tip located at origin
crack_tip = 0;
tau_nuc = 100e6 / mu_SI;        % [Pa], critical resolved shear stress for dislocation nucleation
r_source = 30;     % [m], source position for dislocation nucleation
T = 300; % [K], temperature
% Function to calculate resolved shear stress
tau_interaction = @(ri,rj) mu*b/(2*pi) / (ri-rj);   % dislocation interaction for screw dislocation
tau_imgae = @(r) mu*b/(2*pi) ./ (r-crack_tip);       % image force for screw dislocation

%% Define applied stress field
KappDot = 100e6 / unitSIFrate;
Kapp0 = 0.0;

dt = 100;
Nsteps = 10000;
outputInterval = 1000;

time_curr = dt*linspace(0, Nsteps, Nsteps+1);
Kapp = Kapp0 + KappDot * time_curr;

%% Initiate dislocation configuration
Nd = 2;
disArr = struct('id', {}, 'position', {}, 'velocity', {});

for i = 1: Nd
    disArr(i).id = i; % Unique identifier
    disArr(i).position = r_source*i; % Current position
    disArr(i).velocity = 0.0; % Current velocity
end
%% 

figure
hold on
xlabel('Dislocation Position [b]')
ylabel('Time [b/cs]')

for kInc = 1: Nsteps
tau_applied = @(r) Kapp(kInc) ./ sqrt(2*pi*r);

%% Resolved shear stress
currP = zeros(1, Nd);
tau_int = zeros(1, Nd);

for i = 1: Nd
    ri = disArr(i).position;
    currP(i) = ri;
    for j = 1: Nd
        if i ~= j
            rj = disArr(j).position;
            tau_int(i) = tau_interaction(ri, rj) + tau_int(i);
        end
    end
end

tau_app = tau_applied(currP); 
tau_im = tau_imgae(currP);
rss = tau_app + tau_int + tau_im;

%% Dislocation nucleation
% rss_source = sigma_app(r_source); % Resolved shear stress at source
% if rss_source >= tau_nuc
%     Nd = Nd + 1; % Increment dislocation count
%     disArr(Nd).id = Nd; % Assign new ID
%     disArr(Nd).position = r_source; % Set nucleation position
%     disArr(Nd).velocity = 0; % Initial velocity is zero
% end

%% Mobility law
currV = mobilityLaw_W(rss, T);
newP = currP + currV * dt; % Update positions based on velocities

%% Visualization
if mod(kInc, outputInterval) == 0
    plot(currP, time_curr(kInc), 'o' ,'LineWidth', 2, 'DisplayName', 'Dislocation Positions');
end

end

%% 
grid on
title('Dislocation Dynamics Simulation')
axis([0, 1000, 0, time_curr(end)]);
Kd = 1;