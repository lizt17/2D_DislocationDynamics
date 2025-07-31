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
r_source = 1500;     % [m], source position for dislocation nucleation
T = 300; % [K], temperature
% Function to calculate resolved shear stress
tau_interaction = @(ri,rj) mu*b/(2*pi) / (ri-rj);   % dislocation interaction for screw dislocation
tau_imgae = @(r) mu*b/(4*pi) ./ (r-crack_tip);       % image force for screw dislocation

%% Define applied stress field
KappDot = 0e6 / unitSIFrate;
% Kapp0 = 0.0e6 / unitSIF;
Kapp0 = 0.0;

%% Define output variables
nucleation = false; 
dt = 10;
dxMax = 10;
Nsteps = 100;
outputInterval = 10;
time_curr = 0;
Vmax = 1.0;
time = zeros(1, Nsteps);
back_stress = zeros(1, Nsteps);

%% Initiate dislocation configuration
Nd = 1;
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

    time(kInc) = time_curr;
    Kapp = Kapp0 + KappDot * time_curr;
    tau_applied = @(r) Kapp ./ sqrt(2*pi*r);

%% Dislocation nucleation
if nucleation
    for nd = 1: Nd
        back_stress(kInc) = back_stress(kInc) + tau_interaction(r_source, disArr(nd).position); % Back stress, <0
    end
    rss_source = tau_applied(r_source) + back_stress(kInc) - tau_imgae(r_source); % Resolved shear stress at source
    if rss_source >= tau_nuc
        Nd = Nd + 1; % Increment dislocation count
        disArr(Nd).id = Nd; % Assign new ID
        disArr(Nd).position = r_source; % Set nucleation position
        disArr(Nd).velocity = 0; % Initial velocity is zero
    end
end
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
    rss = tau_app + tau_int - tau_im;


%% Mobility law
    % [currV, athermal] = mobilityLaw_W(rss, T);
    currV = rss .* b / 1;
    athermal = zeros(size(rss));
    newP = currP + currV * dt; % Update positions based on velocities

%% Visualization
    if mod(kInc-1, outputInterval) == 0
        for i = 1:Nd
            if athermal(i)
                plot(currP(i), time_curr, 's', 'LineWidth', 2, 'DisplayName', ['Dislocation ', num2str(i)]); % Square for athermal
            else
                plot(currP(i), time_curr, 'o', 'LineWidth', 2, 'DisplayName', ['Dislocation ', num2str(i)]); % Circle otherwise
            end
        end
    end

%% Update variables for next iteration
    for ndis = 1: Nd
        disArr(ndis).position = newP(ndis); % Update position
        disArr(ndis).velocity = currV(ndis); % store velocity
    end

    if Nd > 0
        x_leadingDis = disArr(Nd).position;
        Vmax = max(currV);
    else
        x_leadingDis = 1000;
        Vmax = 1.0;
    end
    dt = abs(dxMax/Vmax);
    time_curr = time_curr + dt;

end

%% 
grid on
title('Dislocation Dynamics Simulation')
axis([0, r_source*1.5, 0, time_curr]);
Kd = 1;
disp('current SIF [MPa m^0.5] = ')
disp(Kapp*unitSIF/1e6)

%% Analytical solution for single dislocation
% K-field stress
% r_dis = ( 3/2 * time * Kapp0 / sqrt(2*pi) + r_source^(3/2) ).^(2/3);
% plot(r_dis, time, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Analytical Solution');

% image stress
r_dis = sqrt(r_source^2 - mu*b*time/2/pi);
plot(r_dis, time, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Analytical Solution, image stress');