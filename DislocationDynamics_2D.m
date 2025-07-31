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
tau_nuc = 200e6 / mu_SI;        % [Pa], critical resolved shear stress for dislocation nucleation
r_source = 30;     % [m], source position for dislocation nucleation
T = 300; % [K], temperature
% Function to calculate resolved shear stress
tau_interaction = @(ri,rj) mu*b/(2*pi) / (ri-rj);   % dislocation interaction for screw dislocation
tau_imgae = @(r) mu*b/(4*pi) ./ (r-crack_tip);       % image force for screw dislocation

%% Define applied stress field
KappDot = 10e6 / unitSIFrate;
% Kapp0 = 0.0e6 / unitSIF;
Kapp0 = 0.1;

%% Define output variables
nucleation = true; 
dt = 10;
dxMax = 100;
Nsteps = 100000;
outputInterval = 5000;
time_curr = 0;
Vmax = 1.0;
time = zeros(1, Nsteps);
back_stress = zeros(1, Nsteps);
rss_dis1 = zeros(1, Nsteps);
output_rss_source = zeros(1, Nsteps);

%% Initiate dislocation configuration
% Nd always represents the number of dislocations
% nextDisID-1 is total number of dislocations include annihilated ones
Nd = 1;
nextDisID = Nd + 1; 
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
    for nd = 1: Nd
        back_stress(kInc) = back_stress(kInc) + tau_interaction(r_source, disArr(nd).position); % Back stress, <0
        if(back_stress(kInc) < -0.01 || back_stress(kInc) == Inf) 
            back_stress(kInc) = -0.01;
        end
    end
    rss_source = tau_applied(r_source) + back_stress(kInc) - tau_imgae(r_source); % Resolved shear stress at source
    output_rss_source(kInc) = rss_source;
if nucleation
    if rss_source >= tau_nuc
        Nd = Nd + 1; % Increment dislocation count
        disArr(Nd).id = nextDisID; % Assign new ID
        disArr(Nd).position = r_source; % Set nucleation position
        disArr(Nd).velocity = 0; % Initial velocity is zero
        nextDisID = nextDisID + 1; % Increment next dislocation ID
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
    rss_dis1(kInc) = rss(1); % Store resolved shear stress for first dislocation

%% Mobility law
    [currV, athermal] = mobilityLaw_W(rss, T);

%% Move dislocations
    newP = currP + currV * dt; % Update positions based on velocities
    outBoundary = find(newP < crack_tip);
    newP(outBoundary) = crack_tip;

%% Visualization
    if mod(kInc-1, outputInterval) == 0
        for i = 1:Nd
            if athermal(i)
                plot(newP(i), time_curr, 's', 'LineWidth', 2, 'DisplayName', ['Dislocation ', num2str(i)]); % Square for athermal
                disp(['Dislocation ', num2str(i), ' is athermal at time ', num2str(time_curr)]);
            else
                plot(newP(i), time_curr, 'o', 'LineWidth', 2, 'DisplayName', ['Dislocation ', num2str(i)]); % Circle otherwise
            end
        end
    end

%% Update variables for next iteration
    for ndis = 1: Nd
        disArr(ndis).position = newP(ndis); % Update position
        disArr(ndis).velocity = currV(ndis); % store velocity
    end

    if Nd > 0
        x_leadingDis = disArr(1).position;
        Vmax = max(abs(currV));
    else
        x_leadingDis = 1000;
        Vmax = 1.0;
    end
    
%% Dislocation annihilation
    disArr(outBoundary) = []; % Remove dislocations that moved out of boundary
    Nd = Nd - length(outBoundary);

    dt = min(dxMax/Vmax, 1e6);
    time_curr = time_curr + dt;

end

%% 
grid on
title('Dislocation Dynamics Simulation')
axis([0, x_leadingDis*1.5, 0, time_curr]);
Kd = 1;
disp('current SIF [MPa m^0.5] = ')
disp(Kapp*unitSIF/1e6)

figure
hold on
plot(time, rss_dis1, 'LineWidth', 2, 'DisplayName', 'Resolved Shear Stress of Dislocation 1');
plot(time, back_stress, 'LineWidth', 2, 'DisplayName', 'Back Stress');
plot(time, output_rss_source, 'LineWidth', 2, 'DisplayName', 'Stress on source')
axis([0, time(end), -0.01, 0.01])
xlabel('Time [b/cs]');
ylabel('Stress [\mu]')
grid on
legend('Location', 'best')