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
tau_nuc = 1000e6 / mu_SI;        % [Pa], critical resolved shear stress for dislocation nucleation
r_source = 50;     % [b], source position for dislocation nucleation
T = 1300; % [K], temperature
% Function to calculate resolved shear stress
tau_interaction = @(ri,rj) mu*b/(2*pi) / (ri-rj);   % dislocation interaction for screw dislocation
tau_imgae = @(r) mu*b/(4*pi) ./ (r-crack_tip);       % image force for screw dislocation

%% Define applied stress field
KappDot = 100e6 / unitSIFrate;
% Kapp0 = 0.0e6 / unitSIF;
Kapp0 = 0.13;
K0nuc = tau_nuc*sqrt(2*pi*r_source) + mu*b/2/sqrt(2*pi*r_source);
disp(['nucleation Kapp0 = ', num2str(K0nuc)])

%% Define output variables
nucleation = true; 
dt = 10;
dxMax = 20;
Nsteps = 5000000;
outputNum = 50;
outputInterval = Nsteps / outputNum;
time_curr = 0;
Vmax = 1.0;
time = zeros(1, Nsteps);
back_stress = zeros(1, Nsteps);
rss_dis1 = zeros(1, Nsteps);
output_rss_source = zeros(1, Nsteps);
outputKD = zeros(1, Nsteps);

%% Initiate dislocation configuration
% Nd always represents the number of dislocations
% nextDisID-1 is total number of dislocations include annihilated ones
Nd = 0;
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
    if kInc == 1
        Kapp = Kapp0 + KappDot * time_curr;
    else
        Kapp = Kapp0 + KappDot * time_curr - outputKD(kInc-1);
    end
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
if Nd > 0
    currP = zeros(1, Nd);
    tau_int = zeros(1, Nd);

    for i = 1: Nd
        if i > disArr(i).id
            disp('Error: i is not less than id. Stopping execution.');
            return;
        end
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
    Kd = mu*b./sqrt(2*pi*currP);
    KD = sum(Kd);
    outputKD(kInc) = KD;
    % tau_app = (Kapp-KD) ./ sqrt(2*pi*currP); 
    tau_im = tau_imgae(currP);
    rss = tau_app + tau_int - tau_im;
    if ~isempty(rss)
        rss_dis1(kInc) = rss(1); % Store resolved shear stress for first dislocation
    end

%% Mobility law
    [currV, athermal] = mobilityLaw_W(rss, T);
    Vmax = max(abs(currV));

%% Move dislocations, time increment dt 
%% Time increment calculation model 1
    % dcurrP = -diff(currP); % Calculate dislocation spacing, 1-2, 2-3, 3-4, ...
    % dcurrV = -diff(currV); % Calculate velocity difference, 1-2, 2-3, 3-4, ...
    % dtGuess = - dcurrP ./ dcurrV; 
    % dtGuess(dtGuess <= 0) = dxMax / Vmax; 
    % if isempty(dtGuess)
    %     dtGuess = dxMax / Vmax; % Default time step if no valid dtGuess
    % end
    % dt = 0.5*min(dtGuess); % Use the minimum time step to ensure stability
%% Time increment calculation model 2
    dcurrP = -diff(currP); % Calculate dislocation spacing, 1-2, 2-3, 3-4, ...
    if isempty(dcurrP(dcurrP > 0))
        dt = dxMax / Vmax; % Default time step if no dislocations
    else
        dt = min(dcurrP(dcurrP > 0)) / Vmax; 
    end
    dt = min(dt,1e4);
    newP = currP + currV * dt; % Update positions based on velocities
    outBoundary = find(newP <= crack_tip);
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
    end
    
%% Dislocation annihilation
    disArr(outBoundary) = []; % Remove dislocations that moved out of boundary
    Nd = Nd - length(outBoundary);
else
    dt = 1e3;
end % end of no dislocations case
    time_curr = time_curr + dt;

end

%% 
grid on
title('Dislocation Dynamics Simulation')
disp('current SIF [MPa m^0.5] = ')
disp(Kapp*unitSIF/1e6)
disp(['Annihilated dislocations = ', num2str(nextDisID - Nd - 1)]);
% axis([0, x_leadingDis*1.5, 0, time_curr]);

figure
hold on
plot(time, rss_dis1, 'LineWidth', 2, 'DisplayName', 'Resolved Shear Stress of Dislocation 1');
plot(time, back_stress, 'LineWidth', 2, 'DisplayName', 'Back Stress');
plot(time, output_rss_source, 'LineWidth', 2, 'DisplayName', 'Stress on source')
% axis([0, time(end), -0.01, 0.01])
xlabel('Time [b/cs]');
ylabel('Stress [\mu]')
grid on
legend('Location', 'best')

figure
hold on
plot(time, outputKD *unitSIF/1e6, 'LineWidth', 2, 'DisplayName', 'Shielding');
plot(time, (Kapp0 + KappDot * time - outputKD)*unitSIF/1e6, 'LineWidth', 2, 'DisplayName', 'K_{tip} SIF');
xlabel('Time [b/cs]');
ylabel('SIF [MPa m^{0.5}]');
legend('Location', 'best')
grid on