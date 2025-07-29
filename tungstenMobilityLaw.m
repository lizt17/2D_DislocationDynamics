clear
kB = 8.6173303e-5;
Tm = 3695;
a_SI = 3.16e-10;
b_SI = a_SI * sqrt(3)/2;
mu_SI = 161e9;
h_SI = a_SI * sqrt(2/3);
w_SI = a_SI * 25;

%% Uint
a = a_SI;
b = b_SI;
h = h_SI;
w = w_SI;

%% dislocation mobility parameters
dH0 = 1.63;     % [eV]
p = 0.86;
q = 1.69;
T0 = 0.8*Tm;
tauP_SI = 2.03e9;
a0 = 1.5;
Bk = 8.3e-5;
L = 1e-7;   % length of the dislocation segment

%% 
dGkp = @(tau, T) dH0*( (1-(tau/tauP_SI/a0).^p).^q - T/T0);

B_screw = @(tau, T) a*( 2*a * exp(dGkp(tau, T)./(2*kB*T)) + L)*Bk/(2*h*L);

vs = @(tau, T) tau*b./B_screw(tau, T) .* exp(-dGkp(tau, T)./(2*kB*T));

figure
hold on
tau_applied = linspace(100e6, 2000e6, 100);
plot(tau_applied/1e6, vs(tau_applied, 300), 'LineWidth', 2, 'DisplayName', '300K')
plot(tau_applied/1e6, vs(tau_applied, 250), 'LineWidth', 2, 'DisplayName', '250K')
plot(tau_applied/1e6, vs(tau_applied, 200), 'LineWidth', 2, 'DisplayName', '200K')
legend('Location', 'best')
grid on
xlabel('applied stress [MPa]')
ylabel('screw dislocation velocity [m/s]')
% title('Screw Dislocation Velocity vs Applied Stress at 300K')
axis([100 2e3 0 0.001])

T_range = 200:50:1000; % Temperature range from 200K to 1000K
tau_critical = zeros(size(T_range)); % Preallocate array for critical stress

for i = 1:length(T_range)
    T = T_range(i);
    tau_applied = linspace(100e6, 2000e6, 1000); % Fine stress range
    v = vs(tau_applied, T);
    idx = find(v >= 1e-3, 1); % Find the first index where velocity exceeds 1e-3
    if ~isempty(idx)
        tau_critical(i) = tau_applied(idx); % Record the corresponding stress
    else
        tau_critical(i) = NaN; % If no such stress is found, set as NaN
    end
end

figure
plot(T_range, tau_critical / 1e6, '-o', 'LineWidth', 2)
grid on
xlabel('Temperature [K]')
ylabel('Critical Stress [MPa]')
title('Critical Stress vs Temperature for Velocity Threshold of 1e-3 m/s')

% Additional critical stress calculations for velocity thresholds 1e-2 and 1e-1
tau_critical_1e_2 = zeros(size(T_range)); % Preallocate array for critical stress (1e-2)
tau_critical_1e_1 = zeros(size(T_range)); % Preallocate array for critical stress (1e-1)

for i = 1:length(T_range)
    T = T_range(i);
    tau_applied = linspace(100e6, 2000e6, 1000); % Fine stress range
    v = vs(tau_applied, T);
    
    % For velocity threshold 1e-2
    idx_1e_2 = find(v >= 1e-2, 1); % Find the first index where velocity exceeds 1e-2
    if ~isempty(idx_1e_2)
        tau_critical_1e_2(i) = tau_applied(idx_1e_2); % Record the corresponding stress
    else
        tau_critical_1e_2(i) = NaN; % If no such stress is found, set as NaN
    end
    
    % For velocity threshold 1e-1
    idx_1e_1 = find(v >= 1e-1, 1); % Find the first index where velocity exceeds 1e-1
    if ~isempty(idx_1e_1)
        tau_critical_1e_1(i) = tau_applied(idx_1e_1); % Record the corresponding stress
    else
        tau_critical_1e_1(i) = NaN; % If no such stress is found, set as NaN
    end
end

% Plot all critical stress curves on the same figure
figure
hold on
plot(T_range, tau_critical / 1e6, '-o', 'LineWidth', 2, 'DisplayName', 'v >= 1e-3 m/s')
plot(T_range, tau_critical_1e_2 / 1e6, '-s', 'LineWidth', 2, 'DisplayName', 'v >= 1e-2 m/s')
plot(T_range, tau_critical_1e_1 / 1e6, '-^', 'LineWidth', 2, 'DisplayName', 'v >= 1e-1 m/s')
grid on
xlabel('Temperature [K]')
ylabel('Critical Stress [MPa]')
title('Critical Stress vs Temperature for Different Velocity Thresholds')
legend('Location', 'best')
hold off