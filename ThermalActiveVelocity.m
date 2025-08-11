% Plot thermal activation velocity
clear

load('matPara_W.mat');

Kapp_SI = 1.6e6; % [Pa sqrt(m)]

Kapp = Kapp_SI / unitSIF; % Non-dimensional SIF

data800K = load('output_vel_800K.txt'); % Load the output file
data400K = load('output_vel_400K.txt'); % Load the output file
data77K = load('output_vel_77K.txt'); % Load the output file
R800K = data800K(:, 1); % Distance from the crack tip
V800K = data800K(:, 2); % Velocity
athermal800K = data800K(:, 3); % Thermal activation term
dG_kinkpair800K = data800K(:, 4); % Kink pair energy

R400K = data400K(:, 1); % Distance from the crack tip
V400K = data400K(:, 2); % Velocity
athermal400K = data400K(:, 3); % Thermal activation term
dG_kinkpair400K = data400K(:, 4); % Kink pair energy

R77K = data77K(:, 1); % Distance from the crack tip
V77K = data77K(:, 2); % Velocity
athermal77K = data77K(:, 3); % Thermal activation term
dG_kinkpair77K = data77K(:, 4); % Kink pair energy

figure;
hold on;
plot(R800K(athermal800K==0), V800K(athermal800K==0), 'LineWidth', 2, 'DisplayName', 'Athermal Velocity 800K');
plot(R800K(athermal800K==1), V800K(athermal800K==1), 'LineWidth', 2, 'DisplayName', 'Thermal Active Velocity 800K');

plot(R400K(athermal400K==0), V400K(athermal400K==0), 'LineWidth', 2, 'DisplayName', 'Athermal Velocity 400K');
plot(R400K(athermal400K==1), V400K(athermal400K==1), 'LineWidth', 2, 'DisplayName', 'Thermal Active Velocity 400K');

plot(R77K(athermal77K==0), V77K(athermal77K==0), 'LineWidth', 2, 'DisplayName', 'Athermal Velocity 77K');
plot(R77K(athermal77K==1), V77K(athermal77K==1), 'LineWidth', 2, 'DisplayName', 'Thermal Active Velocity 77K');
% plot(R77K, dG_kinkpair77K, 'LineWidth', 2, 'DisplayName', 'Kink Pair Energy 77K');
xlabel('Distance from Crack Tip / b');
ylabel('Velocity /cs');
legend show;
axis([0 2000 0 0.1]);
title('Thermal Activation Velocity');
grid on;