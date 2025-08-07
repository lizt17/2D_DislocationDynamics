clear

data77K = load('output_77K\outputVars.csv');
data300K = load('output_300K\outputVars.csv');
data1800K = load('output_1800K\outputVars.csv');
data800K = load('output_800K\outputVars.csv');
data500K = load('output_500K\outputVars.csv');
data600K = load('output_600K\outputVars.csv');
data600K_tauf1 = load('outputVars/outputVars_600K.csv');
data400K_tauf1 = load('outputVars/outputVars_400K.csv');

Kapp77K = data77K(:,4);
Ktip77K = data77K(:,5);
Kapp300K = data300K(:,4);
Ktip300K = data300K(:,5);
Kapp1800K = data1800K(:,4);
Ktip1800K = data1800K(:,5);
Kapp800K = data800K(:,4);
Ktip800K = data800K(:,5);
Kapp500K = data500K(:,4);
Ktip500K = data500K(:,5);
Kapp600K = data600K(:,4);
Ktip600K = data600K(:,5);
Kapp600K_tauf1 = data600K_tauf1(:,4);
Ktip600K_tauf1 = data600K_tauf1(:,5);
Kapp400K_tauf1 = data400K_tauf1(:,4);
Ktip400K_tauf1 = data400K_tauf1(:,5);
leadingDis800K = data800K(:,8);
leadingDis77K = data77K(:,8);
leadingDis300K = data300K(:,8);
leadingDis1800K = data1800K(:,8);
leadingDis500K = data500K(:,8);
leadingDis600K = data600K(:,8);
leadingDis600K_tauf1 = data600K_tauf1(:,8);
leadingDis400K_tauf1 = data400K_tauf1(:,8);

time = data77K(:,1);

figure
hold on
plot(time, Kapp77K - Ktip77K, 'LineWidth', 2, 'DisplayName', 'Shielding, 77K');
plot(time, Kapp300K - Ktip300K, 'LineWidth', 2, 'DisplayName', 'Shielding, 300K');
plot(time, Kapp500K - Ktip500K, 'LineWidth', 2, 'DisplayName', 'Shielding, 500K');
plot(time, Kapp600K - Ktip600K, 'LineWidth', 2, 'DisplayName', 'Shielding, 600K');
plot(time, Kapp800K - Ktip800K, 'LineWidth', 2, 'DisplayName', 'Shielding, 800K');
plot(time, Kapp1800K - Ktip1800K, 'LineWidth', 2, 'DisplayName', 'Shielding, 1800K');
xlabel('Time [b/cs]', 'FontWeight', 'bold');
ylabel(' $K_{D} \left[ MPa \sqrt{m} \right]$');
legend('Location', 'best');
grid on;

figure
hold on
plot(leadingDis77K, time, 's', 'LineWidth', 2, 'DisplayName', 'plastic zone length, 77K');
plot(leadingDis300K, time, 's', 'LineWidth', 2, 'DisplayName', 'plastic zone length, 300K');
plot(leadingDis500K, time, 's', 'LineWidth', 2, 'DisplayName', 'plastic zone length, 500K');
plot(leadingDis600K, time, 's', 'LineWidth', 2, 'DisplayName', 'plastic zone length, 600K');
plot(leadingDis800K, time, 's', 'LineWidth', 2, 'DisplayName', 'plastic zone length, 800K');
plot(leadingDis1800K, time, 's', 'LineWidth', 2, 'DisplayName', 'plastic zone length, 1800K');

ylabel('Time [b/cs]');
xlabel('Leading dis [b]', 'FontWeight', 'bold');
legend('Location', 'best');
grid on;

figure
hold on
plot(leadingDis600K, time, 's', 'LineWidth', 2, 'DisplayName', 'Leading dis 600K');
plot(leadingDis600K_tauf1, time, 'o', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', 'Leading dis 600K tauf1');
plot(leadingDis400K_tauf1, time, 'o', 'LineWidth', 1, 'MarkerSize', 4, 'DisplayName', 'Leading dis 400K tauf1');
ylabel('Time [b/cs]');
xlabel('Leading dis [b]', 'FontWeight', 'bold');
legend('Location', 'best');
grid on;