clear

outoutDir = "output/";
data = load(outoutDir+"outputVars.csv");

time = data(:,1);
Nd = data(:,2);
rss_source = data(:,3);
Kapp = data(:,4);
Ktip = data(:,5);
back_stress = data(:,6);

figure(1)
hold on
plot(time, rss_source, 'LineWidth', 2, 'DisplayName', 'RSS Source');
plot(time, back_stress, 'LineWidth', 2, 'DisplayName', 'Back Stress');
xlabel('Time');
ylabel('Stress (\mu)');
legend('Location', 'best');
grid on;

figure(2)
hold on
plot(time, Kapp, 'LineWidth', 2, 'DisplayName', 'Kapp')
plot(time, Ktip, 'LineWidth', 2, 'DisplayName', 'Ktip')
xlabel('Time');
ylabel('SIF (MPa \sqrt{m})');
legend('Location', 'best');
grid on;