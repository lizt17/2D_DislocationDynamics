clear
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');

outoutDir = "output/";
data = load(outoutDir + "outputVars.csv");

time = data(:,1);
Nd = data(:,2);
rss_source = data(:,3);
Kapp = data(:,4);
Ktip = data(:,5);
back_stress = data(:,6);
dx = data(:,7);
leadingDis = data(:,8);
leadingDisV = data(:,9);

figure
hold on
plot(Kapp, rss_source, 'LineWidth', 2, 'DisplayName', 'RSS Source');
plot(Kapp, back_stress, 'LineWidth', 2, 'DisplayName', 'Back Stress');
xlabel('Kapp');
ylabel('Stress ($\mu$)');
legend('Location', 'best');
grid on;

figure
hold on
plot([0; Kapp], [0; Kapp], 'LineWidth', 2, 'DisplayName', 'Kapp')
plot([0; Kapp], [0; Ktip], 'LineWidth', 2, 'DisplayName', 'Ktip')
xlabel('Kapp');
ylabel('SIF (MPa $\sqrt{m}$)');
legend('Location', 'best');
grid on;

figure
hold on
plot(Nd, Kapp-Ktip, 'LineWidth', 2, 'DisplayName', 'Number of dislocations')
xlabel('Number of dislocations');
ylabel('KD');
legend('Location', 'best');
grid on;

figure
yyaxis left
plot(Kapp, leadingDis, 'LineWidth', 2, 'DisplayName', 'leading dislcoation position')
ylabel('position')
yyaxis right
plot(Kapp, leadingDisV, 'LineWidth', 2, 'DisplayName', 'leading dislcoation velocity')
legend('Location', 'best');
grid on;

numFrames = 100;

% Initialize a structure array to store all dislocation data
dislocation_data_struct = struct('frameID', {}, 'data', {});

% Loop through files dislocation_0.csv to dislocation_100.csv
% for nframe = 1:numFrames
%     index = nframe-1;
%     filename = outoutDir + "dislocation_" + index + ".csv";
%     if isfile(filename)
%         dislocation_data = load(filename);
%         dislocation_data_struct(nframe).frameID = nframe-1; % Add frameID
%         dislocation_data_struct(nframe).data = dislocation_data; % Store data
%     else
%         warning("File not found: " + filename);
%     end
% end
% 
% figure
% hold on
% Pmax = 100;
% for n = 1:numFrames
%     dis = dislocation_data_struct(n).data;
% 
%     if(~isempty(dis))
%         P = dis(:,3);
%         Pmax = max([Pmax; P]);
%         ids = dis(:,1);
%         scatter(P, time(n) * ones(size(P)), 50, ids, 'filled'); % Use scatter to color by id
%         colorbar; % Add a colorbar to indicate the id values
%         % plot(P, time(n), 'o', 'LineWidth', 2);
%     end
% end
% xlabel('Dislocation Position');
% ylabel('Time [b/cs]');
% axis([0, Pmax, min(time), max(time)])
% % legend('Location', 'best');
% grid on;