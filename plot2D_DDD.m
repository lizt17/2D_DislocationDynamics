clear
set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');

outputDir = "output/";
outputVar = "outputVars/";
data = load(outputVar + "outputVars_400K.csv");

time = data(:,1);
Nd = data(:,2);
rss_source = data(:,3);
Kapp = data(:,4);
Ktip = data(:,5);
back_stress = data(:,6);
dx = data(:,7);
leadingDis = data(:,8);
leadingDisV = data(:,9);
leadingDisRss = data(:,10);

figure
hold on
% plot(Kapp, rss_source, 'LineWidth', 2, 'DisplayName', 'RSS Source');
% plot(Kapp, back_stress, 'LineWidth', 2, 'DisplayName', 'Back Stress');
plot(Kapp, leadingDisRss, 'LineWidth', 2, 'DisplayName', 'RSS dis 0');
plot(Kapp, Kapp*1e6/(161e9*sqrt(0.27e-9))./sqrt(2*pi*leadingDis+100) - 400e6/161e9 - 1./(4*pi*leadingDis), 'LineWidth', 2, 'DisplayName', 'applied stress dis 0')
xlabel('Kapp');
ylabel('Stress ($\mu$)');
legend('Location', 'best');
grid on;

figure
hold on
% plot([0; Kapp], [0; Kapp], 'LineWidth', 2, 'DisplayName', 'Kapp')
% plot([0; Kapp], [0; Ktip], 'LineWidth', 2, 'DisplayName', 'Ktip')
% plot([0; Kapp], [0; Kapp-Ktip], 'LineWidth', 2, 'DisplayName', 'Kd')
% plot([0; Kapp], [0; Nd/max(Nd)], 'o','LineWidth', 1, 'DisplayName', 'Number of dislocations')
plot(time, Kapp, 'LineWidth', 2, 'DisplayName', 'Kapp')
plot(time, Ktip, 'LineWidth', 2, 'DisplayName', 'Ktip')
plot(time, Kapp-Ktip, 'LineWidth', 2, 'DisplayName', 'Kd')
plot(time, Nd/max(Nd), 'o','LineWidth', 1, 'DisplayName', 'Number of dislocations')
xlabel('Kapp');
ylabel('SIF (MPa $\sqrt{m}$)');
legend('Location', 'best');
grid on;

% figure
% hold on
% plot(Kapp, [0; diff(Kapp-Ktip)./diff(Kapp)], 'LineWidth', 2, 'DisplayName', 'Kapp')
% xlabel('Kapp');
% ylabel('SIF (MPa $\sqrt{m}$)');
% legend('Location', 'best');
% grid on;

% figure
% hold on
% plot(Kapp, Nd, 'LineWidth', 2, 'DisplayName', 'Number of dislocations')
% ylabel('Number of dislocations');
% xlabel('Kapp');
% legend('Location', 'best');
% grid on;

% figure
% hold on
% plot(Kapp-Ktip, Nd, 'LineWidth', 2, 'DisplayName', 'Number of dislocations')
% ylabel('Number of dislocations');
% xlabel('KD');
% legend('Location', 'best');
% grid on;

figure
yyaxis left
plot(time, leadingDis, 'LineWidth', 2, 'DisplayName', 'leading dislcoation position')
ylabel('position')
yyaxis right
plot(time, leadingDisV, 'LineWidth', 2, 'DisplayName', 'leading dislcoation velocity')
legend('Location', 'best');
grid on;

numFrames = 200;

% Initialize a structure array to store all dislocation data
dislocation_data_struct = struct('frameID', {}, 'data', {});

% Loop through files dislocation_0.csv to dislocation_100.csv
for nframe = 1:numFrames
    index = nframe-1;
    filename = outputDir + "dislocation_" + index + ".csv";
    if isfile(filename)
        dislocation_data = load(filename);
        dislocation_data_struct(nframe).frameID = nframe-1; % Add frameID
        dislocation_data_struct(nframe).data = dislocation_data; % Store data
    else
        warning("File not found: " + filename);
    end
end

figure
hold on
Pmax = 100;
for n = 1:numFrames
    dis = dislocation_data_struct(n).data;

    if(~isempty(dis))
        P = dis(:,3);
        Pmax = max([Pmax; P]);
        ids = dis(:,1);
        scatter(P, time(n) * ones(size(P)), 50, ids, 'filled'); % Use scatter to color by id
        colorbar; % Add a colorbar to indicate the id values
        % plot(P, time(n), 'o', 'LineWidth', 2);
    end
end
xlabel('Dislocation Position');
ylabel('Time [b/cs]');
axis([0, Pmax, min(time), max(time)])
% legend('Location', 'best');
grid on;

% figure
% hold on
% disEnd = dislocation_data_struct(end).data;
% plot(disEnd(:,3), disEnd(:,4), '-s' ,'LineWidth', 2, 'DisplayName', 'v_dis at frame 100')
% disEnd = dislocation_data_struct(end-10).data;
% plot(disEnd(:,3), disEnd(:,4), '-s' ,'LineWidth', 2, 'DisplayName', 'v_dis at frame 90')
% disEnd = dislocation_data_struct(end-20).data;
% plot(disEnd(:,3), disEnd(:,4), '-s' ,'LineWidth', 2, 'DisplayName', 'v_dis at frame 80')
% disEnd = dislocation_data_struct(end-30).data;
% plot(disEnd(:,3), disEnd(:,4), '-s' ,'LineWidth', 2, 'DisplayName', 'v_dis at frame 70')
% legend('Location','best')