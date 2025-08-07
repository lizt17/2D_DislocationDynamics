clear
load('matpara_W.mat');

data1 = load('output_200K_r300\outputVars.csv');
data2 = load('output_200K_r300_nuc1500\outputVars.csv');
data3 = load('output_200K_r500\outputVars.csv');
data4 = load('output_200K_r500_nuc1500\outputVars.csv');

Kapp1 = data1(:,4);
Ktip1 = data1(:,5);
Kapp2 = data2(:,4);
Ktip2 = data2(:,5);
Kapp3 = data3(:,4);
Ktip3 = data3(:,5);
Kapp4 = data4(:,4);
Ktip4 = data4(:,5);

figure
hold on
plot([0; 1.5*Kapp1], [0; 1.5*Kapp1], 'LineWidth', 2, 'DisplayName', 'Kapp curve');
plot([0; Kapp1], [0; Ktip1],'--', 'LineWidth', 2, 'DisplayName', 'Ktip r=300, tau_nuc=500 MPa');
plot([0; Kapp2], [0; Ktip2],'--', 'LineWidth', 2, 'DisplayName', 'Ktip r=300, tau_nuc=1500 MPa');
plot([0; Kapp3], [0; Ktip3], 'LineWidth', 2, 'DisplayName', 'Ktip r=500, tau_nuc=500 MPa');
plot([0; Kapp4], [0; Ktip4], 'LineWidth', 2, 'DisplayName', 'Ktip r=500, tau_nuc=1500 MPa');
xlabel('Kapp');
ylabel('Ktip');
grid on;
legend('Location', 'best');

data_noFriction = load('output_rate10_nofriction\outputVars.csv');
data_Friction = load('output_rate10_friction\outputVars.csv');
Kapp_noFriction = data_noFriction(:,4);
Ktip_noFriction = data_noFriction(:,5);
Kapp_Friction = data_Friction(:,4);
Ktip_Friction = data_Friction(:,5);

figure
hold on
plot([0; 1.2*Kapp_noFriction], [0; 1.2*Kapp_noFriction], 'LineWidth', 2, 'DisplayName', 'Kapp curve');
plot([0; Kapp_noFriction], [0; Ktip_noFriction], 'LineWidth', 2, 'DisplayName', 'Ktip no friction');
plot([0; Kapp_Friction], [0; Ktip_Friction], 'LineWidth', 2, 'DisplayName', 'Ktip friction 500 MPa');
xlabel('Kapp');
ylabel('Ktip');
grid on;
legend('Location', 'best');

% figure
% hold on
T = 200;
% vel = @(Kapp, r) mobilityLaw_W(Kapp/sqrt(2*pi*r),T);
% K = linspace(0.35e6, 3.1e6, 100);
% r = linspace(500, 5000, 100);
% fplot3(vel);

figure;
unitSIF = mu_SI * sqrt(b_SI);
Kapp_range = linspace(0.35e6/unitSIF, 3.1e6/unitSIF, 100); % Kapp范围
r_range = linspace(500, 5000, 100);        % r范围

% 创建网格
[Kapp_grid, r_grid] = meshgrid(Kapp_range, r_range);

% 计算每个网格点的velocity值
velocity_grid = zeros(size(Kapp_grid));
for i = 1:length(Kapp_range)
    for j = 1:length(r_range)
        velocity_grid(j, i) = mobilityLaw_W(Kapp_range(i)/sqrt(2*pi*r_range(j)), T);
    end
end

% 绘制三维曲面图
surf(Kapp_grid, r_grid, velocity_grid);
xlabel('Kapp');
ylabel('r');
zlabel('Velocity');
title('3D Surface Plot of Velocity');
colorbar;

figure;
hold on
vel = load('output_vel.txt');
for i = 1:size(vel, 1)
    if vel(i, 3) > 0
        plot(vel(:,1), vel(:,2), 'o', 'LineWidth', 0.5);
    else
        plot(vel(:,1), vel(:,2), 's', 'LineWidth', 0.5);
    end 
end
xlabel('Distance [b]');
ylabel('Velocity [cs]');