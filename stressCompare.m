clear

r_dis = 100;
syz = @(x) 1/2/pi*sqrt(r_dis/x)/(x-r_dis);
syz_screw = @(x) 1/2/pi/(x-r_dis);
KD = -1/sqrt(2*pi*r_dis);
syz_shielding = @(x) KD/sqrt(2*pi*x);
syz_image = @(x) -1/4/pi/x;
syz_sup = @(x) syz_screw(x)+syz_shielding(x) + 1/4/pi/r_dis;
% syz_sup1 = @(x) syz_screw(x)+syz_image(x)+ syz_shielding(x);
dsyz = @(x) syz(x) - syz_sup(x);
figure
hold on
fplot(syz, 'LineWidth', 2', 'DisplayName', 'analytical soultion')
fplot(syz_sup, 'LineWidth', 2, 'DisplayName', 'screw dislocation + shielding')
fplot(dsyz, 'LineWidth', 2, 'DisplayName', 'error')
% fplot(syz_sup1, 'LineWidth', 2, 'DisplayName', 'screw dislocation + image + shielding')
axis([0, 300, -0.01, 0.01])
legend('Location', 'best')

figure
hold on
fplot(@(x) syz(x)-syz_screw(x), 'LineWidth', 2, 'DisplayName', 'crack tip stress')
fplot(syz_shielding, 'LineWidth', 2, 'DisplayName', 'shielding stress')
axis([0, 300, -0.01, 0.01])
legend('Location', 'best')
