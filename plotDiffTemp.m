clear

data77K = load('outputVars_diffTemp/outputVars_77K.csv');
data200K = load('outputVars_diffTemp/outputVars_200K.csv');
data500K = load('outputVars_diffTemp/outputVars_500K.csv');

Kapp77K = data77K(:,4);
Ktip77K = data77K(:,5);
Kapp200K = data200K(:,4);
Ktip200K = data200K(:,5);
Kapp500K = data500K(:,4); % Added for 500K data
Ktip500K = data500K(:,5); % Added for 500K data


figure
hold on
plot([0; 1.8*Kapp77K], [1.6; 1.6*ones(size(Kapp77K))], '--', 'LineWidth', 2, 'DisplayName', 'KG');
plot([0; 1.2*Kapp77K], [0; 1.2*Kapp77K], 'LineWidth', 2, 'DisplayName', 'Kapp curve');
plot([0; Kapp77K], [0; Ktip77K], 'LineWidth', 2, 'DisplayName', 'Shielding, 77K');
plot([0; Kapp200K], [0; Ktip200K], 'LineWidth', 2, 'DisplayName', 'Shielding, 200K');
plot([0; Kapp500K], [0; Ktip500K], 'LineWidth', 2, 'DisplayName', 'Shielding, 500K');
xlabel('Kapp [MPa sqrt(m)]', 'FontWeight', 'bold');
ylabel('Ktip [MPa sqrt(m)]');
legend('Location', 'best');
grid on;
