clear; close; close all;

% Parameters
n = 3; % You can change this value to see the effect
x_r = linspace(0, 1.1, 500);

% Functions
Fa = 1 - exp((x_r).^n - 1) ./ (exp(1) - 1);
Fg = 1 - exp((1 - x_r).^n - 1) ./ (exp(1) - 1);

% Plot
figure;
plot(x_r, Fa, 'b-', 'LineWidth', 1.5); hold on;
plot(x_r, (1-Fa).*Fg, 'r--', 'LineWidth', 1.5);
xlabel('$x_r$', 'Interpreter', 'latex');
ylabel('$F_a(x), F_g(x)$', 'Interpreter', 'latex');
legend({'$F_a(x)$', '$F_g(x)$'}, 'Interpreter', 'latex', 'Location', 'best');
title(['$F_a(x)$ and $F_g(x)$ for $n = ', num2str(n), '$'], 'Interpreter', 'latex');
grid on;
