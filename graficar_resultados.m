function graficar_resultados(r, vf, q_est, n_asignado, e_sim, p)

figure('Position', [100 100 1200 400]);

% ── Panel 1: vf vs r (verificación de Stokes) ───────────────────────────
subplot(1,3,1)
r_teo  = linspace(min(r), max(r), 200);
vf_teo = 2*r_teo.^2 * mean(p.rhoOil - 1.205) * p.g / (9 * mean([p.eta_min p.eta_max]));

plot(r*1e6, vf, 'o', 'MarkerFaceColor', '#4d65b4', 'MarkerEdgeColor', 'none')
hold on
plot(r_teo*1e6, vf_teo, '--k', 'LineWidth', 1.2)
xlabel('Radio $r$ ($\mu$m)', 'Interpreter','latex', 'FontSize', 14);
ylabel('Velocidad terminal $v_f$ (m/s)', 'Interpreter', 'latex', 'FontSize', 14)
title('Verificación de Stokes');
legend('Simulados','Teórica','Location','nw')
grid on

% ── Panel 2: cuantización de la carga ───────────────────────────────────
subplot(1,3,2)
scatter(1:length(q_est), q_est/p.e_exact, 50, n_asignado, 'filled')

% Add this line to change the color scheme:
colormap('turbo') % or colormap(turbo) for a smoother look

yline(1:8, '--', 'Color', [0.7 0.7 0.7])
colorbar;  
clim([1 8]) % This ensures 1 is blue and 8 is red
xlabel('Numero de gota ($n_i$)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Carga normalizada ($\hat{q}/e$)', 'Interpreter', 'latex', 'FontSize', 14)
title('Cuantización de la carga')
grid on

% ── Panel 3: regresión lineal ────────────────────────────────────────────
subplot(1,3,3)
n_rng = linspace(0, max(n_asignado)+0.5, 100);
plot(n_asignado, q_est/1e-19, 'o', 'MarkerFaceColor', '#1ebc73', 'MarkerEdgeColor','none')
hold on
plot(n_rng, e_sim*n_rng/1e-19, '-k', 'LineWidth', 1.5)
xlabel('Numero de carga asignado ($n_i$)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Carga estimada $\hat{q}_i$ ($10^{-19}$ C)', 'Interpreter', 'latex', 'FontSize', 14)
title(sprintf('Regresión: e = %.4f \\times 10^{-19} C', e_sim/1e-19))
legend('Datos','Regresión','Location','nw')
grid on

sgtitle('Detemrinación de la Carga Elemental mediante un Análisis del Experimento de Millikan', 'FontSize', 13)

theme('light')
fontname('CMU Serif')
end