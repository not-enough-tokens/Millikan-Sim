function graficar_rc(rc)
% GRAFICAR_RC  Genera las tres gráficas del circuito RC del experimento
%              de Millikan.
%
%   Entrada:
%     rc — struct devuelta por simular_rc
%
%   Gráficas producidas:
%     1. Voltaje en el capacitor vs. tiempo  (con analítica y τ marcado)
%     2. Campo eléctrico entre placas vs. tiempo
%     3. Distribución de potencial entre placas para instantes representativos

tau    = rc.tau;
p      = rc.params;
t_norm = rc.t / tau;            % tiempo normalizado para claridad

colores_snap = lines(length(rc.t_snap));   % un color por instante

figure('Position', [80 80 1250 420], 'Name', 'Circuito RC — Millikan')

% ── Panel 1: Vc(t) ───────────────────────────────────────────────────────
subplot(1, 3, 1)

% Vc numérica (aquí analítica exacta — Euler con dt << tau es indistinguible)
plot(t_norm, rc.Vc, '-', 'Color', '#2563EB', 'LineWidth', 2)
hold on

% Líneas de referencia
yline(p.V0,        '--k',  'V_0',         'LabelHorizontalAlignment', 'left', ...
      'LineWidth', 0.9, 'FontSize', 9)
yline(p.V0*(1-exp(-1)), ':r', '63.2% V_0', 'LabelHorizontalAlignment', 'left', ...
      'LineWidth', 0.9, 'FontSize', 9)
xline(1, '--r', '\tau', 'LabelVerticalAlignment', 'bottom', ...
      'LineWidth', 0.9, 'FontSize', 10)

% Marcar instantes de snapshot
for k = 1:length(rc.t_snap)
    xk = rc.t_snap(k) / tau;
    yk = rc.Vc(round(rc.t_snap(k) / (rc.t(2) - rc.t(1))) + 1);
    scatter(xk, yk, 50, colores_snap(k,:), 'filled', 'ZData', 1)
end

xlabel('Tiempo normalizado  t / \tau')
ylabel('V_C  (V)')
title('Voltaje en el capacitor')
xlim([0 max(t_norm)])
ylim([0 p.V0 * 1.1])
grid on; box on

% ── Panel 2: E(t) ────────────────────────────────────────────────────────
subplot(1, 3, 2)

plot(t_norm, rc.E * 1e-3, '-', 'Color', '#DC2626', 'LineWidth', 2)   % en kV/m
hold on

E_max = p.V0 / p.d;
yline(E_max * 1e-3, '--k', 'E_{max}', ...
      'LabelHorizontalAlignment', 'left', 'LineWidth', 0.9, 'FontSize', 9)
xline(1, '--r', '\tau', 'LabelVerticalAlignment', 'bottom', ...
      'LineWidth', 0.9, 'FontSize', 10)

% Marcar instantes de snapshot
for k = 1:length(rc.t_snap)
    xk = rc.t_snap(k) / tau;
    idx = round(rc.t_snap(k) / (rc.t(2) - rc.t(1))) + 1;
    yk  = rc.E(min(idx, length(rc.E))) * 1e-3;
    scatter(xk, yk, 50, colores_snap(k,:), 'filled', 'ZData', 1)
end

xlabel('Tiempo normalizado  t / \tau')
ylabel('E  (kV m^{-1})')
title('Campo eléctrico entre placas')
xlim([0 max(t_norm)])
ylim([0 E_max * 1.1e-3])
grid on; box on

% ── Panel 3: distribución de potencial V(x) a distintos tiempos ──────────
subplot(1, 3, 3)

hold on
for k = 1:length(rc.t_snap)
    lbl = sprintf('t = %.1f\\tau', rc.t_snap(k) / tau);
    plot(rc.x * 1e3, rc.V_dist(k,:), ...
         'Color', colores_snap(k,:), 'LineWidth', 2.0, 'DisplayName', lbl)
end

% Línea asintótica V0
plot(rc.x * 1e3, p.V0 * rc.x / p.d, '--k', ...
     'LineWidth', 1, 'DisplayName', 'V_0 (t\rightarrow\infty)')

xlabel('Posición entre placas  x  (mm)')
ylabel('V(x, t)  (V)')
title('Distribución de potencial entre placas')
legend('Location', 'nw', 'FontSize', 8)
xlim([0 p.d * 1e3])
ylim([0 p.V0 * 1.1])
grid on; box on

sgtitle('Simulación Circuito RC — Experimento de Millikan', 'FontSize', 13)

try
    theme('light')
    fontname('CMU Serif')
catch
    % Si el entorno no soporta theme/fontname, continúa sin error
end
end