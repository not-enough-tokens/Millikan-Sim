function graficar_transitorio(r, eta, delta_rho, p, dt, tol_rel, tol_abs, t_max)
% GRAFICAR_TRANSITORIO  Dinámica transitoria v(t) para 3 gotas representativas.
%
%   Selecciona automáticamente las gotas en los percentiles 10, 50 y 90
%   del radio, re-ejecuta euler_integrar para cada una en caída libre
%   (E = 0) y construye una curva continua hasta 6τ de la siguiente forma:
%
%     · Segmento Euler   [0, t_stop]     : datos reales del integrador
%     · Segmento cola    [t_stop, 6τ]    : solución analítica v = vf(1−e^{−t/τ})
%
%   Ambos segmentos comparten el mismo color y grosor, produciendo una
%   curva visualmente continua. La transición es suave porque Euler ya
%   convergió en t_stop: v_Euler(t_stop) ≈ v_analitica(t_stop).
%
%   Entradas
%   --------
%   r         : vector de radios de todas las gotas [m]    (1×N)
%   eta       : vector de viscosidades              [Pa·s]  (1×N)
%   delta_rho : vector de (rhoOil − rhoAir)         [kg/m³] (1×N)
%   p         : struct de parámetros (campos: rhoOil, g)
%   dt        : paso temporal del integrador [s]
%   tol_rel   : tolerancia relativa de parada
%   tol_abs   : tolerancia absoluta de parada
%   t_max     : tiempo máximo de simulación [s]

%% ── Selección de las 3 gotas representativas ────────────────────────────

pcts      = [10, 50, 90];
r_pct     = prctile(r, pcts);
idx_gotas = arrayfun(@(rp) find(abs(r - rp) == min(abs(r - rp)), 1), r_pct);

% Límite común del eje x: 6τ de la gota más grande seleccionada
tau_gotas = arrayfun(@(i) (4/3*pi*r(i)^3*p.rhoOil) / (6*pi*eta(i)*r(i)), idx_gotas);
t_xlim    = 6 * max(tau_gotas);     % [s]

colores = {'#4d65b4', '#ae2334', '#1ebc73'};

%% ── Integración y graficado ─────────────────────────────────────────────

figure('Position', [100 100 820 480]);
hold on;

for k = 1:3
    i     = idx_gotas(k);
    r_i   = r(i);
    m_i   = (4/3) * pi * r_i^3 * p.rhoOil;
    tau_i = m_i / (6 * pi * eta(i) * r_i);
    vf_i  = 2 * r_i^2 * delta_rho(i) * p.g / (9 * eta(i));

    f_caida = @(v) (vf_i - v) / tau_i;     % dv/dt = (vf − v)/τ

    % ── Segmento Euler [0, t_stop] ───────────────────────────────────────
    [t_euler, v_euler] = euler_integrar(f_caida, 0, dt, tol_rel, tol_abs, t_max);
    t_stop = t_euler(end);

    % ── Segmento cola analítica [t_stop, 6τ] ─────────────────────────────
    %   Arranca exactamente en t_stop para que no haya hueco ni solapamiento.
    %   Se usa la solución analítica global v = vf(1 − e^{−t/τ}), que pasa
    %   por el mismo punto que Euler porque ambos comparten vf y τ.
    t_cola = linspace(t_stop, t_xlim, 300);
    v_cola = vf_i * (1 - exp(-t_cola / tau_i));

    % ── Concatenar ambos segmentos ────────────────────────────────────────
    t_total = [t_euler,  t_cola( 2:end)] * 1e6;    % µs; evita duplicar t_stop
    v_total = [v_euler,  v_cola(2:end)];

    % ── Plot en una sola llamada → curva completamente continua ──────────
    color_k = colores{k};
    plot(t_total, v_total, '-', 'Color', color_k, 'LineWidth',   1.8, 'DisplayName', sprintf('r = %.2f µm,  τ = %.1f µs', r_i*1e6, tau_i*1e6))

    % Marcador en t_stop: indica dónde terminó Euler y empieza la cola
    plot(t_stop * 1e6, v_euler(end), 'o', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', color_k, 'MarkerSize', 6, 'HandleVisibility','off')
end

%% ── Formato ─────────────────────────────────────────────────────────────

xlabel('Tiempo $t$ ($\mu$s)', 'Interpreter','latex', 'FontSize', 14);
ylabel('Velocidad $v$ (m/s)', 'Interpreter', 'latex', 'FontSize', 14)
title('Velocidad en función del tiempo durante la caída libre', 'FontSize', 13)

legend('Location', 'southeast', 'FontSize', 11)
xlim([0, t_xlim * 1e6])
grid on
box on

% Anotación: Δt y criterio de parada
annotation('textbox', [0.13 0.80 0.30 0.08], ...
    'String',           {sprintf('\\Deltat = %.0e s', dt), ...
                         '○  Fin de Euler'}, ...
    'FitBoxToText',     'on', ...
    'EdgeColor',        [0.6 0.6 0.6], ...
    'BackgroundColor',  'white', ...
    'FontSize',         10)

theme('light')
fontname('CMU Serif')
end