clear; clc; close all;
rng(67);

%% ── Parámetros físicos ──────────────────────────────────────────────────
p.rhoOil  = 886;                % densidad del aceite [kg/m³]
p.g       = 9.81;               % aceleración gravitatoria [m/s²]
p.d       = 6e-3;               % separación de placas [m]
p.e_exact = 1.602176634e-19;    % carga elemental exacta [C]

% Rangos de parámetros variables
p.eta_min = 1.79e-5;
p.eta_max = 1.83e-5;    % viscosidad [Pa·s]
p.rho_min = 1.18;
p.rho_max = 1.23;       % densidad del aire [kg/m³]

% Parámetros del integrador
dt      = 1e-8;     % paso temporal [s]
tol_rel = 1e-6;
tol_abs = 1e-14;
t_max   = 0.5;

usar_cunningham = true;

%% ── Generación de gotas ─────────────────────────────────────────────────
N       = 30;
r       = (0.5 + 4.5*rand(1,N)) * 1e-6;     % radio [m]
n_elec  = randi(8, 1, N);                   % número de electrones
eta     = p.eta_min + (p.eta_max - p.eta_min)*rand(1,N);
rhoAir  = p.rho_min + (p.rho_max - p.rho_min)*rand(1,N);
q_real  = n_elec * p.e_exact;

% Voltaje: 20% sobre el mínimo de equilibrio por gota
delta_rho = p.rhoOil - rhoAir;
Fg        = (4/3)*pi*r.^3 .* delta_rho * p.g;
V         = 1.2 * Fg .* p.d ./ q_real;

%% ── Simulación gota a gota ──────────────────────────────────────────────
vf    = zeros(1,N);
vs    = zeros(1,N);
r_est = zeros(1,N);
q_est = zeros(1,N);

for i = 1:N
    m_i  = (4/3)*pi*r(i)^3 * p.rhoOil;
    E_i  = V(i) / p.d;

    % Función de fuerza: caída libre (sin campo)
    f_caida = @(v) ( Fg(i) - 6*pi*eta(i)*r(i)*v ) / m_i;

    [~, v_hist] = euler_integrar(f_caida, 0, dt, tol_rel, tol_abs, t_max);
    vf(i)       = abs(v_hist(end));

    % Radio estimado (con o sin corrección de Cunningham)
    r_est(i) = estimar_radio(vf(i), eta(i), delta_rho(i), p.g, usar_cunningham);

    % Función de fuerza: ascenso con campo eléctrico
    Fg_i      = (4/3)*pi*r(i)^3 * delta_rho(i) * p.g;
    f_ascenso = @(v) ( -Fg_i + 6*pi*eta(i)*r(i)*v - q_real(i)*E_i ) / m_i;

    % La gota sube → invertimos signo de velocidad para integrar
    [~, v_hist] = euler_integrar(@(v) -f_ascenso(v), 0, dt, tol_rel, tol_abs, t_max);
    vs(i)       = abs(v_hist(end));

    % Carga estimada
    q_est(i) = 6*pi * eta(i) * r_est(i) * (vf(i) + vs(i)) * p.d / V(i);
end

%% ── Determinación de e ──────────────────────────────────────────────────
[e_sim, sigma_e, R2, n_asignado] = estimar_e(q_est);

error_pct = abs(e_sim - p.e_exact) / p.e_exact * 100;

fprintf('\n── Resultado ──────────────────────────────\n')
fprintf('  e_sim  = (%.4f ± %.4f) × 10⁻¹⁹ C\n', e_sim/1e-19, sigma_e/1e-19)
fprintf('  e_real =  1.6022 × 10⁻¹⁹ C\n')
fprintf('  Error  =  %.2f %%\n', error_pct)
fprintf('  R²     =  %.6f\n', R2)
fprintf('───────────────────────────────────────────\n')

%% ── Gráficas ────────────────────────────────────────────────────────────
graficar_resultados(r, vf, q_est, n_asignado, e_sim, p);