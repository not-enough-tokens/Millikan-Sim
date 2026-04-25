clear; clc; close all;

%% ── Parámetros físicos ──────────────────────────────────────────────────
p.rhoOil  = 886;                % densidad del aceite [kg/m³]
p.g       = 9.81;               % aceleración gravitatoria [m/s²]
p.d       = 6e-3;               % separación de placas [m]
p.e_exact = 1.602176634e-19;    % carga elemental exacta [C]

% Rangos de parámetros variables
p.eta_min = 1.79e-5;  p.eta_max = 1.83e-5;    % viscosidad [Pa·s]
p.rho_min = 1.18;     p.rho_max = 1.23;       % densidad del aire [kg/m³]

% Parámetros del integrador (gotas)
dt      = 1e-8;
tol_rel = 1e-6;
tol_abs = 1e-14;
t_max   = 0.5;

usar_cunningham = true;

%% ── Generación de gotas ─────────────────────────────────────────────────
N       = 30;
r       = (0.5 + 4.5*rand(1,N)) * 1e-6;    % radio [m]
n_elec  = randi(8, 1, N);                   % número de electrones
eta     = p.eta_min + (p.eta_max - p.eta_min)*rand(1,N);
rhoAir  = p.rho_min + (p.rho_max - p.rho_min)*rand(1,N);
q_real  = n_elec * p.e_exact;

% Voltaje: 20 % sobre el mínimo de equilibrio por gota
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

    % Caída libre (sin campo)
    f_caida = @(v) ( Fg(i) - 6*pi*eta(i)*r(i)*v ) / m_i;
    [~, v_hist] = euler_integrar(f_caida, 0, dt, tol_rel, tol_abs, t_max);
    vf(i)       = abs(v_hist(end));

    % Radio estimado
    r_est(i) = estimar_radio(vf(i), eta(i), delta_rho(i), p.g, usar_cunningham);

    % Ascenso con campo eléctrico
    Fg_i      = (4/3)*pi*r(i)^3 * delta_rho(i) * p.g;
    f_ascenso = @(v) ( q_real(i)*E_i - Fg_i - 6*pi*eta(i)*r(i)*v ) / m_i;
    [~, v_hist] = euler_integrar(f_ascenso, 0, dt, tol_rel, tol_abs, t_max);
    vs(i)       = abs(v_hist(end));

    % Carga estimada
    q_est(i) = 6*pi * eta(i) * r_est(i) * (vf(i) + vs(i)) * p.d / V(i);
end

%% ── Determinación de e ──────────────────────────────────────────────────
[e_sim, sigma_e, R2, n_asignado] = estimar_e(q_est);

error_pct = abs(e_sim - p.e_exact) / p.e_exact * 100;

fprintf('\n── Resultado (carga elemental) ────────────────────\n')
fprintf('  e_sim  = (%.4f ± %.4f) × 10⁻¹⁹ C\n', e_sim/1e-19, sigma_e/1e-19)
fprintf('  e_real =  1.6022 × 10⁻¹⁹ C\n')
fprintf('  Error  =  %.2f %%\n', error_pct)
fprintf('  R²     =  %.6f\n', R2)
fprintf('───────────────────────────────────────────────────\n')

%% ── Gráficas — experimento de Millikan ──────────────────────────────────
graficar_resultados(r, vf, q_est, n_asignado, e_sim, p);

%% ── Circuito RC ─────────────────────────────────────────────────────────
% Parámetros del circuito RC en el contexto de Millikan
% ---------------------------------------------------------
% Placas cuadradas de 15 cm × 15 cm → A = 0.0225 m²
%   Condición de planos infinitos: sqrt(A) >> d  (0.15 m >> 6 mm ✓)
%
% Resistencia de carga típica en fuentes de alto voltaje de laboratorio:
%   R = 5 MΩ  (limita la corriente y protege la muestra)
%
% Capacitancia calculada: C = ε₀·A/d
%   ε₀ = 8.854×10⁻¹² F/m
%   C  = 8.854e-12 × 0.0225 / 6e-3 ≈ 33.2 pF
%
% Voltaje de fuente: 500 V (rango típico 300–600 V en Millikan)
%
% Constante de tiempo: τ = R·C ≈ 5e6 × 33.2e-12 ≈ 166 µs

epsilon0 = 8.854187817e-12;    % permitividad del vacío [F/m]

p_rc.A   = 0.15^2;                        % área de las placas [m²]
p_rc.d   = p.d;                           % separación [m] — igual que gotas
p_rc.C   = epsilon0 * p_rc.A / p_rc.d;    % capacitancia [F]
p_rc.R   = 5e6;                           % resistencia [Ohm]
p_rc.V0  = 500;                           % voltaje de la fuente [V]

tau_rc   = p_rc.R * p_rc.C;

% Integrar hasta 5τ para ver la saturación completa
dt_rc    = tau_rc / 2000;          % 2000 pasos por τ — resolución adecuada
t_max_rc = 5 * tau_rc;

rc = simular_rc(p_rc, dt_rc, t_max_rc);

%% ── Gráficas — circuito RC ───────────────────────────────────────────────
graficar_rc(rc);