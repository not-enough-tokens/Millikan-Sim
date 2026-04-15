function [e_sim, sigma_e, R2, n_asignado] = estimar_e(q_est)
% ESTIMAR_E  Determina la carga elemental a partir de un conjunto de cargas.
%
%   Etapa 1 — estimación inicial robusta por diferencias sucesivas con mediana
%   Etapa 2 — refinamiento por regresión lineal forzada al origen (2 pasadas)
%
%   Devuelve: e_sim      (carga elemental estimada [C])
%             sigma_e    (incertidumbre estándar de la pendiente [C])
%             R2         (coeficiente de determinación de la regresión)
%             n_asignado (entero de electrones por gota; 0 = descartada)

q_est = q_est(:);           % garantizar vector columna
N     = length(q_est);

%% ── Etapa 1: estimación inicial por diferencias sucesivas ────────────────
%
%  Al ordenar las cargas, la mayoría de gaps consecutivos corresponde a
%  Δn = 1 → diff ≈ e.  Los gaps de 2e, 3e… aparecen solo cuando algún
%  valor de n no fue muestreado.  La mediana es resistente a esos outliers.
%
%  Umbral adaptativo para separar ruido (gotas con mismo n, diff ≈ 0) de
%  diferencias reales (diff ≥ e):
%
%    umbral = 0.5 · min(q_est)
%
%  Fundamento: la carga más pequeña observable es 1·e, por lo que
%  min(q_est) ≥ e (salvo ruido moderado).  Las diferencias intra-grupo
%  son ruido << e, así que 0.5·min(q_est) ≈ 0.5·e las descarta sin
%  eliminar los gaps reales de ≥ 1·e.
%
%  NOTA — Por qué NO usar 0.20·median(diffs):
%    Con ~4 gotas por valor de n, la mayoría de diffs son intra-grupo (≈0).
%    median(diffs) colapsa a casi cero, el umbral colapsa con él, y los
%    near-zeros entran en diffs_pos → e_init ≈ 0 → n_asignado ~ 160 → error 99%.

q_ord = sort(q_est);
diffs = diff(q_ord);                            % N-1 diferencias consecutivas

umbral_bajo = 0.50 * min(q_est);               % ≈ 0.5·e, robusto ante ruido
diffs_pos   = diffs(diffs > umbral_bajo);       % gaps de ≥ 1 electrón

% Primera estimación como mediana de gaps válidos
e_init = median(diffs_pos);

% Refinamiento: descarta gaps de 2e, 3e… (sesgantes hacia arriba)
% para quedarse sólo con diferencias de exactamente 1 electrón
diffs_1e = diffs_pos(diffs_pos < 1.5 * e_init);
if numel(diffs_1e) >= 3
    e_init = median(diffs_1e);
end

%% ── Etapa 2a: primera regresión ─────────────────────────────────────────

n_asig_1 = round(q_est / e_init);

mask_1   = n_asig_1 > 0;                       % descarta gotas con n=0
n1       = n_asig_1(mask_1);
q1       = q_est(mask_1);

e_sim_1  = (n1' * q1) / (n1' * n1);            % mínimos cuadrados sin intercepto

%% ── Etapa 2b: segunda regresión con e refinado ───────────────────────────
%
%  Reasignar con e_sim_1 (más preciso que e_init) mueve las gotas que
%  quedaron a caballo entre dos enteros en la primera pasada.

n_asig_2 = round(q_est / e_sim_1);

mask_2   = n_asig_2 > 0;
n2       = n_asig_2(mask_2);
q2       = q_est(mask_2);

e_sim    = (n2' * q2) / (n2' * n2);

%% ── Estadísticos de la regresión final ──────────────────────────────────

N2       = numel(n2);
residuos = q2 - e_sim * n2;

sigma_e  = sqrt(sum(residuos.^2) / (N2 - 1)) / sqrt(n2' * n2);

SS_res   = sum(residuos.^2);
SS_tot   = sum((q2 - mean(q2)).^2);
R2       = 1 - SS_res / SS_tot;

%% ── Reconstruir vector de salida al tamaño original ─────────────────────

n_asignado         = zeros(N, 1);
n_asignado(mask_2) = n2;

end