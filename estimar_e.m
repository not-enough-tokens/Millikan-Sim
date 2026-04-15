function [e_sim, sigma_e, R2, n_asignado] = estimar_e(q_est)
% ESTIMAR_E  Determina la carga elemental a partir de un conjunto de cargas.
%
%   Etapa 1 — estimación inicial robusta por diferencias por pares
%   Etapa 2 — refinamiento por regresión lineal forzada al origen
%
%   Devuelve: e_sim     (carga elemental estimada [C])
%             sigma_e   (incertidumbre estándar de la pendiente [C])
%             R2        (coeficiente de determinación de la regresión)
%             n_asignado (número entero de electrones asignado a cada gota)

% --- Etapa 1: diferencias por pares ---
q_ord      = sort(q_est);
N          = length(q_ord);
candidatos = [];

for i = 1:N
    for j = i+1:N
        delta = q_ord(j) - q_ord(i);
        for k = 1:7                         
            candidatos(end+1) = delta / k;  
        end
    end
end

% El histograma de candidatos tiene un pico en e: buscamos ese pico
[conteos, bordes] = histcounts(candidatos, 'BinMethod', 'fd');
[~, idx]          = max(conteos);
e_init            = (bordes(idx) + bordes(idx+1)) / 2;

% --- Etapa 2: regresión lineal forzada al origen ---
n_asignado = round(q_est / e_init);

% Mínimos cuadrados para q̂ = e·n → e = (nᵀq̂)/(nᵀn)
n  = n_asignado(:);     % columna
q  = q_est(:);

e_sim   = (n' * q) / (n' * n);

% Incertidumbre estándar de la pendiente
residuos = q - e_sim * n;
sigma_e  = sqrt(sum(residuos.^2) / (N - 1)) / sqrt(n' * n);

% R²
SS_res = sum(residuos.^2);
SS_tot = sum((q - mean(q)).^2);
R2     = 1 - SS_res / SS_tot;
end