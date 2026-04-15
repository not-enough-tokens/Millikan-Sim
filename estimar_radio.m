function r_est = estimar_radio(vf, eta, delta_rho, g, usar_cunningham)
% ESTIMAR_RADIO  Estima el radio de la gota a partir de su velocidad terminal.
%
%   vf               : velocidad terminal de caída [m/s]
%   eta              : viscosidad del aire [Pa·s]
%   delta_rho        : diferencia de densidades (aceite - aire) [kg/m³]
%   g                : aceleración gravitatoria [m/s²]
%   usar_cunningham  : true/false — incluir corrección de Cunningham
%
%   Sin corrección, la inversión de Stokes tiene forma cerrada.
%   Con corrección, se resuelve implícitamente con fzero.

% Solución cerrada (Stokes puro) — sirve como semilla en ambos casos
r_stokes = sqrt(9 * eta * vf / (2 * delta_rho * g));

if ~usar_cunningham
    r_est = r_stokes;
    return
end

% Corrección de Cunningham: vf_real = vf_Stokes * (1 + A*lambda/r)
% Reordenando: (2r²·Δρ·g / 9η)(1 + A·λ/r) - vf = 0
A      = 0.864;
lambda = 68e-9;     % recorrido libre medio del aire [m]

ecuacion = @(r) (2*r^2 * delta_rho * g / (9*eta)) * (1 + A*lambda/r) - vf;

% fzero busca la raíz partiendo de la solución sin corrección
r_est = fzero(ecuacion, r_stokes);
end