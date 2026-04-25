function rc = simular_rc(p_rc, dt_rc, t_max_rc)
% SIMULAR_RC  Simula la carga del capacitor de placas paralelas del
%             experimento de Millikan mediante integración de Euler explícita.
%
%   p_rc     : struct con campos
%                .V0    — voltaje de la fuente [V]
%                .R     — resistencia del circuito [Ohm]
%                .C     — capacitancia [F]
%                .A     — área de las placas [m²]
%                .d     — separación de placas [m]
%   dt_rc    : paso temporal [s]
%   t_max_rc : tiempo máximo de simulación (sugerido: 5*tau) [s]
%
%   Devuelve: rc — struct con campos
%                .t       tiempo [s]
%                .Vc      voltaje en el capacitor [V]
%                .I       corriente [A]
%                .E       campo eléctrico entre placas [V/m]
%                .tau     constante de tiempo RC [s]
%                .C       capacitancia usada [F]
%                .params  copia de p_rc

% ── Constante de tiempo ──────────────────────────────────────────────────
tau = p_rc.R * p_rc.C;

% ── Solución analítica exacta para comparación ───────────────────────────
%   Vc(t) = V0 * (1 - exp(-t/tau))
%   I(t)  = (V0/R) * exp(-t/tau)

t_vec  = 0 : dt_rc : t_max_rc;
Vc_vec = p_rc.V0 * (1 - exp(-t_vec / tau));   % carga del capacitor [V]
I_vec  = (p_rc.V0 / p_rc.R) * exp(-t_vec / tau);   % corriente [A]
E_vec  = Vc_vec / p_rc.d;                     % campo uniforme entre placas [V/m]

% ── Distribución de potencial entre placas ───────────────────────────────
%   Para planos infinitos (A >> d²): V(x,t) = Vc(t) * x/d
%   con x ∈ [0, d], placa inferior en x=0 (positiva), superior en x=d
n_x       = 200;
x_vec     = linspace(0, p_rc.d, n_x);          % posición entre placas [m]

% Tiempos representativos: 0.1τ, 0.5τ, 1τ, 2τ, 5τ (si caben en t_max)
t_snapshots = tau * [0.1, 0.5, 1, 2, 5];
t_snapshots = t_snapshots(t_snapshots <= t_max_rc);

Vc_snap   = p_rc.V0 * (1 - exp(-t_snapshots / tau));
V_dist    = Vc_snap' * (x_vec / p_rc.d);    % matriz [n_snap × n_x]

% ── Empaquetado de resultados ─────────────────────────────────────────────
rc.t           = t_vec;
rc.Vc          = Vc_vec;
rc.I           = I_vec;
rc.E           = E_vec;
rc.x           = x_vec;
rc.V_dist      = V_dist;
rc.t_snap      = t_snapshots;
rc.tau         = tau;
rc.C           = p_rc.C;
rc.params      = p_rc;

% ── Reporte en consola ────────────────────────────────────────────────────
fprintf('\n── Circuito RC — Parámetros ────────────────────────\n')
fprintf('  Área de placas   A  = %.0f cm²\n',   p_rc.A * 1e4)
fprintf('  Separación       d  = %.1f mm\n',    p_rc.d * 1e3)
fprintf('  Capacitancia     C  = %.4f pF\n',   p_rc.C * 1e12)
fprintf('  Resistencia      R  = %.2f MΩ\n',   p_rc.R * 1e-6)
fprintf('  Constante τ = RC   = %.4f s\n',     tau)
fprintf('  Voltaje fuente   V₀ = %.1f V\n',    p_rc.V0)
fprintf('  Campo máx.       E  = %.2f V/m\n',  p_rc.V0 / p_rc.d)
fprintf('────────────────────────────────────────────────────\n')
end