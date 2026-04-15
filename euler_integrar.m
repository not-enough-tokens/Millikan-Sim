function [t_hist, v_hist] = euler_integrar(f, v0, dt, tol_rel, tol_abs, t_max)
% EULER_INTEGRAR  Resuelve dv/dt = f(v) con el método de Euler explícito.
%
%   f        : función @(v) que devuelve dv/dt para una velocidad v dada
%   v0       : velocidad inicial [m/s]
%   dt       : paso temporal [s]
%   tol_rel  : tolerancia relativa para el criterio de paro
%   tol_abs  : tolerancia absoluta (evita división por cero cerca de v=0)
%   t_max    : tiempo máximo de simulación [s]
%
%   Devuelve vectores de tiempo y velocidad hasta alcanzar la terminal.

t = 0;
v = v0;

t_hist = t;
v_hist = v;

while t < t_max
    v_new = v + dt * f(v);          
    t     = t + dt;

    t_hist(end+1) = t;              
    v_hist(end+1) = v_new;          

    % Criterio mixto: robusto tanto cerca de v=0 como en régimen terminal
    if abs(v_new - v) < tol_abs + tol_rel * abs(v_new)
        break
    end

    v = v_new;
end
end