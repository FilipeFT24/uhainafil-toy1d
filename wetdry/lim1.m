function [g] = lim1(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS:
%--------------------------------------------------------------------------
p              = g.p;
q              = quadGaussLobatto(p+1, 'Domain', [0, 1]);
N              = g.N;
R              = numel(q.Points);
bf             = zeros(R, N);
for i = 1:N
    bf(:, i) = g.BF{1, i}(q.Points);
end
%--------------------------------------------------------------------------
Z_dof          = g.x(:, :, 1);
Hudof          = g.x(:, :, 2);
Zbdof          = g.zb;
H_dof          = Z_dof-Zbdof;
H_m            = meanval(g, H_dof);
Hum            = meanval(g, Hudof);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIMITER:
%--------------------------------------------------------------------------
ptol           = 0;
rtol           = eps;
log            = H_m > rtol;
H_dofw         = H_dof(log, :);
Hudofw         = Hudof(log, :);
Zbdofw         = Zbdof(log, :);
H_mw           = H_m  (log, 1);
Humw           = Hum  (log, 1);
H_quadw        = H_dofw*bf';
min_H          = min(H_quadw, [], 2);
theta          = min(1, (H_mw-ptol)./(H_mw-min_H));
theta          = max(0, theta);
g.x(log, :, 1) = H_mw+theta.*(H_dofw-H_mw)+Zbdofw;
g.x(log, :, 2) = Humw+theta.*(Hudofw-Humw);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end