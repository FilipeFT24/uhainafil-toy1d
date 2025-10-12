function [g] = ZhangAndShuLimiter(g)
%--------------------------------------------------------------------------
is_partial = g.WD == 0;
Z_dof      = g.x(:, :, 1);
Hudof      = g.x(:, :, 2);
Zbdof      = g.zb;
H_dof      = Z_dof-Zbdof;
Hubar      = meanval(g, Hudof);
H_bar      = meanval(g, H_dof);
Hudof_1    = Hudof(is_partial, :);
Zbdof_1    = Zbdof(is_partial, :);
H_dof_1    = H_dof(is_partial, :);
Hubar_1    = Hubar(is_partial, 1);
H_bar_1    = H_bar(is_partial, 1);
%--------------------------------------------------------------------------
H_quad     = H_dof_1*g.bf';
min_H      = min(H_quad, [], 2);
theta      = min(1, (H_bar_1-g.drytol)./(H_bar_1-min_H)); % uhaina uses 0 instead of drytol
theta      = max(0, theta);
%--------------------------------------------------------------------------
% Overwrite zeta and hu:
H_                    = H_bar_1+theta.*(H_dof_1-H_bar_1);
g.x(is_partial, :, 1) = H_+Zbdof_1;
g.x(is_partial, :, 2) = Hubar_1+theta.*(Hudof_1-Hubar_1);
end