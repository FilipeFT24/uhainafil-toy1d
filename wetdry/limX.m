function [g] = limX(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fl    = g.bfl;
fr    = g.bfr;
rtol  = eps;
drytol = g.drytol;
%--------------------------------------------------------------------------
Z_dof = g.x(:, :, 1);
Hudof = g.x(:, :, 2);
Zbdof = g.zbinit;
H_dof = Z_dof-Zbdof;
H_m   = meanval(g, H_dof);
Hum   = meanval(g, Hudof);
H_l   = H_dof*fl';
H_r   = H_dof*fr';
%--------------------------------------------------------------------------
log1  = H_m < rtol | (H_l(:, 1) < rtol & H_r(:, 1) < rtol);
c1    = find(log1);
m1    = size(c1, 1);
%log2  = (H_l(:, 1) < drytol & H_r(:, 1) > drytol) | (H_l(:, 1) > drytol & H_r(:, 1) < drytol);
%c2    = find(log2);
%m2    = size(c2, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through "problematic" cells:
%--------------------------------------------------------------------------
for i = 1:m1
    o            = c1   (i, 1);
    g.x(o, :, 1) = Zbdof(o, :);
    g.x(o, :, 2) = 0;
end
%--------------------------------------------------------------------------
% for i = 1:m2
%     o            = c2 (i, 1);
%     min_H        = min(H_l(o, 1), H_r(o, 1));
%     theta        = min(1, (H_m(o, 1)-0)./(H_m(o, 1)-min_H));
%     theta        = max(0, theta);
% 
%     if theta > 0.01
%         xx = 1;
%     end
% 
%     g.x(o, :, 1) = H_m(o, 1)+theta.*(H_dof(o, :)-H_m(o, 1))+Zbdof(o, :);
%     g.x(o, :, 2) = Hum(o, 1)+theta.*(Hudof(o, :)-Hum(o, 1));
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end