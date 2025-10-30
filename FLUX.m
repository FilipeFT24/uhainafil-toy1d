function [g] = FLUX(g, penParam, t)
%--------------------------------------------------------------------------
g     = SWEc1(g, t, penParam);
g     = SWEf1(g);
g.res = g.Fluxc-g.Fluxf;
%--------------------------------------------------------------------------
% P0 correction:
%{
if g.p > 0
    res_aux = meanval(g, g.res);
    x       = g.fix;
    for v = 1:2
        g.res(x, :, v) = repmat(res_aux(x, v), 1, g.N);
    end
end
%}
%--------------------------------------------------------------------------

c1 = g.Fluxc(:, :, 1);
f1 = g.Fluxf(:, :, 1);
c2 = g.Fluxc(:, :, 2);
f2 = g.Fluxf(:, :, 2);
r1 = g.res  (:, :, 1);
r2 = g.res  (:, :, 2);

% if max(abs())
%     xx = 1;
% end

Z_dof = g.x(:, :, 1);
Hudof = g.x(:, :, 2);
Zbdof = g.zb;
H_dof = Z_dof-Zbdof;

if any(H_dof < 0, 'all')
    xx = 1;
end

end