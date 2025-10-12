function [g] = limX(g)
%--------------------------------------------------------------------------
% dry limiter:
g          = dryLimiter(g);
%--------------------------------------------------------------------------
% wet/dry detection:
Z_dof      = g.x(:, :, 1);
Zbdof      = g.zbinit;
H_dof      = Z_dof-Zbdof;
is_dry     = all(H_dof <= g.drytol, 2);
is_wet     = all(H_dof  > g.drytol, 2);
is_partial =~is_dry & ~is_wet;
if g.p > 0
    g.WD(is_wet    , 1) = 0; % wet
    g.WD(is_partial, 1) = 1; % wet and dry
    g.WD(is_dry    , 1) = 2; % dry
end
% -------------------------------------------------------------------------
% wet/dry:
if g.p > 0
    g = wetDryCellAverage(g);
end
%--------------------------------------------------------------------------
end