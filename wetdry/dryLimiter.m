function [g] = dryLimiter(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
drytol              = g.drytol;
veltol              = g.velcutoff;
Z_dof               = g.x(:, :, 1);
Hudof               = g.x(:, :, 2);
Zbdof               = g.zbinit;
H_dof               = Z_dof-Zbdof;
%--------------------------------------------------------------------------
% isdry_drytol        = H_dof < drytol;
% isdry_veltol        = H_dof <= veltol;
% isdry               = isdry_drytol | isdry_veltol;
% Z_dof(isdry_drytol) = Zbdof(isdry_drytol);
% Hudof(isdry       ) = 0;
%--------------------------------------------------------------------------
H_aux               = Z_dof-Zbdof;
g.x(:, :, 1)        = Z_dof;
g.x(:, :, 2)        = Hudof;
g.zb                = g.zbinit;
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end