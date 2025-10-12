function [M] = manning(g) % manning
%--------------------------------------------------------------------------
% CELL
G      = g.data.G;
nm     = g.data.nm;
bf     = g.bf;
Z_dof  = g.x(:, :, 1);
Hudof  = g.x(:, :, 2);
Zbdof  = g.Zbdof;
H_dof  = Z_dof-Zbdof;
Huquad = Hudof*bf';
H_quad = H_dof*bf';
Mquad  =-G.*nm.^2.*kurganov_desingularise(H_quad.^(7./3), Huquad.^2).*sign(Huquad);
M      = Mquad*g.fkc';
%--------------------------------------------------------------------------
end