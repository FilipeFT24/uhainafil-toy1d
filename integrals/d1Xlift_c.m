function [d1X] = d1Xlift_c(g, Xdof) % continuous
%--------------------------------------------------------------------------
% CELL
K         = g.numE;
R         = g.R_disp;
bfD       = g.bfD_disp;
Xdof_perm = permute(Xdof, [2, 3, 1]);
d1Xquad   = reshape(pagemtimes(bfD, Xdof_perm), [R, K])';
d1X       = d1Xquad*g.fc;
%--------------------------------------------------------------------------
end