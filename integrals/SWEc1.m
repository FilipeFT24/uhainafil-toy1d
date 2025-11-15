function [g] = SWEc1(g, t, penParam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS #1:
%--------------------------------------------------------------------------
K       = g.numE;
N       = g.N;
D       = 1;
bf      = g.bf;
G       = g.data.G;
dispers = g.data.dispers;
drytol  = g.drytol;
veltol  = g.velcutoff;
vellim  = g.vellim;
test    = g.test;
S       = zeros(K, N, D);
if dispers
    g = computephi(g, t, penParam);
    S = g.PHIN;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS #2:
%--------------------------------------------------------------------------
Z_dof    = g.x(:, :, 1);
Hudof    = g.x(:, :, 2);
B_dof    = g.zb;
H_dof    = Z_dof-B_dof;
d1B_dof  = d1Xlift_c(g, B_dof);
%--------------------------------------------------------------------------
Z_quad   = Z_dof*bf';
Huquad   = Hudof*bf';
B_quad   = B_dof*bf';
H_quad   = H_dof*bf';
d1B_quad = d1B_dof*bf';
GZ_quad  = G.*Z_quad;
GZ1quad  = GZ_quad.*d1B_quad;
GZ2quad  = GZ_quad.*(1./2.*Z_quad-B_quad);
%--------------------------------------------------------------------------
if vellim == 1 || (vellim == 2 && wetdry == 0)
    W2quad      = Huquad.^2./H_quad;
else
    W2quad      = H_quad.*kurganov_desingularise(H_quad.^2, Huquad.^2);
    %           =         kurganov_desingularise(H_quad   , Huquad.^2); % DO NOT ATTEMPT THIS!
    log         = H_quad < drytol | H_quad < veltol;
    Huquad(log) = 0;
    W2quad(log) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLUX:
%--------------------------------------------------------------------------
DK          = g.DKc;
fk          = g.fc;
HU          = permute(Huquad , [2, 3, 1]); % hu
W2          = permute(W2quad , [2, 3, 1]); % hu^2
GZ2         = permute(GZ2quad, [2, 3, 1]); % G/2*(z^2-2zb)
F_          =-pagemtimes(DK, [HU, W2+GZ2]);
F_          = permute(F_, [3, 1, 2]);
F_(:, :, 2) = F_(:, :, 2)+GZ1quad*fk-S;  % S is defined w/ -1/3 instead of 1/3
if test == 1 % for SW CONVERGENCE
    for i = 1:1+D
        F_(:, :, i) = F_(:, :, i)-inittype(g.itype, @(x) g.data.S{1, i}(x, t), g.xydc, g.xyqc, fk);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P0 fix: CHECK!
fix           = g.fix;
F_(fix, :, :) = 0;
F_(fix, :, 2) = F_(fix, :, 2)+S(fix, :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g.Fluxc       = F_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end