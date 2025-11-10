function [g] = SWEc1(g, t, penParam)
%--------------------------------------------------------------------------
K        = g.numE;
N        = g.N;
bf       = g.bf;
G        = g.data.G;
drytol   = g.drytol;
veltol   = g.velcutoff;
vellim   = g.vellim;
test     = g.test;
S        = zeros(K, N);
if ismembc(test, [1, 7])
    g = computephi(g, t, penParam);
    S = g.PHIN;
end
%--------------------------------------------------------------------------
Z_dof    = g.x(:, :, 1);
Hudof    = g.x(:, :, 2);
B_dof    = g.zb;
H_dof    = Z_dof-B_dof;
d1B_dof  = d1Xlift_c(g, B_dof);
%--------------------------------------------------------------------------
B_quad   = B_dof*bf';
d1B_quad = d1B_dof*bf';
H_quad   = H_dof*bf';
Huquad   = Hudof*bf';
Z_quad   = Z_dof*bf';
GZ_quad  = G.*Z_quad;
GZ1quad  = GZ_quad.*d1B_quad;
GZ2quad  = GZ_quad.*(1./2.*Z_quad-B_quad);
%--------------------------------------------------------------------------
if vellim == 1 || (vellim == 2 && wetdry == 0)
    W2quad      = H_quad.*(Huquad./H_quad).^2;
else
    W2quad      = H_quad.*kurganov_desingularise(H_quad.^2, Huquad.^2);
    %           =         kurganov_desingularise(H_quad   , Huquad.^2); % DO NOT ATTEMPT THIS!
    %           = Huquad.*kurganov_desingularise(H_quad   , Huquad);
    log         = H_quad < drytol | H_quad < veltol;
    Huquad(log) = 0;
    W2quad(log) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DK            = g.DKc;
fk_           = g.fkc;
HU            = permute(Huquad , [2, 3, 1]);
GZ2           = permute(GZ2quad, [2, 3, 1]);
W2            = permute(W2quad , [2, 3, 1]);
F_            = pagemtimes(DK, [HU, W2+GZ2]);
F_            = permute(F_, [3, 1, 2]);
F_(:, :, 2)   = F_(:, :, 2)-GZ1quad*fk_'+S;
if test == 2
    for i = 1:2
        F_(:, :, i) = F_(:, :, i)+inittype(g.itype, @(x) g.data.S{1, i}(x, t), g.xydc, g.xyqc, g.fi_aux);
    end
end
%--------------------------------------------------------------------------
% P0 fix:
fix           = g.fix;
F_(fix, :, :) = 0;
F_(fix, :, 2) = F_(fix, :, 2)+S(fix, :);
%--------------------------------------------------------------------------
g.Fluxc       = F_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end