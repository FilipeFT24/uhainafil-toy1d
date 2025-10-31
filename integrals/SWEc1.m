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
if ismembc(test, [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    g = computephi(g, t, penParam);
    S = g.PHIN;
end
%--------------------------------------------------------------------------
Z_dof    = g.x(:, :, 1);
Hudof    = g.x(:, :, 2);
Zbdof    = g.zb;
H_dof    = Z_dof-Zbdof;
%--------------------------------------------------------------------------
%{
fl       = g.fl_disp;
fr       = g.fr_disp;
Zl       = Z_dof*fl';
Zr       = Z_dof*fr';
Zbl      = Zbdof*fl';
Zbr      = Zbdof*fr';
Zbjmpl   = hydro_reconstruction0(Zl(2:K  , 1), Zbl(2:K  , 1), Zbr(1:K-1, 1));
Zbjmpr   = hydro_reconstruction0(Zr(1:K-1, 1), Zbr(1:K-1, 1), Zbl(2:K  , 1));
%}
%--------------------------------------------------------------------------
d1Zbdof  = d1Xlift_c(g, Zbdof);
%--------------------------------------------------------------------------
H_quad   = H_dof*bf';
Z_quad   = Z_dof*bf';
Zbquad   = Zbdof*bf';
d1Zbquad = d1Zbdof*bf';
GZ_quad  = G.*Z_quad;
GZ1quad  = GZ_quad.*d1Zbquad;
GZ2quad  = GZ_quad.*(1./2.*Z_quad-Zbquad);
Huquad   = Hudof*bf';
%--------------------------------------------------------------------------
switch vellim
    case 1
        W2quad = H_quad.*(Huquad./H_quad).^2;
    case 2
        W2quad = H_quad.*kurganov_desingularise(H_quad.^2, Huquad.^2);
        %      =         kurganov_desingularise(H_quad   , Huquad.^2); % DO NOT ATTEMPT THIS!
        %      = Huquad.*kurganov_desingularise(H_quad   , Huquad);
    otherwise
        return
end
log         = H_quad < drytol | H_quad < veltol;
Huquad(log) = 0;
W2quad(log) = 0;

if any(H_quad < 0, 'all')
    xx = 1;
end
if max(abs(W2quad), [], 'all') > 1
    xx = 1;
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
%--------------------------------------------------------------------------
% P0 fix:
fix           = g.fix;
F_(fix, :, :) = 0;
F_(fix, :, 2) = F_(fix, :, 2)+S(fix, :);
%--------------------------------------------------------------------------
g.Fluxc       = F_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end