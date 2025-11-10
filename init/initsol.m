function [g] = initsol(g, itype, drytol, velcutoff, vellim, wetdry)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
degreec      = max(3.*g.p-1, 1);
%            = min(degreec, 12);
[Q1, W]      = quadRule1D(degreec);
N            = g.N;
K            = g.numE;
R            = size(W, 2);
bf           = g.basesOnQuad.phi1D    {degreec};
bfd          = g.basesOnQuad.gradPhi1D{degreec};
bfl          = g.basesOnQuad.phi0D;
bfr          = g.basesOnQuad.thetaPhi0D;
Wbf          = W'.*bf;
bfD          = zeros(R, N, K);
WbfD         = zeros(R, N, K);
mass         = zeros(N, N);
Mass         = zeros(N, N , K);
xydc         = zeros(K, N);
xyqc         = zeros(K, R);
for j = 1:N
    mass(:, j) = Wbf'*bf(:, j);
end
%--------------------------------------------------------------------------
for o = 1:K
    bfD (:, :, o) = 1./g.detJ0T(o, 1).*bfd;
    WbfD(:, :, o) = W'.*bfD(:, :, o);
    Mass(:, :, o) = g.detJ0T(o, 1).*mass;
end
%--------------------------------------------------------------------------
F_perquadK  = mass\Wbf';
D_perquadK  = zeros(N, R, K);
for o = 1:K
    D_perquadK(:, :, o) = mass\WbfD(:, :, o)';
end
% WEST NORMAL : LEFT FACE
Fi = zeros(K, N, 2);
for o = 2:K
    Fi(o, :, 1) = bfl/Mass(:, :, o);
end
% EAST NORMAL : RIGHT FACE
for o = 1:K-1
    Fi(o, :, 2) = bfr/Mass(:, :, o);
end
fi_aux = Wbf/mass;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for o = 1:K
    for i = 1:N
        xydc(o, i) = g.detJ0T(o, 1).*g.points(i, 1)+g.coordV0T(o, 1);
    end
    for r = 1:R
        xyqc(o, r) = g.detJ0T(o, 1).*Q1(1, r)+g.coordV0T(o, 1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g.itype     = itype;
g.drytol    = drytol;
g.velcutoff = velcutoff;
g.vellim    = vellim;
g.wetdry    = wetdry;
g.inradius  = g.detJ0T; % in_r = dx
g.R         = R;
g.W         = W;
g.bf        = bf;
g.bfl       = bfl;
g.bfr       = bfr;
g.bfD       = bfD;
g.Wbf       = Wbf;
g.DKc       = D_perquadK;
g.fkc       = F_perquadK;
g.Fi_hyp    = Fi;
g.Fi_hypl   = bfl/Mass(:, :, 1);
g.Fi_hypr   = bfr/Mass(:, :, K);
g.Mass      = Mass;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t           = 0;
X           = zeros(K, N, 2);
X(:, :, 1)  = inittype(itype, @(x) g.data.Z (x, t), xydc, xyqc, fi_aux);
X(:, :, 2)  = inittype(itype, @(x) g.data.HU(x, t), xydc, xyqc, fi_aux);
Zbdof       = inittype(itype, @(x) g.data.B (x, t), xydc, xyqc, fi_aux);
g.xydc      = xydc;
g.xyqc      = xyqc;
g.fi_aux    = fi_aux;
g.x         = X;
g.zbinit    = Zbdof;
g.zb        = Zbdof;
%--------------------------------------------------------------------------
g.fix       = false(K, 1);
if g.data.wetdry
    g = limN(g);
end
%--------------------------------------------------------------------------
% RK:
% switch g.p
%     case {0, 1}
        g.omega = [0, 1];
        g.tlvls = [1]; %#ok<NBRAK2>
%     case 2
%         g.omega = [0, 1; 1./2, 1./2];
%         g.tlvls = [0, 1];
%     case 3
%         g.omega = [0, 1; 3./4, 1./4; 1./3, 2./3];
%         g.tlvls = [0, 1, 1./2];
%     case 4
%         g.omega = [0, 1; 1/3, 2/3; 1/3, 2/3; 0, 1];
%         g.tlvls = [0, 1/2, 1/2, 1];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end