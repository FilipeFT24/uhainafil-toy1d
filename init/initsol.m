function [g] = initsol(g, itype, drytol, velcutoff, vellim, wetdry) % OK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE...
%--------------------------------------------------------------------------
degreec  = max(3.*g.p-1, 1);
%        = min(degreec, 12);
[Q1, W]  = quadRule1D(degreec);
K        = g.numE;
N        = g.N;
R        = size (W, 2);
xydc     = zeros(K, N);
xyqc     = zeros(K, R);
for i = 1:N
    xydc(:, i) = g.detJ0T(:, 1).*g.points(i, 1)+g.coordV0T(:, 1);
end
for r = 1:R
    xyqc(:, r) = g.detJ0T(:, 1).*Q1(1, r)+g.coordV0T(:, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
bf       = g.basesOnQuad.phi1D    {degreec};
bd       = g.basesOnQuad.gradPhi1D{degreec};
det_perm = permute(g.detJ0T, [2, 3, 1]);
Wbf      = W'.*bf;
Wbd      = W'.*bd;
bD       = pagemtimes(bd, 1./det_perm);
mass     = zeros(N, N);
for j = 1:N
    mass(:, j) = Wbf'*bf(:, j);
end
fc       = Wbf/mass;
dc       = Wbd/mass;
dc       = dc';
DKc      = pagemtimes(dc, 1./det_perm);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
bfl      = g.basesOnQuad.phi0D;
bfr      = g.basesOnQuad.thetaPhi0D;
fl       = bfl/mass;
fr       = bfr/mass;
Fi       = zeros(K, N, 2);
for o = 2:K
    Fi(o, :, 1) = fl./g.detJ0T(o, 1); % W (<-)
end
for o = 1:K-1
    Fi(o, :, 2) = fr./g.detJ0T(o, 1); % E (->)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSIGN TO...
%--------------------------------------------------------------------------
g.itype     = itype;
g.drytol    = drytol;
g.velcutoff = velcutoff;
g.vellim    = vellim;
g.wetdry    = wetdry;
g.inradius  = g.detJ0T;
g.xydc      = xydc;
g.xyqc      = xyqc;
%--------------------------------------------------------------------------
g.R         = R;
g.bf        = bf;
g.bD        = bD;
g.fc        = fc;
g.DKc       = DKc;
%--------------------------------------------------------------------------
g.bfl       = bfl;
g.bfr       = bfr;
g.Fi_       = Fi;
g.Fi1       = fl./g.detJ0T(1, 1);
g.Fik       = fr./g.detJ0T(K, 1);
%--------------------------------------------------------------------------
g.W         = W;
g.Wbf       = Wbf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INIT:
%--------------------------------------------------------------------------
t           = 0;
X           = zeros(K, N, 2);
X(:, :, 1)  = inittype(itype, @(x) g.data.Z (x, t), xydc, xyqc, fc);
X(:, :, 2)  = inittype(itype, @(x) g.data.HU(x, t), xydc, xyqc, fc);
B_dof       = inittype(itype, @(x) g.data.B (x, t), xydc, xyqc, fc);
g.x         = X;
g.zbinit    = B_dof;
g.zb        = B_dof;
%--------------------------------------------------------------------------
g.fix      = false(K, 1);
if g.data.wetdry
    g = limN(g);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TIME (RK):
%--------------------------------------------------------------------------
switch g.p
    case 0 % RK1
        g.rka = 0;
        g.rkb = 1;
        g.rkc = 0;
    case 1 % RK2
        g.rka = [0; 1];
        g.rkb = [1./2, 1./2];
        g.rkc = [0; 1];
    case 2 % RK3
        g.rka = [0, 0; 1./2, 0; -1, 2];
        g.rkb = [1./6, 2./3, 1./6];
        g.rkc = [0; 1./2; 1];
    case 3 % RK4
        g.rka = [0, 0, 0; 1./2, 0, 0; 0, 1./2, 0; 0, 0, 1];
        g.rkb = [1./6, 1./3, 1./3, 1./6];
        g.rkc = [0; 1./2; 1./2; 1];
    otherwise
        return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end