function [g] = SWEf1(g, t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS #1:
%--------------------------------------------------------------------------
K      = g.numE;
Kf     = g.nume;
N      = g.N;
bfl    = g.bfl;
bfr    = g.bfr;
G      = g.data.G;
drytol = g.drytol;
veltol = g.velcutoff;
vellim = g.vellim;
wetdry = g.wetdry;
test   = g.test;
Fi_    = g.Fi_;
Fi1    = g.Fi1;
Fik    = g.Fik;
F_     = zeros(K , N, 2);
lambda = zeros(Kf, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS #2:
%--------------------------------------------------------------------------
N_dof  = g.x(:, :, 1);
Hudof  = g.x(:, :, 2);
B_dof  = g.zb;
%--------------------------------------------------------------------------
N_l    = N_dof*bfl';
N_r    = N_dof*bfr';
Hul    = Hudof*bfl';
Hur    = Hudof*bfr';
B_l    = B_dof*bfl';
B_r    = B_dof*bfr';
%--------------------------------------------------------------------------
switch test
    case 1
        xv   = g.coordV0T;
        xlb  = xv(1, 1);
        xrb  = xv(K, 2);
        N_lb = g.data.N (xlb, t);
        N_rb = g.data.N (xrb, t);
        Hulb = g.data.HU(xlb, t);
        Hurb = g.data.HU(xrb, t);
    case 5
        N_lb = N_l(1, 1);
        N_rb = 0.40;
        Hulb =-Hul(1, 1);
        Hurb = 0;
    otherwise
        N_lb = N_l(1, 1);
        Hulb =-Hul(1, 1);
        N_rb = N_r(K, 1);
        Hurb =-Hur(K, 1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAMBDA:
%--------------------------------------------------------------------------
% INT:
lambda(2:Kf-1, 1) = lambdaLF(drytol, veltol, vellim, wetdry, G, ...
    N_l(2:K, 1), N_r(1:K-1, 1), ...
    Hul(2:K, 1), Hur(1:K-1, 1), B_l(2:K, 1), B_r(1:K-1, 1));
% BND:
lambda(1 , 1) = lambdaLF(drytol, veltol, vellim, wetdry, G, N_l(1, 1), N_lb, Hul(1, 1), Hulb, B_l(1, 1), B_l(1, 1));
lambda(Kf, 1) = lambdaLF(drytol, veltol, vellim, wetdry, G, N_r(K, 1), N_rb, Hur(K, 1), Hurb, B_r(K, 1), B_r(K, 1));
%--------------------------------------------------------------------------
LAMBDA    = max(lambda, [], 1);
vollambda = [...
    g.inradius./lambda(1:Kf-1, 1), ...
    g.inradius./lambda(2:Kf  , 1)];
vollambda(isinf(vollambda) | isnan(vollambda)) = 1./eps;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLUX:
%--------------------------------------------------------------------------
% INT:
F_l = hydro_reconstruction2(drytol, veltol, vellim, wetdry, G, ...
    N_l(2:K  , 1), N_r(1:K-1, 1), ...
    Hul(2:K  , 1), Hur(1:K-1, 1), B_l(2:K  , 1), B_r(1:K-1, 1), LAMBDA, -1);
F_r = hydro_reconstruction2(drytol, veltol, vellim, wetdry, G, ...
    N_r(1:K-1, 1), N_l(2:K  , 1), ...
    Hur(1:K-1, 1), Hul(2:K  , 1), B_r(1:K-1, 1), B_l(2:K  , 1), LAMBDA, +1);
%--------------------------------------------------------------------------
% BND:
Fbl = hydro_reconstruction2(drytol, veltol, vellim, wetdry, G, N_l(1, 1), N_lb, Hul(1, 1), Hulb, B_l(1, 1), B_l(1, 1), LAMBDA, -1);
Fbr = hydro_reconstruction2(drytol, veltol, vellim, wetdry, G, N_r(K, 1), N_rb, Hur(K, 1), Hurb, B_r(K, 1), B_r(K, 1), LAMBDA, +1);
%--------------------------------------------------------------------------
for v = 1:2
    F_(2:K  , :, v) = F_(2:K  , :, v)+Fi_(2:K  , :, 1).*F_l(:, v);
    F_(1:K-1, :, v) = F_(1:K-1, :, v)+Fi_(1:K-1, :, 2).*F_r(:, v);
    F_(1    , :, v) = F_(1    , :, v)+Fi1.*Fbl(1, v);
    F_(K    , :, v) = F_(K    , :, v)+Fik.*Fbr(1, v);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P0 fix: CHECK!
c           = find(g.fix);
n           = sum (g.fix, 1);
F_(c, :, :) = 0;
for i = 1:n
    o = c(i, 1);
    for v = 1:2
        F_(o, :, v) = 1./g.detJ0T(o, 1).*(F_l(o-1, v)+F_r(o, v));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g.vollambda = vollambda;
g.Fluxf     = F_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end