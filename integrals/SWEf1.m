function [g] = SWEf1(g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K      = g.numE;
Kf     = g.nume;
N      = g.N;
bfl    = g.bfl;
bfr    = g.bfr;
G      = g.data.G;
drytol = g.drytol;
veltol = g.velcutoff;
rtol   = eps;
Fi_    = g.Fi_hyp;
Fil    = g.Fi_hypl;
Fir    = g.Fi_hypr;
F_     = zeros(K , N, 2);
lambda = zeros(Kf, 1);
%--------------------------------------------------------------------------
Z_dof  = g.x(:, :, 1);
Hudof  = g.x(:, :, 2);
Zbdof  = g.zb;
H_dof  = Z_dof-Zbdof;
%--------------------------------------------------------------------------
Hul    = Hudof*bfl';
Hur    = Hudof*bfr';
Z_l    = Z_dof*bfl';
Z_r    = Z_dof*bfr';
Zbl    = Zbdof*bfl';
Zbr    = Zbdof*bfr';
%--------------------------------------------------------------------------
if g.test ~= 16
    Z_lb = Z_l(1, 1);
    Hulb =-Hul(1, 1);
    Z_rb = Z_r(K, 1);
    Hurb =-Hur(K, 1);
else
    switch g.data.opt
        case 1
            Z_lb = Z_l(1, 1);
            Z_rb = 0.45;
        case 2
            Z_lb = 0.60;
            Z_rb = 0.40;
        otherwise
            return
    end
    Hulb =-Hul(1, 1);
    Hurb = 0;
    %{
    Z_lb = 1;
    Hulb = 10;
    Z_rb = Z_r(K, 1);
    Hurb =-Hur(K, 1);
    %}
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LAMBDA:
%--------------------------------------------------------------------------
% INT:
lambda(2:Kf-1, 1) = lambdaLF(drytol, veltol, G, ...
    Z_l(2:K, 1), Z_r(1:K-1, 1), ...
    Hul(2:K, 1), Hur(1:K-1, 1), Zbl(2:K, 1), Zbr(1:K-1, 1));
% BND:
lambda( 1, 1) = lambdaLF(drytol, veltol, G, Z_l(1, 1), Z_lb, Hul(1, 1), Hulb, Zbl(1, 1), Zbl(1, 1));
lambda(Kf, 1) = lambdaLF(drytol, veltol, G, Z_r(K, 1), Z_rb, Hur(K, 1), Hurb, Zbr(K, 1), Zbr(K, 1));
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
F_l = hydro_reconstruction2(drytol, veltol, rtol, G, ...
    Z_l(2:K  , 1), Z_r(1:K-1, 1), ...
    Hul(2:K  , 1), Hur(1:K-1, 1), Zbl(2:K  , 1), Zbr(1:K-1, 1), LAMBDA, -1);
F_r = hydro_reconstruction2(drytol, veltol, rtol, G, ...
    Z_r(1:K-1, 1), Z_l(2:K  , 1), ...
    Hur(1:K-1, 1), Hul(2:K  , 1), Zbr(1:K-1, 1), Zbl(2:K  , 1), LAMBDA, +1);
%--------------------------------------------------------------------------
% BND:
Fbl = hydro_reconstruction2(drytol, veltol, rtol, G, Z_l(1, 1), Z_lb, Hul(1, 1), Hulb, Zbl(1, 1), Zbl(1, 1), LAMBDA, -1);
Fbr = hydro_reconstruction2(drytol, veltol, rtol, G, Z_r(K, 1), Z_rb, Hur(K, 1), Hurb, Zbr(K, 1), Zbr(K, 1), LAMBDA, +1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTEGRALS:
%--------------------------------------------------------------------------
for v = 1:2
    F_(2:K  , :, v) = F_(2:K  , :, v)+Fi_(2:K  , :, 1).*F_l(:, v);
    F_(1:K-1, :, v) = F_(1:K-1, :, v)+Fi_(1:K-1, :, 2).*F_r(:, v);
    F_(1    , :, v) = F_(1    , :, v)+Fil.*Fbl(1, v);
    F_(K    , :, v) = F_(K    , :, v)+Fir.*Fbr(1, v);
end
%--------------------------------------------------------------------------
% P0 fix:
f           = find(g.fix);
n           = sum (g.fix, 1);
F_(f, :, :) = 0;
for i = 1:n
    o = f(i, 1);
    for v = 1:2
        F_(o, :, v) = 1./g.detJ0T(o, 1).*(F_l(o-1, v)+F_r(o, v));
    end
end
%--------------------------------------------------------------------------
g.vollambda = vollambda;
g.Fluxf     = F_;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end