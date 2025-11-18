function [d2X] = d2Xlift_d(g, X, d1X, Xjmpl, Xjmpr, eq) % discontinuous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
d2c    = g.d2c;
detJ0T = g.detJ0T;
d2X    =-X*d2c'./detJ0T.^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
K      = g.numE;
N      = g.N;
zerosN = zeros(1, N);
bfl    = g.bfl_disp;
bfr    = g.bfr_disp;
dKli   = g.dKli_disp;
dKri   = g.dKri_disp;
fKli   = g.fKli_disp;
fKri   = g.fKri_disp;
d1Xl   = d1X*bfl';
d1Xr   = d1X*bfr';
d1Xav  = 1./2.*(d1Xr(1:K-1, 1)+d1Xl(2:K, 1));
d2X    = d2X-[zerosN; fKli.*d1Xav-dKli.*Xjmpl]+[fKri.*d1Xav-dKri.*Xjmpr; zerosN];
switch eq
    case 1
        fKlb      = g.fKlb_disp;
        fKrb      = g.fKrb_disp;
        d2X(1, :) = d2X(1, :)-fKlb.*d1Xl(1, 1);
        d2X(K, :) = d2X(K, :)+fKrb.*d1Xr(K, 1);
    case 2
        dKlb      = g.dKlb_disp;
        dKrb      = g.dKrb_disp;
        Xbl       = X  (1, :)*bfl';
        Xbr       = X  (K, :)*bfr';
        d2X(1, :) = d2X(1, :)+dKlb.*Xbl;
        d2X(K, :) = d2X(K, :)-dKrb.*Xbr;
    otherwise
        return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%{
function [d2X] = d2Xlift_d(g, X, d1X, Xjmpl, Xjmpr, eq) % discontinuous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
K                 = g.numE;
N                 = g.N;
d2kc              = g.d2kc;
X_dof_perm        = permute(X, [2, 3, 1]);
d2Xc              = reshape(pagemtimes(d2kc, X_dof_perm), [N, K])';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
detJ0T            = g.detJ0T;
penParam          = 0;
bfl               = g.bfl_disp;
bfr               = g.bfr_disp;
dKl               = g.dKl_disp;
dKr               = g.dKr_disp;
fKl               = g.fKl_disp;
fKr               = g.fKr_disp;
Xl                = X*bfl';
Xr                = X*bfr';
d1Xl              = d1X*bfl';
d1Xr              = d1X*bfr';
d1Xav             = 1./2.*(d1Xr(1:K-1, 1)+d1Xl(2:K, 1));
d2X_ibp           = zeros(K, N);
d2X_sym           = zeros(K, N);
d2X_pen           = zeros(K, N);
d2X_ibp(2:K  , :) = d2X_ibp(2:K  , :)-   fKl(2:K  , :).*d1Xav;                   % ibp (l)
d2X_ibp(1:K-1, :) = d2X_ibp(1:K-1, :)+   fKr(1:K-1, :).*d1Xav;                   % ibp (r)
d2X_sym(2:K  , :) = d2X_sym(2:K  , :)-   dKl(2:K  , :).*Xjmpl;                   % sym (l)
d2X_sym(1:K-1, :) = d2X_sym(1:K-1, :)+   dKr(1:K-1, :).*Xjmpr;                   % sym (r)
d2X_pen(2:K  , :) = d2X_pen(2:K  , :)-2.*fKl(2:K  , :).*Xjmpl./detJ0T(2:K  , 1); % pen (l)
d2X_pen(1:K-1, :) = d2X_pen(1:K-1, :)+2.*fKr(1:K-1, :).*Xjmpr./detJ0T(1:K-1, 1); % pen (r)
switch eq
    case 1
        d2X_ibp(1, :) = d2X_ibp(1, :)-   fKl(1, :).*d1Xl(1, 1);             % ibp (l)
        d2X_ibp(K, :) = d2X_ibp(K, :)+   fKr(K, :).*d1Xr(K, 1);             % ibp (r)
    case 2
        d2X_sym(1, :) = d2X_sym(1, :)-   dKl(1, :).*Xl(1, 1);               % sym (l)
        d2X_sym(K, :) = d2X_sym(K, :)+   dKr(K, :).*Xr(K, 1);               % sym (r)
        d2X_pen(1, :) = d2X_pen(1, :)-2.*fKl(1, :).*Xl(1, 1)./detJ0T(1, 1); % pen (l)
        d2X_pen(K, :) = d2X_pen(K, :)+2.*fKr(K, :).*Xr(K, 1)./detJ0T(K, 1); % pen (r)
    otherwise
        return
end
%--------------------------------------------------------------------------
d2X               = d2Xc-d2X_ibp+d2X_sym+penParam.*d2X_pen; % B. Rivière (pag. 5)
d2X               =-d2X;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%}