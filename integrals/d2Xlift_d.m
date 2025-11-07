function [d2X] = d2Xlift_d(g, X, d1X, Xjmpl, Xjmpr, eq) % discontinuous
%--------------------------------------------------------------------------
% CELL
K                 = g.numE;
N                 = g.N;
d2kc              = g.d2kc;
X_dof_perm        = permute(X, [2, 3, 1]);
d2Xc              = reshape(pagemtimes(d2kc, X_dof_perm), [N, K])';
%--------------------------------------------------------------------------
% FACE
theta             =-1;
penParam          = 0;
fl                = g.fl_disp;
fr                = g.fr_disp;
dlK               = g.dlKdisp;
drK               = g.drKdisp;
flK               = g.flKdisp;
frK               = g.frKdisp;
Xl                = X*fl';
Xr                = X*fr';
d1Xl              = d1X*fl';
d1Xr              = d1X*fr';
d1Xav             = 1./2.*(d1Xr(1:K-1, 1)+d1Xl(2:K, 1));
d2X_ibp           = zeros(K, N);
d2X_sym           = zeros(K, N);
d2X_pen           = zeros(K, N);
d2X_ibp(2:K  , :) = d2X_ibp(2:K  , :)-   flK(2:K  , :).*d1Xav;                     % ibp (l)
d2X_ibp(1:K-1, :) = d2X_ibp(1:K-1, :)+   frK(1:K-1, :).*d1Xav;                     % ibp (r)
d2X_sym(2:K  , :) = d2X_sym(2:K  , :)-   dlK(2:K  , :).*Xjmpl;                     % sym (l)
d2X_sym(1:K-1, :) = d2X_sym(1:K-1, :)+   drK(1:K-1, :).*Xjmpr;                     % sym (r)
d2X_pen(2:K  , :) = d2X_pen(2:K  , :)-2.*flK(2:K  , :).*Xjmpl./g.detJ0T(2:K  , 1); % pen (l)
d2X_pen(1:K-1, :) = d2X_pen(1:K-1, :)+2.*frK(1:K-1, :).*Xjmpr./g.detJ0T(1:K-1, 1); % pen (r)
switch eq
    case 1
        d2X_ibp(1, :) = d2X_ibp(1, :)-   flK(1, :).*d1Xl(1, 1);               % ibp (l)
        d2X_ibp(K, :) = d2X_ibp(K, :)+   frK(K, :).*d1Xr(K, 1);               % ibp (r)
    case 2
        d2X_sym(1, :) = d2X_sym(1, :)-   dlK(1, :).*Xl(1, 1);                 % sym (l)
        d2X_sym(K, :) = d2X_sym(K, :)+   drK(K, :).*Xr(K, 1);                 % sym (r)
        d2X_pen(1, :) = d2X_pen(1, :)-2.*flK(1, :).*Xl(1, 1)./g.detJ0T(1, 1); % pen (l)
        d2X_pen(K, :) = d2X_pen(K, :)+2.*frK(K, :).*Xr(K, 1)./g.detJ0T(K, 1); % pen (r)
    otherwise
        return
end
%--------------------------------------------------------------------------
d2X               = d2Xc-d2X_ibp+theta.*d2X_sym+penParam.*d2X_pen;  % B. Rivière (pag. 5)
d2X               =-d2X;
%--------------------------------------------------------------------------
end