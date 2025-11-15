function [d1X] = d1Xlift_d(g, X, Xjmpl, Xjmpr, eq) % discontinuous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
K      = g.numE;
fx     = g.fx_disp;
detJ0T = g.detJ0T;
d1X    = X*fx./detJ0T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
N      = g.N;
zerosN = zeros(1, N);
fKli   = g.fKli_disp;
fKri   = g.fKri_disp;
d1X    = d1X+[zerosN; fKli.*Xjmpl]-[fKri.*Xjmpr; zerosN];
if eq == 2
    bfl       = g.bfl_disp;
    bfr       = g.bfr_disp;
    fKlb      = g.fKlb_disp;
    fKrb      = g.fKrb_disp;
    Xbl       = X  (1, :)*bfl';
    Xbr       = X  (K, :)*bfr';
    d1X(1, :) = d1X(1, :)+fKlb.*Xbl;
    d1X(K, :) = d1X(K, :)-fKrb.*Xbr;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%{
function [d1X] = d1Xlift_d(g, X, Xjmpl, Xjmpr, eq) % discontinuous
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
K             = g.numE;
R             = g.R_disp;
bD            = g.bD_disp;
fc            = g.fc_disp;
X_perm        = permute(X, [2, 3, 1]);
d1Xquad       = reshape(pagemtimes(bD, X_perm), [R, K])';
d1X           = d1Xquad*fc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
N             = g.N;
fKl           = g.fKl_disp;
fKr           = g.fKr_disp;
JMP           = zeros(K, N);
JMP(2:K  , :) = JMP(2:K  , :)-fKl(2:K  , :).*Xjmpl;
JMP(1:K-1, :) = JMP(1:K-1, :)+fKr(1:K-1, :).*Xjmpr;
switch eq
    case 1
        % DO NOTHING
    case 2
        bfl       = g.bfl_disp;
        bfr       = g.bfr_disp;
        Xbl       = X  (1, :)*bfl';
        Xbr       = X  (K, :)*bfr';
        JMP(1, :) = JMP(1, :)-fKl(1, :).*Xbl; % Ue = -Ui
        JMP(K, :) = JMP(K, :)+fKr(K, :).*Xbr; % Ue = -Ui
    otherwise
        return
end
%--------------------------------------------------------------------------
d1X           = d1X-JMP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%}