function [d1X] = d1Xlift_d(g, X, Xjmpl, Xjmpr, eq) % discontinuous
%--------------------------------------------------------------------------
% CELL
K             = g.numE;
R             = g.R_disp;
bD            = g.bD_disp;
X_perm        = permute(X, [2, 3, 1]);
d1Xquad       = reshape(pagemtimes(bD, X_perm), [R, K])';
d1X           = d1Xquad*g.fc_disp;
%--------------------------------------------------------------------------
% FACE
N             = g.N;
theta         =-1;
bfl           = g.bfl_disp;
bfr           = g.bfr_disp;
fKl           = g.fKl_disp;
fKr           = g.fKr_disp;
Xl            = X*bfl';
Xr            = X*bfr';
JMP           = zeros(K, N);
JMP(2:K  , :) = JMP(2:K  , :)-fKl(2:K  , :).*Xjmpl;
JMP(1:K-1, :) = JMP(1:K-1, :)+fKr(1:K-1, :).*Xjmpr;
switch eq
    case 1
        % DO NOTHING
    case 2
        JMP(1, :) = JMP(1, :)-fKl(1, :).*Xl(1, 1); % Ue = -Ui
        JMP(K, :) = JMP(K, :)+fKr(K, :).*Xr(K, 1); % Ue = -Ui
    otherwise
        return
end
%--------------------------------------------------------------------------
d1X           = d1X+theta.*JMP;
%--------------------------------------------------------------------------
end