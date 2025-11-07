function [d1X] = d1Xlift_d(g, X, Xjmpl, Xjmpr, eq) % discontinuous
%--------------------------------------------------------------------------
% CELL
K             = g.numE;
R             = g.R_disp;
bfD           = g.bfD_disp;
X_perm        = permute(X, [2, 3, 1]);
d1Xquad       = reshape(pagemtimes(bfD, X_perm), [R, K])';
d1X           = d1Xquad*g.fc;
%--------------------------------------------------------------------------
% FACE
N             = g.N;
theta         =-1;
fl            = g.fl_disp;
fr            = g.fr_disp;
flK           = g.flKdisp;
frK           = g.frKdisp;
Xl            = X*fl';
Xr            = X*fr';
JMP           = zeros(K, N);
JMP(2:K  , :) = JMP(2:K  , :)-flK(2:K  , :).*Xjmpl;
JMP(1:K-1, :) = JMP(1:K-1, :)+frK(1:K-1, :).*Xjmpr;
switch eq
    case 1
        % DO NOTHING
    case 2
        JMP(1, :) = JMP(1, :)-flK(1, :).*Xl(1, 1); % Ue = -Ui
        JMP(K, :) = JMP(K, :)+frK(K, :).*Xr(K, 1); % Ue = -Ui
    otherwise
        return
end
%--------------------------------------------------------------------------
d1X           = d1X+theta.*JMP;
%--------------------------------------------------------------------------
end