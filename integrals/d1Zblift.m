function [d1X] = d1Zblift(g) % discontinuous
%--------------------------------------------------------------------------
X             = g.zb;
%--------------------------------------------------------------------------
% CELL
K             = g.numE;
R             = g.R;
bfD           = g.bfD;
X_perm        = permute(X, [2, 3, 1]);
d1Xquad       = reshape(pagemtimes(bfD, X_perm), [R, K])';
d1X           = d1Xquad*g.fc;
%--------------------------------------------------------------------------
% FACE
withbnd       = g.data.abslayer == 0;
N             = g.N;
theta         =-1;
fl            = g.fl_disp;
fr            = g.fr_disp;
flK           = g.flKdisp;
frK           = g.frKdisp;
Xl            = X*fl';
Xr            = X*fr';
JMP           = zeros(K, N);


% Xjmpl and Xjmpr need to be reconstructed with hydrorecons. this is still
% a flux


JMP(2:K  , :) = JMP(2:K  , :)-flK(2:K  , :).*Xjmpl;
JMP(1:K-1, :) = JMP(1:K-1, :)+frK(1:K-1, :).*Xjmpr;
if withbnd
    switch eq
        case 1
            % DO NOTHING
        case 2
            JMP(1, :) = JMP(1, :)-flK(1, :).*Xl(1, 1); % Ue = -Ui
            JMP(K, :) = JMP(K, :)+frK(K, :).*Xr(K, 1); % Ue = -Ui
        otherwise
            return
    end
end
%--------------------------------------------------------------------------
d1X           = d1X+theta.*JMP;
%--------------------------------------------------------------------------
end