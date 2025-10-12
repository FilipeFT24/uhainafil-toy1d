function [aux] = geth0(p)
[Q, W]    = quadRule1D(min(max(3.*p-1, 1), 12));
[~, F, ~] = BASIS_1D(p);
N         = p+1;
bf        = zeros(size(Q, 2), N);
for j = 1:N
    bf(:, j) = F{1, j}(Q);
end
Wbf       = W'.*bf;
aux.W     = W;
aux.Wbf   = Wbf;
end