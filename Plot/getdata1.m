function [Fh] = getdata1(fid, p, i)
aux    = geth0(p);
data   = load(fid);
U      = data.U(:, :, :, i);
nf     = size(U, 1);
K      = size(U, 2);
U_perm = permute(U, [2, 3, 1]);
Fh     = zeros(nf, K);
for p = 1:nf
    Fh(p, :) = sum(pagemtimes(U_perm(:, :, p), aux.Wbf'), 2)./sum(aux.W, 2);
end
end