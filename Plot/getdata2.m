function [Fh] = getdata2(fid, p, i, ns)
aux          = geth0(p);
data         = load(fid);
U            = data.U(:, :, :, i);
K            = size(U, 2);
DATA         = permute(U, [2, 3, 1]);
Fh           = zeros(ns, 1+K);
Fh(1, 1)     = 0;
Fh(2, 1)     = data.tp(1, 1);
Fh(3, 1)     = data.tp(1, 2);
Fh(1, 2:K+1) = sum(pagemtimes(DATA   (:, :, 1   ), aux.Wbf'), 2)./sum(aux.W, 2); % t0
Fh(2, 2:K+1) = sum(pagemtimes(data.Up(:, :, 1, 1), aux.Wbf'), 2)./sum(aux.W, 2); % tm
Fh(3, 2:K+1) = sum(pagemtimes(data.Up(:, :, 1, 2), aux.Wbf'), 2)./sum(aux.W, 2); % tM
end