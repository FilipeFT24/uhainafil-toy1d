function [Fh] = getdata2(fid, p, i, ns)
aux          = geth0(p);
fh_data      = load(fid);
data         = fh_data.data(:, :, :, i);
K            = size(data, 2);
DATA         = permute(data, [2, 3, 1]);
Fh           = zeros(ns, 1+K);
Fh(1, 1)     = 0;
Fh(2, 1)     = fh_data.tm;
Fh(3, 1)     = fh_data.tM;
Fh(1, 2:K+1) = sum(pagemtimes(DATA      (:, :, 1), aux.Wbf'), 2)./sum(aux.W, 2); % t0
Fh(2, 2:K+1) = sum(pagemtimes(fh_data.Um(:, :, 1), aux.Wbf'), 2)./sum(aux.W, 2); % tm
Fh(3, 2:K+1) = sum(pagemtimes(fh_data.UM(:, :, 1), aux.Wbf'), 2)./sum(aux.W, 2); % tM
end