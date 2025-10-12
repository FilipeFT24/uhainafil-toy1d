function [Fh] = getdata1(fid, p, i)
aux     = geth0(p);
fh_data = load(fid);
data    = fh_data.data(:, :, :, i);
nf      = size(data, 1);
K       = size(data, 2);
DATA    = permute(data, [2, 3, 1]);
Fh      = zeros(nf, K);
for p = 1:nf
    Fh(p, :) = sum(pagemtimes(DATA(:, :, p), aux.Wbf'), 2)./sum(aux.W, 2);
end
%       = Fh(1:end-1, :);
end