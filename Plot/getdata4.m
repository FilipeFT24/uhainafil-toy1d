function [Fh] = getdata4(fid, p)
aux     = geth0(p);
fh_data = load(fid);
data    = fh_data.phi;
nf      = size(data, 1);
K       = size(data, 2);
DATA    = permute(data, [2, 3, 1]);
Fh      = zeros(nf-1, K);
for p = 2:nf
    Fh(p-1, :) = sum(pagemtimes(DATA(:, :, p), aux.Wbf'), 2)./sum(aux.W, 2);
end
end