function [M] = debugrawpointer(filename, K, N, D, str)
data = readmatrix(filename, 'Delimiter', ' ');
switch str
    case "row" % = tuple format
        M  = zeros(N, D, K);
        ND = N*D;
        for o = 1:K
            for i = 1:N
                M(i, :, o) = data((o-1)*ND+(i-1)*D+1:(o-1)*ND+i*D);
            end
        end
        M  = permute(M, [3, 1, 2]);
    case "col"
        M  = zeros(K, N, D);
        ND = N*D;
        for o = 1:K
            for m = 1:D
                M(o, :, m) = data((o-1)*ND+(m-1)*N+1:(o-1)*ND+m*N);
            end
        end
    otherwise
        return
end
end