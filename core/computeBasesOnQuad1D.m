function [basesOnQuad] = computeBasesOnQuad1D(bf, dbf, N, p)
n                          = numel(p);
basesOnQuad.phi1D          = cell (n, 1);
basesOnQuad.gradPhi1D      = cell (n, 1);
basesOnQuad.phi0D          = zeros(1, N);
basesOnQuad.thetaPhi0D     = zeros(1, N);
basesOnQuad.gradPhi0D      = zeros(1, N);
basesOnQuad.thetaGradPhi0D = zeros(1, N);
for it = 1:n
    [Q1, ~] = quadRule1D(it);
    basesOnQuad.phi1D    {it} = zeros(length(Q1), N);
    basesOnQuad.gradPhi1D{it} = zeros(length(Q1), N);
    for i = 1:N
        basesOnQuad.phi1D    {it}(:, i) =  bf{1, i}(Q1);
        basesOnQuad.gradPhi1D{it}(:, i) = dbf{1, i}(Q1);
    end
end
for i = 1:N
    basesOnQuad.phi0D         (1, i) =  bf{1, i}(0);
    basesOnQuad.thetaPhi0D    (1, i) =  bf{1, i}(1);
    basesOnQuad.gradPhi0D     (1, i) = dbf{1, i}(0);
    basesOnQuad.thetaGradPhi0D(1, i) = dbf{1, i}(1);
end
end