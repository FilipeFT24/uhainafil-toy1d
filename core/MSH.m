function [g] = MSH(Xv, p)
%% MESH
g.coordV     = Xv;
g.numE       = size(g.coordV, 1)-1;
g.nume       = g.numE+1;
g.V0T        = [(1:g.numE).', (2:g.numE+1).'];
g.nuE0T      = repmat([-1, 1], g.numE, 1);
g.coordV0T   = zeros(g.numE, 2);
g.markE0TE0T = cell (1, 2);
for k = 1:2
    g.coordV0T(:, k) = g.coordV(g.V0T(:, k), :);
    g.markE0TE0T {k} = sparse(bsxfun(@eq, g.V0T(:, k), g.V0T(:, 3-k)'));
end
g.detJ0T     = g.coordV0T(:, 2)-g.coordV0T(:, 1);

%% BF
[g.dBF, g.BF, g.points] = BASIS_1D(p);
g.N           = p+1;
g.p           = p;
g.basesOnQuad = computeBasesOnQuad1D(g.BF, g.dBF, g.N, 1:12);
g.h = 1./(g.numE.*g.N);
%   = 1./(g.numE);
end