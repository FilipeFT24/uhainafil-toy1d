function [g] = computephi(g, t, penParam) %#ok<INUSL> 
%--------------------------------------------------------------------------
alpha     = g.data.alpha;
G         = g.data.G;
K         = g.numE;
N         = g.N;
KN        = K*N;
R         = g.R_disp;
bf        = g.bf_disp;
Xi        = g.Xi;
Xj        = g.Xj;
%--------------------------------------------------------------------------
Z_dof     = g.x(:, :, 1);
Hudof     = g.x(:, :, 2);
Zbdof     = g.Zbdof;
H_dof     = Z_dof-Zbdof;
%--------------------------------------------------------------------------
fl        = g.fl_disp;
fr        = g.fr_disp;
Hul       = Hudof*fl';
Hur       = Hudof*fr';
Zl        = Z_dof*fl';
Zr        = Z_dof*fr';
Zbl       = Zbdof*fl';
Zbr       = Zbdof*fr';
[H_jmpl, Hujmpl, Z_jmpl, Zbjmpl] = ...
    hydro_reconstruction1(Hul(2:K  , 1), Hur(1:K-1, 1), Zl(2:K  , 1), Zr(1:K-1, 1), Zbl(2:K  , 1), Zbr(1:K-1, 1));
[H_jmpr, Hujmpr, Z_jmpr, Zbjmpr] = ...
    hydro_reconstruction1(Hur(1:K-1, 1), Hul(2:K  , 1), Zr(1:K-1, 1), Zl(2:K  , 1), Zbr(1:K-1, 1), Zbl(2:K  , 1));
%--------------------------------------------------------------------------
d1H_dof   = d1Xlift_d(g,   H_dof         ,   H_jmpl,   H_jmpr, 1); % eq.1: cont.
d2H_dof   = d2Xlift_d(g,   H_dof, d1H_dof,   H_jmpl,   H_jmpr, 1); % eq.1: cont.
d1Z_dof   = d1Xlift_d(g,   Z_dof         ,   Z_jmpl,   Z_jmpr, 1); % eq.1: cont.
d1Zbdof   = d1Xlift_d(g,   Zbdof         ,   Zbjmpl,   Zbjmpr, 1); % eq.1: cont.
d2Zbdof   = d2Xlift_d(g,   Zbdof, d1Zbdof,   Zbjmpl,   Zbjmpr, 1); % eq.1: cont.
%--------------------------------------------------------------------------
d1Zbl     = d1Zbdof*fl';
d1Zbr     = d1Zbdof*fr';
d1Zbjmpl  = 1./2.*(d1Zbl(2:K  , 1)-d1Zbr(1:K-1, 1));
d1Zbjmpr  = 1./2.*(d1Zbr(1:K-1, 1)-d1Zbl(2:K  , 1));
%--------------------------------------------------------------------------
d3Zbdof   = d2Xlift_d(g, d1Zbdof, d2Zbdof, d1Zbjmpl, d1Zbjmpr, 1); % eq.1: cont.
%--------------------------------------------------------------------------
d1Hudof   = d1Xlift_d(g,   Hudof         ,   Hujmpl,   Hujmpr, 2); % eq.2: mom.
d2Hudof   = d2Xlift_d(g,   Hudof, d1Hudof,   Hujmpl,   Hujmpr, 2); % eq.2: mom.
%--------------------------------------------------------------------------
H_quad    = H_dof*bf';
Huquad    = Hudof*bf';
d1H_quad  = d1H_dof*bf';
d2H_quad  = d2H_dof*bf';
d1Z_quad  = d1Z_dof*bf';
d1Zbquad  = d1Zbdof*bf';
d2Zbquad  = d2Zbdof*bf';
d3Zbquad  = d3Zbdof*bf';
d1Huquad  = d1Hudof*bf';
d2Huquad  = d2Hudof*bf';
W1quad    = H_quad.*d1Huquad-d1H_quad.*Huquad; % w1 = hq'-h'q;
W2quad    = H_quad.^2.*d2Huquad-H_quad.*(Huquad.*d2H_quad+2.*d1H_quad.*d1Huquad)+2.*Huquad.*d1H_quad.^2; % w2 = h^2q''-h(qh''+2h'q')+2h'^2q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LHS
beta            = H_quad.*(ones(K, R)+alpha.*(d1H_quad.*d1Zbquad+1./2.*H_quad.*d2Zbquad+d1Zbquad.^2));
beta_perm       = permute(beta, [3, 4, 2, 1]);
MATc            = reshape(sum(pagemtimes(g.FFKc, beta_perm), 3), [N, N, K]);
aux             = d2Ylift_d(g, H_dof, @(x) -1./3.*x.^3, penParam, MATc);  
MATd            = aux.MATd;
MATo            = aux.MATo;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RHS
GHd1Z_quad      = G.*H_quad.*d1Z_quad;
GHd1Z_quad_perm = permute(GHd1Z_quad, [2, 3, 1]);
GHd1Z           = reshape(pagemtimes(g.F_Kc, GHd1Z_quad_perm), [N, K])';
ghd1z1          = reshape(GHd1Z', [], 1);
g.GHd1ZN        = GHd1Z_quad*g.fc;
%--------------------------------------------------------------------------
Q11_quad        = 2.*H_quad.*(d1H_quad+d1Zbquad./2).*kurganov_desingularise(H_quad.^4, W1quad.^2);
%               = 2.*H_quad.*(d1H_quad+d1Zbquad./2).*(W1quad.^2./H_quad.^4);
Q12_quad        = 4./3.*H_quad.^2.*kurganov_desingularise(H_quad.^5, W1quad.*W2quad);
%               = 4./3.*H_quad.^2.*(W1quad.*W2quad./H_quad.^5);
Q13_quad        = Huquad.*d2Zbquad.*kurganov_desingularise(H_quad.^2, W1quad);
%               = Huquad.*d2Zbquad.*(W1quad./H_quad.^2);
Q14_quad        = (d1H_quad.*d2Zbquad+H_quad.*d3Zbquad./2).*kurganov_desingularise(H_quad.^2, Huquad.^2);
%               = (d1H_quad.*d2Zbquad+H_quad.*d3Zbquad./2).*(Huquad.^2./H_quad.^2);
HQ1_quad        = H_quad.*(Q11_quad+Q12_quad+Q13_quad+Q14_quad);
HQ1_quad_perm   = permute(HQ1_quad, [2, 3, 1]);
HQ1             = reshape(pagemtimes(g.F_Kc, HQ1_quad_perm), [N, K])';
hq1             = reshape(HQ1', [], 1);
g.HQ1N          = HQ1_quad*g.fc;
%--------------------------------------------------------------------------
b               =-ghd1z1./alpha-hq1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
n         = 100;
penParami = logspace(-3, 6, n);
psi_ha1   = reshape(inittype(g.inittype, g.data.PSI_H, g.BF, g.points, g.xydc, g.t)', [], 1);
ec        = zeros(K, n);
tc        = zeros(K, n);
e2        = zeros(1, n);
t2        = zeros(1, n);
cA        = zeros(1, n);
for i = 1:n
    A        = sparse(Xi, Xj, reshape(cat(3, MATd+penParami(1, i).*aux.PENd, MATo+penParami(1, i).*aux.PENo), [], 1), KN, KN);   
    psih1    = A\b;
    e1       = psi_ha1-psi_h1;
    t1       = A*psi_ha1-b;
    eN       = reshape(e1, [N, K])';
    tN       = reshape(t1, [N, K])';
    e_quad   = eN*bf';
    t_quad   = tN*bf';
    ec(:, i) = (e_quad.^2)*g.W_disp'.*g.detJ0T(:, 1);
    tc(:, i) = (t_quad.^2)*g.W_disp'.*g.detJ0T(:, 1);
    e2(1, i) = sqrt(sum(ec(:, i), 1));
    t2(1, i) = sqrt(sum(tc(:, i), 1));
    cA(1, i) = condest(A);
end
figure('Color', 'w', 'Renderer', 'painters');
subplot(1, 2, 1);
plot(penParami, cA, '-b');
set(gca, ...
    'Box', 'on', ...
    'Clipping', 'on', ...
    'Layer', 'top', ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'XScale', 'log', ...
    'YScale', 'log');
subplot(1, 2, 2);
plot(penParami, e2, '-b');
set(gca, ...
    'Box', 'on', ...
    'Clipping', 'on', ...
    'Layer', 'top', ...
    'TickLabelInterpreter', 'latex', ...
    'XMinorTick', 'on', ...
    'YMinorTick', 'on', ...
    'XScale', 'log', ...
    'YScale', 'log');
%}
%--------------------------------------------------------------------------
A          = sparse(Xi, Xj, reshape(cat(3, MATd, MATo), [], 1), KN, KN);
psi_h1     = A\b;
psi_hN     = reshape(psi_h1, [N, K])';
psi_h_quad = psi_hN*bf';
psi_quad   = psi_h_quad.*H_quad;
g.PSI_H1   = psi_h1;
g.PSI_HN   = psi_hN;
g.PSIN     = psi_quad*g.fc;
g.PHIN     = g.PSIN+g.GHd1ZN;
%{
psi_ha1    = reshape(inittype(g.inittype, g.data.PSI_H, g.BF, g.points, g.xydc, g.t)', [], 1);
e1         = psi_ha1-psi_h1;
t1         = A*psi_ha1-b;
eN         = reshape(e1, [N, K])';
tN         = reshape(t1, [N, K])';
e_quad     = eN*bf';
t_quad     = tN*bf';
ec         = (e_quad.^2)*g.W_disp'.*g.detJ0T(:, 1);
tc         = (t_quad.^2)*g.W_disp'.*g.detJ0T(:, 1);
e2         = sqrt(sum(ec, 1));
t2         = sqrt(sum(tc, 1));
cA         = condest(A);
%}
end