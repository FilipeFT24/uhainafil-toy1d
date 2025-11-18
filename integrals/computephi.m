function [g] = computephi(g, t, penParam) %#ok<INUSL> 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS #1:
%--------------------------------------------------------------------------
K        = g.numE;
N        = g.N;
R        = g.R_disp;
bf       = g.bf_disp;
bfl      = g.bfl_disp;
bfr      = g.bfr_disp;
alpha    = g.data.alpha;
G        = g.data.G;
detJ0T   = g.detJ0T;
fc       = g.fc_disp;
Wbf      = g.Wbf_disp;
drytol   = g.drytol;
veltol   = g.velcutoff;
vellim   = g.vellim;
wetdry   = g.wetdry;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS #2:
%--------------------------------------------------------------------------
N_dof    = g.x(:, :, 1);
Hudof    = g.x(:, :, 2);
B_dof    = g.zb;
H_dof    = N_dof-B_dof;
%--------------------------------------------------------------------------
Nl       = N_dof*bfl';
Nr       = N_dof*bfr';
Hul      = Hudof*bfl';
Hur      = Hudof*bfr';
Bl       = B_dof*bfl';
Br       = B_dof*bfr';
[N_jmpl, Hujmpl, B_jmpl, H_jmpl] = hydro_reconstruction1(drytol, veltol, vellim, wetdry, Nl(2:K  , 1), Nr(1:K-1, 1), Hul(2:K  , 1), Hur(1:K-1, 1), Bl(2:K  , 1), Br(1:K-1, 1));
[N_jmpr, Hujmpr, B_jmpr, H_jmpr] = hydro_reconstruction1(drytol, veltol, vellim, wetdry, Nr(1:K-1, 1), Nl(2:K  , 1), Hur(1:K-1, 1), Hul(2:K  , 1), Br(1:K-1, 1), Bl(2:K  , 1));
%--------------------------------------------------------------------------
d1H_dof  = d1Xlift_d(g, H_dof         , H_jmpl, H_jmpr, 1); % eq.1: cont.
d2H_dof  = d2Xlift_d(g, H_dof, d1H_dof, H_jmpl, H_jmpr, 1); % eq.1: cont.
d1N_dof  = d1Xlift_d(g, N_dof         , N_jmpl, N_jmpr, 1); % eq.1: cont.
d1Hudof  = d1Xlift_d(g, Hudof         , Hujmpl, Hujmpr, 2); % eq.2: mom.
d2Hudof  = d2Xlift_d(g, Hudof, d1Hudof, Hujmpl, Hujmpr, 2); % eq.2: mom.
%--------------------------------------------------------------------------
%{
d1B_dof  = d1Xlift_c(g,   B_dof);
d2B_dof  = d1Xlift_c(g, d1B_dof);
d3B_dof  = d1Xlift_c(g, d2B_dof);
%}
d1B_dof  = d1Xlift_d(g, B_dof         , B_jmpl, B_jmpr, 1);
d2B_dof  = d2Xlift_d(g, B_dof, d1B_dof, B_jmpl, B_jmpr, 1);
d1B_l    = d1B_dof*bfl';
d1B_r    = d1B_dof*bfr';
d1B_jmpl = 1./2.*(d1B_l(2:K, 1)-d1B_r(1:K-1, 1));
d1B_jmpr =-d1B_jmpl;
d3B_dof  = d2Xlift_d(g, d1B_dof, d2B_dof, d1B_jmpl, d1B_jmpr, 1);
%--------------------------------------------------------------------------
H_quad   = H_dof  *bf';
Huquad   = Hudof  *bf';
d1N_quad = d1N_dof*bf';
d1Huquad = d1Hudof*bf';
d1H_quad = d1H_dof*bf';
d1B_quad = d1B_dof*bf';
d2Huquad = d2Hudof*bf';
d2H_quad = d2H_dof*bf';
d2B_quad = d2B_dof*bf';
d3B_quad = d3B_dof*bf';
W1quad   = H_quad.*d1Huquad-d1H_quad.*Huquad; % h^2u' = hq'-h'q;
W2quad   = H_quad.^2.*d2Huquad-H_quad.*(d2H_quad.*Huquad+2.*d1H_quad.*d1Huquad)+2.*d1H_quad.^2.*Huquad; % h^3u'' = h^2q''-h(h''q+2h'q')+2h'^2q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LHS:
%--------------------------------------------------------------------------
beta            = H_quad.*(ones(K, R)+alpha.*(d1H_quad.*d1B_quad+1./2.*H_quad.*d2B_quad+d1B_quad.^2));
beta_perm       = permute(beta, [3, 4, 2, 1]);
BETA            = pagemtimes(g.FFKc, beta_perm);
aux             = d2Ylift_d(g, H_dof, BETA, penParam);  
MATd            = aux.MATd;
MATo            = aux.MATo;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RHS:
%--------------------------------------------------------------------------
GHd1N_quad      = G.*H_quad.*d1N_quad;
GHd1N           = GHd1N_quad*Wbf.*detJ0T;
%{
GHd1N_quad_perm = permute(GHd1N_quad, [2, 3, 1]);
GHd1N           = reshape(pagemtimes(g.F_Kc, GHd1N_quad_perm), [N, K])';
%}
ghd1n1          = reshape(GHd1N', [], 1);
%               = g.MASS_disp*reshape(inittype(g.itype, @(x) g.data.HYD(x, t), g.xydc_disp, g.xyqc_disp, fc)', [], 1);
g.GHd1NN        = GHd1N_quad*fc;
%--------------------------------------------------------------------------
if vellim == 1 || (vellim == 2 && wetdry == 0)
    a11 = W1quad.^2./H_quad.^4;
    a12 = W1quad.*W2quad./H_quad.^5;
    a13 = W1quad./H_quad.^2;
    a14 = Huquad.^2./H_quad.^2;
else
    a11 = kurganov_desingularise(H_quad.^4, W1quad.^2);
    a12 = kurganov_desingularise(H_quad.^5, W1quad.*W2quad);
    a13 = kurganov_desingularise(H_quad.^2, W1quad);
    a14 = kurganov_desingularise(H_quad.^2, Huquad.^2);
end
Q11_quad        = 2.*H_quad.*(d1H_quad+d1B_quad./2).*a11;
Q12_quad        = 4./3.*H_quad.^2.*a12;
Q13_quad        = Huquad.*d2B_quad.*a13;
Q14_quad        = (d1N_quad.*d2B_quad+H_quad.*d3B_quad./2).*a14;
HQ1_quad        = H_quad.*(Q11_quad+Q12_quad+Q13_quad+Q14_quad);
HQ1             = HQ1_quad*Wbf.*detJ0T;
%{
HQ1_quad_perm   = permute(HQ1_quad, [2, 3, 1]);
HQ1             = reshape(pagemtimes(g.F_Kc, HQ1_quad_perm), [N, K])';
%}
hq1             = reshape(HQ1', [], 1);
%               = g.MASS_disp*reshape(inittype(g.itype, @(x) g.data.HQ1(x, t), g.xydc_disp, g.xyqc_disp, fc)', [], 1);
g.HQ1N          = HQ1_quad*fc;
%--------------------------------------------------------------------------
b               = ghd1n1./alpha+hq1;
%               = g.MASS_disp*reshape(inittype(g.itype, @(x) g.data.RHS(x, t), g.xydc_disp, g.xyqc_disp, fc)', [], 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
val     = reshape(cat(3, MATd, MATo), [], 1);
val1    = val(:);
A       = g.A(val1);
%{
KN      = K*N;
Xi      = g.Xi;
Xj      = g.Xj;
A       = sparse(Xi, Xj, val(:), KN, KN);
%}
P1      = A\b;
PN      = reshape(P1, [N, K])';
P_quad  = PN*bf';
HP_quad = H_quad.*P_quad;
g.PN    = PN;
g.HPN   = HP_quad*fc;
g.DISP  = g.HPN-g.GHd1NN;
%{
Pa1     = reshape(inittype(g.itype, @(x) g.data.P(x, t), g.xydc_disp, g.xyqc_disp, fc)', [], 1);
e1      = Pa1-P1;
t1      = A*Pa1-b;
eN      = reshape(e1, [N, K])';
tN      = reshape(t1, [N, K])';
e_quad  = eN*bf';
t_quad  = tN*bf';
ec      = (e_quad.^2)*g.W_disp'.*detJ0T;
tc      = (t_quad.^2)*g.W_disp'.*detJ0T;
e2      = sqrt(sum(ec, 1));
t2      = sqrt(sum(tc, 1));
cA      = condest(A);
%}
%--------------------------------------------------------------------------
%{
n       = 100;
peni    = logspace(-3, 6, n);
psi_ha1 = reshape(inittype(g.itype, @(x) g.data.PSI_H(x, t), g.xydc_disp, g.xyqc_disp, fc)', [], 1);
ec      = zeros(K, n);
tc      = zeros(K, n);
e2      = zeros(1, n);
t2      = zeros(1, n);
cA      = zeros(1, n);
for i = 1:n
    A        = sparse(Xi, Xj, reshape(cat(3, MATd+peni(1, i).*aux.PENd, MATo+peni(1, i).*aux.PENo), [], 1), KN, KN);   
    psih1    = A\b;
    e1       = psi_ha1-psi_h1;
    t1       = A*psi_ha1-b;
    eN       = reshape(e1, [N, K])';
    tN       = reshape(t1, [N, K])';
    e_quad   = eN*bf';
    t_quad   = tN*bf';
    ec(:, i) = (e_quad.^2)*g.W_disp'.*detJ0T;
    tc(:, i) = (t_quad.^2)*g.W_disp'.*detJ0T;
    e2(1, i) = sqrt(sum(ec(:, i), 1));
    t2(1, i) = sqrt(sum(tc(:, i), 1));
    cA(1, i) = condest(A);
end
figure('Color', 'w', 'Renderer', 'painters');
subplot(1, 2, 1);
plot(peni, cA, '-b');
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
plot(peni, e2, '-b');
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
end