function [g] = initsoldisp(g, theta, penParam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE...
%--------------------------------------------------------------------------
degreec = max(4.*g.p-2, 1);
%       = min(degreec, 12);
[Q1, W] = quadRule1D(degreec);
K       = g.numE;
N       = g.N;
R       = size (W, 2);
xydc    = zeros(K, N);
xyqc    = zeros(K, R);
for i = 1:N
    xydc(:, i) = g.detJ0T(:, 1).*g.points(i, 1)+g.coordV0T(:, 1);
end
for r = 1:R
    xyqc(:, r) = g.detJ0T(:, 1).*Q1(1, r)+g.coordV0T(:, 1);
end
%--------------------------------------------------------------------------
KN  = K*N;
NN  = N*N;
KNN = K*NN;
O   = K-1;
ONN = O*NN;
Xid = zeros(KNN, 1);
Xjd = zeros(KNN, 1);
Xio = zeros(ONN, 2);
Xjo = zeros(ONN, 2);
for o = 1:K
    for i = 1:N
        Xid((o-1)*NN+(i-1)*N+(1:N), 1) = (o-1)*N+(1:N);
        Xjd((o-1)*NN+(i-1)*N+(1:N), 1) = (o-1)*N+(i);
        if o ~= 1
            Xio((o-2)*NN+(i-1)*N+(1:N), 1) = (o-1)*N+(1:N);
            Xjo((o-2)*NN+(i-1)*N+(1:N), 1) = (o-1)*N+(i-N);
        end
        if o ~= K
            Xio((o-1)*NN+(i-1)*N+(1:N), 2) = (o-1)*N+(1:N);
            Xjo((o-1)*NN+(i-1)*N+(1:N), 2) = (o-1)*N+(i+N);
        end
    end
end
Xi = [Xid; Xio(:, 1); Xio(:, 2)];
Xj = [Xjd; Xjo(:, 1); Xjo(:, 2)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
bf        = g.basesOnQuad.phi1D    {degreec};
bd        = g.basesOnQuad.gradPhi1D{degreec};
det_perm3 = permute(g.detJ0T, [2, 3, 1]);
det_perm4 = permute(g.detJ0T, [2, 3, 4, 1]);
Wbf       = W'.*bf;
Wbd       = W'.*bd;
bD        = pagemtimes( bd, 1./det_perm3);
WbD       = pagemtimes(Wbd, 1./det_perm3);
mass      = zeros(N, N);
for j = 1:N
    mass(:, j) = Wbf'*bf(:, j);
end
Mass      = repmat(mass, [1, 1, K]);
Mass      = pagemtimes(Mass, det_perm3);
%--------------------------------------------------------------------------
d2kN      = zeros(R, N, N);
ffN       = zeros(R, N, N);
D2Kc      = zeros(N, N, R, K);
FFKc      = zeros(N, N, R, K);
F_Kc      = zeros(N, R, K);
for j = 1:N
    ffN(:, :, j) = Wbf.*bf(:, j);
end
ffR = permute(ffN, [2, 3, 1]);
for o = 1:K
    for j = 1:N
        d2kN(:, :, j) = WbD(:, :, o).*bD(:, j, o);
    end
    d2kR             = permute(d2kN, [2, 3, 1]);
    D2Kc(:, :, :, o) = d2kR;
    FFKc(:, :, :, o) = ffR;
    F_Kc(:, :, o)    = Wbf';
end
D2Kc      = pagemtimes(D2Kc, det_perm4);
FFKc      = pagemtimes(FFKc, det_perm4);
F_Kc      = pagemtimes(F_Kc, det_perm3);
%--------------------------------------------------------------------------
fc_disp   = Wbf/mass;
fx_disp   = bd'*fc_disp;
d2N       = zeros(R, N, N);
for j = 1:N
    d2N(:, :, j) = Wbd.*bd(:, j);
end
d2R       = permute(d2N, [2, 3, 1]);
d2c       = mass\sum(d2R, 3);
%{
d2kc      = zeros(N, N, K);
for o = 1:K
    d2kc(:, :, o) = Mass(:, :, o)\sum(D2Kc(:, :, :, o), 3);
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
bfl       = g.basesOnQuad.phi0D;
bfr       = g.basesOnQuad.thetaPhi0D;
bdl       = g.basesOnQuad.gradPhi0D;
bdr       = g.basesOnQuad.thetaGradPhi0D;
dl        = bdl/mass;
dr        = bdr/mass;
fl        = bfl/mass;
fr        = bfr/mass;
bDl       = repmat(bdl, [K, 1])./g.detJ0T(:, 1);
bDr       = repmat(bdr, [K, 1])./g.detJ0T(:, 1);
%--------------------------------------------------------------------------
flfl      = zeros(N, N);
flfr      = zeros(N, N);
frfl      = zeros(N, N);
frfr      = zeros(N, N);
Peni      = zeros(N, N, K  , 2);
FiDi      = zeros(N, N, K  , 2);
DiFi      = zeros(N, N, K  , 2);
Pene      = zeros(N, N, K-1, 2);
FiDe      = zeros(N, N, K-1, 2);
DiFe      = zeros(N, N, K-1, 2);
FiDibl    = zeros(N, N);
FiDibr    = zeros(N, N);
for j = 1:N
    flfl(:, j) = bfl.*bfl(:, j); % W (<-)
    flfr(:, j) = bfl.*bfr(:, j); % W (<-)
    frfl(:, j) = bfr.*bfl(:, j); % E (->)
    frfr(:, j) = bfr.*bfr(:, j); % E (->)
end
%--------------------------------------------------------------------------
% INT: W (<-)
for o = 2:K
    DiFi(:, :, o  , 1) =-1./2.*bDl(o, :)'.*bfl;
    DiFe(:, :, o-1, 1) =-1./2.*bDl(o, :)'.*bfr;
    for j = 1:N
        FiDi(:, j, o  , 1) =-1./2.*bfl.*bDl(o  , j);
        FiDe(:, j, o-1, 1) =-1./2.*bfl.*bDr(o-1, j);
    end
    Peni(:, :, o  , 1) = flfl./g.detJ0T(o, 1);
    Pene(:, :, o-1, 1) = flfr./g.detJ0T(o, 1);
end
%--------------------------------------------------------------------------
% INT: E (->)
for o = 1:K-1
    DiFi(:, :, o  , 2) = 1./2.*bDr(o, :)'.*bfr;
    DiFe(:, :, o  , 2) = 1./2.*bDr(o, :)'.*bfl;
    for j = 1:N
        FiDi(:, j, o  , 2) = 1./2.*bfr.*bDr(o  , j);
        FiDe(:, j, o  , 2) = 1./2.*bfr.*bDl(o+1, j);
    end
    Peni(:, :, o  , 2) = frfr./g.detJ0T(o, 1);
    Pene(:, :, o  , 2) = frfl./g.detJ0T(o, 1);
end
%--------------------------------------------------------------------------
% BND: W (<-) and E(->)
DiFibl    =-1./2.*bDl(1, :)'.*bfl;
DiFibr    = 1./2.*bDr(K, :)'.*bfr;
for j = 1:N
    FiDibl(:, j) =-1./2.*bfl.*bDl(1, j);
    FiDibr(:, j) = 1./2.*bfr.*bDr(K, j);
end
Penibl    = flfl./g.detJ0T(1, 1);
Penibr    = frfr./g.detJ0T(K, 1);
%--------------------------------------------------------------------------
dKl       = repmat(dl, [K, 1])./g.detJ0T.^2; % bdl./detJ0T/Mass % W (<-)
dKr       = repmat(dr, [K, 1])./g.detJ0T.^2; % bdr./detJ0T/Mass % E (->)
fKl       = repmat(fl, [K, 1])./g.detJ0T;    % bfl/Mass         % W (<-)
fKr       = repmat(fr, [K, 1])./g.detJ0T;    % bfr/Mass         % E (->)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSIGN TO...
%--------------------------------------------------------------------------
g.Xid       = Xid;
g.Xjd       = Xjd;
g.Xi        = Xi;
g.Xj        = Xj;
g.xydc_disp = xydc;
g.xyqc_disp = xyqc;
g.R_disp    = R;
g.Wbf_disp  = Wbf;
g.MASS_disp = sparse(Xid, Xjd, reshape(Mass, [], 1), KN, KN); % UNUSED
%--------------------------------------------------------------------------
g.bf_disp   = bf;
g.bD_disp   = bD;
g.D2Kc      = D2Kc;
g.FFKc      = FFKc;
g.F_Kc      = F_Kc;
g.fc_disp   = fc_disp;
g.fx_disp   = fx_disp;
%{
g.d2kc      = d2kc;
%}
g.d2c       = d2c;
%--------------------------------------------------------------------------
g.bfl_disp  = bfl;
g.bfr_disp  = bfr;
DiFi1       = DiFi(:, :, :, 1).*theta;
DiFi2       = DiFi(:, :, :, 2).*theta;
FiDi1       = FiDi(:, :, :, 1);
FiDi2       = FiDi(:, :, :, 2);
DiFibl      = DiFibl.*theta;
DiFibr      = DiFibr.*theta;
Penibl      = Penibl.*penParam;
Penibr      = Penibr.*penParam;
g.auxi1     = FiDi1-DiFi1;
g.auxi2     = FiDi2-DiFi2;
g.auxbl     = 2.*FiDibl-DiFibl-Penibl;
g.auxbr     = 2.*FiDibr-DiFibr-Penibr;
%{
g.DiFi1     = DiFi(:, :, :, 1);
g.DiFi2     = DiFi(:, :, :, 2);
g.FiDi1     = FiDi(:, :, :, 1);
g.FiDi2     = FiDi(:, :, :, 2);
g.DiFil     = DiFibl;
g.DiFir     = DiFibr;
g.FiDil     = FiDibl;
g.FiDir     = FiDibr;
g.Penil     = Penibl;
g.Penir     = Penibr;
%}
g.DiFe1     = DiFe(:, :, :, 1).*theta;
g.DiFe2     = DiFe(:, :, :, 2).*theta;
g.FiDe1     = FiDe(:, :, :, 1);
g.FiDe2     = FiDe(:, :, :, 2);
g.Peni1     = Peni(:, :, :, 1).*penParam;
g.Peni2     = Peni(:, :, :, 2).*penParam;
g.Pene1     = Pene(:, :, :, 1).*penParam;
g.Pene2     = Pene(:, :, :, 2).*penParam;
g.A         = spinit(Xi, Xj, [KN, KN]);
%--------------------------------------------------------------------------
%{
g.dKl_disp  = dKl;
g.dKr_disp  = dKr;
g.fKl_disp  = fKl;
g.fKr_disp  = fKr;
%}
g.dKli_disp = dKl(2:K  , :);
g.dKri_disp = dKr(1:K-1, :);
g.dKlb_disp = dKl(1    , :);
g.dKrb_disp = dKr(K    , :);
g.fKli_disp = fKl(2:K  , :);
g.fKri_disp = fKr(1:K-1, :);
g.fKlb_disp = fKl(1    , :);
g.fKrb_disp = fKr(K    , :);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end