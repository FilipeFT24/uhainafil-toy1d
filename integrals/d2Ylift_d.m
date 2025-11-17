function [aux] = d2Ylift_d(g, kdof, BETA, penParam) %#ok<INUSD> 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS #1:
%--------------------------------------------------------------------------
alpha         = g.data.alpha;
wetdry        = g.data.wetdry; %#ok<NASGU>
kdof          = kdof.*alpha;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
K             = g.numE;
N             = g.N;
bf            = g.bf_disp;
D2Kc          = g.D2Kc;
kquad         = kdof*bf';
kquad         = 1./3.*kquad.^3;
kquad_perm    = permute(kquad, [3, 4, 2, 1]);
MATd          = reshape(sum(BETA+pagemtimes(D2Kc, kquad_perm), 3), [N, N, K]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
auxi1         = g.auxi1;
auxi2         = g.auxi2;
auxbl         = g.auxbl;
auxbr         = g.auxbr;
FiDe1         = g.FiDe1;
FiDe2         = g.FiDe2;
DiFe1         = g.DiFe1;
DiFe2         = g.DiFe2;
Peni1         = g.Peni1;
Peni2         = g.Peni2;
Pene1         = g.Pene1;
Pene2         = g.Pene2;
bfl           = g.bfl_disp;
bfr           = g.bfr_disp;
kl            = kdof*bfl';
kr            = kdof*bfr';
kl            =-1./3.*kl.^3;
kr            =-1./3.*kr.^3;
kl_perm       = permute(kl, [3, 2, 1]);
kr_perm       = permute(kr, [3, 2, 1]);
kloperm       = kl_perm(1, 1, 2:K);
kroperm       = kr_perm(1, 1, 1:K-1);
kav           = abs(1./2.*(kloperm+kroperm));
kavl          = cat(3, 0, kav);
kavr          = cat(3, kav, 0);
%--------------------------------------------------------------------------
MATd          = MATd+...
    pagemtimes(auxi1, kl_perm)+...
    pagemtimes(auxi2, kr_perm)-...
    pagemtimes(Peni1, kavl)-...
    pagemtimes(Peni2, kavr);
% 0-Dirichlet condition (p = 0): see section 4.2.1. of F. Marche (p. 303)
%                                see section 4.5.1. of FESTUNG paper.
MATd(:, :, 1) = MATd(:, :, 1)+auxbl.*kl(1, 1); % sign = previous MATd!
MATd(:, :, K) = MATd(:, :, K)+auxbr.*kr(K, 1); % =
%--------------------------------------------------------------------------
MATol         = pagemtimes(FiDe1, kroperm)+pagemtimes(DiFe1, kloperm)+pagemtimes(Pene1, kav);
MATou         = pagemtimes(FiDe2, kloperm)+pagemtimes(DiFe2, kroperm)+pagemtimes(Pene2, kav);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
if wetdry
    % CORRECT WET/DRY (only if specified):
    dry_              = g.WD > 1;
    ndry              = sum(dry_, 1);
    dryu              = dry_;
    dryl              = dry_;
    if ~dry_(1, 1)
        dryl = circshift(dryl, -1);
    end
    MATd (:, :, dry_) = repmat(eye  (N), 1, 1, ndry);
    MATou(:, :, dryu) = repmat(zeros(N), 1, 1, ndry);
    MATol(:, :, dryl) = repmat(zeros(N), 1, 1, ndry);
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
MATo          = cat(3, MATol, MATou);
aux.MATd      = MATd;
aux.MATo      = MATo;
%{
%--------------------------------------------------------------------------
PENd          = zeros(N, N, K);
PENd          = PENd-...
    pagemtimes(Peni1, kavl)-...
    pagemtimes(Peni2, kavr);
PENd(:, :, 1) = PENd(:, :, 1)-Penil.*kl(1, 1));
PENd(:, :, K) = PENd(:, :, K)-Penir.*kr(K, 1));
%--------------------------------------------------------------------------
PENol         = pagemtimes(Pene1, kav);
PENou         = pagemtimes(Pene2, kav);
PENo          = cat(3, PENol, PENou);
aux.MATd      = alpha.*(MATd-PENd);
aux.MATo      = alpha.*(MATo-PENo);
aux.PENd      = alpha.*(PENd)./penParam;
aux.PENo      = alpha.*(PENo)./penParam;
%--------------------------------------------------------------------------
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end