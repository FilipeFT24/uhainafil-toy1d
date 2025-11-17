function [aux] = d2Ylift_d(g, kdof, BETA, penParam) %#ok<INUSD> 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS #1:
%--------------------------------------------------------------------------
alpha      = g.data.alpha;
wetdry     = g.data.wetdry; %#ok<NASGU>
kdof       = kdof.*alpha;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CELL:
%--------------------------------------------------------------------------
K          = g.numE;
N          = g.N;
bf         = g.bf_disp;
D2Kc       = g.D2Kc;
kquad      = kdof*bf';
kquad      = 1./3.*kquad.^3;
kquad_perm = permute(kquad, [3, 4, 2, 1]);
MATd       = reshape(sum(BETA+pagemtimes(D2Kc, kquad_perm), 3), [N, N, K]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FACE:
%--------------------------------------------------------------------------
auxi1      = g.auxi1;
auxi2      = g.auxi2;
auxbl      = g.auxbl;
auxbr      = g.auxbr;
FiDe1      = g.FiDe1;
FiDe2      = g.FiDe2;
DiFe1      = g.DiFe1;
DiFe2      = g.DiFe2;
Peni1      = g.Peni1;
Peni2      = g.Peni2;
Pene1      = g.Pene1;
Pene2      = g.Pene2;
bfl        = g.bfl_disp;
bfr        = g.bfr_disp;
kquadl     = kdof*bfl';
kquadr     = kdof*bfr';
kquadl     =-1./3.*kquadl.^3;
kquadr     =-1./3.*kquadr.^3;
kl_perm    = permute(kquadl, [3, 2, 1]);
kr_perm    = permute(kquadr, [3, 2, 1]);
kav        = abs(1./2.*(kl_perm(1, 1, 2:K)+kr_perm(1, 1, 1:K-1)));
kavl       = cat(3, 0, kav);
kavr       = cat(3, kav, 0);
%--------------------------------------------------------------------------
MATd       = MATd+...
    pagemtimes(auxi1, kl_perm)+...
    pagemtimes(auxi2, kr_perm)-...
    pagemtimes(Peni1, kavl)-...
    pagemtimes(Peni2, kavr);
MATol      = pagemtimes(FiDe1, kr_perm(1, 1, 1:K-1))+pagemtimes(DiFe1, kl_perm(1, 1, 2:K  ))+pagemtimes(Pene1, kav);
MATou      = pagemtimes(FiDe2, kl_perm(1, 1, 2:K  ))+pagemtimes(DiFe2, kr_perm(1, 1, 1:K-1))+pagemtimes(Pene2, kav);
%--------------------------------------------------------------------------
% 0-Dirichlet condition (p = 0): see section 4.2.1. of F. Marche (p. 303)
%                                see section 4.5.1. of FESTUNG paper.
MATd(:, :, 1) = MATd(:, :, 1)+auxbl.*kl_perm(1, 1, 1); % sign = previous MATd!
MATd(:, :, K) = MATd(:, :, K)+auxbr.*kr_perm(1, 1, K); % =
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
MATo     = cat(3, MATol, MATou);
aux.MATd = MATd;
aux.MATo = MATo;
%{
PENd     = zeros(N, N, K);
PENol    = zeros(N, N, K-1);
PENou    = zeros(N, N, K-1);
PENd     = PENd-...
    pagemtimes(Peni1, kavl)-...
    pagemtimes(Peni2, kavr);
PENol    = PENol+pagemtimes(Pene1, kav);
PENou    = PENou+pagemtimes(Pene2, kav);
PENo     = cat(3, PENol, PENou);
if withbnd
    PENd(:, :, 1) = PENd(:, :, 1)-pagemtimes(Penil, kl_perm(1, 1, 1));
    PENd(:, :, K) = PENd(:, :, K)-pagemtimes(Penir, kr_perm(1, 1, K));
end
aux.MATd = alpha.*(MATd-PENd);
aux.MATo = alpha.*(MATo-PENo);
aux.PENd = alpha.*(PENd)./penParam;
aux.PENo = alpha.*(PENo)./penParam;
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end