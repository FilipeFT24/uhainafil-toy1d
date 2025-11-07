function [aux] = d2Ylift_d(g, kdof, func, penParam, MATc)
%--------------------------------------------------------------------------
alpha      = g.data.alpha;
%--------------------------------------------------------------------------
% CELL
K          = g.numE;
N          = g.N;
bf         = g.bf_disp;
D2Kc       = g.D2Kc;
kquad      = kdof*bf';
kquad      = func(kquad);
kquad_perm = permute(kquad, [3, 4, 2, 1]);
MATd       =-reshape(sum(pagemtimes(D2Kc, kquad_perm), 3), [N, N, K]);
%--------------------------------------------------------------------------
% FACE
theta      =-1;
MATol      = zeros(N, N, K-1);
MATou      = zeros(N, N, K-1);
FiDi1      = g.FiDi1;
FiDi2      = g.FiDi2;
FiDe1      = g.FiDe1;
FiDe2      = g.FiDe2;
FiDil      = g.FiDil;
FiDir      = g.FiDir;
DiFi1      = g.DiFi1.*theta;
DiFi2      = g.DiFi2.*theta;
DiFe1      = g.DiFe1.*theta;
DiFe2      = g.DiFe2.*theta;
DiFil      = g.DiFil.*theta;
DiFir      = g.DiFir.*theta;
Peni1      = g.Peni1.*penParam;
Peni2      = g.Peni2.*penParam;
Pene1      = g.Pene1.*penParam;
Pene2      = g.Pene2.*penParam;
Penil      = g.Penil.*penParam;
Penir      = g.Penir.*penParam;
fl         = g.fl_disp;
fr         = g.fr_disp;
kquadl     = kdof*fl';
kquadr     = kdof*fr';
kquadl     = func(kquadl);
kquadr     = func(kquadr);
kl_perm    = permute(kquadl, [3, 2, 1]);
kr_perm    = permute(kquadr, [3, 2, 1]);
kav        = abs(1./2.*(kl_perm(1, 1, 2:K)+kr_perm(1, 1, 1:K-1)));
kavl       = cat(3, 0, kav);
kavr       = cat(3, kav, 0);
%--------------------------------------------------------------------------
MATd       = MATd+...
    pagemtimes(FiDi1-DiFi1, kl_perm)+...
    pagemtimes(FiDi2-DiFi2, kr_perm)-...
    pagemtimes(Peni1, kavl)-...
    pagemtimes(Peni2, kavr);
%--------------------------------------------------------------------------
MATol      = MATol+pagemtimes(FiDe1, kr_perm(1, 1, 1:K-1))+pagemtimes(DiFe1, kl_perm(1, 1, 2:K  ))+pagemtimes(Pene1, kav);
MATou      = MATou+pagemtimes(FiDe2, kl_perm(1, 1, 2:K  ))+pagemtimes(DiFe2, kr_perm(1, 1, 1:K-1))+pagemtimes(Pene2, kav);
%--------------------------------------------------------------------------
% 0-Dirichlet condition (p = 0): see section 4.2.1. of F. Marche (p. 303)
%                                see section 4.5.1. of FESTUNG paper.
MATd(:, :, 1) = MATd(:, :, 1)+(2.*FiDil-DiFil-Penil).*kl_perm(1, 1, 1); % sign = previous MATd!
MATd(:, :, K) = MATd(:, :, K)+(2.*FiDir-DiFir-Penir).*kr_perm(1, 1, K); % =
MATd          = MATd+MATc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if g.data.wetdry
%     % CORRECT WET/DRY (only if specified):
%     dry_              = g.WD > 1;
%     ndry              = sum(dry_, 1);
%     dryu              = dry_;
%     dryl              = dry_;
%     if ~dry_(1, 1)
%         dryl = circshift(dryl, -1);
%     end
%     MATd (:, :, dry_) = repmat(eye  (N), 1, 1, ndry);
%     MATou(:, :, dryu) = repmat(zeros(N), 1, 1, ndry);
%     MATol(:, :, dryl) = repmat(zeros(N), 1, 1, ndry);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
MATo     = cat(3, MATol, MATou);
aux.MATd = alpha.*MATd;
aux.MATo = alpha.*MATo;
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