function [g] = initsoldisp(g)
%--------------------------------------------------------------------------
degreec = max(4.*g.p-2, 1);
degreec = min(degreec, 12);
[~, W]  = quadRule1D(degreec);
N       = g.N;
K       = g.numE;
R       = size(W, 2);
bf      = g.basesOnQuad.phi1D    {degreec};
bfd     = g.basesOnQuad.gradPhi1D{degreec};
bfl     = g.basesOnQuad.phi0D;
bfr     = g.basesOnQuad.thetaPhi0D;
gradbfl = g.basesOnQuad.gradPhi0D;
gradbfr = g.basesOnQuad.thetaGradPhi0D;
Wbf     = W'.*bf;
bfD     = zeros(R, N, K);
WbfD    = zeros(R, N, K);
bfDl    = zeros(K, N);
bfDr    = zeros(K, N);
mass    = zeros(N, N);
Mass    = zeros(N, N, K);
for j = 1:N
    mass(:, j) = Wbf'*bf(:, j);
end
%mass = eye(N);
%--------------------------------------------------------------------------
for o = 1:K
    bfD (:, :, o) = 1./g.detJ0T(o, 1).*bfd;
    WbfD(:, :, o) = W'.*bfD(:, :, o);
    bfDl(o, :)    = 1./g.detJ0T(o, 1).*gradbfl;
    bfDr(o, :)    = 1./g.detJ0T(o, 1).*gradbfr;
    Mass(:, :, o) = g.detJ0T(o, 1).*mass;
end

%--------------------------------------------------------------------------
d2         = zeros(R, N, N);
ff         = zeros(R, N, N);
D2perquadK = zeros(N , N, R, K);
FFperquadK = zeros(N , N, R, K);
F_perquadK = zeros(N ,    R, K);
for j = 1:N
    ff(:, :, j) = Wbf.*bf(:, j);
end
ff_ = permute(ff, [2, 3, 1]);
f__ = Wbf';
for o = 1:K
    for j = 1:N
        d2(:, :, j) = WbfD(:, :, o).*bfD(:, j, o);
    end
    d2_ = permute(d2, [2, 3, 1]);
    %
    D2perquadK(:, :, :, o) = d2_.*g.detJ0T(o, 1);
    FFperquadK(:, :, :, o) = ff_.*g.detJ0T(o, 1);
    F_perquadK(:,    :, o) = f__.*g.detJ0T(o, 1);
end
%--------------------------------------------------------------------------
flfl  = zeros(N, N);
flfr  = zeros(N, N);
frfl  = zeros(N, N);
frfr  = zeros(N, N);
Peni  = zeros(N, N, K  , 2);
FiDi  = zeros(N, N, K  , 2);
DiFi  = zeros(N, N, K  , 2);
Pene  = zeros(N, N, K-1, 2);
FiDe  = zeros(N, N, K-1, 2);
DiFe  = zeros(N, N, K-1, 2);

% WEST NORMAL : LEFT FACE
for j = 1:N
    flfl(:, j) = bfl.*bfl(:, j);
    flfr(:, j) = bfl.*bfr(:, j);
end
for o = 2:K
    DiFi(:, :, o  , 1) = -1./2.*bfDl(o, :)'.*bfl;
    DiFe(:, :, o-1, 1) = -1./2.*bfDl(o, :)'.*bfr;
    for j = 1:N
        FiDi(:, j, o  , 1) = -1./2.*bfl.*bfDl(o  , j);
        FiDe(:, j, o-1, 1) = -1./2.*bfl.*bfDr(o-1, j);
    end
    Peni(:, :, o  , 1) = flfl./g.detJ0T(o, 1);
    Pene(:, :, o-1, 1) = flfr./g.detJ0T(o, 1);
end
% EAST NORMAL : RIGHT FACE
for j = 1:N
    frfl(:, j) = bfr.*bfl(:, j);
    frfr(:, j) = bfr.*bfr(:, j);
end
for o = 1:K-1
    DiFi(:, :, o, 2) =  1./2.*bfDr(o, :)'.*bfr;
    DiFe(:, :, o, 2) =  1./2.*bfDr(o, :)'.*bfl;
    for j = 1:N
        FiDi(:, j, o, 2) =  1./2.*bfr.*bfDr(o  , j);
        FiDe(:, j, o, 2) =  1./2.*bfr.*bfDl(o+1, j);
    end
    Peni(:, :, o, 2) = frfr./g.detJ0T(o, 1);
    Pene(:, :, o, 2) = frfl./g.detJ0T(o, 1);
end
% BOUNDARY
DiFibl = -1./2.*bfDl(1, :)'.*bfl;
DiFibr =  1./2.*bfDr(K, :)'.*bfr;
FiDibl =  zeros(N, N);
FiDibr =  zeros(N, N);
for j = 1:N
    FiDibl(:, j) = -1./2.*bfl.*bfDl(1, j);
    FiDibr(:, j) =  1./2.*bfr.*bfDr(K, j);
end
Penibl =  flfl./g.detJ0T(1, 1);
Penibr =  frfr./g.detJ0T(K, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN  = N*N;
KN  = K*N;
KNN = KN*N;
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



dlK = zeros(K, N);
drK = zeros(K, N);
flK = zeros(K, N);
frK = zeros(K, N);
for o = 1:K
    dlK(o, :) = bfDl(o, :)/Mass(:, :, o);
    drK(o, :) = bfDr(o, :)/Mass(:, :, o);
    flK(o, :) = bfl       /Mass(:, :, o);
    frK(o, :) = bfr       /Mass(:, :, o);
end



%--------------------------------------------------------------------------
g.R_disp    = R;
g.W_disp    = W;
%--------------------------------------------------------------------------
g.bf_disp   = bf;
g.bfd_disp   = bfd;
g.bfD_disp  = bfD;
g.D2Kc      = D2perquadK;
g.FFKc      = FFperquadK;
g.F_Kc      = F_perquadK;
d2kc        = zeros(N, N, K);
for o = 1:K
    d2kc(:, :, o) = Mass(:, :, o)\sum(D2perquadK(:, :, :, o), 3);
end
g.d2kc      = d2kc;
g.fc        = Wbf/mass;
g.Wbf_disp = Wbf;

%--------------------------------------------------------------------------
g.fl_disp   = bfl;
g.fr_disp   = bfr;

g.flkdisp   = bfl/mass;
g.frkdisp   = bfr/mass;

g.flKdisp   = flK;
g.frKdisp   = frK;
g.dlKdisp   = dlK;
g.drKdisp   = drK;


g.DiFi1     = DiFi(:, :, :, 1);
g.DiFi2     = DiFi(:, :, :, 2);
g.DiFe1     = DiFe(:, :, :, 1);
g.DiFe2     = DiFe(:, :, :, 2);
g.FiDi1     = FiDi(:, :, :, 1);
g.FiDi2     = FiDi(:, :, :, 2);
g.FiDe1     = FiDe(:, :, :, 1);
g.FiDe2     = FiDe(:, :, :, 2);
g.Peni1     = Peni(:, :, :, 1);
g.Peni2     = Peni(:, :, :, 2);
g.Pene1     = Pene(:, :, :, 1);
g.Pene2     = Pene(:, :, :, 2);
g.DiFil     = DiFibl;
g.DiFir     = DiFibr;
g.FiDil     = FiDibl;
g.FiDir     = FiDibr;
g.Penil     = Penibl;
g.Penir     = Penibr;
%--------------------------------------------------------------------------
g.Xid       = Xid;
g.Xjd       = Xjd;
g.Xi        = Xi;
g.Xj        = Xj;

g.MASS_disp = sparse(Xid, Xjd, reshape(Mass, [], 1), KN, KN);
g.mass_disp = mass;
g.bfDl_disp = bfDl;
g.bfDr_disp = bfDr;
end