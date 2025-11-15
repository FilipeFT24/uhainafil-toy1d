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
d2N       = zeros(R, N, N);
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
        d2N(:, :, j) = WbD(:, :, o).*bD(:, j, o);
    end
    D2Kc(:, :, :, o) = permute(d2N, [2, 3, 1]);
    FFKc(:, :, :, o) = ffR;
    F_Kc(:, :, o)    = Wbf';
end
D2Kc      = pagemtimes(D2Kc, det_perm4);
FFKc      = pagemtimes(FFKc, det_perm4);
F_Kc      = pagemtimes(F_Kc, det_perm3);
%--------------------------------------------------------------------------
fc_disp   = Wbf/mass;
d2kc      = zeros(N, N, K);
for o = 1:K
    d2kc(:, :, o) = Mass(:, :, o)\sum(D2Kc(:, :, :, o), 3);
end
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
dKl       = repmat(dl, [K, 1])./g.detJ0T.^2; % bdl./detJ0T/Mass
dKr       = repmat(dr, [K, 1])./g.detJ0T.^2; % bdr./detJ0T/Mass
fKl       = repmat(fl, [K, 1])./g.detJ0T;    % bfl/Mass
fKr       = repmat(fr, [K, 1])./g.detJ0T;    % bfr/Mass
%--------------------------------------------------------------------------
flfl      = zeros(N, N);
flfr      = zeros(N, N);
frfl      = zeros(N, N);
frfr      = zeros(N, N);
for j = 1:N
    flfl(:, j) = bfl.*bfl(:, j); % W (<-)
    flfr(:, j) = bfl.*bfr(:, j); % W (<-)
    frfl(:, j) = bfr.*bfl(:, j); % E (->)
    frfr(:, j) = bfr.*bfr(:, j); % E (->)
end



bfDl    = zeros(K, N);
bfDr    = zeros(K, N);

%--------------------------------------------------------------------------
for o = 1:K
    bfDl(o, :)    = 1./g.detJ0T(o, 1).*bdl;
    bfDr(o, :)    = 1./g.detJ0T(o, 1).*bdr;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASSIGN TO...
%--------------------------------------------------------------------------
g.Xid       = Xid;
g.Xjd       = Xjd;
g.Xi        = Xi;
g.Xj        = Xj;
g.xydc_disp = xydc;
g.xyqc_disp = xyqc;
%--------------------------------------------------------------------------
g.R_disp    = R;
g.bf_disp   = bf;
g.bD_disp   = bD;
g.D2Kc      = D2Kc;
g.FFKc      = FFKc;
g.F_Kc      = F_Kc;
g.fc_disp   = fc_disp;
g.d2kc      = d2kc;
%--------------------------------------------------------------------------
g.bfl_disp  = bfl;
g.bfr_disp  = bfr;
g.dKl_disp  = dKl;
g.dKr_disp  = dKr;
g.fKl_disp  = fKl;
g.fKr_disp  = fKr;






%--------------------------------------------------------------------------

Peni  = zeros(N, N, K  , 2);
FiDi  = zeros(N, N, K  , 2);
DiFi  = zeros(N, N, K  , 2);
Pene  = zeros(N, N, K-1, 2);
FiDe  = zeros(N, N, K-1, 2);
DiFe  = zeros(N, N, K-1, 2);

% WEST NORMAL : LEFT FACE
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


g.MASS_disp = sparse(Xid, Xjd, reshape(Mass, [], 1), KN, KN);
g.mass_disp = mass;
end