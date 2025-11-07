function [F] = ...
    hydro_reconstruction2(drytol, veltol, vellim, wetdry, G, zi, ze, hui, hue, zbi, zbe, LAMBDA, n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYDRO RECONSTRUCTION:
%--------------------------------------------------------------------------
hi        = zi-zbi;
he        = ze-zbe;
zbmax     = max(zbi, zbe);
delta     = max(0, zbmax-zi);
zbtilde   = zbmax-delta;
h_tildei  = max(0, zi-zbmax); % =NEW hi
h_tildee  = max(0, ze-zbmax); % =NEW he
z_tildei  = h_tildei+zbtilde; % =NEW zi = free-surface elev.
z_tildee  = h_tildee+zbtilde; % =NEW ze = free-surface elev.
%--------------------------------------------------------------------------
if vellim == 1 || (vellim == 2 && wetdry == 0)
    ui1       = hui./hi;
    ue1       = hue./he;
    ui2       = ui1.^2;
    ue2       = ue1.^2;
else
    ui1       = kurganov_desingularise(hi, hui);
    ue1       = kurganov_desingularise(he, hue);
    ui2       = kurganov_desingularise(hi.^2, hui.^2);
    ue2       = kurganov_desingularise(he.^2, hue.^2);
    logi      = hi < drytol | hi < veltol;
    loge      = he < drytol | he < veltol;
    ui1(logi) = 0;
    ue1(loge) = 0;
    ui2(logi) = 0;
    ue2(loge) = 0;
end
hutildei  = h_tildei.*ui1;
hutildee  = h_tildee.*ue1;
huutildei = h_tildei.*ui2;
huutildee = h_tildee.*ue2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLUX:
%--------------------------------------------------------------------------
w0i   = z_tildei;
w0e   = z_tildee;
w1i   = hutildei;
w1e   = hutildee;
w2i   = huutildei;
w2e   = huutildee;
gzi   = G.*z_tildei.*(1./2.*z_tildei-zbtilde);
gze   = G.*z_tildee.*(1./2.*z_tildee-zbtilde);
%
HU    = 1./2.*(w1i+w1e);
GZ    = 1./2.*(gzi+gze);
W2    = 1./2.*(w2i+w2e);
dZ    =-1./2.*(w0e-w0i);
dHUx  =-1./2.*(w1e-w1i);
gzdzb = G.*z_tildei.*(zbtilde-zbi);
Fx    = [HU, W2+GZ+gzdzb];
Fr    = [dZ, dHUx].*LAMBDA;
F     = Fx.*n+Fr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end