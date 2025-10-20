function [F] = ...
    hydro_reconstruction2(drytol, veltol, rtol, G, zi, ze, hui, hue, zbi, zbe, LAMBDA, n) %#ok<INUSL> 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYDRO RECONSTRUCTION:
hi        = zi-zbi;
he        = ze-zbe;
zbmax     = max(zbi, zbe);
delta     = max(0, zbmax-zi);
zbtilde   = zbmax-delta;


a1 = zi-zbmax;
a2 = ze-zbmax;

h_tildei  = max(0, zi-zbmax); % =NEW hi
h_tildee  = max(0, ze-zbmax); % =NEW he
z_tildei  = h_tildei+zbtilde; % =NEW zi = free-surface elev.
z_tildee  = h_tildee+zbtilde; % =NEW ze = free-surface elev.
%--------------------------------------------------------------------------
%
ui        = hui./hi;
ue        = hue./he;
hutildei  = h_tildei.*ui;
hutildee  = h_tildee.*ue;
huutildei = hutildei.*ui;
huutildee = hutildee.*ue;
%{
hutildei  = h_tildei.*kurganov_desingularise(hi, hui);
hutildee  = h_tildee.*kurganov_desingularise(he, hue);
huutildei = h_tildei.*kurganov_desingularise(hi.^2, hui.^2);
huutildee = h_tildee.*kurganov_desingularise(he.^2, hue.^2);
%
%         = hutildei.*kurganov_desingularise(h_tildei, hutildei);
%         = hutildee.*kurganov_desingularise(h_tildee, hutildee);
%         = h_tildei.*kurganov_desingularise(h_tildei.^2, hutildei.^2);
%         = h_tildee.*kurganov_desingularise(h_tildee.^2, hutildee.^2);
%
hutildei (hi < drytol | hi < veltol) = 0;
hutildee (he < drytol | he < veltol) = 0;
huutildei(hi < drytol | hi < veltol) = 0;
huutildee(he < drytol | he < veltol) = 0;
%}
hutildei (hi < rtol) = 0;
hutildee (he < rtol) = 0;
huutildei(hi < rtol) = 0;
huutildee(he < rtol) = 0;
%--------------------------------------------------------------------------
gzi       = G.*z_tildei.*(1./2.*z_tildei-zbtilde);
gze       = G.*z_tildee.*(1./2.*z_tildee-zbtilde);
w0i       = z_tildei;
w0e       = z_tildee;
w1i       = hutildei;
w1e       = hutildee;
w11i      = huutildei;
w11e      = huutildee;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLUX:
HU        = 1./2.*(w1i+w1e);
GZ        = 1./2.*(gzi+gze);
W11       = 1./2.*(w11i+w11e);
dZ        =-1./2.*(w0e-w0i);
dHUx      =-1./2.*(w1e-w1i);
gzdzb     = G.*z_tildei.*(zbtilde-zbi);
Fx        = [HU, W11+GZ+gzdzb];
Fr        = [dZ, dHUx].*LAMBDA;
F         = Fx.*n+Fr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end