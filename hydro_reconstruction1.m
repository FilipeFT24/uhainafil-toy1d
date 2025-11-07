function [h_jmp, hujmp, z_jmp, zbjmp] = ...
    hydro_reconstruction1(drytol, veltol, vellim, wetdry, hui, hue, zi, ze, zbi, zbe)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYDRO RECONSTRUCTION:
hi       = zi-zbi;
he       = ze-zbe;
zbmax    = max(zbi, zbe);
delta    = max(0, zbmax-zi);
zbtilde  = zbmax-delta;
%--------------------------------------------------------------------------
h_tildei = max(0, zi-zbmax); % =NEW hi
h_tildee = max(0, ze-zbmax); % =NEW he
z_tildei = h_tildei+zbtilde; % =NEW zi = free-surface elev.
z_tildee = h_tildee+zbtilde; % =NEW ze = free-surface elev.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS:
%{
h_jmp    = 1./2.*(hi-he);
hujmp    = 1./2.*(hui-hue);
z_jmp    = 1./2.*(zi-ze);
%}
if vellim == 1 || (vellim == 2 && wetdry == 0)
    ui       = hui./hi;
    ue       = hue./he;
else
    ui       = kurganov_desingularise(hi, hui);
    ue       = kurganov_desingularise(he, hue);
    logi     = hi < drytol | hi < veltol;
    loge     = he < drytol | he < veltol;
    ui(logi) = 0;
    ue(loge) = 0;
end
hutildei = h_tildei.*ui;
hutildee = h_tildee.*ue;
h_jmp    = hi-h_tildei;
zbjmp    = zbi-zbtilde;
h_jmp    = 1./2.*(h_tildei-h_tildee)+h_jmp;
hujmp    = 1./2.*(hutildei-hutildee);
z_jmp    = 1./2.*(z_tildei-z_tildee);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end