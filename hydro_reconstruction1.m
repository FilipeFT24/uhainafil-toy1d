function [z_jmp, hujmp, b_jmp, h_jmp] = ...
    hydro_reconstruction1(drytol, veltol, vellim, wetdry, zi, ze, hui, hue, bi, be)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYDRO RECONSTRUCTION:
%--------------------------------------------------------------------------
hi       = zi-bi;
he       = ze-be;
bmax     = max(bi, be);
delta    = max(0, bmax-zi);
zbtilde  = bmax-delta;
%--------------------------------------------------------------------------
h_tildei = max(0, zi-bmax);  % =NEW hi
h_tildee = max(0, ze-bmax);  % =NEW he
z_tildei = h_tildei+zbtilde; % =NEW zi = free-surface elev.
z_tildee = h_tildee+zbtilde; % =NEW ze = free-surface elev.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARS:
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------
b_jmp    = bi-zbtilde;
h_jmp    = hi-h_tildei;
h_jmp    = 1./2.*(h_tildei-h_tildee)+h_jmp;
hujmp    = 1./2.*(hutildei-hutildee);
z_jmp    = 1./2.*(z_tildei-z_tildee);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end