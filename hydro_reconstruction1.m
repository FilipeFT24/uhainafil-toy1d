function [h_jmp, hujmp, z_jmp, zbjmp] = ...
    hydro_reconstruction1(hui, hue, zi, ze, zbi, zbe)
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
hutildei = h_tildei.*kurganov_desingularise(hi, hui);
hutildee = h_tildee.*kurganov_desingularise(he, hue);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
h_jmp    = 1./2.*(hi-he);
hujmp    = 1./2.*(hui-hue);
z_jmp    = 1./2.*(zi-ze);
%}
%--------------------------------------------------------------------------
h_jmp    = hi-h_tildei;
zbjmp    = zbi-zbtilde;
%--------------------------------------------------------------------------
h_jmp    = 1./2.*(h_tildei-h_tildee)+h_jmp;
hujmp    = 1./2.*(hutildei-hutildee);
z_jmp    = 1./2.*(z_tildei-z_tildee);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end