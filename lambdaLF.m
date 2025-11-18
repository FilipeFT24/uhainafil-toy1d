function [lambda] = lambdaLF(drytol, veltol, vellim, wetdry, G, zi, ze, hui, hue, zbi, zbe)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LIMIT VELOCITY:
hi = zi-zbi;
he = ze-zbe;
if vellim == 1 || (vellim == 2 && wetdry == 0)
    ui = hui./hi;
    ue = hue./he;
else
    ui = kurganov_desingularise(hi, hui);
    ue = kurganov_desingularise(he, hue);
end
ui(hi < drytol | ui < veltol) = 0;
ue(he < drytol | ue < veltol) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYDRO RECONSTRUCTION:
zbmax    = max(zbi, zbe);
h_tildei = max(0, zi-zbmax); % =NEW hi
h_tildee = max(0, ze-zbmax); % =NEW he
lambdai  = abs(ui)+sqrt(G.*max(0, h_tildei));
lambdae  = abs(ue)+sqrt(G.*max(0, h_tildee));
lambda   = max(lambdai, lambdae);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end