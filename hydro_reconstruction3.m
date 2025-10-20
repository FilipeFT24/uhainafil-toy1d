function [h_tildei, h_tildee] = ...
    hydro_reconstruction3(zi, ze, zbi, zbe)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HYDRO RECONSTRUCTION:
zbmax    = max(zbi, zbe);
h_tildei = zi-zbmax;
h_tildee = ze-zbmax;
%--------------------------------------------------------------------------
end