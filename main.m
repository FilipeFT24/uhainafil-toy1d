%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;
close all;
%--------------------------------------------------------------------------
dirs = [...
    "core/", ...
    "debug/", ...
    "init/", ...
    "integrals/", ...
    "quad/", ...
    "wetdry/", ...
    "Plot/", ...
    "Plot/S0/T0 - Convergence (Green-Nagdhi)/", ...
    "Plot/S0/T0 - Convergence (SW)/", ...
    "Plot/S0/T1 - Wave generation/", ...
    "Plot/S1/T0 - Convergence/", ...
    "Plot/S1/T1 - Lake at rest/", ...
    "Plot/S1/T2 - Drying of a lake/", ...
    "Plot/S1/T3 - Oscillating lake/", ...
    "Plot/S2/T1 - Head-on-colission/", ...
    "Plot/S2/T2 - Wall/", ...
    "Plot/S2/T3 - Overtaking colission/", ...
    "Plot/S2/T4 - Breakup of a Gaussian hump/", ...
    "Plot/S2/T5 - Dispersive dam-break/", ...
    "Plot/S3/T1 - Grilli/", ...
    "Plot/S3/T2 - Reversibility check/", ...
    "Plot/S3/T3 - Reflection of shoaling waves/", ...
    "Plot/S3/T4 - Synolakis/", ...
    "Plot/S3/T5 - Submerged bar/", ...
    "Plot/S3/T6 - Wave overtopping over a seawall/"];
for i = 1:numel(dirs)
    folder = dirs(1, i);
    if ~isfolder(folder)
        mkdir(folder);
    end
    addpath(folder);
end
%--------------------------------------------------------------------------
fn = fieldnames(get(groot, 'factory'));
fc = find(contains(fn, 'Interpreter'));
for i = 1:size(fc, 1)
    set(groot, strrep(fn{fc(i, 1)}, 'factory', 'default'), 'latex');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p         = 2;
test      = 1;
data      = setdata(test, 0);
g         = MSH(data.xv, p);
drytol    = 1.0e-02;
velcutoff = drytol;
vellim    = 1;
itype     = 0;     % 0: uhaina: interpolation
                   % 1: my:     projection
CFL       = 0.05;  % CFL./(2.*p+1);
penParam  = 1000;
run       = 1;
%--------------------------------------------------------------------------
g.CFL     = CFL;
g.t       = 0;
g.nit     = 0;
g.test    = test;
g.data    = data;
g         = initsol    (g, itype, drytol, velcutoff, vellim, data.wetdry);
g         = initsoldisp(g);
%--------------------------------------------------------------------------
if run
    obj = PlotA2(g);
    while g.t < g.data.tend
        [g, obj] = ADVt(g, obj, penParam);
    end
else
    switch test %#ok<UNRCH>
        case 1
            PlotT01C(g, [0, 0]);
        case 2
            PlotT02C(g, [0, 0]);
        case 3
            PlotT03C(g, [0, 0]);
        case 4
            PlotT04C(g, [0, 0]);
        case 5
            PlotT05C(g, [0, 0]);
        case 6
            PlotT06C(g, [0, 0]);
        case {7, 8}
            PlotT07C(g, [0, 0]);
        case 9
            PlotT09C(g, [0, 0]);
        case 10
            PlotT10C(g, [0, 0]);
        case 11
            PlotT11C(g, [0, 0]);
        case 12
            PlotT12C(g, [0, 0]);
        case 13
            PlotT13C(g, [0, 0]);
        case 14
            PlotT14C(g, [0, 0]);
        case 15
            PlotT15C(g, [0, 0]);
        case 16
        case 17
        otherwise
            return
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%