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
p         = 1;
test      = 17;
data      = setdata(test, 0);
g         = MSH(data.xv, p);
drytol    = 1.0e-02;
velcutoff = drytol; %change
vellim    = 1;
itype     = 0;     % 0: uhaina: interpolation
                   % 1: my:     projection
CFL       = 0.03;  % CFL./(2.*p+1);
penParam  = 1000;
run       = 1;
%--------------------------------------------------------------------------
g.CFL     = CFL;
g.t       = 0;
g.nit     = 0;
g.test    = test;
g.data    = data;
g         = initsol    (g, itype, drytol, velcutoff, vellim);
g         = initsoldisp(g);
%--------------------------------------------------------------------------
if run
    switch test
        case {1, 2, 3, 4, 5}
            obj = PlotA1(g);
        case {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17}
            obj = PlotA2(g);
        otherwise
            return
    end
    while g.t < g.data.tend
        [g, obj] = ADVt(g, obj, penParam);
    end
else
    switch test %#ok<UNRCH>
        case 1
            PlotT01C(g, [0, 0]); % DONE
        case 2
        case 3
            PlotT03C(g, [0, 0]);
        case 4
            PlotT04C(g, [0, 0]); % DONE
        case 5
        case {6, 7}
            PlotT06C(g, [0, 0]); % DONE
        case 8
            PlotT08C(g, [0, 0]); % DONE
        case 9
            PlotT09C(g, [0, 0]); % DONE
        case 10
            PlotT10C(g, [0, 0]); % DONE
        case 11
            PlotT11C(g, [0, 0]); % DONE
        case 12
        case 13
            PlotT13C(g, [0, 0]); % DONE
        case 14
        case 15
            PlotT15C(g, [0, 0]); % DONE
        otherwise
            return
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%