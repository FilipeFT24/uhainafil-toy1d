function [aux] = PlotA2(g)
%--------------------------------------------------------------------------
test = g.test;
switch test
    case 1
        %------------------------------------------------------------------
        aux.path        = "Plot/S0/T0 - Convergence (SW)/DATA/S00";
        aux.xlim        = [0.0, 1.0];
        aux.xtick       = 0.0:0.2:1.0;
        aux.xticklabel  = ["0.0", "0.2", "0.4", "0.6", "0.8", "1.0"];
        aux.ylim        = [0, 5];
        aux.ytick       = 0:1:5;
        aux.yticklabel  = ["0", "1", "2", "3", "4", "5"];
        aux.xtitle      = "$x$";
        aux.ytitle      = "$\zeta$";
        aux.legendlabel = ["$\zeta_{a}^{\prime}$", "$\zeta_{h}^{\prime}$"];
        aux.numcols     = 1;
        %------------------------------------------------------------------
    case 2
        %------------------------------------------------------------------
        aux.path        = "Plot/S0/T0 - Convergence (GN)/DATA/S00";
        aux.xlim        = [-50, 50];
        aux.xtick       =-50:25:50;
        aux.xticklabel  = ["-50", "-25", "0", "25", "50"];
        aux.ylim        = [-0.01, 0.25];
        aux.ytick       = [-0.01, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25];
        aux.yticklabel  = ["-0.01", "0.00", "0.05", "0.10", "0.15", "0.20", "0.25"];
        aux.xtitle      = "$x$";
        aux.ytitle      = "$\zeta^{\prime}$";
        aux.legendlabel = ["$\zeta_{a}^{\prime}$", "$\zeta_{h}^{\prime}$"];
        aux.numcols     = 1;
        %------------------------------------------------------------------
    case {7, 8}
        %------------------------------------------------------------------
        opt = g.data.opt;
        switch opt
            case 1
                str = "0.21";
            case 2
                str = "0.96";
            otherwise
                return
        end
        %------------------------------------------------------------------
        switch test
            case 7
                aux.path       = sprintf("Plot/S2/T1 - Head-on-colission/DATA/%s/Hoc", str);
                aux.xlim       = [-100, 100];
                aux.xtick      =-100:50:100;
                aux.xticklabel = ["-100", "-50", "0", "50", "100"];
            case 8
                aux.path       = sprintf("Plot/S2/T2 - Wall/DATA/%s/Wall", str);
                aux.xlim       = [-100, 0];
                aux.xtick      =-100:25:0;
                aux.xticklabel = ["-100", "-75", "-50", "-25", "0"];
            otherwise
                return
        end
        switch opt
            case 1
                aux.ylim       = [-0.05, 0.45];
                aux.ytick      = [-0.05, 0.00:0.10:0.40, 0.45];
                aux.yticklabel = ["-0.05", "0.00", "0.10", "0.20", "0.30", "0.40", "0.45"];
            case 2
                aux.ylim       = [-0.45, 2.65];
                aux.ytick      = [-0.45, 0.00:0.50:2.50, 2.65];
                aux.yticklabel = ["-0.45", "0.00", "0.50", "1.00", "1.50", "2.00", "2.50", "2.65"];
            otherwise
                return
        end
        aux.xtitle      = "$x$";
        aux.ytitle      = "";
        aux.legendlabel = ["$\zeta_{a}^{\prime}$", "$\zeta_{h}^{\prime}$"];
        aux.numcols     = 1;
        %------------------------------------------------------------------
    case 9
        %------------------------------------------------------------------
        aux.path        = "Plot/S2/T3 - Overtaking colission/DATA/Oc";
        aux.xlim        = [-200, 200];
        aux.xtick       =-200:50:200;
        aux.xticklabel  = ["-200", "-150", "-100", "-50", "0", "50", "100", "150", "200"];
        aux.ylim        = [-0.10, 1.15];
        aux.ytick       = [-0.10, 0.00:0.25:1.00, 1.15];
        aux.yticklabel  = ["-0.10", "0.00", "0.25", "0.50", "0.75", "1.00", "1.15"];
        aux.xtitle      = "$x$";
        aux.ytitle      = "$\zeta^{\prime}$";
        aux.legendlabel = ["$\zeta_{a}^{\prime}$", "$\zeta_{h}^{\prime}$"];
        aux.numcols     = 1;
        %------------------------------------------------------------------
    case 10
        %------------------------------------------------------------------
        opt = g.data.opt;
        switch opt
            case 1
                str = "01";
            case 2
                str = "20";
            otherwise
                return
        end
        aux.path       = sprintf("Plot/S2/T4 - Breakup of a Gaussian hump/DATA/%s/Bgh", str);
        aux.xlim       = [-100, 100];
        aux.xtick      =-100:50:100;
        aux.xticklabel = ["-100", "-50", "0", "50", "100"];
        %------------------------------------------------------------------
        switch g.data.opt
            case 1
                aux.ylim       = [-0.25, 1.65];
                aux.ytick      = [-0.25, 0.00:0.50:1.50, 1.65];
                aux.yticklabel = ["-0.25", "0.00", "0.50", "1.00", "1.50", "1.65"];
            case 2
                aux.ylim       = [-0.45, 2.10];
                aux.ytick      = [-0.45, 0.00:0.50:2.10];
                aux.yticklabel = ["-0.45", "0.00", "0.50", "1.00", "1.50", "2.00", "2.10"];
            otherwise
                return
        end
        aux.xtitle      = "$x$";
        aux.ytitle      = "$\zeta^{\prime}$";
        %------------------------------------------------------------------
    case 11
        %------------------------------------------------------------------
        aux.path        = "Plot/S2/T5 - Dispersive dam-break/DATA/Ddb";
        %------------------------------------------------------------------
        aux.xlim        = [-300, 300];
        aux.xtick       =-300:100:300;
        aux.xticklabel  = ["-300", "-200", "-100", "0", "100", "200", "300"];
        aux.ylim        = [0.95, 2.00];
        aux.ytick       = [0.95, 1.00:0.25:2.00];
        aux.yticklabel  = ["0.95", "1.00", "1.25", "1.50", "1.75", "2.00"];
        aux.xtitle      = "$x$";
        aux.ytitle      = "$h$";
        %------------------------------------------------------------------
    case 12
        %------------------------------------------------------------------
        bathy           = g.data.bathy;
        xv              = g.data.xv;
        h0              = g.data.h0;
        aux.xlim        = [xv(1, 1), xv(K+1, 1)]./h0;
        aux.xtick       =-50:10:30;
        aux.xticklabel  = ["-50", "-40", "-30", "-20", "-10", "0", "10", "20", "30"];
        aux.xtitle      = "$x^{\prime}$";
        %------------------------------------------------------------------
        switch g.data.opt
            case 1
                aux.path       = sprintf("Plot/S3/T1 - Grilli/DATA/0.10/C%d/Grilli", bathy);
                aux.ylim       = [-0.01, 0.20];
                aux.ytick      = [-0.01, 0.00:0.05:0.20];
                aux.yticklabel = ["-0.01", "0.00", "0.05", "0.10", "0.15", "0.20"];
                aux.ytitle     = "$\zeta^{\prime}$";
            case 2
                aux.path       = sprintf("Plot/S3/T1 - Grilli/DATA/0.15/C%d/Grilli", bathy);
                aux.ylim       = [-0.01, 0.30];
                aux.ytick      = [-0.01, 0.00:0.05:0.30];
                aux.yticklabel = ["-0.01", "0.00", "0.05", "0.10", "0.15", "0.20", "0.25", "0.30"];
                aux.ytitle     = "$\zeta^{\prime}$";
            case 3
                aux.path       = sprintf("Plot/S3/T1 - Grilli/DATA/0.20/C%d/Grilli", bathy);
                aux.ylim       = [-0.02, 0.40];
                aux.ytick      = [-0.02, 0.00:0.10:0.40];
                aux.yticklabel = ["-0.02", "0.00", "0.10", "0.20", "0.30", "0.40"];
                aux.ytitle     = "$\zeta^{\prime}$";
            case 4
                aux.path       = sprintf("Plot/S3/T1 - Grilli/DATA/0.25/C%d/Grilli", bathy);
                aux.ylim       = [-0.02, 0.50];
                aux.ytick      = [-0.02, 0.00:0.10:0.50];
                aux.yticklabel = ["-0.02", "0.00", "0.10", "0.20", "0.30", "0.40", "0.50"];
                aux.ytitle     = "$\zeta^{\prime}$";
            otherwise
                return
        end
        %------------------------------------------------------------------
    otherwise
        return
end
%--------------------------------------------------------------------------
end