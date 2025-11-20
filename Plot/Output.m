classdef Output < handle
    properties
        K
        N
        V
        nf
        xc
        xv
        t
        U
        it
        dt
        fid
        Ph
        isg
        isp
        bf
        ng
        xg
        Ug
        tp
        vp
        Up
    end
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = Output(g, fid, Ph)
            %--------------------------------------------------------------
            K       = g.numE;
            N       = g.N;
            D       = 1;
            V       = 1+D;
            t       = 0;
            nf      = g.data.nf;
            test    = g.test;
            obj.K   = g.numE;
            obj.N   = g.N;
            obj.V   = V;
            obj.nf  = nf;
            obj.xc  = mean(g.coordV0T, 2);
            obj.xv  = g.coordV;
            obj.t   = zeros(nf+1, 1);
            obj.U   = zeros(nf+1, K, N, V);
            %--------------------------------------------------------------
            obj.it  = 0;
            obj.dt  = g.data.tend./nf;
            obj.fid = fid;
            obj.Ph  = Ph;
            obj.isg = ismembc(test, [12, 14]);
            obj.isp = ismembc(test, [7, 8, 9, 14]);
            %--------------------------------------------------------------
            if obj.isg
                xg           = g.data.xg;
                ng           = size(xg, 1);
                obj.bf       = g.BF;
                obj.xg       = xg;
                obj.ng       = ng;
                obj.Ug       = zeros(nf, 1+ng);
            end
            if obj.isp
                obj.tp       = zeros(1, 2);
                obj.vp       = zeros(1, 2);
                obj.vp(1, 1) = realmax;
                obj.vp(1, 2) = realmin;
                obj.Up       = zeros(K, N, V, 2);
                obj.Write2(t, g.x);
            end
            obj.Write1(t, g.x);
            %--------------------------------------------------------------
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SOL.
        function Write1(obj, t, U)
            %--------------------------------------------------------------
            obj.it                 = obj.it+1;
            obj.U(obj.it, :, :, :) = U;
            obj.t(obj.it, 1)       = t;
            %--------------------------------------------------------------
            if obj.isg
                Fg       = zeros(1, obj.ng+1);
                Fg(1, 1) = t;
                for k = 1:obj.ng
                    xgk = obj.xg(k, 1);
                    o   = sum(xgk > obj.xv, 1);
                    xo  = (xgk-obj.xv(o, 1))./(obj.xv(o+1, 1)-obj.xv(o, 1));
                    for j = 1:obj.N
                        Fg(1, k+1) = Fg(1, k+1)+obj.bf{1, j}(xo)*U(o, j, 1);
                    end
                end
                obj.Ug(obj.it, :) = Fg;
            end
            %--------------------------------------------------------------
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PEAKS
        function Write2(obj, t, U)
            %--------------------------------------------------------------
            vmin = min(U(:, :, 1), [], 'all');
            vmax = max(U(:, :, 1), [], 'all');
            %--------------------------------------------------------------
            if vmin < obj.vp(1, 1)
                obj.tp(1, 1)       = t;
                obj.vp(1, 1)       = vmin;
                obj.Up(:, :, :, 1) = U;
            end
            if vmax > obj.vp(1, 2)
                obj.tp(1, 2)       = t;
                obj.vp(1, 2)       = vmax;
                obj.Up(:, :, :, 2) = U;
            end
            %--------------------------------------------------------------
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end