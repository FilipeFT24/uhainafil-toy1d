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
            obj.t   = zeros(nf, 1);
            obj.U   = zeros(nf, K, N, V);
            %--------------------------------------------------------------
            obj.it  = 1;
            obj.dt  = g.data.tend./nf;
            obj.fid = fid;
            obj.Ph  = Ph;
            obj.isg = ismembc(test, [12, 14, 15, 16, 17]);
            obj.isp = ismembc(test, [7, 8, 9]);
            %--------------------------------------------------------------
            if obj.isg
%                 obj.bf       = g.BF;
%                 xg           = g.data.xg;
%                 ng           = size (xg, 2);
%                 obj.Ug       = zeros(nf, 1+ng);
%                 obj.Write1(t, g.x);
            end
            if obj.isp
                obj.tp       = zeros(1, 2);
                obj.vp       = zeros(1, 2);
                obj.vp(1, 1) = realmax;
                obj.vp(1, 2) = realmin;
                obj.Up       = zeros(K, N, V, 2);
                obj.Write2(t, g.x);
            end
            %--------------------------------------------------------------
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % SOL.
        function Write1(obj, t, U)
            %--------------------------------------------------------------
            obj.U(obj.it, :, :, :) = U;
            obj.t(obj.it, 1)       = t;
            obj.it                 = obj.it+1;
            %--------------------------------------------------------------
            if obj.isg
                xx = 1;
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