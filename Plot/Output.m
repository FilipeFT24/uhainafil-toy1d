classdef Output < handle
    properties
        i
        Ph
        fid
        ndt
        nf
        V
        K
        N
        n
        xc
        xv
        xg
        BF
        data
        phi
        datag
        wd
        time
        vm
        vM
        tm
        tM
        Um
        UM
    end
    methods
        function [obj] = Output(Ph, fid, ndt, nf, K, N, n, xc, xv, xg, BF)
            obj.i     = 1;
            obj.Ph    = Ph;
            obj.fid   = fid;
            obj.ndt   = ndt;
            obj.nf    = nf;
            obj.V     = 2;
            obj.K     = K;
            obj.N     = N;
            obj.n     = n;
            obj.xc    = xc;
            obj.xv    = xv;
            obj.xg    = xg;
            obj.BF    = BF;
            obj.time  = zeros(nf, 1);
            obj.data  = zeros(nf, K, N, obj.V);
            obj.phi   = zeros(nf, K, N);
            obj.datag = zeros(nf, 1+n);
            obj.wd    = zeros(nf, K);
            obj.vm    = realmax;
            obj.vM    = realmin;
            obj.tm    = 0;
            obj.tM    = 0;
            obj.Um    = zeros(K, N, obj.V);
            obj.UM    = zeros(K, N, obj.V);
        end
        function Write1(obj, t, U, wd)
            %--------------------------------------------------------------
            dry = wd == 2;
            vmi = min(U(~dry, :, 1), [], 'all');
            vMi = max(U(~dry, :, 1), [], 'all');
            %--------------------------------------------------------------
            if vmi < obj.vm
                obj.tm = t;
                obj.vm = vmi;
                obj.Um = U;
            end
            %--------------------------------------------------------------
            if vMi > obj.vM
                obj.tM = t;
                obj.vM = vMi;
                obj.UM = U;
            end
            %--------------------------------------------------------------
        end
        function Write2(obj, t, U, phi, wd)
            %--------------------------------------------------------------
            if obj.n ~= 0
                Fg       = zeros(1, obj.n+1);
                Fg(1, 1) = t;
                for k = 1:obj.n
                    xgk = obj.xg(k, 1);
                    o   = sum(xgk > obj.xv, 1);
                    xo  = (xgk-obj.xv(o, 1))./(obj.xv(o+1, 1)-obj.xv(o, 1));
                    for j = 1:obj.N
                        Fg(1, k+1) = Fg(1, k+1)+obj.BF{1, j}(xo)*U(o, j, 1);
                    end
                end
                obj.datag(obj.i, :) = Fg;
            end
            %--------------------------------------------------------------
            obj.data(obj.i, :, :, :) = U;
            obj.phi (obj.i, :, :)    = phi;
            obj.wd  (obj.i, :)       = wd;
            obj.time(obj.i, 1)       = t;
            %--------------------------------------------------------------
            obj.i = obj.i+1;
            %--------------------------------------------------------------
        end
    end
end