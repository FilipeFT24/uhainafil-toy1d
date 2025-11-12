function [obj] = PlotB(g, obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h0   = g.data.h0;
tend = g.data.tend;
test = g.test;
Za   = g.data.Z(obj.xc, g.t);
H_   = sum(pagemtimes(g.x(:, :, 1)-g.zb, g.Wbf'), 2)./sum(g.W, 2);
Z_   = sum(pagemtimes(g.x(:, :, 1)     , g.Wbf'), 2)./sum(g.W, 2);
Za   = Za./h0;
Z_   = Z_./h0;
%--------------------------------------------------------------------------
if mod(g.nit, 100) == 0
    if ismembc(test, [1, 2, 6, 7, 8, 9, 10])
        set(obj.Ph{1, 1}, 'YData', Za);
        set(obj.Ph{1, 2}, 'YData', Z_);
    elseif ismembc(test, [3, 4, 5, 12, 13, 14, 15, 16, 17])
        set(obj.Ph{1, 2}, 'YData', Z_);
    elseif test == 11
        set(obj.Ph{1, 2}, 'YData', H_);
    end
    drawnow limitrate;
end
%--------------------------------------------------------------------------
if ismembc(test, [1, 2])
    %----------------------------------------------------------------------
    if g.t+eps >= tend
        ua(:, :, 1) = inittype(g.itype, @(x) g.data.Z (x, g.t), g.xydc, g.xyqc, g.fi_aux);
        ua(:, :, 2) = inittype(g.itype, @(x) g.data.HU(x, g.t), g.xydc, g.xyqc, g.fi_aux);
        uh          = g.x;
        K           = g.numE;
        D           = 1;
        V           = 1+D;
        bf          = g.bf;
        ec          = zeros(K, V);
        e2          = zeros(1, V);
        for i = 1:V
            eN       = ua(:, :, i)-uh(:, :, i);
            eq       = eN*bf';
            ec(:, i) = eq*g.W';
            e2(1, i) = sqrt(sum(((eq.^2)*g.W').*g.detJ0T, 1));
        end
        save(obj.fid, 'uh', 'ec');
    end
    %----------------------------------------------------------------------
else
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
end




% obj.Write1(g.t, g.x, 0);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if g.t > (obj.i-1)*obj.ndt || any(g.t == g.data.tk, 2)
%     %----------------------------------------------------------------------
%     obj.Write2(g.t, g.x, nan, 0);
%     %----------------------------------------------------------------------
%     if obj.i > obj.nf
%         if ismembc(test, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15])
%             tm    = obj.tm;
%             tM    = obj.tM;
%             vm    = obj.vm;
%             vM    = obj.vM;
%             Um    = obj.Um;
%             UM    = obj.UM;
%             time  = obj.time;
%             data  = obj.data;
%             datag = obj.datag;
%             if obj.n ~= 0
%                 save(obj.fid, 'time', 'data', 'datag', 'tm', 'vm', 'Um', 'tM', 'vM', 'UM');
%             else
%                 save(obj.fid, 'time', 'data', 'tm', 'vm', 'Um', 'tM', 'vM', 'UM');
%             end
%         elseif ismembc(test, [12, 13])
%             time  = obj.time;
%             data  = obj.data;
%             wd    = obj.wd;
%             save(obj.fid, 'time', 'data', 'wd');
%         end
%     end
%     %----------------------------------------------------------------------
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end