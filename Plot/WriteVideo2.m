function [] = WriteVideo2(aux, export)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if export
    v           = VideoWriter(aux.path_mp4, 'MPEG-4');
    v.FrameRate = 30;
    v.Quality   = 100;
    fid         = sprintf("%s/1.pdf", aux.path_png);
    open(v);
    exportgraphics(aux.fig, fid, 'ContentType', 'Vector', 'Resolution', 600);
    %exportgraphics(aux.fig, fid);
    %writeVideo(v, imread(fid));
    %delete(fid);
    disp(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 2:aux.nf
    %----------------------------------------------------------------------
    for j = 1:1%size(aux.Fh1, 3)
        set(aux.P1{1, j}, 'YData', aux.Fh1(i, :, j)');
        logj = aux.nk(j, :) == i;
        if any(logj, 2)
            set(aux.P2{j, logj}, 'Visible', 'on');
        end
    end
    %----------------------------------------------------------------------
    if export
        fid = sprintf("%s/%d.pdf", aux.path_png, i);
        exportgraphics(aux.fig, fid, 'ContentType', 'Vector', 'Resolution', 600);
        %exportgraphics(aux.fig, fid);
        %writeVideo(v, imread(fid));
        %delete(fid);
        disp(i);
    else
        pause(1./(aux.nf-1));
    end
    %----------------------------------------------------------------------
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if export
    close(v);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end