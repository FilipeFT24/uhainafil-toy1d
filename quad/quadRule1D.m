function [Q, W] = quadRule1D(p)
%{
P        = min(floor(p./2).*2, 10);
fid_w    = join(["FACE_W", num2str(P), ".txt"],'');
fid_x    = join(["FACE_X", num2str(P), ".txt"],'');
QPointsW = readFile(fid_w);
QPointsX = readFile(fid_x);
Q        = QPointsX;
W        = QPointsW;
%}
q        = quadGaussLegendre(p, 'Domain', [0, 1]);
Q        = q.Points';
W        = q.Weights';
end
