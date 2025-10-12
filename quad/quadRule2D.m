function [Q1, Q2, W] = quadRule2D(p)
P        = min(p, 12);
fid_w    = join(["CELL_W", num2str(P), ".txt"],'');
fid_x    = join(["CELL_X", num2str(P), ".txt"],'');
QPointsW = readFile(fid_w);
QPointsX = readFile(fid_x);
Q1       = QPointsX(:, 1)';
Q2       = QPointsX(:, 2)';
W        = QPointsW;
end