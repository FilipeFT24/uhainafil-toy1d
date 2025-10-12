function [G] = Grilli2()
G       = cell(1, 6);
G{1, 1} = readmatrix('debug/Grilli/DATA2/g0.txt');
G{1, 2} = readmatrix('debug/Grilli/DATA2/g1.txt');
G{1, 3} = readmatrix('debug/Grilli/DATA2/g3.txt');
G{1, 4} = readmatrix('debug/Grilli/DATA2/g5.txt');
G{1, 5} = readmatrix('debug/Grilli/DATA2/g7.txt');
G{1, 6} = readmatrix('debug/Grilli/DATA2/g9.txt');
end