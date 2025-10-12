function [M] = readFile(filename)
fullData = readmatrix(filename, 'Delimiter', ' ');
line1    = fullData(1, :); line1 = line1(~isnan(line1));
line2    = fullData(2, :);
switch numel(line1)
    case 1
        M = line2;
    case 2
        M = reshape(line2, [line1(1), line1(2)]);
    otherwise
        return;
end
end