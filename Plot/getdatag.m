function [Fh] = getdatag(fid)
data = load(fid);
Fh   = data.Ug;
end