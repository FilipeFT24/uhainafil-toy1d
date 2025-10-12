function [Fh] = getdatag(fid)
fh_data = load(fid);
Fh      = fh_data.datag;
end