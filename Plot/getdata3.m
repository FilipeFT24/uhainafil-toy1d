function [Fh] = getdata3(fid)
fh_data = load(fid);
Fh      = fh_data.t;
end