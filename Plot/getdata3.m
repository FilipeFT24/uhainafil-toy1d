function [Fh] = getdata3(fid)
fh_data = load(fid);
Fh      = fh_data.time;
%       = Fh(1:end-1, 1);
end