%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Parse XS File
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f_out = xs_parse_file(fname)

f_out = []; k = 1;
fid = fopen(fname);
tline = fgetl(fid);
while ischar(tline)
    f_temp = regexp(tline,'([^ ,:]*)','tokens');
    for i=1:length(f_temp)
        f_temp{i} = f_temp{i}{:};
        if is_str_numeric(f_temp{i})
            f_temp{i} = str2num(f_temp{i});
        end
    end
    f_out{k} = f_temp;
    tline = fgetl(fid);
    k = k + 1;
end
f_out = f_out';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = is_str_numeric(s)
[f_in_cell, pos] = textscan(s, '%f');
if pos == length(s)
    f = f_in_cell{:};
else
    f = [];
end
out = ~isempty(f);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%