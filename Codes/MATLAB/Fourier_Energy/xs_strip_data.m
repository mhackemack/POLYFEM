%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Strip XS Data
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = xs_strip_data(data, f_in)
disp('Retrieving XS Data.')
out.e_bounds = get_group_boundaries(data,f_in);
out.mats = get_struct_array(data); counter = 1;
% Loop through and get 1G data
for i=1:length(data.enums_1G)
    disp([' -> Stripping MT ',num2str(data.enums_1G(i))])
    out.mats(counter).MT = data.enums_1G(i);
    out.mats(counter).mat = get_1G_data(data, data.enums_1G(i), f_in);
    counter = counter + 1;
end
% Loop through and get scattering data
for i=1:length(data.scatt_enums)
    disp([' -> Stripping MT ',num2str(data.scatt_enums(i))])
    out.mats(counter).MT = data.scatt_enums(i);
    out.mats(counter).mat = zeros(data.num_groups,data.num_groups,data.iscat+1);
    for j=0:data.iscat
        out.mats(counter).mat(:,:,j+1) = get_scattering_data(data, data.scatt_enums(i), j, f_in);
    end
    counter = counter + 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Function List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_struct_array(data)
t_struct = struct('MT',0,'mat',[]);
num = length(data.scatt_enums) + length(data.enums_1G);
for i=1:num
    out(i) = t_struct;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_group_boundaries(data,f_in)
out = zeros(data.num_groups+1,1);
% Find beginning and end
beg_found = false;
beg_num = 0; end_num = 0;
for i=1:length(f_in)
    if beg_found && isempty(f_in{i})
        end_num = i-1;
        break
    elseif isempty(f_in{i})
        continue;
    end
    % Check if beginning
    if strcmpi(f_in{i}{1}, 'Group') && strcmpi(f_in{i}{2}, 'boundaries') && ~beg_found
        beg_num = i+1;
        beg_found = true;
    end
end
% Strip Energy Data
counter = 1;
for i=beg_num:end_num
    ivals = f_in{i};
    for j=1:length(ivals)
        out(counter) = ivals{j};
        counter = counter + 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_1G_data(data, MT_num, f_in)
out = zeros(data.num_groups,1);
% Find beginning and end
beg_found = false;
counter = 1;
for i=1:length(f_in)
    if ~beg_found
        if isempty(f_in{i}), continue; end
        if strcmpi(f_in{i}{1}, 'MT') && f_in{i}{2} == MT_num
            beg_found = true;
        end
    else
        if isempty(f_in{i})
            break
        elseif strcmpi(f_in{i}{1}, 'MT')
            break
        else
            ivals = f_in{i};
            for j=1:length(ivals)
                out(counter) = ivals{j};
                counter = counter + 1;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_scattering_data(data, MT_num, iscat, f_in)
out = zeros(data.num_groups);
% Find beginning and end
beg_found = false;
% sink = 0; g_first = 0; g_last = 0;
% counter = 1;
for i=1:length(f_in)
    if ~beg_found
        if isempty(f_in{i}), continue; end
        if strcmpi(f_in{i}{1}, 'MT') && f_in{i}{2} == MT_num
            if f_in{i}{4} == iscat, beg_found = true; end
        end
    else
        if isempty(f_in{i})
            break
        elseif strcmpi(f_in{i}{1}, 'MT')
            break
        elseif strcmpi(f_in{i}{1}, 'Sink')
            sink = f_in{i}{4}+1;
            g_first = f_in{i}{5}+1;
            g_last = f_in{i}{6}+1;
            dg = g_first:g_last;
            s_counter = 1;
        else
            ivals = f_in{i};
            for j=1:length(ivals)
                out(sink,dg(s_counter)) = ivals{j};
                s_counter = s_counter + 1;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%