function data = set_PA_BC_flags(data, geom_in)
% Boundary Condition Flags
if data.problem.Dimension == 1
    data.Neutronics.Transport.BCFlags = [geom_in.xmin_bound_type,geom_in.xmax_bound_type];
elseif data.problem.Dimension == 2
    data.Neutronics.Transport.BCFlags = [geom_in.xmin_bound_type,geom_in.xmax_bound_type,geom_in.ymin_bound_type,geom_in.ymax_bound_type];
elseif data.problem.Dimension == 3
    data.Neutronics.Transport.BCFlags = [geom_in.xmin_bound_type,geom_in.xmax_bound_type,geom_in.ymin_bound_type,geom_in.ymax_bound_type,geom_in.zmin_bound_type,geom_in.zmax_bound_type];
end
% Boundary Condition Values
data.Neutronics.Transport.BCVals = cell(2*data.problem.Dimension,1);
% xmin vals
data.Neutronics.Transport.BCVals{1} = geom_in.xmin_val;
% xmax vals
data.Neutronics.Transport.BCVals{2} = geom_in.xmax_val;
if data.problem.Dimension > 1
    % ymin vals
    data.Neutronics.Transport.BCVals{3} = geom_in.ymin_val;
    % ymax vals
    data.Neutronics.Transport.BCVals{4} = geom_in.ymax_val;
end
if data.problem.Dimension > 2
    % ymin vals
    data.Neutronics.Transport.BCVals{5} = geom_in.zmin_val;
    % ymax vals
    data.Neutronics.Transport.BCVals{6} = geom_in.zmax_val;
end
