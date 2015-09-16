function [data, geometry] = process_input_data(data, geometry)
global glob
% Process General Inputs
% ------------------------------------------------------------------------------
if isstruct(geometry) && isfield(geometry, 'geometry')
    geometry = geometry.geometry;
end
if geometry.Dimension ~= data.problem.Dimension
    error('Dimensionality does match between input and geometry.')
end
if data.problem.NumberMaterials ~= max(geometry.MatID)
    data.problem.NumberMaterials = max(geometry.MatID);
end
% Determine FEM Elementary Matrix Booleans
% ------------------------------------------------------------------------------
if strcmp(data.Neutronics.transportMethod, 'Diffusion')
    data.Neutronics.FEMVolumeBools = [1,1,0];
    if strcmp(data.Neutronics.FEMType, 'CFEM')
        data.Neutronics.FEMSurfaceBools = [1,0,0,0];
    elseif strcmp(data.Neutronics.FEMType, 'DFEM')
        data.Neutronics.FEMSurfaceBools = [1,1,0,0];
    end
elseif strcmp(data.Neutronics.transportMethod, 'Transport')
    % With DSA
    if data.Neutronics.Transport.performDSA
        data.Neutronics.FEMVolumeBools = [1,1,1];
        if strcmp(data.Neutronics.FEMType, 'CFEM')
            data.Neutronics.FEMSurfaceBools = [1,0,0,0];
        elseif strcmp(data.Neutronics.FEMType, 'DFEM')
            if strcmp(data.Neutronics.Transport.DSAType,'MIP') || ...
               strcmp(data.Neutronics.Transport.DSAType,'IP')
                data.Neutronics.FEMSurfaceBools = [1,1,0,0];
            else
                data.Neutronics.FEMSurfaceBools = [1,1,1,1];
            end
        end
    % Without DSA
    else
        data.Neutronics.FEMVolumeBools = [1,0,1];
        if strcmp(data.Neutronics.FEMType, 'CFEM')
            data.Neutronics.FEMSurfaceBools = [1,0,0,0];
        elseif strcmp(data.Neutronics.FEMType, 'DFEM')
            data.Neutronics.FEMSurfaceBools = [1,0,0,0];
        end
    end
end
% Determine DoF Type
% ------------------------------------------------------------------------------
if strcmp(lower(data.Neutronics.SpatialMethod), 'lagrange')
    data.Neutronics.DoFType = 1;
else
    data.Neutronics.DoFType = 2;
end
% Error Check DoF/FEM Types with Geometry Constraints
% ------------------------------------------------------------------------------
if data.Neutronics.DoFType == 1
    if geometry.Dimension == 2
        if strcmp(geometry.MeshType, 'Triangle')
            % do nothing
        elseif strcmp(geometry.MeshType, 'Quadrilateral')
%             if data.problem.refineMesh
%                 error('Cannot use Lagrange elements for quads with refinement since it turns cells in polygons.')
%             end
        else
            error('Lagrange elements not supported on polygons.')
        end
    elseif geometry.Dimension == 3
        if strcmp(geometry.MeshType, 'Tetrahedron')
            % do nothing
        elseif strcmp(geometry.MeshType, 'Hexahedron')
%             if data.problem.refineMesh
%                 error('Cannot use Lagrange elements for hexes with refinement since it turns cells in polyhedra.')
%             end
        else
            error('Lagrange elements not supported on polyhedra.')
        end
    end
else
    if strcmp(lower(data.Neutronics.SpatialMethod), 'serendipity')
        if strcmp(geometry.MeshType, 'Polygon') || strcmp(geometry.MeshType, 'Polyhedron')
            error('Standard Serendipity elements do not work on polygons/polyhedra.')
        end
    end
end
% Modify Global struct if necessary
% ------------------------------------------------------------------------------
if ~isfield(glob, 'print_info')
    glob.print_info = true;
end
% Check Output Argument Parameters
% ------------------------------------------------------------------------------
if ~isfield(data.problem, 'plotSolution')
    data.problem.plotSolution = 0;
end
if ~isfield(data.problem, 'saveSolution')
    data.problem.saveSolution = 0;
end
if ~isfield(data.problem, 'saveVTKSolution')
    data.problem.saveVTKSolution = 0;
end
if data.problem.saveSolution || data.problem.saveVTKSolution
    if ~isfield(data.problem,'Path')
        error('Output directory missing for solution storage.');
    end
    if ~isfield(data.problem,'Name')
        error('Output name missing for solution storage');
    end
    d_name = ['outputs/',data.problem.Path];
    if ~isequal(exist(d_name, 'dir'),7),mkdir(d_name); end
end
% ------------------------------------------------------------------------------
% Process AMR Data Structures - THIS ONE IS LAST!!!!!!!
% ------------------------------------------------------------------------------
if data.problem.refineMesh
    if ~isfield(data.problem,'refinementLevels')
        error('# of refinement levels required.')
    elseif data.problem.refinementLevels < 1
        data.problem.refineMesh = false;
        return
    end
    if ~isfield(data.problem, 'AMRIrregularity')
        data.problem.AMRIrregularity = inf;
    else
        if ~isnumeric(data.problem.AMRIrregularity)
            error('Cannot determine AMR regularity');
        else
            data.problem.AMRIrregularity = round(data.problem.AMRIrregularity);
        end
    end
    if ~isfield(data.problem, 'projectSolution')
        data.problem.projectSolution = 0;
    end
    if data.problem.refinementTolerance < 0 && abs(data.problem.refinementTolerance) > 1e-13
        error('Refinement tolerance needs to be between 0 and 1.');
    end
    if data.problem.refinementTolerance > 1 && abs(data.problem.refinementTolerance - 1) > 1e-13
        error('Refinement tolerance needs to be between 0 and 1.');
    end
    if ~isfield(data.problem, 'refinementType')
        error('Need to specify a refinement type.')
    end
end
% ------------------------------------------------------------------------------
% Process AMR Data Structures - THIS ONE IS LAST!!!!!!!
% ------------------------------------------------------------------------------

