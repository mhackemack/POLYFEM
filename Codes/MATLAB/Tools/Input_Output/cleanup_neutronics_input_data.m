function data = cleanup_neutronics_input_data(data, geometry)

global glob

% General Neutronics Input
% ------------------------
if strcmp(data.Neutronics, 'transportMethod')
    error('Must specify a transport methodology - specify in field "transportMethod".')
end
if strcmp(data.Neutronics, 'FEMDegree')
    warning('FEM Degree not specified - setting default to 1.')
    data.Neutronics.FEMDegree = 1;
end
if strcmp(data.Neutronics, 'numberEnergyGroups')
    warning('Number of energy groups not specified - setting default to 1.')
    data.Neutronics.numberEnergyGroups = 1;
end
if ~isfield(data.Neutronics, 'IP_Constant')
    data.Neutronics.IP_Constant = 4.0;
end
data = get_spatial_method(data);
% Specific Transport Method Input
% -------------------------------
% Diffusion
% ------------------------------------------------------------------------------
if strcmp(data.Neutronics.transportMethod, 'Diffusion')
    if ~isfield(data.Neutronics, 'Diffusion')
        error('No "Diffusion" structure specified.')
    end
    data.Neutronics.Diffusion.Dimension = data.problem.Dimension;
    data.Neutronics.Diffusion.fluxMoments = 0;
    nm = data.problem.NumberMaterials;
    ng = data.Neutronics.numberEnergyGroups;
    % Method of Manufactured Solutions (MMS)
    if isfield(data.Neutronics.Diffusion, 'MMS')
        if data.Neutronics.Diffusion.MMS
            if data.Neutronics.numberEnergyGroups == 1
                if iscell(data.Neutronics.Diffusion.ExtSource)
                    if ~isa(data.Neutronics.Diffusion.ExtSource{1}, 'function_handle')
                        error('RHS forcing function is not a "function_handle".')
                    end
                else
                    if isa(data.Neutronics.Diffusion.ExtSource, 'function_handle')
                        tff = data.Neutronics.Diffusion.ExtSource;
                        data.Neutronics.Diffusion.ExtSource = cell(1);
                        data.Neutronics.Diffusion.ExtSource{1} = tff;
                    else
                        error('RHS forcing function is not a "function_handle".')
                    end
                end
            else
                if iscell(data.Neutronics.Diffusion.ExtSource)
                    for g=1:data.Neutronics.numberEnergyGroups
                        if ~isa(data.Neutronics.Diffusion.ExtSource{g}, 'function_handle')
                            error(['RHS forcing function for energy group ', num2str(g), ' is not a "function_handle".'])
                        end
                    end
                else
                    error('Function handles for RHS forcing function must lie in a cell structure.')
                end
            end
        end
    else
        data.Neutronics.Diffusion.MMS = false;
    end
    % Check if Multiple Iterations are Required
    % -----------------------------------------
    it_bool = false;
    for m=1:nm
        for g=1:ng
            if abs(data.Neutronics.Diffusion.FissionXS(m,g)) > glob.small
                it_bool = true;
            end
            for gg=1:ng
                if abs(data.Neutronics.Diffusion.ScatteringXS(m,g,gg)) > glob.small
                    it_bool = true;
                end
            end
        end
    end
    data.Neutronics.MultipleIterations = it_bool;
    % Implement Periodic Flux Structure
    % ---------------------------------
    data.Neutronics.Diffusion.HasPeriodicBoundary = false;
    for i=1:length(data.Neutronics.Diffusion.BCFlags)
        if data.Neutronics.Diffusion.BCFlags(i) == glob.Periodic
            data.Neutronics.Diffusion.HasPeriodicBoundary = true;
        end
    end
% Transport 
% ------------------------------------------------------------------------------
elseif strcmp(data.Neutronics.transportMethod, 'Transport')
    data.Neutronics.Transport.Dimension = data.problem.Dimension;
    if ~isfield(data.Neutronics, 'Transport')
        error('No "Transport" structure specified.')
    end
    if ~isfield(data.Neutronics.Transport, 'PnOrder')
        error('Number of flux moment levels not specified.')
    end
    if ~isfield(data.Neutronics.Transport, 'performSweeps')
        data.Neutronics.Transport.performSweeps = false;
    end
    if ~isfield(data.Neutronics.Transport, 'AngleAggregation')
        data.Neutronics.Transport.AngleAggregation = 'auto';
    end
    if ~isfield(data.Neutronics.Transport, 'transportType')
        data.Neutronics.Transport.transportType = 'upwind';
    end
    if strcmp(data.Neutronics.Transport.transportType, 'hybrid') && ...
       data.Neutronics.Transport.performSweeps
        error('Cannot perform sweeps with hybrid transport method.');
    end
    if strcmp(data.Neutronics.Transport.AngleAggregation, 'all') && ...
       data.Neutronics.Transport.performSweeps
       error('Total angle set collapsing cannot occur with sweep operations.')
    end
    if strcmp(data.Neutronics.Transport.transportType, 'hybrid')
        if ~isfield(data.Neutronics.Transport, 'FluxStabilization')
            data.Neutronics.Transport.FluxStabilization = 2.0;
        end
        if ~isfield(data.Neutronics.Transport, 'CurrentStabilization')
            data.Neutronics.Transport.CurrentStabilization = 1.0;
        end
        if ~strcmp(data.Neutronics.Transport.AngleAggregation, 'all')
            error('Full Angle collapse required for hybrid transport.')
        end
        if ~isfield(data.Neutronics.Transport, 'StabilizationMethod')
            error('Need to specify stabilization - will not attempt to guess.')
        end
        if ~strcmp(data.Neutronics.Transport.StabilizationMethod, 'EGDG') && ...
           ~strcmp(data.Neutronics.Transport.StabilizationMethod, 'LDG') && ...
           ~strcmp(data.Neutronics.Transport.StabilizationMethod, 'upwind') && ...
           ~strcmp(data.Neutronics.Transport.StabilizationMethod, 'modified')
            error('Unknown hybrid stabilization method.')
        end
        if strcmp(data.Neutronics.Transport.StabilizationMethod, 'upwind')
            data.Neutronics.Transport.StabilizationType = 0;
        elseif strcmp(data.Neutronics.Transport.StabilizationMethod, 'EGDG')
            data.Neutronics.Transport.StabilizationType = 1;
        elseif strcmp(data.Neutronics.Transport.StabilizationMethod, 'LDG')
            data.Neutronics.Transport.StabilizationType = 2;
        elseif strcmp(data.Neutronics.Transport.StabilizationMethod, 'modified')
            data.Neutronics.Transport.StabilizationType = 3;
        end
    end
    % Create Angular Quadrature Set
    % -----------------------------
    if ~isfield(data.Neutronics.Transport, 'QuadType')
        error('Quadrature type not specified.')
    end
    if strcmp(data.Neutronics.Transport.QuadType,'LS') || strcmp(data.Neutronics.Transport.QuadType,'GLC')
        if ~isfield(data.Neutronics.Transport, 'SnLevels')
            warning('Number of Sn Levels not specified - setting default to 2.')
            data.Neutronics.Transport.SnLevels = 2;
        end
    elseif strcmp(data.Neutronics.Transport.QuadType,'Manual') || strcmp(data.Neutronics.Transport.QuadType,'MANUAL') || strcmp(data.Neutronics.Transport.QuadType,'manual')
%         error('Manual quadrature sets not yet implemented.')
    else
        if ~isfield(data.Neutronics.Transport, 'PolarLevels')
            warning('Number of Polar Levels not specified - setting default to 1.')
            data.Neutronics.Transport.PolarLevels = 1;
        end
        if ~isfield(data.Neutronics.Transport, 'AzimuthalLevels')
            warning('Number of Azimuthal Levels not specified - setting default to 1.')
            data.Neutronics.Transport.AzimuthalLevels = 1;
        end
    end
    data.Neutronics.Transport = get_angular_quadrature(data.Neutronics.Transport, data.Neutronics.Transport.Dimension);
    data.Neutronics.TotalFluxMoments = data.Neutronics.Transport.TotalFluxMoments;
    % Check DSA Inputs
    % ----------------
    if ~isfield(data.Neutronics.Transport, 'performDSA')
        data.Neutronics.Transport.performDSA = 0;
    else
        if data.Neutronics.Transport.performDSA
            ntf = length(data.Neutronics.Transport.BCFlags);
            nm = data.problem.NumberMaterials;
            ng = data.Neutronics.numberEnergyGroups;
            data.Neutronics.Diffusion.fluxMoments = 0;
            data.Neutronics.Diffusion.MMS = false;
            % Allocate Cross-Sections
            data.Neutronics.Diffusion.DiffXS = zeros(nm,ng);
            data.Neutronics.Diffusion.TotalXS = zeros(nm,ng);
            data.Neutronics.Diffusion.AbsorbXS = zeros(nm,ng);
            data.Neutronics.Diffusion.ScatteringXS = zeros(nm,ng,ng);
            data.Neutronics.Diffusion.FissionXS = zeros(nm,ng);
            data.Neutronics.Diffusion.FissSpec = zeros(nm,ng);
            data.Neutronics.Diffusion.ExtSource = zeros(nm,ng);
            % Set Total Cross-Sections
            data.Neutronics.Diffusion.TotalXS = data.Neutronics.Transport.TotalXS;
            % Set Diffusion Coefficients
            for m=1:nm
                for g=1:ng
                    if data.Neutronics.Transport.fluxMoments == 0
                        data.Neutronics.Diffusion.DiffXS(m,g) = 1/(3*data.Neutronics.Transport.TotalXS(m,g));
                    else
                        data.Neutronics.Diffusion.DiffXS(m,g) = 1/(3*(data.Neutronics.Transport.TotalXS(m,g) - data.Neutronics.Transport.ScatteringXS(m,g,g,2)));
                    end
                end
            end
            % Set Scattering Cross-Sections
            for m=1:data.problem.NumberMaterials
                for g=1:data.Neutronics.numberEnergyGroups
                    for gg=1:data.Neutronics.numberEnergyGroups
                        data.Neutronics.Diffusion.ScatteringXS(m,gg,g) = data.Neutronics.Transport.ScatteringXS(m,gg,g,1);
                    end
                end
            end
            % Set Absorption Cross-Sections
            for m=1:data.problem.NumberMaterials
                for g=1:data.Neutronics.numberEnergyGroups
                    data.Neutronics.Diffusion.AbsorbXS(m,g) = data.Neutronics.Transport.TotalXS(m,g) - sum(data.Neutronics.Diffusion.ScatteringXS(m,g,:));
                end
            end
            % Set Bondary Conditions
            data.Neutronics.Diffusion.BCFlags = zeros(ntf);
            data.Neutronics.Diffusion.BCVals = zeros(ntf,ng);
            for f=1:ntf
                for g=1:ng
                    if (data.Neutronics.Transport.BCFlags(f) == glob.Vacuum || ...
                        data.Neutronics.Transport.BCFlags(f)  == glob.IncidentIsotropic || ...
                        data.Neutronics.Transport.BCFlags(f)  == glob.IncidentCurrent)

                        data.Neutronics.Diffusion.BCFlags(f) = glob.Dirichlet;
                        data.Neutronics.Diffusion.BCVals(f,g) = 0.0;
                    elseif data.Neutronics.Transport.BCFlags(f) == glob.Reflecting
                        data.Neutronics.Diffusion.BCFlags(f) = glob.Neumann;
                        data.Neutronics.Diffusion.BCVals(f,g) = 0.0;
                    elseif data.Neutronics.Transport.BCFlags(f) == glob.Periodic
                        data.Neutronics.Diffusion.BCFlags(f) = glob.Periodic;
                        data.Neutronics.Diffusion.BCVals(f,g) = 0.0;
                    end
                end
            end
        end
    end
    % Implement Reflecting/Periodic Angular Flux Structure
    % ----------------------------------------------------
    data.Neutronics.Transport.HasReflectingBoundary = false;
    data.Neutronics.Transport.HasOpposingReflectingBoundary = false;
    data.Neutronics.Transport.HasPeriodicBoundary = false;
    data.Neutronics.Transport.HasBeamBoundary = false;
    data.Neutronics.Transport.ReflectingFluxes = [];
    data.Neutronics.Transport.PeriodicFluxes = [];
    data.Neutronics.Transport.ReflectingFluxesOld = [];
    data.Neutronics.Transport.PeriodicFluxesOld = [];
    data.Neutronics.Transport.OutgoingCurrents = [];
    data.Neutronics.Transport.IncomingCurrents = [];
    data.Neutronics.Transport.OutgoingCurrentsOld = [];
    data.Neutronics.Transport.IncomingCurrentsOld = [];
    for i=1:length(data.Neutronics.Transport.BCFlags)
        if data.Neutronics.Transport.BCFlags(i) == glob.Reflecting
            data.Neutronics.Transport.HasReflectingBoundary = true;
        end
        if data.Neutronics.Transport.BCFlags(i) == glob.Periodic
            data.Neutronics.Transport.HasPeriodicBoundary = true;
        end
        if data.Neutronics.Transport.BCFlags(i) == glob.IncidentBeam
            data.Neutronics.Transport.HasBeamBoundary = true;
        end
    end
    % Method of Manufactured Solutions (MMS)
    % --------------------------------------
    if isfield(data.Neutronics.Transport, 'MMS')
        if data.Neutronics.Transport.MMS
            if iscell(data.Neutronics.Transport.ExtSource)
                for g=1:data.Neutronics.numberEnergyGroups
                    if ~isa(data.Neutronics.Transport.ExtSource{g}, 'function_handle')
                        error(['RHS forcing function for energy group ', num2str(g), ' is not a "function_handle".'])
                    end
                end
            else
                error('Function handles for RHS forcing function must lie in a cell structure.')
            end
        end
    else
        data.Neutronics.Transport.MMS = false;
    end
    % Check if Multiple Iterations are Required
    % -----------------------------------------
    it_bool = false;
    nm = data.problem.NumberMaterials;
    ng = data.Neutronics.numberEnergyGroups;
    nmom = data.Neutronics.Transport.TotalFluxMoments;
    % Check Scattering/Fission XS
    for m=1:nm
        for g=1:ng
            if abs(data.Neutronics.Transport.FissionXS(m,g)) > glob.small
                it_bool = true;
            end
            for gg=1:ng
                for k=1:nmom
                    if abs(data.Neutronics.Transport.ScatteringXS(m,g,gg,k)) > glob.small
                        it_bool = true;
                    end
                end
            end
        end
    end
    % Check Reflecting/Preiodic Boundaries
    for i=1:length(data.Neutronics.Transport.BCFlags)
        if data.Neutronics.Transport.BCFlags(i) == glob.Reflecting || ...
           data.Neutronics.Transport.BCFlags(i) == glob.Periodic
            it_bool = true;
        end
    end
    data.Neutronics.MultipleIterations = it_bool;
else
    error('Cannot determine problem type - options are "Diffusion" and "Transport".');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data] = get_spatial_method(data)
if ~isfield(data.Neutronics, 'SpatialMethod')
    % 1st Order
    if data.Neutronics.FEMDegree == 1
        warning('Spatial discretization method not specified - setting default to "PWLD".')
        data.Neutronics.SpatialMethod = 'PWLD';
    elseif data.Neutronics.FEMDegree == 2
        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%