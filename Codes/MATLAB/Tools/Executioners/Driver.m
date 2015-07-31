%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Driver Class
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB class that will hold all data for problem execution
%                   and perform all necessary solution operations. This was
%                   moved to an object environment to give more control to the
%                   solver system. The solvers used to be a simple list of
%                   functors that were beginning to lose effectiveness for
%                   large-scale problems. It became necessary to construct
%                   personalizable SI/CG/GMRES solvers for solution operations.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Input Space:    1) Data structure - contains all problem information
%                   2) Geometry Class
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Driver < handle
    % Run-Time Properties
    properties (Access = public)
        IsComplete = false
        MeshRefinementTimes
        DoFGenerationTimes
        FEGenerationTimes
    end
    % Class Objects
    properties (Access = public)
        mesh = []
        DoF  = []
        FE   = []
    end
    % Mesh/DoFHandler/FEHandler Properties
    properties (Access = public)
        % Mesh - All these are self explanatory
        Dimension
        TotalVertices
        TotalCells
        TotalFaces
        NumberMaterials
        MeshType
        OriginalMeshType
        % DoFHandler
        DoFType             % Lagrange or Serendipity
        TotalDoFs
        %FEHandler
        FEMDegree           % 1, 2, 3,...
        FEMType             % 1 = CFEM, 2 = DGFEM, 3 = WGFEM
        FEMName             % Basis Function Name   
        FEMVolumeBools
        FEMSurfaceBools
    end
    % Physical/Problem Properties
    properties (Access = public)
        ProblemType
        TransportMethod     % Transport/Diffusion
        EnergyGroups
        ScatteringOrder     % Pn truncation order
        TotalFluxMoments    % Determined from ScatteringOrder
        DiffusionXS
        TransportXS
        BCFlags
        BCValues
    end
    % Solution Properties
    properties (Access = public)
        Solution
        CellVertexNumbers
        MMSErrors
        OuterTimers
        InnerTimers
        L0_Error
        L2_Error
    end
    % Input/Output Properties
    properties (Access = public)
        PlotSolution
        SaveSolution
        SaveVTKSolution
    end
    % Solver Properties
    properties (Access = public)
        MaxIterations
        
    end
    % AMR Properties
    properties (Access = public)
        ActiveAMRBool = false   % Switches on once AMR begins (if at all)
        CurrentRefLvl = 0       % 0 corresponds to no refinement
        AMRProblem              % simple boolean switch
        AMRIrregularity         % Maximum AMR irregularity for the meshes
        RefinementLevels        % # of refinement levels for problem
        RefinementTolerances    % 0 <= tol <= 1
        ThresholdTypes          % e.g. error, # cells, etc.
        ErrorIndicators         % e.g. jump-based, Kelley, etc.
        ProjectSolution         % Bootstrapping Boolean
    end
    % MMS Properties
    properties (Access = public)
        MMSBool
        MMSQuadOrder
        ExactSolution
        MMSSourceFunction
    end
    % Private Properties for internal Driver operations
    properties (Access = private)
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                          Construction Routines
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Constructor
        function obj = Driver (varargin)
            n = nargin;
            global glob
            if isempty(glob), glob.print_info = true; end
            % Input Checking
            if n == 0
                return
            elseif n ~= 2
                error('Incorrect # of input arguments.')
            elseif n == 2
                if isstruct(varargin{1})
                    data = varargin{1};
                else
                    error('First argument is input data structure.');
                end
                if      isa(varargin{2}, 'CartesianGeometry') || ...
                        isa(varargin{2}, 'GeneralGeometry') || ...
                        isa(varargin{2}, 'AMRGeometry')
                    obj.mesh = varargin{2};
                    obj.gather_mesh_information();
                else
                    error('Second input argument needs to be a mesh class.');
                end
                clear varargin;
            end
            % Gather DoF/FE Information
            obj.DoFType   = data.Neutronics.DoFType;
            obj.FEMDegree = data.Neutronics.FEMDegree;
            obj.FEMType   = data.Neutronics.FEMType;
            obj.FEMName   = data.Neutronics.SpatialMethod;
            obj.construct_DoF_FE_Handlers();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Gather Mesh Information Routine
        function gather_mesh_information(obj)
            obj.Dimension        = obj.mesh.Dimension;
            obj.TotalVertices    = obj.mesh.TotalVertices;
            obj.TotalCells       = obj.mesh.TotalCells;
            obj.TotalFaces       = obj.mesh.TotalFaces;
            obj.NumberMaterials  = obj.mesh.TotalFaces;
            obj.MeshType         = obj.mesh.TotalFaces;
            obj.OriginalMeshType = obj.mesh.TotalFaces;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DoFHandler/FEHandler Construction Routine
        function construct_DoF_FE_Handlers(obj)
            % DoFHandler
            obj.DoF = DoFHandler(obj.mesh, obj.FEMDegree, obj.FEMType, obj.DoFType);
            obj.TotalDoFs = obj.DoF.TotalDoFs;
            % FEHandler
            obj.FE = FEHandler(obj.mesh,obj.DoF,obj.FEMName,obj.FEMVolumeBools,obj.FEMSurfaceBools,obj.MMSBool,obj.MMSQuadOrder);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                           Execution Routines
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function execute(obj)
            % Perform Initial Calculation
            obj.calculate_starting_solution();
            % Perform Mesh Refinement Operations
            obj.calculate_starting_solution();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function calculate_starting_solution(obj)
            % Exit immediately if no mesh refinements
            if ~obj.AMRProblem, return; end
            % Loop through refinement levels
            for iter = 1:obj.RefinementLevels
                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function perform_refinement_calculations(obj)
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function perform_refinement_step(obj)
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %                              AMR Routines
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

