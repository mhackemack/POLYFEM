%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outputs = calculate_eigenspectrums(data, inputs)
% Retrieve Function Handles
% -------------------------
[p_func, b_func] = get_function_calls(data);
% Allocate Memory
% ---------------
outputs = cell(data.Neutronics.Transport.NumSnLevels, inputs.TotalMeshes);
% Loop through Input Space and Calculate EigenSpectrums
% -----------------------------------------------------
nlevels = data.Neutronics.Transport.NumSnLevels;
disp('-> Computing EigenSpectrums.'); %rev_str = [];
% Loop through quadratures
for q=1:nlevels
    % Loop through meshes
    for m=1:inputs.TotalMeshes
        % Print mesh/quad combo to screen to keep from going insane...
        msg = sprintf('   -> Computing spectrum for Mesh %d of %d and Quadrature %d of %d',m,inputs.TotalMeshes,q,nlevels);
        disp(msg)
%         fprintf([rev_str,msg]);
%         rev_str = repmat(sprintf('\b'), 1, length(msg));
        % Collect input
        [m_in, phase] = combine_input_set(data, inputs, m, q);
        outputs{q,m} = p_func(b_func, m_in, phase);
    end
end
% fprintf(rev_str);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Auxilliary Function Cals
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p_func, b_func] = get_function_calls(data)
% Problem Type Function Call
if strcmp(data.Type, 'Search')
    p_func = @searcher_func;
else
    p_func = @grid_func;
end
% Matrix Build Function Call
TM = data.Neutronics.TransportMethod;
TT = data.Neutronics.Transport.transportType;
DM = data.Neutronics.DSAType;
b_func = str2func(['func_build_',TM,'_',TT,'_',DM]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = searcher_func(b_func, m_in, p_in)
pnum = p_in.TotalPhases;
pdim = p_in.NumberPhasePerDim;
dim = m_in.mesh.Dimension;
geom_dims = get_problem_dimensions(m_in.mesh);
pmin = zeros(1,dim); pmax = 2*pi*ones(1,dim)./(geom_dims(:,2)');
% Allocate memory
out.Eigen.List = zeros(pnum,1);
out.Eigen.Grid = zeros(pdim*ones(1,dim));
out.Eigen.Max = 0;
out.Search.LamList = zeros(pnum,dim);
out.Search.Evals = zeros(pnum,1);
% Loop through all phases
rev_str = [];
for p=1:pnum
    msg = sprintf('      -> Computing spectrum for Phase %d of %d',p,pnum);
    fprintf([rev_str,msg]);
    rev_str = repmat(sprintf('\b'), 1, length(msg));
    lam = p_in.WNList(p,:);
    [x,fv,~,output] = fminsearchbnd(@(x) func_search(x,b_func,m_in), lam,pmin,pmax);
    % Assign Output Values
    out.Eigen.List(p) = max(abs(fv));
    out.Search.LamList(p,:) = x;
    out.Search.Evals(p) = output.funcCount;
end
fprintf(rev_str);
% Finalize Computations
if dim ~= 1, out.Eigen.Grid = reshape(out.Eigen.List, pdim*ones(1,dim)); end
out.Eigen.Max = max(out.Eigen.List);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = grid_func(b_func, m_in, p_in)
pnum = p_in.TotalPhases;
pdim = p_in.NumberPhasePerDim;
dim = m_in.mesh.Dimension;
% Allocate memory
out.Eigen.List = zeros(pnum,1);
out.Eigen.Grid = zeros(pdim*ones(1,dim));
out.Eigen.Max = 0;
% Loop through all phases
rev_str = [];
for p=1:pnum
    msg = sprintf('      -> Computing spectrum for Phase %d of %d',p,pnum);
    fprintf([rev_str,msg]);
    rev_str = repmat(sprintf('\b'), 1, length(msg));
    lam = p_in.WNList(p,:);
    e_spect = b_func(lam, m_in);
    out.Eigen.List(p) = max(abs(eig( e_spect )));
end
fprintf(rev_str);
% Finalize Computations
if dim ~= 1, out.Eigen.Grid = reshape(out.Eigen.List, pdim*ones(1,dim)); end
out.Eigen.Max = max(out.Eigen.List);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_problem_dimensions( g )
if g.Dimension == 1
    out = [g.minX, g.maxX];
elseif g.Dimension == 2
    out = [g.minX, g.maxX; g.minY, g.maxY];
else
    out = [g.minX, g.maxX; g.minY, g.maxY; g.minZ, g.maxZ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%