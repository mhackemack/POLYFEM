%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Solve CFEM DSA Acceleration Problem
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = perform_cfem_DSA(ndat, solvdat, mesh, DoF, FE, x, A)
global glob
nout = nargout;
ndg = DoF.TotalDoFs;
ndof = ndat.numberEnergyGroups * ndg;
if nargin < 5 || isempty(x)
    x = ones(ndof,1);
else
    x = cell_to_vector(x, DoF);
end
if length(x) > glob.maxSparse
    rhs = get_rhs(x, ndat, mesh, DoF, FE);
    [x,flag,relres,iter] = pcg(@(x) get_Ax(x, ndat, mesh, DoF, FE),rhs,solvdat.relativeTolerance,solvdat.maxIterations,[],[],x);
else
    if nargin < 6 || ~exist('A')
        [A,rhs] = get_global_matrices(x, ndat, mesh, DoF, FE);
    else
        if isempty(A)
            [A,rhs] = get_global_matrices(x, ndat, mesh, DoF, FE);
        else
            rhs = get_rhs(x, ndat, mesh, DoF, FE);
        end
    end
    % Compute Diffusion Solution
    x = A\rhs;
end
% Compute Diffusion Solution
x = A\rhs;
x = cleanup_small_vals(x);
x = vector_to_cell(x,DoF);
% Outputs
varargout{1} = x;
if nout == 2
    if exist('A', 'var')
        varargout{2} = A; 
    else
        varargout{2} = []; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              Function Listing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L,rhs] = get_global_matrices(x, ndat, mesh, DoF, FE)
global glob
n = length(x);
ng = ndat.numberEnergyGroups;
ndg = DoF.TotalDoFs;
% Allocate Memory
if n > glob.maxMatrix
    [L, rhs] = get_sparse_matrices(x, ndat, mesh, DoF, FE);
    return
else
    L = zeros(n,n);
    rhs = zeros(n,1);
end
rev_str = [];
% Loop through cells
for tcell=1:mesh.TotalCells
    cnodes = DoF.ConnectivityArray{tcell};
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    K = FE.CellStiffnessMatrix{tcell};
    % Loop through energy groups
    for g=1:ndat.numberEnergyGroups
        sg = cnodes + (g-1)*ndg;
        L(sg,sg) = L(sg,sg) + ndat.Diffusion.DiffXS(matID,g)*K + ndat.Diffusion.AbsorbXS(matID,g)*M;
        rhs(sg) = rhs(sg) + ndat.Diffusion.ScatteringXS(matID,g,g) * M * x(sg);
    end
end
% Loop through boundary faces
for f=1:mesh.TotalBoundaryFaces
    face = mesh.BoundaryFaces(f);
    fnodes = DoF.FaceCellNodes{face,1};
    flag = mesh.FaceID(face);
    M = FE.FaceMassMatrix{face,1};
    % Loop through energy groups
    fflag = ndat.Diffusion.BCFlags(flag);
    for g=1:ng
        fgnodes = fnodes + (g-1)*ndg;
        if     (ndat.Transport.BCFlags(fflag) == glob.Vacuum || ...
                ndat.Transport.BCFlags(fflag) == glob.IncidentIsotropic || ...
                ndat.Transport.BCFlags(fflag) == glob.IncidentCurrent || ...
                ndat.Transport.BCFlags(fflag) == glob.IncidentBeam)
            L(fgnodes,:) = 0;
            for i=1:length(fgnodes)
                L(fgnodes(i),fgnodes(i)) = 1;
            end
            rhs(fgnodes) = 0;
        elseif  ndat.Transport.BCFlags(fflag) == glob.Reflecting || ...
                ndat.Transport.BCFlags(fflag) == glob.Periodic
            Jin = ndat.Transport.OutgoingCurrents{f}(:,g) - ndat.Transport.OutgoingCurrentsOld{f}(:,g);
            rhs(fgnodes) = rhs(fgnodes) + M*Jin;
        end
%         switch(gflag)
%             case(glob.Dirichlet)
%                 L(fgnodes,:) = 0;
%                 for i=1:length(fgnodes)
%                     L(fgnodes(i),fgnodes(i)) = 1;
%                 end
%                 rhs(fgnodes) = ndat.Diffusion.BCVals{flag,g};
%             case(glob.Neumann)
%                 rhs(fgnodes) = rhs(fgnodes) - ndat.Diffusion.BCVals{flag,g}*F;
%             case(glob.Robin)
%                 L(fgnodes,fgnodes) = L(fgnodes,fgnodes) + 0.5*M;
%                 rhs(fgnodes) = rhs(fgnodes) + 2*ndat.Diffusion.BCVals{flag,g}*F;
%         end
    end
end
if ~issparse(L), L = sparse(L); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mat, rhs] = get_sparse_matrices(v, ndat, mesh, DoF, FE)
global glob
n = length(v);
ng = ndat.numberEnergyGroups;
ndg = DoF.TotalDoFs;

% Allocate Memory
rhs = zeros(n,1);
I = [];
J = [];
TMAT = [];

% Command Line Variables
rev_str = [];

% Loop through cells
for tcell=1:mesh.TotalCells
    % Print Current Cell Information
    if glob.print_info
        msg = sprintf('   -> Building Cell: %d of %d',tcell,mesh.TotalCells);
        fprintf([rev_str,msg]);
        rev_str = repmat(sprintf('\b'), 1, length(msg));
    end
    
    cnodes = DoF.ConnectivityArray{tcell}; ncnodes = length(cnodes);
    onesnodes = ones(ncnodes,1);
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    K = FE.CellStiffnessMatrix{tcell};
    % Loop through energy groups
    for g=1:ndat.numberEnergyGroups
        sg = cnodes + (g-1)*ndg;
        rows = onesnodes*sg;
        cols = (onesnodes*sg)';
        tmat = (ndat.Diffusion.DiffXS(matID,g)*K + ndat.Diffusion.TotalXS(matID,g)*M);
        I = [I;rows(:)];
        J = [J;cols(:)];
        TMAT = [TMAT,tmat(:)];
        tvec =  zeros(ncnodes,1);
        % Loop through right-hand-side energy groups
        for gg=1:ndat.numberEnergyGroups
            sgg = cnodes + (g-1)*ndg;
            colsgg = (onesnodes*sgg)';
            I = [I;rows(:)];
            J = [J;colsgg(:)];
            xsec = ndat.Diffusion.ScatteringXS(matID,gg,g,1);
            tmat = - xsec * M;
            TMAT = [TMAT,tmat(:)];
            tvec = tvec + xsec*M*v(sgg);
        end
        rhs(sg) = rhs(sg) + tvec;
    end
end
mat = sparse(I,J,TMAT,n,n);
dirch_nums = [];
% Loop through boundary faces
for f=1:mesh.TotalBoundaryFaces
    % Print Current Face Information
    if glob.print_info
        msg = sprintf('   -> Building boundary face: %d of %d',f,mesh.TotalBoundaryFaces);
        fprintf([rev_str,msg]);
        rev_str = repmat(sprintf('\b'), 1, length(msg));
    end
    
    face = mesh.BoundaryFaces(f);
    fnodes = DoF.getFaceCellNodes(face,1);
    flag = mesh.FaceID(face);
    M = FE.FaceMassMatrix{face,1};
    F = FE.FaceFunctionMatrix{face,1};
    % Loop through energy groups
    gflag = ndat.Diffusion.BCFlags(flag);
    for g=1:ng
        fgnodes = fnodes + (g-1)*ndg;
        switch(gflag)
            case(glob.Dirichlet)
                dirch_nums = [dirch_nums,fgnodes];
                rhs(fgnodes) = ndat.Diffusion.BCVals{flag,g}*ones(length(fgnodes),1);
            case(glob.Neumann)
                rhs(fgnodes) = rhs(fgnodes) - ndat.Diffusion.BCVals{flag,g}*F;
            case(glob.Robin)
                mat(fgnodes,fgnodes) = mat(fgnodes,fgnodes) + 0.5*M;
                rhs(fgnodes) = rhs(fgnodes) + 2*ndat.Diffusion.BCVals{flag,g}*F;
        end
    end
end
if ~isempty(dirch_nums)
    dirch_nums = unique(dirch_nums);
    mat(dirch_nums,:) = 0;
    diag_id = (dirch_nums-1)*(n+1)+1;
    mat(diag_id) = 1;
end
if glob.print_info, fprintf(rev_str); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhs = get_rhs(v, ndat, mesh, DoF, FE)
global glob
n = length(v);
ng = ndat.numberEnergyGroups;
ndg = DoF.TotalDoFs;
% Allocate Memory
rhs = zeros(n,1);
rev_str = [];
% Loop through cells
for tcell=1:mesh.TotalCells
    % Print Current Cell Information
    if glob.print_info
        msg = sprintf('   -> Building Cell: %d of %d',tcell,mesh.TotalCells);
        fprintf([rev_str,msg]);
        rev_str = repmat(sprintf('\b'), 1, length(msg));
    end
    
    cnodes = DoF.ConnectivityArray{tcell}; ncnodes = length(cnodes);
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    % Loop through energy groups
    for g=1:ndat.numberEnergyGroups
        tvec =  zeros(ncnodes,1);
        % Loop through right-hand-side energy groups
        for gg=1:ndat.numberEnergyGroups
            sgg = cnodes + (g-1)*ndg;
            xsec = ndat.Diffusion.ScatteringXS(matID,gg,g,1);
            tvec = tvec + xsec*M*v(sgg);
        end
        sg = cnodes + (g-1)*ndg;
        rhs(sg) = rhs(sg) + tvec;
    end
end
% Loop through boundary faces
for f=1:mesh.TotalBoundaryFaces
    % Print Current Face Information
    if glob.print_info
        msg = sprintf('   -> Building boundary face: %d of %d',f,mesh.TotalBoundaryFaces);
        fprintf([rev_str,msg]);
        rev_str = repmat(sprintf('\b'), 1, length(msg));
    end
    
    face = mesh.BoundaryFaces(f);
    fnodes = DoF.getFaceCellNodes(face,1);
    flag = mesh.FaceID(face);
    F = FE.FaceFunctionMatrix{face,1};
    % Loop through energy groups
    for g=1:ng
        fgnodes = fnodes + (g-1)*ndg;
        gflag = ndat.Diffusion.BCFlags(flag,g);
        switch(gflag)
            case(glob.Dirichlet)
                rhs(fgnodes) = ndat.Diffusion.BCVals{flag,g}*ones(length(fgnodes),1);
            case(glob.Neumann)
                rhs(fgnodes) = rhs(fgnodes) - ndat.Diffusion.BCVals{flag,g}*F;
            case(glob.Robin)
                rhs(fgnodes) = rhs(fgnodes) + 2*ndat.Diffusion.BCVals{flag,g}*F;
        end
    end
end
if glob.print_info, fprintf(rev_str); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_Ax(v, ndat, mesh, DoF, FE)
global glob
n = length(v);
ng = ndat.numberEnergyGroups;
ndg = DoF.TotalDoFs;
% Allocate Memory
out = zeros(n,1);
% Loop through cells
for tcell=1:mesh.TotalCells
    cnodes = DoF.ConnectivityArray{tcell};
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    K = FE.CellStiffnessMatrix{tcell};
    % Loop through energy groups
    for g=1:ndat.numberEnergyGroups
        sg = cnodes + (g-1)*ndg;
        out(sg) = out(sg) + (ndat.Diffusion.DiffXS(matID,g)*K + ndat.Diffusion.TotalXS(matID,g)*M)*x(sg);
        for gg=1:ndat.numberEnergyGroups
            sgg = cnodes + (gg-1)*ndg;
            out(sg) = out(sg) - ndat.Diffusion.ScatteringXS(matID,gg,g,1)*M*x(sgg);
        end
    end
end
% Loop through boundary faces
for f=1:mesh.TotalBoundaryFaces
    face = mesh.BoundaryFaces(f);
    fnodes = DoF.FaceCellNodes{face,1};
    flag = mesh.FaceID(face);
    M = FE.FaceMassMatrix{face,1};
    % Loop through energy groups
    gflag = ndat.Diffusion.BCFlags(flag);
    for g=1:ng
        fgnodes = fnodes + (g-1)*ndg;
        switch(gflag)
            case(glob.Dirichlet)
                out(fgnodes) = v(fgnodes);
            case(glob.Robin)
                out(fgnodes) = out(fgnodes) + 0.5*M*v(fgnodes);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

