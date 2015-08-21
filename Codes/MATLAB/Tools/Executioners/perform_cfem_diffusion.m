%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Solve CFEM Diffusion Problem
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = perform_cfem_diffusion(ndat, solvdat, mesh, DoF, FE, x, A)
global glob
ndg = DoF.TotalDoFs;
ndof = ndat.numberEnergyGroups * ndg;
if nargin < 5 || isempty(x)
    x = ones(ndof,1);
else
    x = cell_to_vector(x, DoF);
end
if length(x) > glob.maxSparse
    rhs = get_rhs(x, ndat, mesh, DoF, FE);
    [x,~,~,~] = pcg(@(x) get_Ax(x, ndat, mesh, DoF, FE),rhs,solvdat.relativeTolerance,solvdat.maxIterations,[],[],x);
else
    if nargin < 6 || ~exist('A', 'var')
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
x = vector_to_cell(x,DoF);
% Outputs
nout = nargout;
varargout{1} = ndat;
varargout{2} = x;
if nout == 3
    if exist('A', 'var')
        varargout{3} = A; 
    else
        varargout{3} = []; 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                              Function Listing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L,rhs] = get_global_matrices(v, ndat, mesh, DoF, FE)
global glob
n = length(v);
ng = ndat.numberEnergyGroups;
ndg = DoF.TotalDoFs;
% Allocate Memory
if n > glob.maxMatrix
    [L, rhs] = get_sparse_matrices(v, ndat, mesh, DoF, FE);
    return
else
    L = zeros(n,n);
    rhs = zeros(n,1);
end
rev_str = [];
% Loop through cells
for tcell=1:mesh.TotalCells
    % Print Current Cell Information
%     if glob.print_info
%         msg = sprintf('   -> Building Cell: %d of %d',tcell,mesh.TotalCells);
%         fprintf([rev_str,msg]);
%         rev_str = repmat(sprintf('\b'), 1, length(msg));
%     end
    
    cnodes = DoF.ConnectivityArray{tcell}; ncnodes = length(cnodes);
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    K = FE.CellStiffnessMatrix{tcell};
    % Loop through energy groups
    for g=1:ndat.numberEnergyGroups
        sg = cnodes + (g-1)*ndg;
        L(sg,sg) = (L(sg,sg) + ndat.Diffusion.DiffXS(matID,g)*K + ndat.Diffusion.TotalXS(matID,g)*M);
        if ndat.Diffusion.MMS
            gfunc = ndat.Diffusion.ExtSource{g};
            qx = FE.CellQuadNodes{tcell};
            qw = FE.CellQuadWeights{tcell};
            cb = FE.CellBasisValues{tcell};
            tvec = cb'*(qw.*gfunc(qx));
        else
            tvec = ndat.Diffusion.ExtSource(matID,g)*M*ones(ncnodes,1);
        end
        % Loop through right-hand-side energy groups
        for gg=1:ndat.numberEnergyGroups
            sgg = cnodes + (gg-1)*ndg;
            xsec = ndat.Diffusion.FissSpec(matID,g)/ndat.keff*ndat.Diffusion.FissionXS(matID,gg)*ndat.Transport.NuBar(matID,gg) + ndat.Diffusion.ScatteringXS(matID,gg,g);
            tvec = tvec + xsec*M*v(sgg);
        end
    end
    rhs(sg) = rhs(sg) + tvec;
end
% Loop through boundary faces
for f=1:mesh.TotalBoundaryFaces
    % Print Current Face Information
%     if glob.print_info
%         msg = sprintf('   -> Building boundary face: %d of %d',f,mesh.TotalBoundaryFaces);
%         fprintf([rev_str,msg]);
%         rev_str = repmat(sprintf('\b'), 1, length(msg));
%     end
    
    face = mesh.BoundaryFaces(f);
    fnodes = DoF.FaceCellNodes{face,1};
    flag = mesh.FaceID(face);
    M = FE.FaceMassMatrix{face,1};
    F = FE.FaceFunctionMatrix{face,1};
    % Loop through energy groups
    gflag = ndat.Diffusion.BCFlags(flag);
    for g=1:ng
        fgnodes = fnodes + (g-1)*ndg;
        switch(gflag)
            case(glob.Dirichlet)
                L(fgnodes,:) = 0;
                for i=1:length(fgnodes)
                    L(fgnodes(i),fgnodes(i)) = 1;
                end
                rhs(fgnodes) = ndat.Diffusion.BCVals(flag,g);
            case(glob.Neumann)
                rhs(fgnodes) = rhs(fgnodes) - ndat.Diffusion.BCVals(flag,g)*F;
            case(glob.Robin)
                L(fgnodes,fgnodes) = L(fgnodes,fgnodes) + 0.5*M;
                rhs(fgnodes) = rhs(fgnodes) + 2*ndat.Diffusion.BCVals(flag,g)*F;
            case(glob.Function)
                L(fgnodes,:) = 0;
                for i=1:length(fgnodes)
                    L(fgnodes(i),fgnodes(i)) = 1;
                end
                fxn = DoF.NodeLocations(fnodes,:);
                rhs(fgnodes) = ndat.Diffusion.BCVals{flag,g}(fxn);
        end
    end
end
if ~issparse(L), L = sparse(L); end
% if glob.print_info, fprintf(rev_str); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [L, rhs] = get_sparse_matrices(v, ndat, mesh, DoF, FE)
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
%     if glob.print_info
%         msg = sprintf('   -> Building Cell: %d of %d',tcell,mesh.TotalCells);
%         fprintf([rev_str,msg]);
%         rev_str = repmat(sprintf('\b'), 1, length(msg));
%     end
    
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
        I = [I;rows(:)]; J = [J;cols(:)]; TMAT = [TMAT;tmat(:)];
        if ndat.Diffusion.MMS
            gfunc = ndat.Diffusion.ExtSource{g};
            qx = FE.CellQuadNodes{tcell};
            qw = FE.CellQuadWeights{tcell};
            cb = FE.CellBasisValues{tcell};
            tvec = cb'*(qw.*gfunc(qx));
        else
            tvec = ndat.Diffusion.ExtSource(matID,g)*M*ones(ncnodes,1);
        end
        % Loop through right-hand-side energy groups
        for gg=1:ndat.numberEnergyGroups
            sgg = cnodes + (gg-1)*ndg;
            xsec = ndat.Diffusion.FissSpec(matID,g)/ndat.keff*ndat.Diffusion.FissionXS(matID,gg)*ndat.Transport.NuBar(matID,gg) + ndat.Diffusion.ScatteringXS(matID,gg,g);
            tvec = tvec + xsec*M*v(sgg);
        end
    end
    rhs(sg) = rhs(sg) + tvec;
end
L = sparse(I,J,TMAT,n,n);
dirch_nums = [];
% Loop through boundary faces
for f=1:mesh.TotalBoundaryFaces
    % Print Current Face Information
%     if glob.print_info
%         msg = sprintf('   -> Building boundary face: %d of %d',f,mesh.TotalBoundaryFaces);
%         fprintf([rev_str,msg]);
%         rev_str = repmat(sprintf('\b'), 1, length(msg));
%     end
    
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
                rhs(fgnodes) = ndat.Diffusion.BCVals(flag,g)*ones(length(fgnodes),1);
            case(glob.Neumann)
                rhs(fgnodes) = rhs(fgnodes) - ndat.Diffusion.BCVals(flag,g)*F;
            case(glob.Robin)
                L(fgnodes,fgnodes) = L(fgnodes,fgnodes) + 0.5*M;
                rhs(fgnodes) = rhs(fgnodes) + 2*ndat.Diffusion.BCVals(flag,g)*F;
            case(glob.Function)
                dirch_nums = [dirch_nums,fgnodes];
                fxn = DoF.NodeLocations(fnodes,:);
                rhs(fgnodes) = ndat.Diffusion.BCVals{flag,g}(fxn);
        end
    end
end
if ~isempty(dirch_nums)
    dirch_nums = unique(dirch_nums);
    L(dirch_nums,:) = 0;
    diag_id = (dirch_nums-1)*(n+1)+1;
    L(diag_id) = 1;
end
% if glob.print_info, fprintf(rev_str); end
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
%     if glob.print_info
%         msg = sprintf('   -> Building Cell: %d of %d',tcell,mesh.TotalCells);
%         fprintf([rev_str,msg]);
%         rev_str = repmat(sprintf('\b'), 1, length(msg));
%     end
    
    cnodes = DoF.ConnectivityArray{tcell}; ncnodes = length(cnodes);
    matID = mesh.MatID(tcell);
    M = FE.CellMassMatrix{tcell};
    % Loop through energy groups
    for g=1:ndat.numberEnergyGroups
        sg = cnodes + (g-1)*ndg;
        if ndat.Diffusion.MMS
            gfunc = ndat.Diffusion.ExtSource{g};
            qx = FE.CellQuadNodes{tcell};
            qw = FE.CellQuadWeights{tcell};
            cb = FE.CellBasisValues{tcell};
            tvec = cb'*(qw.*gfunc(qx));
        else
            tvec = ndat.Diffusion.ExtSource(matID,g)*M*ones(ncnodes,1);
        end
        % Loop through right-hand-side energy groups
        for gg=1:ndat.numberEnergyGroups
            sgg = cnodes + (gg-1)*ndg;
            xsec = ndat.Diffusion.FissSpec(matID,g)/ndat.keff*ndat.Diffusion.FissionXS(matID,gg)*ndat.Transport.NuBar(matID,gg) + ndat.Diffusion.ScatteringXS(matID,gg,g);
            tvec = tvec + xsec*M*v(sgg);
        end
    end
    rhs(sg) = rhs(sg) + tvec;
end
% Loop through boundary faces
for f=1:mesh.TotalBoundaryFaces
    % Print Current Face Information
%     if glob.print_info
%         msg = sprintf('   -> Building boundary face: %d of %d',f,mesh.TotalBoundaryFaces);
%         fprintf([rev_str,msg]);
%         rev_str = repmat(sprintf('\b'), 1, length(msg));
%     end
    
    face = mesh.BoundaryFaces(f);
    fnodes = DoF.getFaceCellNodes(face,1);
    flag = mesh.FaceID(face);
    F = FE.FaceFunctionMatrix{face,1};
    % Loop through energy groups
    gflag = ndat.Diffusion.BCFlags(flag);
    for g=1:ng
        fgnodes = fnodes + (g-1)*ndg;
        switch(gflag)
            case(glob.Dirichlet)
                rhs(fgnodes) = ndat.Diffusion.BCVals(flag,g)*ones(length(fgnodes),1);
            case(glob.Neumann)
                rhs(fgnodes) = rhs(fgnodes) - ndat.Diffusion.BCVals(flag,g)*F;
            case(glob.Robin)
                rhs(fgnodes) = rhs(fgnodes) + 2*ndat.Diffusion.BCVals(flag,g)*F;
            case(glob.Function)
                fxn = DoF.NodeLocations(fnodes,:);
                rhs(fgnodes) = ndat.Diffusion.BCVals{flag,g}(fxn);
        end
    end
end
% if glob.print_info, fprintf(rev_str); end
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
        out(sg) = out(sg) + (ndat.Diffusion.DiffXS(matID,g)*K + ndat.Diffusion.TotalXS(matID,g)*M)*v(sg);
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
            case(glob.Function)
                out(fgnodes) = v(fgnodes);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%