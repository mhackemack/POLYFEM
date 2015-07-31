%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Piecewise Linear Discontinuous (PWLD) Gauss
%                     Quadrature  Generator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the mass, stiffness, and
%                   gradient matrices for an element's volume using 
%                   the PWPD DGFEM basis functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = PWLD_quad_gen(varargin)
if nargin == 0
    error('--- No inputs specified. ---')
else
    % Collect Input Arguments
    % -----------------------
    nverts = varargin{5};
    verts = varargin{1}(1:nverts,:);
    q_ord = varargin{2};
    faces = varargin{3};
    % Perform Quick Error Checks
    % --------------------------
    if nargout < 2
        error('Too few output arguments - need at least 2.')
    elseif nargout < 3
        mass_out = 0;
        grad_out = 0;
    elseif nargout == 3
        mass_out = 1;
        grad_out = 0;
    elseif nargout == 4
        mass_out = 1;
        grad_out = 1;
    end
    % Prepare Vertices and Dimensional Space
    % --------------------------------------
    [mv,nv] = size(verts); 
    if nv > mv, verts = verts'; end
    [nv,dim] = size(verts);
    rcenter = mean(verts);
    % 1D Generation
    % -------------
    if dim == 1, error('Choosing not to do PWLD in 1D -- is this just LD???'), end
    % 2D Generation
    % -------------
    if dim == 2
        [qx, qw] = Quad_On_Triangle(q_ord); nqx = length(qw);
        nfaces = 1; ffaces{1} = 1:nv;
        total_sides = nv;
    end
    % 3D Generation
    % -------------
    if dim == 3
        [qx, qw] = Quad_On_Tetra(q_ord); nqx = length(qw);
        nfaces = length(faces);
        ffaces = faces;
        c = 0;
        for i=1:nfaces
            c = c + length(faces{i});
        end
        total_sides = c;
    end
    % Get Reference Node Values
    % -------------------------
    db_ref = get_basis_grads(dim);
    if dim == 2
        b_ref = [ones(nqx,1)-qx(:,1)-qx(:,2),qx(:,1),qx(:,2)];
    else
        b_ref = [ones(nqx,1)-qx(:,1)-qx(:,2)-qx(:,3),qx(:,1),qx(:,2),qx(:,3)];
    end
    % Allocate Memory
    % ---------------
    qw = qw / (ones(1,nqx)*qw);
    ntot = nqx*total_sides;
    side_vol = zeros(total_sides,1);
    quad_nodes_out = zeros(ntot, dim);
    quad_wts_out = zeros(ntot, 1);
    if mass_out
        basis_vals_out = zeros(ntot, nv);
    end
    if grad_out
        basis_grads_out = zeros(ntot, nv, dim);
    end
    % Loop through faces
    % ------------------
    tot_vol = 0;
    ccc = 0;
    for f=1:nfaces
        ff = ffaces{f};
        ne = length(ff);
        bb = 1/ne;
        fcenter = mean(verts(ff,:));
        % Loop through Edges on Face
        % --------------------------
        for e=1:ne
            ccc = ccc + 1;
            if e == ne
                ee = ff(1);
            else
                ee = ff(e+1);
            end
            eee = [ff(e),ee];
            % Edge Triangle/Tetrahedron Information
            % -------------------------------------
            if dim==2
                lverts = [verts(eee,:);rcenter];
            else
                lverts = [verts(eee,:);fcenter;rcenter];
            end
            J = zeros(dim,dim);
            llverts = lverts';
            for d=1:dim
                J(:,d) = llverts(:,d+1) - llverts(:,1);
            end
            % Get other jacobian information. This explicit generation
            % is done for speed for large problems. The elementary matrix
            % generation is taking longer than mesh/DoF generation and system
            % solve times combined...
            if dim==2
                detJ = J(1,1)*J(2,2)-J(2,1)*J(1,2);
%                 invJ = [J(2,2),-J(1,2);-J(2,1),J(1,1)]/detJ;
                svol = detJ/2;
            else
                detJ = J(1,1)*(J(3,3)*J(2,2)-J(3,2)*J(2,3)) - J(2,1)*(J(3,3)*J(1,2)-J(3,2)*J(1,3)) + J(3,1)*(J(2,3)*J(1,2)-J(2,2)*J(1,3));
%                 invJ = [(J(3,3)*J(2,2)-J(3,2)*J(2,3)),-(J(3,3)*J(1,2)-J(3,2)*J(1,3)), (J(2,3)*J(1,2)-J(2,2)*J(1,3));...
%                        -(J(3,3)*J(2,1)-J(3,1)*J(2,3)), (J(3,3)*J(1,1)-J(3,1)*J(1,3)),-(J(2,3)*J(2,1)-J(2,1)*J(1,3));...
%                         (J(3,2)*J(2,1)-J(3,1)*J(2,2)),-(J(3,2)*J(1,1)-J(3,1)*J(1,2)), (J(2,2)*J(1,1)-J(2,1)*J(1,2))]/detJ;
                svol = detJ/6;
            end
            tot_vol = tot_vol + svol;
            side_vol(ccc) = svol;
            % Loop through quad nodes within side volume
            % ------------------------------------------
            ii = (ccc-1)*nqx+1:ccc*nqx;
            for q=1:nqx
                quad_nodes_out(ii(q),:) = lverts(1,:) + (J*qx(q,:)')';
            end
            if mass_out
                tbr = b_ref;
                for q=1:nqx
                    basis_vals_out(ii(q),eee) = tbr(q,1:2);
                    basis_vals_out(ii(q),:) = basis_vals_out(ii(q),:) + tbr(q,end)/nv;
                    if dim == 3
                        basis_vals_out(ii(q),ff) = basis_vals_out(ii(q),ff) + tbr(q,end-1)*bb;
                    end
                end
            end
            if grad_out
                for q=1:nqx
                    
                end
            end
        end
    end
end

% Cleanup Quadrature Weights based on cell measure
% ------------------------------------------------
for i=1:total_sides
    quad_wts_out((i-1)*nqx+1:i*nqx) = qw*side_vol(i);
end
% Set Outputs
varargout{1} = quad_nodes_out;
varargout{2} = quad_wts_out;
if mass_out
    for i=1:total_sides
        basis_vals_out(i,:) = basis_vals_out(i,:) / sum(basis_vals_out(i,:));
    end
    varargout{3} = basis_vals_out;
end
if grad_out
    if mass_out
        varargout{4} = basis_grads_out;
    else
        varargout{3} = basis_grads_out;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                             Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_basis_grads(dim)
if dim == 2
    out = [    -1    -1
                1     0
                0     1];
else
    out = [    -1    -1    -1
                1     0     0
                0     1     0
                0     0     1];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
