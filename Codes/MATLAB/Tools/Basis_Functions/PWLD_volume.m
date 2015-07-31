%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Piecewise Linear Discontinuous (PWLD) Basis Function 
%                   Generator - Volume Integrals
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the mass, stiffness, and
%                   gradient matrices for an element's volume using 
%                   the PWLD DGFEM basis functions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes - 2D:     1) 'faces' input not needed (ignored)
%                   2) 'verts' input is in the form (npts x ndim)
%                   3) 'verts' coordinates need to be CCW
%
%   Notes - 3D:     1) 'verts' and 'faces' input both needed
%                   2) 'faces' holds the vertex numberings of 'verts'
%                   3) 'faces' can take either cell or array structure 
%                      - if array structure: form (nfaces x npts_face)
%                                            where npts_face is constant
%                      - array structure can only be used if each face
%                        has the same number of vertices
%                   4) Vertices in 'verts' do not need any proper ordering
%                   5) Vertices on each face need to be in CCW order
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout] = PWLD_volume(varargin)

if nargin == 0
    error('--- No inputs specified. ---')
else
    % Collect Input Arguments
    % -----------------------
    nverts = varargin{5};
    verts = varargin{1}(1:nverts,:);
    faces = varargin{2};
    flags = varargin{3};
    % Prepare Vertices and Dimensional Space
    % --------------------------------------
    [mv,nv] = size(verts); 
    if nv > mv, verts = verts'; end
    [nv,dim] = size(verts);
    rcenter = mean(verts);
    % Allocate Matrix Memory
    % ----------------------
    a = 1/nv;
    M = zeros(nv,nv);       % mass matrix
    K = zeros(nv,nv);       % stiffness matrix
    G = cell(dim,1);        % gradient matrix
    for i=1:dim
        G{i} = zeros(nv,nv);
    end
    % 1D Generation
    % -------------
    if dim == 1, error('Choosing not to do PWLD in 1D -- is this just LD???'), end
    % 2D Generation
    % -------------
    if dim == 2
        ffaces{1} = 1:nv;
%         if nargin ~= 3
%             flags = [1,1,1];
%         end
    end
    % 3D Generation
    % -------------
    if dim == 3
        % Check face structure
        % --------------------
        if nargin < 2
%             error('--- No face input specified. ---')
            faces = convhulln(verts);
        else
            if isempty(faces)
%                 error('--- No face input specified. ---')
                faces = convhulln(verts);
            end
        end
        if iscell(faces)
            ffaces = faces;
        else
            ffaces = cell(length(faces),1);
            for i=1:size(faces,1)
                ffaces{i} = faces(i,:);
            end
        end
    end
    
    % Small error checks
    if sum(flags) ~= nargout
        error('Insufficient output variables.')
    end
    % get ref vals
    if dim==2
        m = [2,1,1;1,2,1;1,1,2]./12;
    elseif dim==3
        m = [2,1,1,1;1,2,1,1;1,1,2,1;1,1,1,2]./20;
    end
    db = get_basis_grads(dim);
    dbt = db';
    if flags(3) == 1
        if dim == 2
            b = ones(1,dim+1)/6;
        else
            b = ones(1,dim+1)/24;
        end
    end
    J = zeros(dim,dim);
    
    % Loop through Faces
    % ------------------
    for f=1:length(ffaces)
        ff = ffaces{f};
        ne = length(ff);
        bb = 1/ne;
        fcenter = mean(verts(ff,:));
        % Loop through Edges on Face
        % --------------------------
        for e=1:ne
            if e == ne
                ee = ff(1);
            else
                ee = ff(e+1);
            end
            eee = [ff(e),ee];
            % Edge Triangle/Tetrahedron Information
            % -------------------------------------
            if dim==2
                lverts = [verts(eee,:);rcenter]';
            else
                lverts = [verts(eee,:);fcenter;rcenter]';
            end
            for d=1:dim
                J(:,d) = lverts(:,d+1) - lverts(:,1);
            end
            % Get other jacobian information. This explicit generation
            % is done for speed for large problems. The elementary matrix
            % generation is taking longer than mesh/DoF generation and system
            % solve times combined...
            if dim==2
                detJ = J(1,1)*J(2,2)-J(2,1)*J(1,2);
                invJ = [J(2,2),-J(1,2);-J(2,1),J(1,1)]/detJ;
                svol = detJ/2;
            else
                detJ = J(1,1)*(J(3,3)*J(2,2)-J(3,2)*J(2,3)) - J(2,1)*(J(3,3)*J(1,2)-J(3,2)*J(1,3)) + J(3,1)*(J(2,3)*J(1,2)-J(2,2)*J(1,3));
                invJ = [(J(3,3)*J(2,2)-J(3,2)*J(2,3)),-(J(3,3)*J(1,2)-J(3,2)*J(1,3)), (J(2,3)*J(1,2)-J(2,2)*J(1,3));...
                       -(J(3,3)*J(2,1)-J(3,1)*J(2,3)), (J(3,3)*J(1,1)-J(3,1)*J(1,3)),-(J(2,3)*J(1,1)-J(2,1)*J(1,3));...
                        (J(3,2)*J(2,1)-J(3,1)*J(2,2)),-(J(3,2)*J(1,1)-J(3,1)*J(1,2)), (J(2,2)*J(1,1)-J(2,1)*J(1,2))]/detJ;
                svol = detJ/6;
            end
            % Mass Matrix Contributions
            % -------------------------
            if flags(1)
                MM = zeros(nv,nv);
                mm = svol*m;
                MM(eee,eee) = mm(1:2,1:2);
                MM = MM + a*a*mm(end,end);
                MM(eee(1),:) = MM(eee(1),:) + a*mm(1,end);
                MM(eee(2),:) = MM(eee(2),:) + a*mm(2,end);
                MM(:,eee(1)) = MM(:,eee(1)) + a*mm(end,1);
                MM(:,eee(2)) = MM(:,eee(2)) + a*mm(end,2);
                if dim == 3
                    MM(ff,ff) = MM(ff,ff) + bb*bb*mm(end-1,end-1);
                    MM(ff,:) = MM(ff,:) + a*bb*mm(end-1,end);
                    MM(:,ff) = MM(:,ff) + a*bb*mm(end,end-1);
                    MM(eee(1),ff) = MM(eee(1),ff) + bb*mm(1,end-1);
                    MM(eee(2),ff) = MM(eee(2),ff) + bb*mm(2,end-1);
                    MM(ff,eee(1)) = MM(ff,eee(1)) + bb*mm(end-1,1);
                    MM(ff,eee(2)) = MM(ff,eee(2)) + bb*mm(end-1,2);
                end
                M = M + MM;
            end
            % Stiffness Matrix Contributions
            % ------------------------------
            if flags(2)
                s = svol*(db*(invJ*invJ')*dbt);
                KK = zeros(nv,nv);
                KK(eee,eee) = s(1:2,1:2);
                KK = KK + a*a*s(end,end);
                KK(eee(1),:) = KK(eee(1),:) + a*s(1,end);
                KK(eee(2),:) = KK(eee(2),:) + a*s(2,end);
                KK(:,eee(1)) = KK(:,eee(1)) + a*s(end,1);
                KK(:,eee(2)) = KK(:,eee(2)) + a*s(end,2);
                if dim == 3
                    KK(ff,ff) = KK(ff,ff) + bb*bb*s(end-1,end-1);
                    KK(ff,:) = KK(ff,:) + a*bb*s(end-1,end);
                    KK(:,ff) = KK(:,ff) + a*bb*s(end,end-1);
                    KK(eee(1),ff) = KK(eee(1),ff) + bb*s(1,end-1);
                    KK(eee(2),ff) = KK(eee(2),ff) + bb*s(2,end-1);
                    KK(ff,eee(1)) = KK(ff,eee(1)) + bb*s(end-1,1);
                    KK(ff,eee(2)) = KK(ff,eee(2)) + bb*s(end-1,2);
                end
                K = K + KK;
            end
            % Gradient Matrix Contributions
            % -----------------------------
            if flags(3)
                c = db*invJ*detJ;
                for j=1:dim
                    g = (c(:,j)*b)';
                    GG = zeros(nv,nv);
                    GG(eee,eee) = g(1:2,1:2);
                    GG = GG + a*a*g(end,end);
                    GG(eee(1),:) = GG(eee(1),:) + a*g(1,end);
                    GG(eee(2),:) = GG(eee(2),:) + a*g(2,end);
                    GG(:,eee(1)) = GG(:,eee(1)) + a*g(end,1);
                    GG(:,eee(2)) = GG(:,eee(2)) + a*g(end,2);
                    if dim == 3
                        GG(ff,ff) = GG(ff,ff) + bb*bb*g(end-1,end-1);
                        GG(ff,:) = GG(ff,:) + a*bb*g(end-1,end);
                        GG(:,ff) = GG(:,ff) + a*bb*g(end,end-1);
                        GG(eee(1),ff) = GG(eee(1),ff) + bb*g(1,end-1);
                        GG(eee(2),ff) = GG(eee(2),ff) + bb*g(2,end-1);
                        GG(ff,eee(1)) = GG(ff,eee(1)) + bb*g(end-1,1);
                        GG(ff,eee(2)) = GG(ff,eee(2)) + bb*g(end-1,2);
                    end
                    G{j} = G{j} + GG;
                end
            end
        end
    end
    % Set Outputs
    counter = 1;
    if flags(1) == 1
        varargout{counter} = M;
        counter = counter + 1;
    end
    if flags(2) == 1
        varargout{counter} = K;
        counter = counter + 1;
    end
    if flags(3) == 1
        varargout{counter} = G;
    end
end


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                             Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J, invJ, detJ, vol] = get_jacobian(verts)
dim = size(verts,2);
J = zeros(dim,dim);
for i=1:dim
    J(:,i) = verts(i+1,:)' - verts(1,:)';
end
invJ = inv(J);
detJ = abs(det(J));
if dim == 2
    vol = detJ/2;
else
    vol = detJ/6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_ref_mass_matrix(dim)
if dim==2
    out = [2,1,1;1,2,1;1,1,2]./12;
elseif dim==3
    out = [2,1,1,1;1,2,1,1;1,1,2,1;1,1,1,2]./20;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_local_stiffness_matrix(dim,invJ,detJ)
% R = lens.^2/(4*A);
% out = [  2*R(1), R(3) - R(1) - R(2), R(2) - R(1) - R(3);...
%          R(3) - R(1) - R(2), 2*R(2), R(1) - R(2) - R(3);...
%          R(2) - R(1) - R(3), R(1) - R(2) - R(3), 2*R(3)    ]./2;
db = get_basis_grads(dim);
if dim == 2
    out = (db*(invJ*invJ')*db').*detJ/2;
else
    out = (db*(invJ*invJ')*db').*detJ/6;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_local_gradient_term(dim,invJ,detJ)
out = cell(dim,1);
% a = -([lens,lens].*vecs)./6;
if dim == 2
    b = ones(1,dim+1)./6;
else
    b = ones(1,dim+1)./24;
end
db = get_basis_grads(dim);
c = db*invJ*detJ;
for i=1:dim
    out{i} = (c(:,i)*b)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_basis_grads(dim)
% out = [-ones(1,dim);diag(ones(dim,1))];
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
function [lens,vecs] = get_side_lengths(verts)
[nv,dim] = size(verts);
lens = zeros(nv,1);
vecs = zeros(nv,dim);
if dim == 2
    dd = verts(3,:)-verts(2,:); lens(1) = norm(dd); vecs(1,:) = [dd(2),-dd(1)]./lens(1);
    dd = verts(1,:)-verts(3,:); lens(2) = norm(dd); vecs(2,:) = [dd(2),-dd(1)]./lens(2);
    dd = verts(2,:)-verts(1,:); lens(3) = norm(dd); vecs(3,:) = [dd(2),-dd(1)]./lens(3);
elseif dim == 3
%     u = cross([verts(1,:)-verts(2,:), 0]',[verts(3,:)-verts(2,:), 0]')'; vecs(1,:) = u./norm(u);
%     u = cross([verts(3,:)-verts(2,:), 0]',[verts(1,:)-verts(2,:), 0]')'; vecs(2,:) = u./norm(u);
%     u = cross([verts(1,:)-verts(3,:), 0]',[verts(2,:)-verts(3,:), 0]')'; vecs(3,:) = u./norm(u);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = matrix_contribution(dim,nv,mat,v,fv)
a = 1/nv;
nnv = length(v);
out = zeros(nv,nv);
out(v,v) = mat(1:nnv,1:nnv);
for i=1:length(v)
    out(v(i),:) = out(v(i),:) + a*mat(i,end);
    out(:,v(i)) = out(:,v(i)) + a*mat(end,i);
end
out = out + a*a*mat(end,end);
if dim == 3
    b = 1/length(fv);
    out(fv,:) = out(fv,:) + a*b*mat(end-1,end);
    out(:,fv) = out(:,fv) + a*b*mat(end,end-1);
    out(fv,fv) = out(fv,fv) + b*b*mat(end-1,end-1);
    for i=1:length(v)
        out(v(i),fv) = out(v(i),fv) + b*mat(i,end-1);
        out(fv,v(i)) = out(fv,v(i)) + b*mat(end-1,i);
    end
end
