%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Piecewise Linear Discontinuous (PWLD) Basis Function 
%                   Generator - Surface Integrals
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the mass, stiffness, and
%                   gradient matrices for an element's surface using 
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
function varargout = PWLD_surface(varargin)
if nargin == 0
    error('--- No inputs specified. ---')
end
% Collect Input Arguments
% -----------------------
nverts = varargin{5};
verts = varargin{1}(1:nverts,:);
sverts = varargin{2}; nsverts = length(sverts);
flags = varargin{3};
% Prepare Vertices and Dimensional Space
% --------------------------------------
[nv,dim] = size(verts);
rcenter = mean(verts);
scenter = zeros(nv, dim);
ns = zeros(nsverts, 1);
for f=1:nsverts
    ns(f) = length(sverts{f});
    scenter(f,:) = mean(verts(sverts{f},:));
end
higher_types = false;
if length(flags) > 2
    if nargout > 2, higher_types = true; end
end
% Small error checks
if dim == 1, error('Choosing not to do PWLD in 1D -- is this just LD???'); end
if sum(flags) ~= nargout, error('Insufficient output variables.'); end
% Allocate Matrix Memory
% ----------------------
J = zeros(dim,dim);
M = cell(nsverts, 1);
if flags(2)
    db = get_basis_grads(dim);
    znv = zeros(nv, nv);
    G = cell(nsverts, 1);
end
if higher_types
    db = get_basis_grads(dim);
    dbt = db';
    G3 = cell(nsverts, 1);
    G4 = cell(nsverts, 1);
end
for f=1:nsverts
    M{f} = zeros(ns(f), ns(f));
    if flags(2)
        G{f} = cell(dim, 1);
        for d=1:dim
            G{f}{d} = znv;
        end
    end
    if higher_types
        G3{f} = znv;
        G4{f} = znv;
    end
end
% Dimension == 2
if dim == 2
    for f=1:nsverts
        fv = verts(sverts{f},:);
        dv = diff(fv);
        len = norm(dv);
        M{f} = len/6*[2,1;1,2];
        % Build Jacobidan if necessary
        if flags(2) || higher_types
            tverts = [fv;rcenter]';
            for d=1:dim
                J(:,d) = tverts(:,d+1) - tverts(:,1);
            end
            detJ = J(1,1)*J(2,2)-J(2,1)*J(1,2);
            invJ = [J(2,2),-J(1,2);-J(2,1),J(1,1)]/detJ;
        end
        if flags(2)
            b = [1/2,1/2,0]*len;
            c = db*invJ;
            for d=1:dim
                g = c(:,d)*b;
                G{f}{d} = matrix_contribution(dim,nv,g,sverts{f},sverts{f});
            end
        end
        if higher_types
            if flags(3)
                n = [dv(2); -dv(1)];
                Jv = invJ*n*len;
                c = db*(Jv*Jv')*dbt;
                G3{f} = matrix_contribution(dim,nv,c,sverts{f},sverts{f});
            end
            if flags(4)
                JJ = invJ*invJ';
                c = db*JJ*dbt*len*len;
                G4{f} = matrix_contribution(dim,nv,c,sverts{f},sverts{f});
            end
        end
    end
end
% Dimension == 3
if dim == 3
    aa = 1/nv;
    for f=1:nsverts
        ssverts = sverts{f};
%         fv = verts(ssverts,:);
        mm = [2,1,1,0;1,2,1,0;1,1,2,0;0,0,0,0]./12;
        bb = 1/ns(f);
        for i=1:ns(f)
            if i==ns(f)
                ii = [ssverts(end),ssverts(1)];
                iii = [ns(f),1];
            else
                ii = ssverts(i:i+1);
                iii = [i,i+1];
            end
            fverts = [verts(ssverts(iii),:);scenter(f,:)];
            area = heron_3D(fverts);
%             area = 0.5*norm(cross(fverts(2,:) - fverts(1,:), fverts(3,:) - fverts(1,:)));
            if flags(1)
                tm = area*mm;
                MM = M{f};
                MM(iii,iii) =  MM(iii,iii) + tm(1:2,1:2);
                MM          =  MM + bb*bb*tm(3,3);
                MM(iii(1),:) = MM(iii(1),:) + bb*tm(1,3);
                MM(iii(2),:) = MM(iii(2),:) + bb*tm(2,3);
                MM(:,iii(1)) = MM(:,iii(1)) + bb*tm(3,1);
                MM(:,iii(2)) = MM(:,iii(2)) + bb*tm(3,2);
                M{f} = MM;
            end
            if flags(2) || higher_types
                tverts = [verts(ii,:);scenter(f,:);rcenter]';
                for d=1:dim
                    J(:,d) = tverts(:,d+1) - tverts(:,1);
                end
                detJ = J(1,1)*(J(3,3)*J(2,2)-J(3,2)*J(2,3)) - J(2,1)*(J(3,3)*J(1,2)-J(3,2)*J(1,3)) + J(3,1)*(J(2,3)*J(1,2)-J(2,2)*J(1,3));
                invJ = [(J(3,3)*J(2,2)-J(3,2)*J(2,3)),-(J(3,3)*J(1,2)-J(3,2)*J(1,3)), (J(2,3)*J(1,2)-J(2,2)*J(1,3));...
                       -(J(3,3)*J(2,1)-J(3,1)*J(2,3)), (J(3,3)*J(1,1)-J(3,1)*J(1,3)),-(J(2,3)*J(1,1)-J(2,1)*J(1,3));...
                        (J(3,2)*J(2,1)-J(3,1)*J(2,2)),-(J(3,2)*J(1,1)-J(3,1)*J(1,2)), (J(2,2)*J(1,1)-J(2,1)*J(1,2))]/detJ;
            end
            if flags(2)
                b = [1/6,1/6,1/6,0]*2*area;
                c = db*invJ;
                for d=1:dim
                    g = c(:,d)*b;
                    GG = G{f}{d};
                    GG(ii,ii) = GG(ii,ii) + g(1:2,1:2);
                    GG = GG + aa*aa*g(4,4);
                    GG(ssverts,ssverts) = GG(ssverts,ssverts) + bb*bb*g(3,3);
                    GG(ii(1),ssverts)  = GG(ii(1),ssverts)  + bb*g(1,3);
                    GG(ii(2),ssverts)  = GG(ii(2),ssverts)  + bb*g(2,3);
                    GG(ssverts,ii(1))  = GG(ssverts,ii(1))  + bb*g(3,1);
                    GG(ssverts,ii(2))  = GG(ssverts,ii(2))  + bb*g(3,2);
                    GG(:,ssverts) = GG(:,ssverts) + aa*bb*g(4,3);
                    GG(ssverts,:) = GG(ssverts,:) + aa*bb*g(3,4);
                    G{f}{d} = GG;
                end
            end
            if higher_types
                if flags(3)
                    v1 = verts(ii(1),:) - scenter(f,:);
                    v2 = verts(ii(2),:) - scenter(f,:);
                    n = cross(v1,v2);
                    Jv = invJ*n*2*area;
                    g = db*(Jv*Jv')*dbt;
                    GG = G3{f};
                    GG(ii,ii) = GG(ii,ii) + g(1:2,1:2);
                    GG = GG + aa*aa*g(4,4);
                    GG(ssverts,ssverts) = GG(ssverts,ssverts) + bb*bb*g(3,3);
                    GG(ii(1),ssverts)  = GG(ii(1),ssverts)  + bb*g(1,3);
                    GG(ii(2),ssverts)  = GG(ii(2),ssverts)  + bb*g(2,3);
                    GG(ssverts,ii(1))  = GG(ssverts,ii(1))  + bb*g(3,1);
                    GG(ssverts,ii(2))  = GG(ssverts,ii(2))  + bb*g(3,2);
                    GG(:,ssverts) = GG(:,ssverts) + aa*bb*g(4,3);
                    GG(ssverts,:) = GG(ssverts,:) + aa*bb*g(3,4);
                    G3{f} = G4{f} + GG;
                end
                if flags(4)
                    JJ = invJ*invJ';
                    g = db*JJ*dbt*area;
                    GG = G4{f};
                    GG(ii,ii) = GG(ii,ii) + g(1:2,1:2);
                    GG = GG + aa*aa*g(4,4);
                    GG(ssverts,ssverts) = GG(ssverts,ssverts) + bb*bb*g(3,3);
                    GG(ii(1),ssverts)  = GG(ii(1),ssverts)  + bb*g(1,3);
                    GG(ii(2),ssverts)  = GG(ii(2),ssverts)  + bb*g(2,3);
                    GG(ssverts,ii(1))  = GG(ssverts,ii(1))  + bb*g(3,1);
                    GG(ssverts,ii(2))  = GG(ssverts,ii(2))  + bb*g(3,2);
                    GG(:,ssverts) = GG(:,ssverts) + aa*bb*g(4,3);
                    GG(ssverts,:) = GG(ssverts,:) + aa*bb*g(3,4);
                    G4{f} = G4{f} + GG;
                end
            end
        end
    end
end
% Set Matrix Output Arguments
% ---------------------------
counter = 1;
if flags(1) == 1
    varargout{counter} = M;
    counter = counter + 1;
end
if flags(2) == 1
    varargout{counter} = G;
    counter = counter + 1;
end
if higher_types
    if flags(3) == 1
        varargout{counter} = G3;
        counter = counter + 1;
    end
    if flags(4) == 1
        varargout{counter} = G4;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                             Function List
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = matrix_contribution(dim,nv,mat,v,fv)
a = 1/nv;
out = zeros(nv,nv);
out(v,v) = mat(1:length(v),1:length(v));
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = heron_3D( fv )
a = norm(fv(2,:) - fv(1,:));
b = norm(fv(3,:) - fv(2,:));
c = norm(fv(1,:) - fv(3,:));
s = (a+b+c)/2;
out = sqrt(s*(s-a)*(s-b)*(s-c));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%