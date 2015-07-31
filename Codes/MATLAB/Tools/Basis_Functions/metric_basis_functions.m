%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Metric Basis Function Generator
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%
%   Description:    MATLAB script to produce the basis function and gradient
%                   values using Metric Coordinates.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = metric_basis_functions( varargin )
nout = nargout;
grad_bool = false;
if nout > 1, grad_bool = true; end
% Collect Input Arguments
% -----------------------
verts = varargin{1};
qx = varargin{2};
faces = varargin{3};
% Prepare Vertices and Dimensional Space
% --------------------------------------
[nverts, dim] = size(verts);
% Quick Error Checking
% --------------------

% Build Basis Function Sets
% -------------------------
[bout, gout] = getBFs(dim, nverts, verts, faces, qx, grad_bool);
% Set Output Arguments
% --------------------
varargout{1} = bout;
if grad_bool, varargout{2} = gout; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function List
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rv, sv] = get_vertex_differences( verts, qx )
[nv, dim] = size(verts);
nqx = size(qx, 1);
rv = zeros(nv, dim, nqx);
sv = zeros(nv, nqx);
zn = ones(nv, 1);
zones = ones(dim,1);
for q=1:nqx
    rv(:,:,q) = zn*qx(q,:) - verts;
    rrv = rv(:,:,q).*rv(:,:,q);
    sv(:,q) = sqrt(rrv*zones);
end
sv = sv';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ef, efder] = get_face_functions(verts, nverts, qx, dim, faces)
nqx = size(qx, 1); nf = length(faces);
rv = verts(1:nverts,:);
ef = zeros(nqx, nf);
efder = zeros(nf, dim, nqx);
if dim == 2
    for q=1:nqx
        tefder = efder(:,:,q);
        drv = [rv(:,1) - qx(q,1), rv(:,2) - qx(q,2)];
        for i=1:nf
            ii = [i, mod(i,nverts)+1];
            norm1 = norm(drv(ii(1),:));
            norm2 = norm(drv(ii(2),:));
            ef(q,i) =  norm1 + norm2 - norm(diff(verts(ii,:)));
            tefder(i,:) = -drv(ii(1),:)/norm1 - drv(ii(2),:)/norm2;
        end
        efder(:,:,q) = tefder;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [adjf,nadjf] = get_face_dependencies(nverts, faces)
adjf = cell(nverts, 1);
nadjf = cell(nverts, 1);
nf = length(faces);
% Loop through faces and determine adjacencies for vertices
for f=1:nf
    fv = faces{f};
    for i=1:length(fv)
        adjf{fv(i)} = [adjf{fv(i)}, f];
    end
end
% Loop through vertices and determine missing faces
for i=1:nverts
    totf = 1:nf; totf(adjf{i}) = [];
    nadjf{i} = totf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = get_wraparound_verts(dim, nverts)
if dim == 2
    out = zeros(nverts, 2);
    for i=1:nverts
        if i==1
            out(i,:) = [nverts,2];
        else
            out(i,:) = [i-1,mod(i,nverts)+1];
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bout, gout] = getBFs(dim, nverts, verts, faces, qx, grad_bool)
nqx = size(qx,1);
gout = zeros(nverts, dim, nqx);
bout = get_single_BF(dim, nverts, verts, faces, qx);
if grad_bool
    b = sqrt(eps);
    tqx = b*(1+qx(:,1));
    tqy = b*(1+qx(:,2));
    % x-dimension
    tq = [qx(:,1)+tqx,qx(:,2)];
    tbx = get_single_BF(dim, nverts, verts, faces, tq);
    % y-dimension
    tq = [qx(:,1),qx(:,2)+tqy];
    tby = get_single_BF(dim, nverts, verts, faces, tq);
    % Loop through quadrature
    for q=1:nqx
        gout(:,1,q) = (tbx(q,:) - bout(q,:))'/tqx(q);
        gout(:,2,q) = (tby(q,:) - bout(q,:))'/tqx(q);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bout = get_single_BF(dim, nverts, verts, faces, qx)
% Get Problem Preliminaries
% -------------------------
[adjf, nadjf] = get_face_dependencies(nverts, faces);
wa_verts = get_wraparound_verts(dim, nverts);
[fbool, w_face] = are_points_on_face(verts, faces, qx);
[rv, sv]      = get_vertex_differences(verts, qx);
[ef, efder]   = get_face_functions(verts, nverts, qx, dim, faces);
nqx = size(qx,1);
% Allocate Memory Space
% ---------------------
bout = zeros(nqx, nverts);
% gout = zeros(nverts, dim, nqx);
k = ones(nqx,nverts); s = ones(nqx,nverts); %sgrad = zeros(nverts,dim,nqx);
% Build Basis Function Sets and Gradients
% ---------------------------------------
for q=1:nqx
    % Loop through vertices and build temporary variables
    for i=1:nverts
        nadf = nadjf{i}; nnadf = length(nadf);
        r = zeros(nnadf);
        % Loop through non-adjacent faces
        for j=1:nnadf
            r(j) = ef(q,nadf(j));
            s(q,i) = s(q,i)*ef(q,nadf(j));
        end
%         % Build sgrad terms
%         for j=1:nnadf
%             rsum = 0;
%             for m=1:nnadf
%                 if j==m, continue; end
%                 rsum = rsum + r(m);
%             end
%             sgrad(i,:,q) = sgrad(i,:,q) + rsum*efder(nadf(j),:,q);
%         end
    end
    
%     % Loop through vertices and build temporary variables
%     for i=1:nverts
%         % get adjacent nodes
%         ii = [i-1,i,mod(i,nverts)+1];
%         if ii(1)==0, ii(1) = nverts; end
%         vv = verts(ii,:);
%         % calculate k value
%         AA = 0.5*abs(det([vv(:,1),vv(:,2),[1;1;1]]));
%         if AA > 1e-13, k(q,i) = AA; end
%         % calculate s value
%         nadf = nadjf{i};
%         % Loop through non-adjacent faces
%         for j=1:length(nadf)
%             fv = faces{nadf(j)};
%             fvm = wa_verts(fv(1),1);
%             fvp = wa_verts(fv(2),2);
%             l2  = sv(q,fv(1));
%             l3  = sv(q,fv(2));
%             l13 = norm(diff(verts([fvm,fv(2)],:)));
%             l24 = norm(diff(verts([fvp,fv(1)],:)));
%             l12 = norm(diff(verts([fvm,fv(1)],:)));
%             l34 = norm(diff(verts([fvp,fv(2)],:)));
%             l23 = norm(diff(verts(fv,:)));
%             c1  = (l12*l12 + l23*l23 - l13*l13)/(2*l12*l23);
%             c2  = (l34*l34 + l23*l23 - l24*l24)/(2*l34*l23);
%             s1  = abs(det([verts([fvm,fv],1),verts([fvm,fv],2),[1;1;1]]))/(l13*l23);
%             s2  = abs(det([verts([fv,fvp],1),verts([fv,fvp],2),[1;1;1]]))/(l23*l24);
% %             s1  = abs(det([verts([fvm,fv],1),verts([fvm,fv],2),[1;1;1]]))/(l12*l23);
% %             s2  = abs(det([verts([fv,fvp],1),verts([fv,fvp],2),[1;1;1]]))/(l34*l23);
%             if s1 < 1e-13, s1=1;end;
%             if s2 < 1e-13, s2=1;end
%             A = l3*l3 + l23*l23 - 2*l23*l3*c2 - l2*l2;
%             B = l2*l2 + l23*l23 - 2*l23*l2*c1 - l3*l3;
%             Ca = l23/2*(A*s1 + B*s2);
%             Cb = A + B;
%             Cc = l23*(A*c1 + B*c2);
%             % apply segement function
%             num   = Ca * ((l2 + l3)^2 - l23*l23);
%             denom = (l2 + l3)*Cb - Cc;
%             t = num/denom;
%             s(q,i) = s(q,i)*t;
%         end
%     end
    
    % Construct basis functions
    sumks = sum(k(q,:).*s(q,:));
    bout(q,:) = k(q,:).*s(q,:) / sumks;
%     % Construct Gradient Terms
%     R = zeros(nverts,dim);
%     for i=1:nverts
%         R(i,:) = sgrad(i,:,q) / (k(q,i)*s(q,i));
%     end
%     for i=1:nverts
%         gout(i,:,q) = bout(q,i) * (R(i,:) - bout(q,:)*R);
%     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

