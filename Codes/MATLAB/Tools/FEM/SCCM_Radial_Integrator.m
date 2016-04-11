%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Polygon-to-disk Mapper
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2016
%
%   Description:    
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Inputs:         1) vertices of the global polygon
%                   2) number of gauss integration points
%
%   Outputs:        1) 
%                   2) 
%                   3) 
%                   4) 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [qx,qw,J,detJ] = SCCM_Radial_Integrator(verts, ng)

t=0; nr = ng;
%compute polygon centroid
temp = polygeom(verts(:,1),verts(:,2));
a = complex(temp(2),temp(3)) ; %centroid of the polygon in complex form
%this will be used as a conformal center for mapping
%construct polygon from the coordinates
Poly = polygon(verts(:,1),verts(:,2));
%map the polygon on to a disk using SCCM
opt.TraceSolution = 'off'; opt.Tolerance = 1e-10;
opt.SolverMethod = 'trust'; opt.InitialGuess = [];
map = diskmap(Poly,opt);
%map the polygon on to a unit disk using SCCM with user centroid of the
%polygon as conformal center
map = center(map,a);
qx = []; qw = []; J = []; detJ = [];
for ir = 1:nr %loop over number of radial locations
    delr = 1/nr;
    r = (ir^2-ir+(1/3))*delr/(ir-0.5); %compute centroid of each sector
    ntheta = round(2*pi/(pi/(ng*(ir-0.5)*delr^2))); %compute number of circles
    deltheta = 2*pi/ntheta; %angle of each sector
    for ith = 1:ntheta %loop over number of circles
        t=t+1;
        theta(ith) = deltheta*(ith-0.5);
        x(t) = r*cos(theta(ith));
        y(t) = r*sin(theta(ith));
        we(t) = (ir-0.5)*deltheta*delr*delr; %compute the weight
        temp = complex(x(t),y(t)); %change it to complex form
        df = evaldiff(map,temp); %compute the derivative of the conformal map
        a = real(df); b = imag(df); %separate real and imaginary part
        jac = [a b;-b a]; %construct the jacobian
        temp2 = eval(map,temp); %use inverse map to locate the points on polygon
        qx = [qx;real(temp2),imag(temp2)];
        qw = [qw;we(t)];
        J(:,:,t) = jac;
        detJ = [detJ;det(jac)];
    end %end loop over number of circles
end %end loop over number of radial locations