%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Draw 3D Angular Quadrature
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB function to plot scalar solutions in 2D.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function draw_quadrature(data)
% Get Angular Quadrature
xa = data.Neutronics.Transport.AngularDirections;
wa = data.Neutronics.Transport.AngularWeights;
[x,w] = get_1st_octant(xa, wa); na = length(w);
w = w / max(w);
% Scatter-plot angles by weight
C = zeros(na,3);
scatter3(x(:,1),x(:,2),x(:,3),100*w,C,'filled');
% Draw 1st octant
n=51;
[az,phi] = meshgrid(linspace(0,pi/2,n));
[xx,yy,zz] = sph2cart(az,phi,1);
C=zeros(n,n,3); C(:,:,3)=1;
hold on
surf(xx,yy,zz,C,'LineStyle','none');
alpha(.5)
view([135,28])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Auxially Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,w] = get_1st_octant(xa, wa)
na = length(wa);
if size(xa,2) ~= 3, error('Give me a 3D angle set...'); end
x = []; w = [];
for m=1:na
    if xa(m,1) > 0 && xa(m,2) > 0 && xa(m,3) > 0
        x = [x;xa(m,:)];
        w = [w,abs(wa(m))];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%