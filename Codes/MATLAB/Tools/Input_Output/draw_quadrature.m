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
xa = data.Neutronics.Transport.AngularDirections; dim = size(xa,2);
wa = data.Neutronics.Transport.AngularWeights;
% Plot 2D quadrature set
if dim == 2
    % Get unit circle information
    R = 1; ang = 0:2*pi/1000:2*pi;
    xcirc = R*cos(ang);
    ycirc = R*sin(ang);
    % Plot quadrature points
    figure(1); hFig = figure(1); ax = gca;
    set(hFig,'Position',[1009 334 561 498]);
    nx = length(wa);
    scatter(xa(:,1), xa(:,2), 120*wa, zeros(nx, 3), 'filled');
    xlabel('\mu');  ax.XTick = (-1:.2:1);
    ylabel('\eta'); ax.YTick = (-1:.2:1);
    set(gca,'fontsize',12);
    % Plot unit circle
    hold on
    plot(xcirc,ycirc,'k')
    plot([-1,1],[0,0],'--k')
    plot([0,0],[-1,1],'--k')
    box on;
    axis square
    hold off
elseif dim == 3
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
end
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