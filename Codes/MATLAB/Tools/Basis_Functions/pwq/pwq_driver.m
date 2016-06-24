clear all; close all; clc;

%
%  4-------3
%  |       |
%  |   c   |
%  |       |
%  1-------2
%
pt{1}=[0; 0];
pt{2}=[1; 0];
pt{3}=[1; 1];
pt{4}=[0; 1];
% (pt1+pt2+pt3+pt)/4;
ptc=[0;0];
for k=1:length(pt)
    ptc=ptc+pt{k}/4;
end

triangle{1}=[pt{1} pt{2} ptc];
triangle{2}=[pt{2} pt{3} ptc];
triangle{3}=[pt{3} pt{4} ptc];
triangle{4}=[pt{4} pt{1} ptc];
for k=1:length(triangle)
    sides{k}.x = triangle{k}(1,:);
    sides{k}.y = triangle{k}(2,:);
end

b{1} = @(x,y,sides) t1(x,y,sides)+0.25*tc(x,y,sides);
b{2} = @(x,y,sides) t2(x,y,sides)+0.25*tc(x,y,sides);
b{3} = @(x,y,sides) t3(x,y,sides)+0.25*tc(x,y,sides);
b{4} = @(x,y,sides) t4(x,y,sides)+0.25*tc(x,y,sides);

n=20;
xx=linspace(0,1,n); yy=xx;
[X,Y]=meshgrid(xx,yy);

% for k=1:length(sides)
%     figure(k);
%     plot_basis_function(b{k},X,Y,sides);
% end

% create quadratic space
for i=1:4
    for j=1:4
        mu{i,j} = @(x,y,sides) b{i}(x,y,sides)*b{j}(x,y,sides);
    end
end

figure(111);
test_same = @(x,y,sides) mu{1,3}(x,y,sides)-mu{2,4}(x,y,sides);
plot_basis_function(test_same,X,Y,sides);


% k=0;
% for i=1:4
%     for j=i:4
% %         if (i==2 && j==4)
% %             continue;
% %         end
%         k=k+1;
%         [i j k]
%         figure(10+k);
%         plot_basis_function(mu{i,j},X,Y,sides);
%     end
% end

% test unity
disp('Testing unity.')
test_unity = @(x,y,sides) -1 + 0*x + 0*y;
for i=1:4
    for j=1:4
        % because j goes from 1 to 4 and not i to 4, the factor 2 is not needed
        test_unity = @(x,y,sides) test_unity(x,y,sides) + mu{i,j}(x,y,sides);
    end
end
figure(200);
plot_basis_function(test_unity,X,Y,sides);

% test linearity
disp('Testing x-linearity.')
test_linx = @(x,y,sides) -1*x+0*y;
for i=1:4
    for j=1:4
        % because j goes from 1 to 4 and not i to 4, the factor 2 is not needed
        xcoef = (pt{i}(1)+pt{j}(1))/2;
        test_linx = @(x,y,sides) test_linx(x,y,sides) + xcoef*mu{i,j}(x,y,sides);
    end
end
figure(300);
plot_basis_function(test_linx,X,Y,sides);

disp('Testing y-linearity.')
test_liny = @(x,y,sides) 0*x-1*y;
for i=1:4
    for j=1:4
        % because j goes from 1 to 4 and not i to 4, the factor 2 is not needed
        ycoef = (pt{i}(2)+pt{j}(2))/2;
        test_liny = @(x,y,sides) test_liny(x,y,sides) + ycoef*mu{i,j}(x,y,sides);
    end
end
figure(400);
plot_basis_function(test_liny,X,Y,sides);

% test quadratic
disp('Testing xx-quadratic.')
test_quad_xx = @(x,y,sides) -x.^2;
for i=1:4
    for j=1:4
        % because j goes from 1 to 4 and not i to 4, the factor 2 is not needed
        xxcoef = (pt{i}(1)*pt{j}(1) + pt{j}(1)*pt{i}(1))/2;
        test_quad_xx = @(x,y,sides) test_quad_xx(x,y,sides) + xxcoef*mu{i,j}(x,y,sides);
    end
end
figure(500); plot_basis_function(test_quad_xx,X,Y,sides);

% test quadratic
disp('Testing xy-quadratic.')
test_quad_xy = @(x,y,sides) -x.*y;
for i=1:4
    for j=1:4
        % because j goes from 1 to 4 and not i to 4, the factor 2 is not needed
        xycoef = (pt{i}(1)*pt{j}(2) + pt{j}(1)*pt{i}(2))/2;
        test_quad_xy = @(x,y,sides) test_quad_xy(x,y,sides) + xycoef*mu{i,j}(x,y,sides);
    end
end
figure(600); plot_basis_function(test_quad_xy,X,Y,sides);

% test quadratic
disp('Testing yy-quadratic.')
test_quad_yy = @(x,y,sides) -y.^2;
for i=1:4
    for j=1:4
        % because j goes from 1 to 4 and not i to 4, the factor 2 is not needed
        yycoef = (pt{i}(2)*pt{j}(2) + pt{j}(2)*pt{i}(2))/2;
        test_quad_yy = @(x,y,sides) test_quad_yy(x,y,sides) + yycoef*mu{i,j}(x,y,sides);
    end
end
figure(700); plot_basis_function(test_quad_yy,X,Y,sides);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%                     Serendipity Conversion and Testing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build A matrix
% --------------
Ap = [-1,0;0,-1;-1,0;0,-1;.5*ones(4,2)];
A = [eye(8),Ap];
% Build B matrix
% --------------
BB = diag(-1*ones(4,1)) - diag(ones((4-1),1),-1);
BB(1,4) = -1;
B = [eye(4),BB;zeros(4),2*eye(4)];
% Build Serendipity Functions
% ---------------------------
