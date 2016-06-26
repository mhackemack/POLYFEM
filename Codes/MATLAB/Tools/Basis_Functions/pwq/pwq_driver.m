clear all; close all; clc;

plot_linear = false;
plot_quadratic = false;
plot_tests = false;
plot_serendipity = false;
plot_serendipity_tests = true;
plot_fem_sol = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                    geometry
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                     Linear basis functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b{1} = @(x,y,sides) t1(x,y,sides)+0.25*tc(x,y,sides);
b{2} = @(x,y,sides) t2(x,y,sides)+0.25*tc(x,y,sides);
b{3} = @(x,y,sides) t3(x,y,sides)+0.25*tc(x,y,sides);
b{4} = @(x,y,sides) t4(x,y,sides)+0.25*tc(x,y,sides);

n=11;
xx=linspace(0,1,n); yy=xx;
[X,Y]=meshgrid(xx,yy);

% create plots
if plot_linear
    for k=1:length(sides)
        figure(k);
        plot_basis_function(b{k},X,Y,sides);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                     Quadratic basis functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create quadratic space
for i=1:4
    for j=1:4
        mu{i,j} = @(x,y,sides) b{i}(x,y,sides)*b{j}(x,y,sides);
    end
end

% create plots
if plot_tests
    figure(111);
    test_same = @(x,y,sides) mu{1,3}(x,y,sides)-mu{2,4}(x,y,sides);
    plot_basis_function(test_same,X,Y,sides);
end

% create plots
if plot_quadratic
    figure(111);
    test_same = @(x,y,sides) mu{1,3}(x,y,sides)-mu{2,4}(x,y,sides);
    plotk=0;
    for i=1:4
        for j=i:4
            %         if (i==2 && j==4)
            %             continue;
            %         end
            k=k+1;
            [i j k]
            figure(10+k);
            plot_basis_function(mu{i,j},X,Y,sides);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                     Testing Quadratic
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create plots
if plot_tests
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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                     Serendipity Conversion and Testing
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build A matrix
% V: (1,1), (2,2), (3,3), (4,4)
% E: (1,2), (2,3), (3,4), (4,1)
% D: (1,3) then (2,4) for diagonal terms
% ------------------------------------------------------------------------------
Ap = [-1,0;0,-1;-1,0;0,-1;0.5*ones(4,2)];
% Ap = [-1,0;0,-1;-1,0;0,-1;ones(4,2)];
A = [eye(8),Ap];
% Build B matrix
% ------------------------------------------------------------------------------
BB = diag(-1*ones(4,1)) - diag(ones((4-1),1),-1);
BB(1,4) = -1;
B = [eye(4),BB;zeros(4),2*eye(4)];
% Build Serendipity Functions
% ------------------------------------------------------------------------------
% First, consolidate repeating quadratic functions
% i,j ordering for i=1:4, j=1:4 is
%      i     j     k     Mike's row index in A
%      1     1     1       1
%      1     2     2       5
%      1     3     3       9
%      1     4     4       8
%      2     1     5
%      2     2     6       2
%      2     3     7       6
%      2     4     8      10
%      3     1     9
%      3     2    10
%      3     3    11       3
%      3     4    12       7
%      4     1    13
%      4     2    14
%      4     3    15
%      4     4    16       4
% usein_A=zeros(16,2);
% usein_A(:,1) = 1:16; % k from above
% usein_A(:,2) = [ 1 5 9 8 0 2 6 10 0 0 3 7 0 0 0 4];

% i,j ordering for i=1:4, j=i:4 is
%      i     j     k   Mike's row index in A
%      1     1     1   1
%      1     2     2   5
%      1     3     3   9
%      1     4     4   8
%      2     2     5   2
%      2     3     6   6
%      2     4     7  10
%      3     3     8   3
%      3     4     9   7
%      4     4    10   4
usein_A=zeros(10,2);
usein_A(:,1) = 1:10; % k from above
usein_A(:,2) = [ 1 5 9 8 2 6 10 3 7 4];

disp('Creating serendipity.')
for ixi=1:8
    % create zero function first
    xi{ixi} = @(x,y,sides) 0*x+0*y;
    % grab row of A but permute columns
    coef = A(ixi,usein_A(:,2));
    k=0;
    for i=1:4
        for j=i:4
            k=k+1;
            % because j goes from i to 4 and not 1 to 4, the factor 2 is needed
            fac=1;
            if abs(i-j) == 1 % bottom, right, and top faces
                fac=2;
            elseif i==1 && j==4 % left face
                fac=2;
            elseif abs(i-j) > 1
                fac=2;
            end
%             if abs(i-j) == 1 % bottom, right, and top faces
%                 fac=2;
%             elseif i==1 && j==4 % left face
%                 fac=2;
%             elseif abs(i-j) > 1
%                 fac=1;
%             end
            % only do [xi]=A[mu] when the entries in A are not zero
            if abs(coef(k)) > eps
                xi{ixi}= @(x,y,sides) xi{ixi}(x,y,sides) + fac*coef(k)*mu{i,j}(x,y,sides);
            end
        end
    end
end

if plot_serendipity
    % test serendipity
    disp('Plotting serendipity.')
    for ixi=1:8
        fprintf('%d/8 \n',ixi);
        figure(800+ixi);
        plot_basis_function(xi{ixi},X,Y,sides);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                     Testing serendipity
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create plots
if plot_serendipity_tests
    % test unity
    disp('Testing serendipity unity.')
    test_unity = @(x,y,sides) -1 + 0*x + 0*y;
    for ixi=1:8
        test_unity = @(x,y,sides) test_unity(x,y,sides) + xi{ixi}(x,y,sides);
    end
    figure(920);
    plot_basis_function(test_unity,X,Y,sides);
    % test x-linearity
    disp('Testing serendipity x-linearity.')
    test_xlin = @(x,y,sides) -1*x + 0*y;
    for ixi=1:8
        if ixi <=4
            xcoef = pt{ixi}(1);
        else
            xcoef = (pt{ixi-4}(1)+pt{mod(ixi-4,4)+1}(1))/2;
        end
        test_xlin = @(x,y,sides) test_xlin(x,y,sides) + xcoef*xi{ixi}(x,y,sides);
    end
    figure(921);
    plot_basis_function(test_xlin,X,Y,sides);
    % test y-linearity
    disp('Testing serendipity y-linearity.')
    test_ylin = @(x,y,sides) 0*x + -1*y;
    for ixi=1:8
        if ixi <=4
            ycoef = pt{ixi}(2);
        else
            ycoef = (pt{ixi-4}(2)+pt{mod(ixi-4,4)+1}(2))/2;
        end
        test_ylin = @(x,y,sides) test_ylin(x,y,sides) + ycoef*xi{ixi}(x,y,sides);
    end
    figure(922);
    plot_basis_function(test_ylin,X,Y,sides);
    % test xy-quadratic
    disp('Testing serendipity xy-quadratic.')
    test_xy = @(x,y,sides) -1*x.*y;
    for ixi=1:8
        if ixi <=4
            coef = pt{ixi}(1)*pt{ixi}(2);
        else
            coef = (pt{ixi-4}(1)+pt{mod(ixi-4,4)+1}(1))/2*(pt{ixi-4}(2)+pt{mod(ixi-4,4)+1}(2))/2;
        end
        test_xy = @(x,y,sides) test_xy(x,y,sides) + coef*xi{ixi}(x,y,sides);
    end
    figure(923);
    plot_basis_function(test_xy,X,Y,sides);
    % test xx-quadratic
    disp('Testing serendipity xx-quadratic.')
    test_xx = @(x,y,sides) -1*x.*x;
    for ixi=1:8
        if ixi <=4
            coef = pt{ixi}(1)*pt{ixi}(1);
        else
            coef = pt{ixi-4}(1)*pt{mod(ixi-4,4)+1}(1);
        end
        test_xx = @(x,y,sides) test_xx(x,y,sides) + coef*xi{ixi}(x,y,sides);
    end
    figure(924);
    plot_basis_function(test_xx,X,Y,sides);
    % test yy-quadratic
    disp('Testing serendipity yy-quadratic.')
    test_yy = @(x,y,sides) -1*y.*y;
    for ixi=1:8
        if ixi <=4
            coef = pt{ixi}(2)*pt{ixi}(2);
        else
            coef = pt{ixi-4}(2)*pt{mod(ixi-4,4)+1}(2);
        end
        test_yy = @(x,y,sides) test_yy(x,y,sides) + coef*xi{ixi}(x,y,sides);
    end
    figure(925);
    plot_basis_function(test_yy,X,Y,sides);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                     Create a FEM solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4 vertices followed by 4 mid-edges
sol_coef=[ 1 0 0 1 0.2 0 0.2 1];
fem_sol_serendipity = @(x,y,sides) 0*x+0*y;
for ixi=1:8
    fem_sol_serendipity = @(x,y,sides) fem_sol_serendipity(x,y,sides) + sol_coef(ixi)*xi{ixi}(x,y,sides);
end

if plot_fem_sol
    disp('Testing fem_sol_serendipity.')
    figure(900);
    plot_basis_function(fem_sol_serendipity,X,Y,sides);
end



% coef ordering based on k using i=1:4, j=i:4
sol_coef=[ 1 0.2 0.2 1 0 0 0.2 0 0.2 1];
fem_sol_quad = @(x,y,sides) 0*x+0*y;
k=0;
for i=1:4
    for j=i:4
        k=k+1;
        % because j goes from i to 4 and not 1 to 4, the factor 2 is needed
        fac=1;
        if i~=j
            fac=2;
        end
        fem_sol_quad = @(x,y,sides) fem_sol_quad(x,y,sides) + fac*sol_coef(k)*mu{i,j}(x,y,sides);
    end
end

if plot_fem_sol
    disp('Testing fem quad solution.')
    figure(901);
    plot_basis_function(fem_sol_quad,X,Y,sides);
end