function [A, P] = polyareaperimeter(xin)

x = xin(:,1);
y = xin(:,2);
n = size(xin,1);
% number of vertices
% [ x, ns ] = shiftdim( x );
% [ y, ns ] = shiftdim( y );
% temporarily shift data to mean of vertices for improved accuracy
xm = mean(x);
ym = mean(y);
x = x - xm*ones(n,1); 
y = y - ym*ones(n,1); 
% delta x and delta y
dx = x( [ 2:n 1 ] ) - x;
dy = y( [ 2:n 1 ] ) - y;
% area and perimeter
A = sum( y.*dx - x.*dy )/2;
P = sum( sqrt( dx.*dx +dy.*dy ) ); 

return