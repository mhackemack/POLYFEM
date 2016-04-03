function [V,F] = RegularPolygon(N, rad)
if N<3, error('Cannot form a polygon with N<3.'); end
% generate the vertices.
theta = 0:2*pi/N:(2*pi-2*pi/N);
Pos = (rad*exp(-1i*theta))';

V = [real(Pos),imag(Pos)];
F = cell(N,1);
for i=1:N
    F{i} = [i,mod(i,N)+1];
end