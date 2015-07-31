%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Triangle Area (Heron's Formula)
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%
%   Description:    MATLAB function to calculate a triangle's area by use of
%                   Heron's Formula.
%                   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Notes:   
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = triangle_area(verts)
v1 = verts(1,:);
v2 = verts(2,:);
v3 = verts(3,:);
a = norm(v2-v1);
b = norm(v3-v1);
c = norm(v2-v3);
s = (a+b+c)/2;
out = sqrt(s*(s-a)*(s-b)*(s-c));