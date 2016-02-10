function out = ExactSol_IncidentLeftFace_2D_45deg(xx,~)
x = xx(:,1); y = xx(:,2);
val = cos(pi/4);
out = exp(-x./val);
out(x>y) = 0.;