function out = ExactSol_IncidentLeftFace_2D_45degDown_manual(xx,~)
x = xx(:,1); y = xx(:,2);
val = cos(pi/4);
out = (1/val)*exp(-x./val);
out(y>1-x) = 0.;