function out = ExactSol_IncidentIsotropicLeftFace(xx, ~)
x = xx(:,1);
out = exp(-x) - x.*expint(x);
out(isnan(out)) = 1.0;