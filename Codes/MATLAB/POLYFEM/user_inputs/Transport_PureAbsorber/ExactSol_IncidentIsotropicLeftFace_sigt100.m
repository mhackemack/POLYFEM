function out = ExactSol_IncidentIsotropicLeftFace_sigt100(xx, ~)
sigt = 100;
x = sigt*xx(:,1);
out = exp(-x) - x.*expint(x);
out(isnan(out)) = 1.0;