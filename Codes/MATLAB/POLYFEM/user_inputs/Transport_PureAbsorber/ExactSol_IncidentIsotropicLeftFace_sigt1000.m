function out = ExactSol_IncidentIsotropicLeftFace_sigt1000(xx, ~)
sigt = 1000;
x = sigt*xx(:,1);
out = exp(-x) - x.*expint(x);
out(isnan(out)) = 1.0;