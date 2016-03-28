function out = ExactSol_IncidentIsotropicLeftFace_sigt10(xx, ~)
sigt = 10;
x = sigt*xx(:,1);
out = exp(-x) - x.*expint(x);
out(isnan(out)) = 1.0;