function out = ExactSol_IncidentLeftTopFace_2D_45degDown_LS4_sigt10(xx,~)
x = xx(:,1); y = xx(:,2); nx = length(x);
val = cos(pi/4); valdenom = 3.500211745815407e-01;
sigt = 10;
out = zeros(nx,1);
for i=1:nx
    if y(i)<(1-x(i))
        out(i) = (1/val)*exp(-sigt*x(i)/valdenom);
    else
        out(i) = (1/val)*exp(-sigt*(1-y(i))/valdenom);
    end
end