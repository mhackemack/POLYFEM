function out = get_orthogonal_length(dim, cVol, nv, fA, cSA)

if dim == 1
    out = cVol;
elseif dim == 2
    if nv == 3
        out = 2*cVol/fA;
    elseif nv == 4
        out = cVol/fA;
    elseif nv > 4 && mod(nv,2) == 0
        out = 4*cVol/cSA;
    elseif nv > 4 && mod(nv,2) ~= 0
        out = 2*cVol/cSA + sqrt((2*cVol)/(nv*sin(2*pi/nv)));
    end
elseif dim == 3
    if nv == 4
        out = 3*cVol/fA;
    elseif nv == 8
        out = cVol/fA;
    else
        out = 6*cVol/cSA;
    end
end
