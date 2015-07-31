%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Execute LD Sweep Chunk
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    Invert Transport Operator by Sweeping - upwind scheme is
%                   strongly enforced.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flux = exec_func_transport_sweep_LD(ndat, mesh, DoF, FE, x, m, groups)
% Process Input Space
% -------------------
global glob
dim = mesh.Dimension;
ndof = DoF.TotalDoFs;
nldof = FE.LDNumDoFs;
ng = length(groups);
angs = ndat.Transport.AngleSets{m}; na = length(angs);
angdirs = ndat.Transport.AngularDirections;
angNorm = ndat.Transport.AngQuadNorm;
m2d = ndat.Transport.moment_to_discrete;
% Get Sweep Information
% ---------------------
sweep = ndat.Transport.Sweeping;
CellOrder = sweep.CellSweepOrder{m};
USFaces = sweep.UpstreamFaces{m};
