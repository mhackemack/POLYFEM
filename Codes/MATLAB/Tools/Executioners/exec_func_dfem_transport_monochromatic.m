%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Execution Functor - DFEM Transport (upwind)
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = exec_func_dfem_transport_monochromatic(ndat, solvdat, mesh, DoF, FE, x, A)
global glob
% Quick Error Checking
% ------------------------------------------------------------------------------
groups = 1:ndat.numberEnergyGroups;
ang_sets = ndat.Transport.AngleSets; nas = length(ang_sets);
if nas ~= 1, error('Only 1 angle set allowed.'); end
rev_str = [];
% Setup Solution Space
% ------------------------------------------------------------------------------
[ndat, flux_out] = setup_solution_space(ndat, mesh, DoF);
DSA_iterations = 0; DSA_time = 0;
% Collect Matrix and RHS and compute angular fluxes
if glob.print_info
    msg = '   Performing Transport Calculation.';
    fprintf([rev_str,msg]);
    rev_str = repmat(sprintf('\b'), 1, length(msg));
end
L = exec_func_LHS_dfem_transport_monochromatic(ndat, mesh, DoF, FE, ang_sets{1}, groups);
rhs = exec_func_RHS_dfem_transport_monochromatic(x, ndat, mesh, DoF, FE, ang_sets{1}, groups);
y = L\rhs;
% Postprocess angular flux solutions
flux_out = add_to_flux(y, ndat, DoF, ang_sets{1}, groups, flux_out);
ndat = compute_partial_boundary_currents(y, ndat, mesh, DoF, ang_sets{1}, groups);
% ndat = add_reflecting_angular_fluxes(y, ndat, mesh, DoF, ang_sets{1}, groups);
% Perform DSA Update
% ------------------------------------------------------------------------------
if ndat.Transport.performDSA
    % Print DSA Calculation Information
    if glob.print_info
        msg = '   Performing DSA Correction Calculation.';
        fprintf([rev_str,msg]);
        rev_str = repmat(sprintf('\b'), 1, length(msg));
    end
    [flux_out, A, DSA_iterations, DSA_time] = perform_DSA_step(ndat, solvdat, mesh, DoF, FE, flux_out, x, A);
end
% Set Outputs
% ------------------------------------------------------------------------------
nout = nargout;
varargout{1} = ndat;
varargout{2} = flux_out;
if nout > 2, varargout{3} = A; end
if nout > 3, varargout{4} = DSA_iterations; end
if nout > 4, varargout{5} = DSA_time; end
% Clear Command Line Text
if glob.print_info, fprintf(rev_str); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ndat, flux_out] = setup_solution_space(ndat, mesh, DoF)
if ndat.Transport.HasReflectingBoundary || ndat.Transport.HasPeriodicBoundary
    ndat.Transport.ReflectingFluxesOld = ndat.Transport.ReflectingFluxes;
    ndat.Transport.PeriodicFluxesOld = ndat.Transport.PeriodicFluxes;
end
ndat.Transport.OutgoingCurrentsOld = ndat.Transport.OutgoingCurrents;
ndat.Transport.IncomingCurrentsOld = ndat.Transport.IncomingCurrents;
ndat = zero_partial_currents(ndat, mesh, DoF);
% Set zero flux moments
ndof = DoF.TotalDoFs;
flux_out = cell(ndat.numberEnergyGroups, ndat.Transport.TotalFluxMoments);
for g=1:ndat.numberEnergyGroups
    for m=1:ndat.Transport.TotalFluxMoments
        flux_out{g,m} = zeros(ndof, 1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flux = add_to_flux(psi, ndat, DoF, angs, groups, flux)
ndof = DoF.TotalDoFs; dofs = (1:ndof)';
ng = length(groups);
na = length(angs);
q_offset = (1:na)*ndof - ndof;
g_offset = (1:ng)*ndof*na - ndof*na;
for q=1:na
    tq = angs(q);
    for g=1:ng
        grp = groups(g);
        qgdofs = dofs + g_offset(g) + q_offset(q);
        for m=1:ndat.Transport.TotalFluxMoments
            d2m = ndat.Transport.discrete_to_moment(m,tq);
            flux{grp,m} = flux{grp,m} + d2m*psi(qgdofs);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ndat = add_reflecting_angular_fluxes(psi, ndat, mesh, DoF, angs, groups)
if ~ndat.Transport.HasReflectingBoundary, return; end
global glob
ndof = DoF.TotalDoFs;
ng = length(groups);
na = length(angs);
q_offset = (1:na)*ndof - ndof;
g_offset = (1:ng)*ndof*na - ndof*na;
for f=1:mesh.TotalBoundaryFaces
    ff = mesh.BoundaryFaces(f);
    fflag = mesh.FaceID(ff);
    fnorm = mesh.FaceNormal(ff,:)';
    for q=1:na
        tq = angs(q);
        adir = ndat.Transport.AngularDirections(tq,:);
        afdot = adir * fnorm;
        if ndat.Transport.BCFlags(fflag) == glob.Reflecting && afdot > 0
            fnodes = DoF.FaceCellNodes{ff,1};
            for g=1:ng
                grp = groups(g);
                fnqg = fnodes + g_offset(g) + q_offset(q);
                ndat.Transport.ReflectingFluxes{ff}{tq}(:,grp) = psi(fnqg);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ndat = compute_partial_boundary_currents(x, ndat, mesh, DoF, angs, groups)
ndof = DoF.TotalDoFs;
ng = length(groups);
na = length(angs);
q_offset = (1:na)*ndof - ndof;
g_offset = (1:ng)*ndof*na - ndof*na;
for f=1:mesh.TotalBoundaryFaces
    ff = mesh.BoundaryFaces(f);
    fnorm = mesh.FaceNormal(ff,:)';
    fnodes = DoF.FaceCellNodes{ff,1};
    for q=1:na
        tq = angs(q);
        adir = ndat.Transport.AngularDirections(tq,:);
        afdot = adir * fnorm;
        wt = ndat.Transport.AngularWeights(tq);
        for g=1:ng
            fnqg = fnodes + g_offset(g) + q_offset(q);
            if afdot > 0
                tv = ndat.Transport.OutgoingCurrents{ff};
                tv(:,g) = tv(:,g) + wt * afdot * x(fnqg);
                ndat.Transport.OutgoingCurrents{ff}(:,g) = tv(:,g);
            elseif afdot < 0
                tv = ndat.Transport.IncomingCurrents{ff};
                tv(:,g) = tv(:,g) + wt * afdot * x(fnqg);
                ndat.Transport.IncomingCurrents{ff}(:,g) = tv(:,g);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ndat = zero_partial_currents(ndat, mesh, DoF)
ng = ndat.numberEnergyGroups;
ndat.Transport.OutgoingCurrents = cell(mesh.TotalFaces, 1);
ndat.Transport.IncomingCurrents = cell(mesh.TotalFaces, 1);
for f=1:mesh.TotalBoundaryFaces
    ff = mesh.BoundaryFaces(f);
    fdn = length(DoF.FaceCellNodes{ff,1});
    ndat.Transport.OutgoingCurrents{ff} = zeros(fdn,ng);
    ndat.Transport.IncomingCurrents{ff} = zeros(fdn,ng);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, A, DSA_it, DSA_time] = perform_DSA_step(ndat, solvdat, mesh, DoF, FE, x, x0, A)
% Form scalar flux differences
for g=1:ndat.numberEnergyGroups
    x0{g} = x{g,1} - x0{g,1};
end
% Switch based on DSA type
if strcmp(ndat.Transport.DSAType, 'MIP')
    [x0, A, DSA_it, DSA_time] = perform_MIP_DSA(ndat, solvdat, mesh, DoF, FE, x0, A);
elseif strcmp(ndat.Transport.DSAType, 'IP')
    [x0, A, DSA_it, DSA_time] = perform_IP_DSA(ndat, solvdat, mesh, DoF, FE, x0, A);
elseif strcmp(ndat.Transport.DSAType, 'M4S')
    [x0, A, DSA_it, DSA_time] = perform_M4S_DSA(ndat, solvdat, mesh, DoF, FE, x0, A);
end
% Add DSA correction
for g=1:ndat.numberEnergyGroups
    x{g,1} = x{g,1} + x0{g};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%