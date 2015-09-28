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
function varargout = exec_func_dfem_transport_upwind_Rev1(data,xsid,qid,groups,mesh,DoF,FE)
global glob
% Setup Solution Space
% ------------------------------------------------------------------------------
[data, flux_out] = setup_solution_space(data, mesh, DoF);
% Loop through Angle Sets
% ------------------------------------------------------------------------------
ang_sets = data.Quadrature.AngleSets; nas = length(ang_sets);
rev_str = [];
for m=1:nas
    % Print Angle Set Iteration
    if glob.print_info
        msg = sprintf('   Calculating Flux for Angle Set: %d of %d',m,nas);
        fprintf([rev_str,msg]);
        rev_str = repmat(sprintf('\b'), 1, length(msg));
    end
    % Collect Matrix and RHS and compute angular fluxes
    L = exec_func_LHS_dfem_transport_upwind_Rev1(data, mesh, DoF, FE, ang_sets{m}, groups);
    rhs = exec_func_RHS_dfem_transport_Rev1(x, data, mesh, DoF, FE, ang_sets{m}, groups);
    y = L\rhs;
    % Postprocess angular flux solutions
    flux_out = add_to_flux(y, data, DoF, ang_sets{m}, groups, flux_out);
%     ndat = compute_partial_boundary_currents(y, ndat, mesh, DoF, ang_sets{m}, groups);
    data = add_reflecting_angular_fluxes(y, data, mesh, DoF, ang_sets{m}, groups);
end
% Set Outputs
% ------------------------------------------------------------------------------
varargout{1} = data;
varargout{2} = flux_out;
% Clear Command Line Text
if glob.print_info, fprintf(rev_str); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [data, flux_out] = setup_solution_space(data, mesh, DoF)
if data.Fluxes.HasReflectingBoundary || data.Fluxes.HasPeriodicBoundary
    data.Fluxes.ReflectingFluxesOld = data.Fluxes.ReflectingFluxes;
    data.Fluxes.PeriodicFluxesOld = data.Fluxes.PeriodicFluxes;
end
% data.Fluxes.OutgoingCurrentsOld = data.Fluxes.OutgoingCurrents;
% data.Fluxes.IncomingCurrentsOld = data.Fluxes.IncomingCurrents;
% data = zero_partial_currents(data, mesh, DoF);
% Set zero flux moments
ndof = DoF.TotalDoFs;
flux_out = cell(data.Groups.numberEnergyGroups, data.Fluxes.TotalFluxMoments);
for g=1:data.Groups.numberEnergyGroups
    for m=1:data.Fluxes.TotalFluxMoments
        flux_out{g,m} = zeros(ndof, 1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function flux = add_to_flux(psi, data, DoF, angs, groups, flux)
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
        for m=1:data.Fluxes.TotalFluxMoments
            d2m = data.Quadrature.discrete_to_moment(m,tq);
            flux{grp,m} = flux{grp,m} + d2m*psi(qgdofs);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = add_reflecting_angular_fluxes(psi, data, mesh, DoF, angs, groups)
if ~data.Quadrature.HasReflectingBoundary, return; end
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
        adir = data.Quadrature.AngularDirections(tq,:);
        afdot = adir * fnorm;
        if data.XS.BCFlags(fflag) == glob.Reflecting && afdot > 0
            fnodes = DoF.FaceCellNodes{ff,1};
            for g=1:ng
                grp = groups(g);
                fnqg = fnodes + g_offset(g) + q_offset(q);
                data.Fluxes.ReflectingFluxes{ff}{tq}(:,grp) = psi(fnqg);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = compute_partial_boundary_currents(x, data, mesh, DoF, angs, groups)
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
        adir = data.Quadrature.AngularDirections(tq,:);
        afdot = adir * fnorm;
        wt = data.Quadrature.AngularWeights(tq);
        for g=1:ng
            fnqg = fnodes + g_offset(g) + q_offset(q);
            if afdot > 0
                tv = data.Fluxes.OutgoingCurrents{ff};
                tv(:,g) = tv(:,g) + wt * afdot * x(fnqg);
                data.Fluxes.OutgoingCurrents{ff}(:,g) = tv(:,g);
            elseif afdot < 0
                tv = data.Fluxes.IncomingCurrents{ff};
                tv(:,g) = tv(:,g) + wt * afdot * x(fnqg);
                data.Fluxes.IncomingCurrents{ff}(:,g) = tv(:,g);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = zero_partial_currents(data, mesh, DoF)
ng = data.Groups.numberEnergyGroups;
data.Fluxes.OutgoingCurrents = cell(mesh.TotalFaces, 1);
data.Fluxes.IncomingCurrents = cell(mesh.TotalFaces, 1);
for f=1:mesh.TotalBoundaryFaces
    ff = mesh.BoundaryFaces(f);
    fdn = length(DoF.FaceCellNodes{ff,1});
    data.Fluxes.OutgoingCurrents{ff} = zeros(fdn,ng);
    data.Fluxes.IncomingCurrents{ff} = zeros(fdn,ng);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%