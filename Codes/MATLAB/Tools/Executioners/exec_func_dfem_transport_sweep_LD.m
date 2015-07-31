%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Execution Functor - LDFEM Transport with Sweeps
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    Invert Transport Operator by Sweeping - upwind scheme is
%                   strongly enforced.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function varargout = exec_func_dfem_transport_sweep_LD(ndat, solvdat, mesh, DoF, FE, x, A)
global glob
% Setup Solution Space
