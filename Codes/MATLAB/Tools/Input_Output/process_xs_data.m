%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Check Cross Sections
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2015
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function XS = process_xs_data(XS)
% Set the XS booleans
nxs = length(XS);
tot_bools    = false(nxs,1);
diff_bools   = false(nxs,1);
abs_bools    = false(nxs,1);
fiss_bools   = false(nxs,1);
nu_bools     = false(nxs,1);
chi_bools    = false(nxs,1);
scatt_bools  = false(nxs,1);
src_bools    = false(nxs,1);
% Loop through all cross section structs
for i=1:nxs
% Set Default Booleans in struct
% ------------------------------------------------------------------------------
XS(i).HasFission = false;
XS(i).HasScattering = false;
XS(i).HasExtSource = false;
% Check Total Cross Section
% ------------------------------------------------------------------------------
if ~isempty(XS(i).TotalXS)
    if nnz(XS(i).TotalXS) > 0, tot_bools(i) = true; end
end
% Check Diffusion Coefficients - we will set this automatically if needed
% ------------------------------------------------------------------------------
if ~isempty(XS(i).DiffXS)
    % Check for any non-zero entries
    if all(XS(i).DiffXS(:))
        diff_bools(i) = true;
    else
        diff_bools(i) = false;
    end
else
    if tot_bools(i)
        XS(i).DiffXS = 1/3/XS(i).TotalXS;
        diff_bools(i) = true;
    else
        diff_bools(i) = false;
    end
end
% Check Absorption Cross Section
% ------------------------------------------------------------------------------
if ~isempty(XS(i).AbsorbXS)
    if nnz(XS(i).AbsorbXS) > 0, abs_bools(i) = true; end
end
% Check Fission Cross Section
% ------------------------------------------------------------------------------
if ~isempty(XS(i).FissionXS)
    if nnz(XS(i).FissionXS) > 0, fiss_bools(i) = true; end
end
% Check NuBar Cross Section
% ------------------------------------------------------------------------------
if ~isempty(XS(i).NuBar)
    if nnz(XS(i).NuBar) > 0, nu_bools(i) = true; end
end
% Check Fission Emission Spectrum (Chi)
% ------------------------------------------------------------------------------
if ~isempty(XS(i).FissSpec)
    if nnz(XS(i).FissSpec) > 0, chi_bools(i) = true; end
end
% Check External Source
% ------------------------------------------------------------------------------
if ~isempty(XS(i).ExtSource)
    % Check if an array - corresponds to non-MMS values
    if isnumeric(XS(i).ExtSource)
        if nnz(XS(i).ExtSource) > 0, src_bools(i) = true; end
    end
    % Check cell array for function handles - corresponds to MMS
    if iscell(XS(i).ExtSource)
        for j=1:size(XS(i).ExtSource,1)
            for k=1:size(XS(i).ExtSource,2)
                if ~isa(XS(i).ExtSource{j,k},'function_handle')
                    error('External Source is not a function handle for XS %d with index (%d,%d).',i,j,k);
                end
            end
        end
    end
end
% Set Boolean Flags in 'data' struct
% ------------------------------------------------------------------------------
if fiss_bools(i) && chi_bools(i) && nu_bools(i), XS(i).HasFission = true; end
if scatt_bools(i), XS(i).HasScattering = true; end
if src_bools(i), XS(i).HasExtSource = true; end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%