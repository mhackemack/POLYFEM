%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Jacobi Energy Collapsing
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
function data = collapse_jacobi_xs(data,xsin,xsout,accid)
global glob
% Retrieve some problem information
% ------------------------------------------------------------------------------
nm = data.problem.NumberMaterials;
grps = data.Acceleration.Info(accid).Groups; ngrps = length(grps);
txs = data.XS(xsin).TotalXS(:,grps);
sxs = data.XS(xsin).ScatteringXS(:,grps,grps,1);
% Generate Eigenshape for energy collapse
% ------------------------------------------------------------------------------
data.Acceleration.Info(accid).ErrorShape = zeros(nm,ngrps); y = cell(nm,1);
for m=1:nm
    A = diag(txs(m,:))\squeeze(sxs(m,:,:));
    [y{m},~,~] = power_method(A,ones(ngrps,1),2000,1e-15);
    y{m} = y{m} / sum(y{m});
    data.Acceleration.Info(accid).ErrorShape(m,:) = y{m};
end
% Populate Acceleration XS
% ------------------------------------------------------------------------------
if data.Acceleration.Info(accid).AccelerationType == glob.Accel_WGS_DSA
    data.XS(xsout).DiffXS = zeros(nm,1);
    data.XS(xsout).AbsorbXS = zeros(nm,1);
    data.XS(xsout).ScatteringXS = data.XS(xsin).ScatteringXS;
    for m=1:nm
        if data.Transport.PnOrder == 0
            DC = 1./(3*txs(m,:));
        elseif data.Transport.PnOrder > 0
            tsxs = sum(squeeze(sxs(m,:,:,2)),2)';
            DC = 1./(3*(txs(m,:) - tsxs));
        end
        data.XS(xsout).TotalXS(m)  = txs(m,:)*y{m};
        data.XS(xsout).DiffXS(m) = DC*y{m};
        data.XS(xsout).AbsorbXS(m) = txs(m,:)*y{m} - sum(squeeze(sxs(m,:,:))*y{m});
    end
elseif data.Acceleration.Info(accid).AccelerationType == glob.Accel_WGS_TSA
    
else
    error('Incorrect acceleration type for this function call.');
end
% Set Boundary Conditions
% ------------------------------------------------------------------------------
data.XS(xsout).BCFlags = data.XS(xsin).BCFlags; 
data.XS(xsout).BCVals  = data.XS(xsin).BCVals;
return