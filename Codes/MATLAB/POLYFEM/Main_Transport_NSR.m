%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Title:          Spectral Radius Transport Run Script
%
%   Author:         Michael W. Hackemack
%   Institution:    Texas A&M University
%   Year:           2014
%   
%   Description:    
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Note(s):        
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clear Project Space
% ------------------------------------------------------------------------------
if exist('pbool', 'var')
    clearvars -except pbool
else
    clear; pbool = false;
end
clc; close all; format long e
if ~pbool, fpath = get_path(); addpath(fpath); pbool = true; end
% Populate global space
% ------------------------------------------------------------------------------
global glob
glob = get_globals('Office');
glob.print_info = true;
addpath([glob.input_path,'Transport_NSR']);
% Begin user input section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bf, quad, bc
% bf_name = {'WACHSPRESS','MV'};
% bf_name = {'WACHSPRESS','MV','MAXENT'};
bf_name = {'PWLD'};
fdeg = [1];
q_type = 'LS'; sn_levels = [2,4];
bc_type = 'Vacuum';
% geometry
dim = 3; m_type = 'tet';
% dx_num_start = 2; L = 1;
dx_num_start = 21; L = 1;
dx_start = linspace(0,L,dx_num_start);
ar = 10;
% xs
c = 0.9999;
mfp_lower = 1; mfp_upper = 41;
mfp_min = 0; mfp_max = 3;
% mfp_lower = 0; mfp_upper = 41;
% mfp_min = -1; mfp_max = 3;
mfp_vals = logspace(mfp_min, mfp_max, mfp_upper);
% DSA
diff_type = 'M4S'; C_IP = [4];
% End user input section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load some initial data
% ----------------------
geom = load_geometry_input(dim, m_type, dx_start, ar, 1);

% Loop through problem space and execute test suite
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mfp_tot = mfp_lower + mfp_upper;
sn_num = length(sn_levels);
C_num = length(C_IP);
txs = 0;
t_max_dim = (max(geom.CellVolume))^(1/dim);
dname = 'outputs/Transport_NSR/';
dname = [dname, diff_type, '_', bc_type, '_', m_type, '/'];
data.Neutronics.Transport.DSAType = diff_type;
% MIP and IP DSA schemes
% ------------------------------------------------------------------------------
if strcmpi(data.Neutronics.Transport.DSAType,'mip') || strcmpi(data.Neutronics.Transport.DSAType,'ip')
% Allocate Memory Space
SI_err_L2   = cell(sn_num, C_num, mfp_tot);
SI_err_inf  = cell(sn_num, C_num, mfp_tot);
SI_norm_L2  = cell(sn_num, C_num, mfp_tot);
SI_norm_inf = cell(sn_num, C_num, mfp_tot);
SI_iters = zeros(sn_num, C_num, mfp_tot);
NSR_err  = zeros(sn_num, C_num, mfp_tot);
NSR_norm = zeros(sn_num, C_num, mfp_tot);
% Loop through finite element orders
for f=1:length(fdeg)
    disp(['-> Finite Element Degree: ',num2str(f),' of ', num2str(length(fdeg))])
    for b=1:length(bf_name)
        disp(['  -> Basis Function: ',num2str(b),' of ', num2str(length(bf_name))])
        data = load_user_input(dim, bf_name{b}, fdeg(f), bc_type);
        data.Neutronics.Transport.DSAType = diff_type;
        % Loop through quadrature
        rev_str = [];
        for m=1:sn_num
            disp(['    -> Quadrature Set: ',num2str(m),' of ', num2str(sn_num)])
            data = load_quad_input( data, q_type, sn_levels(m) );
            % Loop through IP constants
            for i=1:C_num
                disp(['      -> IP Constant: ',num2str(i),' of ', num2str(C_num)])
                mfp = zeros(mfp_tot, 1);
                data.Neutronics.IP_Constant = C_IP(i);
                %         if dim == 3, data.Neutronics.IP_Constant = 1e-3*data.Neutronics.IP_Constant; end
                dx_num = dx_num_start;
                geom = load_geometry_input(dim, m_type, dx_start, ar, 1);
                tc = 0;
                % First loop through upper mfp values in reverse order
                for j=mfp_upper:-1:1
                    tc = tc + 1;
                    msg = sprintf('        -> MFP Number: %d of %d',tc,mfp_tot);
                    fprintf([rev_str,msg]);
                    rev_str = repmat(sprintf('\b'), 1, length(msg));
                    
                    jj = mfp_lower + j;
                    mfp(jj) = mfp_vals(j);
                    txs = mfp_vals(j) / t_max_dim;
                    data = load_xs_input( data, txs, c);
                    [data, geom] = process_input_data(data, geom);
                    data = cleanup_neutronics_input_data(data, geom);
                    % Execute problem
                    [data, sol, ~, ~, ~] = execute_problem(data, geom);
                    % Collect statistics
                    SI_iters(m,i,jj) = sol.iter;
                    SI_err_L2{m,i,jj} = sol.error_L2;
                    SI_err_inf{m,i,jj} = sol.error_inf;
                    SI_norm_L2{m,i,jj} = sol.norm_L2;
                    SI_norm_inf{m,i,jj} = sol.norm_inf;
                end
                % Next loop through mesh refinement steps
                jj = 1;
                for j=mfp_lower:-1:1
                    jj = jj + 1;
                    tc = tc + 1;
                    msg = sprintf('      -> MFP Number: %d of %d',tc,mfp_tot);
                    fprintf([rev_str,msg]);
                    rev_str = repmat(sprintf('\b'), 1, length(msg));
                    
                    dx_num = (dx_num-1)*2+1;
                    dx = linspace(0,L,dx_num);
                    geom = load_geometry_input(dim, m_type, dx, ar, jj);
                    mfp(j) = txs * (max(geom.CellVolume))^(1/dim);
                    [data, geom] = process_input_data(data, geom);
                    data = cleanup_neutronics_input_data(data, geom);
                    % Execute problem
                    glob.print_info = true;
                    [data, sol, ~, ~, ~] = execute_problem(data, geom);
                    glob.print_info = false;
                    % Collect statistics
                    SI_iters(m,i,j) = sol.iter;
                    SI_err_L2{m,i,j} = sol.error_L2;
                    SI_err_inf{m,i,j} = sol.error_inf;
                    SI_norm_L2{m,i,j} = sol.norm_L2;
                    SI_norm_inf{m,i,j} = sol.norm_inf;
                end
                % Process statistics
                for j=1:mfp_tot
                    it = SI_iters(m,i,j); counter = 0;
                    terr_L2 = 0; t_norm_L2 = 0;
                    if it < 4, continue; end
                    NSR_err(m,i,j) = SI_err_L2{m,i,j}(end)/SI_err_L2{m,i,j}(end-1);
                    NSR_norm(m,i,j) = SI_norm_L2{m,i,j}(end)/SI_norm_L2{m,i,j}(end-1);
                end
                % Save off information
                if ~isequal(exist(dname, 'dir'),7), mkdir(dname); end
                q = sn_levels(m); C = C_IP(i);
                if dim == 1
                    fname = sprintf('%s%d_%s%d_C=%d',bf_name{b},fdeg(f),q_type,q,C);
                else
                    fname = sprintf('%s%d_%s%d_C=%d_AR=%d',bf_name{b},fdeg(f),q_type,q,C,ar);
                end
                mat_out = [mfp, squeeze(NSR_norm(m,i,:)), squeeze(NSR_err(m,i,:))];
                dlmwrite([dname,fname,'.dat'],mat_out,'precision','%14.8e');
            end
        end
    end
end
fprintf(rev_str);
end
% M4S DSA scheme
% ------------------------------------------------------------------------------
if strcmpi(data.Neutronics.Transport.DSAType,'m4s')
% Allocate Memory Space
SI_err_L2   = cell(sn_num, mfp_tot);
SI_err_inf  = cell(sn_num, mfp_tot);
SI_norm_L2  = cell(sn_num, mfp_tot);
SI_norm_inf = cell(sn_num, mfp_tot);
SI_iters = zeros(sn_num, mfp_tot);
NSR_err  = zeros(sn_num, mfp_tot);
NSR_norm = zeros(sn_num, mfp_tot);
% Loop through finite element orders
for f=1:length(fdeg)
    disp(['-> Finite Element Degree: ',num2str(f),' of ', num2str(length(fdeg))])
    for b=1:length(bf_name)
        disp(['  -> Basis Function: ',num2str(b),' of ', num2str(length(bf_name))])
        data = load_user_input(dim, bf_name{b}, fdeg(f), bc_type);
        data.Neutronics.Transport.DSAType = diff_type;
        % Loop through quadrature
        rev_str = [];
        for m=1:sn_num
            disp(['    -> Quadrature Set: ',num2str(m),' of ', num2str(sn_num)])
            data = load_quad_input( data, q_type, sn_levels(m) );
            dx_num = dx_num_start;
            geom = load_geometry_input(dim, m_type, dx_start, ar, 1);
            tc = 0;
            mfp = zeros(mfp_tot, 1);
            % First loop through upper mfp values in reverse order
            for j=mfp_upper:-1:1
                tc = tc + 1;
                msg = sprintf('        -> MFP Number: %d of %d',tc,mfp_tot);
                fprintf([rev_str,msg]);
                rev_str = repmat(sprintf('\b'), 1, length(msg));
                
                jj = mfp_lower + j;
                mfp(jj) = mfp_vals(j);
                txs = mfp_vals(j) / t_max_dim;
                data = load_xs_input( data, txs, c);
                [data, geom] = process_input_data(data, geom);
                data = cleanup_neutronics_input_data(data, geom);
                % Execute problem
                [data, sol, ~, ~, ~] = execute_problem(data, geom);
                % Collect statistics
                SI_iters(m,jj) = sol.iter;
                SI_err_L2{m,jj} = sol.error_L2;
                SI_err_inf{m,jj} = sol.error_inf;
                SI_norm_L2{m,jj} = sol.norm_L2;
                SI_norm_inf{m,jj} = sol.norm_inf;
            end
            % Next loop through mesh refinement steps
            jj = 1;
            for j=mfp_lower:-1:1
                jj = jj + 1;
                tc = tc + 1;
                msg = sprintf('      -> MFP Number: %d of %d',tc,mfp_tot);
                fprintf([rev_str,msg]);
                rev_str = repmat(sprintf('\b'), 1, length(msg));
                
                dx_num = (dx_num-1)*2+1;
                dx = linspace(0,L,dx_num);
                geom = load_geometry_input(dim, m_type, dx, ar, jj);
                mfp(j) = txs * (max(geom.CellVolume))^(1/dim);
                [data, geom] = process_input_data(data, geom);
                data = cleanup_neutronics_input_data(data, geom);
                % Execute problem
                glob.print_info = true;
                [data, sol, ~, ~, ~] = execute_problem(data, geom);
                glob.print_info = false;
                % Collect statistics
                SI_iters(m,j) = sol.iter;
                SI_err_L2{m,j} = sol.error_L2;
                SI_err_inf{m,j} = sol.error_inf;
                SI_norm_L2{m,j} = sol.norm_L2;
                SI_norm_inf{m,j} = sol.norm_inf;
            end
            % Process statistics
            for j=1:mfp_tot
                it = SI_iters(m,j); counter = 0;
                terr_L2 = 0; t_norm_L2 = 0;
                if it < 4, continue; end
                NSR_err(m,j) = SI_err_L2{m,j}(end)/SI_err_L2{m,j}(end-1);
                NSR_norm(m,j) = SI_norm_L2{m,j}(end)/SI_norm_L2{m,j}(end-1);
            end
            % Save off information
            if ~isequal(exist(dname, 'dir'),7), mkdir(dname); end
            q = sn_levels(m);
            if dim == 1
                fname = sprintf('%s%d_%s%d',bf_name{b},fdeg(f),q_type,q);
            else
                fname = sprintf('%s%d_%s%d_AR=%d',bf_name{b},fdeg(f),q_type,q,ar);
            end
            mat_out = [mfp, squeeze(NSR_norm(m,:))', squeeze(NSR_err(m,:))'];
            dlmwrite([dname,fname,'.dat'],mat_out,'precision','%14.8e');
        end
    end
end
fprintf(rev_str);
end