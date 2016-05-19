function out = get_globals( s )
% Numbers
%%%%%%%%%
out.small = 1e-14;
out.large = 1/eps/10;
if strcmp(s, 'Home')
    % Home - more memory
    out.maxMatrix = 2.5e4;
    out.maxSparse = 1.1e5;
else
    % Office
    out.maxMatrix = 1e4;
    out.maxSparse = 6e4;
end
% Constants
%%%%%%%%%%%
out.speed_of_light = 2.99792458e10;
out.planck = 6.62606957e-34;
out.boltzmann_ev = 8.6173324e-5;
out.boltzmann_J = 1.3806488e-23;
out.boltzmann_erg = 1.3806488e-16;
out.elementary_charge = 1.6021892e-19;
out.neutron_mass = 1.6749544e-27;
out.avogadro = 6.022141e23;
out.gas_constant = 8.3144621;

out.amu_to_kg = 1.6605655e-27;
out.amu_to_MeV = 931.5016;
out.MeV_to_J = 1.601892e-13;
out.J_to_MeV = 6.2415065e12;
% Miscellaneous
%%%%%%%%%%%%%%%
out.input_path = 'user_inputs/';
out.output_path = '../Outputs/';
out.geom_path = 'geometry_inputs/precompiled/';
out.geom_raw_tri_tet_path = 'geometry_inputs/raw_tri_tet_files/';
out.geom_raw_dg_path = 'geometry_inputs/raw_dg_files/';
out.xs_dir = '../XSFiles/';
out.directory = [];
out.print_info = true;
out.plot_counter = 0;
% Global Face Condition Types
% ---------------------------
out.Periodic = -1;
out.Interior = 0;
out.Function = 9;
% Boundary Conditions
% -------------------
out.Dirichlet = 1;
out.Neumann = 2;
out.Robin = 3;
% Transport Boundary Conditions
% -----------------------------
out.Vacuum = 1;
out.Reflecting = 2;
out.IncidentIsotropic = 3;
out.IncidentCurrent = 4;
out.IncidentBeam = 5;
% Acceleration Types
% ------------------
out.Accel_WGS_DSA = 1;
out.Accel_WGS_TSA = 2;
out.Accel_AGS_TG  = 3;
out.Accel_AGS_MTG  = 4;
out.Accel_AGS_TTG = 5;
out.Accel_AGS_MTTG = 6;
out.Accel_WGS_MJIA_DSA = 7;
out.Accel_Fission_DSA = 8;
% Acceleration Discretizations
% ----------------------------
out.Accel_DSA_MIP = 1;
out.Accel_DSA_IP  = 2;
out.Accel_DSA_DCF = 3;
out.Accel_DSA_M4S = 4;
out.Accel_TSA_DFEM = 5;
out.Accel_TSA_CFEM = 6;
