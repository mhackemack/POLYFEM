function out = get_globals( s )
% Numbers
% ------------------------------------------------------------------------------
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
% ------------------------------------------------------------------------------
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
% ------------------------------------------------------------------------------
out.input_path = 'inputs/';
out.output_path = 'outputs/';
out.XS_path = 'XSFiles/';
out.print_info = true;

