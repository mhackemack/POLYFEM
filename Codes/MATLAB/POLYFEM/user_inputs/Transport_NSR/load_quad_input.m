function data = load_quad_input( data, q_type, sn_levels )
data.Neutronics.Transport.fluxMoments = 0;
data.Neutronics.Transport.QuadType = q_type;
data.Neutronics.Transport.SnLevels = sn_levels;