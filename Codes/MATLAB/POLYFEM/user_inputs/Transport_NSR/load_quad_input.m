function data = load_quad_input( data, q_type, sn_levels )
data.Neutronics.Transport.PnOrder = 0;
data.Neutronics.Transport.QuadType = q_type;
data.Neutronics.Transport.SnLevels = sn_levels;
data.Neutronics.Transport.AngleAggregation = 'single';