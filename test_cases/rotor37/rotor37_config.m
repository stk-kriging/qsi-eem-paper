% ROTOR37_CONFIG EEM-SMC configuration of the ROTOR37 case

function config = rotor37_config ()

prm = rotor37_struct ();

config = eem_smc_config_base (prm);

% smaller initial design than the default 10 (dim_x + dim_s)
config.pts_init = 5 * (prm.dim_x + prm.dim_s);

config.stop = 0.2;
config.stop_final = 0.2;

end % function
