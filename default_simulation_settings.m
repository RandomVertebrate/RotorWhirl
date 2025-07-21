function simulation_settings = default_simulation_settings()

simulation_settings.num_processor_cores = 8;
simulation_settings.num_assumed_modes = 15;
simulation_settings.campbell_resolution = 10;
simulation_settings.modeshape_plot_rate = 5;
simulation_settings.max_nodal_locations = 4;
simulation_settings.max_num_modeshapes = 6;
simulation_settings.max_speed = 600;
simulation_settings.air_gap = (10:5:60)*1e-6;

end