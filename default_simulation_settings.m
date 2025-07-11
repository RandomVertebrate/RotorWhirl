function simulation_settings = default_simulation_settings()

simulation_settings.num_processor_cores = 4;
simulation_settings.num_assumed_modes = 10;
simulation_settings.campbell_resolution = 10;
simulation_settings.modeshape_plot_rate = 5;
simulation_settings.max_nodal_locations = 4;
simulation_settings.max_num_modeshapes = 9;
simulation_settings.max_speed = 1000;
simulation_settings.air_gap = 20e-6;

end