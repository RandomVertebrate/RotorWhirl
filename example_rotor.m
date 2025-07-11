function rotor_data = example_rotor()

rotor_data.length = 1;
rotor_data.shaft_initial_radius = 5e-3;
rotor_data.lateral_support_locations = [0.15 0.6 0.8];
rotor_data.lateral_support_stiffnesses = [5e6 5e6 5e6];
rotor_data.axial_load_location = 0.8;
rotor_data.shaft_step_locations = [0.3 0.35 0.4 0.45];
rotor_data.shaft_step_values = [20 130 -130 -20]*1e-3;
rotor_data.axial_hole_radius = 18e-3;
rotor_data.axial_hole_ends = [0.325 0.425];
rotor_data.lumped_mass_location = [];
rotor_data.lumped_mass_value = [];

end