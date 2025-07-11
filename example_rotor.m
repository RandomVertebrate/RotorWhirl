function rotor_data = example_rotor()

rotor_data.length = 1;
rotor_data.shaft_initial_radius = 5e-3;
rotor_data.lateral_support_locations = [0.3 0.7];
rotor_data.lateral_support_stiffnesses = [5e6 5e6];
rotor_data.axial_load_location = 0.7;
rotor_data.shaft_step_locations = [0.3 0.35 0.4 0.45];
rotor_data.shaft_step_values = [20 30 -30 -20]*1e-3;
rotor_data.axial_hole_radius = 2;
rotor_data.axial_hole_start = 500;
rotor_data.axial_hole_end = 500;
rotor_data.lumped_mass_location = [];
rotor_data.lumped_mass_value = [];

end