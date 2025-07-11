function bearing_data = example_bearing()

% Updated bearing parameters
bearing_data.C1 = 230.571424;        % Curve-fitted bearing constant
bearing_data.C2 = 14112611380.972;   % Curve-fitted bearing constants
bearing_data.C3 = 2.035349;          % Curve-fitted bearing constants

bearing_data.disc_radius = 50e-3;
bearing_data.disc_thickness = 11e-3;
bearing_data.num_orifices = 6;
bearing_data.orifice_radial_location = 36e-3;
bearing_data.bearing_angle_offset = 0;

end