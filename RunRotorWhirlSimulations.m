function RunRotorWhirlSimulations(rotor_data,bearing_data,material_properties,simulation_settings)

    outputfile = fopen('output.txt', 'w');
    critfile = fopen('critspeeds.txt', 'w');
    
    
    processorcores = simulation_settings.num_processor_cores;
    
    modes = simulation_settings.num_assumed_modes;         % Number of assumed modes (higher = more accurate and slower)
    plotres = simulation_settings.campbell_resolution;       % Resolution of Campbell diagram (increase if diagram looks garbled)
    modespeeds = simulation_settings.modeshape_plot_rate;    % Number of speeds to draw modeshapes at
    maxnodes = simulation_settings.max_nodal_locations;
    max_num_modeshapes = simulation_settings.max_num_modeshapes;
    
    % All physical quantities in SI units
    
    maxspeed = simulation_settings.max_speed;   % Maximum spin speed (Campbell diagram horizontal axis upper limit)
    
    fprintf('Defining system parameters...\n')
    
    l = rotor_data.length;                   % Length of rotor
    sploc = rotor_data.lateral_support_locations;                     % x-coordinates of lateral support bearings

    fp = rotor_data.axial_load_location;                    % Point (x-coordinate) of application of axial load
    
    lmploc = rotor_data.lumped_mass_location;                    % lumped-mass locations

    lmpval = rotor_data.lumped_mass_value;                   % lumped-mass values
    num_lmp = length(lmploc);

    
    % Material Properties
    materialDensity = material_properties.density;
    youngsModulus = material_properties.youngs_modulus;
    shearModulus = material_properties.shear_modulus;

    % DEFINING DISC (gas bearing disc, attached to left end of shaft)
    DiscRadius = bearing_data.disc_radius;
    DiscThickness = bearing_data.disc_thickness;
    DiscMass = materialDensity*pi*DiscRadius^2*DiscThickness;
    DiscJp = 0.5*DiscMass*DiscRadius^2;
    DiscJt = 0.5*DiscJp+DiscMass*DiscThickness^2/12;
    DiscLoc = DiscThickness/2;
    
    % Axial hole details for hollow shaft
    holerad = rotor_data.axial_hole_radius;
    holestart = rotor_data.axial_hole_ends(1);  % x-coordinate of start of hole
    holestop = rotor_data.axial_hole_ends(2);   % x-coordinate of end of hole
    % (holestart == holestop) implies no hole
    
    % DEFINING SHAFT GEOMETRY
    if length(rotor_data.shaft_step_locations)>0
        baseradius = rotor_data.shaft_initial_radius;                   % Radius of shaft at junction with disc
        steploc = rotor_data.shaft_step_locations;                 % Axial locations (distances from shaft-disc junction) of steps in shaft profile
        stepval = rotor_data.shaft_step_values;                   % Values of radius increments at step locations
        % Rotor geometry
        syms xcoord real    
        ORexp = TrueShaftProfile(xcoord, baseradius, steploc, stepval);
        OR2exp = StiffShaftProfile(xcoord, baseradius, steploc, stepval, l, participationslope);
    
        % Radius function used for linear mass density, rotational inertia, and shear stiffness
        outerradius = @(x) subs(ORexp, xcoord, x);
    
        % Radius function used for bending stiffness
        outerradius2 = @(x) subs(OR2exp, xcoord, x);
    
        % Hole radius function (ignore for solid shaft)
        innerradius = @(x) holerad*heaviside(x-holestart) - holerad*heaviside(x-holestop);

    else
        % get radius functions from rotor_data
        outerradius = rotor_data.outer_radius;
        outerradius2 = rotor_data.stiff_radius;
        innerradius = @(x) rotor_data.inner_radius(x).*(heaviside(holestop-x)-heaviside(holestart-x));
    end
    
    participationslope = 1;
    
    % Gas bearing parameters
    C1 = bearing_data.C1;        % Curve-fitted bearing constant
    C2 = bearing_data.C2; % Curve-fitted bearing constants
    C3 = bearing_data.C3;          % Curve-fitted bearing constants
    
    nholes = bearing_data.num_orifices;             % number of holes
    rholes = bearing_data.orifice_radial_location;       % radial location of holes
    phi0 = bearing_data.bearing_angle_offset;               % bearing offset angle
    
    % AIR GAP VALUES TO RUN SIMULATIONS FOR
    gap = simulation_settings.air_gap;            % airgapvalues
    numloadvals = length(gap);
    
    % BEARING STIFFNESSES TO RUN SIMULATIONS FOR (kb1 in column 1, kb2 in column 2)
    kblist = rotor_data.lateral_support_stiffnesses;       %in N/m
    kbmatsize = size(kblist);
    num_kbs = kbmatsize(1);
    
    firstiteration = true;
    
    fprintf(critfile, "Axial Load (N)\tClearance (m)\tTotal Axial Stiffness (N/m)\tLateral Stiffness (N/m)\tCritical Speeds (rad/s)\n");
    
    for kbiter = 1:num_kbs
        kbvals = kblist(kbiter,:);
        for gapiter = 1:numloadvals
            h0=gap(gapiter);
            RotorWhirl_RigidDisc;
        
            critformat = '%.4e';
            numcrits = length(criticalspeeds);
            for i=1:numcrits-1
                critformat = strcat(critformat,'\t%.4e');
            end
            critformat = strcat(critformat,'\n');
            
            mainformat = '';
            for i=1:plotres
                mainformat = strcat(mainformat,'\t%f');
            end
            mainformat = strcat(mainformat,'\n');
            
            fprintf(outputfile, '\nBearing stiffnesses %.e,', kbvals(1));
            for i=2:num_bearings
                fprintf(outputfile, '%.e, ', kbvals(i));
            end
            fprintf(outputfile, ...
                'Load = %f N, Clearance = %.4e m, Total Axial Stiffness = %.4e N/m\n', ...
                Fax, h0, TotalStiffness);
            fprintf(outputfile, 'Spin Speeds (rad/s)');
            fprintf(outputfile, mainformat, speeds);
            
            fprintf(outputfile, 'Natural Frequencies (Hz)');
            for counter = 1:modes
                fprintf(outputfile, mainformat, points(counter,:));
            end
            
            fprintf(outputfile, "Critical speeds in rad/s:\t"+critformat, criticalspeeds(:));
    
            fprintf(critfile, "%.4e\t%.4e\t%.4e\t", [Fax, h0, TotalStiffness, kbvals(1)]);
            fprintf(critfile, critformat, criticalspeeds(:));
            
            saveas(gcf, 'Campbell/Frames/kbval'+string(kbiter)+'Frame'+string(gapiter)+'.png');
            
            if firstiteration
                firstiteration = false;
            end
        end
    end
    
    fprintf(outputfile, '\n\nRotor Length = %f m\n', l);
    fprintf(outputfile, 'Disc radius = %f m\n', DiscRadius);
    fprintf(outputfile, 'Disc thickness = %f m\n', DiscThickness);
    fprintf(outputfile, 'True shaft profile = %s\n', string(outerradius(x)));
    fprintf(outputfile, 'Effective rigidity profile = %s\n', string(outerradius2(x)));
    fprintf(outputfile, 'hole profile = %s from %s to %s\n', string(innerradius(x)), string(holestart), string(holestop));
    fprintf(outputfile, 'Support locations from left side (gas bearing side)');
    for i=1:num_bearings
        fprintf(outputfile, ', %f', sploc(i))
    end
    fprintf(outputfile, '\nBearing constants: C1 = %f N, C2 = %f m^(-C3), C3 = %f\n', C1, C2, C3);
    fprintf(outputfile, 'Number of bearing holes = %d\n', nholes);
    fprintf(outputfile, 'Radial location of bearing holes = %f m\n', rholes);
    
    fclose(outputfile)
    fclose(critfile)
    %fclose(Nodal_loc)
end