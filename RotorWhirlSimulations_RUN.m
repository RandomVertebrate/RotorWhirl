clear all
outputfile = fopen('output.txt', 'w');
critfile = fopen('critspeeds.txt', 'w');

processorcores = 8;

modes = 7;          % Number of assumed modes (higher = more accurate and slower)
plotres = 8;        % Resolution of Campbell diagram (increase if diagram looks garbled)
modespeeds = 3;     % Number of speeds to draw modeshapes at

% All physical quantities in SI units

maxspeed = 3000;   % Maximum spin speed (Campbell diagram horizontal axis upper limit)

fprintf('Defining system parameters...\n')

l = 657.5e-3;                   % Length of rotor
sploc = ...                     % x-coordinates of lateral support bearings
     [50, 330, 587]*1e-3;
fp = 516e-3;                    % Point (x-coordinate) of application of axial load

lmploc = ...                    % lumped-mass locations
    [300 657.5]*1e-3;
lmpval = ...                    % lumped-mass values
    [0.5 0.15];
num_lmp = length(lmploc);

% DEFINING DISC (gas bearing disc, attached to left end of shaft)
DiscRadius = 60e-3;
DiscThickness = 15e-3;
DiscMass = 7850*pi*DiscRadius^2*DiscThickness;
DiscJp = 0.5*DiscMass*DiscRadius^2;
DiscJt = 0.5*DiscJp+DiscMass*DiscThickness^2/12;
DiscLoc = DiscThickness/2;

% Axial hole details for hollow shaft
holerad = 6e-3;
holestart = 25e-3;  % x-coordinate of start of hole
holestop = 175e-3;   % x-coordinate of end of hole
% (holestart == holestop) implies no hole

% DEFINING SHAFT GEOMETRY
baseradius = 12.5e-3;                   % Radius of shaft at junction with disc
steploc = [100 150 350 400 500 520]*1e-3;                 % Axial locations (distances from shaft-disc junction) of steps in shaft profile
stepval = [40 -30 -10 20 50 -60]*1e-3;                   % Values of radius increments at step locations

participationslope = 1;

% Gas bearing parameters
% PCD75
C1 = 290.120806;        % Curve-fitted bearing constant
C2 = 9024693618.143471; % Curve-fitted bearing constants
C3 = 2.004049;          % Curve-fitted bearing constants

% %PCD25
% C1 = 75.520732;        % Curve-fitted bearing constant
% C2 = 9110281087.324850; % Curve-fitted bearing constants
% C3 = 1.998757;          % Curve-fitted bearing constants

nholes = 6;             % number of holes
rholes = 50e-3;         % radial location of holes
phi0 = 0;               % bearing offset angle
mingap = 10e-6;         % minimum bearing gap
maxgap = 100e-6;        % maximum bearing gap

% AIR GAP VALUES TO RUN SIMULATIONS FOR
Airgapvalues = [10:2:70];            % airgapvalues
%Airgapvalues = [70];                % airgapvalues
gap = Airgapvalues*1e-6;
numloadvals = length(gap);

% BEARING STIFFNESSES TO RUN SIMULATIONS FOR (kb1 in column 1, kb2 in column 2)
kblist = [1e5 1e5 1e5];       %in N/m
kbmatsize = size(kblist);
num_kbs = kbmatsize(1);

% Mode shapes plotting
maxnodes = 4
max_num_modeshapes = 6

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
fprintf(outputfile, 'True shaft profile = %s\n', string(ORexp));
fprintf(outputfile, 'Effective rigidity profile = %s\n', string(ORexp));
fprintf(outputfile, 'Support locations from left side (gas bearing side)');
for i=1:num_bearings
    fprintf(outputfile, ', %f', sploc(i))
end
fprintf(outputfile, '\nBearing constants: C1 = %f N, C2 = %f m^(-C3), C3 = %f\n', C1, C2, C3);
fprintf(outputfile, 'Number of bearing holes = %d\n', nholes);
fprintf(outputfile, 'Radial location of bearing holes = %f m\n', rholes);

fclose(outputfile)
fclose(critfile)