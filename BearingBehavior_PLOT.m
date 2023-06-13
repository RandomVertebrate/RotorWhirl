%Aerostatic Thrust Bearing Description
%bearing pad with n orifices equally distributed along a circle of radius r.
%bearing gap h is measured at the center of the circle.
%bearing pad tilts by angle theta about an axis that:
%a) is in the plane of the bearing face
%b) passes through the center of the pad
%c) makes angle phi0 with a radial line extending to a hole.

clear all
syms x C1 C2 C3 n r phi0 k h theta hmicrons thetasecs

phi = 2*pi/n %angle between radial lines extending to successive orifices

%Nonlinear spring behaviour of each orifice (force as a function of gap)
f(x) = C1/(1+C2*x^(C3))

%Total thrust force in (or load capacity of) bearing
F(h,theta, C1, C2, C3) = symsum(f(h+theta*r*sin(k*phi+phi0)),k,0,n-1)

%Resultant moment on bearing pad about tilt axis due to orifice forces
M(h, theta, C1, C2, C3) = symsum(r*sin(k*phi+phi0)*f(h+theta*r*sin(k*phi+phi0)),k,0,n-1)

%constants defining nonlinear spring behavior W=C1/(1+C2*h^C3) for each orifice
%obtained by minimizing square of load capacity error relative to data from CFD or experiment
C1val = 956.428643
C2val = 303648875.241726
C3val = 1.669899

%Parameter value substitutions for plots
nval = 6
rval = 37.5e-3
phi0val = 0
params = [n r phi0]
values = [nval rval phi0val]

%Force as a function of h in microns and theta in seconds
%(instead of meters and radians)
Feval(hmicrons,thetasecs) = subs(F(hmicrons*1e-6, thetasecs*pi/(180*3600), C1val, C2val, C3val), params, values)

%Moment as a function of h in microns and theta in seconds
%(instead of meters and radians)
Meval(hmicrons,thetasecs) = subs(M(hmicrons*1e-6, thetasecs*pi/(180*3600), C1val, C2val, C3val), params, values)

%Outer radius of disc, used to calculate maximum tilt given h
%or minimum h given tilt
Ro = 50e-3
%Alternate expression for Ro: makes minimum orifice gap zero
%Ro = rval*max(cos(phi0val),cos(2*pi/nval-phi0val))

hmicronsChoice = [10 15 20 25 30 35 40] %Select operating clearances for plotting bearing behavior
tiltsecsChoice = [0 8.2506 12.376 16.501 24.752 37.128 49.504 74.255 ] %Select tilt values for plotting bearing behavior
hmicronsmax = 40 %Upper limit of clearance for plots

numhChoices = length(hmicronsChoice)
numtiltChoices = length(tiltsecsChoice)

%Maximum tilt angles (in seconds) corresponding to chosen operating clearances
tiltsecsmax = (hmicronsChoice*1e-6/Ro)*180*3600/pi

%Minimum clearances corresponding to chosen tilt angles
hmicronsmin = (tiltsecsChoice*(pi/(180*3600)))*Ro*1e6

heights = zeros(numhChoices, 1000);
tilts = zeros(numtiltChoices, 1000);

%Generating arrays of heights and tilt angles for horizontal axis data of plots
heights0 = linspace(0, hmicronsmax(1), 1000);
for i = 1:numtiltChoices
    heights(i,:) = linspace(hmicronsmin(i), hmicronsmax, 1000);
end
for i = 1:numhChoices
    tilts(i,:) = linspace(-tiltsecsmax(i), tiltsecsmax(i), 1000);
end

LW = 1 %line weight for plots
%plot styles
linestyles = ["-k" "-g" "--k" "--g" "-.k" "-.g" ":k" ":g"]
pointstyles = ["xg" "xk" "+g" "+k" "^g" "^k" "squareg" "squarek"]

%Read CFD/Experiment data from excel file (for comparison)
DataTable = readtable('FitData.xlsx')
NumDataPts = height(DataTable)
data.height = DataTable.h_microns_
data.tilt = DataTable.Tilt_seconds_
data.load = DataTable.W_kN_*1e3 %convert kN to N
data.moment = DataTable.Moment_Nm_

%"relevant data" is data at the chosen operating clearances or tilt values
%this data can be plotted alongside fitted curves to show the quality of the fit
%gaprelevantdata is data at a relevant gap value (chosen operating clearance)
%tiltrelevantdata is data at a relevant tilt value (chosen operating tilt)
%now defining data structures to hold "relevant" subsets of available data
for i=1:numhChoices
    %gaprelevantdata(i) will hold data at ith chosen operating clearance
    %i.e. at gap = hmicronsChoice(i) microns
    gaprelevantdata(i).height = [];
    gaprelevantdata(i).tilt = [];
    gaprelevantdata(i).load = [];
    gaprelevantdata(i).moment = [];
end
for i=i:numtiltChoices
    %tiltrelevantdata(i) will hold data at ith chosen operating tilt
    %i.e. at tilt = tiltsecsChoice(i) seconds
    tiltrelevantdata(i).height = [];
    tiltrelevantdata(i).tilt = [];
    tiltrelevantdata(i).load = [];
    tiltrelevantdata(i).moment = [];
end

tolerance = 0.05; %tolerance for datapoint being "at" the operating point

%now populating "relevant" data sets
for i=1:numhChoices %iterating through operating points (and relevant data sets)
   j = 0;
   for k=1:NumDataPts %iterating through datapoints
       %if the current data point is at the current operating clearance
       if abs(data.height(k)-hmicronsChoice(i)) < tolerance
           %add current datapoint to current relevant dataset
           j = j + 1;
           gaprelevantdata(i).height(j) = data.height(k);
           gaprelevantdata(i).tilt(j) = data.tilt(k);
           gaprelevantdata(i).load(j) = data.load(k);
           gaprelevantdata(i).moment(j) = data.moment(k);
       end
   end
end

for i=1:numtiltChoices
   j = 0;
   for k=1:NumDataPts %iterating through datapoints
        %if the current data point is at the current operating tilt
        if abs(-data.tilt(k)-tiltsecsChoice(i)) < tolerance
            %add current datapoint to current relevant dataset
            j = j + 1;
            tiltrelevantdata(i).height(j) = data.height(k);
            tiltrelevantdata(i).tilt(j) = data.tilt(k);
            tiltrelevantdata(i).load(j) = data.load(k);
            tiltrelevantdata(i).moment(j) = data.moment(k);
        end
    end
end

%additional relevant data set: zero tilt condition
j = 0;
for k=1:NumDataPts
    if abs(data.tilt(k)) < tolerance
        j = j + 1;
        tiltrelevantdata0tilt.height(j) = data.height(k);
        tiltrelevantdata0tilt.load(j) = data.load(k);
        tiltrelevantdata0tilt.moment(j) = data.moment(k);
    end
end

%display relevant data
for i=1:numhChoices
    disp(strcat('Data at ',num2str(hmicronsChoice(i)),' microns'))
    [gaprelevantdata(i).height(:) gaprelevantdata(i).tilt(:) ...
        gaprelevantdata(i).load(:) gaprelevantdata(i).moment(:)]
end
for i=1:numtiltChoices
    disp(strcat('Data at ',num2str(tiltsecsChoice(i)),' seconds of tilt'))
    [tiltrelevantdata(i).height(:) tiltrelevantdata(i).tilt(:) ...
        tiltrelevantdata(i).load(:) tiltrelevantdata(i).moment(:)]
end

%Plotting gap vs load
subplot(2, 2, 1)
hold on
for i=1:numtiltChoices
    plot(heights(i,:), Feval(heights(i,:), tiltsecsChoice(i)), ...
        linestyles(i), 'LineWidth', LW)
end
for i=1:numtiltChoices
    plot(tiltrelevantdata(i).height, tiltrelevantdata(i).load, pointstyles(i))
end
grid minor
xlabel('Gap (\mu{}m)')
ylabel('Load (N)')
legendentries1 = strings(0);
for i=1:numtiltChoices
    legendentries1(end+1)=strcat('Curve fit, tilt = ', num2str(tiltsecsChoice(i)), ' ''''');
end
for i=1:numtiltChoices
    legendentries1(end+1)=strcat('CFD, tilt = ', num2str(tiltsecsChoice(i)), ' ''''');
end
legend(legendentries1)
hold off

%Plotting tilt vs load
subplot(2, 2, 2)
hold on
for i=1:numhChoices
    plot(tilts(i,:), Feval(hmicronsChoice(i), tilts(i,:)), linestyles(i), 'LineWidth', LW)
end
for i=1:numhChoices
    plot(gaprelevantdata(i).tilt, gaprelevantdata(i).load, pointstyles(i))
end
grid minor
xlabel('Tilt ('''')')
ylabel('Load (N)')
legendentries2 = strings(0);
for i=1:numhChoices
    legendentries2(end+1)=strcat('Curve fit, gap = ', num2str(hmicronsChoice(i)), ' \mu{}m');
end
for i=1:numhChoices
    legendentries2(end+1)=strcat('CFD, gap = ', num2str(hmicronsChoice(i)), ' \mu{}m');
end
legend(legendentries2)
hold off

%Plotting height vs moment
subplot(2, 2, 3)
hold on
for i=1:numtiltChoices
    plot(heights(i,:), -Meval(heights(i,:), tiltsecsChoice(i)), linestyles(i), 'LineWidth', LW)
end
for i=1:numtiltChoices
    plot(tiltrelevantdata(i).height, tiltrelevantdata(i).moment, pointstyles(i))
end
grid minor
xlabel('Gap (\mu{}m)')
ylabel('Moment (Nm)')
legendentries3 = strings(0);
for i=1:numtiltChoices
    legendentries3(end+1)=strcat('Curve fit, tilt = ', num2str(tiltsecsChoice(i)), ' ''''');
end
for i=1:numtiltChoices
    legendentries3(end+1)=strcat('CFD, tilt = ', num2str(tiltsecsChoice(i)), ' ''''');
end
legend(legendentries3)
hold off

%Plotting tilt vs moment
subplot(2, 2, 4)
hold on
for i=1:numhChoices
    plot(tilts(i,:), Meval(hmicronsChoice(i), tilts(i,:)), linestyles(i), 'LineWidth', LW)
end
for i=1:numhChoices
    plot(gaprelevantdata(i).tilt, gaprelevantdata(i).moment, pointstyles(i))
end
grid minor
xlabel('Tilt ('''')')
ylabel('Moment (Nm)')
legendentries4 = strings(0);
for i=1:numhChoices
    legendentries4(end+1)=strcat('Curve fit, gap = ', num2str(hmicronsChoice(i)), ' \mu{}m');
end
for i=1:numhChoices
    legendentries4(end+1)=strcat('CFD, gap = ', num2str(hmicronsChoice(i)), ' \mu{}m');
end
legend(legendentries4)
hold off

