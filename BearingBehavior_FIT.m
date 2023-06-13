%Aerostatic Thrust Bearing Description
%bearing pad with n orifices equally distributed along a circle of radius r.
%bearing gap h is measured at the center of the circle.
%bearing pad tilts by angle theta about an axis that:
%a) is in the plane of the bearing face
%b) passes through the center of the pad
%c) makes angle phi0 with a radial line extending to a hole.

%The mathematical model:
%Each orifice behaves like a nonlinear spring whose force is given by:
%  C1/(1+C2*h^C3)
%where C1 C2 and C3 are constants and where h is the bearing gap at the orifice location.

%This program:
%reads CFD/Experimental data for net bearing load at different tilt and height conditions
%and fits C1 C2 and C3 to this data.

clear all

syms C1 C2 C3 theta h

nholes = 6 %number of evenly distributed holes
r = 0.0375 %radial distance of holes from center of disc
phi0 = 0 %angular coordinate of first hole (tilt axis is at zero)

%Read CFD/Experiment data from excel file
%Excel file has the following columns:
%  h_microns_  Tilt_seconds_  W_kN_  Moment_Nm_
DataTable = readtable('FitData.xlsx')
NumDataPts = height(DataTable)
data.height = DataTable.h_microns_*1e-6; %convert microns to meters
data.tilt = DataTable.Tilt_seconds_*pi/(180*3600); %convert seconds to radians
data.load = DataTable.W_kN_*1e3; %convert kN to N
data.moment = DataTable.Moment_Nm_

%anonymous function on which fminsearch operates. Passes system parameters
%to ForceError(), which actually calculates error.
ErrorFunc = @(C) ForceError(nholes, r, phi0, data, NumDataPts, C(1), C(2), C(3));

%Find values of C1, C2, C3 that minimize error between data and calculated value.
fittolerance = 1e-8
options = optimset('MaxFunEvals', 1e5, 'MaxIter', 1e5, 'Display', 'iter', 'TolX', fittolerance, 'TolFun', fittolerance);
Cvals = fminsearch(ErrorFunc, [1000 1e10 2], options);
fprintf('C1 = %f, C2 = %f, C3 = %f \n', Cvals(1), Cvals(2), Cvals(3))

%calculates total squared error over dataset given C1, C2, C3
function Error = ForceError(nholes, r, phi0, data, NumDataPts, C1, C2, C3)
    
    Error = 0;
    %iterate through datapoints and add squares of differences between data and calculated value
    for i = 1:NumDataPts
        Force = TotalForce(C1, C2, C3, nholes, r, phi0, data.height(i), data.tilt(i));
        Error = Error + (Force - data.load(i))^2;
    end
end

%calculates force given C1 C2 C3 and height and tilt
function Force = TotalForce(C1, C2, C3, nholes, r, phi0, h, theta)
    
    phi = 2*pi/nholes;
    
    Force = 0;
    
    %iterate through holes and add forces
    for i = 0:(nholes-1)
        hi = h + theta*r*sin(phi0+i*phi); %height of ith hole
        Force = Force + C1/(1+C2*hi^C3);
    end
end