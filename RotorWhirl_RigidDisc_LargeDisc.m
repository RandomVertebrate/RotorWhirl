% The following variables must be defined before running this script
% kb1 (Lateral stiffness of bearing 1)
% kb2 (Lateral stiffness of bearing 2)
% FaxPercent (Percentage of maximum axial load applied)
% firstiteration (logical: some computer algebra skipped for subsequent iterations)

processorcores = 48;
delete(gcp('nocreate'));
parpool(processorcores);

if firstiteration
    modes = 15;        % Number of assumed modes (higher = more accurate and slower)
    plotres = 7;      % Resolution of Campbell diagram (increase if diagram looks garbled)

    % All physical quantities in SI units

    maxspeed = 18000;  % Maximum spin speed (Campbell diagram horizontal axis upper limit)

    fprintf('Defining system parameters...\n')

    l = 136.6e-3;           % Length of rotor
    sp1 = 57e-3;            % Distance of first support from left side
    sp2 = 77e-3;            % Distance of second support from left side      
    fp = l;                 % Point (x-coordinate) of application of axial load
    
    % DEFINING DISC (gas bearing disc, attached to left end of shaft)
    DiscRadius = 50.5e-3;
    DiscThickness = 11e-3;
    DiscMass = 7850*pi*DiscRadius^2*DiscThickness;
    DiscJp = 0.5*DiscMass*DiscRadius^2;
    DiscJt = 0.5*DiscJp+DiscMass*DiscThickness^2/12;
    DiscLoc = DiscThickness/2;

    % Axial hole details for hollow shaft
    holerad = 3e-3;
    holestart = 127e-3;  % x-coordinate of start of hole
    holestop = 127e-3;   % x-coordinate of end of hole
    % (holestart == holestop) implies no hole
    
    % DEFINING SHAFT GEOMETRY
    baseradius = 20e-3;                             % Radius of shaft at junction with disc
    steploc = [4 24 91 125.6]*1e-3;                 % Axial locations (distances from shaft-disc junction) of steps in shaft profile
    stepval = [-12 7 -3 -8]*1e-3;                   % Values of radius increments at step locations
    numsteps = length(steploc);
    
    participationslope = 1;
    
    ORexp = sym(baseradius);
    OR2exp = sym(baseradius);
    
    syms xcoord
    
    for i = 1:numsteps
        ORexp = ORexp + stepval(i)*heaviside(xcoord-steploc(i));
    end
    
    for i = 1:numsteps
        if stepval(i)<0
            rampbeg = steploc(i) + stepval(i)/participationslope;
            rampend = steploc(i);
            OR2exp = OR2exp - (xcoord-rampbeg)*participationslope ...
                *heaviside(xcoord-rampbeg)*heaviside(rampend-xcoord) ...
                + stepval(i)*heaviside(xcoord-rampend);
        else
            rampbeg = steploc(i);
            rampend = steploc(i) + stepval(i)/participationslope;
            OR2exp = OR2exp + (xcoord-rampbeg)*participationslope ...
                *heaviside(xcoord-rampbeg)*heaviside(rampend-xcoord) ...
                + stepval(i)*heaviside(xcoord-rampend);
        end
    end

    % Radius function used for linear mass density, rotational inertia, and shear stiffness
    outerradius = @(x) subs(ORexp, xcoord, x);

    % Radius function used for bending stiffness
    outerradius2 = @(x) subs(OR2exp, xcoord, x);

    % Hole radius function (ignore for solid shaft)
    innerradius = @(x) holerad*heaviside(x-holestart) - holerad*heaviside(x-holestop);

    % Drawing rotor geometry...
    xvals = linspace(0, l, 1000);
    holexvals = linspace(0.999*holestart, 1.001*holestop, 1000);
    bearinglocs = [sp1, sp1, sp2, sp2];
    bearingrads = [outerradius(sp1), -outerradius(sp1), outerradius(sp2), -outerradius(sp2)];
    zeroradius = double(outerradius(0));
    zeroholeradius = double(innerradius(0));
    endradius = double(outerradius(l));
    endholeradius = double(innerradius(l));
    figure(1)
    clf
    plot(xvals, outerradius(xvals), 'k')
    hold on
    plot(xvals, outerradius2(xvals), '--k')
    plot(xvals, -outerradius(xvals), 'k')
    plot(xvals, -outerradius2(xvals), '--k')
    fill([-DiscThickness, 0, 0, -DiscThickness], [DiscRadius, DiscRadius, -DiscRadius, -DiscRadius], 'k')
    plot(holexvals, innerradius(holexvals), 'k')
    plot(holexvals, -innerradius(holexvals), 'k')
    plot([0 0], [-zeroradius -zeroholeradius], 'k', [0 0], [zeroradius zeroholeradius], 'k',...
        [l l], [-endradius -endholeradius], 'k', [l l], [endradius endholeradius], 'k')
    plot(bearinglocs, bearingrads, 'ob')
    forcearrow = quiver(fp+10e-3, 0, -10e-3, 0, 'r');
    forcearrow.LineWidth = 2;
    forcearrow.MaxHeadSize = 1;
    ylim([-DiscRadius*1.1, DiscRadius*1.1])
    xlim([-DiscThickness-0.1*(l+DiscThickness), 1.1*(l+DiscThickness)])
    title('Rotor Geometry (meters)')
    hold off
    axis equal
    grid minor
    shg
    % End drawing rotor geometry
end

% Gas bearing parameters
C1 = 925.084065;        % Curve-fitted bearing constant
C2 = 336266293.849243;  % Curve-fitted bearing constants
C3 = 1.683283;          % Curve-fitted bearing constants
nholes = 6;             % number of holes
rholes = 37.5e-3;       % radial location of holes
phi0 = 0;               % bearing offset angle
mingap = 0;             % minimum bearing gap
maxgap = 1e-3;          % maximum bearing gap

syms h theta k

% CALCULATING BEARING BEHAVIOR

% angle between bearing holes
phigap = 2*pi/nholes;
% Total bearing force (load capacity as a function of gap and tilt)
BearingForce = @(h, theta) symsum(C1/(1+C2*(h+rholes*theta*sin(phi0+k*phigap))^(C3)), k, 0, nholes-1);
% Total bearing moment
BearingMoment = @(h, theta) symsum(C1*rholes*sin(phi0+k*phigap)/(1+C2*(h+rholes*theta*sin(phi0+k*phigap))^(C3)), k, 0, nholes-1);

MaxLoad = BearingForce(0,0);

Fax = FaxPercent*MaxLoad/100;                           % Axial load

NetFlatForce = @(h) eval(Fax - BearingForce(h,0));      % Single-argument load function to solve for operating clearance

h0 = lsqnonlin(NetFlatForce, 1e-5, mingap, maxgap);     % Operating clearance for chosen load value

TotalStiffness = subs(diff(NetFlatForce(h), h), h, h0); % Total translational film stiffness at operating clearance

fprintf('\nCurrent load = %f N\nCurrent operating clearance = %f m\nCurrent total film stiffness = %f N/m\nForce residual = %f N\n\n', ...
    Fax, h0, TotalStiffness, NetFlatForce(h0))

syms epsilon x t

% Local linearization of bearing behavior
bfApprox = @(h,theta) subs(taylor(BearingForce(h0+epsilon*(h-h0), epsilon*theta), epsilon, 0, 'Order', 2), epsilon, 1);
bmApprox = @(h,theta) subs(taylor(BearingMoment(h0+epsilon*(h-h0), epsilon*theta), epsilon, 0, 'Order', 2), epsilon, 1);

if firstiteration

    % CALCULATING ROTOR PROPERTIES

    rhoA = @(x) expand(eval(7850*pi*(outerradius(x)^2-innerradius(x)^2)));           % Linear Mass Density

    EI = @(x) expand(eval(2e11*pi*(outerradius2(x)^4-innerradius(x)^4)/4));           % Flexural Rigidity

    KAG = @(x) expand(eval(0.9*pi*(outerradius(x)^2-innerradius(x)^2)*7.6923e10));   % Effective Shear Rigidity

    Jp = @(x) expand(eval((1/4)*rhoA(x)*(outerradius(x)^2+innerradius(x)^2)));       % Axial Mass Moment of Inertia per Unit Length

    syms thetaz(x,t) thetay(x,t) omegas

    % CALCULATING ROTATION MATRICES
    % Epsilon multiplied with small angles for linearization later on

    % First rotation matrix, thetaz about the z axis
    R1 = [cos(epsilon*thetaz), -sin(epsilon*thetaz), 0; ...
        sin(epsilon*thetaz), cos(epsilon*thetaz), 0; ...
        0, 0, 1];
    % Second rotation matrix, thetay about the y axis
    R2 = [cos(epsilon*thetay), 0, sin(epsilon*thetay); ...
        0, 1, 0; ...
        -sin(epsilon*thetay), 0, cos(epsilon*thetay)];
    % Third rotation matrix, omegas*t about the x axis
    R3 = [1, 0, 0; ...
        0, cos(omegas*t), -sin(omegas*t); ...
        0, sin(omegas*t), cos(omegas*t)];
    % Effective rotation matrix (rotations are body-fixed)
    R = R1*R2*R3;

    % FINDING OMEGA FROM R

    % First, finding [Omega cross], i.e. the skew symmetric matrix by which premultiplication is equivalent to a cross product with omega on the left side
    OmegaMatrixFunction = diff(R,t)*transpose(R);
    % Converting the matrix from a symbolic function to an expression so that it can be indexed
    OmegaMatrix = OmegaMatrixFunction(x,t);
    % Pulling out the components of Omega individually to create a column vector
    Omega = simplify(expand([OmegaMatrix(3,2); OmegaMatrix(1,3); OmegaMatrix(2,1)]));

    % CALCULATING J (LENGTH-NORMALIZED MOMENT OF INERTIA TENSOR)

    % Each rotor element is a disc, so transverse inertia is half of the axial inertia
    Jt = @(x) Jp(x)/2;
    % Inertia in the reference configuration
    J0 = @(x) [Jp(x) 0 0; 0 Jt(x) 0; 0 0 Jt(x)];
    % Transform by R to find the current inertia
    J = @(x) R*J0(x)*transpose(R);
    
    DiscJ0 = [DiscJp 0 0; 0 DiscJt 0; 0 0 DiscJt];
    DiscJ = R*DiscJ0*transpose(R);

    fprintf('Defining assumed modes...\n')

    % POLYNOMIAL ASSUMED MODE BASIS

    syms Qshapes

    parfor i = 1:2
        Qshapes(i)=x^(i-1);
    end
    
    parfor i = 3:modes
        Qshapes(i)=expand((x/l)^(i-1));
    end

    % DEFINING GENERALIZED COORDINATES
    A = {};
    for i = 1:4*modes+1
        A{i}(t) = str2sym('q'+string(i)+'(t)');
    end

    % Defining y- and z-direction deflection and y- and z-axis shear as sum of products of space and time-dependent parts
    syms ydef zdef phiy phiz
    ydef = 0; zdef = 0; phiy = 0; phiz = 0;
    for i=1:modes
        ydef = ydef + A{i}(t)*Qshapes(i);
        zdef = zdef + A{i+modes}(t)*Qshapes(i);
        phiy = phiy + A{i+2*modes}(t)*Qshapes(i);
        phiz = phiz + A{i+3*modes}(t)*Qshapes(i);
    end

    fprintf('Finding energy expressions...\n')

    % Rotational Kinetic Energy per Unit Length
    fprintf('Defining Rotational KE...\n')
    KEr = subs(...
        taylor(expand(0.5*((transpose(Omega)*J(x))*Omega)), epsilon, 'Order', 3), ...
        {epsilon, thetaz(x,t), thetay(x,t)}, ...
        {1, diff(ydef, x)-phiz, -diff(zdef, x)-phiy});
    TrDisc = subs(subs(...
        taylor(expand(0.5*((transpose(Omega)*DiscJ(x,t))*Omega)), epsilon, 'Order', 3), ...
        {epsilon, thetaz(x,t), thetay(x,t)}, ...
        {1, diff(ydef, x)-phiz, -diff(zdef, x)-phiy}), x, 0);
    % Translational Kinetic Energy per Unit Length
    fprintf('Defining Translational KE...\n')
    KEt = 0.5*rhoA(x)*(diff(ydef, t)^2+diff(zdef, t)^2);
    % Bending Potential Energy per Unit Length
    fprintf('Defining Bending PE...\n')
    PEb = 0.5*EI(x)*(diff(diff(ydef, x)-phiz, x)^2 + diff(-diff(zdef, x)-phiy, x)^2);
    % Shearing Potential Energy per Unit Length
    fprintf('Defining Shearing PE...\n')
    PEs = 0.5*KAG(x)*(phiy^2+phiz^2);
    % Net change (reduction) in rotor span due to deflection
    fprintf('Defining Axial Shortening...\n')
    deltal = int(0.5*(diff(ydef,x)^2 + diff(zdef,x)^2 - phiy^2 - phiz^2), x, 0, fp);

    syms u(t) % A new generalized coordinate that gives the horizontal location of the bearing side of the rotor, zero when the gap equals h0 and positive when the gap is larger.
end

fprintf('Defining Axial Rigid Body PE...\n')
Vc = Fax*(u(t)-deltal); % Potential energy due to the applied compressive force

if firstiteration
    
    fprintf('Defining Axial Rigid Body KE...\n')
    Tx = 0.5*eval(int(rhoA(x), x, 0, l))*(diff(u(t),t))^2; % Kinetic energy due to axial rigid-body motion of the shaft

    fprintf('Integrating Translational KE...\n')
    Tt = expand(int(KEt, x, 0, l) + Tx); % Translational kinetic energy of the shaft
    
    Tdisc = TrDisc + 0.5*DiscMass*(diff(u(t), t))^2 + ...
        0.5*DiscMass*(subs((diff(ydef, t) - DiscLoc*diff(diff(ydef,x),t))^2, x, 0) + ...
        subs((diff(zdef, t) - DiscLoc*diff(diff(zdef,x),t))^2, x, 0));
    
    fprintf('Integrating Rotational KE... ')
    
    % Parallelized integration of Rotational KE
    KErTerms = children(KEr);
    numKErTerms = length(KErTerms);
    fprintf("%d terms found.\n", numKErTerms)
    Tr = sym(0);
    parfor i = 1:numKErTerms
        Tr = Tr + int(KErTerms(i), x, 0, l)
    end

    fprintf('Adding Total KE''s...\n')
    T = Tt + Tr + Tdisc;

    fprintf('Integrating Bending and Shearing PE...\n')
    Vbs = expand(int(PEb + PEs, x, 0, l)); % Potential energy due to bending and shearing

end
    
% Total potential energy
fprintf('Adding Bearing Energy...\n')
V = simplify(Vbs + Vc + 0.5*(kb1*subs(ydef, x, sp1)^2 + kb1*subs(zdef, x, sp1)^2 + ...
    kb2*subs(ydef, x, sp2)^2 + kb2*subs(zdef, x, sp2)^2));

fprintf('Writing the Lagrangian...\n')
L = T - V; % THE LAGRANGIAN

fprintf('Writing Euler-Lagrange equations...\n')

% DEFINING EULER LAGRANGE EQUATIONS IN FULL SYMBOLIC FORM WITH ARBITRARY TIME FUNCTIONS

NumEqs = 4*modes+1; % The number of equations is equal to the number of generalized coordinates. There is one these for each assumed mode shape considered in each of two directions each for shear and for bending, plus one more, u(t).

eqns = {}; % eqns will contain the full EoM containing arbitrary time functions and their derivatives

syms Ai Aidot uu uudot

% COMPUTER ALGEBRA STRATEGY FOR VARIATIONAL DERIVATIVES:
% 1) Substitute variables in place of functions
% 2) Use diff()
% 3) Substitute functions back in in place of variables

parfor i=1:4*modes % For each generalized coordinate,
    fprintf('Writing equation %d of %d...\n', i, NumEqs);
    % Currently, the generalized coordinate A{i}(t) is an arbitrary symbolic function of time, as is its time-derivative.
    % FIRST substitute a symbolic variable in place of the time derivative of the generalized coordinate into the Lagrangian, THEN substitute another sym variable in place of the coordinate itself.
    Lconst = subs(subs(L, diff(A{i}(t),t), Aidot), A{i}(t), Ai);
    % Differentiate the variable-substituted Lagrangian with respect to the variable representing the generalized coordinate, then substitute the original arbitrary time-functions back in place of the placeholder symvars
    delLdelA = subs(diff(Lconst, Ai), {Ai, Aidot}, {A{i}(t), diff(A{i}(t),t)});
    % Differentiate the variable-substituted Lagrangian with respect to the variable representing the time-derivative of the generalized coordinate, then substitute the original arbitrary time-functions back in place of the placeholder symvars
    delLdelAdot = subs(diff(Lconst, Aidot), {Ai, Aidot}, {A{i}(t), diff(A{i}(t),t)});
    
    % The second and (modes+2)th generalized coordinates give the slope of the shaft at x=0.
    % The (2*modes+1)the and (3*modes+1)th gen. coords. give the shear angle of the shaft at x=0.
    % All of these contribute to the tilt of the gas bearing, so the bearing moment acts along (against) them.
    % For the equations corresponding to these generalized coordinates, the RHS of the Euler-Lagrange equation is the generalized force, i.e. the bearing moment.
    if ((i==2) | (i==modes+2)) | ((i==2*modes+1) | (i==3*modes+1))
        % Actually negative, but the negative sign is built into bmApprox.
        RHS = bmApprox(h0+u(t), A{i}(t));
    else
        RHS = 0;
    end
    
    eqns{i} = (diff(delLdelAdot, t) - delLdelA == RHS);
end

% Last equation, same process as above but now with respect to u(t). The RHS is the generalized force along u(t), i.e. the bearing force.
fprintf('Writing equation %d of %d...\n', NumEqs, NumEqs);
Lconst = subs(subs(L, diff(u(t), t), uudot), u(t), uu);
delLdelu = subs(diff(Lconst, uu), {uu, uudot}, {u(t), diff(u(t), t)});
delLdeludot = subs(diff(Lconst, uudot), {uu, uudot}, {u(t), diff(u(t), t)});
RHS = bfApprox((h0+u(t)), (A{2}(t)^2 + A{modes+2}(t)^2)^(1/2));
eqns{4*modes+1} = (diff(delLdeludot, t) - delLdelu == RHS);

% NOW SUBSITUTING SYMBOLIC VARIABLES IN PLACE OF ARBITRARY TIME FUNCTIONS INTO THE SYSTEM OF EQUATIONS SO THAT MATRICES CAN BE FORMED

sys = {}; %sys will contain the EoM with variable placeholers for arbitrary time functions and their derivatives

syms Aa [1 modes*4]
syms Adot [1 modes*4]
syms Addot [1 modes*4]

parfor i=1:4*modes
    fprintf('Substituting variables into equation %d of %d...\n', i, NumEqs);
    tempeq = eqns{i};
    for j=1:modes*4 % For each equation, i.e., each generalized coordinate
        % Substitute in variables for the second and first time derivatives of the generalized coordinate and the generalized coordinate, in that order, into the corresponding Euler-Lagrange equation
        tempeq = subs(subs(subs(tempeq, diff(A{j}(t), 't', 2), Addot(j)),...
            diff(A{j}(t), t), Adot(j)), A{j}(t), Aa(j));
    end
    sys{i} = tempeq;
end

syms Uu Udot Uddot

% Variable substitutions like above forthe last equation, i.e. the equation in u(t)
fprintf('Substituting variables into equation %d of %d...\n', NumEqs, NumEqs);
tempeq = eqns{4*modes+1};
tempeq = subs(subs(subs(tempeq, diff(u(t), 't', 2), Uddot), diff(u(t), t), Udot), u(t), Uu);
sys{4*modes+1} = tempeq;

fprintf('Converting equations to matrix form...\n')

X = [Aa(:); Uu];           % Variable placeholders for generalized coordinated
Xdot = [Adot(:); Udot];    % Variable placeholders for first time-derivatives of generalized coordinates
Xddot = [Addot(:); Uddot]; % Variable placeholders for second time-derivatives of generalized coordinates

% The mass matrix is the coefficient matrix in the EoM of the second derivatives of the generalized coordinates
% The gyroscopic/damping matrix is the coefficient matrix in the EoM of the first derivatives of the generalized coordinates
% The stiffness matrix is the coefficient matrix in the EoM of the generalized coordinates
M = sym(zeros(NumEqs^2, 1));
G = sym(zeros(NumEqs^2, 1));
K = sym(zeros(NumEqs^2, 1));
syms eq2expr [1 NumEqs]
parfor i = 1:NumEqs
    eq2expr(i) = lhs(sys{i}) - rhs(sys{i});
end
parfor i = 1:NumEqs^2
    cureqnum = ceil(i/NumEqs);
    curvarnum = mod(i,NumEqs) + 1
    M(i) = diff(eq2expr(cureqnum), Xddot(curvarnum));
    G(i) = diff(eq2expr(cureqnum), Xdot(curvarnum));
    K(i) = diff(eq2expr(cureqnum), X(curvarnum));
end
M = reshape(M, NumEqs, NumEqs)';
G = reshape(G, NumEqs, NumEqs)';
K = reshape(K, NumEqs, NumEqs)';

% Stacking matrices to form a generalized eigenvalue problem, thus:
% 
% /  M   G \ / qddot \   / [0]  K  \ / qdot \
% |        | |       | + |         | |      | = 0
% \ [0]  I / \ qdot  /   \ -I  [0] / \  q   /
% 
% So that
% 
%      /  M   G \
% Rr = |        |
%      \ [0]  I /
% 
%      
%      / [0]  K  \
% Ss = |         |
%      \ -I  [0] /

fprintf('Stacking matrices...\n')
Rr = [M G; zeros(NumEqs) eye(NumEqs)];
Ss = [zeros(NumEqs) K; -eye(NumEqs), zeros(NumEqs)];

fprintf('Calculating eigenvalues...\n')

% Equally distributed speeds (omegas values) at which to compute eigenvalues
speeds = linspace(0, maxspeed, plotres);

% A matrix to store eigenvalues. Each column will correspond to a spin speed from speeds
points = zeros(2*NumEqs, plotres);

% Calculating and storing eigenvalues at different spin speeds...
parfor i=1:plotres % for each spin speed
    currentspeed = speeds(i);
    fprintf('Calculating eigenvalues at spin speed %f rad/s\n', speeds(i));
    
%   System governing equations (and so Rr and Ss) should be independent of t
%   But imperfect simplification might lead to t being present in Rr and Ss
%   Arbitrary substitution t=0 is used to deal with this.
    
    RrSubbed = vpa(subs(Rr, [omegas, t], [currentspeed, 0]));   % Substitute spin speed and t into Rr
    SsSubbed = vpa(subs(Ss, [omegas, t], [currentspeed, 0]));   % Substitute spin speed and t into Ss
    EffMat = -RrSubbed\SsSubbed;                                % Find Rr\Ss i.e. inv(Rr)*Ss
    points(:,i) = sort(abs(imag(eig(EffMat))))/(2*pi);          % Imaginary parts of eigenvalues divided by 2*pi are natural frequencies
end

% Discard alternate points in each column because they repeat (one instance for the y-direction and one for the z-direction)
for i=1:NumEqs-1
    points(i,:) = points(2*i-1,:);
end
points(NumEqs:2*NumEqs,:) = zeros(NumEqs+1,plotres);

fprintf('Organizing natural frequencies into whirl lines...\n')

for i=plotres-1:-1:2 % for each column except the first and last (moving right to left: assumption is that the last two entries of each row are in the correct column)
    for j=1:2*NumEqs % for each frequency in the column
        for k=j:2*NumEqs % (slope matching) find and move to the same row the fequency in the column to the left whose difference from the current frequency is closest the difference between the current frequency and the frequency to its immediate right
            if  abs(2*points(j,i)-points(k,i-1)-points(j,i+1)) < abs(2*points(j,i)-points(j,i-1)-points(j,i+1))
                temp = points(k,i-1);
                points(k,i-1) = points(j,i-1);
                points(j,i-1) = temp;
            end
        end
    end
end

% Now the same thing going from left to right instead of right to left
for i=2:plotres-1
    for j=1:2*NumEqs
        for k=j:2*NumEqs
            if abs(2*points(j,i)-points(k,i+1)-points(j,i-1)) < abs(2*points(j,i)-points(j,i+1)-points(j,i-1))
                temp = points(k,i+1);
                points(k,i+1) = points(j,i+1);
                points(j,i+1) = temp;
            end
        end
    end
end

fprintf('Finding intersections...\n')

criticalspeeds = [];

nspeeds = modes;

for i=1:nspeeds
    for j=1:plotres-1
        % If one natural frequency in a line is less than the spin speed and the next is greater, or if the one is greater and the next is less
        if (points(i,j) - speeds(j)/(2*pi)) * (points(i,j+1) - speeds(j+1)/(2*pi)) < 0
            % Find the point of intersection using a linear interpolation formula
            x1 = speeds(j);
            x2 = speeds(j+1);
            y1 = points(i, j);
            y2 = points(i, j+1);
            criticalspeeds(end+1) = (2*pi*(x1*y2 - x2*y1))/(x1 - x2 - 2*pi*y1 + 2*pi*y2);
        end
    end
end

% SHOWING CRITICAL SPEEDS AND PLOTTING CAMPBELL DIAGRAM

criticalspeeds

fprintf('Plotting Campbell diagram...\n')

figure(2)

clf

plot([0,maxspeed], [0, maxspeed/(2*pi)], '-k')
hold on
plot(criticalspeeds, criticalspeeds/(2*pi), 'k^')

for i=1:nspeeds
    if points(i,1) == 0
        continue
    elseif points(i,1)<points(i,plotres)
        plot(speeds, points(i,:), '-r')
    else
        plot(speeds, points(i,:), '-b')
    end
end

grid minor

ylim([0 maxspeed/(2*pi)])
xlim([0 maxspeed])
title('Bearing Stiffnesses '+string(sprintf('%.1e', kb1))+...
    ' N/m, Load = '+string(sprintf('%.1f', Fax))+' N')
xlabel('Spin Speed (rad/s)')
ylabel('Natural Frequency (Hz)')

shg