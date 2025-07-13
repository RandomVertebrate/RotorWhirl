
% Uses many variables defined in RUN file
% To be called by RUN file only

% Start parallel processing pool
delete(gcp('nocreate'));
parpool(processorcores);

if firstiteration
    
    % Drawing rotor geometry...
    xvals = linspace(0, l, 1000);
    holexvals = linspace(0.999*holestart, 1.001*holestop, 1000);
    zeroradius = double(outerradius(0));
    zeroholeradius = double(innerradius(0));
    endradius = double(outerradius(l));
    endholeradius = double(innerradius(l));
    figure()
    clf
    trueprofileplot = plot(xvals, outerradius(xvals), 'k');
    hold on
    stiffprofileplot = plot(xvals, outerradius2(xvals), '--k');
    plot(xvals, -outerradius(xvals), 'k')
    plot(xvals, -outerradius2(xvals), '--k')
    discplot = fill([-DiscThickness, 0, 0, -DiscThickness], ...
        [DiscRadius, DiscRadius, -DiscRadius, -DiscRadius], 'k');
    plot(holexvals, innerradius(holexvals), 'k')
    plot(holexvals, -innerradius(holexvals), 'k')
    plot([0 0], [-zeroradius -zeroholeradius], 'k', ...
        [0 0], [zeroradius zeroholeradius], 'k', ...
        [l l], [-endradius -endholeradius], 'k', ...
        [l l], [endradius endholeradius], 'k')
    splocplot = plot([sploc, sploc], [outerradius(sploc), -outerradius(sploc)], 'ob');
    lmpplot = scatter(lmploc, lmploc*0, 1000*DiscRadius*lmpval/max(lmpval), "filled");
    forcearrow = quiver(fp+10e-3, 0, -l/10, 0, 'r');
    forcearrow.LineWidth = 2;
    forcearrow.MaxHeadSize = 1;
    ylim([-DiscRadius*1.1, DiscRadius*1.1])
    xlim([-DiscThickness-0.1*(l+DiscThickness), 1.1*(l+DiscThickness)])
    title('Rotor Geometry (Meters)')
    legend([trueprofileplot; stiffprofileplot; discplot; splocplot; lmpplot; forcearrow], ...
           ["True shaft profile"; "Effective rigidity profile"; ...
           "Aerostatic bearing disc"; "Lateral support bearing"; ...
           "Lumped Mass"; "Axial load"])
    hold off
    axis equal
    grid minor
    shg
    % End drawing rotor geometry
end

syms h theta real
syms k

% CALCULATING BEARING BEHAVIOR

% angle between bearing holes
phigap = 2*pi/nholes;
% Total bearing force (load capacity as a function of gap and tilt)
BearingForce = @(h, theta) symsum(C1/(1+C2*(h+rholes*theta*sin(phi0+k*phigap))^(C3)), k, 0, nholes-1);
% Total bearing moment
BearingMoment = @(h, theta) symsum(C1*rholes*sin(phi0+k*phigap)/(1+C2*(h+rholes*theta*sin(phi0+k*phigap))^(C3)), k, 0, nholes-1);

Fax = subs(BearingForce(h,0), h, h0);
NetFlatForce = @(h) eval(Fax - BearingForce(h,0));      % Single-argument load function to solve for operating clearance

TotalStiffness = subs(diff(NetFlatForce(h), h), h, h0); % Total translational film stiffness at operating clearance

fprintf('\nCurrent load = %f N\nCurrent operating clearance = %f m\nCurrent total film stiffness = %f N/m\nForce residual = %f N\n\n', ...
    Fax, h0, TotalStiffness, NetFlatForce(h0))

syms epsilon
syms x t real

% Local linearization of bearing behavior
bfApprox = @(h,theta) subs(taylor(BearingForce(h0+epsilon*(h-h0), epsilon*theta), epsilon, 0, 'Order', 2), epsilon, 1);
bmApprox = @(h,theta) subs(taylor(BearingMoment(h0+epsilon*(h-h0), epsilon*theta), epsilon, 0, 'Order', 2), epsilon, 1);

if firstiteration

    % CALCULATING ROTOR PROPERTIES

    rhoA = @(x) expand(eval(materialDensity*pi*(outerradius(x)^2-innerradius(x)^2)));           % Linear Mass Density

    EI = @(x) expand(eval(youngsModulus*pi*(outerradius2(x)^4-innerradius(x)^4)/4));           % Flexural Rigidity

    KAG = @(x) expand(eval(0.9*pi*(outerradius(x)^2-innerradius(x)^2)*shearModulus));   % Effective Shear Rigidity

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
    PEb = 0.5*EI(x)*expand(diff(diff(ydef, x)-phiz, x)^2) + 0.5*EI(x)*expand(diff(-diff(zdef, x)-phiy, x)^2);
    % Shearing Potential Energy per Unit Length
    fprintf('Defining Shearing PE...\n')
    PEs = 0.5*KAG(x)*expand(phiy^2) + 0.5*KAG(x)*expand(phiz^2);
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
    
    % Kinetic Energy due to lumped masses
    fprintf('Defining Kinetic Energy of Lumped Masses...\n')
    Tlmp = sym(0);
    for i=1:num_lmp
        Tlmp = Tlmp + 0.5*lmpval(i)*(diff(u(t), t))^2 + ...
        0.5*lmpval(i)*(subs((diff(ydef, t))^2+(diff(zdef, t))^2, x, lmploc(i)));
    end

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
    T = Tt + Tr + Tdisc + Tlmp;

    fprintf('Integrating Bending and Shearing PE... ')
    % Parallelized integration of bending and shearing energy
    PEbsTerms = children(PEb + PEs);
    numPEbsTerms = length(PEbsTerms);
    fprintf("%d terms found.\n", numPEbsTerms)
    Vbs = sym(0);
    parfor i = 1:numPEbsTerms
        Vbs = Vbs + int(PEbsTerms(i), x, 0, l)
    end
end
    
% Total potential energy
V = simplify(Vbs + Vc);
fprintf('Adding Bearing Energy...\n')
num_bearings = length(sploc);
for i = 1:num_bearings
    V = V + 0.5*(...
        kbvals(i)*subs(ydef, x, sploc(i))^2 + ...
        kbvals(i)*subs(zdef, x, sploc(i))^2);
end

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
        % Substitute in variables for the second and first time derivatives of
        % the generalized coordinate and the generalized coordinate, in that
        % order, into the corresponding Euler-Lagrange equation
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

X = [Aa(:); Uu];           % Variable placeholders for generalized coordinates
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
    fprintf('Calculating eigenvalues at spin speed %f rad/s ...\n', speeds(i));
    
    %   System governing equations (and so Rr and Ss) should be independent of t
    %   But imperfect simplification might lead to t being present in Rr and Ss
    %   Arbitrary substitution t=0 is used to deal with this.
    
    RrSubbed = vpa(subs(Rr, [omegas, t], [currentspeed, 0]));   % Substitute spin speed and t into Rr
    SsSubbed = vpa(subs(Ss, [omegas, t], [currentspeed, 0]));   % Substitute spin speed and t into Ss
    EffMat = -RrSubbed\SsSubbed;                                % Find Rr\Ss i.e. inv(Rr)*Ss
    try
        [EigVec, EigVal] = eig(EffMat);
        points(:,i) = sort(abs(imag(diag(EigVal))))/(2*pi);     % Imaginary parts of eigenvalues divided by 2*pi are natural frequencies
    catch
        points(:,i) = NaN;
        fprintf('Eigenvalue calculation failed at spin speed %f rad/s\n', speeds(i));
    end    
end

% Fill in any NaN data in points (where eig failed) using neighboring datapoints
for i = 1:plotres-1                             % For each nonfinal column in points
    nanlist = isnan(points(:,i));               % Find NaN entries in column
    points(nanlist,i) = points(nanlist,i+1);    % Fill in NaN values from next column
end
% Fill in NaN data in last column from previous column
nanlist = isnan(points(:,end));
points(nanlist,end) = points(nanlist,end-1);

% Serial loop to allow global access to EigVec, EigVal
for i = 1:modespeeds
    currentspeed = (i-1)*maxspeed/modespeeds;
    RrSubbed = vpa(subs(Rr, [omegas, t], [currentspeed, 0]));   % Substitute spin speed and t into Rr
    SsSubbed = vpa(subs(Ss, [omegas, t], [currentspeed, 0]));   % Substitute spin speed and t into Ss
    EffMat = -RrSubbed\SsSubbed;                                % Find Rr\Ss i.e. inv(Rr)*Ss
    try
        [EigVec, EigVal] = eig(EffMat);
        fprintf('Drawing modeshapes at spin speed %f rad/s\n', speeds(i));
        modeshape_filename = 'Modeshapes/kb'+string(sprintf('%.1e', kbvals(1)))+...
            'Fax'+string(sprintf('%.1f', Fax))+'omega'+string(round(currentspeed))+...
            '.svg'
        mode_shapes(l,EigVec,EigVal,Qshapes,outerradius,currentspeed,Fax,h0,maxnodes,max_num_modeshapes,DiscRadius,modeshape_filename)
    catch
        fprintf('Drawing modeshapes failed at spin speed %f rad/s\n', speeds(i));
    end
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

if firstiteration
    CampbellFig = figure();
else
    figure(CampbellFig)
end

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
title('Air gap = '+string(sprintf('%d', round(h0*1e6)))+...
    ' \mu{m}, Load = '+string(sprintf('%.1f', Fax))+' N')
xlabel('Spin Speed (rad/s)')
ylabel('Natural Frequency (Hz)')

shg