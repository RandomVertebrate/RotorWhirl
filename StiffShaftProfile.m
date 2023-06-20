function ORexp = StiffShaftProfile(xcoord, baseradius, steploc, stepval, l, slope)

% This function finds and returns a symbolic expression for a stepped
% shaft's effective profile as per the 45-degree rule (or any other angle).
% xcoord is the symbolic variable in terms of which to return the output
% (axial spatial coordinate)
% baseradius is the radius of the shaft at the left end
% steploc is a row vector of step locations (x-coordinates)
% steval is a row vector of step values. Positive values are increases in
% radius and negative values are decreases.
% l is the length of the shaft.

margin = max(stepval)*1e-8;                 % Numerical margin for floating-point comparisons

numsteps = length(steploc);

stepwidth = [steploc(1) diff(steploc) l-steploc(end)];

% Will hold all points on the true shaft profile immediately after a step
% down, from where a north-westerly sloped line might emanate
rampdownendpts = zeros(2, 0);

% Will hold all points on the true shaft profile immediately before a step
% up, from where a north-easterly sloped line might emanate
rampupstartpts = zeros(2, 0);

% Finding all points along the true shaft profile
% (Including [0; 0] and [0, l] to "close" the profile)
truepts = [0, 0; 0, baseradius];
for i = 1:numsteps
    truepts(:, end+1) = truepts(:, end) + [stepwidth(i); 0];
    truepts(:, end+1) = truepts(:, end) + [0; stepval(i)];
end
truepts(:, end+1) = truepts(:, end) + [stepwidth(end); 0];
truepts(:, end+1) = [l, 0];
ncnr = size(truepts);
numtruepts = ncnr(2);             % Number of points on the true shaft profile

% Each column will hold the end points [x1; y1; x2; y2] of a line segment
linesegs = zeros(4, 0);
% Adding true-profile segments to linesegs
for i = 1:numtruepts-1
    linesegs(:, end+1) = [truepts(:, i); truepts(:, i+1)];
end

% allpts will hold all true-profile points except [0; 0] and [0; l]
% AND all possible intersections between upramps, downramps and the true
% profile
allpts = truepts(:, 2:end-1);

% Collecting upramp start points and downramp end points
for i = 1:numsteps
    if stepval(i) > 0
        rampupstartpts(:, end+1) = truepts(:, 2*i+1);
    else
        rampdownendpts(:, end+1) = truepts(:, 2*i+2);
    end
end
ncnr = size(rampupstartpts);
numupramps = ncnr(2);
ncnr = size(rampdownendpts);
numdownramps = ncnr(2);

% Finding intersections of upramps with horizontal and vertical lines
for i = 1:numupramps             % Iterating through upramps
    % Horizontal
    for j = 1:numsteps+1         % Iterating through all horizontal lines
        % Potential ramp height if the current upramp intersected with the current horizontal line
        rampheight = truepts(2, 2*j) - rampupstartpts(2, i);
        % Potential ramp width
        rampwidth = rampheight/slope;
        % Potential crossing point
        crossing = rampupstartpts(:, i) + [rampwidth; rampheight];
        % IF the crossing point lies within the limits of the horizontal line segment
        if crossing(1) > truepts(1, 2*j) & crossing(1) <= truepts(1, 2*j+1)
            % Add the crossing point to allpts
            allpts(:, end+1) = crossing;
            % Add the upramp line segment to linesegs
            linesegs(:, end+1) = [rampupstartpts(:, i); crossing];
        end
    end
    % Vertical
    % An upramp will never intersect the first vertical line
    for j = 2:numsteps+2        % Iterating through all vertical lines except the first
        % Potential ramp width if current upramp intersects with the current vertical line
        rampwidth = truepts(1, 2*j-1) - rampupstartpts(1, i);
        % Potential ramp height
        rampheight = rampwidth*slope;
        % Potential crossing point
        crossing = rampupstartpts(:, i) + [rampwidth; rampheight];
        % IF the crossing point lies within the limits of the vertical line segment
        if (crossing(2) > truepts(2, 2*j-1) & crossing(2) <= truepts(2, 2*j)) ...
                | (crossing(2) <= truepts(2, 2*j-1) & crossing(2) > truepts(2, 2*j))
            % Add the crossing point to allpts
            allpts(:, end+1) = crossing;
            % Add the upramp line segment to linesegs
            linesegs(:, end+1) = [rampupstartpts(:, i); crossing];
        end
    end
end

% Finding intersections of downramps with horizontal and vertical lines
for i = 1:numdownramps          % Iterating through downramps
    % Horizontal
    for j = 1:numsteps          % Iterating through horizontal lines
        % Potential ramp height
        rampheight = truepts(2, 2*j) - rampdownendpts(2, i);
        % Potential ramp width
        rampwidth = rampheight/slope;
        % Potential crossing point
        crossing = rampdownendpts(:, i) + [-rampwidth; rampheight];
        % IF crossing within horizontal line segment limits
        if crossing(1) > truepts(1, 2*j) & crossing(1) <= truepts(1, 2*j+1)
            % Add crossing to allpts
            allpts(:, end+1) = crossing;
            % Add downramp to linesegs
            linesegs(:, end+1) = [rampdownendpts(:, i); crossing];
        end
    end
    % Vertical
    % A downramp will never intersect the last vertical line
    for j = 1:numsteps+1        % Iterating through all vertical lines except the last 
        % Potential ramp width
        rampwidth = rampdownendpts(1, i) - truepts(1, 2*j-1);
        % Potential ramp height
        rampheight = rampwidth*slope;
        % Potential crossing
        crossing = rampdownendpts(:, i) + [-rampwidth; rampheight];
        % IF crossing within vertical line segment limits
        if (crossing(2) > truepts(2, 2*j-1) & crossing(2) <= truepts(2, 2*j)) ...
                | (crossing(2) <= truepts(2, 2*j-1) & crossing(2) > truepts(2, 2*j))
            % Add crossing to allpts
            allpts(:, end+1) = crossing;
            % Add downramp to linesegs
            linesegs(:, end+1) = [rampdownendpts(:, i); crossing];
        end
    end
end

% Finding intersections between upramps and downramps
for i = 1:numupramps                 % Iterating through upramps
    for j = 1:numdownramps           % Iterating through downramps
        % If ramps can intersect, i.e. downramp ends after upramp starts
        if rampupstartpts(1, i) < rampdownendpts(1, j)
            % Calculating crossing
            vertsep = rampdownendpts(2, j) - rampupstartpts(2, i);
            horzsep = rampdownendpts(1, j) - rampupstartpts(1, i);
            rampwidth = (vertsep + slope*horzsep)/(2*slope);
            rampheight = rampwidth*slope;
            crossing = rampupstartpts(:, i) + [rampwidth; rampheight];
            % Add crossing to allpts
            allpts(:, end+1) = crossing;
            % Add upramp to linesegs
            linesegs(:, end+1) = [rampupstartpts(:, i); crossing];
            % Add downramp to linesegs
            linesegs(:, end+1) = [crossing; rampdownendpts(:, j)];
        end
    end
end

% Total number of points (excluding [0; 0] and [l; 0])
ncnr = size(allpts);
totnumpts = ncnr(2);

%Total number of line segments
ncnr = size(linesegs);
numlinesegs = ncnr(2);

% This will hold indices (in allpts) of points to be deleted
deletelist = [];

% If multiple points have the same x-coordinate, only keep the smallest y-coordinate
for i = 1:totnumpts              % For each point
    for j = 1:totnumpts          % Look at every other point
        heightdiff = allpts(2, j) - allpts(2, i);
        widthdiff = allpts(1, j) - allpts(1, i);
        % If the other point is directly overhead, mark it to be deleted
        if abs(widthdiff) < margin & heightdiff > margin
            % (Only mark to be deleted if not already marked)
            if ~ismember(j, deletelist)
                deletelist(end+1) = j;
            end
        end
    end
end

% If any point is above any line segment, delete that point.
% (Only the lowest line segments will be kept)
for i = 1:totnumpts              % For each point
    for j = 1:numlinesegs        % For each line segment
        % If the point is within the horizontal limits of the line segment
        if (allpts(1, i) > linesegs(1, j) & allpts(1, i) <= linesegs(3, j))...
                | (allpts(1, i) > linesegs(3, j) & allpts(1, i) <= linesegs(1, j))
            % And if the piint's y-coordinate is greater than obtained by extending the line segment to the horizontal location of the point
            if allpts(2, i) > margin + linesegs(2, j) + (linesegs(4, j) - linesegs(2, j))*(allpts(1, i) - linesegs(1, j))/(linesegs(3, j) - linesegs(1, j))
                % And if the point is not already marked to be deleted
                if ~ismember(i, deletelist)
                    % Mark the point to be deleted
                    deletelist(end+1) = i;
                end
            end
        end
    end
end

% Number of points to be deleted
numdel = length(deletelist);

% DELETING by copying points from allpts to somepts and omitting those marked to be deleted
somepts = [];
for i = 1:totnumpts
    if ~ismember(i, deletelist)
        somepts(:, end+1) = allpts(:, i);
    end
end

% Sorting remaining points by x-coordinate
[xvals, order] = sort(somepts(1, :));
yvals = somepts(2, order);

% Removing redundant points
finalpoints = uniquetol([xvals' yvals'], margin, 'ByRows', true);
xvals = finalpoints(:, 1)';
yvals = finalpoints(:, 2)';

% Plotting result
figure()
plot(truepts(1, :), truepts(2, :), "--k")
axis image
hold on
plot(rampupstartpts(1, :), rampupstartpts(2, :), "*k")
plot(rampdownendpts(1, :), rampdownendpts(2, :), "^k")
plot(allpts(1, numsteps*2+1:end), allpts(2, numsteps*2+1:end), "+k")
i = 2*numsteps+2;
plot([linesegs(1, i), linesegs(3, i)], [linesegs(2, i), linesegs(4, i)], ':k')
for i = 2*numsteps+3:numlinesegs
    plot([linesegs(1, i), linesegs(3, i)], [linesegs(2, i), linesegs(4, i)], ':k', 'HandleVisibility', 'off')
end
plot(xvals, yvals, '-k', 'LineWidth', 1.5)
grid minor
title("Effective Shaft Profile")
xlabel("x")
ylabel("Radius")
legend(["True Profile"; "Upramp Start"; "Downramp End"; "Intersections"; "Ramps"; "Effective Profile"])
shg

% Building symbolic output
ORexp = sym(yvals(1));
numterms = length(xvals);
for i = 1:numterms-1
    ORexp = ORexp + heaviside(xcoord-xvals(i))*heaviside(xvals(i+1)-xcoord) ...
        *((xcoord-xvals(i))*(yvals(i+1)-yvals(i))/(xvals(i+1)-xvals(i))) ...
        + heaviside(xcoord-xvals(i+1))*(yvals(i+1)-yvals(i));
end
ORexp = simplify(ORexp, 'Steps', 10);

end