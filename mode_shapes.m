function modeshapes(rotorlength, EigVec_matrix, EigVal_matrix, Qshapes, outerradius, currentspeed, filename)

% Assumes spatial coordinate x
syms x

num_assumed_modes = length(Qshapes)

outfig = figure();

% Number of modes to discard (highest modes will be inaccurate)
%discard_modeshapes = floor(num_assumed_modes/2)
%discard_modeshapes = 0
max_modeshapes = 6;

% Populate displacement modeshape function coefficient arrays
for i=1:8*num_assumed_modes+2
    idx_range = 4*(num_assumed_modes+2):4*(num_assumed_modes+2)+num_assumed_modes-1;
    shapecoeffs(:,i) = EigVec_matrix(idx_range,i);
    Fnat_array(i) = abs(imag(EigVal_matrix(i,i)));
end

for i=1:2*num_assumed_modes+2
    shapecoeffs(:, i) = shapecoeffs(:, i)/norm(shapecoeffs(:,i));
end

% delete NaN modes
realmodes = ~isnan(shapecoeffs(1,:))
shapecoeffs = shapecoeffs(:,realmodes)
Fnat_array = Fnat_array(realmodes)
% only keep modes with unique frequencies
[Fnat_array, uq_idx, ~] = uniquetol(double(Fnat_array));
shapecoeffs = shapecoeffs(:,uq_idx);
num_modeshapes = length(uq_idx)
num_modeshapes == length(Fnat_array);
    
shapecoeffs = imag(shapecoeffs);        % Use imaginary parts only
szcoeff = size(shapecoeffs);
num_modeshapes == szcoeff(2)

% Symbolic modeshapes
for i=1:num_modeshapes
    modeshape(i) = Qshapes*shapecoeffs(:,i);
end

modeshape = transpose(simplify(modeshape));     % Convert to column vector

% sort modes by frequency
[Fnat_array, mode_order] = sort(Fnat_array);
modeshape = modeshape(mode_order);

% Discretize, normalize, plot and save modeshapes
xvals = linspace(0,rotorlength,100);
shaftradius = outerradius(xvals);
mdsps = subs(modeshape,x,xvals);
if num_modeshapes > max_modeshapes
    num_modeshapes = max_modeshapes;
end
% normalizing...
shaftradius = shaftradius/max(shaftradius);
for i=1:num_modeshapes
    scale = max(abs(mdsps(i,:)));
    if scale ~= 0
        mdsps(i,:) = mdsps(i,:)/scale;
    end
end
% plotting...
combined_plot = tiledlayout(num_modeshapes, 1)
for i=1:num_modeshapes
    nexttile;
    plot(xvals, mdsps(i,:), "LineWidth", 3.0)
    hold on
    plot(xvals, mdsps(i,:)+shaftradius,"--k")
    plot(xvals, mdsps(i,:)-shaftradius,"--k")
    title(string(round(Fnat_array(i)/(2*pi)))+" Hz")
    xlim([0 rotorlength])
    yticks(0)
    grid on
    box off
end
title(combined_plot, "Modeshapes at "+string(round(currentspeed,1))+" rad/s")
posvec = get(outfig, "Position");
posvec(3) = 300;
posvec(4) = 100*(num_modeshapes);
set(outfig, "Position", posvec)
saveas(outfig, filename)
close(outfig)

end
