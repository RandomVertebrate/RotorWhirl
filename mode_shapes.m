function mode_shapes(l, EigVec_matrix, EigVal_matrix, Qshapes, outerradius, currentspeed, maxnodes, maxmodes, DiscRadius, filename)

syms x
n = double(maxmodes);                         % no of modes shapes to plot
outfig = figure();

modes = length(Qshapes);
EigVec = EigVec_matrix;
EigVal = EigVal_matrix;
EigVal_matrix = imag(EigVal);
EigVec_imag = imag(EigVec);
EigVec_real = real(EigVec);
mat_s = size(EigVec_imag(:,1));
mat_size = mat_s(1);
EigValue_vec = diag(EigVal_matrix,0);
EigValue_vec_sort = sort(EigValue_vec);
Vec_size = mat_size/2;
coeff_def = zeros(modes,n);
coeff_shr = zeros(modes,n);
E_imag = zeros(mat_size,1);
E_real = zeros(mat_size,1);

for i = Vec_size+1:Vec_size+n
    ind = find(round(EigValue_vec,6) == round(EigValue_vec_sort(i),6));
    size_index = length(ind);
    for j = 1:size_index
        aa = 0;
        if aa == 1
        else
            index = ind(j);
            E_imag = EigVec_imag(:,index);
            E_real = EigVec_real(:,index);
            if E_imag(Vec_size+modes) ~= 0
                aa = 1;
            else
                aa = 0;
                index = min(ind);
            end
        end
    end
    E_imag = EigVec_imag(:,index);
    E_real = EigVec_real(:,index);
    if E_imag(Vec_size+modes) ~= 0
        coeff_def(:,i-Vec_size) = E_imag(Vec_size+modes:Vec_size+2*modes-1);
        coeff_shr(:,i-Vec_size) = E_imag(Vec_size:Vec_size+modes-1);
    else
        coeff_def(:,i-Vec_size) = E_real(Vec_size+modes:Vec_size+2*modes-1);
        coeff_shr(:,i-Vec_size) = E_real(Vec_size:Vec_size+modes-1);
    end
end

npts = 500;
scgap = 10;

x_vals = linspace(0, l, npts);

shaftprofile = double(outerradius(x_vals));
shaftprofile(1) = DiscRadius;
shaftdia = max(shaftprofile);

if maxmodes>8
    combined_plot = tiledlayout(ceil(maxmodes/3),3);
else
    combined_plot = tiledlayout('flow');
end
combined_plot.TileSpacing = 'compact';
combined_plot.Padding = 'compact';

for i = 1:n
    
    def_shape = 0;
    shr_shape = 0;

    for j = 1:modes
        def_shape = def_shape+coeff_def(j,i)*Qshapes(j);
        shr_shape = shr_shape+coeff_shr(j,i)*Qshapes(j);
    end

    y(i) = def_shape;
    y_va = double(subs(def_shape, x, x_vals));
    sh_va = double(subs(shr_shape, x, x_vals));
    y_va_max = max(abs(y_va));
    sh_va_max = max(abs(sh_va));

    if y_va_max == 0
        y_val = y_va;
    else
        y_val = shaftdia*y_va./y_va_max;   
    end
    if sh_va_max == 0
        sh_val = sh_va;
    else
        sh_val = shaftdia*sh_va./sh_va_max;   
    end

    spin_axis_angle = (diff(y_val)./diff(x_vals))+sh_val(1:end-1);
    maxspangle = max(abs(spin_axis_angle));
    meanspanglediff = mean(abs(diff(spin_axis_angle)./diff(x_vals(1:end-1))))/maxspangle;
    if maxspangle ~= 0
        anglescale = (meanspanglediff^2 + meanspanglediff)/(meanspanglediff^2 + 1);
        spin_axis_angle = anglescale*(pi/4)*spin_axis_angle/maxspangle;
    end
    spin_axis_angle(end+1) = spin_axis_angle(end);

    profile_xshift = -shaftprofile.*sin(spin_axis_angle);
    profile_yshift = shaftprofile.*cos(spin_axis_angle);

    zero_crossings = diff(sign(y_val)) ~= 0;
    zero_x_vals = x_vals(zero_crossings);
    
    num_nodes = length(zero_x_vals);
    if num_nodes <= maxnodes % && y_va_max~=0
        modeshape_ax = nexttile;
        plot(x_vals, y_val, "LineWidth", 1.0)
        hold on
        plot(x_vals, 0*x_vals,":")

        plot([x_vals(1)+profile_xshift(1),x_vals(1)-profile_xshift(1)], ...
            [y_val(1)+profile_yshift(1),y_val(1)-profile_yshift(1)], ...
            "-","Color","k","LineWidth",1.5)
        for cr_sc = 1:scgap:npts
            plot([x_vals(cr_sc)+profile_xshift(cr_sc),x_vals(cr_sc)-profile_xshift(cr_sc)], ...
                [y_val(cr_sc)+profile_yshift(cr_sc),y_val(cr_sc)-profile_yshift(cr_sc)], ...
                "-","Color",double([0.5 0.5 0.5]))
        end
        plot([x_vals(end)+profile_xshift(end),x_vals(end)-profile_xshift(end)], ...
            [y_val(end)+profile_yshift(end),y_val(end)-profile_yshift(end)], ...
            "-","Color","k")

        plot(x_vals+profile_xshift,y_val+profile_yshift,'-k')
        plot(x_vals-profile_xshift,y_val-profile_yshift,'-k')
            
        scatter(zero_x_vals, zeros(size(zero_x_vals)), 'ro', 'LineWidth', 1); % Zero-crossing points
        title(string(round(EigValue_vec_sort(Vec_size+i)/(2*pi)))+" Hz")
        yticks(0)
        hold off
        grid on
        box on
        axis equal
        xlim([-l*0.1 l*1.1])
        ylim([-2*shaftdia 2*shaftdia])
    end
end

title(combined_plot, "Modeshapes at "+string(round(currentspeed,1))+" rad/s")
saveas(outfig, filename)
close(outfig)

end