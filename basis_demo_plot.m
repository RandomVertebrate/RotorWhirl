clear all
syms x k h theta C1 C2 C3

f(x) = C1/(1+C2*x^(C3))
slope_f(x) = diff(f(x), x)

minheight = 20e-6;
maxheight = 100e-6;
npts = 100;
gap_values = linspace(minheight, maxheight, npts);

plot_threshold = [5 60];                   % limits (in % of C1) on average function value

C1vals = 1e3*linspace(0.5, 1, 3)
C2vals = 1e6*2.^(linspace(0, 40, 8))
C3vals = [linspace(1, 1.75, 5), linspace(2, 4, 3)]

figure()
set(gcf, "DefaultAxesLineStyleOrder", ["-", "--", "-.", ":"])
legend_entries = [];
hold on

for C1val = C1vals
    for C2val = C2vals
        for C3val = C3vals
            
            thisline = subs(f(gap_values), [C1, C2, C3], [C1val, C2val, C3val]);
            mean_val = mean(thisline);
    
            if mean_val>C1val*plot_threshold(1)/100 && mean_val<C1val*plot_threshold(2)/100
                plot(gap_values, thisline, "linewidth", 1)
                legend_entries = [legend_entries; ...
                    "C1 = "+num2str(C1val, "%.2e")+", C2 = "+ ... 
                     num2str(C2val, "%.2e")+", C3 = "+num2str(C3val, "%.2f")];
            else
                C3val
            end
        end
    end
end

xlim([minheight, maxheight])
xlabel("$h$", "Interpreter", "latex")
ylabel("$\frac{C_1}{C_2h^{C_3}+1}$", "Interpreter", "latex")

legend(legend_entries, "Location", "eastoutside")