data = readmatrix("critspeeds.txt")';

horz_axis = data(2, 1:end-1)

data = rowData2lines(data(5:end, 1:end-1));
ncrits = size(data,1)

% Export processed/aligned data
outfile = "critspeeds_processed.xls";
writematrix(transpose(["Gap (m)";"Critical Speeds"]),outfile,"Range","A1")
writematrix(transpose(horz_axis),outfile,"Range","A2")
writematrix(transpose(data),outfile,"Range","B2")

tl = tiledlayout(1,2);

% Plot critical speeds against gap
nexttile
plot(horz_axis*1e6, data, 'r')
axis tight
grid minor

nexttile
stackedplot(horz_axis*1e6, transpose(data), 'r', "DisplayLabels", repelem("",ncrits))
grid on

xlabel(tl,"Clearance (\mu{m})")
ylabel(tl,"Critical Speed (rad/s)")

set(gcf,"Position",[750 350 600 460])