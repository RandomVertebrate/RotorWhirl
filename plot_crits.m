data = readmatrix("critspeeds.txt")';

horz_axis = data(2, 1:end-1)

data = rowData2lines(data(5:end, 1:end-1));

plot(horz_axis, data, 'k')
ylabel("Critical Speeds (rad/s)")
xlabel("Clearance (m)")
grid minor
xlim([min(horz_axis) max(horz_axis)])