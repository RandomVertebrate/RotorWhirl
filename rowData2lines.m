function points = rowData2lines(points)

% This function attempts tp rearrange every column of a matrix so that the
% rows make sense when plotted against something evenly spaced.

% Since there is no inherent ordering to eigenvalues, study of the change
% in eigenvalues as a function of spin speed, load, gap, or any other
% system parameter requires mapping between individual eigenvalues from
% different sets calculated independently. One simple way to do this is
% through effective minimization of the change in slope of
% parameter-vs-eigenvalue lines at every point. Assuming equally
% distributed parameter values, this reduces to choosing a mapping that
% minimizes the second difference of the data. Since an exhaustive search
% is silly, this assumes the values in the first and last columns start off
% in place.

% Data must be along rows, i.e. each row taken to be a line

[num_rows, num_columns] = size(points);

for i=num_columns-1:-1:2 % for each column except the first and last (moving right to left: assumption is that the last two entries of each row are in the correct column)
    for j=1:2*num_rows % for each value in the column

        % (slope matching) find and move to the same row the value in the 
        % column to the left whose difference from the current value is 
        % closest to the difference between the current value and the 
        % value to its immediate right
        for k=j:num_rows
            if  abs(2*points(j,i)-points(k,i-1)-points(j,i+1)) < abs(2*points(j,i)-points(j,i-1)-points(j,i+1))
                temp = points(k,i-1);
                points(k,i-1) = points(j,i-1);
                points(j,i-1) = temp;
            end
        end
    end
end

% Now the same thing going from left to right instead of right to left
for i=2:num_columns-1
    for j=1:num_rows
        for k=j:num_rows
            if abs(2*points(j,i)-points(k,i+1)-points(j,i-1)) < abs(2*points(j,i)-points(j,i+1)-points(j,i-1))
                temp = points(k,i+1);
                points(k,i+1) = points(j,i+1);
                points(j,i+1) = temp;
            end
        end
    end
end

end