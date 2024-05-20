% constraint function keeping radius greater than native radius
function output = dx_constraint(input)
%Read in parameter values to evaluate from inpfile
B = importdata('Native_in_');
native_value = B.data(1,1);

%Calculate the objective function value for that parameter set

if length(input) == 1
    cons = -0.1;
else
    cons = -0.1;
    norm = input;
    norm(end) = [];
    dx_input = diff(input)./norm;
    for i = 1:length(dx_input)
        if dx_input(i) < 0 && norm(i) < native_value
            cons = cons + 1;
        end
    end
end


%Write the objective function value to file Jfile
% Exporting data to .dat
% Creating file to be written to
filename = fopen('cons.txt','w');
%Writing data to file
fprintf(filename, '%f %f\n', cons);
%Creating output
output = cons;
%Closing file
fclose(filename);
