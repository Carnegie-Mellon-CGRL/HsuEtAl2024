%objective function for radius and compliance only - no inflammation
function J = dual_objective(rad_input, comp_input)
%Calculate objective and write to file
%in outfile = "sse.txt"

%Read in parameter values to evaluate from inpfile
R = importdata('Native_in_');
native_rad = R.data(1,1);

native_comp = 0.064;

%Calculate the objective function value for that parameter set
GnR_outputs_changed = ((rad_input-native_rad)/native_rad).^2;
rad_function = sqrt((1/length(rad_input))*(sum(GnR_outputs_changed)));

comp_output_changed = ((comp_input-native_comp)/native_comp).^2;
comp_function = sqrt((1/length(comp_input))*(sum(comp_output_changed)));


objective_function = (rad_function/(0.201615159212792*2) + comp_function/(0.823762611015949*2));

if (isnan(objective_function))
    objective_function = 2;
end

% Creating file to be written to
filename = fopen('sse.txt','w');
%Writing data to file
fprintf(filename, '%f\n', objective_function);
%Creating output
J = objective_function;
%Closing file
fclose(filename);
end