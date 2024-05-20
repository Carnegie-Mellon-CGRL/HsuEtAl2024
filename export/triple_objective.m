%objective function for root mean square error radius and compliance
function J = triple_objective(rad_input, comp_input, infl_input)
%Calculate objective and write to file
%in outfile = "sse.txt"

%Read in parameter values to evaluate from inpfile
R = importdata('Native_in_');
native_rad = R.data(1,1);

native_comp = 0.064;

exp = importdata('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/master/GnR_out_experimental_');
exp_infl = exp(:,6) + exp(:,7);

%Calculate the objective function value for that parameter set
%need to normalize here?
GnR_outputs_changed = ((rad_input-native_rad)/native_rad).^2;
rad_function = sqrt((1/length(rad_input))*(sum(GnR_outputs_changed)));

comp_output_changed = ((comp_input-native_comp)/native_comp).^2;
comp_function = sqrt((1/length(comp_input))*(sum(comp_output_changed)));

infl_output_changed = sum((infl_input/max(exp_infl)).^2);
infl_function = sqrt((1/length(infl_input))*(infl_output_changed));

Rad = rad_function/(0.201615159212792*3);
Comp = comp_function/(0.823762611015949*3);
Infl = infl_function/(0.379384971295092*3);

objective_function = (Rad + Comp + Infl);

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