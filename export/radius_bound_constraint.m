% constraint function that caps the rate of change and the flat radius

%return a negative value if the constraint isn't violated, otherwise return a positive number equal to how many points violated the constraint
function output = radius_bound_constraint(input)
    %Read in parameter values to evaluate from inpfile
    B = importdata('Native_in_');
    native_value = B.data(1,1);
    
    %Calculate the objective function value for that parameter set
    
    if length(input) == 1
        cons = input;
    else
        cons = -0.1;
        norm = input;
        for i = 1:length(norm)
            if norm(i) > native_value * 1.1 | norm(i) < native_value * 0.9
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
    