function output = sphere(x)
%Read in parameter values to evaluate from inpfile
A = importdata(x);
%Calculate the objective function value for that parameter set
sse = (A(1) - 5)^2 + (A(2) - 5)^2 + (A(3) - 5)^2;

%Write the objective function value to file Jfile
% Exporting data to .dat
% Creating file to be written to
fileName = fopen('sse.txt','w');
%Writing data to file
fprintf(fileName, '%f %f\n', sse);
%Creating output
output = sse;
%Closing file
fclose(fileName);
