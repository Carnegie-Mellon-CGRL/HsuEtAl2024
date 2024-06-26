clear all; clc;

%ouput directory
out_dir = 'sphere_test';

%code that copies smf.inp to a new file with a different name,
%e.g smf_clone.inp
copyfile smf.inp smf_clone.inp

%code that deletes these files if they already exist
files = {'xnew.dat' 'temp.dat' 'c.txt' 'sse.txt' 'cons.txt' 'Jnew.dat' 'Hnew.dat' 'stopfile' 'Jhist.dat' 'Hhist.dat' 'BestJHist.dat' 'x_init_lhs.dat' 'xhist.dat' 'filter.txt' 'theta_hist'};
for k = 1 : length(files)
    if exist(files{k}, 'file')==2
        delete(files{k});
    end
end
ptsfile="xnew.dat";
temppts="temp.dat";
inpfile="c.txt";
outfile="sse.txt";
consfile="cons.txt";
Jfile="Jnew.dat";
Hfile="Hnew.dat";
endfile="stopfile";

while (~isfile(endfile))

    %Run the current optimization step
    command = "./krig";
    system(command);
    
    if (isfile(Jfile))
        command = strcat("rm ", Jfile);
        system(command);
    end

    if (isfile(Hfile))
        command = strcat("rm ", Hfile);
        system(command);
    end

    command = strcat("cp ", ptsfile, " ", temppts);
    system(command);

    pts_to_eval = load(temppts);
    
    %Run the objective and constraint function for each trial pt
    for i = 1:length(pts_to_eval(:,1))
        
        fid = fopen(inpfile, 'w');
        fprintf(fid, "%f ", pts_to_eval(i,:));
        fclose(fid);
        
        %To run Matlab objective function
        output = sphere(inpfile);
        
        %Run without constraint, i.e. every point is a feasible
        %point since H is written to 0
        fid = fopen(consfile, 'w');
        fprintf(fid, "%f ", 0.0);
        fclose(fid);

        if (isfile(Jfile))
          command = strcat("cat ", Jfile, " ", outfile, " > ", "1_", Jfile);
          system(command);
          command = strcat("mv ", "1_", Jfile, " ", Jfile);
          system(command);
        else
          command = strcat("cp ", outfile, " ", Jfile);
          system(command);
        end

        if (isfile(Hfile))
          command = strcat("cat ", Hfile, " ", consfile, " > ", "1_", Hfile);
          system(command);
          command = strcat("mv ", "1_", Hfile, " ", Hfile);
          system(command);
        else
          command = strcat("cp ", consfile, " ", Hfile);
          system(command);
        end

    end

command = strcat("rm ", temppts);
system(command);

end

%copies smf.inp to a new file with a different name, e.g smf_end.inp
copyfile smf.inp smf_end.inp

%code that renames smf_clone.inp to smf.inp
movefile('smf_clone.inp','smf.inp');

%code that saves all files created in a new folder
%Remember to update folder name before every run
mkdir(out_dir);
files_end = {'c.txt' 'sse.txt' 'cons.txt' 'Jnew.dat' 'Hnew.dat' 'stopfile' 'Jhist.dat' 'Hhist.dat' 'BestJHist.dat' 'x_init_lhs.dat' 'xhist.dat' 'filter.txt' 'theta_hist' 'smf_end.inp'};
for k = 1 : length(files_end)
    movefile(files_end{k}, out_dir);
end