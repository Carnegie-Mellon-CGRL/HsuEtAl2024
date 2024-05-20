clear; clc;
%This script runs a full optimization.  The only modification that needs to be made is for the objective and constraint selection on line ~100
%Run this using the bash script "master" to run multiple trials at once

%make a default structure for Scaffold_in
Scaffold_default = read_Scaffold_in("Scaffold_in_default");
Scaffold_temp = Scaffold_default;


%deletes these files if they already exist
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

%until the stopfile exists, loop
outdir = "results";
%copy smf input
dim = 5;
system(strcat("python3 build_smf.py ", num2str(dim)));
%output directory
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
    


%Adding in scaling factors for 0 to 1 SMF bounds
% changed bounds for both polymers to enclose initial values for both polymers
vfrac1_up_bound = 0.3;
vfrac1_low_bound = 0.025;

%Adding in scaling factors for fiber diameter
fdiam1_up_bound = 20; %originally 6.416; raising to 8 - 16 for knitted diameter
fdiam1_low_bound = 2;

%Adding in scaling factors for k1
k1_up_bound = 0.1;
k1_low_bound = 0.003;

%Adding in scaling factors for zeta1
zeta1_up_bound = 392;
zeta1_low_bound = 50;

%Adding in scaling factors for cp1_1
cp1_1_up_bound = 3200000000;
cp1_1_low_bound = 450000000; %lowering bound by 2 orders of magnitude for compliance matching - originally 800 MPa
%Run the objective and constraint function for each trial pt
for i = 1:length(pts_to_eval(:,1))
    
    
    vfrac1_scaled = (vfrac1_up_bound - vfrac1_low_bound) * pts_to_eval(i,1) + vfrac1_low_bound;
    fdiam1_scaled = (fdiam1_up_bound - fdiam1_low_bound) * pts_to_eval(i,2) + fdiam1_low_bound;
    k1_scaled = (k1_up_bound - k1_low_bound) * pts_to_eval(i,3) + k1_low_bound;
    zeta1_scaled = (zeta1_up_bound - zeta1_low_bound) * pts_to_eval(i,4) + zeta1_low_bound;
    cp1_1_scaled = exp((log(cp1_1_up_bound) - log(cp1_1_low_bound)) * pts_to_eval(i,5) + log(cp1_1_low_bound));

    %Write points to eval to appropriate parameter location in the
    %current Scaffold
    Scaffold_temp.epsilon_p1_0 = vfrac1_scaled;
    Scaffold_temp.fd_p_1 = fdiam1_scaled;
    Scaffold_temp.k_p1 = k1_scaled;
    Scaffold_temp.zeta_p1 = zeta1_scaled;
    Scaffold_temp.c1_p1 = cp1_1_scaled;

    %Write new structure to file for input to gnr code
    write_Scaffold_in("Scaffold_in_", Scaffold_temp);
    command = "/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/master/gnr -m 361";
    system(command);

    filename = fopen('sse.txt','w');
    fprintf(filename, '%f\n', 2);
    fclose(filename);
    fid = fopen(consfile, 'w');
    fprintf(fid, "%f ", 1.0);
    fclose(fid);
    %check if gnr_out and exp_out exist
    if (isfile("GnR_out_") && isfile("Exp_out_"))
      GnR_outputs = importdata("GnR_out_");
      if (height(GnR_outputs) == 361)
        GnR_comp = importdata("Exp_out_");
        GnR_radius = GnR_outputs(:,1);
        GnR_radius = GnR_radius*1000; %Converting to mm
        GnR_infl = GnR_outputs(:,6) + GnR_outputs(:,7);
        J = triple_objective(GnR_radius, GnR_comp, GnR_infl);  %change objective here
        C = radius_min_cap_constraint(GnR_radius); %change constraint here
      end

      %we are done with both _out_ files, so wipe them here
      command = "rm GnR_out_";
      system(command);
      command = "rm Exp_out_";
      system(command);
    end
    
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

%rename SMF file to smf_end
system("mv smf.inp smf_end.inp");

%re-call gnr on optimized point to get the correct output files - this should be a copy of the code in the loop above
filt = readlines('filter.txt');
opt_params = str2double(split(filt(1)));

vfrac1_scaled = (vfrac1_up_bound - vfrac1_low_bound) * opt_params(1) + vfrac1_low_bound;
fdiam1_scaled = (fdiam1_up_bound - fdiam1_low_bound) * opt_params(2) + fdiam1_low_bound;
k1_scaled = (k1_up_bound - k1_low_bound) * opt_params(3) + k1_low_bound;
zeta1_scaled = (zeta1_up_bound - zeta1_low_bound) * opt_params(4) + zeta1_low_bound;
cp1_1_scaled = exp((log(cp1_1_up_bound) - log(cp1_1_low_bound)) * opt_params(5) + log(cp1_1_low_bound));



%Write points to eval to appropriate parameter location in the
%current Scaffold
Scaffold_temp.epsilon_p1_0 = vfrac1_scaled;
Scaffold_temp.fd_p_1 = fdiam1_scaled;
Scaffold_temp.k_p1 = k1_scaled;
Scaffold_temp.zeta_p1 = zeta1_scaled;
Scaffold_temp.c1_p1 = cp1_1_scaled;


%Write new structure to file for input to gnr code
write_Scaffold_in("Scaffold_in_", Scaffold_temp);
command = "/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/master/gnr -m 361";
system(command);



%save all files created to the output directory
command = strcat("mkdir ", outdir)
system(command);
files_end = {'c.txt' 'sse.txt' 'cons.txt' 'Jnew.dat' 'Hnew.dat' 'stopfile' 'Jhist.dat' 'Hhist.dat' 'BestJHist.dat' 'x_init_lhs.dat' 'xhist.dat' 'filter.txt' 'theta_hist' 'smf_end.inp' 'xnew.dat' 'Exp_out_' 'GnR_out_' 'Scaffold_in_'};
for k = 1 : length(files_end)
  if exist(files_end{k}, 'file')==2
      movefile(files_end{k}, outdir);
  end
end

