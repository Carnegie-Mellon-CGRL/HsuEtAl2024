%Matlab script to plot radius, thickness, compliance, and inflammation
%Ensure compliance has been calculated, or comment out the compliance graphing section

clear; clc;
curr = pwd;

len = 361;
outfile = fopen('output.txt', 'w');
objective = ' for 5 variable optimization with triple objective and radius between 100% and 110% of native';
objs = "J(R)";
%Output:
trials = 10;

for i = 1:trials
    out = importdata(strcat(pwd, '/', string(i), '/results/GnR_out_'));
    expout = importdata(strcat(pwd, '/', string(i), '/results/Exp_out_'));
    filt = readlines(strcat(pwd, '/', string(i), '/results/filter.txt'));
    fprintf(outfile, "Trial %d: %s\n", i, (filt(1)));
    
    rad = out(:,1);
    comp = expout(:,1);
    infl = out(:,6) + out(:,7);
    rad = rad*1000;
    opt_obj = triple_objective(rad, comp, infl);
    fprintf(outfile, "Objective: %f\n", opt_obj);
    system("rm sse.txt");
end

gnr_exp = importdata('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/master/Gnr_out_experimental_');
rad_exp = gnr_exp(:,1);
infl_exp = gnr_exp(:,6) + gnr_exp(:,7);
rad_exp = rad_exp*1000;

expout_exp = importdata('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/master/Exp_out_experimental_1day');
comp_exp = expout_exp(:,1);

exp_obj = triple_objective(rad_exp, comp_exp, infl_exp);
fprintf(outfile, "Unoptimized objective: %.6f\n", exp_obj);
system("rm sse.txt");


%Graphing radius:

figure;
hold on;
xlabel("Time (days)");
ylabel("Radius (mm)");
title(strcat("Radius", objective));

    %plotting the trials
for i = 1:trials
    out = importdata(strcat(pwd, '/', string(i), '/results/GnR_out_'));
    rad = out(:,1);
    rad = rad*1000;
    plot(1:length(rad), rad, 'DisplayName', strcat(objs, string(i)), 'linewidth', 1.25);
end


    %plotting experimental and native
plot(1:length(rad_exp), rad_exp, 'DisplayName', 'Exp', 'color', '#888', 'linewidth', 1.25);
plot(1:length(rad_exp), 8.573 + zeros(length(rad_exp), 1), 'DisplayName', 'native', 'color', 'black', 'linewidth', 1.25);
legend('Location', 'southeast');
saveas(gcf, 'radius.pdf');


%plot compliance
figure;
hold on;
xlabel("Time (days)");
ylabel("Compliance");
title(strcat("Compliance", objective));

for i = 1:trials
    expout = importdata(strcat(pwd, '/', string(i), '/Exp_out_'));
    comp = expout(:,1);
    plot((0:(length(comp))-1).*5, comp, 'DisplayName', strcat(objs, string(i)), 'linewidth', 1.25);
end


plot((1:length(comp_exp)).*5, comp_exp, 'DisplayName', 'Exp', 'color', '#888', 'linewidth', 1.25);
plot((1:length(comp_exp)).*5, 0.06 + zeros(length(comp_exp), 1), 'DisplayName', 'Native', 'color', 'black', 'linewidth', 1.25);
legend('Location', 'southeast');
saveas(gcf, 'compliance.pdf');

%plot thickness
figure;
hold on;
xlabel("Time (days)");
ylabel("Thickness (mm)");
title(strcat("Thickness", objective));

for i = 1:trials
    out = importdata(strcat(pwd, '/', string(i), '/results/GnR_out_'));
    th = out(:,2);
    th = th*1000;
    plot(1:length(th), th, 'DisplayName', strcat(objs, string(i)), 'linewidth', 1.25);
end
legend;

th_exp = gnr_exp(:,2);
th_exp = th_exp*1000;
plot(1:length(th_exp), th_exp, 'DisplayName', 'Exp', 'color', '#888', 'linewidth', 1.25);
plot(1:length(th_exp), 0.7 + zeros(length(th_exp), 1), 'DisplayName', 'native', 'color', 'black', 'linewidth', 1.25);
saveas(gcf, 'thickness.pdf');

%plot inflammation


figure;
hold on;
xlabel("Time (days)");
ylabel("Inflammation (kg/m^3)");
title(strcat("Inflammation", objective));

for i = 1:trials
    out = importdata(strcat(pwd, '/', string(i), '/results/GnR_out_'));
    inflam1 = out(:,6);
    inflam2 = out(:,7);
    plot(1:length(inflam1), inflam1+inflam2, 'DisplayName', strcat(objs, string(i)), 'linewidth', 1.25);
end
legend;
infl1 = gnr_exp(:,6);
infl2 = gnr_exp(:,7);
plot(1:length(infl1), infl1+infl2, 'DisplayName', 'Exp', 'color', '#888', 'linewidth', 1.25);

saveas(gcf, 'inflammation.pdf');

system("mkdir results");
system("mv radius.pdf compliance.pdf thickness.pdf inflammation.pdf results");