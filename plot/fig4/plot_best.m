%Matlab script to plot radius, thickness, compliance, and inflammation for best trial.  run the daily compliance GNR script before running this to get better compliance plotting.

clear; clc;

opt1 = 10;
opt2 = 7;
opt3 = 1;

obj1 = "Floor";
obj2 = "Ceiling";
obj3 = "Range";


%importing
out1 = importdata(strcat('/Users/ryanhsu/Documents/GitHub/SMF/GnR/final_paper/fig4/floor/', string(opt1), '/results/GnR_out_'));
expout1 = importdata(strcat('/Users/ryanhsu/Documents/GitHub/SMF/GnR/final_paper/fig4/floor/', string(opt1), '/Exp_out_'));

out2 = importdata(strcat('/Users/ryanhsu/Documents/GitHub/SMF/GnR/final_paper/fig4/ceiling/', string(opt2), '/results/GnR_out_'));
expout2 = importdata(strcat('/Users/ryanhsu/Documents/GitHub/SMF/GnR/final_paper/fig4/ceiling/', string(opt2), '/Exp_out_'));

out3 = importdata(strcat('/Users/ryanhsu/Documents/GitHub/SMF/GnR/final_paper/fig4/range/', string(opt3), '/results/GnR_out_'));
expout3 = importdata(strcat('/Users/ryanhsu/Documents/GitHub/SMF/GnR/final_paper/fig4/range/', string(opt3), '/Exp_out_'));


rad1 = out1(:,1)*1000;
comp1 = expout1(:,1);
infl1 = out1(:,6) + out1(:,7);
th1 = out1(:,2)*1000;

rad2 = out2(:,1)*1000;
comp2 = expout2(:,1);
infl2 = out2(:,6) + out2(:,7);
th2 = out2(:,2)*1000;

rad3 = out3(:,1)*1000;
comp3 = expout3(:,1);
infl3 = out3(:,6) + out3(:,7);
th3 = out3(:,2)*1000;

gnr_exp = importdata('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/master/Gnr_out_experimental_');
rad_exp = gnr_exp(:,1);
infl_exp = gnr_exp(:,6) + gnr_exp(:,7);
rad_exp = rad_exp*1000;
th_exp = gnr_exp(:,2)*1000;

expout_exp = importdata('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/master/Exp_out_experimental_1day');
comp_exp = expout_exp(:,1);

%HERE
linewidth = 5;

for i = 1:4
    fig = figure;
    hold on;
    fig.Position(3:4) = [800,800];
    fig.PaperSize = fig.Position(3:4)./72;
    xlabel("Time (days)");

    switch i
        case 1
            yvar = "rad";
            ylabel("Radius (mm)");
        case 2
            yvar = "comp";
            ylabel("Compliance");
        case 3
            yvar = "th";
            ylabel("Thickness (mm)");
        case 4
            yvar = "infl";
            ylabel("Inflammation (kg/m^3)");
    end
    for j = 1:3
        switch j
            case 1
                style = '--';
            case 2
                style = '-.';
            case 3
                style = ':';
        end
        xvar = eval(strcat(yvar,int2str(j)));
        plot(1:length(xvar), xvar, 'DisplayName', eval(strcat("obj",int2str(j))), 'linestyle', style, 'color', 'black', 'LineWidth', linewidth);
    end
    %HERE
    set(gca, 'FontSize', 35, 'FontWeight', 'bold', 'LineWidth', linewidth*.7);
    switch i
        case 1
            plot(1:361, rad_exp, 'DisplayName', 'Exp', 'color', '#888', 'LineWidth', linewidth);
            plot(1:361, 8.573 + zeros(361, 1), 'DisplayName', 'Native', 'color', 'black', 'LineWidth', linewidth);
            plot(1:361, 8.573*0.9 + zeros(361, 1), 'color', 'black', 'LineWidth', 1, 'LineStyle','--', 'HandleVisibility','off');
            plot(1:361, 8.573*1.1 + zeros(361, 1), 'color', 'black', 'LineWidth', 1, 'LineStyle','--', 'HandleVisibility', 'off');
            [~, hobj, ~, ~] = legend('Location','southeast');
            set(hobj, 'LineWidth', 2.5);
        case 2
            plot(1:361, comp_exp, 'color', '#888', 'LineWidth', linewidth, 'DisplayName', 'Exp');
            plot(1:361, 0.064 + zeros(1,361), 'color', 'black', 'LineWidth', linewidth, 'DisplayName', 'Native');
            [~, hobj, ~, ~] = legend('Location','northeast');
            set(hobj, 'LineWidth', 2.5);
        case 3
            plot(1:361, th_exp, 'color', '#888', 'LineWidth', linewidth, 'DisplayName', 'Exp');
            plot(1:361, 0.7 + zeros(1,361), 'color', 'black', 'LineWidth', linewidth, 'DisplayName', 'Native');
            [~, hobj, ~, ~] = legend('Location','northeast');
            set(hobj, 'LineWidth', 2.5);
        case 4
            plot(1:361, infl_exp, 'color', '#888', 'LineWidth', linewidth, 'DisplayName', 'Exp');
            plot(1:361, zeros(1,361), 'color', 'black', 'LineWidth', linewidth, 'DisplayName', 'Native');
            [~, hobj, ~, ~] = legend('Location','northeast');
            set(hobj, 'LineWidth', 2.5);
    end
    
    xlim([0 361]);
    currTicks = yticks;
    if length(currTicks) > 6
        newTicks = currTicks(1:2:end);
        newTicks(length(newTicks)+1) = newTicks(end) + (newTicks(end) - newTicks(end-1));
        yticks(newTicks);
    end
    saveas(gcf, strcat(yvar,".pdf"));
end   

system("mkdir best_results");
system("mv *.pdf best_results");
