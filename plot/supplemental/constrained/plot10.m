%Plotting all 10 trials - radius, compliance, thickness, inflammation - for
% >100% (floor), 90-110% (range), 100-110% (ceiling).  

clear; clc;
%Importing Experimental Data
gnr_exp = importdata('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/master/Gnr_out_experimental_');
rad_exp = gnr_exp(:,1);
infl_exp = gnr_exp(:,6) + gnr_exp(:,7);
rad_exp = rad_exp*1000;
th_exp = gnr_exp(:,2)*1000;

expout_exp = importdata('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/master/Exp_out_experimental_1day');
comp_exp = expout_exp(:,1);



%ADJUST LINEWIDTH AND COLOR HERE
linewidth = 2;
cmap = colormap(turbo(10));

%Loop through 'floor', 'range', 'ceiling
for k = 1:3
    switch k
        case 1
            f = '+minrad/';
            name = 'floor';
        case 2
            f = '+minrad+radcap/';
            name = 'ceiling';
        case 3
            f = '+rad_range/';
            name = 'range';
    end

    %Loop through radius, compliance, inflammation, thickness
    for i = 1:4
        fig = figure;
        hold on;
        fig.Position(3:4) = [800,800];
        fig.PaperSize = fig.Position(3:4)./72;
        xlabel("Time (days)");
        %FONTS HERE
        set(gca, 'FontSize', 20, 'FontWeight', 'bold', 'LineWidth', linewidth*.7);
    
        %Loop through trials
        for j = 1:10
            color = cmap(j,:);
            infile = importdata(strcat('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/fig4/inflam', f, num2str(j), '/results/GnR_out_'));
            
            switch i
                case 1
                    yvar = "rad";
                    ylabel("Radius (mm)");
                    xvar = infile(:,1).*1000; %radius
                case 2
                    yvar = "comp";
                    ylabel("Compliance");
                    infile = importdata(strcat('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/fig4/inflam', f, num2str(j), '/results/Exp_out_'));
                    xvar = infile(:,1); %compliance
                case 3
                    yvar = "th";
                    ylabel("Thickness (mm)");
                    xvar = infile(:,2).*1000; %thickness
                case 4
                    yvar = "infl";
                    ylabel("Inflammation (kg/m^3)");
                    xvar = infile(:,6) + infile(:,7); %inflammation
            end
            if i == 2
                plot((1:length(xvar)).*5-5, xvar, 'DisplayName', int2str(j), 'color', color, 'LineWidth', linewidth);
            else
                plot(1:length(xvar), xvar, 'DisplayName', int2str(j), 'color', color, 'LineWidth', linewidth);
            end
        end
            

        switch i
            case 1
                plot(1:361, rad_exp, 'DisplayName', 'Exp', 'color', '#888', 'LineWidth', linewidth);
                plot(1:361, 8.573 + zeros(361, 1), 'DisplayName', 'Native', 'color', 'black', 'LineWidth', linewidth);
                [~, hobj, ~, ~] = legend('Location','southeast');
                set(hobj, 'LineWidth', 2.5);
                %Plotting constraint
                switch k
                    case 1
                        plot(1:361, 8.573*1.1 + zeros(361, 1),'HandleVisibility','off', 'color', 'red', 'LineWidth', 1, 'LineStyle','--');
                        plot(1:361, 8.573*0.9 + zeros(361, 1),'HandleVisibility','off', 'color', 'red', 'LineWidth', 1, 'LineStyle','--');
                    case 2
                        plot(1:361, 8.573*1.1 + zeros(361, 1),'HandleVisibility','off', 'color', 'black', 'LineWidth', 1, 'LineStyle','--');
                        plot(1:361, 8.573*0.9 + zeros(361, 1),'HandleVisibility','off', 'color', 'red', 'LineWidth', 1, 'LineStyle','--');
                    case 3
                        plot(1:361, 8.573*1.1 + zeros(361, 1),'HandleVisibility','off', 'color', 'black', 'LineWidth', 1, 'LineStyle','--');
                        plot(1:361, 8.573*0.9 + zeros(361, 1),'HandleVisibility','off', 'color', 'black', 'LineWidth', 1, 'LineStyle','--');
                end
                ylim([5.5 10]);
            case 2
                plot(1:361, comp_exp, 'DisplayName', 'Exp', 'color', '#888', 'LineWidth', linewidth);
                plot(1:361, 0.064 + zeros(1,361), 'DisplayName', 'Native', 'color', 'black', 'LineWidth', linewidth);
                [~, hobj, ~, ~] = legend('Location','northwest');
                set(hobj, 'LineWidth', 2.5);
            case 3
                plot(1:361, th_exp, 'DisplayName', 'Exp', 'color', '#888', 'LineWidth', linewidth);
                plot(1:361, 0.7 + zeros(1,361), 'DisplayName', 'Native', 'color', 'black', 'LineWidth', linewidth);
                [~, hobj, ~, ~] = legend('Location','northeast');
                set(hobj, 'LineWidth', 2.5);
            case 4
                plot(1:361, infl_exp, 'DisplayName', 'Exp', 'color', '#888', 'LineWidth', linewidth);
                plot(1:361, zeros(1,361), 'DisplayName', 'Native', 'color', 'black', 'LineWidth', linewidth);
                [~, hobj, ~, ~] = legend('Location','northeast');
                set(hobj, 'LineWidth', 2.5);
        end
        
        %Axes
        xlim([0 361]);
        currTicks = yticks;
        if length(currTicks) > 6
            newTicks = currTicks(1:2:end);
            newTicks(length(newTicks)+1) = newTicks(end) + (newTicks(end) - newTicks(end-1));
            yticks(newTicks);
        end
        saveas(gcf, strcat(yvar, '_', name, ".pdf"));
    end   
end
system("mkdir 10_results");
system("mv -f *.pdf 10_results");