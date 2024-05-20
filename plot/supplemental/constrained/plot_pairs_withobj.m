%Plotting all 10 trials - pairwise parameter plots - for
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

%Importing Trial Data
values = zeros(10, 6, 3);

for k = 1:3
    for i = 1:10
        switch k
            case 1
                f = '+minrad/';
            case 2
                f = '+minrad+radcap/';
            case 3
                f = '+rad_range/';
        end    
    filt = importdata(strcat('/Users/ryanhsu/Documents/GitHub/SMF/GnR/paper/fig4/inflam', f, num2str(i), '/results/filter.txt'));
    filt_nums = filt(1,:);
        for j = 1:6
            values(i,j,k) = filt_nums(j);
        end
    end
end

%scaling the parameters to their actual values
values(:,1,:) = values(:,1,:).*(0.3-0.025)+0.025; %volume fraction (%)
values(:,2,:) = values(:,2,:).*(20-2)+2; %fiber diameter (µm)
values(:,3,:) = values(:,3,:).*(0.1-0.003)+0.003; %degradation rate (1/days)
values(:,4,:) = values(:,4,:).*(392-50)+50; %degradation shape (days)
values(:,5,:) = values(:,5,:).*(3200000000-480000000)+480000000; %shear modulus (Pa)


%ADJUST LINEWIDTH AND COLOR HERE
linewidth = 2;
markersize = 50;



%Loop through 'floor', 'range', 'ceiling
for i = 1:3
    fig = figure;
    hold on;
    fig.Position(3:4) = [800,800];
    fig.PaperSize = fig.Position(3:4)./72;
    %FONTS HERE
    set(gca, 'FontSize', 20, 'FontWeight', 'bold', 'LineWidth', linewidth*.7);

    cmap = colormap(hot);
    caxis([min(min(values(:,6,:))), max(max(values(:,6,:)))]);
    colorbar;
    switch i
        case 1
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
                scatter(values(:,1,k), values(:,2,k), markersize*3, values(:,6,k), 'o', "filled");
            end
            xlabel("Volume Fraction (%)");
            ylabel("Fiber Diameter (µm)");
            xlim([0.025, 0.3]);
            ylim([2, 20]);
            title("Fiber Diameter vs. Volume Fraction");
            type = 'fdiam_vfrac';
        case 2
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
                scatter(values(:,1,k), values(:,5,k), markersize*3, values(:,6,k), 'o', "filled");
            end
            xlabel("Volume Fraction (%)");
            ylabel("Shear Modulus (Pa)");
            xlim([0.025, 0.3]);
            ylim([480000000, 3200000000]);
            title("Shear Modulus vs. Volume Fraction");
            type = 'modulus_vfrac';
        case 3
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
                scatter(values(:,3,k), values(:,4,k), markersize*3, values(:,6,k), 'o', "filled");
            end
            xlabel("Degradation Rate (1/days)");
            ylabel("Degradation Shape (days)");
            xlim([0.003, 0.1]);
            ylim([50, 392]);
            title("Degradation Shape vs. Degradation Rate");
            type = 'shape_rate';
    end
    %Axes
    currTicks = yticks;
    if length(currTicks) > 6
        newTicks = currTicks(1:2:end);
        newTicks(length(newTicks)+1) = newTicks(end) + (newTicks(end) - newTicks(end-1));
        yticks(newTicks);
    end

    currTicks = xticks;
    if length(currTicks) > 6
        newTicks = currTicks(1:2:end);
        newTicks(length(newTicks)+1) = newTicks(end) + (newTicks(end) - newTicks(end-1));
        xticks(newTicks);
    end
    saveas(gcf,strcat(type,'.pdf'));
end   

system("mkdir pair_obj_results");
system("mv -f *.pdf pair_obj_results");






