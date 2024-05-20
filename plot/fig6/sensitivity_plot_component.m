%definitions

vfrac_up_bound = 0.3;
vfrac_low_bound = 0.025;

fdiam_up_bound = 20;
fdiam_low_bound = 2.0;

k1_up_bound = 0.1;
k1_low_bound = 0.003;

zeta1_up_bound = 392.0;
zeta1_low_bound = 50.0;

cp1_1_up_bound = 3200000000.0;
cp1_1_low_bound = 450000000.0;

%fetch optimized values

dirnum = 7;
dirstr = num2str(dirnum);
inp_path = fullfile('/Users/ryanhsu/Documents/GitHub/SMF/GnR/final_paper/fig4/ceiling', dirstr, 'results', 'filter.txt');
fid = fopen(inp_path, 'r');
line = fgetl(fid);
fclose(fid);
sline = split(line);

optvals = [str2double(sline(1)), str2double(sline(2)), str2double(sline(3)), str2double(sline(4)), str2double(sline(5))];
optvals(5) = (exp((log(cp1_1_up_bound) - log(cp1_1_low_bound)) * optvals(5) + log(cp1_1_low_bound)) - cp1_1_low_bound) / (cp1_1_up_bound - cp1_1_low_bound);


results = zeros(5,9,4);
constr = zeros(5,9,4);
%un-normalize parameters
param = zeros(5,9);
for i = 1:5
    param(i,:) = [0, 0.05, 0.1, max(optvals(i) - 0.1, 0), optvals(i), min(optvals(i) + 0.1, 1), 0.9, 0.95, 1];
end
factors = [vfrac_up_bound - vfrac_low_bound, fdiam_up_bound - fdiam_low_bound, k1_up_bound - k1_low_bound, zeta1_up_bound - zeta1_low_bound, cp1_1_up_bound - cp1_1_low_bound];
offset = [vfrac_low_bound, fdiam_low_bound, k1_low_bound, zeta1_low_bound, cp1_1_low_bound];

for i = 1:9
    param(:,i) = reshape((reshape(param(:,i), [1,5]) .* factors + offset), [5,1]);
end

%fetch sensitivity results and plot
for i = 1:5
    for j = 1:9
    obj = importdata(strcat(pwd, '/sensitivity/', num2str(i), '_', num2str(j), '/sse.txt'));
    cons = importdata(strcat(pwd, '/sensitivity/', num2str(i), '_', num2str(j), '/cons.txt'));
    results(i,j,:) = obj;
    constr(i,j,:) = cons;
    end
end

xsize = 20;
dotsize = 45;
linewidth = 5;

for i = 1:5
    for k = 1:3
        fig = figure;
        fig.Position(3:4) = [800,800];
        fig.PaperSize = fig.Position(3:4)./72;
        set(gca, 'FontSize', 35, 'FontWeight', 'bold');
        ylim([0 1]);
        yticks([0 0.2 0.4 0.6 0.8 1])
        hold on;
        for j = 1:9
            if results(i,j,1) == 2
                plot(param(i,j), 0, 'X', 'Color', 'Red', 'MarkerSize', xsize, 'LineWidth', linewidth);
            elseif constr(i,j,1) > 0
                plot(param(i,j), results(i,j,k), '.', 'Color', 'Red', 'MarkerSize', dotsize);
            elseif j ~= 5
                plot(param(i,j), results(i,j,k), '.', 'Color', 'Black', 'MarkerSize', dotsize);
            end
        end
        plot(param(i,5), results(i,5,k), '.', 'Color', 'Blue', 'MarkerSize', dotsize);

        name = "";


        switch k
        case 1
            ylabel("Radius Component");
            name = strcat(name, "rad");
        case 2
            ylabel("Compliance Component");
            name = strcat(name, "comp");
        case 3
            ylabel("Inflammation Component");
            name = strcat(name, "infl");
        end
        currTicks = xticks;
        if length(currTicks) > 6
            newTicks = currTicks(1:2:end);
            newTicks(length(newTicks)+1) = newTicks(end) + (newTicks(end) - newTicks(end-1));
            xlim([0 newTicks(length(newTicks))]);
            xticks(newTicks);
        end
        
        if i == 1
            xlabel("Volume Fraction");
            name = strcat(name, "_vfrac.pdf");
            %sgtitle("Volume Fraction Sensitivity");
        elseif i == 2
            xlabel("Fiber Diameter (µm)");
            name = strcat(name, "_fdiam.pdf");
            %sgtitle("Fiber Diameter Sensitivity");
        elseif i == 3
            xlabel("Degradation Rate (1/days)");
            name = strcat(name, "_rate.pdf");
            %sgtitle("Degradation Rate Sensitivity");
        elseif i == 4
            xlabel("Degradation Shape (days)");
            name = strcat(name, "_shape.pdf");
            %sgtitle("Degradation Shape Sensitivity");
        elseif i == 5
            xlabel("Shear Modulus (Pa)");
            name = strcat(name, "_modulus.pdf");
            %sgtitle("Shear Modulus Sensitivity");
        end
        saveas(gcf,name);
    end
end

system("mkdir -p sensitivity_plots");
system("mv *.pdf sensitivity_plots");
% 
% for i = 1:5
%     figure;
%     subplot(1,3,1);
%     hold on;
%     ylim([0 1]);
%     fig = gcf;
%     fig.Position(3:4) = [1200, 500];
% 
%     subplot(1,3,2);
%     hold on;
%     ylim([0 1]);
%     fig = gcf;
%     fig.Position(3:4) = [1200, 500];
%     
%     subplot(1,3,3);
%     hold on;
%     ylim([0 1]);
%     fig = gcf;
%     fig.Position(3:4) = [1200, 500];
% 
%     for j = 1:9
%         if results(i,j,1) == 2
%             subplot(1,3,1);
%             plot(param(i,j), 0, 'X', 'Color', 'Red', 'MarkerSize', 10);
%             subplot(1,3,2);
%             plot(param(i,j), 0, 'X', 'Color', 'Red', 'MarkerSize', 10);
%             subplot(1,3,3);
%             plot(param(i,j), 0, 'X', 'Color', 'Red', 'MarkerSize', 10);
%         elseif constr(i,j,1) > 0
%             subplot(1,3,1);
%             plot(param(i,j), results(i,j,1), 'o', 'Color', 'Red', 'MarkerSize', 5);
%             subplot(1,3,2);
%             plot(param(i,j), results(i,j,2), 'o', 'Color', 'Red', 'MarkerSize', 5);
%             subplot(1,3,3);
%             plot(param(i,j), results(i,j,3), 'o', 'Color', 'Red', 'MarkerSize', 5);
%         else
%             subplot(1,3,1);
%             plot(param(i,j), results(i,j,1), '.', 'Color', 'Black', 'MarkerSize', 10);
%             subplot(1,3,2);
%             plot(param(i,j), results(i,j,2), '.', 'Color', 'Black', 'MarkerSize', 10);
%             subplot(1,3,3);
%             plot(param(i,j), results(i,j,3), '.', 'Color', 'Black', 'MarkerSize', 10);
%         end
%     end
%     for k = 1:3
%         subplot(1,3,k);
%         plot(param(i,5), results(i,5,k), '.', 'Color', 'Blue', 'MarkerSize', 10);
%         switch k
%         case 1
%             ylabel("Radius Component of Objective");
%         case 2
%             ylabel("Compliance Component of Objective");
%         case 3
%             ylabel("Inflammation Component of Objective");
%         end
%         if i == 1
%             xlabel("Volume Fraction");
%             title("Volume Fraction Sensitivity");
%             saveas(gcf,"vfrac_part_sens.pdf");
%         elseif i == 2
%             xlabel("Fiber Diameter (µm)");
%             title("Fiber Diameter Sensitivity");
%             saveas(gcf,"fdiam_part_sens.pdf");
%         elseif i == 3
%             xlabel("Degradation Rate (1/days)");
%             title("Degradation Rate Sensitivity");
%             saveas(gcf,"rate_part_sens.pdf");
%         elseif i == 4
%             xlabel("Degradation Shape (days)");
%             title("Degradation Shape Sensitivity");
%             saveas(gcf,"shape_part_sens.pdf");
%         elseif i == 5
%             xlabel("Shear Modulus (Pa)");
%             title("Shear Modulus Sensitivity");
%             saveas(gcf,"shear_part_sens.pdf");
%         end
%     end
% end
% 
% system("mkdir -p sensitivity_plots");
% system("mv vfrac_part_sens.pdf fdiam_part_sens.pdf rate_part_sens.pdf shape_part_sens.pdf shear_part_sens.pdf sensitivity_plots");