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

xsize = 20;
dotsize = 45;
linewidth = 5;

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

for i = 1:5
    figure;
    hold on;
    fig = gcf;
    fig.Position(3:4) = [800, 800];
    fig.PaperSize = fig.Position(3:4)./72;
    set(gca, 'FontSize', 28, 'FontWeight', 'bold')
    ylim([0 2]);
    yticks([0 0.4 0.8 1.2 1.6 2])
    for j = 1:9
        if results(i,j,4) == 2
            plot(param(i,j), 0, 'X', 'Color', 'Black', 'MarkerSize', xsize, 'LineWidth', linewidth);
        elseif constr(i,j,4) > 0
            plot(param(i,j), results(i,j,4), '.', 'Color', 'Red', 'MarkerSize', dotsize);
        else
            plot(param(i,j), results(i,j,4), '.', 'Color', 'Black', 'MarkerSize', dotsize);
        end
            
    end
    plot(param(i,5), results(i,5,4), '.', 'Color', 'Blue', 'MarkerSize', dotsize);

    currTicks = xticks;
    if length(currTicks) > 6
        newTicks = currTicks(1:2:end);
        newTicks(length(newTicks)+1) = newTicks(end) + (newTicks(end) - newTicks(end-1));
        xlim([0 newTicks(length(newTicks))]);
        xticks(newTicks);
    end
    if i == 1
        xlabel("Volume Fraction");
        ylabel("Triple Objective");
        %title("Volume Fraction Sensitivity");
        saveas(gcf,"vfrac_total_sens.pdf");
    elseif i == 2
        xlabel("Fiber Diameter (µm)");
        ylabel("Triple Objective");
        %title("Fiber Diameter Sensitivity");
        saveas(gcf,"fdiam_total_sens.pdf");
    elseif i == 3
        xlabel("Degradation Rate (1/days)");
        ylabel("Triple Objective");
        %title("Degradation Rate Sensitivity");
        saveas(gcf,"rate_total_sens.pdf");
    elseif i == 4
        xlabel("Degradation Shape (days)");
        ylabel("Triple Objective");
        %title("Degradation Shape Sensitivity");
        saveas(gcf,"shape_total_sens.pdf");
    elseif i == 5
        xlabel("Shear Modulus (Pa)");
        ylabel("Triple Objective");
        %title("Shear Modulus Sensitivity");
        saveas(gcf,"shear_total_sens.pdf");
    end
end

system("mkdir -p sensitivity_plots");
system("mv *.pdf sensitivity_plots");


