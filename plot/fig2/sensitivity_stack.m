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
factors = [vfrac_up_bound - vfrac_low_bound, fdiam_up_bound - fdiam_low_bound, k1_up_bound - k1_low_bound, zeta1_up_bound - zeta1_low_bound, cp1_1_up_bound - cp1_1_low_bound];
offset = [vfrac_low_bound, fdiam_low_bound, k1_low_bound, zeta1_low_bound, cp1_1_low_bound];
optvals = [0.15 8 0.0125 275 2000000000];
optvals = (optvals - offset)./factors;


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
    fig = figure;
    fig.Position(3:4) = [800,800];
    fig.PaperSize = fig.Position(3:4)./72;
    set(gca, 'FontSize', 35, 'FontWeight', 'bold');
    ylim([0 1]);
    yticks([0 0.2 0.4 0.6 0.8 1])
    hold on;
    cmap = colormap(hsv(3));
    for k = 1:3
        switch k
            case 1
                name = "Radius";
            case 2
                name = "Compliance";
            case 3
                name = "Inflammation";
        end
        for j = 1:9
            if results(i,j,1) == 2
                plot(param(i,j), 0, 'X', 'Color', 'Black', 'MarkerSize', xsize, 'LineWidth', linewidth, 'HandleVisibility','off');
            elseif j ~= 5
                plot(param(i,j), results(i,j,k), '.', 'Color', cmap(k,:), 'MarkerSize', dotsize, 'HandleVisibility','off');
            end
        end
        plot(param(i,5), results(i,5,k), '.', 'Color', 'k', 'MarkerSize', dotsize,'HandleVisibility','off');
        plot(param(i,5), results(i,5,k), '.', 'Color', cmap(k,:), 'MarkerSize', dotsize*.75,'DisplayName', name);


        switch k
        case 1
            name = strcat(name, "rad");
        case 2
            name = strcat(name, "comp");
        case 3 
            name = strcat(name, "infl");
        end
        ylabel("Objective Value");
        currTicks = xticks;
        if length(currTicks) > 6
            newTicks = currTicks(1:2:end);
            newTicks(length(newTicks)+1) = newTicks(end) + (newTicks(end) - newTicks(end-1));
            xlim([0 newTicks(length(newTicks))]);
            xticks(newTicks);
        end
        [~, hobj, ~, ~] = legend('Location','northeast');
        set(hobj, 'LineWidth', 2.5);
        if i == 1
            xlabel("Volume Fraction");
            name = "vfrac.pdf";
        elseif i == 2
            xlabel("Fiber Diameter (µm)");
            name = "fdiam.pdf";
        elseif i == 3
            xlabel("Degradation Rate (1/days)");
            name = "rate.pdf";
        elseif i == 4
            xlabel("Degradation Shape (days)");
            name = "shape.pdf" ;
        elseif i == 5
            xlabel("Shear Modulus (Pa)");
            name = "modulus.pdf";
        end
        saveas(gcf,name);
    end
end

system("mkdir -p sensitivity_stack_plots");
system("mv *.pdf sensitivity_stack_plots");