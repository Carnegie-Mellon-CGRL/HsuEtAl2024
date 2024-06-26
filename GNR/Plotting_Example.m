clc; clear all; clf;

% colors = parula(5);

% colors = bone(6);

colors = [
          0.0 0.0 0.0
          1.0 0.0 0.0
          0.0 0.0 1.0
          ];
      
styles = ['-', ":", ":"];
      

data_files = {
              'GnR_out_'
              };
          
%Can plot multiple files with outputs from different runs (uncomment below)
%Each run overwrites GnR_out_ unless it is given a non-default name
% data_files = {
%               'GnR_out_',...
%               'GnR_out_PGA',...
%               'GnR_out_PGA'
%               };


%Parameterized inputs
n_par = length(data_files);
time = 0:360;

for i = 1:n_par
    
    data = load(data_files{i});
    
    a = data(:,1);
    h = data(:,2);
    rhoR = data(:,3);
    rhoR_p1 = data(:,4);
    rhoR_p2 = data(:,5);
    rhoR_i_m = data(:,6);
    rhoR_i_c = data(:,7);

    figure(1)
    subplot(2,2,1)
    hold on
    plot(time / 7, a/a(1), styles(i), 'LineWidth', 2.0, 'Color', colors(i,:))
    xlabel('Time (weeks)'); ylabel('Inner Radius (-)')
    xlim([0 53])
    set(gca, 'FontSize', 14, 'LineWidth', 2.0, 'FontWeight', 'bold')

    subplot(2,2,2)
    hold on
    plot(time / 7, h/h(1), styles(i), 'LineWidth', 2.0, 'Color', colors(i,:))
    xlabel('Time (weeks)'); ylabel('Thickness (-)')
    xlim([0 53])
    set(gca, 'FontSize', 14, 'LineWidth', 2.0, 'FontWeight', 'bold')

    subplot(2,2,3)
    hold on
    plot(time / 7, rhoR_p1, styles(i), 'LineWidth', 2.0, 'Color', colors(i,:))
    xlabel('Time (weeks)'); ylabel('rho^{p}_{sp} (kg/m^3)')
    xlim([0 53])
    set(gca, 'FontSize', 14, 'LineWidth', 2.0, 'FontWeight', 'bold')

    subplot(2,2,4)
    hold on
    plot(time / 7, rhoR_p2, styles(i), 'LineWidth', 2.0, 'Color', colors(i,:))
    xlabel('Time (weeks)'); ylabel('rho^{p}_{kn} (kg/m^3)')
    xlim([0 53])
    set(gca, 'FontSize', 14, 'LineWidth', 2.0, 'FontWeight', 'bold')

end


