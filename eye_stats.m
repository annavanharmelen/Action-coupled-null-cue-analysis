%% Script for doing stats on saccade and gaze bias data.
% So run those scripts first.
% by Anna, 04-07-2023
%% Saccade bias data - stats
statcfg.xax = saccade.time;
statcfg.npermutations = 10000;
statcfg.clusterStatEvalaluationAlpha = 0.05;
statcfg.nsub = size(pp2do, 2);
statcfg.statMethod = 'montecarlo';
%statcfg.statMethod = 'analytic';

% ft_size = 26;
timeframe = [951:2451]; %0 - 1500 ms post-cue
% data_cond1 = d3(:,2,timeframe);
% data_cond2 = d3(:,4,timeframe);
% data_press = (d3(:,11,timeframe) + d3(:,12,timeframe) ./ 2);
% data_notpress = (d3(:,13,timeframe) + d3(:,14,timeframe) ./ 2);

% data_old = squeeze(d3(3,:,4,timeframe));
% data_new = squeeze(d3(1,:,5,timeframe));
% null_data = zeros(size(data_new));

con_above_median = d3(:,16,timeframe);
con_below_median = d3(:,12,timeframe);
incon_above_median = d3(:,17,timeframe);
incon_below_median = d3(:,13,timeframe);

% stat1 = frevede_ftclusterstat1D(statcfg, data_cond1, null_data)
% stat2 = frevede_ftclusterstat1D(statcfg, data_cond2, null_data)
% stat_comp = frevede_ftclusterstat1D(statcfg, data_cond1, data_cond2)

% stat_press_comp = frevede_ftclusterstat1D(statcfg, data_press, data_notpress);
% stat_press = frevede_ftclusterstat1D(statcfg, data_press, null_data);
% stat_notpress = frevede_ftclusterstat1D(statcfg, data_notpress, null_data);

% stat_old = frevede_ftclusterstat1D(statcfg, data_old, null_data);
% stat_new = frevede_ftclusterstat1D(statcfg, data_new, null_data);
% stat_no_comp = frevede_ftclusterstat1D(statcfg, data_old, data_new);

% stat_con = frevede_ftclusterstat1D(statcfg, con_above_median, con_below_median);
% stat_incon = frevede_ftclusterstat1D(statcfg, incon_above_median, incon_below_median);

stat_fast_vs_slow_con_vs_incon = frevede_ftclusterstat1D(statcfg, con_below_median - con_above_median, incon_below_median - incon_above_median);
%% Saccade bias data - plot only effect
mask_xxx = double(stat_comp.mask);
mask_xxx(mask_xxx==0) = nan; % nan data that is not part of mark

figure; hold on;
ylimit = [-0.3, 0.3];
p1 = frevede_errorbarplot(saccade.time, squeeze(d3(:,2,:)), colours(2,:), 'se');
p2 = frevede_errorbarplot(saccade.time, squeeze(d3(:,4,:)), colours(1,:), 'se');
xlim(xlimtoplot);
p1.LineWidth = 1.5;
p2.LineWidth = 1.5;
legend([p1, p2], saccade.label([2,4]));
plot(xlim, [0,0], '--', 'LineWidth',2, 'Color', [0.6, 0.6, 0.6]);
plot([0,0], ylimit, '--', 'LineWidth',2, 'Color', [0.6, 0.6, 0.6]);

xlim(xlimtoplot);
sig = plot(saccade.time(timeframe), mask_xxx*-0.11, 'Color', 'k', 'LineWidth', 4); % verticaloffset for positioning of the "significance line"
ylim([-0.2, 0.2]);
ylabel('Rate effect (delta Hz)');
xlabel('Time (ms)');
% set(gcf,'position',[0,0, 1800,900])
% fontsize(ft_size*1.5,"points")
legend([p1,p2], saccade.label([4,6]));

%% Saccade bias - old vs. new
mask_comp = double(stat_no_comp.mask);
mask_comp(mask_comp==0) = nan; % nan data that is not part of mark

mask_new = double(stat_new.mask);
mask_new(mask_new==0) = nan; % nan data that is not part of mark

mask_old = double(stat_old.mask);
mask_old(mask_old==0) = nan; % nan data that is not part of mark

figure; hold on;
ylimit = [-0.3, 0.3];
p1 = frevede_errorbarplot(saccade.time, squeeze(d3(1,:,5,:)), 'k', 'both');
p2 = frevede_errorbarplot(saccade.time, squeeze(d3(3,:,4,:)), 'c', 'both');
xlim(xlimtoplot);
p1.LineWidth = 1.5;
p2.LineWidth = 1.5;
yline(0,'--', 'LineWidth',2, 'Color', [0.6, 0.6, 0.6])
xline(0,'--', 'LineWidth',2, 'Color', [0.6, 0.6, 0.6])

xlim(xlimtoplot);
sig1 = plot(saccade.time(timeframe), mask_comp*-0.20, 'Color', 'r', 'LineWidth', 4); % verticaloffset for positioning of the "significance line"
sig2 = plot(saccade.time(timeframe), mask_new*-0.22, 'Color', 'k', 'LineWidth', 4); % verticaloffset for positioning of the "significance line"
sig3 = plot(saccade.time(timeframe), mask_old*-0.24, 'Color', 'c', 'LineWidth', 4); % verticaloffset for positioning of the "significance line"
ylim([-0.3, 0.5]);
ylabel('Rate effect (delta Hz)');
xlabel('Time (ms)');
% set(gcf,'position',[0,0, 1800,900])
% fontsize(ft_size*1.5,"points")
legend([p1, p2], {'new (task)', 'old (no task)'});


%% Saccade bias - effect of pressing
mask_comp = double(stat_press_comp.mask);
mask_comp(mask_comp==0) = nan; % nan data that is not part of mark

mask_press = double(stat_press.mask);
mask_press(mask_press==0) = nan; % nan data that is not part of mark

mask_notpress = double(stat_notpress.mask);
mask_notpress(mask_notpress==0) = nan; % nan data that is not part of mark

figure; hold on;
ylimit = [-0.3, 0.3];
p1 = frevede_errorbarplot(saccade.time, squeeze((d3(:,11,:) + d3(:,12,:)) ./ 2), 'b', 'both');
p2 = frevede_errorbarplot(saccade.time, squeeze((d3(:,13,:) + d3(:,14,:)) ./ 2), 'r', 'both');
xlim(xlimtoplot);
p1.LineWidth = 1.5;
p2.LineWidth = 1.5;
plot(xlim, [0,0], '--', 'LineWidth',2, 'Color', [0.6, 0.6, 0.6]);
plot([0,0], ylimit, '--', 'LineWidth',2, 'Color', [0.6, 0.6, 0.6]);

xlim(xlimtoplot);
sig1 = plot(saccade.time(timeframe), mask_comp*-0.20, 'Color', 'k', 'LineWidth', 4); % verticaloffset for positioning of the "significance line"
sig2 = plot(saccade.time(timeframe), mask_press*-0.22, 'Color', 'b', 'LineWidth', 4); % verticaloffset for positioning of the "significance line"
sig3 = plot(saccade.time(timeframe), mask_notpress*-0.24, 'Color', 'r', 'LineWidth', 4); % verticaloffset for positioning of the "significance line"
ylim([-0.3, 0.5]);
ylabel('Rate effect (delta Hz)');
xlabel('Time (ms)');
% set(gcf,'position',[0,0, 1800,900])
% fontsize(ft_size*1.5,"points")
legend([p1, p2], {'should have pressed', 'should not have pressed'});

