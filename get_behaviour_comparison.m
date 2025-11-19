%% Script for comparing the behavioural results from the action-coupled cue to a non-action one

%% start clean
clear; clc; close all;

%% parameters
display_percentageok = 1;
plot_individuals = 0;
plot_averages = 0;

pp2do = [[1:9,11:26];[1:8, 10:16, 18:27];[1:25];[3:10,12:17,19:23,25:30]];

%% load and aggregate the data from each study
congruency_labels = {"match", "match", "match", };

for s = [1:4]
    disp([newline(), 'getting data from study ', num2str(s)]);
    p = 0;

    for pp = pp2do(s,:)
        p = p+1;
    
        % get participant data
        param = getAllSubjParam(s, pp);
    
        % load
        disp(['getting data from participant ', param.subjName]);
        behdata = readtable(param.log);
    
        % update signed error to stay within -90/+90
        behdata.signed_difference(behdata.signed_difference>90) = behdata.signed_difference(behdata.signed_difference>90)-180;
        behdata.signed_difference(behdata.signed_difference<-90) = behdata.signed_difference(behdata.signed_difference<-90)+180;
        
        % check ok trials, just based on decision time, because this one is unlimited.
        oktrials = abs(zscore(behdata.idle_reaction_time_in_ms))<=3; 
        percentageok(s,p) = mean(oktrials)*100;

        %display percentage ok
        if display_percentageok
        fprintf('%s has %.2f%% OK trials ', param.subjName, percentageok(s,p))
        fprintf('and an average score of %.2f \n', mean(behdata.performance))
        end

        % trial selections
        congruent_trials = ismember(behdata.trial_condition, {'congruent'});
        incongruent_trials = ismember(behdata.trial_condition, {'incongruent'});
        
        if s == 1 || s == 4
            neutral_trials = ismember(behdata.trial_condition, {'neutral'});
        else
            cc_trials = ismember(behdata.block_type, {'colour_probe'})&ismember(behdata.cue_form, {'colour_cue'});
            oktrials = oktrials&cc_trials;
        end
   
        % extract data of interest
        overall_dt(s,p) = mean(behdata.idle_reaction_time_in_ms(oktrials));
        overall_error(s,p) = mean(behdata.absolute_difference(oktrials));
        overall_perf(s,p) =  mean(behdata.performance(oktrials));
    
        congruency_dt(s,p,1) = mean(behdata.idle_reaction_time_in_ms(congruent_trials&oktrials));
        congruency_dt(s,p,2) = mean(behdata.idle_reaction_time_in_ms(incongruent_trials&oktrials));
    
        congruency_er(s,p,1) = mean(behdata.absolute_difference(congruent_trials&oktrials));
        congruency_er(s,p,2) = mean(behdata.absolute_difference(incongruent_trials&oktrials));
        
        if s == 1 || s == 4
            congruency_dt(s,p,3) = mean(behdata.idle_reaction_time_in_ms(neutral_trials&oktrials));
            congruency_er(s,p,3) = mean(behdata.absolute_difference(neutral_trials&oktrials));
        end

    end

    % calculate aggregates
    con_dt_effect(:,s) = congruency_dt(s,:,2) - congruency_dt(s,:,1);
    con_er_effect(:,s) = congruency_er(s,:,2) - congruency_er(s,:,1);

end


if plot_averages
%% plot pure congruecy data for all studies
figure;

% study 1
subplot(4,2,1)
hold on 
b = bar(mean(squeeze(congruency_dt(1,:,1:3))));
errorbar([1:3], mean(squeeze(congruency_dt(1,:,1:3))), std(squeeze(congruency_dt(1,:,1:3))) ./ sqrt(size(pp2do(1,:), 2)), 'LineStyle', 'none', 'LineWidth', 1.5)
xticks([1,2,3]);
xticklabels(congruency_labels);
title(['dt as function of congruency']);
% add individuals
plot([1:3], squeeze(congruency_dt(1,:,1:3)), 'Color', [0, 0, 0, 0.25]);

subplot(4,2,2)
hold on 
b = bar(mean(squeeze(congruency_er(1,:,1:3))));
errorbar([1:3], mean(squeeze(congruency_er(1,:,1:3))), std(squeeze(congruency_er(1,:,1:3))) ./ sqrt(size(pp2do(1,:), 2)), 'LineStyle', 'none', 'LineWidth', 1.5)
xticks([1,2,3]);
xticklabels(congruency_labels);
title(['er as function of congruency']);
% add individuals
plot([1:3], squeeze(congruency_er(1,:,1:3)), 'Color', [0, 0, 0, 0.25]);

% study 2
subplot(4,2,3)
hold on 
b = bar(mean(squeeze(congruency_dt(2,:,1:2))));
errorbar([1:2], mean(squeeze(congruency_dt(2,:,1:2))), std(squeeze(congruency_dt(2,:,1:2))) ./ sqrt(size(pp2do(2,:), 2)), 'LineStyle', 'none', 'LineWidth', 1.5)
xticks([1,2]);
xticklabels(congruency_labels);
% add individuals
plot([1:2], squeeze(congruency_dt(2,:,1:2)), 'Color', [0, 0, 0, 0.25]);

subplot(4,2,4)
hold on 
b = bar(mean(squeeze(congruency_er(2,:,1:2))));
errorbar([1:2], mean(squeeze(congruency_er(2,:,1:2))), std(squeeze(congruency_er(2,:,1:2))) ./ sqrt(size(pp2do(2,:), 2)), 'LineStyle', 'none', 'LineWidth', 1.5)
xticks([1,2]);
xticklabels(congruency_labels);
% add individuals
plot([1:2], squeeze(congruency_er(2,:,1:2)), 'Color', [0, 0, 0, 0.25]);

% study 3
subplot(4,2,5)
hold on 
b = bar(mean(squeeze(congruency_dt(3,:,1:2))));
errorbar([1:2], mean(squeeze(congruency_dt(3,:,1:2))), std(squeeze(congruency_dt(3,:,1:2))) ./ sqrt(size(pp2do(3,:), 2)), 'LineStyle', 'none', 'LineWidth', 1.5)
xticks([1,2]);
xticklabels(congruency_labels);
% add individuals
plot([1:2], squeeze(congruency_dt(3,:,1:2)), 'Color', [0, 0, 0, 0.25]);

subplot(4,2,6)
hold on 
b = bar(mean(squeeze(congruency_er(3,:,1:2))));
errorbar([1:2], mean(squeeze(congruency_er(3,:,1:2))), std(squeeze(congruency_er(3,:,1:2))) ./ sqrt(size(pp2do(3,:), 2)), 'LineStyle', 'none', 'LineWidth', 1.5)
xticks([1,2]);
xticklabels(congruency_labels);
% add individuals
plot([1:2], squeeze(congruency_er(3,:,1:2)), 'Color', [0, 0, 0, 0.25]);

% study 1
subplot(4,2,7)
hold on 
b = bar(mean(squeeze(congruency_dt(4,:,1:3))));
errorbar([1:3], mean(squeeze(congruency_dt(4,:,1:3))), std(squeeze(congruency_dt(4,:,1:3))) ./ sqrt(size(pp2do(1,:), 2)), 'LineStyle', 'none', 'LineWidth', 1.5)
xticks([1,2,3]);
xticklabels(congruency_labels);
title(['dt as function of congruency']);
% add individuals
plot([1:3], squeeze(congruency_dt(4,:,1:3)), 'Color', [0, 0, 0, 0.25]);

subplot(4,2,8)
hold on 
b = bar(mean(squeeze(congruency_er(4,:,1:3))));
errorbar([1:3], mean(squeeze(congruency_er(4,:,1:3))), std(squeeze(congruency_er(4,:,1:3))) ./ sqrt(size(pp2do(1,:), 2)), 'LineStyle', 'none', 'LineWidth', 1.5)
xticks([1,2,3]);
xticklabels(congruency_labels);
title(['er as function of congruency']);
% add individuals
plot([1:3], squeeze(congruency_er(4,:,1:3)), 'Color', [0, 0, 0, 0.25]);


%% plot congruency effects between studies
% bar
figure;
subplot(2,2,1)
hold on
b = bar(mean(con_dt_effect));
errorbar([1:4], mean(con_dt_effect), std(con_dt_effect) ./ sqrt(size(pp2do,2)), 'LineWidth', 1.5);
xticks([1,2,3]);
xticklabels({"Study 1", "Study 2", "Study 3", "Study 4"});
title('dt congruency effect')

subplot(2,2,2)
hold on
b = bar(mean(con_er_effect));
errorbar([1:4], mean(con_er_effect), std(con_er_effect) ./ sqrt(size(pp2do,2)), 'LineWidth', 1.5);
xticks([1,2,3]);
xticklabels({"Study 1", "Study 2", "Study 3"});
title('er congruency effect')

% scatter with individuals
subplot(2,2,3)
scatter([ones(25,1), ones(25,1)*2, ones(25,1)*3, ones(25,1)*4], con_dt_effect);
xlimit([0, 4]);
xticks([1,2,3,4]);
xticklabels({"Study 1", "Study 2", "Study 3", "Study 4"});

subplot(2,2,4)
scatter([ones(25,1), ones(25,1)*2, ones(25,1)*3, ones(25,1)*4], con_er_effect);
xlimit([0, 4])
xticks([1,2,3,4]);
xticklabels({"Study 1", "Study 2", "Study 3", "Study 4"});

%% congruency effect dt x error mega scatter
figure;
scatter(con_dt_effect, con_er_effect);
xlabel('congruency effect on DECISION TIME');
ylabel('congruency effect on ERROR');
legend({"Study 1", "Study 2", "Study 3", "Study 4"});

%% congruency effect saccade x behaviour mega scatter
% only works if you run GA_saccadeBias_comparison before this
figure;
subplot(1,2,1)
hold on
scatter(squeeze(avg_saccade_effect(1,:,1)), con_dt_effect(:,1), marker="diamond", markerFaceColor=get_colour("pink",""), MarkerEdgeColor='none');
scatter(squeeze(avg_saccade_effect(1,:,2)), con_dt_effect(:,2), marker="diamond", markerFaceColor=[0.5, 0.5, 0.5], MarkerEdgeColor='none');
scatter(squeeze(avg_saccade_effect(1,:,3)), con_dt_effect(:,3), marker="o", markerFaceColor=[0.5, 0.5, 0.5], MarkerEdgeColor='none');
scatter(squeeze(avg_saccade_effect(1,:,4)), con_dt_effect(:,4), marker="o", markerFaceColor=get_colour("pink",""), MarkerEdgeColor='none');
ylabel('congruency effect of DECISION TIME');
xlabel('congruency effect of SACCADE BIAS');
title('DT effect x saccade bias');
b_dt = reshape(avg_saccade_effect(1,:,:), [100,1])\reshape(con_dt_effect, [100,1]);
y_dt = b_dt*reshape(avg_saccade_effect(1,:,:), [100,1]);
plot(reshape(avg_saccade_effect(1,:,:), [100,1]),y_dt, LineWidth=1.0, Color=[0,0,0]);
% h = lsline();
% set(h(4),'color',[17,113,190]/255);
% set(h(3),'color',[17,113,190]/255);
% set(h(2),'color',[221,84,0]/255);
% set(h(1),'color',[237,177,32]/255);
legend({"Study 1", "Study 2", "Study 3", "Study 4"});

subplot(1,2,2)
hold on
scatter(squeeze(avg_saccade_effect(1,:,1)), con_er_effect(:,1), marker="diamond", markerFaceColor=get_colour("pink",""), MarkerEdgeColor='none');
scatter(squeeze(avg_saccade_effect(1,:,2)), con_er_effect(:,2), marker="diamond", markerFaceColor=[0.5, 0.5, 0.5], MarkerEdgeColor='none');
scatter(squeeze(avg_saccade_effect(1,:,3)), con_er_effect(:,3), marker="o", markerFaceColor=[0.5, 0.5, 0.5], MarkerEdgeColor='none');
scatter(squeeze(avg_saccade_effect(1,:,4)), con_er_effect(:,4), marker="o", markerFaceColor=get_colour("pink",""), MarkerEdgeColor='none');
ylabel('congruency effect of ERROR');
xlabel('congruency effect of SACCADE BIAS');
title('Error effect x saccade bias');
b_er = reshape(avg_saccade_effect(1,:,:), [100,1])\reshape(con_er_effect, [100,1]);
y_er = b_er*reshape(avg_saccade_effect(1,:,:), [100,1]);
plot(reshape(avg_saccade_effect(1,:,:), [100,1]),y_er, LineWidth=1.0, Color=[0,0,0]);
% h = lsline();
% set(h(3),'color',[17,113,190]/255);
% set(h(2),'color',[221,84,0]/255);
% set(h(1),'color',[237,177,32]/255);
legend({"Study 1", "Study 2", "Study 3", "Study 4"});

[r,pval] = corr(con_dt_effect(:,1), avg_saccade_effect(1,:,1)');
if pval <= 0.05
    disp('study 1 has corr for dt');
    disp(['r= ', num2str(r), 'p= ', num2str(pval)]);
end
[r,pval] = corr(con_dt_effect(:,2), avg_saccade_effect(1,:,2)');
if pval <= 0.05
    disp('study 2 has corr for dt');
    disp(['r= ', num2str(r), 'p= ', num2str(pval)]);
end
[r,pval] = corr(con_dt_effect(:,3), avg_saccade_effect(1,:,3)');
if pval <= 0.05
    disp('study 3 has corr for dt');
    disp(['r= ', num2str(r), 'p= ', num2str(pval)]);
end
[r,pval] = corr(con_dt_effect(:,4), avg_saccade_effect(1,:,4)');
if pval <= 0.05
    disp('study 4 has corr for dt');
    disp(['r= ', num2str(r), 'p= ', num2str(pval)]);
end

[r,pval] = corr(con_er_effect(:,1), avg_saccade_effect(1,:,1)');
if pval <= 0.05
    disp('study 1 has corr for er');
    disp(['r= ', num2str(r), 'p= ', num2str(pval)]);
end
[r,pval] = corr(con_er_effect(:,2), avg_saccade_effect(1,:,2)');
if pval <= 0.05
    disp('study 2 has corr for er');
    disp(['r= ', num2str(r), 'p= ', num2str(pval)]);
end
[r,pval] = corr(con_er_effect(:,3), avg_saccade_effect(1,:,3)');
if pval <= 0.05
    disp('study 3 has corr for er');
    disp(['r= ', num2str(r), 'p= ', num2str(pval)]);
end
[r,pval] = corr(con_er_effect(:,4), avg_saccade_effect(1,:,4)');
if pval <= 0.05
    disp('study 4 has corr for er');
    disp(['r= ', num2str(r), 'p= ', num2str(pval)]);
end

%% congruency effect saccade x behaviour less-mega scatter
% only works if you run GA_saccadeBias_comparison before this
figure;
subplot(1,2,1)
hold on
scatter(reshape(con_dt_effect(:,[1,4]), [50,1]), reshape(squeeze(avg_saccade_effect(1,:,[1,4])), [50,1]));
scatter(reshape(con_dt_effect(:,2:3), [50,1]), reshape(squeeze(avg_saccade_effect(1,:,2:3)), [50,1]));
xlabel('congruency effect of DECISION TIME');
ylabel('congruency effect of SACCADE BIAS');
title('DT effect x saccade bias');
h = lsline();
set(h(1),'color',[17,113,190]/255);
set(h(2),'color',[221,84,0]/255);
legend({"Study 1 + 4", "Study 2 + 3"});

subplot(1,2,2)
hold on
scatter(reshape(con_er_effect(:,[1,4]), [50,1]), reshape(squeeze(avg_saccade_effect(1,:,[1,4])), [50,1]));
scatter(reshape(con_er_effect(:,2:3), [50,1]), reshape(squeeze(avg_saccade_effect(1,:,2:3)), [50,1]));
xlabel('congruency effect of DECISION TIME');
ylabel('congruency effect of SACCADE BIAS');
title('ER effect x saccade bias');
h = lsline();
set(h(2),'color',[17,113,190]/255);
set(h(1),'color',[221,84,0]/255);
legend({"Study 1 + 4", "Study 2 + 3"});

[r,pval] = corr(con_dt_effect(:,1), squeeze(avg_saccade_effect(1,:,1))');
if pval <= 0.05
    disp('study 1 has corr for dt');
    disp(['r= ', num2str(r), 'p= ', num2str(pval)]);
end
[r,pval] = corr(reshape(con_dt_effect(:,2:3), [50,1]), reshape(squeeze(avg_saccade_effect(1,:,2:3)), [50,1]));
if pval <= 0.05
    disp('study 2/3 has corr for dt');
    disp(['r= ', num2str(r), 'p= ', num2str(pval)]);
end

[r,pval] = corr(con_er_effect(:,1), squeeze(avg_saccade_effect(1,:,1))');
if pval <= 0.05
    disp('study 1 has corr for er');
    disp(['r= ', num2str(r), 'p= ', num2str(pval)]);
end
[r,pval] = corr(reshape(con_dt_effect(:,2:3), [50,1]), reshape(squeeze(avg_saccade_effect(1,:,2:3)), [50,1]));
if pval <= 0.05
    disp('study 2/3 has corr for er');
    disp(['r= ', num2str(r), 'p= ', num2str(pval)]);
end

%% main figure for poster?
dt_ylim = [200 1400];
dt_yticks = [200, 600, 1000, 1400];
er_ylim = [5 35];
er_yticks = [5, 15, 25, 35];
xlimit = [0.35 2.65];
bar_width = 0.64;


figure;
tL = subplot(2,3,1);
hold on 
b = bar(mean(squeeze(mean(congruency_dt([1,4],:,1:2), 2))), 'FaceColor', 'flat', 'LineStyle', 'none', 'BarWidth', bar_width);
b.CData(1,:) = get_colour("blue", "");
b.CData(2,:) = get_colour("red", "");
errorbar([1], mean(squeeze(mean(congruency_dt([1,4],:,1), 2))), std(congruency_dt([1,4],:,1),[],"all") ./ sqrt(size(pp2do, 2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("blue", "dark"))
errorbar([2], mean(squeeze(mean(congruency_dt([1,4],:,2), 2))), std(congruency_dt([1,4],:,2),[],"all") ./ sqrt(size(pp2do, 2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("red", "dark"))
xticks([1,2]);
xticklabels(congruency_labels(1:2));
xlim(xlimit);
ylim(dt_ylim);
yticks(dt_yticks);
ylabel('Decision time (ms)');
% add individuals
plot([1:2], reshape(squeeze(congruency_dt([1,4],:,1:2)), [50 2]), 'Color', [0, 0, 0, 0.25], 'LineWidth', 0.75);

tM = subplot(2,3,2);
hold on 
b = bar(mean(squeeze(mean(congruency_dt(2:3,:,1:2), 2))), 'FaceColor', 'flat', 'LineStyle', 'none', 'BarWidth', bar_width);
b.CData(1,:) = get_colour("blue", "");
b.CData(2,:) = get_colour("red", "");
errorbar([1], mean(squeeze(mean(congruency_dt(2:3,:,1), 2))), std(congruency_dt(2:3,:,1),[],"all") ./ sqrt(size(pp2do, 2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("blue", "dark"))
errorbar([2], mean(squeeze(mean(congruency_dt(2:3,:,2), 2))), std(congruency_dt(2:3,:,2),[],"all") ./ sqrt(size(pp2do, 2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("red", "dark"))
xticks([1,2]);
xticklabels(congruency_labels(1:2));
xlim(xlimit);
ylim(dt_ylim);
yticks(dt_yticks);
% add individuals
plot([1:2], reshape(squeeze(congruency_dt(2:3,:,1:2)), [50 2]), 'Color', [0, 0, 0, 0.25], 'LineWidth', 0.75);


tR = subplot(2,3,3);
hold on
b = bar([mean(con_dt_effect(:,[1,4]), "all"), mean(con_dt_effect(:,2:3), "all")], 'FaceColor', 'flat', 'LineStyle', 'none', 'BarWidth', bar_width);
b.CData(1,:) = get_colour("pink", "");
b.CData(2,:) = [0.5, 0.5, 0.5];
errorbar([1], mean(con_dt_effect(:,[1,4]), "all"), std(con_dt_effect(:,[1,4]), [], "all") ./ sqrt(size(pp2do,2)), 'LineWidth', 1.5, 'Color', get_colour("pink", "dark"));
errorbar([2], mean(con_dt_effect(:,2:3), "all"), std(con_dt_effect(:,2:3), [], "all") ./ sqrt(size(pp2do,2)), 'LineWidth', 1.5, 'Color', [0,0,0]);
xticks([1,2]);
xticklabels({"action", "no action"});
xlim(xlimit);
ylim([0, 200]);

bL = subplot(2,3,4);
hold on 
b = bar(mean(squeeze(mean(congruency_er([1,4],:,1:2), 2))), 'FaceColor', 'flat', 'LineStyle', 'none', 'BarWidth', bar_width);
b.CData(1,:) = get_colour("blue", "");
b.CData(2,:) = get_colour("red", "");
errorbar([1], mean(squeeze(mean(congruency_er([1,4],:,1), 2))), std(congruency_er([1,4],:,1),[],"all") ./ sqrt(size(pp2do, 2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("blue", "dark"))
errorbar([2], mean(squeeze(mean(congruency_er([1,4],:,2), 2))), std(congruency_er([1,4],:,2),[],"all") ./ sqrt(size(pp2do, 2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("red", "dark"))
xticks([1,2]);
xticklabels(congruency_labels(1:2));
xlim(xlimit);
ylim(er_ylim);
yticks(er_yticks);
ylabel('Reproduction error (°)');
% add individuals
plot([1:2], reshape(squeeze(congruency_er([1,4],:,1:2)), [50 2]), 'Color', [0, 0, 0, 0.25], 'LineWidth', 0.75);


bM = subplot(2,3,5);
hold on 
b = bar(mean(squeeze(mean(congruency_er(2:3,:,1:2), 2))), 'FaceColor', 'flat', 'LineStyle', 'none', 'BarWidth', bar_width);
b.CData(1,:) = get_colour("blue", "");
b.CData(2,:) = get_colour("red", "");
errorbar([1], mean(squeeze(mean(congruency_er(2:3,:,1), 2))), std(congruency_er(2:3,:,1),[],"all") ./ sqrt(size(pp2do, 2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("blue", "dark"))
errorbar([2], mean(squeeze(mean(congruency_er(2:3,:,2), 2))), std(congruency_er(2:3,:,2),[],"all") ./ sqrt(size(pp2do, 2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("red", "dark"))
xticks([1,2]);
xticklabels(congruency_labels(1:2));
xlim(xlimit);
ylim(er_ylim);
yticks(er_yticks);
% add individuals
plot([1:2], reshape(squeeze(congruency_er(2:3,:,1:2)), [50 2]), 'Color', [0, 0, 0, 0.25], 'LineWidth', 0.75);


bR = subplot(2,3,6);
hold on
b = bar([mean(con_er_effect(:,[1,4]), "all"), mean(con_er_effect(:,2:3), "all")], 'FaceColor', 'flat', 'LineStyle', 'none', 'BarWidth', bar_width);
b.CData(1,:) = get_colour("pink", "");
b.CData(2,:) = [0.5, 0.5, 0.5];
errorbar([1], mean(con_er_effect(:,[1,4]), "all"), std(con_er_effect(:,[1,4]), [], "all") ./ sqrt(size(pp2do,2)), 'LineWidth', 1.5, 'Color', get_colour("pink", "dark"));
errorbar([2], mean(con_er_effect(:,2:3), "all"), std(con_er_effect(:,2:3), [], "all") ./ sqrt(size(pp2do,2)), 'LineWidth', 1.5, 'Color', [0,0,0]);
xticks([1,2]);
xticklabels({"action", "no action"});
xlim(xlimit);
ylim([0, 6]);
% xtickangle(30);

% general
set(gcf(), 'Position', [800 300 1200 800]);

axes = {tL, tM, tR, bL, bM, bR};
for i = 1:size(axes,2)
    set(axes{i}, 'Box', 'on');
    set(axes{i}, 'FontSize', [24.8]);
    set(axes{i}, 'FontName', 'Aptos');
    set(axes{i}.XAxis, 'FontSize', 18);
    set(axes{i}, 'LineWidth', 1.33);
end

bL.YLabel.Position = [-0.4435   20.0000   -1.0000];

% print("C:\Users\annav\Documents\Surfdrive\Conferences\ICON 2025\Figures\behaviour_effect_of_action", "-dsvg")
% print("C:\Users\annav\Documents\Surfdrive\Conferences\ICON 2025\Figures\behaviour_effect_of_action", "-dpng")

% stats stats stats
p_values = zeros(2,3);
[h,p_values(1,1),ci,stats] = ttest(congruency_dt(1,:,1), congruency_dt(1,:,2));
[h,p_values(1,2),ci,stats] = ttest(reshape(congruency_dt(2:3,:,1), [50,1]), reshape(congruency_dt(2:3,:,2), [50,1]));
[h,p_values(1,3),ci,stats] = ttest2(con_dt_effect(:,1), reshape(con_dt_effect(:,2:3), [50,1]));
[h,p_values(2,1),ci,stats] = ttest(congruency_er(1,:,1), congruency_er(1,:,2));
[h,p_values(2,2),ci,stats] = ttest(reshape(congruency_er(2:3,:,1), [50,1]), reshape(congruency_er(2:3,:,2), [50,1]));
[h,p_values(2,3),ci,stats] = ttest2(con_er_effect(:,1), reshape(con_er_effect(:,2:3), [50,1]));

end
