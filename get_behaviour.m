clear all
close all
clc

%% set parameters and loops
see_performance = 1;
display_percentageok = 1;
plot_individuals = 0;
plot_averages = 1;

pp2do = setdiff(1:30,[1,2,11,18,24]); 
p = 0;

[bar_size, colours, dark_colours, labels, subplot_size, percentageok] = setBehaviourParam(pp2do);

for pp = pp2do
    p = p+1;
    ppnum(p) = pp;
    figure_nr = 1;

    param = getSubjParam(pp);
    disp(['getting data from ', param.subjName]);
    
    %% load actual behavioural data
    behdata = readtable(param.log);
    
    %% update signed error to stay within -90/+90
    behdata.signed_difference(behdata.signed_difference>90) = behdata.signed_difference(behdata.signed_difference>90)-180;
    behdata.signed_difference(behdata.signed_difference<-90) = behdata.signed_difference(behdata.signed_difference<-90)+180;
    
    %% check ok trials, just based on decision time, because this one is unlimited.
    oktrials = abs(zscore(behdata.idle_reaction_time_in_ms))<=3; 
    % oktrials = abs(zscore(behdata.absolute_difference))<=3; 
    percentageok(p) = mean(oktrials)*100;

    %% display percentage OK
    if display_percentageok
        fprintf('%s has %.2f%% OK trials ', param.subjName, percentageok(p))
        fprintf('and an average score of %.2f \n\n', mean(behdata.performance))
    end

    %% basic data checks, each pp in own subplot
    if plot_individuals
        figure(figure_nr);
        figure_nr = figure_nr+1;
        subplot(subplot_size,subplot_size,p);
        histogram(behdata.response_time_in_ms(oktrials),50);
        title(['response time - pp ', num2str(pp2do(p))]);
        xlim([0 2000]);

        figure(figure_nr);
        figure_nr = figure_nr+1;
        subplot(subplot_size,subplot_size,p);
        histogram(behdata.idle_reaction_time_in_ms(oktrials),50);
        title(['decision time - pp ', num2str(pp2do(p))]);  

        figure(figure_nr);
        figure_nr = figure_nr+1;
        subplot(subplot_size,subplot_size,p);
        histogram(behdata.signed_difference(oktrials),50);
        title(['signed error - pp ', num2str(pp2do(p))]);
        xlim([-100 100]);

        figure(figure_nr);
        figure_nr = figure_nr+1;
        subplot(subplot_size,subplot_size,p);
        histogram(behdata.absolute_difference(oktrials),50);
        title(['error - pp ', num2str(pp2do(p))]);
        xlim([0 100]);
        
    end
    
    %% trial selections
    colour_one_trials = behdata.capture_colour_id == 1;
    colour_two_trials = behdata.capture_colour_id == 2;
    colour_three_trials = behdata.capture_colour_id == 3;
    colour_four_trials = behdata.capture_colour_id == 4;

    respond_trials = ismember(behdata.block_type, {'respond 3'});
    not_respond_trials = ismember(behdata.block_type, {'respond not 3'});
    
    pressed_trials = ismember(behdata.cue_hit, {'True'}) | ismember(behdata.cue_false_alarm, {'True'});
    unpressed_trials = ones(size(pressed_trials)) - pressed_trials;

    hit_trials = ismember(behdata.cue_hit, {'True'});
    false_alarm_trials = ismember(behdata.cue_false_alarm, {'True'});
    miss_trials = (respond_trials & behdata.capture_colour_id == 3 & unpressed_trials) | (not_respond_trials & behdata.capture_colour_id ~= 3 & unpressed_trials);
    correct_rejection_trials = (respond_trials & behdata.capture_colour_id ~= 3 & unpressed_trials) | (not_respond_trials & behdata.capture_colour_id == 3 & unpressed_trials);
    
    congruent_trials = ismember(behdata.trial_condition, {'congruent'});
    incongruent_trials = ismember(behdata.trial_condition, {'incongruent'});
    neutral_trials = ismember(behdata.trial_condition, {'neutral'});

    left_target_trials = ismember(behdata.target_bar, {'left'});
    right_target_trials = ismember(behdata.target_bar, {'right'});
    
    %% save hit/miss/false_alarm/correct_rejection
    % hits, misses, false_alarms, correct_rejections
    scores(p,1) = sum(behdata.cue_hit == "True");
    scores(p,2) = sum(miss_trials);
    scores(p,3) = sum(behdata.cue_false_alarm == "True");
    scores(p,4) = sum(correct_rejection_trials);
    scores(p,5) = sum(scores(p,1:4));


    %% mixture models of target error
    non_target_orientations = behdata.left_orientation;
    non_target_orientations(behdata.target_bar == "left") = behdata.right_orientation(behdata.target_bar == "left");
    [con_B(p,:), ~, W] = mixtureFit(deg2rad(behdata.report_orientation(congruent_trials)), deg2rad(behdata.target_orientation(congruent_trials)), deg2rad(non_target_orientations(congruent_trials)));
    [incon_B(p,:), ~, ~] = mixtureFit(deg2rad(behdata.report_orientation(incongruent_trials)), deg2rad(behdata.target_orientation(incongruent_trials)), deg2rad(non_target_orientations(incongruent_trials)));
    [neu_B(p,:), ~, ~] = mixtureFit(deg2rad(behdata.report_orientation(neutral_trials)), deg2rad(behdata.target_orientation(neutral_trials)), deg2rad(non_target_orientations(neutral_trials)));


    %% extract data of interest
    overall_dt(p,1) = mean(behdata.idle_reaction_time_in_ms(oktrials));
    overall_error(p,1) = mean(behdata.absolute_difference(oktrials));
    overall_perf(p,1) =  mean(behdata.performance(oktrials));
    
    congruency_labels = {"match", "match", "match"};

    congruency_dt(p,1) = mean(behdata.idle_reaction_time_in_ms(congruent_trials&oktrials));
    congruency_dt(p,2) = mean(behdata.idle_reaction_time_in_ms(neutral_trials&oktrials));
    congruency_dt(p,3) = mean(behdata.idle_reaction_time_in_ms(incongruent_trials&oktrials));

    congruency_er(p,1) = mean(behdata.absolute_difference(congruent_trials&oktrials));
    congruency_er(p,2) = mean(behdata.absolute_difference(neutral_trials&oktrials));
    congruency_er(p,3) = mean(behdata.absolute_difference(incongruent_trials&oktrials));

    pressed_labels = {"pressed", "unpressed"};
      
    pressed_dt(p,1) = mean(behdata.idle_reaction_time_in_ms(pressed_trials&oktrials));
    pressed_dt(p,2) = mean(behdata.idle_reaction_time_in_ms(unpressed_trials&oktrials));
   
    pressed_congruency_dt(p,1) = mean(behdata.idle_reaction_time_in_ms(pressed_trials&oktrials&congruent_trials));
    pressed_congruency_dt(p,2) = mean(behdata.idle_reaction_time_in_ms(pressed_trials&oktrials&incongruent_trials));
    pressed_congruency_dt(p,3) = mean(behdata.idle_reaction_time_in_ms(unpressed_trials&oktrials&congruent_trials));
    pressed_congruency_dt(p,4) = mean(behdata.idle_reaction_time_in_ms(unpressed_trials&oktrials&incongruent_trials));
    
    pressed_er(p,1) = mean(behdata.absolute_difference(pressed_trials&oktrials));
    pressed_er(p,2) = mean(behdata.absolute_difference(unpressed_trials&oktrials));
   
    pressed_congruency_er(p,1) = mean(behdata.absolute_difference(pressed_trials&oktrials&congruent_trials));
    pressed_congruency_er(p,2) = mean(behdata.absolute_difference(pressed_trials&oktrials&incongruent_trials));
    pressed_congruency_er(p,3) = mean(behdata.absolute_difference(unpressed_trials&oktrials&congruent_trials));
    pressed_congruency_er(p,4) = mean(behdata.absolute_difference(unpressed_trials&oktrials&incongruent_trials));
    
    cue_response_labels = {"hit", "miss", "false alarm", "correct rejection"};

    cue_response_dt(p,1) = mean(behdata.idle_reaction_time_in_ms(hit_trials&oktrials), "omitnan");
    cue_response_dt(p,2) = mean(behdata.idle_reaction_time_in_ms(miss_trials&oktrials), "omitnan");
    cue_response_dt(p,3) = mean(behdata.idle_reaction_time_in_ms(false_alarm_trials&oktrials), "omitnan");
    cue_response_dt(p,4) = mean(behdata.idle_reaction_time_in_ms(correct_rejection_trials&oktrials), "omitnan");
    
    cue_response_er(p,1) = mean(behdata.absolute_difference(hit_trials&oktrials), "omitnan");
    cue_response_er(p,2) = mean(behdata.absolute_difference(miss_trials&oktrials), "omitnan");
    cue_response_er(p,3) = mean(behdata.absolute_difference(false_alarm_trials&oktrials), "omitnan");
    cue_response_er(p,4) = mean(behdata.absolute_difference(correct_rejection_trials&oktrials), "omitnan");
    
    block_labels = {"respond_3", "respond_not_3"};

    block_dt(p,1) = mean(behdata.idle_reaction_time_in_ms(respond_trials&oktrials));
    block_dt(p,2) = mean(behdata.idle_reaction_time_in_ms(not_respond_trials&oktrials));
    
    block_congruency_dt(p,1) = mean(behdata.idle_reaction_time_in_ms(respond_trials&oktrials&congruent_trials));
    block_congruency_dt(p,2) = mean(behdata.idle_reaction_time_in_ms(respond_trials&oktrials&incongruent_trials));
    block_congruency_dt(p,3) = mean(behdata.idle_reaction_time_in_ms(not_respond_trials&oktrials&congruent_trials));
    block_congruency_dt(p,4) = mean(behdata.idle_reaction_time_in_ms(not_respond_trials&oktrials&incongruent_trials));

    block_er(p,1) = mean(behdata.absolute_difference(respond_trials&oktrials));
    block_er(p,2) = mean(behdata.absolute_difference(not_respond_trials&oktrials));
    
    block_congruency_er(p,1) = mean(behdata.absolute_difference(respond_trials&oktrials&congruent_trials));
    block_congruency_er(p,2) = mean(behdata.absolute_difference(respond_trials&oktrials&incongruent_trials));
    block_congruency_er(p,3) = mean(behdata.absolute_difference(not_respond_trials&oktrials&congruent_trials));
    block_congruency_er(p,4) = mean(behdata.absolute_difference(not_respond_trials&oktrials&incongruent_trials));

    %% calculate aggregates of interest
   pressed_dt_effect(p,1) = pressed_congruency_dt(p,2) - pressed_congruency_dt(p,1);
   pressed_dt_effect(p,2) = pressed_congruency_dt(p,4) - pressed_congruency_dt(p,3);
    
   pressed_er_effect(p,1) = pressed_congruency_er(p,2) - pressed_congruency_er(p,1);
   pressed_er_effect(p,2) = pressed_congruency_er(p,4) - pressed_congruency_er(p,3);

   block_dt_effect(p,1) = block_congruency_dt(p,2) - block_congruency_dt(p,1);
   block_dt_effect(p,2) = block_congruency_dt(p,4) - block_congruency_dt(p,3);

   block_er_effect(p,1) = block_congruency_er(p,2) - block_congruency_er(p,1);
   block_er_effect(p,2) = block_congruency_er(p,4) - block_congruency_er(p,3);

    %% plot individuals
    if plot_individuals
        figure(figure_nr);
        figure_nr = figure_nr+1;
        subplot(5,5,p);
        title(p);
        b = bar(congruency_dt(p,:), 'FaceColor', 'flat', 'LineStyle', 'none');
        b.CData(1,:) = get_colour("blue", "");
        b.CData(2,:) = get_colour("green", "");
        b.CData(3,:) = get_colour("red", "");
    end
end

if plot_averages
%% plot mixture moduling results
    figure;
    tL = subplot(2,2,1);
    hold on
    b = bar([mean(con_B(:,1)), mean(neu_B(:,1)), mean(incon_B(:,1))], 'faceColor', 'flat', 'EdgeColor','none');
    b.CData(1,:) = get_colour("blue", "");
    b.CData(2,:) = get_colour("green", "");
    b.CData(3,:) = get_colour("red", "");
    errorbar(1, mean(con_B(:,1)), std(con_B(:,1)) ./ sqrt(size(pp2do,2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("blue", "dark"));
    errorbar(2, mean(neu_B(:,1)), std(neu_B(:,1)) ./ sqrt(size(pp2do,2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("green", "dark"));
    errorbar(3, mean(incon_B(:,1)), std(incon_B(:,1)) ./ sqrt(size(pp2do,2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("red", "dark"));
    % add individuals
    plot([1,2,3], [con_B(:,1), neu_B(:,1), incon_B(:,1)], 'Color', [0, 0, 0, 0.25], 'LineWidth', 0.75);
    xlim([0.35 3.65]);
    xticks([1,2,3]);
    xticklabels(congruency_labels);
    ylabel('Precision (K)');

    tR = subplot(2,2,2);
    hold on
    b = bar([mean(con_B(:,2)), mean(neu_B(:,2)), mean(incon_B(:,2))], 'faceColor', 'flat', 'EdgeColor','none');
    b.CData(1,:) = get_colour("blue", "");
    b.CData(2,:) = get_colour("green", "");
    b.CData(3,:) = get_colour("red", "");
    errorbar(1, mean(con_B(:,2)), std(con_B(:,2)) ./ sqrt(size(pp2do,2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("blue", "dark"));
    errorbar(2, mean(neu_B(:,2)), std(neu_B(:,2)) ./ sqrt(size(pp2do,2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("green", "dark"));
    errorbar(3, mean(incon_B(:,2)), std(incon_B(:,2)) ./ sqrt(size(pp2do,2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("red", "dark"));
    % add individuals
    plot([1,2,3], [con_B(:,2), neu_B(:,2), incon_B(:,2)], 'Color', [0, 0, 0, 0.25], 'LineWidth', 0.75);
    ylabel('P(target response)');
    ylim([0.5 1]);
    yticks([0.5, 0.75, 1]);
    xlim([0.35 3.65]);
    xticks([1,2,3]);
    xticklabels(congruency_labels);

    bL = subplot(2,2,3);
    hold on
    b = bar([mean(con_B(:,3)), mean(neu_B(:,3)), mean(incon_B(:,3))], 'faceColor', 'flat', 'EdgeColor','none');
    b.CData(1,:) = get_colour("blue", "");
    b.CData(2,:) = get_colour("green", "");
    b.CData(3,:) = get_colour("red", "");
    errorbar(1, mean(con_B(:,3)), std(con_B(:,3)) ./ sqrt(size(pp2do,2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("blue", "dark"));
    errorbar(2, mean(neu_B(:,3)), std(neu_B(:,3)) ./ sqrt(size(pp2do,2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("green", "dark"));
    errorbar(3, mean(incon_B(:,3)), std(incon_B(:,3)) ./ sqrt(size(pp2do,2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("red", "dark"));
    % add individuals
    plot([1,2,3], [con_B(:,3), neu_B(:,3), incon_B(:,3)], 'Color', [0, 0, 0, 0.25], 'LineWidth', 0.75);
    ylabel('P(non-target response)');
    xlim([0.35 3.65]);
    xticks([1,2,3]);
    xticklabels(congruency_labels);

    bR = subplot(2,2,4);
    hold on
    b = bar([mean(con_B(:,4)), mean(neu_B(:,4)), mean(incon_B(:,4))], 'faceColor', 'flat', 'EdgeColor','none');
    b.CData(1,:) = get_colour("blue", "");
    b.CData(2,:) = get_colour("green", "");
    b.CData(3,:) = get_colour("red", "");
    errorbar(1, mean(con_B(:,4)), std(con_B(:,4)) ./ sqrt(size(pp2do,2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("blue", "dark"));
    errorbar(2, mean(neu_B(:,4)), std(neu_B(:,4)) ./ sqrt(size(pp2do,2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("green", "dark"));
    errorbar(3, mean(incon_B(:,4)), std(incon_B(:,4)) ./ sqrt(size(pp2do,2)), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("red", "dark"));
    % add individuals
    plot([1,2,3], [con_B(:,4), neu_B(:,4), incon_B(:,4)], 'Color', [0, 0, 0, 0.25], 'LineWidth', 0.75);
    ylabel('P(uniform response)');
    ylim([0 0.4]);
    yticks([0, 0.2, 0.4]);
    xlim([0.35 3.65]);
    xticks([1,2,3]);
    xticklabels(congruency_labels);

    % general
    set(gcf(), 'Position', [800 300 880 800]);

    axes = {tL, tR, bL, bR};
    for i = 1:size(axes,2)
        set(axes{i}, 'Box', 'on');
        set(axes{i}, 'FontSize', [24.8]);
        set(axes{i}, 'FontName', 'Aptos');
        set(axes{i}.XAxis, 'FontSize', 18);
        set(axes{i}, 'LineWidth', 1.33);
    end

    tL.YLabel.Position = [-0.3717   10.3214   -1.0000];
    bR.YLabel.Position = [-0.6498    0.2000   -1.0000];

    tL.Position = [0.1469    0.5838    0.3177    0.3412];
    bL.Position = [0.1469    0.1100    0.3177    0.3412];
    tR.Position = [0.6530    0.5838    0.3177    0.3412];
    bR.Position = [0.6530    0.1100    0.3177    0.3412];

%% stats over previous figure (mixture modelling)
    % some stats
    n = size(pp2do,2);
    K_tbl = table([1:n,1:n,1:n]', [con_B(:,1);neu_B(:,1);incon_B(:,1)],[ones(n,1);ones(n,1)*2;ones(n,1)*3], 'VariableNames', {'participant_id', 'precision', 'condition'});
    pT_tbl = table([1:n,1:n,1:n]', [con_B(:,2);neu_B(:,2);incon_B(:,2)],[ones(n,1);ones(n,1)*2;ones(n,1)*3], 'VariableNames', {'participant_id', 'pT', 'condition'});
    pNT_tbl = table([1:n,1:n,1:n]', [con_B(:,3);neu_B(:,3);incon_B(:,3)],[ones(n,1);ones(n,1)*2;ones(n,1)*3], 'VariableNames', {'participant_id', 'pNT', 'condition'});
    pU_tbl = table([1:n,1:n,1:n]', [con_B(:,4);neu_B(:,4);incon_B(:,4)],[ones(n,1);ones(n,1)*2;ones(n,1)*3], 'VariableNames', {'participant_id', 'pU', 'condition'});

    K_lme = fitlme(K_tbl,'precision ~ condition + (1|participant_id)')
    anova(K_lme)
    pT_lme = fitlme(pT_tbl,'pT ~ condition + (1|participant_id)')
    anova(pT_lme)
    pNT_lme = fitlme(pNT_tbl,'pNT ~ condition + (1|participant_id)')
    anova(pNT_lme)
    pU_lme = fitlme(pU_tbl,'pU ~ condition + (1|participant_id)')
    anova(pU_lme)

    %post-hoc pairwise comparison
    mm_p_values = zeros(3,4)
    [h,mm_p_values(1,1),ci,stats] = ttest(con_B(:,1), neu_B(:,1));
    [h,mm_p_values(2,1),ci,stats] = ttest(neu_B(:,1), incon_B(:,1));
    [h,mm_p_values(3,1),ci,stats] = ttest(con_B(:,1), incon_B(:,1));

    [h,mm_p_values(1,2),ci,stats] = ttest(con_B(:,2), neu_B(:,2));
    [h,mm_p_values(2,2),ci,stats] = ttest(neu_B(:,2), incon_B(:,2));
    [h,mm_p_values(3,2),ci,stats] = ttest(con_B(:,2), incon_B(:,2));

    [h,mm_p_values(1,3),ci,stats] = ttest(con_B(:,3), neu_B(:,3));
    [h,mm_p_values(2,3),ci,stats] = ttest(neu_B(:,3), incon_B(:,3));
    [h,mm_p_values(3,3),ci,stats] = ttest(con_B(:,3), incon_B(:,3));

    [h,mm_p_values(1,4),ci,stats] = ttest(con_B(:,4), neu_B(:,4));
    [h,mm_p_values(2,4),ci,stats] = ttest(neu_B(:,4), incon_B(:,4));
    [h,mm_p_values(3,4),ci,stats] = ttest(con_B(:,4), incon_B(:,4));

%% all pp plot
    figure; 
    subplot(4,1,1);
    bar(ppnum, overall_dt(:,1));
    title('overall decision time');
    ylim([0 900]);
    xlabel('pp #');

    subplot(4,1,2);
    bar(ppnum, overall_error(:,1));
    title('overall error');
    xlabel('pp #');

    subplot(4,1,3);
    bar(ppnum, overall_perf(:,1));
    title('overall performance');
    xlabel('pp #');
    
    subplot(4,1,4);
    bar(ppnum, percentageok);
    title('percentage ok trials');
    ylim([90 100]);
    xlabel('pp #');

    %% does it work at all?
    %congruent vs. incongruent cue
    
    % decision time
    figure;
    hold on 
    b = bar(mean(congruency_dt), 'FaceColor', colours(3,:), 'LineStyle', 'none');
    errorbar([1:3], mean(congruency_dt), std(congruency_dt) ./ sqrt(p), 'LineStyle', 'none', 'Color', dark_colours(3,:), 'LineWidth', 1.5)
    xticks([1,2,3]);
    xticklabels(congruency_labels);
    title(['dt as function of congruency']);
    % add individuals
    plot([1:3], congruency_dt, 'Color', [0, 0, 0, 0.25]);
    
    % error
    figure;
    hold on 
    b = bar(mean(congruency_er), 'FaceColor', colours(3,:), 'LineStyle', 'none');
    errorbar([1:3], mean(congruency_er), std(congruency_er) ./ sqrt(p), 'LineStyle', 'none', 'Color', dark_colours(3,:), 'LineWidth', 1.5)
    xticks([1,2,3]);
    xticklabels(congruency_labels);
    title(['error as function of congruency']);
    % add individuals
    plot([1:3], congruency_er, 'Color', [0, 0, 0, 0.25]);
    
    %% effect of pressing
    % pressed or not
    % dt
    figure;
    hold on
    b = bar(mean(pressed_dt));
    errorbar([1:2], mean(pressed_dt), std(pressed_dt) ./ sqrt(p))
    xticks([1,2]);
    xticklabels(pressed_labels);
    title(['decision time as function of pressing']);
    % add individuals
    plot([1:2], pressed_dt);

    % er
    figure;
    hold on
    b = bar(mean(pressed_er));
    errorbar([1:2], mean(pressed_er), std(pressed_er) ./ sqrt(p))
    xticks([1,2]);
    xticklabels(pressed_labels);
    title(['error as function of pressing']);
    % add individuals
    plot([1:2], pressed_er);

    % pressed or not x should have pressed or not
    % dt
    figure;
    hold on
    b = bar(mean(cue_response_dt, "omitnan"));
    errorbar([1:4], mean(cue_response_dt, "omitnan"), std(cue_response_dt, "omitmissing") ./ sqrt(p))
    xticks([1:4]);
    xticklabels(cue_response_labels);
    title(['decision time as function of pressing x should have pressed']);
    % add individuals
    plot([1:4], cue_response_dt);

    % er
    figure;
    hold on
    b = bar(mean(cue_response_er, "omitnan"));
    errorbar([1:4], mean(cue_response_er, "omitnan"), std(cue_response_er, "omitmissing") ./ sqrt(p))
    xticks([1:4]);
    xticklabels(cue_response_labels);
    title(['error as function of pressing x should have pressed']);
    % add individuals
    plot([1:4], cue_response_er);

    %% congruency effect x pressing
    % dt
    figure; 
    hold on
    b = bar(mean(pressed_dt_effect));
    errorbar([1:2], mean(pressed_dt_effect), std(pressed_dt_effect) ./ sqrt(p))
    xticks([1:2]);
    xticklabels(pressed_labels);
    title(['congruency effect (dt) as function of press'])
    % add individuals
    plot([1:2], pressed_dt_effect);

    % er
    figure; 
    hold on
    b = bar(mean(pressed_er_effect));
    errorbar([1:2], mean(pressed_er_effect), std(pressed_er_effect) ./ sqrt(p))
    xticks([1:2]);
    xticklabels(pressed_labels);
    title(['congruency effect (er) as function of press'])
    % add individuals
    plot([1:2], pressed_er_effect);

    %% effect of block type
    % respond_3 vs. respond_not_3
    % dt
    figure;
    hold on
    b = bar(mean(block_dt));
    errorbar([1:2], mean(block_dt), std(block_dt) ./ sqrt(p))
    xticks([1,2]);
    xticklabels(block_labels);
    title(['decision time as function of block type']);
    % add individuals
    plot([1:2], block_dt);

    % er
    figure;
    hold on
    b = bar(mean(block_er));
    errorbar([1:2], mean(block_er), std(block_er) ./ sqrt(p))
    xticks([1,2]);
    xticklabels(block_labels);
    title(['error as function of block type']);
    % add individuals
    plot([1:2], block_er);

    %% congruency effect x block type
    % dt
    figure; 
    hold on
    b = bar(mean(block_dt_effect));
    errorbar([1:2], mean(block_dt_effect), std(block_dt_effect) ./ sqrt(p))
    xticks([1:2]);
    xticklabels(block_labels);
    title(['congruency effect (dt) as function of block type'])
    % add individuals
    plot([1:2], block_dt_effect);

    % er
    figure; 
    hold on
    b = bar(mean(block_er_effect));
    errorbar([1:2], mean(block_er_effect), std(block_er_effect) ./ sqrt(p))
    xticks([1:2]);
    xticklabels(block_labels);
    title(['congruency effect (er) as function of block type'])
    % add individuals
    plot([1:2], block_er_effect);

   
    %% main behavioural figure - pT
    % figure; 
    % hold on
    % b = bar([1,2], [mean(pT_effect(:,1)); mean(pT_effect(:,2))], 'LineStyle', 'none');
    % % add errorbars
    % errorbar([1:2], mean(pT_effect(:,1:2)), std(pT_effect(:,1:2)) ./ sqrt(p), "black", 'Color', dark_colours(1,:), 'LineWidth', 1.5);
    % 
    % % add individuals
    % plot([1:2], [pT_effect(:,1:2)]', 'Color', [0, 0, 0, 0.25]);
    % % highlight excluded participants (only works if included at the top)
    % % plot([x(1),x(2)], [congruency_er_effect(9,1:2)]', 'Color', [1, 0, 0, 0.25]);
    % % plot([x(3),x(4)], [congruency_er_effect(9,3:4)]', 'Color', [1, 0, 0, 0.25]);
    % % plot([x(1),x(2)], [congruency_er_effect(17,1:2)]', 'Color', [0, 0, 1, 0.25]);
    % % plot([x(3),x(4)], [congruency_er_effect(17,3:4)]', 'Color', [0, 0, 1, 0.25]);
    % xticks([1,2]);
    % xticklabels(probe_labels);
    % % ylim([-3.5 6.5]);
    % legend("location probe", "colour probe");
    % title(['pT effect - averaged']);
    % 
    % % set(gcf,'position',[0,0, 700,1080])

     %% main behavioural figure - pNT
    % figure; 
    % hold on
    % b = bar([1,2], [mean(pNT_effect(:,1)); mean(pNT_effect(:,2))], 'LineStyle', 'none');
    % % add errorbars
    % errorbar([1,2], mean(pNT_effect(:,1:2)), std(pNT_effect(:,1:2)) ./ sqrt(p), "black", 'Color', dark_colours(1,:), 'LineWidth', 1.5);
    % % add individuals
    % plot([1:2], [pNT_effect(:,1:2)]', 'Color', [0, 0, 0, 0.25]);
    % % highlight excluded participants (only works if included at the top)
    % % plot([x(1),x(2)], [congruency_er_effect(9,1:2)]', 'Color', [1, 0, 0, 0.25]);
    % % plot([x(3),x(4)], [congruency_er_effect(9,3:4)]', 'Color', [1, 0, 0, 0.25]);
    % % plot([x(1),x(2)], [congruency_er_effect(17,1:2)]', 'Color', [0, 0, 1, 0.25]);
    % % plot([x(3),x(4)], [congruency_er_effect(17,3:4)]', 'Color', [0, 0, 1, 0.25]);
    % xticks([1,2]);
    % xticklabels(probe_labels);
    % % ylim([-3.5 6.5]);
    % legend("location cue", "colour cue");
    % title(['pNT effect - averaged']);
    % 
    % % set(gcf,'position',[0,0, 700,1080])

    %% main figure for poster
    % settings
    w = [100 2.5];
    xpos = [1 2];
    xpos_offset = [1.1, 1.9];

    % decision time
    figure;
    tL = subplot(2,2,1);
    hold on 
    b = bar(mean(congruency_dt), 'FaceColor','flat', 'LineStyle', 'none');
    b.CData(1,:) = get_colour("blue", "");
    b.CData(2,:) = get_colour("green", "");
    b.CData(3,:) = get_colour("red", "");
    errorbar(1, mean(congruency_dt(:,1)), std(congruency_dt(:,1)) ./ sqrt(p), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("blue", "dark"))
    errorbar(2, mean(congruency_dt(:,2)), std(congruency_dt(:,2)) ./ sqrt(p), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("green", "dark"))
    errorbar(3, mean(congruency_dt(:,3)), std(congruency_dt(:,3)) ./ sqrt(p), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("red", "dark"))
    xticks([1,2,3]);
    xlim([0.35 3.65]);
    xticklabels(congruency_labels);
    ylabel('Decision time (ms)');
    yticks([0, 750, 1500]);
    % add individuals
    plot([1:3], congruency_dt, 'Color', [0, 0, 0, 0.25], 'LineWidth', 0.75);

    % calculate kernel densities
    [f_c, xi_c] = ksdensity(congruency_dt(:,1) - congruency_dt(:,2));
    [m,c] = min(abs(xi_c - mean(congruency_dt(:,1)-congruency_dt(:,2))));
    [f_i, xi_i] = ksdensity(congruency_dt(:,3) - congruency_dt(:,2));
    [m,i] = min(abs(xi_i - mean(congruency_dt(:,3)-congruency_dt(:,2))));

    tR = subplot(2,2,2);
    hold on 
    patch('XData', [f_c(:)*-w(1) + xpos(1), zeros(numel(xi_c(:)),1)],'yData', [xi_c(:),xi_c(:)],'FaceColor', get_colour("blue", ""), 'EdgeColor', 'none');
    patch('XData', [f_i(:)*w(1) + xpos(2), zeros(numel(xi_i(:)),1)],'yData', [xi_i(:),xi_i(:)],'FaceColor', get_colour("red", ""), 'EdgeColor', 'none');
    plot([xpos(1), xpos(1) + f_c(c)*-w(1)] , [mean(congruency_dt(:,1) - congruency_dt(:,2)), mean(congruency_dt(:,1) - congruency_dt(:,2))], 'Color', get_colour("blue", "dark"), 'LineWidth', 2);
    plot([xpos(2), xpos(2) + f_i(i)*w(1)] , [mean(congruency_dt(:,3) - congruency_dt(:,2)), mean(congruency_dt(:,3) - congruency_dt(:,2))], 'Color', get_colour("red", "dark"), 'LineWidth', 2);
    yline(0, 'LineWidth',1, 'LineStyle','--');
    xlim([0.25 2.75]);
    xticks([1 2]);
    xticklabels(congruency_labels([1,3]));
    yticks([-300, 0, 300]);
    ylim([-450 450]);
    % add individuals
    plot(xpos_offset, [congruency_dt(:,1)-congruency_dt(:,2),congruency_dt(:,3)-congruency_dt(:,2)], 'Color', [0, 0, 0, 0.25], 'LineWidth', 0.75);
    scatter(ones(25,1)*xpos_offset(1), congruency_dt(:,1)-congruency_dt(:,2), 20, get_colour("blue", ""), 'filled');
    scatter(ones(25,1)*xpos_offset(2), congruency_dt(:,3)-congruency_dt(:,2), 20,  get_colour("red", ""), 'filled');
    
    
    % error
    bL = subplot(2,2,3);
    hold on
    b = bar(mean(congruency_er), 'FaceColor','flat','LineStyle', 'none');
    b.CData(1,:) = get_colour("blue", "");
    b.CData(2,:) = get_colour("green", "");
    b.CData(3,:) = get_colour("red", "");
    errorbar(1, mean(congruency_er(:,1)), std(congruency_er(:,1)) ./ sqrt(p), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("blue", "dark"))
    errorbar(2, mean(congruency_er(:,2)), std(congruency_er(:,2)) ./ sqrt(p), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("green", "dark"))
    errorbar(3, mean(congruency_er(:,3)), std(congruency_er(:,3)) ./ sqrt(p), 'LineStyle', 'none', 'LineWidth', 1.5, 'Color', get_colour("red", "dark"))
    xticks([1,2,3]);
    xlim([0.35 3.65]);
    xticklabels(congruency_labels);
    ylabel('Reproduction error (°)');
    % add individuals
    plot([1:3], congruency_er, 'Color', [0, 0, 0, 0.25], 'LineWidth', 0.75);

    % calculate kernel densities
    [f_c, xi_c] = ksdensity(congruency_er(:,1) - congruency_er(:,2));
    [m,c] = min(abs(xi_c - mean(congruency_er(:,1)-congruency_er(:,2))));
    [f_i, xi_i] = ksdensity(congruency_er(:,3) - congruency_er(:,2));
    [m,i] = min(abs(xi_i - mean(congruency_er(:,3)-congruency_er(:,2))));

    bR = subplot(2,2,4);
    hold on 
    patch('XData', [f_c(:)*-w(2) + xpos(1), zeros(numel(xi_c(:)),1)],'yData', [xi_c(:),xi_c(:)],'FaceColor', get_colour("blue", ""), 'EdgeColor', 'none');
    patch('XData', [f_i(:)*w(2) + xpos(2), zeros(numel(xi_i(:)),1)],'yData', [xi_i(:),xi_i(:)],'FaceColor', get_colour("red", ""), 'EdgeColor', 'none');
    plot([xpos(1), xpos(1) + f_c(c)*-w(2)] , [mean(congruency_er(:,1) - congruency_er(:,2)), mean(congruency_er(:,1) - congruency_er(:,2))], 'Color', get_colour("blue", "dark"), 'LineWidth', 2);
    plot([xpos(2), xpos(2) + f_i(i)*w(2)] , [mean(congruency_er(:,3) - congruency_er(:,2)), mean(congruency_er(:,3) - congruency_er(:,2))], 'Color', get_colour("red", "dark"), 'LineWidth', 2);
    yline(0, 'LineWidth',1, 'LineStyle','--');
    xlim([0.25 2.75]);
    xticks([1 2]);
    xticklabels(congruency_labels([1,3]));
    yticks([-10, 0, 10]);
    ylim([-15 15]);
    % add individuals
    plot(xpos_offset, [congruency_er(:,1)-congruency_er(:,2),congruency_er(:,3)-congruency_er(:,2)], 'Color', [0, 0, 0, 0.25], 'LineWidth', 0.75);
    scatter(ones(25,1)*xpos_offset(1), congruency_er(:,1)-congruency_er(:,2), 20, get_colour("blue", ""), 'filled');
    scatter(ones(25,1)*xpos_offset(2), congruency_er(:,3)-congruency_er(:,2), 20,  get_colour("red", ""), 'filled');  
    
    % general
    set(gcf(), 'Position', [800 300 880 800]);

    axes = {tL, tR, bL, bR};
    for i = 1:size(axes,2)
        set(axes{i}, 'Box', 'on');
        set(axes{i}, 'FontSize', [24.8]);
        set(axes{i}, 'FontName', 'Aptos');
        set(axes{i}.XAxis, 'FontSize', 18);
        set(axes{i}, 'LineWidth', 1.33);
    end

    bL.YLabel.Position = [-0.6417   20.0000   -1.0000];

    %% do some stats over previous figure
    n = size(pp2do,2);
    dt_tbl = table([1:n,1:n,1:n]', [congruency_dt(:,1);congruency_dt(:,2);congruency_dt(:,3)],[ones(n,1);ones(n,1)*2;ones(n,1)*3], 'VariableNames', {'participant_id', 'decision_time', 'condition'});
    er_tbl = table([1:n,1:n,1:n]', [congruency_er(:,1);congruency_er(:,2);congruency_er(:,3)],[ones(n,1);ones(n,1)*2;ones(n,1)*3], 'VariableNames', {'participant_id', 'error', 'condition'});

    dt_lme = fitlme(dt_tbl,'decision_time ~ condition + (1|participant_id)')
    anova(dt_lme)
    er_lme = fitlme(er_tbl,'error ~ condition + (1|participant_id)')
    anova(er_lme)

    %post-hoc pairwise comparison
    p_values = zeros(3,2)
    [h,p_values(1,1),ci,stats] = ttest(congruency_dt(:,1), congruency_dt(:,2));
    [h,p_values(2,1),ci,stats] = ttest(congruency_dt(:,2), congruency_dt(:,3));
    [h,p_values(3,1),ci,stats] = ttest(congruency_dt(:,1), congruency_dt(:,3));

    [h,p_values(1,2),ci,stats] = ttest(congruency_er(:,1), congruency_er(:,2));
    [h,p_values(2,2),ci,stats] = ttest(congruency_er(:,2), congruency_er(:,3));
    [h,p_values(3,2),ci,stats] = ttest(congruency_er(:,1), congruency_er(:,3));
   
end

%% see performance per block (for fun for participants)
if see_performance
    x = (1:16);
    block_performance = zeros(1, length(x));
    block_performance_std = zeros(1, length(x));
    block_speed = zeros(1, length(x));
    block_speed_std = zeros(1, length(x));
    
    for i = x
        block_performance(i) = mean(behdata.performance(behdata.block == i));
        block_performance_std(i) = std(behdata.performance(behdata.block == i));
        block_speed(i) = mean(behdata.idle_reaction_time_in_ms(behdata.block == i));
        block_speed_std(i) = std(behdata.idle_reaction_time_in_ms(behdata.block == i));
    end
    
    figure;
    hold on
    plot(block_performance)
    %errorbar(block_performance, block_performance_std)
    ylim([50 100])
    xlim([1 16])
    ylabel('Performance score')
    yyaxis right
    plot(block_speed)
    %errorbar(block_speed, block_speed_std)
    ylim([100 2000])
    ylabel('Reaction time (ms)')
    xlabel('Block number')
    xticks(x)
end
