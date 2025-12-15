%% Script for comparing the saccade results from the action-coupled cue to a non-action one

%% start clean
clear; clc; close all;
    
%% parameters
pp2do           = [[1:9,11:26];[1:8, 10:16, 18:27];[1:25];[3:10,12:17,19:23,25:30]];

oneOrTwoD       = 1;        oneOrTwoD_options = {'_1D','_2D'};
nsmooth         = 200;
plotSinglePps   = 0;
plotGAs         = 0;
plotFigures     = 0;
xlimtoplot      = [-100 1500];

%% load and aggregate the data from each study
% this all only works so easily because all the data was time-locked to the
% cue in the same manner (from 1s before to 2.5s after)
for s = [1:4]
    disp([newline(), 'getting data from study ', num2str(s)]);
    p = 0;

    for pp = pp2do(s,:)
        p = p+1;
    
        % get participant data
        param = getAllSubjParam(s, pp);
    
        % load
        disp(['getting data from participant ', param.subjName]);
        load([param.path, '\saved_data\saccadeEffects', oneOrTwoD_options{oneOrTwoD} '__', param.subjName], 'saccade','saccadesize');
           
        % smooth?
        if nsmooth > 0
            for x1 = 1:size(saccade.toward,1)
                saccade.toward(x1,:)  = smoothdata(squeeze(saccade.toward(x1,:)), 'gaussian', nsmooth);
                saccade.away(x1,:)    = smoothdata(squeeze(saccade.away(x1,:)), 'gaussian', nsmooth);
                saccade.effect(x1,:)  = smoothdata(squeeze(saccade.effect(x1,:)), 'gaussian', nsmooth);
            end
            % also smooth saccadesize data over time.
            for x1 = 1:size(saccadesize.toward,1)
                for x2 = 1:size(saccadesize.toward,2)
                    saccadesize.toward(x1,x2,:) = smoothdata(squeeze(saccadesize.toward(x1,x2,:)), 'gaussian', nsmooth);
                    saccadesize.away(x1,x2,:)   = smoothdata(squeeze(saccadesize.away(x1,x2,:)), 'gaussian', nsmooth);
                    saccadesize.effect(x1,x2,:) = smoothdata(squeeze(saccadesize.effect(x1,x2,:)), 'gaussian', nsmooth);
                end
            end
        end
    
        % put into matrix, with study as first and participant as second dimension
        d1(s,p,1:size(saccade.toward,1),:) = saccade.toward;
        d2(s,p,1:size(saccade.away,1),:) = saccade.away;
        d3(s,p,1:size(saccade.effect,1),:) = saccade.effect;
    
        d4(s,p,1:size(saccadesize.toward,1),:,:) = saccadesize.toward;
        d5(s,p,1:size(saccadesize.away,1),:,:) = saccadesize.away;
        d6(s,p,1:size(saccadesize.effect,1),:,:) = saccadesize.effect;

        % save average saccade bias per pp
        if s == 1
            avg_saccade_effect(1,p,s) = mean(saccade.effect(5,saccade.time>=200 & saccade.time<=600));
            avg_saccade_effect(2,p,s) = mean(saccade.effect(5,:));
        elseif s == 4
            avg_saccade_effect(1,p,s) = mean(saccade.effect(7,saccade.time>=200 & saccade.time<=600));
            avg_saccade_effect(2,p,s) = mean(saccade.effect(7,:));
        else
            avg_saccade_effect(1,p,s) = mean(saccade.effect(4,saccade.time>=200 & saccade.time<=600));
            avg_saccade_effect(2,p,s) = mean(saccade.effect(4,:));
        end

        % save average saccade bias per pp
        if s == 1
            max_saccade_effect(1,p,s) = max(saccade.effect(5,saccade.time>=200 & saccade.time<=600));
            max_saccade_effect(2,p,s) = max(saccade.effect(5,:));
        elseif s == 4
            max_saccade_effect(1,p,s) = max(saccade.effect(7,saccade.time>=200 & saccade.time<=600));
            max_saccade_effect(2,p,s) = max(saccade.effect(7,:));
        else
            max_saccade_effect(1,p,s) = max(saccade.effect(4,saccade.time>=200 & saccade.time<=600));
            max_saccade_effect(2,p,s) = max(saccade.effect(4,:));
        end
    end
end

% make GA for the saccadesize fieldtrip structure data (update code below,
% because this averages over studies, not pp's)
% saccadesize.toward = squeeze(mean(d4));
% saccadesize.away   = squeeze(mean(d5));
% saccadesize.effect = squeeze(mean(d6));

if plotFigures
    %% Compare all 4 studies
    figure;
    hold on
    p1 = frevede_errorbarplot(saccade.time, squeeze(d3(1,:,5,:)), [1,0,0], 'se');
    p2 = frevede_errorbarplot(saccade.time, squeeze(d3(2,:,4,:)), [0,0,1], 'se');
    p3 = frevede_errorbarplot(saccade.time, squeeze(d3(3,:,4,:)), [0,1,1], 'se');
    p4 = frevede_errorbarplot(saccade.time, squeeze(d3(4,:,7,:)), [1,0,1], 'se');
    xlim(xlimtoplot);
    legend([p1, p2, p3, p4], {'study 1', 'study 2', 'study 3', 'study 4'}, 'EdgeColor', 'none', 'AutoUpdate','off', 'FontSize', 25.4);

    %% Effect of having a secondary task for the cue
    action_data = [squeeze(d3(1,:,5,:)); squeeze(d3(4,:,7,:))];
    no_action_data = [squeeze(d3(2,:,4,:)); squeeze(d3(3,:,4,:))];

    figure;
    hold on
    p1 = frevede_errorbarplot(saccade.time, action_data, get_colour("pink",""), 'se');
    p2 = frevede_errorbarplot(saccade.time, no_action_data, [0.5,0.5,0.5], 'se');
    p1.LineWidth = 3;
    p2.LineWidth = 3;
    xline(0, 'LineWidth',3, 'Color',[107, 107, 107]/255, 'LineStyle','--');
    yline(0, 'LineWidth',3, 'Color',[107, 107, 107]/255, 'LineStyle','--');
    xlim(xlimtoplot);
    ylim([-0.2 0.4]);
    yticks([-0.2 0 0.2 0.4]);
    ylabel('Saccade bias (ΔHz)');
    xlabel('Time after cue (ms)');
    legend([p1, p2], {'cue action required', 'no action required'}, 'EdgeColor', 'none', 'AutoUpdate','off', 'FontSize', 25.4);
    set(gca(), 'FontSize', [25.4]);
    set(gca(), 'FontName', 'Aptos');
    set(gca(), 'Box', 'on');
    set(gca(), 'LineWidth', 1.49);
    set(gcf(), 'Position', [500 500 850 800])

    % add stats
    statcfg.xax = saccade.time;
    statcfg.npermutations = 10000;
    statcfg.clusterStatEvalaluationAlpha = 0.05;
    statcfg.nsub = 50; %for within pp
    statcfg.nsub1 = 50; %for between pp
    statcfg.nsub2 = 50; %for between pp
    statcfg.statMethod = 'montecarlo';

    timeframe = [951:2451]; %0 - 1500 ms post-cue

    stat_action = frevede_ftclusterstat1D(statcfg, action_data(:,timeframe), zeros(size(action_data(:,timeframe))));
    stat_no_action = frevede_ftclusterstat1D(statcfg, no_action_data(:,timeframe), zeros(size(no_action_data(:,timeframe))));
    stat_compare = frevede_ftclusterstat1D_indep(statcfg, action_data(:,timeframe), no_action_data(:,timeframe));

    mask_action = double(stat_action.mask); mask_action(mask_action==0) = nan; % nan data that is not part of mark
    mask_no_action = double(stat_no_action.mask); mask_no_action(mask_no_action==0) = nan;
    mask_compare = double(stat_compare.mask); mask_compare(mask_compare==0) = nan;

    sig1 = plot(saccade.time(timeframe), mask_action*-0.16, 'Color', get_colour("pink",""), 'LineWidth', 4); %shown already in first figure
    sig2 = plot(saccade.time(timeframe), mask_no_action*-0.15, 'Color', [0.5,0.5,0.5], 'LineWidth', 4); % not significant
    sig3 = plot(saccade.time(timeframe), mask_compare*-0.14, 'Color', 'k', 'LineWidth', 4);

    % print("C:\Users\annav\Documents\Surfdrive\Conferences\ICON 2025\Figures\saccade_effect_of_action", "-dsvg")
    % print("C:\Users\annav\Documents\Surfdrive\Conferences\ICON 2025\Figures\saccade_effect_of_action", "-dpng")

    %% Interaction between microsaccade bias and RT (combined E1 and E4)
    figure;
    hold on;
    p1 = frevede_errorbarplot(saccade.time, squeeze(d3(1,:,12,:)), 'r', 'se'); % fast congruent
    p2 = frevede_errorbarplot(saccade.time, squeeze(d3(1,:,16,:)), 'b', 'se'); % slow congruent
    p3 = frevede_errorbarplot(saccade.time, squeeze(d3(1,:,13,:)), 'r', 'se'); % fast incongruent
    p4 = frevede_errorbarplot(saccade.time, squeeze(d3(1,:,17,:)), 'b', 'se'); % slow incongruent
    legend([p1,p2,p3,p4], {'congruent fast', 'congruent slow', 'incongruent fast', 'incongruent slow'});
    title('E1');
    xlim(xlimtoplot);

    figure;
    hold on;
    p1 = frevede_errorbarplot(saccade.time, squeeze(d3(4,:,14,:)), 'r', 'se'); % fast congruent
    p2 = frevede_errorbarplot(saccade.time, squeeze(d3(4,:,18,:)), 'b', 'se'); % slow congruent
    p3 = frevede_errorbarplot(saccade.time, squeeze(d3(4,:,15,:)), 'r', 'se'); % fast incongruent
    p4 = frevede_errorbarplot(saccade.time, squeeze(d3(4,:,19,:)), 'b', 'se'); % slow incongruent
    legend([p1,p2,p3,p4], {'congruent fast', 'congruent slow', 'incongruent fast', 'incongruent slow'});
    title('E2');
    xlim(xlimtoplot);

    figure;
    hold on;
    p1 = frevede_errorbarplot(saccade.time, [squeeze(d3(1,:,12,:));squeeze(d3(4,:,14,:))], 'r', 'se'); % fast congruent
    p2 = frevede_errorbarplot(saccade.time, [squeeze(d3(1,:,16,:));squeeze(d3(4,:,18,:))], 'b', 'se'); % slow congruent
    p3 = frevede_errorbarplot(saccade.time, [squeeze(d3(1,:,13,:));squeeze(d3(4,:,15,:))], 'm', 'se'); % fast incongruent
    p4 = frevede_errorbarplot(saccade.time, [squeeze(d3(1,:,17,:));squeeze(d3(4,:,19,:))], 'c', 'se'); % slow incongruent
    legend([p1,p2,p3,p4], {'congruent fast', 'congruent slow', 'incongruent fast', 'incongruent slow'});
    title('E1 + E2 combined');
    xlim(xlimtoplot);

    high_capture_data = [squeeze(d3(1,:,12,:));squeeze(d3(4,:,14,:));squeeze(d3(1,:,17,:));squeeze(d3(4,:,19,:))];
    low_capture_data = [squeeze(d3(1,:,16,:));squeeze(d3(4,:,18,:));squeeze(d3(1,:,13,:));squeeze(d3(4,:,15,:))];
    figure;
    hold on;
    p1 = frevede_errorbarplot(saccade.time, high_capture_data, 'r', 'se'); % fast congruent + slow incongruent
    p2 = frevede_errorbarplot(saccade.time, low_capture_data, 'b', 'se'); % slow congruent + fast incongruent
    legend([p1,p2], {'congruent fast + incongruent slow', 'congruent slow + incongruent fast'}, 'AutoUpdate','off');
    title('E1 + E2 combined');
    xlim(xlimtoplot);
    % add stats
    statcfg.xax = saccade.time;
    statcfg.npermutations = 10000;
    statcfg.clusterStatEvalaluationAlpha = 0.05;
    statcfg.nsub = size(pp2do, 2) * 2;
    statcfg.statMethod = 'montecarlo';

    timeframe = [951:2451]; %0 - 1500 ms post-cue

    stat_high_capture = frevede_ftclusterstat1D(statcfg, high_capture_data(:,timeframe), zeros(size(high_capture_data(:,timeframe))));
    stat_low_capture = frevede_ftclusterstat1D(statcfg, low_capture_data(:,timeframe), zeros(size(low_capture_data(:,timeframe))));
    stat_comp_capture = frevede_ftclusterstat1D(statcfg, high_capture_data(:,timeframe), low_capture_data(:,timeframe));

    mask_high_cap = double(stat_high_capture.mask); mask_high_cap(mask_high_cap==0) = nan; % nan data that is not part of mark
    mask_low_cap = double(stat_low_capture.mask); mask_low_cap(mask_low_cap==0) = nan;
    mask_comp_cap = double(stat_comp_capture.mask); mask_comp_cap(mask_comp_cap==0) = nan;

    sig1 = plot(saccade.time(timeframe), mask_high_cap*-0.16, 'r', 'LineWidth', 3);
    sig2 = plot(saccade.time(timeframe), mask_low_cap*-0.165, 'b', 'LineWidth', 3);
    sig3 = plot(saccade.time(timeframe), mask_comp_cap*-0.14, 'k', 'LineWidth', 6); %not significant, p=0.138 for first cluster (during peak)

    fast_data = [squeeze(d3(1,:,12,:));squeeze(d3(4,:,14,:))]-[squeeze(d3(1,:,13,:));squeeze(d3(4,:,15,:))];
    slow_data = [squeeze(d3(1,:,16,:));squeeze(d3(4,:,18,:))]-[squeeze(d3(1,:,17,:));squeeze(d3(4,:,19,:))];
    figure;
    % this one makes no sense no?
    hold on;
    p1 = frevede_errorbarplot(saccade.time, fast_data, 'r', 'se'); % fast congruent - fast incongruent
    p2 = frevede_errorbarplot(saccade.time, slow_data, 'b', 'se'); % slow congruent - slow incongruent
    legend([p1,p2], {'congruent fast - incongruent fast', 'congruent slow - incongruent slow'}, 'AutoUpdate','off');
    title('E1 + E2 combined');
    xlim(xlimtoplot);
    % add stats
    statcfg.xax = saccade.time;
    statcfg.npermutations = 10000;
    statcfg.clusterStatEvalaluationAlpha = 0.05;
    statcfg.nsub = size(pp2do, 2) * 2;
    statcfg.statMethod = 'montecarlo';

    timeframe = [951:2451]; %0 - 1500 ms post-cue

    stat_fast = frevede_ftclusterstat1D(statcfg, fast_data(:,timeframe), zeros(size(fast_data(:,timeframe))));
    stat_slow = frevede_ftclusterstat1D(statcfg, slow_data(:,timeframe), zeros(size(slow_data(:,timeframe))));
    stat_comp_speed = frevede_ftclusterstat1D(statcfg, fast_data(:,timeframe), slow_data(:,timeframe));

    mask_fast = double(stat_fast.mask); mask_fast(mask_fast==0) = nan; % nan data that is not part of mark
    mask_slow = double(stat_slow.mask); mask_slow(mask_slow==0) = nan;
    mask_comp_speed = double(stat_comp_speed.mask); mask_comp_speed(mask_comp_speed==0) = nan;

    sig1 = plot(saccade.time(timeframe), mask_fast*-0.16, 'r', 'LineWidth', 3);
    sig2 = plot(saccade.time(timeframe), mask_slow*-0.165, 'b', 'LineWidth', 3);
    sig3 = plot(saccade.time(timeframe), mask_comp_speed*-0.14, 'k', 'LineWidth', 6); %not significant, p=0.138 for first cluster (during peak)

end
