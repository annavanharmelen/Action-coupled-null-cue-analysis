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
    statcfg.nsub = size(pp2do, 2);
    statcfg.statMethod = 'montecarlo';

    timeframe = [951:2451]; %0 - 1500 ms post-cue

    stat_action = frevede_ftclusterstat1D(statcfg, action_data(:,timeframe), zeros(size(action_data(:,timeframe))));
    stat_no_action = frevede_ftclusterstat1D(statcfg, no_action_data(:,timeframe), zeros(size(no_action_data(:,timeframe))));
    stat_compare = frevede_ftclusterstat1D(statcfg, action_data(:,timeframe), no_action_data(:,timeframe));

    mask_action = double(stat_action.mask); mask_action(mask_action==0) = nan; % nan data that is not part of mark
    mask_no_action = double(stat_no_action.mask); mask_no_action(mask_no_action==0) = nan;
    mask_compare = double(stat_compare.mask); mask_compare(mask_compare==0) = nan;

    % sig1 = plot(saccade.time(timeframe), mask_action*-0.16, 'Color', get_colour("pink",""), 'LineWidth', 3); %shown already in first figure
    % sig2 = plot(saccade.time(timeframe), mask_no_action*-0.165, 'Color', [0.5,0.5,0.5], 'LineWidth', 3); % not significant
    sig3 = plot(saccade.time(timeframe), mask_compare*-0.14, 'Color', 'k', 'LineWidth', 6);

    % print("C:\Users\annav\Documents\Surfdrive\Conferences\ICON 2025\Figures\saccade_effect_of_action", "-dsvg")
    % print("C:\Users\annav\Documents\Surfdrive\Conferences\ICON 2025\Figures\saccade_effect_of_action", "-dpng")

end
