
%% Step3b--grand average plots of gaze-shift (saccade) results

%% start clean
clear; clc; close all;
    
%% parameters
pp2do           = [1:9,11:26];
oneOrTwoD       = 1;        oneOrTwoD_options = {'_1D','_2D'};
nsmooth         = 200;
plotSinglePps   = 0;
plotGAs         = 0;
plotFigures     = 0;
xlimtoplot      = [-100 1500];

%% set visual parameters
[bar_size, colours, dark_colours, labels, subplot_size, percentageok] = setBehaviourParam(pp2do);

avg_effect_labels = {'fast congruent', 'slow congruent', 'fast incongruent', 'slow incongruent'};

%% load and aggregate the data from all pp
s = 0;
for pp = pp2do
    s = s+1;

    % get participant data
    param = getSubjParam(pp);

    % load
    disp(['getting data from participant ', param.subjName]);
    load([param.path, '\saved_data\saccadeEffects', oneOrTwoD_options{oneOrTwoD} '__', param.subjName], 'saccade','saccadesize');

    % save averages (saccade effect (capture cue effect and probe cue reaction)
    avg_saccade_effect(s, 1) = mean(saccade.effect(12,saccade.time>=200 & saccade.time<=600));
    avg_saccade_effect(s, 2) = mean(saccade.effect(16,saccade.time>=200 & saccade.time<=600));
    avg_saccade_effect(s, 3) = mean(saccade.effect(13,saccade.time>=200 & saccade.time<=600));
    avg_saccade_effect(s, 4) = mean(saccade.effect(17,saccade.time>=200 & saccade.time<=600));
       
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

    % put into matrix, with pp as first dimension
    d1(s,:,:) = saccade.toward;
    d2(s,:,:) = saccade.away;
    d3(s,:,:) = saccade.effect;

    d4(s,:,:,:) = saccadesize.toward;
    d5(s,:,:,:) = saccadesize.away;
    d6(s,:,:,:) = saccadesize.effect;
end

%% make GA for the saccadesize fieldtrip structure data, to later plot as "time-frequency map" with fieldtrip. For timecourse data, we directly plot from d structures above. 
saccadesize.toward = squeeze(mean(d4));
saccadesize.away   = squeeze(mean(d5));
saccadesize.effect = squeeze(mean(d6));

%% all subs
if plotSinglePps
    % effect
    figure;
    for sp = 1:s
        subplot(subplot_size,subplot_size,sp); hold on;
        plot(saccade.time, squeeze(d3(sp,5,:)));
        plot(xlim, [0,0], '--k');
        xlim(xlimtoplot); ylim([-1.5 1.5]);
        title(pp2do(sp));
    end
    legend(saccade.label);

    % toward vs. away
    figure;
    for sp = 1:s
        subplot(subplot_size,subplot_size,sp); hold on;
        plot(saccade.time, squeeze(d3(sp,:,:)));
        plot(xlim, [0,0], '--k');
        xlim(xlimtoplot); ylim([-1.5 1.5]);
        title(pp2do(sp));
    end
    legend(saccade.label);

    % towardness for all conditions condition - gaze shift effect X saccade size
    figure;
    cfg = [];
    cfg.parameter = 'effect_individual';
    cfg.figure = 'gcf';
    cfg.zlim = [-.1 .1];
    cfg.xlim = xlimtoplot;
    for sp = 1:s
        subplot(subplot_size,subplot_size,sp); hold on;
        saccadesize.effect_individual = squeeze(d6(sp,:,:,:)); % put in data from this pp
        cfg.channel = 2; % congruent
        ft_singleplotTFR(cfg, saccadesize);
        title(pp2do(sp));
    end
    colormap('jet');

    figure;
    cfg = [];
    cfg.parameter = 'effect_individual';
    cfg.figure = 'gcf';
    cfg.zlim = [-.1 .1];
    cfg.xlim = xlimtoplot;
    for sp = 1:s
        subplot(subplot_size,subplot_size,sp); hold on;
        saccadesize.effect_individual = squeeze(d6(sp,:,:,:)); % put in data from this pp
        cfg.channel = 3; % incongruent
        ft_singleplotTFR(cfg, saccadesize);
        title(pp2do(sp));
    end
    colormap('jet');
end

%% plot grand average data patterns of interest, with error bars
if plotGAs
    % right and left cues, per condition
    figure;
    for sp = [1:4]
        subplot(2,4,sp); hold on; title(saccade.label(sp));
        p1 = frevede_errorbarplot(saccade.time, squeeze(d1(:,sp,:)), [1,0,0], 'se');
        p2 = frevede_errorbarplot(saccade.time, squeeze(d2(:,sp,:)), [0,0,1], 'se');
        plot(xlim, [0,0], '--k');
        xlim(xlimtoplot); ylim([0 1]);
    end
    legend([p1, p2], {'toward','away'});
    
    % towardness per condition - gaze shift effect X saccade size
    figure;
    for sp = [1:4]
        subplot(2,4,sp); hold on; title(saccade.label(sp));
        frevede_errorbarplot(saccade.time, squeeze(d3(:,sp,:)), [0,0,0], 'both');
        plot(xlim, [0,0], '--k');
        xlim(xlimtoplot); ylim([-0.3 0.3]);
    end
    legend({'effect'});
    
    %% towardness overlay of all conditions
    ylimit = [-0.3, 0.3];
    figure;
    hold on;
    p1 = frevede_errorbarplot(saccade.time, squeeze(d3(:,1,:)), 'k', 'both');
    plot(xlim, [0,0], '--k');
    plot([0,0], ylimit, '--k');
    legend([p1], saccade.label(1));
    xlim(xlimtoplot);
    ylabel('Rate (Hz)');
    xlabel('Time (ms)');
    
    %% effect congruent - incongruent
    ylimit = [-0.3, 0.3];
    figure;
    hold on;
    p1 = frevede_errorbarplot(saccade.time, squeeze((d3(:,2,:) - d3(:,3,:)) ./ 2), 'k', 'both'); %should be same as cue-match condition (nr. 5)
    plot(xlim, [0,0], '--k');
    plot([0,0], ylimit, '--k');
    legend([p1], saccade.label(1));
    xlim(xlimtoplot);
    ylabel('Rate (Hz)');
    xlabel('Time (ms)');

    %% effect of pressing to cue
    ylimit = [-0.3, 0.3];
    figure;
    hold on;
    p1 = frevede_errorbarplot(saccade.time, squeeze(d3(:,8,:)), 'b', 'both');
    p2 = frevede_errorbarplot(saccade.time, squeeze(d3(:,9,:)), 'r', 'both');
    plot(xlim, [0,0], '--k');
    plot([0,0], ylimit, '--k');
    legend([p1, p2], {'pressed', 'did not press'});
    xlim(xlimtoplot);
    ylabel('Rate (Hz)');
    xlabel('Time (ms)');

    %% effect of having to press to cue
    ylimit = [-0.3, 0.3];

    figure;
    hold on;
    p1 = frevede_errorbarplot(saccade.time, squeeze(d3(:,10,:)), 'b', 'both');
    p2 = frevede_errorbarplot(saccade.time, squeeze(d3(:,11,:)), 'r', 'both');
    plot(xlim, [0,0], '--k');
    plot([0,0], ylimit, '--k');
    legend([p1, p2], {'should have pressed', 'should not have pressed'});
    xlim(xlimtoplot);
    ylabel('Rate (Hz)');
    xlabel('Time (ms)');

    %% effect of block type
    figure;
    hold on;
    p1 = frevede_errorbarplot(saccade.time, squeeze(d3(:,20,:)), [1,0,0], 'both');
    p2 = frevede_errorbarplot(saccade.time, squeeze(d3(:,21,:)), [0,0,1], 'both');
    plot(xlim, [0,0], '--k');
    plot([0,0], ylimit, '--k');
    legend([p1, p2], {'respond_3 block', 'respond_not_3 block'});
    xlim(xlimtoplot);
    ylabel('Rate (Hz)');
    xlabel('Time (ms)');
    title('Effect of block type')

    %% as function of saccade size
    cfg = [];
    cfg.parameter = 'effect';
    cfg.figure = 'gcf';
    cfg.zlim = [-0.1, 0.1];
    cfg.xlim = xlimtoplot;
    cfg.colormap = 'jet';
    % per condition
    figure;
    for chan = [1:4]
        cfg.channel = chan;
        subplot(2,4,chan); ft_singleplotTFR(cfg, saccadesize);
    end
    % cfg.channel = 4;
    % ft_singleplotTFR(cfg, saccadesize);
    ylabel('Saccade size (dva)')
    xlabel('Time (ms)')
    hold on
    plot([0,0], [0, 7], '--k');
    % plot([1500,1500], [0, 7], '--', 'LineWidth',3, 'Color', [0.6, 0.6, 0.6]);
    ylim([0.2 6.8]);

    %% plot effect of performance (median-split) on saccade data
    figure;
    subplot(2,2,1)
    hold on
    c1 = frevede_errorbarplot(saccade.time, squeeze(d3(:,16,:)), [0.5, 0.5, 0.5], 'se');
    c2 = frevede_errorbarplot(saccade.time, squeeze(d3(:,12,:)), 'b', 'se');
    legend([c1, c2], {'above median DT', 'below median DT'});
    xlim(xlimtoplot);
    title('Congruent');
    ylabel('Saccade bias (Hz)');

    subplot(2,2,2)
    hold on
    i1 = frevede_errorbarplot(saccade.time, squeeze(d3(:,17,:)), [0.5, 0.5, 0.5], 'se');
    i2 = frevede_errorbarplot(saccade.time, squeeze(d3(:,13,:)), 'r', 'se');
    legend([i1, i2], {'above median DT', 'below median DT'});
    xlim(xlimtoplot);
    title('Incongruent');

    subplot(2,2,3)
    hold on
    c3 = frevede_errorbarplot(saccade.time, squeeze(d3(:,18,:)), [0.5, 0.5, 0.5], 'se');
    c4 = frevede_errorbarplot(saccade.time, squeeze(d3(:,14,:)), 'b', 'se');
    legend([c3, c4], {'above median ER', 'below median ER'});
    xlim(xlimtoplot);
    ylabel('Saccade bias (Hz)');

    subplot(2,2,4)
    hold on
    i3 = frevede_errorbarplot(saccade.time, squeeze(d3(:,19,:)), [0.5, 0.5, 0.5], 'se');
    i4 = frevede_errorbarplot(saccade.time, squeeze(d3(:,15,:)), 'r', 'se');
    legend([i3, i4], {'above median ER', 'below median ER'});
    xlim(xlimtoplot);

    %% plot fast vs. slow (median-split) on saccade data
    % as timecourse
    figure;
    hold on
    j1 = frevede_errorbarplot(saccade.time, squeeze(d3(:,12,:)) - squeeze(d3(:,16,:)), 'b', 'se');
    j2 = frevede_errorbarplot(saccade.time, squeeze(d3(:,13,:)) - squeeze(d3(:,17,:)), 'r', 'se');
    ylabel('Saccade bias (Hz)');
    xlabel('Time (ms)');
    legend([j1, j2], {'Congruent', 'Incongruent'});
    xlim(xlimtoplot);
    title('Fast vs. slow effect');

    % as bargraph
    figure;
    hold on
    b1 = bar([1,2], mean(avg_saccade_effect(:,1:2)), 0.6);
    b2 = bar([3,4], mean(avg_saccade_effect(:,3:4)), 0.6);
    errorbar(mean(avg_saccade_effect), std(avg_saccade_effect) ./ sqrt(size(pp2do,2)), 'LineStyle', 'None', 'Color', 'k');
    ylabel('Mean saccade rate (ΔHz)');
    xticks([1,2,3,4]);
    xticklabels({'fast', 'slow', 'fast', 'slow'});
    legend({'congruent', 'incongruent'});
    
    % first level = condition, second level = fast or slow
    stats = rm_anova2([avg_saccade_effect(:,1);avg_saccade_effect(:,2);avg_saccade_effect(:,3);avg_saccade_effect(:,4)], repmat([1:25]',4,1), [ones(50,1);ones(50,1)*2], [ones(25,1);ones(25,1)*2;ones(25,1);ones(25,1)*2], {'condition',  'speed'})
    
    [h,p,ci,stats] = ttest(avg_saccade_effect(:,1), avg_saccade_effect(:,2))
    [h,p,ci,stats] = ttest(avg_saccade_effect(:,3), avg_saccade_effect(:,4))
end
%% main figures for poster
if plotFigures
    %% saccade bias post-cue, including size course
    %some settings
    ylim1 = [0, 1];
    ylim2 = [-0.2, 0.4];
    line = 2;
    x_ticks = [0, 500, 1000, 1500];

    cfg = [];
    cfg.parameter = 'effect';
    cfg.figure = 'gcf';
    cfg.zlim = [-0.12, 0.12];
    cfg.xlim = xlimtoplot;
    cfg.colormap = brewermap(1000, 'RdBu');

    statcfg.xax = saccade.time;
    statcfg.npermutations = 10000;
    statcfg.clusterStatEvalaluationAlpha = 0.05;
    statcfg.nsub = size(pp2do, 2);
    statcfg.statMethod = 'montecarlo';
    
    timeframe = [951:2451]; %0 - 1500 ms post-cue
    stat = frevede_ftclusterstat1D(statcfg, squeeze(d3(:,5,timeframe)), zeros(size(squeeze(d3(:,5,timeframe)))));
    mask = double(stat.mask); mask(mask==0) = nan; % nan data that is not part of mark
    
   %%
    % make figure 🎉
    figure;

    % cue-matching
    tL = subplot(1,3,1);
    hold on
    p1 = frevede_errorbarplot(saccade.time, squeeze(d1(:,5,:)), get_colour("blue", ""), 'se');
    p2 = frevede_errorbarplot(saccade.time, squeeze(d2(:,5,:)), get_colour("red", ""), 'se');
    p1.LineWidth = line;
    p2.LineWidth = line;
    legend([p1, p2], {'toward', 'away'}, 'EdgeColor','none', 'AutoUpdate','off', 'FontSize', 17);
    xline(0, 'LineWidth',2, 'Color',[107, 107, 107]/255, 'LineStyle','--');
    xlim(xlimtoplot);
    xticks(x_ticks);
    ylim(ylim1);

    yticks([0, 0.5, 1]);
    ylabel('Saccade rate (Hz)');

    mL = subplot(1,3,2);
    hold on
    p3 = frevede_errorbarplot(saccade.time, squeeze(d3(:,5,:)), get_colour("pink", ""), 'se');
    p3.LineWidth = line;
    yline(0, 'LineWidth',2, 'Color',[107, 107, 107]/255, 'LineStyle','--');
    xline(0, 'LineWidth',2, 'Color',[107, 107, 107]/255, 'LineStyle','--');
    sig = plot(saccade.time(timeframe), mask*-0.15, 'Color', get_colour("pink",""), 'LineWidth', 3);
    xlim(xlimtoplot);
    xticks(x_ticks);
    ylim(ylim2);
    yticks([-0.2 0 0.2 0.4]);
    ylabel('Saccade bias (ΔHz)');
    
    bL = subplot(1,3,3);
    cfg.channel = 5;
    ft_singleplotTFR(cfg, saccadesize);
    yline(6-4/2, 'LineWidth',2, 'Color',[107, 107, 107]/255, 'LineStyle','--');
    xline(0, 'LineWidth',2, 'Color',[107, 107, 107]/255, 'LineStyle','--');
    c= colorbar();
    c.Position = [0.92 0.32 0.012 0.41];
    c.Ticks = [-0.1, 0, 0.1];
    xticks(x_ticks);
    ylim([0.25, 6]);
    ylabel('Saccade size (degrees)');
    title([]);
    
    % general
    set(gcf(), 'Position', [500 500 1450 370])

    axes = {tL, mL, bL};
    for i = 1:size(axes,2)
        xlabel(axes{i}, 'Time after cue (ms)', 'FontName', 'Aptos');
        set(axes{i}, 'Box', 'on');
        set(axes{i}, 'FontSize', [17]);
        set(axes{i}, 'FontName', 'Aptos');
        set(axes{i}, 'LineWidth', 1);
        
    end

    print("C:\Users\annav\Documents\Surfdrive\Conferences\ICON 2025\Figures\saccade_post_cue", "-dsvg")
    print("C:\Users\annav\Documents\Surfdrive\Conferences\ICON 2025\Figures\saccade_post_cue", "-dpng")

    %% effect on attentional selection latency (bias post-probe)
    line = 2;

    figure;
    top = subplot(5,1,[1:4]);
    hold on
    p1 = frevede_errorbarplot(saccade.time, squeeze(d3(:,2,:)), get_colour("blue", ""), 'se');
    p2 = frevede_errorbarplot(saccade.time, squeeze(d3(:,3,:)), get_colour("green", ""), 'se');
    p3 = frevede_errorbarplot(saccade.time, squeeze(d3(:,4,:)), get_colour("red", ""), 'se');
    p1.LineWidth = 3;
    p2.LineWidth = 3;
    p3.LineWidth = 3;
    xline(1500, 'LineWidth', 3, 'Color',[107, 107, 107]/255, 'LineStyle','--');
    yline(0, 'LineWidth',3, 'Color',[107, 107, 107]/255, 'LineStyle','--');
    ylim([-0.2 0.5]);
    yticks([-0.2, 0, 0.2, 0.4]);
    xlim([1400, 2400]);
    xticks([1500, 1800, 2100, 2400]);
    xticklabels({});
    ylabel('Saccade bias (ΔHz)');
    legend([p1, p2, p3], {'match tested', 'match none', 'match untested'}, 'EdgeColor', 'none', 'FontSize', 19.5);

    %get jackknife estimation of latency + stats
    % time = 2650:3050;  % 200 - 600 ms post probe cue onset
    time = 2450:3250;  % 0 - 800 ms post probe cue onset
    nTime = length(time);
    
    % Function to get latency of peak (or half-peak)
    get_halfpeak_latency = @(tc, t) ...
        t(find(mean(tc, 2) >= 0.5 * max(mean(tc, 2)), 1, 'first'));  % Half-peak (onset)
    
    % Initialize matrices to store jackknife latency values
    peakLatencies = zeros(size(pp2do,2), 3);
    halfLatencies = zeros(size(pp2do,2), 3);
    halfDownLatencies = zeros(size(pp2do,2), 3);

    % Single pp loop: calculate peak latencies and half-peak latencies for each participant
    for i = 1:size(pp2do,2)
        
        [~, pp_peakLatencies(i, 1)] = max(squeeze(d3(i,2,time)));
        [~, pp_peakLatencies(i, 2)] = max(squeeze(d3(i,3,time)));
        [~, pp_peakLatencies(i, 3)] = max(squeeze(d3(i,4,time)));
        
                                   % half peak cannot be before full peak!
        [~, pp_halfLatencies(i, 1)] = min(abs(squeeze(d3(i,2,time(1:pp_peakLatencies(i,1)))) - max(squeeze(d3(i,2,time))) / 2));
        [~, pp_halfLatencies(i, 2)] = min(abs(squeeze(d3(i,3,time(1:pp_peakLatencies(i,2)))) - max(squeeze(d3(i,3,time))) / 2));
        [~, pp_halfLatencies(i, 3)] = min(abs(squeeze(d3(i,4,time(1:pp_peakLatencies(i,3)))) - max(squeeze(d3(i,4,time))) / 2));

        [~, pp_halfDownLatencies(i, 1)] = min(abs(squeeze(d3(i,2,time(pp_peakLatencies(i,1):end))) - max(squeeze(d3(i,2,time))) / 2));
        [~, pp_halfDownLatencies(i, 2)] = min(abs(squeeze(d3(i,3,time(pp_peakLatencies(i,2):end))) - max(squeeze(d3(i,3,time))) / 2));
        [~, pp_halfDownLatencies(i, 3)] = min(abs(squeeze(d3(i,4,time(pp_peakLatencies(i,3):end))) - max(squeeze(d3(i,4,time))) / 2));
    end

    pp_halfDownLatencies = pp_halfDownLatencies + pp_peakLatencies;

    
    % Jackknife loop: leave-one-participant-out
    for i = 1:size(pp2do,2)
        idx = setdiff(1:size(pp2do,2), i);  % all participants except i
        
        % Compute latency estimates with one pp removed
        [~, peakLatencies(i, 1)] = max(mean(squeeze(d3(idx,2,time))));
        [~, peakLatencies(i, 2)] = max(mean(squeeze(d3(idx,3,time))));
        [~, peakLatencies(i, 3)] = max(mean(squeeze(d3(idx,4,time))));

                                   % half peak cannot be before full peak!
        [~, halfLatencies(i, 1)] = min(abs(mean(squeeze(d3(idx,2,time(1:peakLatencies(i,1))))) - max(mean(squeeze(d3(idx,2,time)))) / 2));
        [~, halfLatencies(i, 2)] = min(abs(mean(squeeze(d3(idx,3,time(1:peakLatencies(i,2))))) - max(mean(squeeze(d3(idx,3,time)))) / 2));
        [~, halfLatencies(i, 3)] = min(abs(mean(squeeze(d3(idx,4,time(1:peakLatencies(i,3))))) - max(mean(squeeze(d3(idx,4,time)))) / 2));

        [~, halfDownLatencies(i, 1)] = min(abs(mean(squeeze(d3(idx,2,time(peakLatencies(i,1):end)))) - max(mean(squeeze(d3(idx,2,time)))) / 2));
        [~, halfDownLatencies(i, 2)] = min(abs(mean(squeeze(d3(idx,3,time(peakLatencies(i,2):end)))) - max(mean(squeeze(d3(idx,3,time)))) / 2));
        [~, halfDownLatencies(i, 3)] = min(abs(mean(squeeze(d3(idx,4,time(peakLatencies(i,3):end)))) - max(mean(squeeze(d3(idx,4,time)))) / 2));
    end

    halfDownLatencies = halfDownLatencies + peakLatencies;

    % compute stats using jackknife adjusted formula
    % start by estimating individual latencies
    for i = 1:size(pp2do,2)
        halfLatency_ests(i,1) = size(pp2do, 2)*mean(halfLatencies(:,1)) - (size(pp2do, 2) - 1)*halfLatencies(i,1);
        halfLatency_ests(i,2) = size(pp2do, 2)*mean(halfLatencies(:,2)) - (size(pp2do, 2) - 1)*halfLatencies(i,2);
        halfLatency_ests(i,3) = size(pp2do, 2)*mean(halfLatencies(:,3)) - (size(pp2do, 2) - 1)*halfLatencies(i,3);
    
        peakLatency_ests(i,1) = size(pp2do, 2)*mean(peakLatencies(:,1)) - (size(pp2do, 2) - 1)*peakLatencies(i,1);
        peakLatency_ests(i,2) = size(pp2do, 2)*mean(peakLatencies(:,2)) - (size(pp2do, 2) - 1)*peakLatencies(i,2);
        peakLatency_ests(i,3) = size(pp2do, 2)*mean(peakLatencies(:,3)) - (size(pp2do, 2) - 1)*peakLatencies(i,3);

        halfDownLatency_ests(i,1) = size(pp2do, 2)*mean(halfDownLatencies(:,1)) - (size(pp2do, 2) - 1)*halfDownLatencies(i,1);
        halfDownLatency_ests(i,2) = size(pp2do, 2)*mean(halfDownLatencies(:,2)) - (size(pp2do, 2) - 1)*halfDownLatencies(i,2);
        halfDownLatency_ests(i,3) = size(pp2do, 2)*mean(halfDownLatencies(:,3)) - (size(pp2do, 2) - 1)*halfDownLatencies(i,3);
    end

    bottom = subplot(5,1,[5]);
    hold on
    xline(1500, 'LineWidth',3, 'Color',[107, 107, 107]/255, 'LineStyle','--');
    scatter(mean(halfDownLatencies)+saccade.time(time(1)), [1,2,3], 100, [get_colour("blue", "");get_colour("green", "");get_colour("red", "")], 'filled');
    errorbar(mean(halfDownLatencies(:,1))+saccade.time(time(1)),[1], std(halfDownLatency_ests(:,1)) ./ sqrt(size(pp2do,2)), "horizontal", 'LineStyle', 'none', 'Color', get_colour("blue", ""), 'LineWidth', 3);
    errorbar(mean(halfDownLatencies(:,2))+saccade.time(time(1)),[2], std(halfDownLatency_ests(:,2)) ./ sqrt(size(pp2do,2)), "horizontal", 'LineStyle', 'none', 'Color', get_colour("green", ""), 'LineWidth', 3);
    errorbar(mean(halfDownLatencies(:,3))+saccade.time(time(1)),[3], std(halfDownLatency_ests(:,3)) ./ sqrt(size(pp2do,2)), "horizontal", 'LineStyle', 'none', 'Color', get_colour("red", ""), 'LineWidth', 3);
    ylim([0,4]);
    yticks([]);
    xlim([1400, 2400]);
    xticks([1500, 1800, 2100, 2400]);
    xticklabels({'0', '300', '600', '900'});
    xlabel('Time after probe (ms)');
    % set(gca(), 'Position', [0.1300 0.0450 0.7750 0.1243]);


    % general
    set(gcf(), 'Position', [500 500 850 800])

    axes = {top, bottom};
    for i = 1:size(axes,2)
        set(axes{i}, 'Box', 'on');
        set(axes{i}, 'FontSize', [25.4]);
        set(axes{i}, 'FontName', 'Aptos');
        set(axes{i}, 'LineWidth', 1.49);
    end

    print("C:\Users\annav\Documents\Surfdrive\Conferences\ICON 2025\Figures\saccade_post_probe", "-dsvg")
    print("C:\Users\annav\Documents\Surfdrive\Conferences\ICON 2025\Figures\saccade_post_probe", "-dpng")

    
end
