%% Step3-- gaze-shift calculation

%% start clean
clear; clc; close all;

%% parameters
for pp = setdiff(1:26,[2,18,24]);

    oneOrTwoD       = 1; oneOrTwoD_options = {'_1D','_2D'};
    plotResults     = 0;

    %% load epoched data of this participant data
    param = getSubjParam(pp);
    load([param.path, '\epoched_data\eyedata_vidi5','_'  param.subjName], 'eyedata');

    %% add relevant behavioural file data

    %% only keep channels of interest
    cfg = [];
    cfg.channel = {'eyeX','eyeY'}; % only keep x & y axis
    eyedata = ft_selectdata(cfg, eyedata); % select x & y channels

    %% reformat such that all data in single matrix of trial x channel x time
    cfg = [];
    cfg.keeptrials = 'yes';
    tl = ft_timelockanalysis(cfg, eyedata); % realign the data: from trial*time cells into trial*channel*time?

    %% pixel to degree
    [dva_x, dva_y] = frevede_pixel2dva(squeeze(tl.trial(:,1,:)), squeeze(tl.trial(:,2,:)));
    tl.trial(:,1,:) = dva_x;
    tl.trial(:,2,:) = dva_y;

    %% selection vectors for conditions -- this is where it starts to become interesting!
    % cued item location
    targL = ismember(tl.trialinfo(:,1), [21,23,25,27,29,211,213,215,217,219,221,223]);
    targR = ismember(tl.trialinfo(:,1), [22,24,26,28,210,212,214,216,218,220,222,224]);
    
    congruent = ismember(tl.trialinfo(:,1), [21,22,25,26,213,214,217,218]);
    incongruent = ismember(tl.trialinfo(:,1), [23,24,27,28,215,216,219,220]);
    neutral = ismember(tl.trialinfo(:,1), [29,210,211,212,221,222,223,224]);
    neutral_3 = ismember(tl.trialinfo(:,1), [29,210,221,222]);
    neutral_4 = ismember(tl.trialinfo(:,1), [211,212,223,224]);

    captureL = ismember(tl.trialinfo(:,1), [21,24,25,28,213,216,217,220]);
    captureR = ismember(tl.trialinfo(:,1), [22,23,26,27,214,215,218,219]);

    % responded or not
    behdata = readtable(param.log);
    if size(behdata,1) ~= size(eyedata.trial, 2) %quick sanity check just in case
        continue
        throw(MException('get_saccadebias:trial_numbers', "Loaded behavioural data has a different number of trials to the eyedata."));
    end

    pressed = logical(ismember(behdata.cue_hit, {'True'}) + ismember(behdata.cue_false_alarm, {'True'}));
    not_pressed = and(ismember(behdata.cue_hit, {'False'}), ismember(behdata.cue_false_alarm, {'False'}));

    if sum(pressed) + sum(not_pressed) ~= size(eyedata.trial, 2) %quick sanity check just in case
        throw(MException('get_saccadebias:pressed_trials', "Some trials were not recognised as 'pressed' or 'not_pressed'."));
    end

    % should respond or not
    should_respond = ismember(tl.trialinfo(:,1), [29,210,213:220,223,224]);
    should_not_respond = ismember(tl.trialinfo(:,1), [21:28,211,212,221,222]);

    % block type
    respond_neutral = ismember(tl.trialinfo(:,1), [21:29, 210:212]);
    respond_not_neutral = ismember(tl.trialinfo(:,1), [213:224]);

    % median-split on behaviour
    behdata = readtable(getSubjParam(pp).log);
    oktrials = abs(zscore(behdata.idle_reaction_time_in_ms))<=3; 
    bm_dt = (behdata.idle_reaction_time_in_ms < median(behdata.idle_reaction_time_in_ms(oktrials)))&oktrials;
    am_dt = (behdata.idle_reaction_time_in_ms > median(behdata.idle_reaction_time_in_ms(oktrials)))&oktrials;
    bm_er = (behdata.absolute_difference < median(behdata.absolute_difference(oktrials)))&oktrials;
    am_er = (behdata.absolute_difference > median(behdata.absolute_difference(oktrials)))&oktrials;
    
    % channels
    chX = ismember(tl.label, 'eyeX');
    chY = ismember(tl.label, 'eyeY');

    %% get gaze shifts using our custom function
    cfg = [];
    data_input = squeeze(tl.trial);
    time_input = tl.time*1000;

    if oneOrTwoD == 1         [shiftsX, velocity, times]             = PBlab_gazepos2shift_1D(cfg, data_input(:,chX,:), time_input);
    elseif oneOrTwoD == 2     [shiftsX,shiftsY, peakvelocity, times] = PBlab_gazepos2shift_2D(cfg, data_input(:,chX,:), data_input(:,chY,:), time_input);
    end

    %% select usable gaze shifts
    minDisplacement = 0;
    maxDisplacement = 1000;

    if oneOrTwoD == 1     saccadesize = abs(shiftsX);
    elseif oneOrTwoD == 2 saccadesize = abs(shiftsX+shiftsY*1i);
    end
    shiftsL = shiftsX<0 & (saccadesize>minDisplacement & saccadesize<maxDisplacement);
    shiftsR = shiftsX>0 & (saccadesize>minDisplacement & saccadesize<maxDisplacement);

    %% get relevant contrasts out
    saccade = [];
    saccade.time = times;
    saccade.label = {'all', 'con_tar', 'neu_tar', 'incon_tar', ...
        'neu_special_tar', 'neu_no_task_tar', 'cuematch_cue', ...
        'con_cue', 'incon_cue', 'cuematch_press', 'cuematch_notpress', ...
        'cuematch_should_press', 'cuematch_should_notpress',...
        'c_bm_dt', 'i_bm_dt', 'c_bm_er', 'i_bm_er', ...
        'c_am_dt', 'i_am_dt', 'c_am_er', 'i_am_er'};

    for selection = [1:21] % conditions.
        if     selection == 1  sel = ones(size(targL));
        elseif selection == 2  sel = congruent;
        elseif selection == 3  sel = neutral;
        elseif selection == 4  sel = incongruent;
        elseif selection == 5  sel = neutral_3;
        elseif selection == 6  sel = neutral_4;
        elseif selection == 7  sel = any([congruent;incongruent]);
        elseif selection == 8  sel = congruent;
        elseif selection == 9  sel = incongruent;
        elseif selection == 10  sel = any([congruent;incongruent])&pressed;
        elseif selection == 11 sel = any([congruent;incongruent])&not_pressed;
        elseif selection == 12 sel = any([congruent;incongruent])&should_respond;
        elseif selection == 13 sel = any([congruent;incongruent])&should_not_respond;
        elseif selection == 14 sel = congruent&bm_dt;
        elseif selection == 15 sel = incongruent&bm_dt;
        elseif selection == 16 sel = congruent&bm_er;
        elseif selection == 17 sel = incongruent&bm_er;
        elseif selection == 18 sel = congruent&am_dt;
        elseif selection == 19 sel = incongruent&am_dt;
        elseif selection == 20 sel = congruent&am_er;
        elseif selection == 21 sel = incongruent&am_er;
        end

        if selection <= 6
            saccade.toward(selection,:) =  (mean(shiftsL(targL&sel,:)) + mean(shiftsR(targR&sel,:))) ./ 2;
            saccade.away(selection,:)  =   (mean(shiftsL(targR&sel,:)) + mean(shiftsR(targL&sel,:))) ./ 2;

        else
            saccade.toward(selection,:) =  (mean(shiftsL(captureL&sel,:)) + mean(shiftsR(captureR&sel,:))) ./ 2;
            saccade.away(selection,:)  =   (mean(shiftsL(captureR&sel,:)) + mean(shiftsR(captureL&sel,:))) ./ 2;
        end
    
    % add towardness field
    saccade.effect = (saccade.toward - saccade.away);
    
    end
    
    %% smooth and turn to Hz
    integrationwindow = 100; % window over which to integrate saccade counts
    saccade.toward = smoothdata(saccade.toward,2,'movmean',integrationwindow)*1000; % *1000 to get to Hz, given 1000 samples per second.
    saccade.away   = smoothdata(saccade.away,2,  'movmean',integrationwindow)*1000;
    saccade.effect = smoothdata(saccade.effect,2,'movmean',integrationwindow)*1000;

    %% plot
    if plotResults
        figure;
        for sp = 1:4
            subplot(2,2,sp);
            hold on;
            plot(saccade.time, saccade.toward(sp,:), 'r');
            plot(saccade.time, saccade.away(sp,:), 'b');
            title(saccade.label(sp));
            legend({'toward','away'},'autoupdate', 'off');
            plot([0,0], ylim, '--k');
            plot([1500,1500], ylim, '--k');
        end

        figure;
        for sp = 1:4
            subplot(2,2,sp); hold on;
            plot(saccade.time, saccade.effect(sp,:), 'k');
            plot(xlim, [0,0], '--k');
            title(saccade.label(sp));
            legend({'effect'},'autoupdate', 'off');
            plot([0,0], ylim, '--k');
            plot([1500,1500], ylim, '--k');
        end

        figure;
        hold on;
        plot(saccade.time, saccade.effect([1:4],:));
        plot(xlim, [0,0], '--k');
        legend(saccade.label([1:3]),'autoupdate', 'off');
        plot([0,0], ylim, '--k');
        plot([1500,1500], ylim, '--k');
        drawnow;
    end

    %% also get as function of saccade size - identical as above, except with extra loop over saccade size.
    binsize = 0.5;
    halfbin = binsize/2;

    saccadesize = [];
    saccadesize.dimord = 'chan_freq_time';
    saccadesize.freq = halfbin:0.1:7-halfbin; % shift sizes, as if "frequency axis" for time-frequency plot
    saccadesize.time = times;
    saccadesize.label = {'all', 'con_tar', 'neu_tar', 'incon_tar', 'neu_special_tar', 'neu_no_task_tar', 'cuematch_cue'};
    
    cnt = 0;
    for sz = saccadesize.freq;
        cnt = cnt+1;
        shiftsL = [];
        shiftsR = [];
        shiftsL = shiftsX<-sz+halfbin & shiftsX > -sz-halfbin; % left shifts within this range
        shiftsR = shiftsX>sz-halfbin  & shiftsX < sz+halfbin; % right shifts within this range

        for selection = [1:7] % conditions.
            if     selection == 1  sel = ones(size(targL));
            elseif selection == 2  sel = congruent;
            elseif selection == 3  sel = neutral;
            elseif selection == 4  sel = incongruent;
            elseif selection == 5  sel = neutral_3;
            elseif selection == 6  sel = neutral_4;
            elseif selection == 7  sel = any([congruent;incongruent]);
            end

            if selection <= 6
                saccadesize.toward(selection,cnt,:) = (mean(shiftsL(targL&sel,:)) + mean(shiftsR(targR&sel,:))) ./ 2;
                saccadesize.away(selection,cnt,:) =   (mean(shiftsL(targR&sel,:)) + mean(shiftsR(targL&sel,:))) ./ 2;
            else
                saccadesize.toward(selection,cnt,:) =  (mean(shiftsL(captureL&sel,:)) + mean(shiftsR(captureR&sel,:))) ./ 2;
                saccadesize.away(selection,cnt,:)  =   (mean(shiftsL(captureR&sel,:)) + mean(shiftsR(captureL&sel,:))) ./ 2;
            end

        end

    end
    % add towardness field
    saccadesize.effect = (saccadesize.toward - saccadesize.away);

    %% smooth and turn to Hz
    integrationwindow = 100; % window over which to integrate saccade counts
    saccadesize.toward = smoothdata(saccadesize.toward,3,'movmean',integrationwindow)*1000; % *1000 to get to Hz, given 1000 samples per second.
    saccadesize.away   = smoothdata(saccadesize.away,3,  'movmean',integrationwindow)*1000;
    saccadesize.effect = smoothdata(saccadesize.effect,3,'movmean',integrationwindow)*1000;

    if plotResults
        cfg = [];
        cfg.parameter = 'effect';
        cfg.figure = 'gcf';
        %cfg.zlim = [-0.01, 0.01];
        figure;
        for chan = 1:5
            cfg.channel = chan;
            subplot(2,3,chan); ft_singleplotTFR(cfg, saccadesize);
        end
        colormap('jet');
        drawnow;
    end

    %% save
    save([param.path, '\saved_data\saccadeEffects', oneOrTwoD_options{oneOrTwoD} '__', param.subjName], 'saccade','saccadesize');

    %% close loops
end % end pp loop
