function param = getSubjParam(pp)

%% participant-specific notes

%% set path and pp-specific file locations
unique_numbers = [25,44,43,24,59,13,14,63,37,21,62,64,42,71,35,89,61,97,26,78,91,47,70,94,88,46]; %needs to be in the right order

param.path = '\\labsdfs.labs.vu.nl\labsdfs\FGB-ETP-CogPsy-ProactiveBrainLab\core_lab_members\Anna\Data\vidi5 - action-coupled-null-cue\';

if pp < 10
    param.subjName = sprintf('pp0%d', pp);
else
    param.subjName = sprintf('pp%d', pp);
end

log_string = sprintf('data_session_%d.csv', pp);
param.log = [param.path, log_string];

eds_string = sprintf('%d_%d.asc', pp, unique_numbers(pp));
param.eds = [param.path, eds_string];