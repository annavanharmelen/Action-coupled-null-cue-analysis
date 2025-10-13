function param = getAllSubjParam(study, pp)

%% set path and pp-specific file locations, dependent on study
if study == 1
    unique_numbers = [25,44,43,24,59,13,14,63,37,21,62,64,42,71,35,89,61,97,26,78,91,47,70,94,88,46]; %needs to be in the right order
    param.path = '\\labsdfs.labs.vu.nl\labsdfs\FGB-ETP-CogPsy-ProactiveBrainLab\core_lab_members\Anna\Data\vidi5 - action-coupled-null-cue\';

elseif study == 2
    unique_numbers = [89, 57, 75, 12, 92, 95, 90, 79, 47, 10, 85, 53, 99, 34, 74, 55, 94, 68, 73, 97, 50, 36, 93, 35, 61, 84, 56]; %needs to be in the right order
    param.path = '\\labsdfs.labs.vu.nl\labsdfs\FGB-ETP-CogPsy-ProactiveBrainLab\core_lab_members\Anna\Data\vidi3 - location-by-colour null-cue\';

elseif study == 3
    unique_numbers = [95, 37, 92, 74, 49, 52, 99, 39, 96, 42, 90, 13, 72, 44, 85, 24, 76, 71, 50, 47, 32, 59, 19, 93, 88]; %needs to be in the right order
    param.path = '\\labsdfs.labs.vu.nl\labsdfs\FGB-ETP-CogPsy-ProactiveBrainLab\core_lab_members\Anna\Data\vidi3.2 - location-by-colour null-cue\';

end

if pp < 10
    param.subjName = sprintf('pp0%d', pp);
else
    param.subjName = sprintf('pp%d', pp);
end

log_string = sprintf('data_session_%d.csv', pp);
param.log = [param.path, log_string];

eds_string = sprintf('%d_%d.asc', pp, unique_numbers(pp));
param.eds = [param.path, eds_string];