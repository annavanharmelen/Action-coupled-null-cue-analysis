# ACTION-COUPLED NULL-CUE ANALYSIS 

# ***[work in progress]***

Analysis scripts (in MATLAB) for the data acquired from the action-coupled null-cue experiment(in Python). For the experiment, see: [Action-coupled-main](https://github.com/annavanharmelen/Action-coupled-null-cue-experiment).

## Author
Made by Anna van Harmelen in 2025, with scripts from Dr. Freek van Ede.

## Installation
Some of these analysis scripts are dependent on the [Fieldtrip toolbox](https://www.fieldtriptoolbox.org), and were originally built using the 2020.10.23 version of Fieldtrip.

## Configuration
To make sure the scripts run correctly, open the getSubjParam.m file to either...:
- Enter the randomised participant numbers (in order of session number), if your filing system is the same as mine.
- Change the code, so this function can find the data corresponding to each participant.

## Running
The analysis runs in multiple parts, some of which are dependent on each other.

All main behavioural data is analysed in get_behaviour.m. Some extra behavioural data is analysed in extra_behavioural_data.m. If you have additional questionnaire data it is loaded in and analysed in qualtrics_data.m.

All eye-tracking data is analysed by entering the desired participantnumbers into these scripts, and then running them in precisely this sequence:
1. epoch_eye_data.m
2. get_gazePositionBias.m and get_saccadeBias.m
3. GA_gazePositionBias.m and GA_saccadeBias.m

Note that most of these eye-tracking analysis scripts are dependent on the 'frevede_' functions in the repository.
