%%

%% BROAD-NESS - (AUDITORY PATTERN RECOGNITION 2020) - LEONARDO BONETTI

%%

%% *** START UP FUNCTIONS.. (LBPD_startup_D) ***

%starting up some functions for LBPD toolbox.

%starting up some of the functions that I wrote for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl);

%%

%% *** PREPROCESSING ***

%%% OBS!! the preprocessing was computed for (nearly) all the data
%%% collected and not only for the data that was actually analyzed and
%%% reported in this paper.

%%

%% Maxfilter

%OBS! before running maxfilter you need to close matlab, open the terminal and write: 'use anaconda', then open matlab and run maxfilter script

maxfilter_path = '/neuro/bin/util/maxfilter';
project = 'MINDLAB2020_MEG-AuditoryPatternRecognition';
maxDir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2'; %output path

path = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition'; %path with all the subjects folders
jj = dir([path '/0*']); %list all the folders starting with '0' in order to avoid hidden files
for ii = 90:length(jj) %over subjects (STARTING FROM 13, SINCE THE PREVIOUS ONES WERE THE PILOTS (IF THEY DID ALSO THE PROPER EXPERIMENT, THE MAXFILTER COMPUTATION WILL BE DONE IN THE NEXT  SECTION))
    cart = [jj(ii).folder '/' jj(ii).name]; %create a path that combines the folder name and the file name
    pnana = dir([cart '/2*']); %search for folders starting with '2'
    for pp = 1:length(pnana) %loop to explore ad analyze all the folders inside the path above
        cart2 = [pnana(1).folder '/' pnana(pp).name];
        pr = dir([cart2 '/ME*']); %looks for meg folder
        if ~isempty(pr) %if pr is not empty, proceed with subfolders inside the meg path
            pnunu = dir([pr(1).folder '/' pr(1).name '/00*']);
            %if length(pnunu) > 1
            %warning(['subj ' num2str(ii) ' has files nuber = ' num2str(length(pnunu))]) %show a warning message if any subj has more thatn 1 meg sub-folder
            %end
            for dd = 1:length(pnunu)
                if strcmp(pnunu(dd).name(5:6),'re') || strcmp(pnunu(dd).name(5:6),'sa') || strcmp(pnunu(dd).name(5:6),'vi') || strcmp(pnunu(dd).name(5:6),'pd') %checks whether characters 5 to 6 are equal to 're', 'vi, 'sa' or 'pd'; the loop continues if this is true (1) and it stops if this is false (0)
                    %idx2 = strfind(pnunu(1).name,'Mus'); % search for musmelo folder in order to avoid other projects
                    %if ~isempty(idx2)
                    fpath = dir([pnunu(1).folder '/' pnunu(dd).name '/files/*.fif']); % looks for .fif file
                    rawName = ([fpath.folder '/' fpath.name]); %assigns the final path of the .fif file to the rawName path used in the maxfilter command
                    maxfName = ['SUBJ' jj(ii).name '_' fpath.name(1:end-4)]; %define the output name of the maxfilter processing
                    %movement compensation
                    cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -movecomp -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
                    %no movement compensation (to be used if HPI coils did not work properly)
%                     cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
                    system(cmd);
                end
            end
        end
    end
end

%% Maxfilter - Subjects who did the experiment (version 2.0) after the first pilot (version 1.0)

%%% DONE SUBJECTS 2 10 11 7 %%%

%OBS! before running maxfilter you need to close matlab, open the terminal and write: 'use anaconda', then open matlab and run maxfilter script

%settings
maxfilter_path = '/neuro/bin/util/maxfilter';
project = 'MINDLAB2020_MEG-AuditoryPatternRecognition';
maxDir = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2'; %output path

clear rawName maxfName
%path to raw files for these particular subjects
pathraw{1} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0002/20210413_000000/MEG'; %assigns the final path of the .fif file to the rawName path used in the maxfilter command
pathraw{2} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0010/20210416_000000/MEG';
pathraw{3} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0011/20210423_000000/MEG';
pathraw{4} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0007/20210511_000000/MEG';
pathraw{5} = '/raw/sorted/MINDLAB2020_MEG-AuditoryPatternRecognition/0012/20210622_000000/MEG';

for ii = 5:length(pathraw) %over subjects
    list = dir([pathraw{ii} '/0*']);
    for jj = 1:length(list) %over experimental blocks
        if strcmp(list(jj).name(5:6),'re') || strcmp(list(jj).name(5:6),'sa') || strcmp(list(jj).name(5:6),'vi') || strcmp(list(jj).name(5:6),'pd') %checks whether characters 5 to 6 are equal to 're', 'vi, 'sa' or 'pd'; the loop continues if this is true (1) and it stops if this is false (0)
            
            fpath = dir([list(jj).folder '/' list(jj).name '/files/*.fif']); % looks for .fif file
            rawName = [fpath(1).folder '/' fpath(1).name];
            maxfName = ['SUBJ' pathraw{ii}(56:59) '_' fpath.name(1:end-4) '_bis']; %define the output name of the maxfilter processing
            %movement compensation
            cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -movecomp -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
            %no movement compensation
%             cmd = ['submit_to_cluster -q maxfilter.q -n 4 -p ' ,project, ' "',maxfilter_path,' -f ',[rawName],' -o ' [maxDir '/' maxfName '_tsssdsm.fif'] ' -st 4 -corr 0.98 -ds 4 ',' -format float -v | tee ' [maxDir '/log_files/' maxfName '_tsssdsm.log"']];
            system(cmd);
        end
    end
end

%% Converting the .fif files into SPM objects

%OBS! remember to run 'starting up OSL' first

%setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu

%%

fif_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/*.fif'); %creates a list with the .fif files

for ii = 340:341%:length(fif_list) %over the .fif files
    S = []; %structure 'S'                   
    S.dataset = [fif_list(ii).folder '/' fif_list(ii).name];
    D = spm_eeg_convert(S);
%     D = job2cluster(@cluster_spmobject, S); %actual function for conversion
end

%% Removing bad segments using OSLVIEW

%checks data for potential bad segments (periods)
%marking is done by right-clicking in the proximity of the event and click on 'mark event'
%a first click (green dashed label) marks the beginning of a bad period
%a second click indicates the end of a bad period (red)
%this will mean that we are not using about half of the data, but with such bad artefacts this is the best we can do
%we can still obtain good results with what remains
%NB: Push the disk button to save to disk (no prefix will be added, same name is kept)

%OBS! remember to check for bad segments of the signal both at 'megplanar' and 'megmag' channels (you can change the channels in the OSLVIEW interface)
%OBS! remember to mark the trial within the bad segments as 'badtrials' and use the label for removing them from the Averaging (after Epoching) 

spm_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/spmeeg*.mat'); %path to SPM objects

for ii = 1:3%:length(spm_list) %over experimental blocks %OBS!
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    D = oslview(D);
    D.save(); %save the selected bad segments and/or channels in OSLVIEW
    disp(ii)
end

%% AFRICA denoising (part I)

%setting up the cluster
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu

%%

%ICA calculation
spm_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/spmeeg*.mat');

for ii = 339:length(spm_list) %OBS!
    S = [];
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    S.D = D;
    
    jobid = job2cluster(@cluster_africa,S);
%   D = osl_africa(D,'do_ica',true,'do_ident',false,'do_remove',false,'used_maxfilter',true); 
%   D.save();
end

%% AFRICA denoising (part II)

% v = [11 12 19 32];
%visual inspection and removal of artifacted components
%look for EOG and ECG channels (usually the most correlated ones, but check a few more just in case)
spm_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/spmeeg*.mat');

for ii = 339:length(spm_list) %OBS!%38:41
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]);
    D = osl_africa(D,'do_ident','manual','do_remove',false,'artefact_channels',{'EOG','ECG'});
    %hacking the function to manage to get around the OUT OF MEMORY problem..
    S = [];
    S.D = D;
    jobid = job2cluster(@cluster_rembadcomp,S);
%   D.save();
    disp(ii)
end

%% Epoching: one epoch per old/new excerpt (baseline = (-100ms)

prefix_tobeadded = 'e'; %adds this prefix to epoched files
spm_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/spmeeg*.mat');

for ii = 339:length(spm_list) %over .mat files
    D = spm_eeg_load([spm_list(ii).folder '/' spm_list(ii).name]); %load spm_list .mat files
    dummy = D.fname; %OBS! D.fname does not work, so we need to use a 'dummy' variable instead
    %if strcmp(dummy(22:26), 'speed') %checks whether characters 22 to 26 are equal to 'speed'; the loop continues if this is true (1) and it stops if this is false (0)
    events = D.events; %look for triggers
    %takes the correct triggers sent during the recording
    clear trigcor
    count_evval = 0; %???
    for ieve = 1:length(events) %over triggers
        if strcmp(events(ieve).type,'STI101_up') %only triggers at the beginning of each stimuli
            if events(ieve).value ~= 103 && events(ieve).value ~= 104 && events(ieve).value ~= 128 && events(ieve).value ~= 8 && events(ieve).value ~= 132 && events(ieve).value ~= 48 && events(ieve).value ~= 32 && events(ieve).value ~= 64 %discard 104 and 128 for random triggers
                count_evval = count_evval + 1;
                trigcor(count_evval,1) = events(ieve).time; %+ 0.010; %this takes the correct triggers and add 10ms of delay of the sound travelling into the tubes
                %variable with all the triggers we need
            end
        end
    end
    trl_sam = zeros(length(trigcor),3); %prepare the samples matrix with 0's in all its cells
    trl_sec = zeros(length(trigcor),3); %prepare the seconds matrix with 0's in all its cells
    %deftrig = zeros(length(trigcor),1); %this is not useful
    for k = 1:length(trigcor) %over selected triggers
        %deftrig(k,1) = 0.012 + trigcor(k,1); %adding a 0.012 seconds delay to the triggers sent during the experiment (this delay was due to technical reasons related to the stimuli)
        trl_sec(k,1) = trigcor(k,1) - 0.1000; %beginning time-window epoch in s (please note that we computed this operation two times, obtaining two slightly different pre-stimulus times.
        %this was done because for some computations was convenient to have a slightly longer pre-stimulus time
        %remove 1000ms of baseline
        if strcmp(dummy(22:26), 'speed')
            trl_sec(k,2) = trigcor(k,1) + 5.5; %end time-window epoch in seconds
        else
            trl_sec(k,2) = trigcor(k,1) + 4.4; %end time-window epoch in seconds
        end
        trl_sec(k,3) = trl_sec(k,2) - trl_sec(k,1); %range time-windows in seconds
        trl_sam(k,1) = round(trl_sec(k,1) * 250) + 1; %beginning time-window epoch in samples %250Hz per second
        trl_sam(k,2) = round(trl_sec(k,2) * 250) + 1; %end time-window epoch in samples
        trl_sam(k,3) = -25; %sample before the onset of the stimulus (corresponds to 0.100ms)
    end
    dif = trl_sam(:,2) - trl_sam(:, 1); %difference between the end and the beginning of each sample (just to make sure that everything is fine)
    if ~all(dif == dif(1)) %checking if every element of the vector are the same (i.e. the length of the trials is the same; we may have 1 sample of difference sometimes because of different rounding operations..)
        trl_sam(:,2) = trl_sam(:,1) + dif(1);
    end
    %creates the epochinfo structure that is required for the source reconstruction later
    epochinfo.trl = trl_sam;
    epochinfo.time_continuous = D.time;
    %switch the montage to 0 because for some reason OSL people prefer to do the epoching with the not denoised data
    D = D.montage('switch',0);
    %build structure for spm_eeg_epochs
    S = [];
    S.D = D;
    S.trl = trl_sam;
    S.prefix = prefix_tobeadded;
    D = spm_eeg_epochs(S);
    
    %store the epochinfo structure inside the D object
    D.epochinfo = epochinfo;
    D.save();
    %take bad segments registered in OSLVIEW and check if they overlap with the trials. if so, it gives the number of overlapped trials that will be removed later   
    count = 0;
    Bad_trials = zeros(length(trigcor),1);
    for kkk = 1:length(events) %over events
        if strcmp(events(kkk).type,'artefact_OSL')
            for k = 1:length(trl_sec) %over trials
                if events(kkk).time - trl_sec(k,2) < 0 %if end of trial is > than beginning of artifact
                    if trl_sec(k,1) < (events(kkk).time + events(kkk).duration) %if beginning of trial is < than end of artifact
                        Bad_trials(k,1) = 1; %it is a bad trial (stored here)
                        count = count + 1;
                    end
                end                  
            end
        end
    end
    %if bad trials were detected, their indices are stored within D.badtrials field
    disp(spm_list(ii).name);
    if count == 0
        disp('there are no bad trials marked in oslview');
    else
        D = badtrials(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
%         D = conditions(D,find(Bad_trials),1); %get the indices of the badtrials marked as '1' (that means bad)
        epochinfo = D.epochinfo;
        xcv = find(Bad_trials == 1);
        %this should be done only later.. in any case.. not a problem..
        for jhk = 1:length(xcv)
            D = D.conditions(xcv(jhk),'Bad');
            epochinfo.conditionlabels(xcv(jhk)) = {'Bad'};
            disp([num2str(ii) ' - ' num2str(jhk) ' / ' num2str(length(xcv))])
        end
        D.epochinfo = epochinfo;
        D.save(); %saving on disk
        disp('bad trials are ')
        length(D.badtrials)
    end
    D.save();
    disp(ii)
end


%% Defining the conditions - All blocks

%define conditions - 1 epoch for each old/new excerpt (baseline = (-)100ms)

xlsx_dir_behav = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/BehavioralTaskMEG/Version_2/Final_xlsx'; %dir to MEG behavioral results (.xlsx files)
epoch_list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*.mat'); %dir to epoched files

for ii = 339:length(epoch_list) %over epoched data
    D = spm_eeg_load([epoch_list(ii).folder '/' epoch_list(ii).name]);
    dummy = D.fname;
    %barbaric solution.. to build the name to be read for the excel files with the MEG behavioral tasks performance
    if strcmp(dummy(18:23),'recogm')
        dumbloc = 'Block_3.xlsx';
        bl = 3;
    elseif strcmp(dummy(18:23),'recogs')
        dumbloc = 'Block_4.xlsx';
        bl = 4;
    elseif strcmp(dummy(18:23),'sameme')
        dumbloc = 'Block_5.xlsx';
        bl = 5;
    elseif strcmp(dummy(18:23),'visual')
        dumbloc = 'Block_6.xlsx';
        bl = 6;
    elseif strcmp(dummy(18:19),'pd')
        dumbloc = 'Project_PD.xlsx';
        bl = 7;
    end
    if strcmp(dummy((end-14):(end-14)+2),'bis')
        dumls = ['Subj_' dummy(13:16) 'bis_' dumbloc]; %getting subject ID directly from the SPM object to reduce the probability to make mistakes
    else
        dumls = ['Subj_' dummy(13:16) '_' dumbloc];
    end
    [~,~,raw_recog] = xlsread([xlsx_dir_behav '/' dumls]); %excel files
    %picking the current block
    if bl == 3 %block 3 (THIS IS THE ACTUAL BLOCK FOR THIS PAPER)
        for k = 1:length(D.trialonset)
            if raw_recog{(k + 1),3} == 0 %if there was no response
                D = D.conditions(k,'No_response');
            elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 1 %old correct
                D = D.conditions(k,'Old_Correct'); %assign old correct
            elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 2 %old incorrect
                D = D.conditions(k,'Old_Incorrect'); %otherwise assign new correct
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 2 %new t1 correct
                D = D.conditions(k,'New_T1_Correct');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 1 %new t1 incorrect
                D = D.conditions(k,'New_T1_Incorrect');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t2') && raw_recog{(k + 1),3} == 2 %new t2 correct
                D = D.conditions(k,'New_T2_Correct');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t2') && raw_recog{(k + 1),3} == 1 %new t2 incorrect
                D = D.conditions(k,'New_T2_Incorrect');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 2 %new t3 correct
                D = D.conditions(k,'New_T3_Correct');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 1 %new t3 incorrect
                D = D.conditions(k,'New_T3_Incorrect');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t4') && raw_recog{(k + 1),3} == 2 %new t4 correct
                D = D.conditions(k,'New_T4_Correct');
            elseif strcmp(raw_recog{(k + 1),2}(14:15),'t4') && raw_recog{(k + 1),3} == 1 %new t4 incorrect
                D = D.conditions(k,'New_T4_Incorrect');
            end
        end
    elseif bl == 4 %block 4
        for k = 1:length(D.trialonset)
            if raw_recog{(k + 1),3} == 0 %if there was no response
                D = D.conditions(k,'No_response');
            elseif strcmp(raw_recog{(k + 1),2}(6:13),'fast_old') && raw_recog{(k + 1),3} == 1 %old fast correct
                D = D.conditions(k,'Old_Fast_Correct'); %assign old fast correct
            elseif strcmp(raw_recog{(k + 1),2}(6:13),'fast_old') && raw_recog{(k + 1),3} == 2 %old fast incorrect
                D = D.conditions(k,'Old_Fast_Incorrect'); %assign old fast incorrect
            elseif strcmp(raw_recog{(k + 1),2}(6:13),'slow_old') && raw_recog{(k + 1),3} == 1 %old slow correct
                D = D.conditions(k,'Old_Slow_Correct'); %assign old slow correct
            elseif strcmp(raw_recog{(k + 1),2}(6:13),'slow_old') && raw_recog{(k + 1),3} == 2 %old slow incorrect
                D = D.conditions(k,'Old_Slow_Incorrect'); %assign old slow incorrect
            elseif strcmp(raw_recog{(k + 1),2}(8:15),'fast_new') && raw_recog{(k + 1),3} == 2 %new fast correct
                D = D.conditions(k,'New_Fast_Correct'); %assign new fast correct
            elseif strcmp(raw_recog{(k + 1),2}(8:15),'fast_new') && raw_recog{(k + 1),3} == 1 %new fast incorrect
                D = D.conditions(k,'New_Fast_Incorrect'); %assign new fast incorrect
            elseif strcmp(raw_recog{(k + 1),2}(8:15),'slow_new') && raw_recog{(k + 1),3} == 2 %new slow correct
                D = D.conditions(k,'New_Slow_Correct'); %assign new slow correct
            else
                D = D.conditions(k,'New_Slow_Incorrect'); %assign new incorrect
            end
        end
    elseif bl == 5 %block 5
        for k = 1:20 %over the 20 encoding trials
            D = D.conditions(k,'Encoding');
        end
        for k = 1:(length(D.trialonset)-20) %over recognition trials (this would work even if the last trial was out of bonds and thus was not epoched..)
            if raw_recog{(k + 23),3} == 0 %if there was no response
                D = D.conditions(k + 20,'No_response');
            elseif strcmp(raw_recog{(k + 23),2}(9:11),'enc') && raw_recog{(k + 23),3} == 1 %encoding correct
                D = D.conditions(k + 20,'Old_Correct'); %assign encoding correct
            elseif strcmp(raw_recog{(k + 23),2}(9:11),'enc') && raw_recog{(k + 23),3} == 2 %encoding incorrect
                D = D.conditions(k + 20,'Old_Incorrect'); %otherwise assign encoding incorrect
            elseif strcmp(raw_recog{(k + 23),2}(9:11),'rec') && raw_recog{(k + 23),3} == 2 %recognition correct
                D = D.conditions(k + 20,'New_Correct'); %assign recognition correct
            else
                D = D.conditions(k + 20,'New_Incorrect'); %assign new incorrect
            end
        end
    elseif bl == 6 %block 6
        for k = 1:20 %over the 20 encoding trials
            D = D.conditions(k,'Encoding');
        end
        for k = 1:(length(D.trialonset)-20) %over recognition trials (this would work even if the last trial was out of bonds and thus was not epoched..)
            if raw_recog{(k + 23),3} == 0 %if there was no response
                D = D.conditions(k + 20,'No_response');
            elseif strcmp(raw_recog{(k + 23),2}(8:10),'old') && raw_recog{(k + 23),3} == 1 %encoding correct
                D = D.conditions(k + 20,'Old_Correct'); %assign encoding correct
            elseif strcmp(raw_recog{(k + 23),2}(8:10),'old') && raw_recog{(k + 23),3} == 2 %encoding incorrect
                D = D.conditions(k + 20,'Old_Incorrect'); %otherwise assign encoding incorrect
            elseif strcmp(raw_recog{(k + 23),2}(8:10),'new') && raw_recog{(k + 23),3} == 2 %recognition correct
                D = D.conditions(k + 20,'New_Correct'); %assign recognition correct
            else
                D = D.conditions(k + 20,'New_Incorrect'); %assign new incorrect
            end
        end
    elseif bl == 7
        for k = 1:length(D.trialonset) %over trials
            if strcmp(raw_recog{(k + 1),2}(3:7),'porco')
                D = D.conditions(k,'pd');
            else
                D = D.conditions(k,'pm');
            end
        end
    end
    %this is for every block
    if ~isempty(D.badtrials) %overwriting badtrials (if any) on condition labels
        BadTrials = D.badtrials;
        for badcount = 1:length(BadTrials) %over bad trials
            D = D.conditions(BadTrials(badcount),'Bad_trial');
        end
    end
    D = D.montage('switch',1);
    D.epochinfo.conditionlabels = D.conditions; %to add for later use in the source reconstruction
    D.save(); %saving data on disk
    disp(num2str(ii))
end

%%

%% *** SOURCE RECONSTRUCTION (LBPD) ***

%%

%% CREATING 8mm PARCELLATION FOR EASIER INSPECTION IN FSLEYES
%OBS!! This section is done only for better handling of some visualization purposes, but it does not affect any of the beamforming algorithm;
% it is just important not to mix up the MNI coordinates, thus I would recommend to use the following lines

%1) USE load_nii TO LOAD A PREVIOUS NIFTI IMAGE
imag_8mm = load_nii('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_T1_8mm_brain.nii.gz');
Minfo = size(imag_8mm.img); %get info about the size of the original image
M8 = zeros(Minfo(1), Minfo(2), Minfo(3)); %Initialize an empty matrix with the same dimensions as the original .nii image
cc = 0; %set a counter
M1 = imag_8mm.img;
for ii = 1:Minfo(1) %loop across each voxel of every dimension
    for jj = 1:Minfo(2)
        for zz = 1:Minfo(3)
            if M1(ii,jj,zz) ~= 0 %if we have an actual brain voxel
                cc = cc+1;
                M8(ii,jj,zz) = cc;
            end
        end
    end
end
%2) PUT YOUR MATRIX IN THE FIELD ".img"
imag_8mm.img = M8; %assign values to new matrix 
%3) SAVE NIFTI IMAGE USING save_nii
save_nii(imag_8mm,'/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_brain_diy.nii.gz');
%4) USE FSLEYES TO LOOK AT THE FIGURE
%Create parcellation on the 8mm template
for ii = 1:3559 %for each 8mm voxel
    cmd = ['fslmaths /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_brain_diy.nii.nii.gz -thr ' num2str(ii) ' -uthr ' num2str(ii) ' -bin /scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/AAL_80mm_3559ROIs/' num2str(ii) '.nii.gz'];
    system(cmd)
    disp(ii)
end
%5) GET MNI COORDINATES OF THE NEW FIGURE AND SAVE THEM ON DISK
MNI8 = zeros(3559,3);
for mm = 1:3559 %over brain voxel
    path_8mm = ['/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/parcel_80mm_3559ROIs/' num2str(mm) '.nii.gz']; %path for each of the 3559 parcels
    [mni_coord,pkfo] = osl_mnimask2mnicoords(path_8mm);  %getting MNI coordinates
    MNI8(mm,:) = mni_coord; %storing MNI coordinates
end
%saving on disk
save('/scratch7/MINDLAB2017_MEG-LearningBach/DTI_Portis/Templates/MNI152_8mm_coord_dyi.mat', 'MNI8');

%% CONVERSION T1 - DICOM TO NIFTI

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/osl/dicm2nii'); %adds path to the dcm2nii folder in osl
MRIsubj = dir('/projects/MINDLAB2020_MEG-AuditoryPatternRecognition/raw/0*');
MRIoutput = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/MRI_nifti';
MRIout_block{1} = 'Block_3'; MRIout_block{2} = 'Block_4'; MRIout_block{3} = 'Block_5'; MRIout_block{4} = 'Block_6'; MRIout_block{5} = 'Block_7';

for bb = 1:5 %over experimental blocks
    for ii = 1:length(MRIsubj) %over subjects
        asd = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/MRI_nifti/Block_' num2str(bb+2) '/' MRIsubj(ii).name];
        if ~exist(asd,'dir') %checking whether the directory exists
            mkdir(asd); %if not, creating it
        end
        if isempty(dir([asd '/*.nii'])) %if there are no nifti images.. I need to convert them
            flagg = 0;
            MRIMEGdate = dir([MRIsubj(ii).folder '/' MRIsubj(ii).name '/20*']);
            niiFolder = [MRIoutput '/' MRIout_block{bb} '/' MRIsubj(ii).name];
            for jj = 1:length(MRIMEGdate) %over dates of recording
                if ~isempty(dir([MRIMEGdate(jj).folder '/' MRIMEGdate(jj).name '/MR*'])) %if we get an MRI recording
                    MRI2 = dir([MRIMEGdate(jj).folder '/' MRIMEGdate(jj).name '/MR/*fatsat']); %looking for T1
                    if ~isempty(MRI2) %if we have it
                        flagg = 1; %determining that I could convert MRI T1
                        dcmSource = [MRI2(1).folder '/' MRI2(1).name '/files/'];
                        if ii ~= 68 || jj ~= 3 %this is because subject 0068 got two MRIs stored.. but the second one (indexed by jj = 3) is of another subject (0086); in this moment, subject 0086 is perfectly fine, but in subject 0068 there are still the two MRIs (for 0068 (jj = 2) and for 0086 (jj = 3))
                            dicm2nii(dcmSource, niiFolder, '.nii');
                        end
                    end
                end
            end
            if flagg == 0
                warning(['subject ' MRIsubj(ii).name ' has no MRI T1']);
            end
        end
        disp(ii)
    end
end

%% SETTING FOR CLUSTER (PARALLEL COMPUTING)

% clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 1); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% RHINO coregistration

%block to be run RHINO coregistrartion on
block = 6; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd

if block == 3
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 4
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogspeed*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 5
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*samemel*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 6
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*visualpat*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 7
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*pd*_tsssdsm.mat'); %dir to epoched files (encoding)
end

%running rhino
%OBS! check that all MEG data are in the same order and number as MRI nifti files!
a = ['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/MRI_nifti/Block_' num2str(block)]; %set path to MRI subjects' folders
for ii = 35%1:length(list) %OBS! change this depending on atonal vs. major
    S = [];
    S.ii = ii;
    S.D = [list(ii).folder '/' list(ii).name]; %path to major files
    D = spm_eeg_load(S.D);
    if ~isfield(D,'inv') %checking if the coregistration was already run
        dummyname = D.fname;
        if 7 == exist([a '/' dummyname(13:16)],'dir') %if you have the MRI folder
            dummymri = dir([a '/' dummyname(13:16) '/*.nii']); %path to nifti files (ending with .nii)
            if ~isempty(dummymri)
                S.mri = [dummymri(1).folder '/' dummymri(1).name];
                %standard parameters
                S.useheadshape = 1;
                S.use_rhino = 1; %set 1 for rhino, 0 for no rhino
                %         S.forward_meg = 'MEG Local Spheres';
                S.forward_meg = 'Single Shell'; %CHECK WHY IT SEEMS TO WORK ONLY WITH SINGLE SHELL!!
                S.fid.label.nasion = 'Nasion';
                S.fid.label.lpa = 'LPA';
                S.fid.label.rpa = 'RPA';
                jobid = job2cluster(@coregfunc,S); %running with parallel computing
            else
                warning(['subject ' dummyname(13:16) ' does not have the MRI'])
            end
        end
    else
        if isempty(D.inv{1}) %checking whether the coregistration was run but now it is empty..
            warning(['subject ' D.fname ' has an empty rhino..']);
        end
    end
    disp(ii)
end

%% checking (or copying) RHINO

copy_label = 0; % 1 = pasting inv RHINO from epoched data (where it was computed) to continuous data; 0 = simply showing RHINO coregistration
%block to be run RHINO coregistration on
block = 3; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd

if block == 3
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 4
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogspeed*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 5
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*samemel*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 6
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*visualpat*_tsssdsm.mat'); %dir to epoched files (encoding)
elseif block == 7
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*pd*_tsssdsm.mat'); %dir to epoched files (encoding)
end

for ii = 1:length(list)
    D = spm_eeg_load([list(ii).folder '/' list(ii).name]);
    if isfield(D,'inv')
        if copy_label == 0 %simply displaying RHINO coregistration
            if isfield(D,'inv') %checking if the coregistration was already run
                rhino_display(D)
            end
        else %pasting inv RHINO from epoched data (where it was computed) to continuous data
            inv_rhino = D.inv;
            D2 = spm_eeg_load([list(ii).folder '/' list(ii).name(2:end)]);
            D2.inv = inv_rhino;
            D2.save();
        end
    end
    disp(['Block ' num2str(block) ' - Subject ' num2str(ii)])
end

%%

%% BEAMFORMING

%% SETTING FOR CLUSTER (PARALLEL COMPUTING)

% clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 1); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 1); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% FUNCTION FOR SOURCE RECONSTRUCTION

%user settings
clust_l = 1; %1 = using cluster of computers (CFIN-MIB, Aarhus University); 0 = running locally
timek = 1:1026; %time-points
freqq = []; %frequency range (empty [] for broad band)
% freqq = [0.1 1]; %frequency range (empty [] for broad band)
% freqq = [2 8]; %frequency range (empty [] for broad band)
sensl = 1; %1 = magnetometers only; 2 = gradiometers only; 3 = both magnetometers and gradiometers (SUGGESTED 1!)
workingdir2 = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD'; %high-order working directory (a subfolder for each analysis with information about frequency, time and absolute value will be created)
block = 3; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd
invers = 1; %1-4 = different ways (e.g. mean, t-values, etc.) to aggregate trials and then source reconstruct only one trial; 5 for single trial independent source reconstruction

if isempty(freqq)
    absl = 0; % 1 = absolute value of sources; 0 = not
else
    absl = 0;
end
%actual computation
%list of subjects with coregistration (RHINO - OSL/FSL) - epoched
if block == 3
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogminor*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Old_Correct','New_T1_Correct','New_T2_Correct','New_T3_Correct','New_T4_Correct'};
    load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat');
elseif block == 4
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*recogspeed*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Old_Fast_Correct','Old_Slow_Correct','New_Fast_Correct','New_Slow_Correct'};
elseif block == 5
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*samemel*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Encoding','Old_Correct','New_Correct'};
elseif block == 6
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*visualpat*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'Encoding','Old_Correct','New_Correct'};
elseif block == 7
    list = dir ('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/e*pd*_tsssdsm.mat'); %dir to epoched files (encoding)
    condss = {'pd','pm'};
end
if isempty(freqq)
    workingdir = [workingdir2 '/Block_' num2str(block) '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_broadband_invers_' num2str(invers)];
else
    workingdir = [workingdir2 '/Block_' num2str(block) '/Beam_abs_' num2str(absl) '_sens_' num2str(sensl) '_freq_' num2str(freqq(1)) '_' num2str(freqq(2)) '_invers_' num2str(invers)];
end
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing');
if ~exist(workingdir,'dir') %creating working folder if it does not exist
    mkdir(workingdir)
end
for ii = 1:length(list) %over subjects
    S = [];
    if ~isempty(freqq) %if you want to apply the bandpass filter, you need to provide continuous data
        %             disp(['copying continuous data for subj ' num2str(ii)])
        %thus pasting it here
        %             copyfile([list_c(ii).folder '/' list_c(ii).name],[workingdir '/' list_c(ii).name]); %.mat file
        %             copyfile([list_c(ii).folder '/' list_c(ii).name(1:end-3) 'dat'],[workingdir '/' list_c(ii).name(1:end-3) 'dat']); %.dat file
        %and assigning the path to the structure S
        S.norm_megsensors.MEGdata_c = [list(ii).folder '/' list(ii).name(2:end)];
    end
    %copy-pasting epoched files
    %         disp(['copying epoched data for subj ' num2str(ii)])
    %         copyfile([list(ii).folder '/' list(ii).name],[workingdir '/' list(ii).name]); %.mat file
    %         copyfile([list(ii).folder '/' list(ii).name(1:end-3) 'dat'],[workingdir '/' list(ii).name(1:end-3) 'dat']); %.dat file
    
    S.Aarhus_cluster = clust_l; %1 for parallel computing; 0 for local computation
    
    S.norm_megsensors.zscorel_cov = 1; % 1 for zscore normalization; 0 otherwise
    S.norm_megsensors.workdir = workingdir;
    S.norm_megsensors.MEGdata_e = [list(ii).folder '/' list(ii).name];
    S.norm_megsensors.freq = freqq; %frequency range
    S.norm_megsensors.forward = 'Single Shell'; %forward solution (for now better to stick to 'Single Shell')
    
    S.beamfilters.sensl = sensl; %1 = magnetometers; 2 = gradiometers; 3 = both MEG sensors (mag and grad) (SUGGESTED 3!)
    S.beamfilters.maskfname = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz'; % path to brain mask: (e.g. 8mm MNI152-T1: '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_T1_8mm_brain.nii.gz')
    
    S.inversion.znorml = 0; % 1 for inverting MEG data using the zscored normalized one; (SUGGESTED 0 IN BOTH CASES!)
    %                                 0 to normalize the original data with respect to maximum and minimum of the experimental conditions if you have both magnetometers and gradiometers.
    %                                 0 to use original data in the inversion if you have only mag or grad (while e.g. you may have used zscored-data for covariance matrix)
    %
    S.inversion.timef = timek; %data-points to be extracted (e.g. 1:300); leave it empty [] for working on the full length of the epoch
    S.inversion.conditions = condss; %cell with characters for the labels of the experimental conditions (e.g. {'Old_Correct','New_Correct'})
    S.inversion.bc = [1 26]; %extreme time-samples for baseline correction (leave empty [] if you do not want to apply it)
    S.inversion.abs = absl; %1 for absolute values of sources time-series (recommendnded 1!)
    S.inversion.effects = invers;
    
    S.smoothing.spatsmootl = 0; %1 for spatial smoothing; 0 otherwise
    S.smoothing.spat_fwhm = 100; %spatial smoothing fwhm (suggested = 100)
    S.smoothing.tempsmootl = 0; %1 for temporal smoothing; 0 otherwise
    S.smoothing.temp_param = 0.01; %temporal smoothing parameter (suggested = 0.01)
    S.smoothing.tempplot = [1 2030 3269]; %vector with sources indices to be plotted (original vs temporally smoothed timeseries; e.g. [1 2030 3269]). Leave empty [] for not having any plot.
    
    S.nifti = 1; %1 for plotting nifti images of the reconstructed sources of the experimental conditions
    S.out_name = ['SUBJ_' list(ii).name(13:16)]; %name (character) for output nifti images (conditions name is automatically detected and added)
    
    if clust_l ~= 1 %useful  mainly for begugging purposes
        MEG_SR_Beam_LBPD(S);
    else
        jobid = job2cluster(@MEG_SR_Beam_LBPD,S); %running with parallel computing
    end
end

%% SETTING FOR CLUSTER (PARALLEL COMPUTING)

% clusterconfig('scheduler', 'none'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('scheduler', 'cluster'); %If you do not want to submit to the cluster, but simply want to test the script on the hyades computer, you can instead of 'cluster', write 'none'
clusterconfig('long_running', 0); % This is the cue we want to use for the clsuter. There are 3 different cues. Cue 0 is the short one, which should be enough for us
clusterconfig('slot', 3); %slot is memory, and 1 memory slot is about 8 GB. Hence, set to 2 = 16 GB
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing')

%% COMPUTING MEAN ACROSS PARTICIPANTS (USING PARALLEL COMPUTING, PARTICULARLY USEFUL IF YOU HAVE SEVERAL CONTRASTS)

block = 3; % 3 = recogminor; 4 = recogspeed; 5 = samemel; 6 = visualpat; 7 = pd
clust = 1; % 1 = using Aarhus cluster (parallel computing); 0 = run locally
analys_n = 2; %analysis number (in the list indexed below)

%building structure
asd = dir(['/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_' num2str(block) '/Beam*']);
S = [];
% dumm = dir([asd(analys_n).folder '/' asd(analys_n).name '/SUBJ*mat']);
% S.list = dumm(1:10);
S.workingdir = [asd(analys_n).folder '/' asd(analys_n).name '/test']; %path where the data from MEG_SR_Beam_LBPD.m is stored
S.sensl = 1; % 1 = magnetometers only; 2 = gradiometers only; 3 = both magnetometers and gradiometers.
S.plot_nifti = 1; %1 to plot nifti images; 0 otherwise
S.plot_nifti_name = []; %character with name for nifti files (it may be useful if you run separate analysis); Leave empty [] to not  specify any name
% S.contrast = [1 0 0 0 0 0 -1; 0 1 0 0 0 0 -1; 0 0 1 0 0 0 -1; 0 0 0 1 0 0 -1; 0 0 0 0 1 0 -1; 0 0 0 0 0 1 -1; 1 1 1 1 1 1 -1]; %one row per contrast (e.g. having 3 conditions, [1 -1 0; 1 -1 -1; 0 1 -1]; two or more -1 or 1 are interpreted as the mean over them first and then the contrast. Leave empty [] for no contrasts. 
if block == 3
    S.contrast = [1 -1 0 0 0; 1 0 -1 0 0; 1 0 0 -1 0; 1 0 0 0 -1];
elseif block == 4
        S.contrast = [1 0 -1 0; 0 1 0 -1];
elseif block == 5 || block == 6
    S.contrast = [0 1 -1; -1 1 0];
else
    S.contrast = [1 -1];
end
S.effects = 1; %mean over subjects for now
if clust == 1
    S.Aarhus_clust = 1; %1 to use paralle computing (Aarhus University, contact me, Leonardo Bonetti, for more information; leonardo.bonetti@clin.au.dk)
    %actual function
    jobid = job2cluster(@MEG_SR_Stats1_Fast_LBPD,S); %running with parallel computing
else
    S.Aarhus_clust = 0; %1 to use paralle computing (Aarhus University, contact me, Leonardo Bonetti, for more information; leonardo.bonetti@clin.au.dk)
    MEG_SR_Stats1_Fast_LBPD(S)
end

%%

%%

%% *** PCA COMPUTATION ***

%%

% Please, note that the following code is for getting PCA computed on the data averaged across participants.
% To inspect the several tests that we made (including statistical testing) please, refer to the following script provided in the current directory: BroadNess_APR2020_2 

% If you wish to use our code for running PCA, you can safely refer to the function provided in this script: PCA_LBPD
% If you wish to compute more elaborated solutions (as well as for statistical testing), please refer to the procedures and functions described in: BroadNess_APR2020_2 


%% setting up the cluster

addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Cluster_ParallelComputing') %add the path to the function that submits the jobs to the cluster
clusterconfig('scheduler', 'cluster');
clusterconfig('long_running', 1); %there are different queues for the cluster depending on the number and length of the jobs you want to submit 
clusterconfig('slot', 1); %slot in the queu

%% LBPD functions

%1) Add LBPD functions
% starting up some of the functions that I wrote for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl); % add to path LBPD and other OSL functions

%% LOADING AVERAGE OVER SUBJECTS

%AGGREGATED SUBJECTS

%loading data in 3559-voxel space (8mm)
%getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/sources_main_effects.mat');

%actual computation
%adjusting polarity
timex = 45:52;
vect = zeros(3559,1);
for jj = 1:3559 %over brain voxels
    if squeeze(mean(t_val_s(jj,timex,1),2)) > 0 %if the data in voxel jj is positive during N100 time
        vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
    else
        vect(jj,1) = 1;
    end
end
dum = zeros(size(t_val_s,1),size(t_val_s,2),size(t_val_s,3));
for cc = 1:size(dum,3) %over conditions
    for jj = 1:size(t_val_s,1) %over brain voxels
        dum(jj,:,cc) = t_val_s(jj,:,cc) .* vect(jj,1); %reversing (or not)..
        disp(['condition ' num2str(cc) ' - source ' num2str(jj)])
    end
end

% load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband_invers_5/GenEig_AverageOverSubjects.mat');

%% PRINCIPAL COMPONENT ANALYSIS

load('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/time_normal.mat');

S = [];
S.H = dum(:,:,1:5);
S.permnum = 1;
S.fig_l = 1;
S.sign_eig = '0';
S.namenii = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Test_PCA/Bor';
S.time = time;
S.rand_l = 1;
S.onefig = 0;

[ OUT ] = PCA_LBPD( S )

wcoeff = OUT.W;

%% TIME SERIES OF SINGLE SUBJECTS USING PCA WEIGHTS ON AVERAGED SUBJECTS DATA (this requires the two sections above to be run first)

% This is just to have a further confirmation that everything worked fine, but it is not the actual statistical testing

list = dir('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/after_maxfilter_v2/Source_LBPD/Block_3/Beam_abs_0_sens_1_freq_broadband/SUBJ*.mat');
J = zeros(size(wcoeff,2),3,5,length(list));

for ii = 1:length(list)
    load([list(ii).folder '/' list(ii).name])
    t_val_s = OUT.sources_ERFs;
    dum = zeros(size(t_val_s,1),size(t_val_s,2),size(t_val_s,3));
    for cc = 1:size(dum,3) %over conditions
        for jj = 1:size(t_val_s,1) %over brain voxels
            dum(jj,:,cc) = t_val_s(jj,:,cc) .* vect(jj,1); %reversing (or not)..
            %            disp(['condition ' num2str(cc) ' - source ' num2str(jj)])
        end
    end
    for iii = 1:5 %over experimental conditions
        J(:,:,iii,ii) = dum(:,1:size(wcoeff,2),iii)' * wcoeff(:,1:3); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
    end
    disp(ii)
end

%%

%%
