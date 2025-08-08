%%

%% BROAD-NESS - (AUDITORY PATTERN RECOGNITION 2020) - LEONARDO BONETTI - REVISION I - ADVANCED SCIENCE

%%

%% *** START UP FUNCTIONS.. (LBPD_startup_D) ***

%starting up some functions for LBPD toolbox.

%starting up some of the functions that I wrote for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl);

%%

%% CORRELATING TIME SERIES OF BRAIN NETWORKS WITH BEHAVIORAL INFORMATION (ACCURACY AND REACTION TIMES IN MEG TASK)

xlsx_dir_behav = '/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/BehavioralTaskMEG/Version_2/Final_xlsx/'; %dir to MEG behavioral results (.xlsx files)
bb = 3;
list_beh = dir([xlsx_dir_behav 'Subj*' num2str(bb) '.xlsx']); %normal blocks
bl = bb;
if bl == 3
    Block_3 = cell(length(list_beh)+1,13);
elseif bl == 4
    Block_4 = cell(length(list_beh)+1,11);
elseif bl == 5
    Block_5 = cell(length(list_beh)+1,7);
elseif bl == 6
    Block_6 = cell(length(list_beh)+1,7);
end
for ii = 1:length(list_beh) %over subjects for block bb
    %barbaric solution.. to build the name to be read for the excel files with the MEG behavioral tasks performance
    [~,~,raw_recog] = xlsread([list_beh(ii).folder '/' list_beh(ii).name]); %excel files
    %picking the current block
    %legend
    Block_3{1,1} = 'Subject'; Block_3{1,2} = 'OLD_Cor'; Block_3{1,3} = 'OLD_Cor %'; Block_3{1,4} = 'New_T1_Cor'; Block_3{1,5} = 'New_T1_Cor %'; %1st row
    Block_3{1,6} = 'New_T2_Cor'; Block_3{1,7} = 'New_T2_Cor %'; Block_3{1,8} = 'New_T3_Cor'; Block_3{1,9} = 'New_T3_Cor %'; Block_3{1,10} = 'New_T4_Cor'; Block_3{1,11} = 'New_T4_Cor %'; Block_3{1,12} = 'No response'; Block_3{1,13} = 'No response %'; %1st row
    Block_3{1,14} = 'OLD_RT'; Block_3{1,15} = 'NEWT1_RT'; Block_3{1,16} = 'NEWT2_RT'; Block_3{1,17} = 'NEWT3_RT'; Block_3{1,18} = 'NEWT4_RT';
    nr = 0; old = 0; n1 = 0; n2 = 0; n3 = 0; n4 = 0;
    ort = []; n1rt = []; n2rt = []; n3rt = []; n4rt = [];
    for k = 1:135
        if raw_recog{(k + 1),3} == 0 %if there was no response
            nr = nr + 1;
        elseif strcmp(raw_recog{(k + 1),2}(8:9),'ol') && raw_recog{(k + 1),3} == 1 %old correct
            old = old + 1;
            ort = cat(1,ort,raw_recog{(k + 1),4});
        elseif strcmp(raw_recog{(k + 1),2}(14:15),'t1') && raw_recog{(k + 1),3} == 2 %new t1 correct
            n1 = n1 + 1;
            n1rt = cat(1,n1rt,raw_recog{(k + 1),4});
        elseif strcmp(raw_recog{(k + 1),2}(14:15),'t2') && raw_recog{(k + 1),3} == 2 %new t2 correct
            n2 = n2 + 1;
            n2rt = cat(1,n2rt,raw_recog{(k + 1),4});
        elseif strcmp(raw_recog{(k + 1),2}(14:15),'t3') && raw_recog{(k + 1),3} == 2 %new t3 correct
            n3 = n3 + 1;
            n3rt = cat(1,n3rt,raw_recog{(k + 1),4});
        elseif strcmp(raw_recog{(k + 1),2}(14:15),'t4') && raw_recog{(k + 1),3} == 2 %new t4 correct
            n4 = n4 + 1;
            n4rt = cat(1,n4rt,raw_recog{(k + 1),4});
        end
    end
    disp(num2str(['Block ' num2str(bb) ' - Subject ' num2str(ii)]))
    Block_3{ii+1,1} = list_beh(ii).name(6:9); Block_3{ii+1,2} = old; Block_3{ii+1,3} = (old/27)*100; Block_3{ii+1,4} = n1; Block_3{ii+1,5} = (n1/27)*100;
    Block_3{ii+1,6} = n2; Block_3{ii+1,7} = (n2/27)*100; Block_3{ii+1,8} = n3; Block_3{ii+1,9} = (n3/27)*100; Block_3{ii+1,10} = n4; Block_3{ii+1,11} = (n4/27)*100; Block_3{ii+1,12} = nr; Block_3{ii+1,13} = (nr/135)*100;
    Block_3{ii+1,14} = mean(ort); Block_3{ii+1,15} = mean(n1rt); Block_3{ii+1,16} = mean(n2rt); Block_3{ii+1,17} = mean(n3rt); Block_3{ii+1,18} = mean(n4rt);
end
Block_3_t = cell2table(Block_3); %remove the possible empty cell


%preparing information for behavioral analysis
accur_idx = [2:2:10];
RT_idx = [14:18];
reduced_index = find(time >= starting_time & time <= ending_time);
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/figure2_TEST/timeserie_wcoeff_from_maineffect.mat');
dataa = permute(J, [2,1,3,4]);
bumba = cell(2,5,3,3);

%accuracy
rho = zeros(size(dataa,1),776,length(accur_idx)); %PCs x time x experimental conditions
p = zeros(size(dataa,1),776,length(accur_idx)); %PCs x time x experimental conditions
for cc = 1:size(dataa,1) %over PCs
    for ii = 1:length(accur_idx) %over experimental conditions
        for tt = 1:776 %time for beginning to 2.5 seconds
            beha = cell2mat(Block_3(2:end,accur_idx(ii))); %extracting inforation about musical training
            datadum = squeeze(dataa(cc,tt,ii,~isnan(beha))); %extracting neural data for brain network cc and experimental condition ii (making sure there are no NaNs (which is not the case but this way of writing can be useful for the future))
            [rho(cc,tt,ii),p(cc,tt,ii)] = corr(datadum,beha(~isnan(beha))); %correlation
        end
       disp(ii)
   end
end
for PC = 1:2
    for cond = 1:5
        [pthr] = fdr(p(PC,:,cond))
        bumba{PC,cond,1,1} = time(find(p(PC,:,cond)<pthr));
        bumba{PC,cond,1,2} = find(p(PC,:,cond)<pthr);
        bumba{PC,cond,1,3} = cat(1,time(find(p(PC,:,cond)<pthr)),rho(PC,find(p(PC,:,cond)<pthr),cond));
    end
end

%RTs
rho = zeros(size(dataa,1),776,length(RT_idx)); %PCs x time x experimental conditions
p = zeros(size(dataa,1),776,length(RT_idx)); %PCs x time x experimental conditions
for cc = 1:size(dataa,1) %over PCs
    for ii = 1:length(RT_idx) %over experimental conditions
        for tt = 1:776 %time for beginning to 2.5 seconds
            beha = cell2mat(Block_3(2:end,RT_idx(ii))); %extracting inforation about musical training
            datadum = squeeze(dataa(cc,tt,ii,~isnan(beha))); %extracting neural data for brain network cc and experimental condition ii (making sure there are no NaNs (which is not the case but this way of writing can be useful for the future))
            [rho(cc,tt,ii),p(cc,tt,ii)] = corr(datadum,beha(~isnan(beha))); %correlation
        end
       disp(ii)
   end
end
for PC = 1:2
    for cond = 1:5
        [pthr] = fdr(p(PC,:,cond))
        bumba{PC,cond,2,1} = time(find(p(PC,:,cond)<pthr));
        bumba{PC,cond,2,2} = find(p(PC,:,cond)<pthr);
        bumba{PC,cond,2,3} = cat(1,time(find(p(PC,:,cond)<pthr)),rho(PC,find(p(PC,:,cond)<pthr),cond));
    end
end

%musical training
[~,~,raw] = xlsread('/scratch7/MINDLAB2020_MEG-AuditoryPatternRecognition/leonardo/Papers/Systematic_Variations_DCM/NatureCommunications_Revision/Revision_I/MusicalTrainingPracticeWM.xlsx');
beha = cell2mat(raw(2:end,3)); %extracting inforation about musical training
rho = zeros(size(dataa,1),776,length(RT_idx)); %PCs x time x experimental conditions
p = zeros(size(dataa,1),776,length(RT_idx)); %PCs x time x experimental conditions
for cc = 1:size(dataa,1) %over PCs
    for ii = 1:length(RT_idx) %over experimental conditions
        for tt = 1:776 %time for beginning to 2.5 seconds
            datadum = squeeze(dataa(cc,tt,ii,~isnan(beha))); %extracting neural data for brain network cc and experimental condition ii (making sure there are no NaNs (which is not the case but this way of writing can be useful for the future))
            [rho(cc,tt,ii),p(cc,tt,ii)] = corr(datadum,beha(~isnan(beha))); %correlation
        end
       disp(ii)
   end
end
for PC = 1:2
    for cond = 1:5
        [pthr] = fdr(p(PC,:,cond))
        bumba{PC,cond,3,1} = time(find(p(PC,:,cond)<pthr));
        bumba{PC,cond,3,2} = find(p(PC,:,cond)<pthr);
        bumba{PC,cond,3,3} = cat(1,time(find(p(PC,:,cond)<pthr)),rho(PC,find(p(PC,:,cond)<pthr),cond));
    end
end


%% inspecting results


clc
PC = 1;
cond = 5;

bumba = time(find(p(PC,:,cond)<0.05));
bumba = cat(1,bumba,rho(PC,find(p(PC,:,cond)<0.05),cond))

[pthr] = fdr(p(PC,:,cond))
bumba = time(find(p(PC,:,cond)<pthr))
bumba = cat(1,bumba,rho(PC,find(p(PC,:,cond)<pthr),cond))


%%
