%% SouConTACS Analyses
% the script to perform preprocessing and all subsequent analysis on the
% Source-Confidence tACS dataset
clearvars;close all

%%%%%%%%%%%%%%%%%%%%%%
%% GLOBAL VARIABLES %%
%%%%%%%%%%%%%%%%%%%%%%
% list of preprocessing variations to be done: 
% - eeg (EEG)
% - beh (behavioral)
curexperiment.preproc                   = {'eeg'}; % 'eeg' or 'beh'
% list of analyses to be done: 
% - ana (single-subject analysis)
% - pow (spectral analysis (& phase-amplitude coupling, power envelope correlation, phase coherence)
% - ctl (control analysis)
% - beh (behavioral analysis)
% - ga (grand averages)
% - plt (plotting of TFR data)
% - erp (ERP analysis)
% - erpga (ERP grand average) 
curexperiment.analyses                  = {'ana','pow'}; % {'ana','erp','pow','plt','ga','erpga'}; 
% which output to obtain from the fourier transform, the power only or the complex numbers.
% For 'pow' low and high frequencies are treated differently, see below.
% For 'fourier' all frequencies are treated the same, aka no multitapering.
% pow: 
% - time-frequency plotting in analysis level 1 or 2
% - when you just need the time-frequency data
% - getting trial-by-trial frequency data for analysis level 3
% - getting the baseline corrected peak amplitudes from the time-frequency data in level 2
% fourier: 
% - getting frequency data in analysis level 2
% - getting phase amplitude coupling in analysis levels 2 & 3
% - getting power envelope correlation in analysis level 2
% - getting phase coherence in level 2
% This experiment is split up into two papers, depedent on which paper, certain decisions are altered.
% In general, anaylsis level 1 is used for plotting, level 2 is used for the tACS paper analyses, and 3 is used for the eeg paper analyses
curexperiment.expvariant                = 'eeg'; % 'eeg' or 'tacs'

%%%%%%%%%%%%%%%%%%%%%%
curexperiment.name                      = 'SouConTACS'; % current experiment
curexperiment.manualcheck               = 3; % determine if artifact rejection and ICA need to be done checked (both=1, only artifacts=2, only ICA=3, none=4)
curexperiment.dirroot                   = ''; % set current base directory
curexperiment.scriptdir                 = ''; % directory of the script
curexperiment.datafolder_input          = [curexperiment.dirroot,curexperiment.name,'/Raw']; % location of the EEG data
curexperiment.datafolder_inputbehav     = [curexperiment.dirroot,curexperiment.name,'/BehavData']; % location of the behavioral data
curexperiment.datafolder_output         = [curexperiment.dirroot,curexperiment.name,'/Output']; % location of outputfiles
curexperiment.analysis_loc              = fullfile(curexperiment.datafolder_output, sprintf('%s_Analyses',curexperiment.expvariant)); % location of analyses
curexperiment.outputfile                = fullfile(curexperiment.datafolder_output, sprintf('_%s_Stats.mat',curexperiment.name));
curexperiment.extension                 = {'*.eeg','*.set'}; % extension of the EEG files
curexperiment.marker_offset             = 0; % there is no offset in the EEG markers
curexperiment.Nses                      = 4; % number of sessions
curexperiment.stims                     = {'nostim','theta','gamma','sham'}; % the stimulations
curexperiment.eventtype                 = {'Response','Stimulus','trigger'}; % epoch event type
curexperiment.fs_org                    = 10000; % original sampling frequency
curexperiment.fs_ds                     = 1000; % downsampled sampling frequency
curexperiment.elec.lay                  = 'easycap-M1.txt'; % EEG template, which will be loaded as a layout in curexperiment.elecs. This will have a few additional channels that are not used
curexperiment.elec.impref               = 'Fz'; % implicit reference (non-recorded)
curexperiment.elec.newref               = 'EEG'; % desired new/offline reference: average reference
curexperiment.Nextelectrodes            = 0; % number of external electrodes
curexperiment.filtype                   = 'but'; % filter type; I tried 'fir' but it took ages
curexperiment.bs_freq                   = [58 62; 118 122; 178 182]; % notch filter alternative
curexperiment.hp_freq                   = .1; % high-pass filter
curexperiment.lp_freq                   = 100; % low-pass filter
curexperiment.Ntrials_enc               = 200; % number of encoding trials
curexperiment.Ntrials_ret               = 400; % number of retrieval trials
curexperiment.datasets_names            = {'data_enc','data_ret','data_rest'}; % dataset variable names
curexperiment.dataset_name              = {'_Enc','_Ret','_Rest'}; % dataset name for naming files
curexperiment.subject_groups            = 1; % number of subject groups
curexperiment.Nsubs                     = 54; % number of subjects
curexperiment.subjectsWrong             = {}; % participants that switched up old/new responses
curexperiment.do_art                    = false; % do (not) do artifact rejection
curexperiment.do_ica                    = true; % do (not) do ICA
curexperiment.Nanalyses.erp             = 99; %amount of outputfiles
curexperiment.Nanalyses.tp              = 99; % amount of outputfiles
curexperiment.Nanalyses.plt             = 0; % amount of outputfiles
curexperiment.Nanalyses.con             = 999; % amount of outputfiles
curexperiment.erp.basewin               = [-.2 0]; % erp baseline: -200 - 0 ms
curexperiment.erp.time_interest         = [.2 .4;... % erp timewindows: 200 - 400 ms (frontal)
                                           .5 .9]; % 500 - 900 ms (parietal)
curexperiment.pow.curpow                = {'_Total'}; % set the current power type of interest (total, induced, evoked)
curexperiment.pow.basewin               = [-.5 -.25]; % tfr baseline: -500 - -250 ms (because of the sliding time window)
curexperiment.pow.freq_interest         = 1:1:57; % frequencies of interest
curexperiment.pow.timwin                = 0.5.*ones(size(curexperiment.pow.freq_interest)); % length of sliding timewindow
curexperiment.pow.taptype               = 'hanning'; % single taper
curexperiment.pow.low.freq_interest     = 1:1:29; % frequencies of interest
curexperiment.pow.low.timwin            = 0.5.*ones(size(curexperiment.pow.low.freq_interest)); % length of sliding timewindow
curexperiment.pow.low.taptype           = 'hanning'; % single taper
curexperiment.pow.high.freq_interest    = 30:1:57; % frequencies of interest
curexperiment.pow.high.timwin           = 10./curexperiment.pow.high.freq_interest; % length of sliding timewindow, 10 cycles per time window
curexperiment.pow.high.tapsmo           = 0.2.*curexperiment.pow.high.freq_interest; % width of the frequency smoothing, total (symmetrical) smoothing is between 12 and 22.4 Hz
curexperiment.pow.high.taptype          = 'dpss'; % multi taper
if strcmp(curexperiment.expvariant,'eeg')
    curexperiment.pow.TbTFreqs          = [3,7;30,50]; % frequencies of interest in the TbT analysis
    curexperiment.pow.freqnms           = {'theta','gamma'}; % frequencies of interest in the TbT analysis
elseif strcmp(curexperiment.expvariant,'tacs')
    curexperiment.pow.TbTFreqs          = [1,7;45,55]; % frequencies of interest in the TbT analysis
end
curexperiment.pow.plot                  = true; % determine whether to do the single-subject TFR plotting
% as the electrodes were put on different positions on the head in sessions 2-4, we needed to relabel them
curexperiment.oldchans{1}               = {'Fp1';'F3';'F7';'FT9';'FC5';'FC1';'C3';'T7';'TP9';'CP5';'CP1';'Pz';'P3';'P7';'O1';'Oz';'O2';... % session 1
                                            'P4';'P8';'TP10';'CP6';'CP2';'Cz';'C4';'T8';'FT10';'FC6';'FC2';'F4';'F8';'Fp2';'AF7';'AF3';...
                                            'AFz';'F1';'F5';'FT7';'FC3';'FCz';'C1';'C5';'TP7';'CP3';'P1';'P5';'PO7';'PO3';'POz';'PO4';...
                                            'PO8';'P6';'P2';'CPz';'CP4';'TP8';'C6';'C2';'FC4';'FT8';'F6';'F2';'AF4';'AF8'};
curexperiment.oldchans{2}               = {'Fp1';'F3';'F7';'FT9';'FC5';'F5';'C3';'T7';'TP9';'P1';'CP4';'Pz';'P2';'POz';'O1';'Oz';'O2';'P4';'P8';... % session 2-4
                                            'TP10';'CP6';'P6';'F1';'C4';'T8';'FT10';'FC6';'FC3';'AF3';'F8';'FC4';'AF7';'AF3x';'AFz';...
                                            'F1x';'F5x';'FT7';'FC3x';'FCz';'C1';'C5';'TP7';'CP3';'P1x';'P5';'PO7';'PO3';'POzx';'PO4';...
                                            'PO8';'P6x';'P2x';'CPz';'CP4x';'TP8';'C6';'C2';'FC4x';'FT8';'F6';'F2';'AF4';'AF8';'tACS'};
curexperiment.newchans{1}               = {'Fp1';'F3';'F7';'FT9';'FC5';'FC1';'C3';'T7';'TP9';'CP5';'CP1';'Pz';'P3';'P7';'O1';'Oz';'O2';... % session 1
                                            'P4';'P8';'TP10';'CP6';'CP2';'Cz';'C4';'T8';'FT10';'FC6';'FC2';'F4';'F8';'Fp2';'AF7';'AF3';...
                                            'AFz';'F1';'F5';'FT7';'FC3';'FCz';'C1';'C5';'TP7';'CP3';'P1';'P5';'PO7';'PO3';'POz';'PO4';...
                                            'PO8';'P6';'P2';'CPz';'CP4';'TP8';'C6';'C2';'FC4';'FT8';'F6';'F2';'AF4';'AF8'};
curexperiment.newchans{2}               = {'Fp1';'F3';'F7';'FT9';'FC5';'F5';'C3';'T7';'TP9';'P1';'CP4';'Pz';'P2';'POz';'O1';'Oz';'O2';'P4';'P8';... % session 2-4
                                            'TP10';'CP6';'P6';'F1';'C4';'T8';'FT10';'FC6';'FC3';'AF3';'F8';'FC4';'tACS'};
if strcmp(curexperiment.expvariant,'eeg')
    curexperiment.chngrp.frontal            = { % grouping of channels
                                                'F5','F3','F1',... %left 
                                                'Fz',... %central
                                                'F6','F4','F2'}; %right
    curexperiment.chngrp.parietal           = { % grouping of channels
                                                'P5','P3','P1',... %left
                                                'Pz',... %central
                                                'P6','P4','P2'};
elseif strcmp(curexperiment.expvariant,'tacs')
    curexperiment.chngrp.frontal           = { % grouping of channels
                                                'AF4','AFz','AF8','F2','F4'};
    curexperiment.chngrp.parietal            = { % grouping of channels
                                                'P5','P3','P7','CP5','PO7'};
end
% make a table with the original encoding EEG markers
description.enc                         = {'Stimulus Encoding Pleasant','Stimulus Encoding Place'... % encoding
                                            'Response Unsuccessful','Response Part Successful','Response Successful','No Response','Response Onset'}; % encoding
original_marker.enc                     = {'S  2','S  3',... % encoding
                                            'R  1','R  2','R  3','R  9','R  8'}'; % encoding
count_without_practice.enc              = {100,100,...
                                            [],[],[],[],[]}';
cur_count.enc                           = zeros(length(description.enc),1);
curexperiment.original_markers.enc      = table(original_marker.enc,count_without_practice.enc,cur_count.enc,'RowNames',description.enc); curexperiment.original_markers.enc.Properties.VariableNames = {'original_marker' 'count_without_practice' 'cur_count'}; % the original encoding EEG markers
% make a table with the original retrieval & rest EEG markers
description.ret                         = {'Start/End RestEEG','Eyes Open','Eyes Closed',... % restEEG
                                            'Stimulus Retrieval New','Stimulus Retrieval Pleasant','Stimulus Retrieval Place'... % retrieval
                                            'Response Very Sure Old','Response Bit Sure Old','Response Unsure ON','Response Bit Sure New','Response Very Sure New','No Response',...; % retrieval
                                            'Response Very Sure Pleasant','Response Bit Sure Pleasant','Response Unsure Source','Response Bit Sure Place','Response Very Sure Place'}; % retrieval
original_marker.ret                     = {'R  6','R  7','R  8',... % restEEG
                                            'S  1','S  2','S  3'... % retrieval
                                            'R  1','R  2','R  3','R  4','R  5','R  9',... % retrieval
                                            'R 11','R 12','R 13','R 14','R 15'}'; % retrieval
count_without_practice.ret              = {4,4,4,... % restEEG
                                            200,100,100,... % retrieval
                                            [],[],[],[],[],[],...
                                            [],[],[],[],[]}'; % retrieval
cur_count.ret                           = zeros(length(description.ret),1);
curexperiment.original_markers.ret      = table(original_marker.ret,count_without_practice.ret,cur_count.ret,'RowNames',description.ret); curexperiment.original_markers.ret.Properties.VariableNames = {'original_marker' 'count_without_practice' 'cur_count'}; % the original retrieval and rest EEG markers
clear original_marker count_without_practice cur_count description
curexperiment.markers.encret            = {'S  1','S  2','S  3'}; % stimulus onset markers
curexperiment.stimtime(1,:)             = [.7 (1.5+1/curexperiment.fs_ds)]; % encoding epoch length
curexperiment.stimtime(3,:)             = [0 (60+1/curexperiment.fs_ds)]; % rest epoch length
% Conditions
% Encoding
% x... : (stimulus) 2=pleasant, 3=place
% .x.. : (enc resp) 1=unsuccessful, 2=partial successful, 3=successful, 9=no resp
% ..x. : (ON resp)  1=vs old, 2=bs old, 3=not sure 4=bs new, 5=vs new, 9=no resp
% ...x : (SO resp)  1=vs pleas, 2=bs pleas, 3=not sure 4=bs place, 5=vs place, 9=no resp
curexperiment.data1l1_name{1}           = 'sHit'; % item hit, source hit
curexperiment.data1.l1.condition1       = str2double(strrep(join(string([combvec(2,[1,2,3,9],[1,2],[1,2])';combvec(3,[1,2,3,9],[1,2],[4,5])'])),' ',''));                                         
curexperiment.data1l1_name{2}           = 'sMiss'; % item hit, source miss
curexperiment.data1.l1.condition2       = str2double(strrep(join(string([combvec(2,[1,2,3,9],[1,2],[3,4,5])';combvec(3,[1,2,3,9],[1,2],[1,2])'])),' ',''));                                         
curexperiment.data1l1_name{3}           = 'iHit'; % item hit, all source
curexperiment.data1.l1.condition3       = str2double(strrep(join(string(combvec([2,3],[1,2,3,9],[1,2],[1,2,3,4,5,9])')),' ',''));                     
curexperiment.data1l1_name{4}           = 'iMiss'; % item miss
curexperiment.data1.l1.condition4       = str2double(strrep(join(string(combvec([2,3],[1,2,3,9],[3,4,5],[1,2,3,4,5,9])')),' ',''));  
curexperiment.data1l1_name{5}           = 'sHC'; % source high-confidence
curexperiment.data1.l1.condition5       = str2double(strrep(join(string(combvec([2,3],[1,2,3,9],[1,2,3,4,5,9],[1,5])')),' ',''));                                               
curexperiment.data1l1_name{6}           = 'sLC'; % source low-confidence
curexperiment.data1.l1.condition6       = str2double(strrep(join(string(combvec([2,3],[1,2,3,9],[1,2,3,4,5,9],[2,3,4])')),' ',''));                                                                                      
curexperiment.data1l1_name{7}           = 'iHC'; % item high-confidence
curexperiment.data1.l1.condition7       = str2double(strrep(join(string(combvec([2,3],[1,2,3,9],[1,5],[1,2,3,4,5,9])')),' ',''));                   
curexperiment.data1l1_name{8}           = 'iLC'; % item low-confidence
curexperiment.data1.l1.condition8       = str2double(strrep(join(string(combvec([2,3],[1,2,3,9],[2,3,4],[1,2,3,4,5,9])')),' ','')); 
% Combinations
cmbs.c1 = {'Old'};cmbs.c2 = {'iCor'};cmbs.c3 = {'iHC','iLC'};cmbs.c4 = {'sCor','sInc'};cmbs.c5 = {'sHC','sLC'};
[cmbs.c1,cmbs.c2,cmbs.c3,cmbs.c4,cmbs.c5] = ndgrid(cmbs.c1,cmbs.c2,cmbs.c3,cmbs.c4,cmbs.c5);
M = [cmbs.c1(:),cmbs.c2(:),cmbs.c3(:),cmbs.c4(:),cmbs.c5(:)];
cmbs.c1 = {'Old'};cmbs.c2 = {'iInc'};cmbs.c3 = {'iHC','iLC'};cmbs.c4 = {''};cmbs.c5 = {''};
[cmbs.c1,cmbs.c2,cmbs.c3,cmbs.c4,cmbs.c5] = ndgrid(cmbs.c1,cmbs.c2,cmbs.c3,cmbs.c4,cmbs.c5);
M2 = [cmbs.c1(:),cmbs.c2(:),cmbs.c3(:),cmbs.c4(:),cmbs.c5(:)];
M = [M;M2];
Mlabel = arrayfun(@(x)strjoin(M(x,:),'_'),(1:size(M,1))','un',0);
Mlabel = strip(Mlabel,'right','_');
clear cmbs
for i=1:length(Mlabel)
    curexperiment.data1l2_name{i}  = Mlabel{i};
    cnd.c1_1        = 2;
    cnd.c1_2        = 3;
    cnd.c2          = [1,2,3,9];
    if all(contains(M(i,2:3),{'iCor','iHC'}))
        cnd.c3      = 1;
    elseif all(contains(M(i,2:3),{'iCor','iLC'}))
        cnd.c3      = 2;
    elseif all(contains(M(i,2:3),{'iInc','iLC'}))
        cnd.c3      = [3,4]; % guesses coded as incorrect
    elseif all(contains(M(i,2:3),{'iInc','iHC'}))
        cnd.c3      = 5;
    end
    if all(contains(M(i,4:5),{'sCor','sHC'}))
        cnd.c4_1    = 1;                cnd.c4_2    = 5;
    elseif all(contains(M(i,4:5),{'sCor','sLC'}))
        cnd.c4_1    = 2;                cnd.c4_2    = 4;
    elseif all(contains(M(i,4:5),{'sInc','sLC'}))
        cnd.c4_1    = [3,4];            cnd.c4_2    = [3,2]; % guesses coded as incorrect
    elseif all(contains(M(i,4:5),{'sInc','sHC'}))
        cnd.c4_1    = 5;                cnd.c4_2    = 1;
    else
        cnd.c4_1    = [1,2,3,4,5,9];    cnd.c4_2    = [1,2,3,4,5,9];
    end
    curexperiment.data1.l2.(['condition' num2str(i)]) = str2double(strrep(join(string([...
                                                            combvec(cnd.c1_1,cnd.c2,cnd.c3,cnd.c4_1)';...
                                                            combvec(cnd.c1_2,cnd.c2,cnd.c3,cnd.c4_2)'])),' ','')); 
    clear cnd
end
curexperiment.data1l3_name{1}           = 'All'; % all trials
curexperiment.data1.l3.condition1       = str2double(strrep(join(string(combvec([2,3],[1,2,3,9],[1,2,3,4,5,9],[1,2,3,4,5,9])')),' ',''));                     
% Retrieval
% x .  .   : (stimulus) 1=new, 2=pleasant, 3=place
% . x  .   : (response) 1=vs old, 2=bs old, 3=not sure 4=bs new, 5=vs new, 9=no resp
% . . (x)  : (response) 1=vs pleas, 2=bs pleas, 3=not sure 4=bs place, 5=vs place, 9=no resp
curexperiment.data2l1_name{1}           = 'sHit'; % item hit & source hit
curexperiment.data2.l1.condition1       = str2double(strrep(join(string([combvec(2,[1,2,3],[1,2])';combvec(3,[1,2,3],[4,5])'])),' ',''));                   
curexperiment.data2l1_name{2}           = 'sMiss'; % item hit & source miss
curexperiment.data2.l1.condition2       = str2double(strrep(join(string([combvec(2,[1,2,3],[3,4,5])';combvec(3,[1,2,3],[1,2])'])),' ',''));
curexperiment.data2l1_name{3}           = 'iHit'; % item hit (all source (non) responses)
curexperiment.data2.l1.condition3       = str2double(strrep(join(string(combvec([2,3],[1,2],[1,2,3,4,5,9])')),' ',''));  
curexperiment.data2l1_name{4}           = 'iMiss'; % item miss
curexperiment.data2.l1.condition4       = str2double(strrep(join(string(combvec([2,3],[3,4,5])')),' ',''));
curexperiment.data2l1_name{5}           = 'FA'; % (item) false alarm
curexperiment.data2.l1.condition5       = str2double(strrep(join(string(combvec(1,[1,2,3],[1,2,3,4,5,9])')),' ',''));
curexperiment.data2l1_name{6}           = 'CR'; % (item) correct rejection
curexperiment.data2.l1.condition6       = str2double(strrep(join(string(combvec(1,[4,5])')),' ',''));                      
curexperiment.data2l1_name{7}           = 'sHC'; % HC source
curexperiment.data2.l1.condition7       = str2double(strrep(join(string(combvec([2,3],[1,2,3],[1,5])')),' ','')); 
curexperiment.data2l1_name{8}           = 'sLC'; % LC source
curexperiment.data2.l1.condition8       = str2double(strrep(join(string(combvec([2,3],[1,2,3],[2,3,4])')),' ',''));  
curexperiment.data2l1_name{9}           = 'iHC'; % HC item 
curexperiment.data2.l1.condition9       = str2double(strrep(strrep(join(string([combvec([2,3],[1,5],[1,2,3,4,5,9])'; combvec(1,5,0)'])),' ',''),'0','')); 
curexperiment.data2l1_name{10}          = 'iLC'; % LC item
curexperiment.data2.l1.condition10      = str2double(strrep(strrep(join(string([combvec([2,3],[2,3,4],[1,2,3,4,5,9])'; combvec(1,[3,4],0)'])),' ',''),'0','')); 
curexperiment.data2l1_name{11}          = 'iHitHC'; % Hit HC item 
curexperiment.data2.l1.condition11      = str2double(strrep(join(string(combvec([2,3],1,[1,2,3,4,5,9])')),' ','')); 
curexperiment.data2l1_name{12}          = 'iHitLC'; % Hit LC item
curexperiment.data2.l1.condition12      = str2double(strrep(join(string(combvec([2,3],[2,3],[1,2,3,4,5,9])')),' ','')); 
curexperiment.data2l1_name{13}          = 'iCRHC'; % CR HC item 
curexperiment.data2.l1.condition13      = str2double(strrep(join(string(combvec(1,5)')),' ','')); 
curexperiment.data2l1_name{14}          = 'iCRLC'; % CR LC item
curexperiment.data2.l1.condition14      = str2double(strrep(join(string(combvec(1,[3,4])')),' ','')); 

% Combinations
cmbs.c1 = {'Old'};cmbs.c2 = {'iCor'};cmbs.c3 = {'iHC','iLC'};cmbs.c4 = {'sCor','sInc'};cmbs.c5 = {'sHC','sLC'};
[cmbs.c1,cmbs.c2,cmbs.c3,cmbs.c4,cmbs.c5] = ndgrid(cmbs.c1,cmbs.c2,cmbs.c3,cmbs.c4,cmbs.c5);
M = [cmbs.c1(:),cmbs.c2(:),cmbs.c3(:),cmbs.c4(:),cmbs.c5(:)];
cmbs.c1 = {'Old'};cmbs.c2 = {'iInc'};cmbs.c3 = {'iHC','iLC'};cmbs.c4 = {''};cmbs.c5 = {''};
[cmbs.c1,cmbs.c2,cmbs.c3,cmbs.c4,cmbs.c5] = ndgrid(cmbs.c1,cmbs.c2,cmbs.c3,cmbs.c4,cmbs.c5);
M2 = [cmbs.c1(:),cmbs.c2(:),cmbs.c3(:),cmbs.c4(:),cmbs.c5(:)];
cmbs.c1 = {'New'};cmbs.c2 = {'iCor'};cmbs.c3 = {'iHC','iLC'};cmbs.c4 = {''};cmbs.c5 = {''};
[cmbs.c1,cmbs.c2,cmbs.c3,cmbs.c4,cmbs.c5] = ndgrid(cmbs.c1,cmbs.c2,cmbs.c3,cmbs.c4,cmbs.c5);
M3 = [cmbs.c1(:),cmbs.c2(:),cmbs.c3(:),cmbs.c4(:),cmbs.c5(:)];
cmbs.c1 = {'New'};cmbs.c2 = {'iInc'};cmbs.c3 = {'iHC','iLC'};cmbs.c4 = {'sInc'};cmbs.c5 = {'sHC','sLC'};
[cmbs.c1,cmbs.c2,cmbs.c3,cmbs.c4,cmbs.c5] = ndgrid(cmbs.c1,cmbs.c2,cmbs.c3,cmbs.c4,cmbs.c5);
M4 = [cmbs.c1(:),cmbs.c2(:),cmbs.c3(:),cmbs.c4(:),cmbs.c5(:)];
M = [M;M2;M3;M4];
Mlabel = arrayfun(@(x)strjoin(M(x,:),'_'),(1:size(M,1))','un',0);
Mlabel = strip(Mlabel,'right','_');
for i=1:length(Mlabel)
    curexperiment.data2l2_name{i}  = Mlabel{i};
    cnd.c1_1        = 2;
    cnd.c1_2        = 3;
    cnd.c1_3        = 1;
    if contains(M(i,1),'Old') && contains(M(i,2),'iCor')
        if all(contains(M(i,3),{'iHC'}))
            cnd.c2      = 1;
        elseif all(contains(M(i,3),{'iLC'}))
            cnd.c2      = 2;
        end
        if all(contains(M(i,4:5),{'sCor','sHC'}))
            cnd.c3_1    = 1;        cnd.c3_2    = 5;
        elseif all(contains(M(i,4:5),{'sCor','sLC'}))
            cnd.c3_1    = 2;        cnd.c3_2    = 4;
        elseif all(contains(M(i,4:5),{'sInc','sLC'}))
            cnd.c3_1    = [3,4];    cnd.c3_2    = [3,2]; % guesses coded as incorrect
        elseif all(contains(M(i,4:5),{'sInc','sHC'}))
            cnd.c3_1    = 5;        cnd.c3_2    = 1;
        end
        curexperiment.data2.l2.(['condition' num2str(i)]) = str2double(strrep(join(string([...
                                                                combvec(cnd.c1_1,cnd.c2,cnd.c3_1)';...
                                                                combvec(cnd.c1_2,cnd.c2,cnd.c3_2)'])),' ',''));  
    elseif contains(M(i,1),'Old') && contains(M(i,2),'iInc')
        if all(contains(M(i,3),{'iLC'}))
            cnd.c2      = [3,4]; % guesses coded as incorrect
        elseif all(contains(M(i,3),{'iHC'}))
            cnd.c2      = 5;
        end
        curexperiment.data2.l2.(['condition' num2str(i)]) = str2double(strrep(join(string([...
                                                                combvec(cnd.c1_1,cnd.c2)';...
                                                                combvec(cnd.c1_2,cnd.c2)'])),' ','')); 
    elseif contains(M(i,1),'New') && contains(M(i,2),'iCor')
        if all(contains(M(i,3),{'iHC'}))
            cnd.c2      = 5;
        elseif all(contains(M(i,3),{'iLC'}))
            cnd.c2      = 4;
        end
        curexperiment.data2.l2.(['condition' num2str(i)]) = str2double(strrep(join(string(...
                                                                combvec(cnd.c1_3,cnd.c2)')),' ',''));
    elseif contains(M(i,1),'New') && contains(M(i,2),'iInc')
        if all(contains(M(i,3),{'iLC'}))
            cnd.c2      = [3,2]; % guesses coded as incorrect
        elseif all(contains(M(i,3),{'iHC'}))
            cnd.c2      = 1;
        end
        if all(contains(M(i,5),{'sLC'}))
            cnd.c3      = [2,3,4]; % guesses coded as incorrect
        elseif all(contains(M(i,5),{'sHC'}))
            cnd.c3      = [1,5];
        end
        curexperiment.data2.l2.(['condition' num2str(i)]) = str2double(strrep(join(string(...
                                                                combvec(cnd.c1_3,cnd.c2,cnd.c3)')),' ',''));
        % as I was silly and coded new & not sure as 13, I need to adjust for this
        curexperiment.data2.l2.(['condition' num2str(i)])(ismember(curexperiment.data2.l2.(['condition' num2str(i)]),[131,132,133,134,135])) = 13;
        curexperiment.data2.l2.(['condition' num2str(i)]) = unique(curexperiment.data2.l2.(['condition' num2str(i)]));
        clear cnd
    end
end
clear Mlabel M M2 M3 M4
curexperiment.data2l3_name{1}           = 'All'; % all responses (except no resp ON)
curexperiment.data2.l3.condition1       = str2double(strrep([join(string(combvec([1,2,3],[3,4,5])'));join(string(combvec([1,2,3],[1,2],[1,2,3,4,5,9])'))],' ',''));                  
% restEEG
% x. : (timepoint) 1=pre stimulation, 2=post stimulation
% .x : (eyes) 7=open, 8=closed
curexperiment.data3l1_name{1}           = 'AllPre'; % prestim
curexperiment.data3.l1.condition1       = str2double(strrep(join(string(combvec(1,[7,8])')),' ',''));    
curexperiment.data3l1_name{2}           = 'AllPost'; % poststim
curexperiment.data3.l1.condition2       = str2double(strrep(join(string(combvec(2,[7,8])')),' ',''));                      
curexperiment.data3l2_name{1}           = 'EyesOpenPre'; % eyes open prestim
curexperiment.data3.l2.condition1       = str2double(strrep(join(string(combvec(1,7)')),' ','')); 
curexperiment.data3l2_name{2}           = 'EyesClosedPre'; % eyes closed prestim
curexperiment.data3.l2.condition2       = str2double(strrep(join(string(combvec(1,8)')),' ','')); 
curexperiment.data3l2_name{3}           = 'EyesOpenPost'; % eyes open prestim
curexperiment.data3.l2.condition3       = str2double(strrep(join(string(combvec(2,7)')),' ','')); 
curexperiment.data3l2_name{4}           = 'EyesClosedPost'; % eyes closed prestim
curexperiment.data3.l2.condition4       = str2double(strrep(join(string(combvec(2,8)')),' ',''));                         
curexperiment.data3l3_name{1}           = 'All'; % prestim & poststim
curexperiment.data3.l3.condition1       = str2double(strrep(join(string(combvec([1,2],[7,8])')),' ',''));                
% levels of processing
curexperiment.data1.levels              = length(fieldnames(curexperiment.data1)); % encoding analysis levels
curexperiment.data1.level_name{1}       = '_SubMem';
curexperiment.data1.level_name{2}       = '_SubAllConds';
curexperiment.data1.level_name{3}       = '_All';
curexperiment.data2.levels              = length(fieldnames(curexperiment.data2)); % retrieval analysis levels
curexperiment.data2.level_name{1}       = '_RetMem';
curexperiment.data2.level_name{2}       = '_AllConds';
curexperiment.data2.level_name{3}       = '_All';
curexperiment.data3.levels              = length(fieldnames(curexperiment.data3)); % rest analysis levels
curexperiment.data3.level_name{1}       = '_PrePost';
curexperiment.data3.level_name{2}       = '_PrePostOpenClosed';
curexperiment.data3.level_name{3}       = '_All'; 
curexperiment.Nelectrodes               = [64 32]; % number of EEG electrodes in session 1 and session 2-4

%%%%%%%%%%%%
%% SET-UP %%
%%%%%%%%%%%%
% set the script folder as the current directory
cd(curexperiment.scriptdir)
% add the Fieldtripfolder to the path
addpath([curexperiment.dirroot 'Fieldtrip/Fieldtrip/fieldtrip-20200128'])
ft_defaults % sets the defaults and configure the minimal required path settings
% make the analysis location if it doesn't exist
if ~exist([curexperiment.analysis_loc filesep 'Plots' filesep 'Archive' filesep 'SingleSubs'], 'dir')
    mkdir([curexperiment.analysis_loc filesep 'Plots' filesep 'Archive' filesep 'SingleSubs']);
end
% locate the raw EEG files on the disk
df = [];
for i = 1:length(curexperiment.extension)
    df = [df; dir(fullfile(curexperiment.datafolder_input, '**', '*', curexperiment.extension{i}))];
end
files = sort({df.name});
% match the encoding and retrieval files
for i = 1:curexperiment.Nsubs
    pattern = num2str(i + 100) + "_[0-9]\.eeg";  % Match exactly one character after "_"
    matchIdx = ~cellfun('isempty', regexp(files, pattern, 'once'));  % Apply regex
    curexperiment.files(i, :) = {files{matchIdx}};  % Store matched filenames
end
clear files df
% create the EEG layout template
try
    if isfield(curexperiment,'elec')
        curexperiment.elecs          = ft_read_sens(curexperiment.elec.lay);
    end
catch
    fprintf('\nNo EEG template created\n')
end
% detect the channel neighbours
try
    cfg               = [];
    cfg.method        = 'triangulation'; %distance
    cfg.feedback      = 'no';
    curexperiment.neighbours        = ft_prepare_neighbours(cfg,curexperiment.elecs);
catch
    fprintf('\nno channel neighbours determined\n')
end
% add the external colormap
addpath('BrewerMap') 
% adjust warning messages to only be printed once per timeout period
ft_warning timeout 60
ft_warning once

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEHAVIORAL PREPROCESSING %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(contains(curexperiment.preproc,'beh'))
    % load in the data
    load(curexperiment.outputfile)
    %% PREPROCESSING
    % loop over the files
    txtfiles = sort({dir(fullfile(curexperiment.datafolder_inputbehav,'*txt*')).name});
    for t=1:length(txtfiles)
        % get the task name
        cur_task = txtfiles{t}(1:find(ismember(txtfiles{t},'_'))-1);
        % get the participant number
        ind = find(ismember(txtfiles{t},'-'));
        cur_ppn = txtfiles{t}(ind(1)+1:ind(2)-1);
        % get the session number
        cur_ses = txtfiles{t}(ind(2)+1);
        % see if the data already exists, if not get the data
        if ~contains(['session_',cur_ses], fieldnames(Data.(cur_task).(['ppn_', cur_ppn])) )
            fprintf('\n########\nLoading... %s task of participant %s during session %s\n########\n',cur_task,cur_ppn,cur_ses)
            % load in the text file
            data_txt = importdata(fullfile(curexperiment.datafolder_inputbehav, txtfiles{t}));
            data_tbl = table;
            % find the header variables in between the header start and end
            hdrs = find(contains(data_txt,'Header'));
            for i=hdrs(1)+1:hdrs(2)-1
                item    = split(data_txt{i},': ');
                item{1} = erase(item{1},'.'); % remove the . from the variable name to eliminate future table confusion
                if isempty(item{2})
                    item{2} = 'None';
                elseif all(ismember(item{2}, '0123456789'))
                    % add variable to table (as number)
                    item{2} = str2double(item{2});
                end
                % add header variable to the table
                data_tbl.(item{1}) = item{2}; 
            end
            % get the other variables in the trials
            trls = find(contains(data_txt,'LogFrame Start'));
            trle = find(contains(data_txt,'LogFrame End'));
            % loop over the trials
            for trli=1:length(trls)
                % get the trialnumber
                data_tbl.trial(trli,1)=trli;
                % loop over the variables
                for i=trls(trli)+1:trle(trli)-1
                    item    = split(data_txt{i},': ');
                    item{1} = strtrim(erase(item{1},'.')); % remove the . and trailing whitespade from the variable name to eliminate future table confusion
                    if isempty(item{2})
                        item{2} = 'None';
                    elseif all(ismember(item{2}, '-0123456789.,/ '))
                        % add variable to table (as number)
                        item{2} = str2double(item{2});
                    end
                    % add trial variable to the table
                    data_tbl.(item{1})(trli,1:length(item{2})) = item{2}; 
                    clear item
                end
            end
            Data.(cur_task).(['ppn_' cur_ppn]).(['session_' cur_ses]) = Data_temp;
            clear data_tbl data_txt trl* txtfiles
            fprintf('Done with ppn %s ses %s',cur_ppn, cur_ses)
        end
    end
    %% FURTHER PREPROCESSING
    % loop over participants
    for f=101:curexperiment.Nsubs+100
        % loop over the sessions
        for s=1:curexperiment.Nses
            fprintf('\n########\nBusy.. participant %d, session %d\n########\n',f,s)
            % DATA CREATION FOR TRIAL-BY-TRIAL ANALYSIS
            %% ENCODING (STUDY)
            cur_task = 'Kahn7Study';
            % store the currently relevant data in a temp variable
            cur_data = Data.(cur_task).(['ppn_' num2str(f)]).(['session_' num2str(s)]);
            % get participant data
            ppn.Subject    = repmat(cur_data.Subject(1),curexperiment.Ntrials_enc,1);
            ppn.Age        = repmat(cur_data.Age(1),curexperiment.Ntrials_enc,1);
            ppn.Sex        = repmat({cur_data.Sex(1)},curexperiment.Ntrials_enc,1);
            ppn.Session    = repmat(cur_data.Session(1),curexperiment.Ntrials_enc,1);
            % keep only experimental trials
            cur_data = cur_data(contains(cellstr(cur_data.Procedure),'StudyProc'),:);
            % adjust the values where needed
            cur_data.Rating(cur_data.Rating==78|cur_data.Rating==9) = NaN; % 9/78 = no response, change to NaN
            cur_data.EncodeDecisionRT(cur_data.EncodeDecisionRT==0)=NaN;
            % get the variables of interest for this dataset
            ppn.EncTrial   = [1:curexperiment.Ntrials_enc]';
            ppn.EncWord    = cellstr(cur_data.Word);
            ppn.SourceTask = cellstr(cur_data.Task);
            ppn.EncResp    = cur_data.Rating(:,1); %1 = unsuccessful, 2 = partially successful, 3 = successful, NaN = no response
            ppn.EncRT      = cur_data.EncodeDecisionRT;
            % save the encoding data
            idx1 =((f-100-1)*curexperiment.Nses*curexperiment.Ntrials_enc)+s*curexperiment.Ntrials_enc-(curexperiment.Ntrials_enc-1);
            idx2 = idx1+curexperiment.Ntrials_enc-1;
            Data.ENC(idx1:idx2,:) = struct2table(ppn);
            clear cur_data idx*
            %% RETRIEVAL (TEST)
            cur_task = 'Kahn7Test';
            % store the currently relevant data in a temp variable
            cur_data = Data.(cur_task).(['ppn_' num2str(f)]).(['session_' num2str(s)]);
            % get participant data
            ppn.Subject    = repmat(ppn.Subject(1),curexperiment.Ntrials_ret,1);
            ppn.Age        = repmat(ppn.Age(1),curexperiment.Ntrials_ret,1);
            ppn.Sex        = repmat({ppn.Sex(1)},curexperiment.Ntrials_ret,1);
            ppn.Session    = repmat(ppn.Session(1),curexperiment.Ntrials_ret,1);
            if strcmp(cur_data.Stimulation(1,:),'NoStim')
                ppn.Stimulation = repmat({'None'},curexperiment.Ntrials_ret,1);
            elseif strcmp(cur_data.Stimulation(1,:),'B') % 4 Hz
                ppn.Stimulation = repmat({'Theta'},curexperiment.Ntrials_ret,1);
            elseif strcmp(cur_data.Stimulation(1,:),'C') % 50 Hz
                ppn.Stimulation = repmat({'Gamma'},curexperiment.Ntrials_ret,1);
            elseif strcmp(cur_data.Stimulation(1,:),'D') % sham
                ppn.Stimulation = repmat({'Sham'},curexperiment.Ntrials_ret,1);
            end
            % keep only experimental trials
            cur_data = cur_data(contains(cellstr(cur_data.Procedure),'TestProc'),:);
            % DATA CORRECTIOM
            % adjust for a task coding mess up in the first 18 participants
            cur_data.TestSourceRT(contains(cellstr(cur_data.TestWordRESP(:,1)),{'{','i','l','N'}))=NaN;
            cur_data.TestSourceRESP(contains(cellstr(cur_data.TestWordRESP(:,1)),{'{','i','l','N'}))=NaN;
            % if there is a reponse before 800 ms, use that and overwrite the later response
            for i=1:size(cur_data,1)
                cur_data.tmp_preTestWordRESP(i,:) = pad(cur_data.preTestWordRESP(i,:),7);
                cur_data.tmp_TestWordRESP(i,:) = pad(cur_data.TestWordRESP(i,:),7);
            end
            cur_data.preTestWordRESP = cur_data.tmp_preTestWordRESP;
            cur_data.TestWordRESP = cur_data.tmp_TestWordRESP;
            for i=1:size(cur_data,1) 
                if ~contains(cur_data.preTestWordRESP(i,:),'None')
                    cur_data.TestWordRESP(i,:)=cur_data.preTestWordRESP(i,:);
                    cur_data.TestWordRT(i)=cur_data.preTestWordRT(i);
                elseif ~ismember(cur_data.TestWordRESP(i),'None')
                    cur_data.TestWordRT(i)=cur_data.TestWordRTTime(i)-cur_data.preTestWordOnsetTime(i); % add the 800 ms + delay to the RT
                end
                clear preWord
            end
            if f==117 && s==2 % adjust for the missing data
                cur_data.Word(end:curexperiment.Ntrials_ret)='';
                cur_data.Task(end:curexperiment.Ntrials_ret)='';
                cur_data.TestWordRT(end:curexperiment.Ntrials_ret)=NaN;
                cur_data.TestSourceRT(end:curexperiment.Ntrials_ret)=NaN;
                cur_data.Condition(end:curexperiment.Ntrials_ret)='';
            end
            % make the RT=0 to NaN
            cur_data.TestWordRT(cur_data.TestWordRT==0)=NaN;
            cur_data.TestSourceRT(cur_data.TestSourceRT==0)=NaN;
            % recode the responses to ratings
            %             item        source
            % 1 = h     = OldHC     = PleasHC
            % 2 = u     = OldLC     = PleasLC
            % 3 = space = Guess     = Guess
            % 4 = i     = NewLC     = PlaceLC
            % 5 = l     = NewHC     = PlaceHC
            cur_data.ItemRESP = NaN(size(cur_data,1),1);
            cur_data.ItemRESP(cur_data.TestWordRESP(:,1)=='h')=1;
            cur_data.ItemRESP(cur_data.TestWordRESP(:,1)=='u')=2;
            try cur_data.ItemRESP(all(cur_data.TestWordRESP(:,1:7)=='{SPACE}',2),:)=3;catch;end  
            cur_data.ItemRESP(cur_data.TestWordRESP(:,1)=='i')=4;
            cur_data.ItemRESP(cur_data.TestWordRESP(:,1)=='l')=5;
            try cur_data.ItemRESP(all(cur_data.TestWordRESP(:,1:4)=='None',2),:)=NaN;catch;end
            cur_data.SourceRESP = NaN(size(cur_data,1),1);
            cur_data.SourceRESP(cur_data.TestSourceRESP(:,1)=='h')=1;
            cur_data.SourceRESP(cur_data.TestSourceRESP(:,1)=='u')=2;
            try cur_data.SourceRESP(all(cur_data.TestSourceRESP(:,1:7)=='{SPACE}',2),:)=3;catch;end
            cur_data.SourceRESP(cur_data.TestSourceRESP(:,1)=='i')=4;
            cur_data.SourceRESP(cur_data.TestSourceRESP(:,1)=='l')=5;
            % get the variables of interest for this dataset
            ppn.RetTrial   = [1:curexperiment.Ntrials_ret]';
            ppn.RetWord    = cellstr(cur_data.Word);
            ppn.Source     = cellstr(cur_data.Task);
            ppn.OldNew     = cellstr(cur_data.Condition);
            ppn.ONResp     = cur_data.ItemRESP;
            ppn.ONRT       = cur_data.TestWordRT;
            ppn.SourceResp = cur_data.SourceRESP;
            ppn.SourceRT   = cur_data.TestSourceRT;
            %% COMBINE ENCODING AND RETRIEVAL TRIALS
            % copy the encoding variables we want to add
            ppn.tmpEncTrial = ppn.EncTrial;
            ppn.tmpEncResp  = ppn.EncResp;
            ppn.tmpEncRT    = ppn.EncRT;
            % create new empty versions of those variables
            ppn.EncTrial = NaN(curexperiment.Ntrials_ret,1);
            ppn.EncResp  = NaN(curexperiment.Ntrials_ret,1);
            ppn.EncRT    = NaN(curexperiment.Ntrials_ret,1);
            % loop over the encoding trials to match them to the retrieval trials
            for t=1:curexperiment.Ntrials_enc
                ppn.EncTrial(matches(ppn.RetWord,ppn.EncWord{t})) = ppn.tmpEncTrial(t);
                ppn.EncResp(matches(ppn.RetWord,ppn.EncWord{t})) = ppn.tmpEncResp(t);
                ppn.EncRT(matches(ppn.RetWord,ppn.EncWord{t})) = ppn.tmpEncRT(t);
            end
            % remove the fields we don't need anymore
            ppn.Word = ppn.RetWord;
            ppn = rmfield(ppn,{'EncWord','RetWord','SourceTask','tmpEncTrial','tmpEncResp','tmpEncRT'});
            % reorder the fieldnames
            ppn = orderfields(ppn,{'Subject','Age','Sex','Session','Stimulation','Word','RetTrial','EncTrial','OldNew','Source','EncResp','EncRT','ONResp','ONRT','SourceResp','SourceRT'});
            % save the data
            idx1 =((f-100-1)*curexperiment.Nses*curexperiment.Ntrials_ret)+s*curexperiment.Ntrials_ret-(curexperiment.Ntrials_ret-1);
            idx2 = idx1+curexperiment.Ntrials_ret-1;
            try
                Data.RET(idx1:idx2,:) = struct2table(ppn);
            catch
                Data = rmfield(Data,'RET');
                Data.RET(idx1:idx2,:) = struct2table(ppn);
            end
            clear cur_data idx* ppn
        end
    end
    % save the data
    save(curexperiment.outputfile,'Data')
    outputfile2 = strcat(curexperiment.outputfile(1:end-4),'_',string(datetime('today')),'.mat');
    save(outputfile2,'Data')
end

%%%%%%%%%%%%%%%%%%%%%%%
%% EEG PREPROCESSING %%
%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n###################\n')
fprintf('## PREPROCESSING ##\n')
fprintf('###################\n')
if any(contains(curexperiment.preproc, 'eeg'))
    % loop over subjects
    for f=101:curexperiment.Nsubs+100
        % loop over sessions
        for s=1%:curexperiment.Nses
            %% SUBJECT DATA
            subjectdata.nr        = num2str(f); % subject number
            subjectdata.ses       = num2str(s); % subject session
            subjectdata.dir_main  = fullfile(curexperiment.datafolder_output, sprintf('%s',curexperiment.name,'_', subjectdata.nr)); % subject folder
            subjectdata.dir_ses   = fullfile(curexperiment.datafolder_output, sprintf('%s',curexperiment.name,'_', subjectdata.nr),sprintf('Session_%s',subjectdata.ses)); % subject folder
            subjectdata.data      = curexperiment.files(contains(curexperiment.files,[num2str(f),'_',num2str(s)]));
            % if the subject folder does not exist make it
            if ~exist(subjectdata.dir_ses, 'dir')
              mkdir(subjectdata.dir_ses);
            end
            % check if preprocessing is already done for this participant
            if sum(contains({dir(subjectdata.dir_main).name},'PreProcessed.mat')) ~= 999 %length(curexperiment.dataset_name)
                fprintf('\n\nSUBJECT: %s\nSESSION: %s\n',subjectdata.nr, subjectdata.ses)
                %% LOAD DATA
                % check if we need to load the raw data or if that is already done for this participant
                if sum(contains({dir(subjectdata.dir_ses).name},'RawDataMarkers.mat')) ~= 1
                    fprintf('\nREADING DATA\n\n');
                    cfg             = [];
                    % don't load in the channels we don't need
                    cfg.channel     = {'all', '-EXG1-1','-EXG2-1','-EXG3-1','-EXG4-1','-EXG5-1','-EXG6-1','-EXG6-1','-EXG7-1','-EXG8-1',...
                                       '-66','-67','-68','-69','-70','-71','-72','-tACS'};
                    % loop over the datasets of this participant
                    for d=1:length(subjectdata.data)
                        % some datasets are too big to fit into memory so we preprocess one channel at a time, downsample and past together.
                        % https://www.fieldtriptoolbox.org/faq/how_can_i_preprocess_a_dataset_that_is_too_large_to_fit_into_memory/
                        % loop over all channels that need to be included in the current dataset
                        if strcmp(subjectdata.ses,'1')
                            chanN = 1; %64 channels (session 1)
                        else
                            chanN = 2; %32 channels (session 2,3,4)
                        end
                        if contains(subjectdata.data{d},'.eeg')
                            chn_strt = 1;
                            % load in the data channel by channel, skipping channels have already been processed
                            srtt0=tic;
                            for i=chn_strt:length(curexperiment.newchans{chanN})
                                srtt=tic;
                                fprintf('channel %d of %d (%s)\n',i,curexperiment.Nelectrodes(chanN),curexperiment.dataset_name{d}(2:end))
                                cfgp            = [];
                                cfgp.channel    = curexperiment.newchans{chanN}(i);
                                cfgp.dataset    = [curexperiment.datafolder_input,filesep,subjectdata.data{d}];
                                datp            = ft_preprocessing(cfgp);
                                cfgr            = [];
                                cfgr.resamplefs = curexperiment.fs_ds;
                                datr{i}         = ft_resampledata(cfgr, datp);
                                clear datp cfgp cfgr
                                endt=toc(srtt);
                                fprintf('Elapsed time for this channel is %d minutes and %d seconds\n', floor(endt/60), round(rem(endt,60)));
                                endt0=toc(srtt0);
                                fprintf('Elapsed time for this dataset is %d hour(s), %d minute(s) and %d second(s)\n', floor(endt0/3600), floor(endt0/60), round(rem(endt0,60)));
                            end
                            clear srtt* endt* chanN chn_strt
                            cfg = [];
                            data_org(d) = ft_appenddata(cfg, datr{:}); 
                            clear datr
                            save([subjectdata.dir_ses filesep subjectdata.nr '_RawData.mat'],'data_org');
                        elseif contains(subjectdata.data{d},'.set')
                            % load in the data that was merged earlier (see SouConTACS_CombineFiles.m) and resample
                            % https://www.fieldtriptoolbox.org/getting_started/eeglab/
                            fprintf('Loading..\n')
                            cd(strcat(curexperiment.dirroot, 'EEGLAB/eeglab2024_0'))
                            eeglab
                            EEG = pop_fileio(subjectdata.data{d}, 'dataformat','auto');
                            FT_EEG = eeglab2fieldtrip(EEG,'raw','none');
                            cd(strcat(curexperiment.dirroot, 'Fieldtrip'))
                            warning('off', 'MATLAB:rmpath:DirNotFound');
                            close all
                            rmpath(genpath(strcat(curexperiment.dirroot, 'EEGLAB/eeglab2024_0')))
                            cfgp            = [];
                            cfgp.channel    = curexperiment.newchans{chanN};
                            FT_EEG          = ft_preprocessing(cfgp,FT_EEG);
                            FT_EEG          = rmfield(FT_EEG,'sampleinfo');
                            data_org(d)     = FT_EEG;
                        else
                            error('What is the EEG data extension?')
                        end
                    end
                    if strcmp(curexperiment.name,'SouConTACS') && strcmp(subjectdata.ses,'1') && (strcmp(subjectdata.nr,'101') || strcmp(subjectdata.nr,'102'))
                        % switch from actiChamp layout to ActiCapSnap layout
                        for i=1:length(data_org)
                            data_org(i).label = {'Fp1';'F3';'F7';'FT9';'FC5';'FC1';'C3';'T7';'TP9';'CP5';'CP1';'Pz';'P3';'P7';'O1';'Oz';'O2';'P4';'P8';'TP10';'CP6';'CP2';'Cz';...
                                'C4';'T8';'FT10';'FC6';'FC2';'F4';'F8';'Fp2';'AF7';'AF3';'AFz';'F1';'F5';'FT7';'FC3';'C1';'C5';'TP7';'CP3';'P1';'P5';'PO7';'PO3';'POz';...
                                'PO4';'PO8';'P6';'P2';'CPz';'CP4';'TP8';'C6';'C2';'FC4';'FT8';'F6';'AF8';'AF4';'F2';'FCz'};
                        end
                    end
                    % save the data
                    save([subjectdata.dir_ses filesep subjectdata.nr '_RawData.mat'],'data_org');
                    % "define" the trials to adjust the markers
                    % trl = NxM, m1 = sample indices of start of the trial, m2 = sample indices of the end of the trial, m3 = offet of the trigger with respect to the trial.
                    cfg=[];
                    cfg.trialdef.eventtype  = curexperiment.eventtype;
                    cfg.trialfun            = 'ft_trialfun_general';
                    for d=1:size(subjectdata.data,1)
                        cfg.dataset     = subjectdata.data{d};
                        data_markers(d) = ft_definetrial(cfg);
                    end
                    % resample the trialinfo
                    for d=1:size(subjectdata.data,1)
                        if contains(subjectdata.data{d},'.eeg')
                            for i=2:length(data_markers(d).event)
                                data_markers(d).event(i).sample = data_markers(d).event(i).sample / (curexperiment.fs_org/curexperiment.fs_ds);
                            end
                            data_markers(d).trl(:,1:2) = data_markers(d).trl(:,1:2) ./ (curexperiment.fs_org/curexperiment.fs_ds);
                        elseif contains(subjectdata.data{d},'.set')
                            data_markers(d).event=rmfield(data_markers(d).event,'urevent');
                        end
                    end
                    save([subjectdata.dir_ses filesep subjectdata.nr '_RawDataMarkers.mat'],'data_markers');
                end
                
                %% ADD TRIALNUMBERS
                % check if we need to add trialnumbers to the data or if that is already done for this participant
                if sum(contains({dir(subjectdata.dir_ses).name},'AlterMarkers.mat')) ~= 1
                    % check if we need to load in data
                    if exist('data_markers','var')~=1
                        load(fullfile(subjectdata.dir_ses, strcat(subjectdata.nr,'_RawDataMarkers.mat')));
                    end
                    sprintf('\nADDING TRIALNUMBERS\n');
                    cfg                = [];
                    cfg.markers        = curexperiment.markers;
                    cfg.curexperiment  = curexperiment;
                    cfg.offset         = curexperiment.marker_offset;
                    if exist(curexperiment.outputfile,'file')
                       load(curexperiment.outputfile)
                    end
                    % call trial count function
                    for i=1:size(data_markers,2)
                        cfg.data           = data_markers(i);
                        [data_markers(i),Data.trls.(['session_' subjectdata.ses]).(['ppn' subjectdata.nr]).(['data' num2str(i)])]= trlnr(cfg);
                    end
                    save(curexperiment.outputfile,'Data')
    
                    %% ALTER MARKERS
                    % run a script which adjust the EEG markers
                    fprintf('\nALTERING MARKERS\n\n');
                    % loop over the datasets
                    for d=1:length(curexperiment.dataset_name)
                        % add a column for the original EEG markers
                        tmp=cell(size(data_markers(d).event));
                        [data_markers(d).event(:).original_marker]=deal(tmp{:});
                        clear tmp
                        i=1;
                        while i < length(data_markers(d).event)+1
                           % exclude all non-marker rows
                           if ismember(cellstr(data_markers(d).event(i).type),cellstr(curexperiment.eventtype)) == 0
                              % delete non-marker row
                               data_markers(d).event(i) = [];
                           else
                               % add a column with the original markers
                                data_markers(d).event(i).original_marker = data_markers(d).event(i).value;
                               % go to the next row
                                i = i +1;         
                           end
                        end
                        clear i
                        % COUNT ALL MARKERS
                        % for readability use a temp variable
                        tmp = curexperiment.datasets_names{d};
                        curexperiment.original_markers.(tmp).cur_count = [];
                        curexperiment.original_markers.(tmp).cur_count(1)=0;
                        for i=1:length(data_markers(d).event)     
                            % count the values
                            curexperiment.original_markers.(tmp).cur_count(strcmp(data_markers(d).event(i).original_marker, [curexperiment.original_markers.(tmp).original_marker])) = ...
                                curexperiment.original_markers.(tmp).cur_count(ismember(curexperiment.original_markers.(tmp).original_marker, data_markers(d).event(i).original_marker)) +1;
                        end
                        clear i
                        count_error = {};
                        c=1;
                        % check if all markers are there
                        for i=1:size(curexperiment.original_markers.(tmp),1) % skip practice count
                            % check if there is a predefined count value
                            if ~isempty(curexperiment.original_markers.(tmp).count_without_practice{i})
                                % check whether the predefined count value matches the actual countin this participant
                                if ~logical(curexperiment.original_markers.(tmp).count_without_practice{i} == curexperiment.original_markers.(tmp).cur_count(i))
                                    fprintf(2,[sprintf('ERROR count %s', curexperiment.original_markers.(tmp).Properties.RowNames{i}) char(10)]);
                                    count_error{c} = curexperiment.original_markers.(tmp).Properties.RowNames{i};
                                    c = c+1;
                                end
                            end
                        end
                        clear c
                        if isempty(count_error)
                            fprintf('\nNo (remaining) trialcount inconsistencies\n')
                        else
                            close all
                            figure(1);
                            t = annotation('textbox');
                            t.String = sprintf('Check the errors\nClick on me when you are done\nor abort the script');
                            t.FontSize = 20;
                            t.Position = [0.25 0.35 0.5 0.35];
                            waitforbuttonpress
                            fprintf('ok\n')
                            close figure 1
                        end
                        clear count_error t
                        % MAKE THE MARKERS NUMERIC
                        % get the indices of the stimulus onsets
                        ind = zeros(curexperiment.(strcat('Ntrials_', tmp)),1);
                        cnt=1;
                        for i=1:length(data_markers(d).event) 
                            if ~isempty(data_markers(d).event(i).trlnr)
                                ind(cnt) = i;
                                cnt=cnt+1;
                            end
                        end
                        clear cnt
                        % remove the zeros in case trial(s) is/are missing
                        ind=nonzeros(ind);
                        ind(end+1)=length(data_markers(d).event)+1;
                        % clear all the char event.values
                        [data_markers(d).event.value]=deal([]);
                        % SUBSEQUENT MEMORY
                        if d==1
                            % get all the encoding words
                            w_enc = Data.Kahn7Study.(['ppn_' subjectdata.nr]).(['session_' subjectdata.ses]).Word...
                                (all(ismember(Data.Kahn7Study.(['ppn_' subjectdata.nr]).(['session_' subjectdata.ses]).Procedure(:,1:9),'StudyProc'),2),:);
                            % get all the retrieval words
                            w_ret = Data.Kahn7Test.(['ppn_' subjectdata.nr]).(['session_' subjectdata.ses]).Word; 
                            % get the indices of the encoding words in the retrieval list
                            [~,wind]=ismember(w_enc,w_ret(:,1:size(w_enc,2)),'rows');
                            % adjust the response ratings (1,2,3,4,5,9)
                            ONRating = Data.Kahn7Test.(['ppn_' subjectdata.nr]).(['session_' subjectdata.ses]).ONRating;
                            ONRating(ONRating==0)=9; % no resp: 0->9
                            SORating = Data.Kahn7Test.(['ppn_' subjectdata.nr]).(['session_' subjectdata.ses]).SORating;
                            SORating = SORating-10; % 11->1, 12->2, etc.
                            SORating(SORating==-10)=9; % no resp: 0->9
                            % loop over the encoding words to get the subsequent memory info
                            for w=1:length(w_enc)
                                r_enc(w,1) = str2double(strcat(num2str(ONRating(wind(w))),num2str(SORating(wind(w)))));
                            end
                            clear w SORating ONRating w_enc w_ret wind
                        end
                        % loop over the stimulus onset indices to find the behavioral info
                        cnt = 1;
                        for i=ind(1:end-1)'
                            % find the following markers:
                            % (S 1), S 2, S 3 
                            % and change them in numeric values
                            if d==1 % encoding
                                % ENCODING TASK + SUBSEQUENT MEMORY
                                if strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(ind(cnt+1)-1).original_marker,'R  1') % unsuccesful
                                    data_markers(d).event(i).value=2100+r_enc(cnt);
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(ind(cnt+1)-1).original_marker,'R  2') % partly succesful
                                    data_markers(d).event(i).value=2200+r_enc(cnt);
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(ind(cnt+1)-1).original_marker,'R  3') % succesful
                                    data_markers(d).event(i).value=2300+r_enc(cnt);
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && (strcmp(data_markers(d).event(ind(cnt+1)-1).original_marker,'R  9')...
                                    || strcmp(data_markers(d).event(ind(cnt+1)-2).original_marker,'R  9')) % no response
                                    data_markers(d).event(i).value=2900+r_enc(cnt);
                    
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(ind(cnt+1)-1).original_marker,'R  1') % unsuccesful
                                    data_markers(d).event(i).value=3100+r_enc(cnt);
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(ind(cnt+1)-1).original_marker,'R  2') % partly succesful
                                    data_markers(d).event(i).value=3200+r_enc(cnt);
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(ind(cnt+1)-1).original_marker,'R  3') % succesful
                                    data_markers(d).event(i).value=3300+r_enc(cnt);
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && (strcmp(data_markers(d).event(ind(cnt+1)-1).original_marker,'R  9')...
                                    || strcmp(data_markers(d).event(ind(cnt+1)-2).original_marker,'R  9')) % no response
                                    data_markers(d).event(i).value=3900+r_enc(cnt);
                                end
                                cnt=cnt+1;
                            elseif d==2 % retrieval
                                % NEW
                                if strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 11') % very sure pleasant
                                    data_markers(d).event(i).value=111;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 12') % bit sure pleasant
                                    data_markers(d).event(i).value=112;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 13') % not sure
                                    data_markers(d).event(i).value=113;  
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 14') % bit sure place
                                    data_markers(d).event(i).value=114;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 15') % very sure place
                                    data_markers(d).event(i).value=115;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R  9') % no response
                                    data_markers(d).event(i).value=119;
                    
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 11') % very sure pleasant
                                    data_markers(d).event(i).value=121;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 12') % bit sure pleasant
                                    data_markers(d).event(i).value=122;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 13') % not sure
                                    data_markers(d).event(i).value=123;  
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 14') % bit sure place
                                    data_markers(d).event(i).value=124;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 15') % very sure place
                                    data_markers(d).event(i).value=125;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R  9') % no response
                                    data_markers(d).event(i).value=129;
                    
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  3') % not sure
                                    data_markers(d).event(i).value=13;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  4') % bit sure new
                                    data_markers(d).event(i).value=14;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  1')... % new
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  5') % very sure new
                                    data_markers(d).event(i).value=15;
                                % OLD-PLEASANT
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 11') % very sure pleasant
                                    data_markers(d).event(i).value=211;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 12') % bit sure pleasant
                                    data_markers(d).event(i).value=212;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 13') % not sure
                                    data_markers(d).event(i).value=213;  
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 14') % bit sure place
                                    data_markers(d).event(i).value=214;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 15') % very sure place
                                    data_markers(d).event(i).value=215;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R  9') % no response
                                    data_markers(d).event(i).value=219;
                    
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 11') % very sure pleasant
                                    data_markers(d).event(i).value=221;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 12') % bit sure pleasant
                                    data_markers(d).event(i).value=222;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 13') % not sure
                                    data_markers(d).event(i).value=223;  
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 14') % bit sure place
                                    data_markers(d).event(i).value=224;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 15') % very sure place
                                    data_markers(d).event(i).value=225;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R  9') % no response
                                    data_markers(d).event(i).value=229;
                    
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  3') % not sure
                                    data_markers(d).event(i).value=23;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  4') % bit sure new
                                    data_markers(d).event(i).value=24;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  2')... % pleasant
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  5') % very sure new
                                    data_markers(d).event(i).value=25;
                                %OLD-PLACE
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 11') % very sure pleasant
                                    data_markers(d).event(i).value=311;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 12') % bit sure pleasant
                                    data_markers(d).event(i).value=312;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 13') % not sure
                                    data_markers(d).event(i).value=313;  
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 14') % bit sure place
                                    data_markers(d).event(i).value=314;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 15') % very sure place
                                    data_markers(d).event(i).value=315;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  1')... % very sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R  9') % very sure place
                                    data_markers(d).event(i).value=319;
                    
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 11') % very sure pleasant
                                    data_markers(d).event(i).value=321;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 12') % bit sure pleasant
                                    data_markers(d).event(i).value=322;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 13') % not sure
                                    data_markers(d).event(i).value=323;  
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 14') % bit sure place
                                    data_markers(d).event(i).value=324;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R 15') % very sure place
                                    data_markers(d).event(i).value=325;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  2')... % bit sure old
                                    && strcmp(data_markers(d).event(i+2).original_marker,'R  9') % no response
                                    data_markers(d).event(i).value=329;
                    
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  3') % not sure
                                    data_markers(d).event(i).value=33;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  4') % bit sure new
                                    data_markers(d).event(i).value=34;
                                elseif strcmp(data_markers(d).event(i).original_marker,'S  3')... % place
                                    && strcmp(data_markers(d).event(i+1).original_marker,'R  5') % very sure new
                                    data_markers(d).event(i).value=35;
                                end
                                %NO RESPONSE
                                % add no response marker
                                if isempty(data_markers(d).event(i).value)
                                    data_markers(d).event(i).value=99;
                                end
                            end
                        end
                        clear cnt r_enc
                        % numeric markers RestEEG
                        cnt=0; %count onset/offset restEEG markers
                        if d==2
                            for i=1:length(data_markers(d).event)
                                if strcmp(data_markers(d).event(i).original_marker,'R  6') % start/end restEEG
                                    data_markers(d).event(i).value=6;
                                    cnt=cnt+1;
                                elseif strcmp(data_markers(d).event(i).original_marker,'R  7') % eyes open
                                    if cnt<2
                                        data_markers(d).event(i).value=17;
                                    else
                                        data_markers(d).event(i).value=27;
                                    end
                                elseif strcmp(data_markers(d).event(i).original_marker,'R  8') % eyes closed
                                    if cnt<2
                                        data_markers(d).event(i).value=18;
                                    else
                                        data_markers(d).event(i).value=28;
                                    end
                                end
                            end
                        end
                        clear i ind
                        % WRAPPING UP
                        % remove the 'original_marker field'
                        data_markers(d).event = rmfield(data_markers(d).event,'original_marker');
                        % create an event file
                        event{d} = data_markers(d).event;
                        % save the trl count
                        if exist(curexperiment.outputfile,'file')
                           load(curexperiment.outputfile)
                        end
                        if f==1
                            cur_trlcnt = curexperiment.original_markers.(tmp)(:,[1 3]);
                            evalc(sprintf('cur_trlcnt.Properties.VariableNames(end) = {''ppn%s_ses%s''};',subjectdata.nr,subjectdata.ses));
                            Data.trlcnt_pre.(tmp) = cur_trlcnt;
                            clear cur_trlcnt
                        else
                            Data.trlcnt_pre.(tmp).(strcat('ppn',subjectdata.nr,'_ses',subjectdata.ses)) = curexperiment.original_markers.(tmp).cur_count;
                        end
                        save(curexperiment.outputfile,'Data')
                        clear tmp ind
                    end
                    % save the data with altered markers
                    save([subjectdata.dir_ses filesep subjectdata.nr '_RawData_AlterMarkers.mat'],'data_markers'); 
                    save([subjectdata.dir_ses filesep subjectdata.nr '_Raw_AlterMarkers_Events.mat'],'event');   
                end
            
                %% RE-REFERENCING
                % check if we need to rereference the data or if that is already done for this participant
                if sum(contains({dir(subjectdata.dir_ses).name},'Rereferenced.mat')) ~= 1
                    sprintf('\nRE-REFERENCING\n');
                    % check if we need to load in data
                    if exist('data_org','var')~=1
                        load(fullfile(subjectdata.dir_ses, strcat(subjectdata.nr,'_RawData.mat')));
                    end
                    % rereference the EEG data
                    cfg             = [];
                    cfg.reref       = 'yes';
                    cfg.channel     = 'EEG';
                    cfg.refchannel  = curexperiment.elec.newref; % the average of these channels is used as the new reference
                    cfg.implicitref = curexperiment.elec.impref; % the implicit (non-recorded) reference channel is added to the data representation
                    % rereference the data
                    data_reref      = ft_preprocessing(cfg,data_org);
                    % save the rereferenced data
                    fprintf('\nSaving rereferenced data\n')
                    save([subjectdata.dir_ses filesep subjectdata.nr '_Rereferenced.mat'], 'data_reref');
                end
        
                %% FILTERING
                % check if we need to filter the data or if that is already done for this participant
                if sum(contains({dir(subjectdata.dir_ses).name},'Filtered.mat')) ~= 1
                    fprintf('\nFILTERING\n\n');
                    % check if we need to load in data
                    if exist('data_reref','var')~=1
                        load(fullfile(subjectdata.dir_ses, strcat(subjectdata.nr,'_Rereferenced.mat')));
                    end
                    % https://mailman.science.ru.nl/pipermail/fieldtrip/2011-October/017300.html
                    % https://eeglab.org/others/Firfilt_FAQ.html#q-for-granger-causality-analysis-what-filter-should-be-used-11212020-updated
                    cfg             = [];
                    cfg.bsfilter    = 'yes';
                    cfg.bsfreq      = curexperiment.bs_freq; % alternative to notch filter
                    cfg.bsfilttype  = curexperiment.filtype;
                    cfg.bsfiltord   = 2; % default=4 Butterworth and filter order of 2 based on: The effect of filtering on Granger causality based multivariate causality measures 
                    cfg.hpfilter    = 'yes';
                    cfg.hpfreq      = curexperiment.hp_freq;
                    cfg.hpfilttype  = curexperiment.filtype;
                    cfg.hpfiltord   = 2; % default=4 Butterworth and filter order of 2 based on: The effect of filtering on Granger causality based multivariate causality measures
                    cfg.lpfilter    = 'yes';
                    cfg.lpfreq      = curexperiment.lp_freq;
                    cfg.lpfilttype  = curexperiment.filtype;
                    cfg.lpfiltord   = 2; % default=4 Butterworth and filter order of 2 based on: The effect of filtering on Granger causality based multivariate causality measures
                    % filter the data
                    data_filter     = ft_preprocessing(cfg, data_reref);
                    % save the filtered data
                    save([subjectdata.dir_ses filesep subjectdata.nr '_Filtered.mat'], 'data_filter');
                end
        
                %% EPOCHING
                % check if we need to epoch the data or if that is already done for this participant
                if sum(contains({dir(subjectdata.dir_ses).name},'_Trials.mat')) ~= length(curexperiment.dataset_name)
                    fprintf('DEFINING TRIALS\n')
                    % check if we need to load in data
                    if exist('event','var')~=1
                        load(fullfile(subjectdata.dir_ses, strcat(subjectdata.nr,'_Raw_AlterMarkers_Events.mat')));
                    end
                    if exist('data_filter','var')~=1
                        load(fullfile(subjectdata.dir_ses, strcat(subjectdata.nr,'_Filtered.mat')));
                    end
                    % ENCODING
                    cfg                         = [];
                    cfg.dataset                 = [curexperiment.datafolder_input filesep subjectdata.data{contains(subjectdata.data,'Study')}];
                    cfg.trialdef.eventtype      = curexperiment.eventtype;
                    cfg.trialdef.eventvalue     = unique(cell2mat(struct2cell(curexperiment.data1.l3)));
                    cfg.trialdef.prestim        = curexperiment.stimtime(1,1)/(curexperiment.fs_org/curexperiment.fs_ds); % in seconds. adjust for the downsampling
                    cfg.trialdef.poststim       = curexperiment.stimtime(1,2)/(curexperiment.fs_org/curexperiment.fs_ds); % in seconds. adjust for the downsampling
                    cfg.event                   = event{1};
                    cfg.trialfun                = 'ft_trialfun_general';
                    % define the trials
                    trials_enc                  = ft_definetrial(cfg);
                    % add the trialcount
                    for i=1:length(trials_enc.trl)
                        % find the trlnr associated with the current epoch
                        trials_enc.trl(i,5)=event{1}( [event{1}.sample] == trials_enc.trl(i,1)+curexperiment.stimtime(1,1)*curexperiment.fs_ds ).trlnr;
                    end
                    % use the trl info to make the epochs
                    cfg = [];
                    cfg.trl = round(trials_enc.trl);
                    data_enc = ft_redefinetrial(cfg,data_filter(1));
                    data_enc.cfg.previous = []; % clear the history
                    % save trlinfo
                    save([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{1} '_Trials.mat'], curexperiment.datasets_names{1},'-v7.3'); % first dataset
                    % RETRIEVAL
                    event{2} = event{2}(~cellfun(@isempty,{event{2}.value}));
                    cnt=1;
                    for i=1:length(event{2})
                        % find the stimulus markers and select the samples between the markers
                        if ismember(event{2}(i).value,unique(cell2mat(struct2cell(curexperiment.data2.l3))))
                            trials_ret.trl(cnt,1)=event{2}(i).sample-(.500*curexperiment.fs_ds); %start epoch is -500 ms (fixations is 650-850 ms)
                            trials_ret.trl(cnt,2)=event{2}(i+1).sample-(.505*curexperiment.fs_ds); % end epoch 505 ms before the next stimulus onset (so minimally 5 ms gap between trials)
                            if (trials_ret.trl(cnt,2)-event{2}(i).sample) > (2*curexperiment.fs_ds) % if the trial is longer than 2 seconds, alter the length to 2 sec
                                trials_ret.trl(cnt,2)=event{2}(i).sample+2*curexperiment.fs_ds;
                            end
                            trials_ret.trl(cnt,3)=-.500*curexperiment.fs_ds; % epoch offset
                            trials_ret.trl(cnt,4)=event{2}(i).value; % marker
                            trials_ret.trl(cnt,5)=event{2}(i).trlnr; % marker
                            cnt=cnt+1;
                        end
                    end
                    clear cnt
                    % use the trl info to make the epochs
                    cfg = [];
                    cfg.trl = round(trials_ret.trl);
                    data_ret = ft_redefinetrial(cfg,data_filter(2));
                    data_ret.cfg.previous = []; % clear the history
                    % save trlinfo
                    save([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{2} '_Trials.mat'], curexperiment.datasets_names{2},'-v7.3'); % second dataset
                    % REST
                    cnt=1;
                    for i=1:length(event{2})
                        % find the rest markers and select the samples between the markers
                        if ismember(event{2}(i).value,unique(cell2mat(struct2cell(curexperiment.data3.l2))))
                            trials_rest.trl(cnt,1)=event{2}(i).sample; % start epoch
                            trials_rest.trl(cnt,2)=event{2}(i+1).sample-1; % end epoch
                            trials_rest.trl(cnt,3)=0; % epoch offset
                            trials_rest.trl(cnt,4)=event{2}(i).value; % marker
                            cnt=cnt+1;
                        end
                    end
                    cfg = [];
                    cfg.trl = round(trials_rest.trl);
                    data_rest = ft_redefinetrial(cfg,data_filter(2));
                    % cut the long epochs in 2 second epochs
                    cfg = [];
                    cfg.length = 2; % in seconds. adjust for the downsampling
                    data_rest = ft_redefinetrial(cfg,data_rest);
                    data_rest.cfg.previous = []; % clear the history
                    % save the data
                    save([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{3} '_Trials.mat'], curexperiment.datasets_names{3},'-v7.3'); % third dataset
                end
        
                %% TRIAL COUNT
                trl_cnt = false;
                if trl_cnt
                    fprintf('\nTRIAL COUNT PRE')
                    % load the data
                    load(curexperiment.outputfile)
                    % loop through the datasets
                    for d=1:length(curexperiment.datasets_names)
                        % load in the data
                        load([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_Trials.mat']);
                        % store the data of the current dataset
                        cur_data = eval(curexperiment.datasets_names{d});
                        % get all possible markers
                        markers = curexperiment.(['data',num2str(d)]).l3.condition1; % level: all
                        % count the markers
                        try
                            for i=1:length(markers)
                              count(i) = sum(cur_data.trialinfo(:,1)==markers(i));
                            end
                        catch
                            fprintf('no markers')
                            count(i) = 0;
                        end
                        fldnms = fieldnames(curexperiment.(['data',num2str(d)]).l2); % level: Conf/Open-Closed
                        % link count to conditions
                        for i=1:length(fldnms)
                            conmarkers = curexperiment.(['data',num2str(d)]).l2.(fldnms{i});
                            try
                                concount(i) = sum(count(ismember(markers,conmarkers)));
                            catch
                                fprintf('no markers in this condition')
                                concount(i)=0;
                            end
                        end
                        % wrapping up
                        connames = char(curexperiment.(['data',num2str(d),'l2_name']));
                        count_tab = array2table(concount,'VariableNames',cellstr(connames));
                        ppn_tab = table(str2double(subjectdata.nr),'VariableNames',{'Participant'});
                        cur_tab = [ppn_tab,count_tab];
                        Data.TrlCount.PrePreProc.(curexperiment.dataset_name{d}(2:end))(str2double(subjectdata.nr)-100,:) = cur_tab;
                        % clear the variables
                        clear markers count fldnms conmarkers connames concount *_tab cur_data
                    end
                    % save the data
                    save(curexperiment.outputfile,'Data')
                end
        
                %% ARTIFACT REJECTION (ROUND 1)
                % check if we need to do non-ocular artifact rejection or if that is already done for this participant
                if sum(contains({dir(subjectdata.dir_ses).name},'ArtiRemoved.mat')) ~= length(curexperiment.dataset_name)
                    % check if we need to load in data
                    if exist('event','var')~=1
                        load(fullfile(subjectdata.dir_ses, strcat(subjectdata.nr,'_Raw_AlterMarkers_Events.mat')));
                    end
                    % copy the events from dataset 2 (ret) to dataset 3 (rest)
                    event{3} = event{2};
                    % loop over the datasets
                    for d=1:length(curexperiment.datasets_names)
                        % check if we need to load in data
                        if exist(curexperiment.datasets_names{d},'var')~=1
                            load([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_Trials.mat']);
                        end
                        % store the data in a temporary variable
                        cur_data = eval(curexperiment.datasets_names{d});
                        % get the events of stimulus onset to plot them
                        event_arti = [];
                        for i=1:length(event{d})
                            if ismember(event{d}(i).value, curexperiment.(['data', num2str(d)]).l3.condition1)
                                event_arti = [event_arti event{d}(i)];
                            end
                        end
                        % check if there is already artifact data
                        cfg = [];
                        if exist([subjectdata.dir_ses filesep subjectdata.nr, curexperiment.dataset_name{d}, '_Artifacts.mat'],'file')==2
                            load([subjectdata.dir_ses filesep subjectdata.nr, curexperiment.dataset_name{d}, '_Artifacts.mat']);
                            cfg.artfctdef.visual.artifact = artifacts;
                        end
                        % if there is no artifact data yet select the artifacts
                        if (isfield(cfg,'artfctdef') == 0 || isempty(cfg.artfctdef.visual.artifact))... # no prior artifact data
                                && (contains(Data.Kahn7Test.(['ppn_',subjectdata.nr]).(['session_',subjectdata.ses]).Stimulation(1),{'N','D'}))... # stimulation conditions either none (N) or sham (D)
                                && curexperiment.do_art
                            cfg.viewmode        = 'vertical';
                            cfg.ploteventlabels = 'colorvalue';
                            cfg.renderer        = 'zbuffer'; %'opengl' 'painters';
                            cfg.continuous      = 'yes';
                            cfg.blocksize       = 15;
                            cfg.event           = event_arti;
                            cfg.channel         = {'all', '-tACS'};
                            cfg.colorgroups     = 'allblack';
                            cfg.ylim            = [-80 80];
                            cfg.fontsize        = 5;
                            fprintf('\nSUBJECT: %s\n', subjectdata.nr);
                            fprintf('SESSION: %s\n', subjectdata.ses);
                            fprintf('DATASET: %s\n', curexperiment.dataset_name{d}(2:end));
                            % show the data to select the artifacts
                            cfg                  = ft_databrowser(cfg,cur_data);
                            artifacts            = cfg.artfctdef.visual.artifact;
                        end
                        % load in the outputfile
                        load(curexperiment.outputfile)
                        % save the artifacts
                        Data.AR_badseg.(['ppn',subjectdata.nr]).(['ses',subjectdata.ses]).(curexperiment.dataset_name{d}(2:end)).(['date_',num2str(year(datetime)),num2str(month(datetime)),num2str(day(datetime))]) = artifacts;
                        save(curexperiment.outputfile,'Data')
                        % reject the artifacts
                        cfg.artfctdef.reject        = 'complete';
                        try
                            cleandata               = ft_rejectartifact(cfg,cur_data);
                        catch
                            % if there is an error, there is probably no data left after
                            % artifact rejection
                            cleandata            = cur_data;
                            cleandata.trial      = {};
                            cleandata.time       = {};
                            cleandata.trialinfo  = [];
                            cleandata.sampleinfo = [];
                        end
                        % save the artifact data
                        evalc(sprintf('%s = cleandata', curexperiment.datasets_names{d}));
                        save([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_ArtiRemoved.mat'], curexperiment.datasets_names{d});
                        save([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_Artifacts.mat'],'artifacts');
                        save([curexperiment.datafolder_input filesep 'Artifacts' filesep subjectdata.nr '_' subjectdata.ses curexperiment.dataset_name{d} '_Artifacts_' ['date_',num2str(year(datetime)),num2str(month(datetime)),num2str(day(datetime))] '.mat'],'artifacts');
                        clear artifacts cleandata cur_data
                    end
                end
        
                %% ICA (OCULAR ARTIFACT CORRECTION)
                % check if we need to do non-ocular artifact rejection or if that is already done for this participant
                if sum(contains({dir(subjectdata.dir_ses).name},'PostICA.mat')) ~= length(curexperiment.dataset_name)
                    if curexperiment.do_ica
                        % check if we need to load in data
                        if exist('event','var')~=1
                            load(fullfile(subjectdata.dir_ses, strcat(subjectdata.nr,'_Raw_AlterMarkers_Events.mat')));
                        end
                        % loop over the datasets
                        for d=1:length(curexperiment.datasets_names)
                            % check if we need to load in data
                            if exist(curexperiment.datasets_names{d},'var')~=1
                                load([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_ArtiRemoved.mat']);
                            end
                            % store the data in a temporary variable
                            cur_data = eval(curexperiment.datasets_names{d});
                            if ~isempty(cur_data.trial)
                                % EEGLAB ICA
                                % In order to use ICLabel to automatically select ICA components,
                                % we have to convert our Fieldtrip data to EEGLAB data.
                                % https://eeglab.org/others/EEGLAB_and_Fieldtrip.html
                                % adjust the data
                                tempdata = cur_data;
                                tempdata.hdr.nTrials=length(cur_data.trial); % add number of trials
                                tempdata.hdr.nSamplesPre = -cur_data.time{1,1}(1)*cur_data.hdr.Fs; % adjust the samples to account for the time window before stimulus onset
                                tempdata.hdr.label=cur_data.hdr.label';
                                if d==2 % make the trials identical in size by making them equal to the largest trial (adding NaNs and making them longer did not work, making them smaller made all trials smaller after component rejection)
                                    for i=1:length(cur_data.trial)
                                        trlsz(i) = size(cur_data.trial{i},2);
                                    end
                                    trlszmi = min(trlsz);
                                    trlszma = max(trlsz);
                                    for i=1:length(cur_data.trial)
                                        %tempdata.trial{i}(:,trlszmi+1:end)=[];
                                        tempdata.trial{i}(:,end+1:trlszma)=0;
                                    end
                                end
                                tempdata.trial = cat(3,tempdata.trial{:}); % reshape the data structure (chan x time x trial)
                                % adjust the events
                                if d==1 || d==2
                                    e=d;
                                    event{e} = event{e}(~cellfun(@isempty,{event{e}.trlnr}));
                                    event{e} = event{e}(cur_data.trialinfo(:,2));
                                elseif d==3
                                    e=2;
                                    event{e} = [];
                                end
                                % convert the Fieldtrip data to EEGLAB data
                                addpath(fullfile(curexperiment.dirroot,'EEGLAB','eeglab2024_0'))
                                eeglab nogui
                                EEG = pop_fileio(tempdata.hdr,tempdata.trial,event{e});
                                EEG.setname = [subjectdata.nr '_' subjectdata.ses curexperiment.dataset_name{d} '_preICA']; 
                                clear tempdata
                                % add channel location info
                                EEG=pop_chanedit(EEG, 'lookup',[curexperiment.dirroot 'EEGLAB/eeglab2024_0/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp']);
                                % save ICA EEGLAB data for SIFT
                                save([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_predataICA_EEGLAB.mat'], 'EEG');
                                % 1.5 Hz FILTER (to improve ICA, temporary filter)
                                % https://sccn.ucsd.edu/wiki/Makoto%27s_preprocessing_pipeline#Adjust_data_rank_for_ICA_.2805.2F17.2F2019_updated.29
                                EEG_org = EEG;
                                EEG = pop_eegfiltnew(EEG, 'locutoff',1.5);
                                % PERFORM ICA
                                % check the rank of the data
                                % https://sccn.ucsd.edu/pipermail/eeglablist/2018/013802.html
                                % https://sccn.ucsd.edu/pipermail/eeglablist/2018/013885.html
                                % https://sccn.ucsd.edu/pipermail/eeglablist/2010/003339.html
                                % https://eeglab.org/tutorials/06_RejectArtifacts/RunICA.html#which-ica-algorithm
                                % https://sccn.ucsd.edu/pipermail/eeglablist/2016/011049.html
                                % load in the outputfile
                                load(curexperiment.outputfile)
                                Data.ICA_rank.(['ses',subjectdata.ses]).(curexperiment.dataset_name{d}(2:end))(str2double(subjectdata.nr(2:end)),1) = rank(double(EEG.data([1:EEG.nbchan],:)'));Data.ICA_rank.(['ses',subjectdata.ses]).(curexperiment.dataset_name{d}(2:end))(str2double(subjectdata.nr(2:end)),1) = sum(eig(cov(double(EEG.data([1:EEG.nbchan],:)'))) > 1E-7);
                                if Data.ICA_rank.(['ses',subjectdata.ses]).(curexperiment.dataset_name{d}(2:end))(str2double(subjectdata.nr(2:end)),1)>63
                                    Ncomp = 63; % reducing the number of components by 1, due to re-referencing.
                                else
                                    Ncomp = Data.ICA_rank.(['ses',subjectdata.ses]).(curexperiment.dataset_name{d}(2:end))(str2double(subjectdata.nr(2:end)),1); % reducing the number of components based on the rank of the data.
                                end
                                save(curexperiment.outputfile,'Data')
                                % ICA decomposition
                                EEG = pop_runica(EEG, 'icatype', 'runica','interupt','off','pca', Ncomp, 'stop', 1e-7); % lower criterion to return cleaner compositions for high-density data
                                EEG = eeg_checkset(EEG,'ica');
                                clear cur_rank
                                % make two copied of the EEG data, one with the max data and
                                % one with the min data. This is only needed if d==2. if d~=2,
                                % EEG = EEG_long = EEG_short.
                                % copy the ICA details back to the "unfiltered" data
                                EEG_org.icaweights = EEG.icaweights;
                                EEG_org.icasphere = EEG.icasphere;
                                EEG_org.icawinv = EEG.icawinv;
                                EEG_org.icachansind = EEG.icachansind;
                                EEG_long  = EEG_org; % to be transferred back to FT data
                                EEG_short = EEG; % to be used for IClabel, which does not like zeros
                                if d==2
                                    EEG_short.data(:,trlszmi+1:end,:)=[];
                                end
                                % IClabel
                                EEG_short = pop_iclabel(EEG_short, 'default');
                                EEG_short = eeg_checkset(EEG_short);
                                % plot the ICA components
                                figure;pop_viewprops(EEG_short, 0, [1:35], {}, {}, 1, 'ICLabel' );
                                set(findall(gcf,'-property','FontSize'),'FontSize',5)
                                pause(0.05)
                                exportgraphics(gcf,[subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_ICs.png'],'ContentType','image');
                                % select components based on the IClabels
                                % https://sccn.ucsd.edu/wiki/Makoto%27s_useful_EEGLAB_code#How_to_perform_component_selection_using_ICLabel.28.29_and_dipole_information_.2802.2F20.2F2020_updated.29
                                % https://www.researchgate.net/post/EEG_IC_component_rejection
                                c.muscleIdx = find(EEG_short.etc.ic_classification.ICLabel.classifications(:,2) >= 0.9); % > 90% muscle artifacts
                                c.eyeIdx = find(EEG_short.etc.ic_classification.ICLabel.classifications(:,3) >= 0.8);    % > 80% eye artifacts
                                c.heartIdx = find(EEG_short.etc.ic_classification.ICLabel.classifications(:,4) >= 0.8);  % > 80% heart artifacts
                                c.chanIdx = find(EEG_short.etc.ic_classification.ICLabel.classifications(:,6) >= 0.9); % > 90% channel noise artifacts
                                % number of components to check for additional artifacts
                                if strcmp(subjectdata.nr,'108')||...
                                        strcmp(subjectdata.nr,'131')||...
                                        strcmp(subjectdata.nr,'145')
                                    artN = 15; % check additional components for these participants
                                elseif strcmp(subjectdata.nr,'144')
                                    artN = 20; % check additional components for this participant
                                else
                                    artN = 10;
                                end
                                [~,I] = max(EEG_short.etc.ic_classification.ICLabel.classifications(1:artN,:),[],2);   
                                I = find((I==2 | I==3 | I==4 | I==6 | I==7)); % for the first artN components, check if these are mainly channel noise, muscle, eye, heart, or other artifacts
                                c.etcIdx = I(~(ismember(I,[c.muscleIdx;c.eyeIdx;c.heartIdx;c.chanIdx])))
                                components = [c.muscleIdx; c.eyeIdx; c.heartIdx; c.chanIdx; c.etcIdx];
                                % save the excluded components
                                % load in the outputfile
                                load(curexperiment.outputfile)
                                Data.ICA_totcomp.(['ppn',subjectdata.nr]).(['ses',subjectdata.ses]).(curexperiment.dataset_name{d}(2:end)) = size(EEG.icaweights,1);
                                Data.ICA_badcomp.(['ppn',subjectdata.nr]).(['ses',subjectdata.ses]).(curexperiment.dataset_name{d}(2:end)) = c;
                                save(curexperiment.outputfile,'Data')
                                % transfer back to fieldtrip data
                                ic_data = eeglab2fieldtrip(EEG_long,'comp','none');
                                ic_data = rmfield(ic_data,'elec');
                                ic_data.fsample = cur_data.fsample;
                                ic_data.trialinfo = cur_data.trialinfo;
                                ic_data.sampleinfo = cur_data.sampleinfo;
                                if d==2
                                    % remove the zero padding from the data
                                    for i=1:length(cur_data.trial)
                                        ic_data.trial{i}(:,trlsz(i)+1:end)=[];
                                        ic_data.time{i}(:,trlsz(i)+1:end)=[];
                                    end
                                    clear trl*
                                end
                                save([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_dataICA.mat'], 'ic_data');
                                % remove the components
                                cfg = [];
                                cfg.component  = components;
                                data_iccleaned = ft_rejectcomponent(cfg,ic_data);
                                data_iccleaned.label = data_iccleaned.label';
%                                 % EOG CONTROL ANALYSES
%                                 if d==1 || d ==2
%                                     % load in EEGLAB as otherwise epoch2continuous will not run
%                                     eeglab;close
%                                     % make just eye component data to compare blinks between conditions
%                                     EEG = pop_subcomp(EEG_long, c.eyeIdx, 0, 1);
%                                     save([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_ICA_EOG.mat'], 'EEG');
%                                     % loop over the conditions
%                                     fldnms = fieldnames(curexperiment.(['data' num2str(d)]).l2);
%                                     % make an empty table
%                                     cnds = cellfun(@(x) x(2:end), curexperiment.(['data' num2str(d) 'l2_name']), 'UniformOutput', 0);
%                                     tmp = table('RowNames',cnds);
%                                     tmp.Trials          = double.empty(length(cnds),0); % The number of trials in this condition
%                                     tmp.Channel         = cell.empty(length(cnds),0); % The channel with the best signal
%                                     tmp.BlinksTotal     = double.empty(length(cnds),0); % The number of potential blinks detected in this signal
%                                     tmp.BlinksTotalGood = double.empty(length(cnds),0); % Total number of good blinks in the data set.
%                                     tmp.BlinksMin       = double.empty(length(cnds),0); % The total number of blinks per minute
%                                     tmp.BlinksMinGood   = double.empty(length(cnds),0); % The total number of good blinks per minute
%                                     for cond = 1:length(fldnms)
%                                         try
%                                             % select only the epochs from that condition
%                                             EEG2 = pop_selectevent( EEG, 'type',curexperiment.(['data' num2str(d)]).l2.([fldnms{cond}]) ,'deleteevents','off','deleteepochs','on','invertepochs','off');
%                                             tmp.Trials(cond) = EEG2.trials;
%                                             % make the data continous again for blink detection
%                                             EEG2 = epoch2continuous(EEG2);
%                                             % check the blinks per condition
%                                             % https://vislab.github.io/EEG-Blinks/
%                                             [EEG2, ~, blinks, ~, ~, blinkStatistics, ~] = pop_blinker(EEG2,struct());
%                                             close
%                                             % store the information in the table
%                                             tmp.Channel{cond} = blinkStatistics.usedLabel;
%                                             tmp.BlinksTotal(cond) = blinkStatistics.numberBlinks;
%                                             tmp.BlinksTotalGood(cond) = blinkStatistics.numberGoodBlinks;
%                                             tmp.BlinksMin(cond) = blinkStatistics.blinksPerMin(1);
%                                             tmp.BlinksMinGood(cond) = blinkStatistics.blinksPerMin(5);
%                                         catch
%                                             % store blanks in the table
%                                             tmp.Channel{cond} = '';
%                                             tmp.BlinksTotal(cond) = NaN;
%                                             tmp.BlinksTotalGood(cond) = NaN;
%                                             tmp.BlinksMin(cond) = NaN;
%                                             tmp.BlinksMinGood(cond) = NaN;
%                                         end
%                                     end
%                                     % save the blink variables
%                                     load(curexperiment.outputfile)
%                                     Data.Blinks.(curexperiment.dataset_name{d}(2:end)).(['ppn_' subjectdata.nr]) = tmp;
%                                     save(curexperiment.outputfile,'Data')
%                                 end
                                % close eeglab
                                close all
                                rmpath(genpath(fullfile(curexperiment.dirroot,'EEGLAB','eeglab2024_0')))
                                clear e cond ALL* fldnms CURRENT* EEG eeglab* ERP eyeIdx globalvars heartIdx LASTCOM Ncomp plotset PLUGINLIST STUDY c
                            else
                                % no trials left after artifact rejection
                                components = [];
                                ic_data = [];
                                data_iccleaned = cur_data;
                            end
                            % SAVE DATA
                            data_iccleaned.hdr = cur_data.hdr;
                            data_iccleaned.trl = cur_data.cfg.trl;
                            evalc(sprintf('%s = data_iccleaned', curexperiment.datasets_names{d}));
                            save([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_Components.mat'],'components');
                            save([curexperiment.datafolder_input filesep 'ICA' filesep subjectdata.nr curexperiment.dataset_name{d} '_Components_' date '.mat'],'components');
                            save([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_dataICA.mat'], 'ic_data');
                            save([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_PostICA.mat'], curexperiment.datasets_names{d});
                            % plot the channel variance after ICA
                            ft_rejectvisual(cfg,data_iccleaned)
                            pause(.5) % pause half a second before saving the figure
                            exportgraphics(gcf,[subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_PostICATrials.png'],'ContentType','image');
                            clear components data_iccleaned ic_data respICA
                            close all
                        end
                    else
                        % SKIP ICA
                        fprintf('\nSKIP ICA\n')
                        % loop through the datasets
                        for d=1:length(curexperiment.datasets)
                            % check if we need to load in data
                            if exist(curexperiment.datasets_names{d},'var')~=1
                                load([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_ArtiRemoved.mat']);
                            end
                            components = [];
                            ic_data = [];
                            fprintf('\nSaving data\n')
                            save([subjectdata.subjectdir filesep subjectdata.nr curexperiment.dataset_name{d} '_Components.mat'],'components');
                            save([curexperiment.datafolder_input filesep 'ICA' filesep subjectdata.nr curexperiment.dataset_name{d} '_Components_' date '.mat'],'components');
                            save([subjectdata.subjectdir filesep subjectdata.nr curexperiment.dataset_name{d} '_dataICA.mat'], 'ic_data');
                            save([subjectdata.subjectdir filesep subjectdata.nr curexperiment.dataset_name{d} '_PostICA.mat'], curexperiment.datasets_names{d});
                            clear ic_data components
                        end
                    end
                    clear event
                end

                %% ARTIFACT REJECTION (ROUND 2)
                % check if we need to do non-ocular artifact rejection or if that is already done for this participant
                if sum(contains({dir(subjectdata.dir_ses).name},'ArtiRemovedFin.mat')) ~= 999 %length(curexperiment.dataset_name)
                    % load in the outputfile
                    load(curexperiment.outputfile)
                    % check if we need to load in data
                    if exist('event','var')~=1
                        load(fullfile(subjectdata.dir_ses, strcat(subjectdata.nr,'_Raw_AlterMarkers_Events.mat')));
                    end
                    % copy the events from dataset 2 (ret) to dataset 3 (rest)
                    event{3} = event{2};
                    % loop over the datasets
                    for d=1:length(curexperiment.datasets_names)
                        % check if we need to load in data
                        if exist(curexperiment.datasets_names{d},'var')~=1
                            load([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_PostICA.mat']);
                        end
                        % store the data in a temporary variable
                        cur_data = eval(curexperiment.datasets_names{d});
                        % get the events of stimulus onset to plot them
                        event_arti = [];
                        for i=1:length(event{d})
                            if ismember(event{d}(i).value, curexperiment.(['data', num2str(d)]).l3.condition1)
                                event_arti = [event_arti event{d}(i)];
                            end
                        end
                        % check if there is already artifact data
                        cfg = [];
                        if exist([subjectdata.dir_ses filesep subjectdata.nr, curexperiment.dataset_name{d}, '_ArtifactsFin.mat'],'file')==2
                            load([subjectdata.dir_ses filesep subjectdata.nr, curexperiment.dataset_name{d}, '_ArtifactsFin.mat']);
                            cfg.artfctdef.visual.artifact = artifacts;
                        end
                        % (optional) channel color groups
                        colorgrps = true;
                        if colorgrps
                            chans_highl = {'F7','AF8'};
                            cfg.colorgroups     = ismember(cur_data.label,chans_highl)+1;
                            cfg.linecolor = [.5 .5 .5;0 0 0; 1 0 0]; % has to be at least 3 due to fieldtrip code, but we only use 2
                        else
                            cfg.colorgroups     = 'allblack';
                        end
                        % if there is no artifact data yet select the artifacts
                        if (isfield(cfg,'artfctdef') == 0 || isempty(cfg.artfctdef.visual.artifact))... # no prior artifact data
                                && (contains(Data.Kahn7Test.(['ppn_',subjectdata.nr]).(['session_',subjectdata.ses]).Stimulation(1),{'N','D'}))... # stimulation conditions either none (N) or sham (D)
                                && curexperiment.do_art
                            cfg.viewmode        = 'vertical';
                            cfg.ploteventlabels = 'colorvalue';
                            cfg.renderer        = 'zbuffer'; %'opengl' 'painters';
                            %cfg.continuous      = 'yes';
                            %cfg.blocksize       = 15;
                            cfg.event           = event_arti;
                            cfg.channel         = {'all', '-tACS'};
                            cfg.ylim            = [-30 30]; %[-80 80];
                            cfg.fontsize        = 5;
                            fprintf('\nSUBJECT: %s\n', subjectdata.nr);
                            fprintf('SESSION: %s\n', subjectdata.ses);
                            fprintf('DATASET: %s\n', curexperiment.dataset_name{d}(2:end));
                            % show the data to select the artifacts
                            cfg                  = ft_databrowser(cfg,cur_data);
                            artifacts            = cfg.artfctdef.visual.artifact;
                        end
                        % save the artifacts
                        Data.AR_badseg2.(['ppn',subjectdata.nr]).(['ses',subjectdata.ses]).(curexperiment.dataset_name{d}(2:end)).(['date_',num2str(year(datetime)),num2str(month(datetime)),num2str(day(datetime))]) = artifacts;
                        save(curexperiment.outputfile,'Data')
                        % for some reason some trials seem to be duplicated, so this removes the duplicates (quick fix)
                        if any(diff(cur_data.sampleinfo)<0)
                            idx1 = find(diff(cur_data.sampleinfo)<0, 1, 'first')+1;
                            smp_min = min(cur_data.sampleinfo(idx1:end,:),[],'all');
                            smp_max = max(cur_data.sampleinfo(idx1:end,:),[],'all');
                            cfg.artfctdef.visual.artifact(end+1,:) = [smp_min smp_max];
                            clear smp*
                        end
                        % reject the artifacts
                        cfg.artfctdef.reject        = 'complete';
                        try
                            cleandata               = ft_rejectartifact(cfg,cur_data);
                        catch
                            % if there is an error, there is probably no data left after
                            % artifact rejection
                            cleandata            = cur_data;
                            cleandata.trial      = {};
                            cleandata.time       = {};
                            cleandata.trialinfo  = [];
                            cleandata.sampleinfo = [];
                        end
                        % save the artifact data
                        evalc(sprintf('%s = cleandata', curexperiment.datasets_names{d}));
                        save([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_ArtiRemovedFin.mat'], curexperiment.datasets_names{d});
                        save([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_ArtifactsFin.mat'],'artifacts');
                        save([curexperiment.datafolder_input filesep 'Artifacts' filesep subjectdata.nr '_' subjectdata.ses curexperiment.dataset_name{d} '_ArtifactsFin_' ['date_',num2str(year(datetime)),num2str(month(datetime)),num2str(day(datetime))] '.mat'],'artifacts');
                        clear artifacts cleandata cur_data
                    end
                end
        
                %% TRIAL COUNT
                trl_cnt = true;
                if trl_cnt
                    fprintf('\nTRIAL COUNT POST')
                    % load the data
                    load(curexperiment.outputfile)
                    % loop through the datasets
                    for d=1:length(curexperiment.datasets_names)
                        % load in the data
                        load([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_ArtiRemovedFin.mat']);
                        % store the data of the current dataset
                        cur_data = eval(curexperiment.datasets_names{d});
                        % get all possible markers
                        markers = curexperiment.(['data',num2str(d)]).l3.condition1; % level: all
                        % count the markers
                        try
                            for i=1:length(markers)
                              count(i) = sum(cur_data.trialinfo(:,1)==markers(i));
                            end
                        catch
                            fprintf('no markers')
                            count(i) = 0;
                        end
                        fldnms = fieldnames(curexperiment.(['data',num2str(d)]).l2); % level: Conf/Open-Closed
                        % link count to conditions
                        for i=1:length(fldnms)
                            conmarkers = curexperiment.(['data',num2str(d)]).l2.(fldnms{i});
                            try
                                concount(i) = sum(count(ismember(markers,conmarkers)));
                            catch
                                fprintf('no markers in this condition')
                                concount(i)=0;
                            end
                        end
                        % wrapping up
                        connames = char(curexperiment.(['data',num2str(d),'l2_name']));
                        count_tab = array2table(concount,'VariableNames',cellstr(connames));
                        ppn_tab = table(str2double(subjectdata.nr),'VariableNames',{'Participant'});
                        cur_tab = [ppn_tab,count_tab];
                        Data.TrlCount.PostPreProc.(curexperiment.dataset_name{d}(2:end))(str2double(subjectdata.nr)-100,:) = cur_tab;
                        % clear the variables
                        clear markers count fldnms conmarkers connames concount *_tab cur_data
                    end
                    % save the data
                    save(curexperiment.outputfile,'Data')
                end
        
                %% END OF PREPROCESSING
                if sum(contains({dir(subjectdata.dir_main).name},'PreProcessed.mat')) ~= length(curexperiment.dataset_name)
                    fprintf('\nEND OF PREPROCESSING\n');
                    % loop over datasets
                    for d=1:length(curexperiment.datasets_names)
                        % load in data
                        load([subjectdata.dir_ses filesep subjectdata.nr curexperiment.dataset_name{d} '_ArtifactsFin.mat']);
                        fprintf('\nSAVING PREPROCESSED %s\n',curexperiment.dataset_name{d}(2:end));
                        eval(sprintf('%s.cfg.previous=[];',curexperiment.datasets_names{d}));
                        % save the preprocessed data
                        save([subjectdata.dir_main filesep subjectdata.nr curexperiment.dataset_name{d} '_Session' subjectdata.ses '_PreProcessed.mat'], curexperiment.datasets_names{d});
                    end
                    % save data
                    if exist('Data','var')==1
                        save(curexperiment.outputfile,'Data')
                        save(strcat(curexperiment.outputfile(1:end-4),'_',string(datetime('today')),'.mat'),'Data')
                    end
                end
                clear *dat data_* filesm PreProc* preproc* subjectdata temp i  *cur_data event
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EEG SUBJECT-LEVEL ANALYSES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n############################\n')
fprintf('## SUBJECT-LEVEL ANALYSES ##\n')
fprintf('############################\n')
% loop over subjects
if any(contains(curexperiment.analyses,'ana'))
    for f=101:curexperiment.Nsubs+100
        % loop over sessions
        for s=1%:curexperiment.Nses
            %% SUBJECT DATA
            subjectdata.nr        = num2str(f); % subject number
            subjectdata.ses       = num2str(s); % subject session
            subjectdata.dir_main  = fullfile(curexperiment.datafolder_output, sprintf('%s',curexperiment.name,'_', subjectdata.nr)); % subject folder
            subjectdata.dir_ses   = fullfile(curexperiment.datafolder_output, sprintf('%s',curexperiment.name,'_', subjectdata.nr),sprintf('Session_%s',subjectdata.ses)); % subject folder
            subjectdata.data      = curexperiment.files(contains(curexperiment.files,[num2str(f),'_',num2str(s)]));
            fprintf('\n\nSUBJECT: %s\nSESSION: %s\n',subjectdata.nr, subjectdata.ses)
            %% LOAD THE DATA
            fprintf('\nLOADING FULLY PREPROCESSED DATA\n');
            for d=1:length(curexperiment.dataset_name)
                load(fullfile(subjectdata.dir_main, [subjectdata.nr,curexperiment.dataset_name{d},'_Session',subjectdata.ses,'_PreProcessed.mat']));
            end
            %% STORE INFO ON RET TRIAL DURATION
            % load Data
            load(curexperiment.outputfile)
            % get the ret times to get the minimun trial end, the mean trial end and the standard deviation of the trial end
            ret_times = cell2mat(cellfun(@(x) x(end), data_ret.time,'UniformOutput',false));
            % if the table does not exist yet, make it
            if ~any(contains(fieldnames(Data),'RetTrialInfo'))
                sz = [curexperiment.Nsubs 3];
                varTypes = ["double","double","double"];
                varNames = ["MinDur","MeanDur","StdDur"];
                Data.RetTrialInfo = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);
                Data.RetTrialInfo{:,:}(Data.RetTrialInfo{:,:}==0) = nan;
            end
            Data.RetTrialInfo(f-100,1:3) = {min(ret_times),mean(ret_times,'omitnan'),std(ret_times,'omitnan')};
            % save data
            save(curexperiment.outputfile,'Data')
            save(strcat(curexperiment.outputfile(1:end-4),'_',string(datetime('today')),'.mat'),'Data')

            %% EVENT-RELATED POTENTIALS
            % check if we need to run the script that does ERP analyses on ubject level
            if any(contains('erp',curexperiment.analyses)) && sum(~cellfun('isempty',regexp({dir(curexperiment.analysis_loc).name},[subjectdata.nr '.*ERP\.mat']))) < 10 % change this number to the amount of files that you expect to get out of this analysis for each participant
                fprintf('\n## ERPs ##\n')
                % load Data
                load(curexperiment.outputfile)
                % loop through the datasets
                for d=1:length(curexperiment.dataset_name)-1 % exclude Rest EEG
                    fprintf('\nCurrent dataset: %s\n', curexperiment.dataset_name{d}(2:end))
                    % store the data of the current dataset
                    cur_data = eval(curexperiment.datasets_names{d});
                    % loop over levels of processing
                    for l=3%1:curexperiment.(['data',num2str(d)]).levels      
                        %% ERP ANALYSES
                        % Select trials based upon condition
                        for t=1:length(fieldnames(curexperiment.(['data', num2str(d)]).(['l', num2str(l)])))
                            cfg = [];
                            % find the trails corresponding to the current condition
                            cfg.trials = find(ismember(cur_data.trialinfo(:,1),curexperiment.(['data', num2str(d)]).(['l', num2str(l)]).(['condition',num2str(t)])));
                            if length(cfg.trials) > 5 % only include if there are more than 5 trial(s) in the condition. I lowered this because the low number of LC responses.
                                % ERP analysis per condition
                                cfg.channel            = 'EEG';     
                                cfg.vartrllength       = 2; % accept variable trial length
                                data_ERP(t)            = ft_timelockanalysis(cfg,cur_data);
                                % perform baseline correction
                                cfg.baseline           = curexperiment.erp.basewin;
                                data_ERP_norm(t)       = ft_timelockbaseline(cfg,data_ERP(t));
                %                 % NO BASELINE 
                %                 data_ERP_norm(t) = data_ERP(t);
                                % save the data
                                data_cond              = data_ERP_norm(t);
                                data_cond.cfg.previous = [];
                                save([curexperiment.analysis_loc filesep subjectdata.nr '_' subjectdata.ses curexperiment.dataset_name{d} curexperiment.(['data', num2str(d)]).level_name{l} '_' curexperiment.(['data',num2str(d),'l',num2str(l),'_name']){t} '_ERP'],'-v7.3','data_cond'); % added the '-v7.3' because of the size of the files
                            end   
                        end  
                        clear data_ERP* curcnd curdat t data_cond
                    end
                end
                clear d i
            end

            %% TIME-FREQUENCY ANALYSES
            % check if we need to run the script that does time-frequency analyses on subject level
            if any(contains('pow',curexperiment.analyses)) && sum(~cellfun('isempty',regexp({dir(curexperiment.analysis_loc).name},[subjectdata.nr '.*Total_TF\.mat']))) < 10 % change this number to the amount of files that you expect to get out of this analysis for each participant 
                % Do within subject time-frequency analyses. If the level is "all" compute
                % the TRF on a trial-by-trial basis, otherwise, just average over conditions.
                fprintf('\n## TRFs+ ##\n')
                % load Data
                load(curexperiment.outputfile)
                % loop over the datasets
                for d=1:length(curexperiment.dataset_name)-1 % skip resting-state EEG
                    fprintf('\nCurrent dataset: %s\n', curexperiment.dataset_name{d}(2:end))
                    % store the data of the current dataset
                    cur_data = eval(curexperiment.datasets_names{d});
                    % loop over levels of processing
                    for l=1:curexperiment.(['data',num2str(d)]).levels 
                        %% TRIAL SELECTION
                        % Select trials based upon condition
                        for t=1:length(fieldnames(curexperiment.(['data', num2str(d)]).(['l', num2str(l)])))
                            cfg = [];
                            % find the trails corresponding to the current condition
                            cfg.trials = find(ismember(cur_data.trialinfo(:,1),curexperiment.(['data', num2str(d)]).(['l', num2str(l)]).(['condition',num2str(t)])));
                            if length(cfg.trials) > 5 % only include if there are more than 5 trial(s) in the condition. I lowered this because the low number of LC responses.
                                data(t) = ft_selectdata(cfg,cur_data);   
                            else
                                if exist('data','var')
                                    data(t).fsample   =[];
                                    data(t).trialinfo =[];
                                end
                            end
                        end
                        % handle empty data
                        if ~exist('data','var')
                            data.fsample   =[];
                            data.trialinfo =[];
                        end
                        %% TOTAL POWER
                        % frequency analyses per condition
                        for i=1:length(data)
                            fprintf('\nCurrent frequency analysis: %s - Condition %s\n',curexperiment.dataset_name{d}(2:end),num2str(i))
                            if not(isempty(data(i).trialinfo))
                                cfg                 = [];
                                cfg.channel         = 'EEG';
                                cfg.padtype         = 'zero';
                                cfg.pad             = 5;
                                % time-frequency analysis, output: FOURIER
                                cfg.output          = 'fourier';
                                cfg.method          = 'mtmfft';
                                cfg.keeptrials      = 'yes';
                                cfg.taper           = curexperiment.pow.taptype;
                                cfg.keeptapers      = 'yes';
                                cfg.foi             = curexperiment.pow.freq_interest; 
                                % cut the data between 300 and 800 ms
                                cfg_d               = [];
                                cfg_d.toilim        = [.3 .8];
                                data_trm(i)         = ft_redefinetrial(cfg_d, data(i));
                                data_freq_fou(i)    = ft_freqanalysis(cfg, data_trm(i)); 
                                
                                % time-frequency analysis, output: POWER
                                cfg.output          = 'pow';
                                cfg.method          = 'mtmconvol';
                                cfg.keeptapers      = 'no';
                                % if all trials are included, we will get the data per trial
                                if strcmp(curexperiment.(['data', num2str(d)]).level_name{l},'_All')
                                    cfg.keeptrials  = 'yes';
                                % otherwise we will average over trials
                                else
                                    cfg.keeptrials  = 'no';
                                end
                                % get the maximum duration of trial length
                                t_min               = min(cellfun(@min, data(i).time));
                                t_max               = max(cellfun(@max, data(i).time));
                                temp_pad            = (cfg.pad-(abs(t_min)+t_max))-abs(-1-t_min); % the difference between the rounded and actual trial time minus prepad. Now everything has the same timepoints.
                                cfg.toi             = t_min-abs(-1-t_min):.01:t_max+temp_pad; % add the temp_pad duration before the start and after the end of the trail.
                                clear temp_pad t_max t_min
                                % LOW FREQUENCIES
                                cfg.taper           = curexperiment.pow.low.taptype;
                                cfg.foi             = curexperiment.pow.low.freq_interest; 
                                cfg.t_ftimwin       = curexperiment.pow.low.timwin; 
                                data_freqLow_pow(i) = ft_freqanalysis(cfg, data(i));    
                                % HIGH FREQUENCIES
                                cfg.taper           = curexperiment.pow.high.taptype;
                                cfg.foi             = curexperiment.pow.high.freq_interest; 
                                % if all trials are included, we will get the data per trial
                                cfg.t_ftimwin       = curexperiment.pow.high.timwin; 
                                cfg.tapsmofrq       = curexperiment.pow.high.tapsmo;
                                data_freqHigh_pow(i) = ft_freqanalysis(cfg, data(i));    

                                % PASTE THE LOW AND HIGH POWER DATA TOGETHER
                                cfg                 = [];
                                cfg.appenddim       = 'freq';
                                data_freq_pow(i)    = ft_appendfreq(cfg,data_freqLow_pow(i),data_freqHigh_pow(i));

                                % GET THE AVERAGE AND PEAK POWER PER CONDITION
                                if l==2
                                    cur_freqdata = data_freq_pow(i);
                                    % perform baseline correction
                                    cfg                  = [];
                                    cfg.baseline         = curexperiment.pow.basewin;
                                    cfg.baselinetype     = 'relative'; % Expressing, for each frequency, the raw power values as the relative increase or decrease with respect to the power in the baseline interval. This means active period/baseline. Note that the relative baseline is expressed as a ratio; i.e. no change is represented by 1.
                                    cur_freqdata_bl      = ft_freqbaseline(cfg, cur_freqdata);
                                    flnms = fieldnames(curexperiment.chngrp);
                                    for c=1:length(flnms)
                                        % get the channels of interest
                                        cur_chan = ismember(cur_freqdata(find(~cellfun(@isempty,{cur_freqdata.label}),1)).label,curexperiment.chngrp.(flnms{c}));
                                        for fr=1:length(curexperiment.pow.TbTFreqs)
                                            % get the frequencies of interest
                                            cur_freqs = curexperiment.pow.TbTFreqs;
                                            cur_freq  = find(cur_freqdata.freq >= cur_freqs(fr,1) & cur_freqdata.freq <= cur_freqs(fr,2));
                                            cur_time  = cur_freqdata.time >= .3 & cur_freqdata.time <= .8;
                                            % get the (relative) peak frequency
                                            % & the average power over trials, channels and frequencies
                                            if fr == 1
                                                pow_t.(flnms{c})(i)  = mean(cur_freqdata_bl.powspctrm(cur_chan,cur_freq,cur_time),'all','omitnan');
                                                [M,I]                = max(mean(cur_freqdata_bl.powspctrm(cur_chan,cur_freq,cur_time),[1,3],'omitnan'));
                                                peak_t.(flnms{c})(i) = cur_freq(I);
                                            elseif fr == 2
                                                pow_g.(flnms{c})(i)  = mean(cur_freqdata_bl.powspctrm(cur_chan,cur_freq,cur_time),'all','omitnan');
                                                [M,I]                = max(mean(cur_freqdata_bl.powspctrm(cur_chan,cur_freq,cur_time),[1,3],'omitnan'));
                                                peak_g.(flnms{c})(i) = cur_freq(I);
                                            end
                                        end
                                        
                                    end
                                    clear cur_freqdata*
                                end
                            else
                                fprintf('Not enough trials, so no FT analyses done')
                                if exist('data_freq_pow','var')
                                    data_freq_pow(i).freq = [];
                                end
                                if exist('data_freq_fou','var')
                                    data_freq_fou(i).freq = [];
                                end
                                flnms = fieldnames(curexperiment.chngrp);
                                for c=1:length(flnms)
                                    pow_t.(flnms{c})(i)  = NaN;
                                    pow_g.(flnms{c})(i)  = NaN;
                                    peak_t.(flnms{c})(i)  = NaN;
                                    peak_g.(flnms{c})(i)  = NaN;
                                end
                            end
                        end
                        clear data tmp_data_freq cur_freqs i
                        %% DIFFERENCE PLOTS (SINGLEPLOTS & TOPOPLOTS)
                        try
                            % only make plots in certain cases
                            if d==2 && l==1 
                                curconname = curexperiment.(['data' num2str(d) 'l' num2str(l) '_name']);
                                % CALCULATE THE DIFFERENCE
                                cur_diffs = {'iHitHC-iHitLC', 'iCRHC-iCRLC', 'iHit-CR'};
                                for dif = 1:length(cur_diffs)
                                    cur_diff = contains(strcat('_',curconname,'_'),strcat('_',split(cur_diffs{dif},'-'),'_'));
                                    limt = [1.2, .1];
                                    % store only the datasets that are currently needed
                                    temp_GA = data_freq_pow(cur_diff);
                                    % compute the difference
                                    cfg = [];
                                    cfg.parameter = 'powspctrm';
                                    cfg.operation = 'subtract'; %'x1 - x2';
                                    TF_diff = ft_math(cfg,temp_GA(1),temp_GA(2));
                                    clear temp_GA
                                    % SINGLEPLOT (difference plot)
                                    % no baseline correction for the difference plots
                                    for fr=1:size(curexperiment.pow.TbTFreqs,2)
                                        cngrps = fieldnames(curexperiment.chngrp);
                                        for cn=1:length(cngrps)
                                            % plot both the raw values and the t-values
                                            cur_data = TF_diff;
                                            % compute the mean over the relevant channels
                                            pl.timewin = [-.2 1];
                                            if fr==1
                                                pl.freqwin = [3 30];
                                                cur_freq = 'theta';
                                                lim = limt(1);
                                                tcks = 2;
                                                sig_area = [.3,curexperiment.pow.TbTFreqs(1,1),.5,curexperiment.pow.TbTFreqs(1,2)-curexperiment.pow.TbTFreqs(1,1)];
                                            elseif fr==2
                                                pl.freqwin = [30 55];
                                                cur_freq = 'gamma';
                                                lim = limt(2);
                                                tcks = 4;
                                                sig_area = [.3,curexperiment.pow.TbTFreqs(2,1),.5,curexperiment.pow.TbTFreqs(2,2)-curexperiment.pow.TbTFreqs(2,1)];
                                            end
                                            pl.time = logical(cur_data.time >= pl.timewin(1) & cur_data.time <= pl.timewin(2));
                                            pl.freq = logical(cur_data.freq >= pl.freqwin(1) & cur_data.freq <= pl.freqwin(2));
                                            pl.chan = ismember(cur_data.label,curexperiment.chngrp.(cngrps{cn}));
                                            meanpow = squeeze(mean(cur_data.powspctrm(pl.chan,pl.freq,pl.time),1,'omitnan'));
                                            pl.time = cur_data.time(pl.time);
                                            pl.freq = cur_data.freq(pl.freq);
                                            % interpolate to make the plot smooth
                                            pl.time_inp = linspace(pl.timewin(1),pl.timewin(2),512);
                                            pl.freq_inp = linspace(pl.freqwin(1),pl.freqwin(2),512);
                                            [org_time,org_freq] = meshgrid(pl.time,pl.freq);
                                            [inp_time,inp_freq] = meshgrid(pl.time_inp,pl.freq_inp);
                                            meanpow = interp2(org_time,org_freq,meanpow,inp_time,inp_freq,'spline');
                                            % start the figure
                                            fig = figure();
                                            imagesc(pl.time,pl.freq,meanpow);
                                            xlim(pl.timewin);
                                            axis xy
                                            % add the labels
                                            xlabel('Time (s)')
                                            ylabel('Frequency (Hz)')
                                            % adjust the ticks
                                            yticks(pl.freq(1):tcks:pl.freq(end))
                                            % make the colorbar
                                            colormap(brewermap(256,'*RdBu'));
                                            h = colorbar();
                                            h.Location = 'southoutside';
                                            clim([-lim lim+lim*.01]);
                                            ylabel(h,'Power difference (µV)')
                                            % line at t=0
                                            hold on;
                                            xline(0,'--k','LineWidth',2);
                                            % add a rectangle to show the important area
                                            rectangle('Position',sig_area,'LineWidth',2)
                                            % change the font and plot size
                                            fontsize(fig, 25,"points")
                                            set(gcf,'Position', [1,1,1000,1000])
                                            pause(.2)
                                            %set(gcf,'Position', get(0, 'Screensize')) % make full screen
                                            exportgraphics(fig,[curexperiment.analysis_loc filesep 'Plots' filesep 'SingleSubs' filesep subjectdata.nr '_Session' subjectdata.ses curexperiment.dataset_name{d} '_TF_DiffSing' curexperiment.(['data' num2str(d)]).level_name{l} '_' cur_freq '_' cngrps{cn} '_' cur_diffs{dif} '.png']);
                                            exportgraphics(fig,[curexperiment.analysis_loc filesep 'Plots' filesep 'Archive' filesep 'SingleSubs' filesep subjectdata.nr '_Session' subjectdata.ses curexperiment.dataset_name{d} '_TF_DiffSing' curexperiment.(['data' num2str(d)]).level_name{l} '_' cur_freq '_' cngrps{cn} '_' cur_diffs{dif} '_' date '.png']);
                                        end
                                        clear pl meanpow cur_freq
                                    end
                                    close all
                                    % TOPOPLOT (difference plot)
                                    cfg               = [];
                                    cfg.xlim          = [.3 .8];
                                    cfg.colorbar      = 'yes';
                                    cfg.interactive   = 'no';
                                    cfg.layout        = curexperiment.elec.lay;
                                    cfg.renderer      = 'painters';                               
                                    cfg.marker        = 'off';
                                    cfg.highlight     = 'labels';
                                    %cfg.fontsize      = 25;
                                    for fr=1:size(curexperiment.pow.TbTFreqs,2)
                                        % plot both the raw values and the t-values
                                        cur_data = TF_diff;
                                        if fr==1
                                            cfg.ylim = curexperiment.pow.TbTFreqs(1,:); % theta range
                                            cur_freq = 'theta';
                                            cfg.zlim = [-5 5];
                                        elseif fr==2
                                            cfg.ylim = curexperiment.pow.TbTFreqs(2,:); % gamma range
                                            cur_freq = 'gamma';
                                            cfg.zlim = [-.1 .1];
                                        end
                                        figure;ft_topoplotTFR(cfg,cur_data);
                                        %colorbar('Ticks',cfg.zlim, 'FontSize',18)
                                        %uistack(r,'bottom')
                                        colormap(brewermap(256,'*RdBu'));
                                        h = colorbar();
                                        ylabel(h,'Power difference (µV)')
                                        title({curexperiment.dataset_name{d}(2:end);cur_diffs{dif};cur_freq}, 'Interpreter', 'none');
                                        %set(gca,'fontsize',20)
                                        %set(gcf, 'Position', get(0, 'Screensize')) % make full screen
                                        pause(0.05) % in seconds %to save correctly
                                        exportgraphics(gcf,[curexperiment.analysis_loc filesep 'Plots' filesep 'SingleSubs' filesep subjectdata.nr '_Session' subjectdata.ses curexperiment.dataset_name{d} '_TF_DiffTopo' curexperiment.(['data' num2str(d)]).level_name{l} '_' cur_freq '_' cur_diffs{dif} '.png']);
                                        exportgraphics(gcf,[curexperiment.analysis_loc filesep 'Plots' filesep 'Archive' filesep 'SingleSubs' filesep subjectdata.nr '_Session' subjectdata.ses curexperiment.dataset_name{d} '_TF_DiffTopo' curexperiment.(['data' num2str(d)]).level_name{l} '_' cur_freq '_' cur_diffs{dif} '_' date '.png']);
                                        clear cur_freq
                                    end
                                    clear TF_diff cur_diff
                                    close all
                                end
                            end
                        catch
                            fprintf('\nCould not make difference plot\n\n')
                            end
                        %% PHASE-AMPLITUDE COUPLING
                        % We have to keep in mind that because every Fourier 
                        % transform only gives complex number, every trial only
                        % has one phase (per frequency) and not all 20 bins. 
                        if l==2 || l==3
                            flnms = fieldnames(curexperiment.chngrp);
                            if exist('data_freq_fou','var')
                                cfg             = [];
                                cfg.freqlow     = curexperiment.pow.TbTFreqs(1,:);
                                cfg.freqhigh    = curexperiment.pow.TbTFreqs(2,:);
                                cfg.method      = 'pac';
                                if l==2
                                    cfg.keeptrials  = 'no';
                                elseif l==3
                                    cfg.keeptrials  = 'yes';
                                end
                                for c=1:length(flnms)
                                    cfg.channel     = curexperiment.chngrp.(flnms{c});
                                    for i=1:length(data_freq_fou)
                                        fprintf('\nCurrent pac analysis: %s - Condition %s\n',curexperiment.dataset_name{d}(2:end),num2str(i))
                                        if not(isempty(data_freq_fou(i).freq))      
                                            crossfreq       = ft_crossfrequencyanalysis(cfg, data_freq_fou(i));
                                            if l==2
                                                % average over channels, frequencies, and phase
                                                pac.(flnms{c})(i)  = mean(crossfreq.crsspctrm,'all','omitnan');    
                                            elseif l==3
                                                % average over channels, frequencies, and phase
                                                pac.(flnms{c})  = squeeze(mean(crossfreq.crsspctrm,[2:5],'omitnan'));
                                            end
                                            clear crossfreq
                                        else
                                            pac.(flnms{c})(i)  = NaN;
                                        end
                                    end
                                end
                            else
                                if l==3
                                    if d==1
                                        trlN = curexperiment.Ntrials_enc;
                                    elseif d==2
                                        trlN = curexperiment.Ntrials_ret;
                                    end
                                    for i=1:length(flnms)
                                        pac.(flnms{i})  = NaN([trlN,1]);
                                    end
                                end
                            end
                        end
                        %% POWER ENVELOPE CORRELATION
                        % between frontal and parietal channels
                        % https://mailman.science.ru.nl/pipermail/fieldtrip/2013-March/006367.html
                        if l==2 
                            flnms = fieldnames(curexperiment.chngrp);
                            if exist('data_freq_fou','var')
                                for i=1:length(data_freq_fou)
                                    if not(isempty(data_freq_fou(i).freq))
                                        fprintf('\nCurrent pocr analysis: %s - Condition %s\n',curexperiment.dataset_name{d}(2:end),num2str(i))
                                        cur_freqdata = data_freq_fou(i);
                                        % combine the frontal channels and the parietal channels
                                        for c=1:length(flnms)
                                            cfg              = [];
                                            cfg.channel      = curexperiment.chngrp.(flnms{c});
                                            cfg.avgoverchan  = 'yes';
                                            %cfg.nanmean      = 'yes';
                                            data_frq.(flnms{c})  = ft_selectdata(cfg,cur_freqdata);
                                            data_frq.(flnms{c}).label = flnms{c};
                                        end
                                        % append the frontal and parietal averages 
                                        data_frq.cmb               = cur_freqdata;
                                        data_frq.cmb.fourierspctrm = [data_frq.(flnms{1}).fourierspctrm ,data_frq.(flnms{2}).fourierspctrm];
                                        data_frq.cmb.label         = {data_frq.(flnms{1}).label; data_frq.(flnms{2}).label};
                                        % power envelope correlation
                                        cfg             = [];
                                        cfg.method      = 'powcorr';
                                        cfg.keeptrials  = 'yes';
                                        powcorr         = ft_connectivityanalysis(cfg, data_frq.cmb);
                                        for fr=1:length(curexperiment.pow.TbTFreqs)
                                            % get the frequencies of interest
                                            cur_freqs = curexperiment.pow.TbTFreqs;
                                            cur_freq = cur_freqdata(find(~cellfun(@isempty,{powcorr.freq}),1)).freq >=cur_freqs(fr,1) & powcorr(find(~cellfun(@isempty,{powcorr.freq}),1)).freq <=cur_freqs(fr,2);
                                            % average and save
                                            if fr==1
                                                poc_t(i)  = mean(powcorr.powcorrspctrm(1,2,cur_freq),'all','omitnan'); 
                                            elseif fr==2
                                                poc_g(i)  = mean(powcorr.powcorrspctrm(1,2,cur_freq),'all','omitnan'); 
                                            end
                                        end
                                        clear powcorr
                                    else
                                        poc_t(i)  = NaN;
                                        poc_g(i)  = NaN;
                                    end
                                end
                            end
                        end
                        %% PHASE COHERENCE
                        % between frontal and parietal channels
                        % https://mailman.science.ru.nl/pipermail/fieldtrip/2013-March/006367.html
                        if l==2
                            flnms = fieldnames(curexperiment.chngrp);
                            if exist('data_freq_fou','var')
                                for i=1:length(data_freq_fou)
                                    if not(isempty(data_freq_fou(i).freq))
                                        fprintf('\nCurrent phco analysis: %s - Condition %s\n',curexperiment.dataset_name{d}(2:end),num2str(i))
                                        cur_freqdata = data_freq_fou(i);
                                        % combine the frontal channels and the parietal channels
                                        for c=1:length(flnms)
                                            cfg              = [];
                                            cfg.channel      = curexperiment.chngrp.(flnms{c});
                                            cfg.avgoverchan  = 'yes';
                                            %cfg.nanmean      = 'yes';
                                            data_frq.(flnms{c})  = ft_selectdata(cfg,cur_freqdata);
                                            data_frq.(flnms{c}).label = flnms{c};
                                        end
                                        % append the frontal and parietal averages 
                                        data_frq.cmb               = cur_freqdata;
                                        data_frq.cmb.fourierspctrm = [data_frq.(flnms{1}).fourierspctrm ,data_frq.(flnms{2}).fourierspctrm];
                                        data_frq.cmb.label         = {data_frq.(flnms{1}).label; data_frq.(flnms{2}).label};
                                        % power envelope correlation
                                        cfg             = [];
                                        cfg.method      = 'ppc';
                                        cfg.keeptrials  = 'yes';
                                        phacoh          = ft_connectivityanalysis(cfg, data_frq.cmb);
                                        for fr=1:length(curexperiment.pow.TbTFreqs)
                                            % get the frequencies of interest
                                            cur_freqs = curexperiment.pow.TbTFreqs;
                                            cur_freq = cur_freqdata(find(~cellfun(@isempty,{phacoh.freq}),1)).freq >=cur_freqs(fr,1) & phacoh(find(~cellfun(@isempty,{phacoh.freq}),1)).freq <=cur_freqs(fr,2);
                                            % average and save
                                            if fr==1
                                                phc_t(i)  = mean(phacoh.ppcspctrm(1,2,cur_freq),'all','omitnan'); 
                                            elseif fr==2
                                                phc_g(i)  = mean(phacoh.ppcspctrm(1,2,cur_freq),'all','omitnan'); 
                                            end
                                        end
                                        clear phacoh
                                    else
                                        phc_t(i)  = NaN;
                                        phc_g(i)  = NaN;
                                    end
                                end
                            end
                        end
                        clear data_frq
                        %% SAVE DATA 
                        % save l=2 & l=3 in Data for the stats later
                        %% BY-CONDITIONS
                        if l==2
                            cur_d = curexperiment.dataset_name{d}(2:end);
                            if d==1
                                % all possible behavioral combinations
                                [x1,x2,x3,x4,x5] = ndgrid({'Old','New'},{'iCor','iInc'},{'iHC','iLC'},{'sCor','sInc','NaN'},{'sHC','sLC','NaN'});
                                S.Subject       = repmat({subjectdata.nr},[length(x1(:)),1]);
                                S.OldNew        = x1(:); 
                                S.ItemAcc       = x2(:); 
                                S.ItemConf      = x3(:); 
                                S.SourceAcc     = x4(:); 
                                S.SourceConf    = x5(:); 
                                S.label         = repmat({''},[length(x1(:)),1]);
                                S.Enc_pow_theta_frontal     = NaN([length(x1(:)),1]);
                                S.Enc_pow_theta_parietal    = NaN([length(x1(:)),1]);
                                S.Enc_pow_gamma_frontal     = NaN([length(x1(:)),1]);
                                S.Enc_pow_gamma_parietal    = NaN([length(x1(:)),1]);
                                S.Enc_pek_theta_frontal     = NaN([length(x1(:)),1]);
                                S.Enc_pek_theta_parietal    = NaN([length(x1(:)),1]);
                                S.Enc_pek_gamma_frontal     = NaN([length(x1(:)),1]);
                                S.Enc_pek_gamma_parietal    = NaN([length(x1(:)),1]);
                                S.Enc_pac_frontal           = NaN([length(x1(:)),1]);
                                S.Enc_pac_parietal          = NaN([length(x1(:)),1]);
                                S.Enc_poc_theta             = NaN([length(x1(:)),1]);
                                S.Enc_poc_gamma             = NaN([length(x1(:)),1]);
                                S.Enc_phc_theta             = NaN([length(x1(:)),1]);
                                S.Enc_phc_gamma             = NaN([length(x1(:)),1]);
                                S.Ret_pow_theta_frontal     = NaN([length(x1(:)),1]);
                                S.Ret_pow_theta_parietal    = NaN([length(x1(:)),1]);
                                S.Ret_pow_gamma_frontal     = NaN([length(x1(:)),1]);
                                S.Ret_pow_gamma_parietal    = NaN([length(x1(:)),1]);
                                S.Ret_pek_theta_frontal     = NaN([length(x1(:)),1]);
                                S.Ret_pek_theta_parietal    = NaN([length(x1(:)),1]);
                                S.Ret_pek_gamma_frontal     = NaN([length(x1(:)),1]);
                                S.Ret_pek_gamma_parietal    = NaN([length(x1(:)),1]);
                                S.Ret_pac_frontal           = NaN([length(x1(:)),1]);
                                S.Ret_pac_parietal          = NaN([length(x1(:)),1]);
                                S.Ret_poc_theta             = NaN([length(x1(:)),1]);
                                S.Ret_poc_gamma             = NaN([length(x1(:)),1]);
                                S.Ret_phc_theta             = NaN([length(x1(:)),1]);
                                S.Ret_phc_gamma             = NaN([length(x1(:)),1]);
                                clear x*
                            end
                            % match the EEG variables to the behavioral ones
                            if ~(d==1 && strcmp(subjectdata.nr,'153'))
                                for i=1:length(curexperiment.(['data', num2str(d), 'l', num2str(l), '_name']))
                                    cur_lab_o = curexperiment.(['data', num2str(d), 'l', num2str(l), '_name']){i};
                                    cur_lab = split(cur_lab_o,'_');
                                    if length(cur_lab) == 3 % some items won't have source info, so we will mark these as NaN here
                                        cur_lab(4:5) = {'NaN'};
                                    end
                                    ind1 = ismember(S.OldNew,cur_lab(1));
                                    ind2 = ismember(S.ItemAcc,cur_lab(2));
                                    ind3 = ismember(S.ItemConf,cur_lab(3));
                                    ind4 = ismember(S.SourceAcc,cur_lab(4));
                                    ind5 = ismember(S.SourceConf,cur_lab(5));
                                    ind  = all([ind1,ind2,ind3,ind4,ind5],2);
                                    if sum(ind) >1
                                        error('why is there more than one match?')
                                    end
                                    % paste the variables in the right row
                                    S.label(ind)                            = {cur_lab_o};
                                    S.([cur_d,'_pow_theta_frontal'])(ind)   = pow_t.frontal(i);
                                    S.([cur_d,'_pow_theta_parietal'])(ind)  = pow_t.parietal(i);
                                    S.([cur_d,'_pow_gamma_frontal'])(ind)   = pow_g.frontal(i);
                                    S.([cur_d,'_pow_gamma_parietal'])(ind)  = pow_g.parietal(i);
                                    S.([cur_d,'_pek_theta_frontal'])(ind)   = peak_t.frontal(i);
                                    S.([cur_d,'_pek_theta_parietal'])(ind)  = peak_t.parietal(i);
                                    S.([cur_d,'_pek_gamma_frontal'])(ind)   = peak_g.frontal(i);
                                    S.([cur_d,'_pek_gamma_parietal'])(ind)  = peak_g.parietal(i);
                                    S.([cur_d,'_pac_frontal'])(ind)         = pac.frontal(i);
                                    S.([cur_d,'_pac_parietal'])(ind)        = pac.parietal(i);
                                    S.([cur_d,'_poc_theta'])(ind)           = poc_t(i);
                                    S.([cur_d,'_poc_gamma'])(ind)           = poc_g(i);
                                    S.([cur_d,'_phc_theta'])(ind)           = phc_t(i);
                                    S.([cur_d,'_phc_gamma'])(ind)           = phc_g(i);
                                end
                            end
                            clear p*
                        end
                        clear ind*
                        %% TRIAL-BY-TRIAL
                        if l==3
                            % PAC
                            if d==1
                                trlN = curexperiment.Ntrials_enc;
                            elseif d==2
                                trlN = curexperiment.Ntrials_ret;
                            end
                            % deal with empty data
                            if ~exist('data_freq_fou','var')
                                data_freq_fou.freq = 999;
                                data_freq_fou.time = 999;
                                data_freq_fou.label = {'nochan'};
                                data_freq_fou.powspctrm = NaN;
                                data_freq_fou.trialinfo = [NaN([trlN,1]) [1:trlN]'];
                            end
                            % save the variables to a table
                            T                                                                   = table;
                            T.Subject                                                           = repmat(str2double(subjectdata.nr),trlN,1);
                            T.Session                                                           = repmat(str2double(subjectdata.ses),trlN,1);
                            T.([curexperiment.dataset_name{d}(2:end) 'Trial'])                  = [1:trlN]';
                            % get the trialnumbers
                            trls = data_freq_fou.trialinfo(:,2);
                            for i=1:length(flnms)
                                T.([curexperiment.dataset_name{d}(2:end) '_' flnms{i} 'PAC'])   = NaN([trlN,1]);
                                T.([curexperiment.dataset_name{d}(2:end) '_' flnms{i} 'PAC'])(trls,:)   = pac.(flnms{i});
                            end
                            % save data in Data for the stats later
                            Data.(['PAC_' (curexperiment.dataset_name{d}(2:end))])((f-100-1)*trlN+1:(f-100)*trlN,:) = T; 
                            clear T clear flnms p*
                            % POWER
                            % deal with empty data
                            if ~exist('data_freq_pow','var')
                                data_freq_pow.freq = 999;
                                data_freq_pow.time = 999;
                                data_freq_pow.label = {'nochan'};
                                data_freq_pow.powspctrm = NaN;
                                data_freq_pow.trialinfo = [999 999];
                            end
                            if d==1
                                trlN = curexperiment.Ntrials_enc;
                            elseif d==2
                                trlN = curexperiment.Ntrials_ret;
                            end
                            % save the variables to a table
                            T = table;
                            T.Subject = repmat(str2double(subjectdata.nr),trlN,1);
                            T.Session = repmat(str2double(subjectdata.ses),trlN,1);
                            T.([curexperiment.dataset_name{d}(2:end) 'Trial'])  = [1:trlN]';
                            % freqs 
                            cur_freqs = curexperiment.pow.TbTFreqs;
                            cur_freq_t = data_freq_pow(find(~cellfun(@isempty,{data_freq_pow.freq}),1)).freq >=cur_freqs(1,1) & data_freq_pow(find(~cellfun(@isempty,{data_freq_pow.freq}),1)).freq <=cur_freqs(1,2);
                            cur_freq_g = data_freq_pow(find(~cellfun(@isempty,{data_freq_pow.freq}),1)).freq >=cur_freqs(2,1) & data_freq_pow(find(~cellfun(@isempty,{data_freq_pow.freq}),1)).freq <=cur_freqs(2,2);
                            % time 300-800 ms
                            cur_time = (data_freq_pow(find(~cellfun(@isempty,{data_freq_pow.time}),1)).time >=.3 & data_freq_pow(find(~cellfun(@isempty,{data_freq_pow.time}),1)).time <=.8);
                            % channels
                            cur_chan_f = ismember(data_freq_pow(find(~cellfun(@isempty,{data_freq_pow.label}),1)).label,curexperiment.chngrp.frontal);
                            cur_chan_p = ismember(data_freq_pow(find(~cellfun(@isempty,{data_freq_pow.label}),1)).label,curexperiment.chngrp.parietal);
                            T.([curexperiment.dataset_name{d}(2:end) '_ThetaFro']) = NaN([trlN,length(cur_freqs(1,1):cur_freqs(1,2))]);
                            T.([curexperiment.dataset_name{d}(2:end) '_GammaFro']) = NaN([trlN,length(cur_freqs(2,1):cur_freqs(2,2))]);
                            T.([curexperiment.dataset_name{d}(2:end) '_ThetaPar']) = NaN([trlN,length(cur_freqs(1,1):cur_freqs(1,2))]);
                            T.([curexperiment.dataset_name{d}(2:end) '_GammaPar']) = NaN([trlN,length(cur_freqs(2,1):cur_freqs(2,2))]);
                            for i=1:length(data_freq_pow)
                                if length(data_freq_pow(i).time)>1
                                    % get the trialnumbers
                                    trls = data_freq_pow(i).trialinfo(:,2);
                                    % add the power to the corresponding trials
                                    T.([curexperiment.dataset_name{d}(2:end) '_ThetaFro'])(trls,:,:) = squeeze(mean(data_freq_pow(i).powspctrm(:,cur_chan_f,cur_freq_t,cur_time),[2,4],'omitnan'));
                                    T.([curexperiment.dataset_name{d}(2:end) '_GammaFro'])(trls,:,:) = squeeze(mean(data_freq_pow(i).powspctrm(:,cur_chan_f,cur_freq_g,cur_time),[2,4],'omitnan'));
                                    T.([curexperiment.dataset_name{d}(2:end) '_ThetaPar'])(trls,:,:) = squeeze(mean(data_freq_pow(i).powspctrm(:,cur_chan_p,cur_freq_t,cur_time),[2,4],'omitnan'));
                                    T.([curexperiment.dataset_name{d}(2:end) '_GammaPar'])(trls,:,:) = squeeze(mean(data_freq_pow(i).powspctrm(:,cur_chan_p,cur_freq_g,cur_time),[2,4],'omitnan'));
                                    clear trls
                                end
                            end
                            % save data in Data for the stats later
                            Data.(['TF_' (curexperiment.dataset_name{d}(2:end))])((f-100-1)*trlN+1:(f-100)*trlN,:) = T; 
                            clear T
                        else
                            % save the data for l<3
                            if exist('data_freq_pow','var')
                                for i=1:length(data_freq_pow)
                                    if ~isempty(data_freq_pow(i).label)
                                        data_cond = data_freq_pow(i);
%                                         if ~contains('powspctrm',fieldnames(data_cond))
%                                             data_cond.powspctrm = abs(data_cond.fourierspctrm).^2;
%                                         end
                                        data_cond.cfg.previous = [];
                                        save([curexperiment.analysis_loc filesep subjectdata.nr '_' subjectdata.ses curexperiment.dataset_name{d} curexperiment.(['data', num2str(d)]).level_name{l} '_' curexperiment.(['data',num2str(d),'l',num2str(l),'_name']){i} '_Total_TF'],'-v7.3','data_cond'); % added the '-v7.3' because of the size of the files
                                    end
                                end
                            end
                        end
                        clear cur_freq* cur_time cur_chan* trlN
                        clear data_cond data_freq* curdat* curcond data
                    end
                end
                if l==2
                    T = struct2table(S);
                    % save data in Data for the stats later
                    Data.TF_byCond((f-100-1)*size(T,1)+1:(f-100)*size(T,1),:) = T;
                end
                clear T S
                save(curexperiment.outputfile,'Data')
                save(strcat(curexperiment.outputfile(1:end-4),'_',string(datetime('today')),'.mat'),'Data')
                clear d fr l do_plt
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EEG GROUP-LEVEL ANALYSES %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\n##########################\n')
fprintf('## GROUP-LEVEL ANALYSES ##\n')
fprintf('##########################\n')
if any(contains(curexperiment.analyses,'ga'))
    % OBTAIN THE DATA FOR THE TRIAL BY TRIAL ANALYSIS FOR R
    % get the data
    load(curexperiment.outputfile)
    if any(contains(curexperiment.analyses,'erpga'))
        %% ERP
        % find the files in the analyses folder
        ERPdir          = fullfile(curexperiment.analysis_loc, sprintf('*ERP*'));
        ERPdf           = dir(ERPdir);
        ERPfiles        = {ERPdf.name};
        ERPfiles        = ERPfiles(~contains(ERPfiles,'GrandAverage'));
        clear ERPdir  
        % loop over sessions
        for ses=1%:curexperiment.Nses
            ERPfiles_ses = ERPfiles(contains(ERPfiles,['_' num2str(ses) '_']));
            % loop over datasets
            for d=2%1:length(curexperiment.datasets_names)-1 % exclude Rest EEG
                % find the files from the current dataset
                ERPfiles_dat = ERPfiles_ses(contains(ERPfiles_ses,curexperiment.dataset_name{d}));
                for l=1%:curexperiment.(['data' num2str(d)]).levels
                    curlevname = curexperiment.(['data' num2str(d)]).level_name{l}(2:end);
                    fprintf('\n%s\n',curlevname)
                    ERPfiles_lev = ERPfiles_dat(contains(ERPfiles_dat,[curexperiment.(['data' num2str(d)]).level_name{l}, '_']));
                    % define the current dataset
                    curdat = curexperiment.(['data' num2str(d)]).(['l' num2str(l)]);
                    curconname = curexperiment.(['data' num2str(d) 'l' num2str(l) '_name']);
                    curconname = strrep(curconname,'_','');
                    % loop over conditions
                    for c=1:length(fieldnames(curdat))
                        % create an array to hold the inputfiles for the different conditions
                        curmatfiles = ERPfiles_lev(contains(ERPfiles_lev,strcat(curconname{c},'_')));
                        for cf=1:length(curmatfiles)
                            inputfiles(c,cf) = fullfile(curexperiment.analysis_loc, curmatfiles(cf));
                        end
                    end
                    clear c curdat curfiles cf
                    cfg = [];
                    cfg.keepindividual = 'no';
                    % calculate the grand average per condition
                    % the timewindow is equal to that of the shortest of the inputs
                    for i=1:size(inputfiles,1)
                        cfg.inputfile = inputfiles(i, ~cellfun('isempty',inputfiles(i,:)));
                        fprintf('\nERP analysis %s %s\n',curexperiment.dataset_name{d}(2:end), curconname{i});
                        GrandAverage_ERP(i) = ft_timelockgrandaverage(cfg);
                    end
                    % PLOTTING
                    cfg                  = [];
                    % I changed the FT_singleplotER function to make this
                    % work (the cfg.parameter needs to be changed to 'avg'
                    % after the preprocessing).
                    cfg.preproc.lpfilter = 'yes';
                    cfg.preproc.lpfreq   = 30;
                    cfg.xlim             = [-.2 1];
                    % https://i.pinimg.com/originals/eb/0c/a4/eb0ca44520b566e9a11b607c0334b07f.gif
                    if l==1
                       cfg.linecolor  = [
                                        000 102 015;...   % green (dark)    - sHit
                                        000 204 051;...   % green (light)   - sMiss
                                        255 235 000;...   % yellow          - iHit
                                        204 000 000;...   % red             - iMiss
                                        255 153 153;...   % pink            - FA
                                        255 102 000;...   % orange          - CR
                                        102 000 204;...   % purple (dark)   - sHC
                                        153 102 255;...   % purple (light)  - sLC
                                        000 051 255;...   % blue (dark)     - iHC
                                        000 153 255 ...   % blue (light)    - iLC
                                        ]./255; 
                    else
                    end
                    cfg.fontsize    = 20;
                    cfg.linewidth   = 3;
                    for ch=1:length(fieldnames(curexperiment.chngrp))
                        cnnames = fieldnames(curexperiment.chngrp);
                        cfg.channel = curexperiment.chngrp.(cnnames{ch});
                        %curdata=sprintf('GrandAverage_ERP(%d),',1:length(GrandAverage_ERP));
                        slctn = find(contains(curconname, {'iHit','CR'}));
                        curdata=sprintf('GrandAverage_ERP(%d),',slctn);
                        eval(sprintf('ft_singleplotER(cfg,%s)',curdata(1:end-1)))
                        clear curdata
                        hold on
                        ylim 'auto'
                        cfg.ylim = get(gca,'YLim');
                        plot([0,0],cfg.ylim,'k-','LineWidth',1.5); % plot straight line at stimulus onset
                        hold on
                        cur_leg = legend(curconname{slctn},'Location','southoutside','Orientation','horizontal');
                        set(cur_leg, 'FontSize',cfg.fontsize);
                        clear r
                        set(gcf, 'Position', get(0, 'Screensize')) % make full screen
                        title({curexperiment.dataset_name{d}(2:end);cnnames{ch};});
                        pause(0.05) %in seconds %to save correctly
                        exportgraphics(gcf,[curexperiment.analysis_loc filesep 'Plots' filesep 'Session' num2str(ses) curexperiment.dataset_name{d} '_' cnnames{ch} curexperiment.(['data' num2str(d)]).level_name{l} '_ERP_Singleplot_GrandAverage.png']);
                        exportgraphics(gcf,[curexperiment.analysis_loc filesep 'Plots' filesep 'Archive' filesep 'Session' num2str(ses) curexperiment.dataset_name{d} '_' cnnames{ch} curexperiment.(['data' num2str(d)]).level_name{l} '_ERP_Singleplot_GrandAverage_' date '.png']);
                        close all
                        clear curdata cur_leg
                    end 
                    clear inputfiles i
                    % save the data
                    for i=1:length(GrandAverage_ERP)
                        data_cond = GrandAverage_ERP(i);
                        data_cond.cfg.previous = []; % clear previous
                        fprintf('\nSaving %s %s\n',curexperiment.dataset_name{d}(2:end), curconname{i});
                        save([curexperiment.analysis_loc filesep '_Session' num2str(ses) curexperiment.dataset_name{d} '_' curlevname curexperiment.(['data' num2str(d) 'l' num2str(l) '_name']){i} '_ERP_GrandAverage'],'data_cond');
                    end
                    clear data_cond i GrandAverage_ERP curconname
                end
                clear ERPfiles_lev
            end
            clear ERPfiles_dat
        end
        clear ERPfiles_ses ERPfiles curmatfiles d g i ERP* ses
    end

    %% POWER
    % make the average of the oscillatory data (over frequencies(, or more))
    % make a placeholder for the average values
    Data.TF_Ret_avg = table;
    Data.TF_Ret_avg = Data.TF_Ret(:,1:3);
    % do the averages over the frequencies, separate for brain regions
    for v=4:length(Data.TF_Ret.Properties.VariableNames)
        Data.TF_Ret_avg.(Data.TF_Ret.Properties.VariableNames{v}) = mean(Data.TF_Ret.(Data.TF_Ret.Properties.VariableNames{v}),2,'omitnan');
    end
    % add the new variables to the RET table as a placeholder
    flnms = fieldnames(Data.TF_Ret_avg);
    pow_flnms = flnms(contains(flnms,{'Theta','Gamma'}));
    for i=1:length(pow_flnms)
        Data.RET.(pow_flnms{i}) = NaN(size(Data.RET,1),1);
    end
    % add the oscillatory data to the behavioral trial data
    for i=1:size(Data.TF_Ret_avg,1)
        % find the corresponding trial in RET
        idx = all([Data.RET.Subject Data.RET.Session Data.RET.RetTrial]==[Data.TF_Ret_avg.Subject(i) Data.TF_Ret_avg.Session(i) Data.TF_Ret_avg.RetTrial(i)],2);
        % sanity check
        if sum(idx)>1
            error('There should be just a single trial that matches, but there are more. Check what is happening')
        end
        Data.RET(idx,pow_flnms) = Data.TF_Ret_avg(i,pow_flnms);
    end
    %% PHASE AMPLITUDE COUPLING
    % add the new variables to the RET table as a placeholder
    flnms_enc = fieldnames(Data.PAC_Enc);
    pac_flnms_enc = flnms_enc(contains(flnms_enc,{'PAC'}));
    flnms_ret = fieldnames(Data.PAC_Ret);
    pac_flnms_ret = flnms_ret(contains(flnms_ret,{'PAC'}));
    pac_flnms = [pac_flnms_enc;pac_flnms_ret];
    for i=1:length(pac_flnms)
        Data.RET.(pac_flnms{i}) = NaN(size(Data.RET,1),1);
    end
    % add the pac data to the behavioral trial data
    for i=1:size(Data.RET,1)
        % find the corresponding retrieval trial
        idx_ret = all([Data.RET.Subject(i) Data.RET.Session(i) Data.RET.RetTrial(i)]==[Data.PAC_Ret.Subject Data.PAC_Ret.Session Data.PAC_Ret.RetTrial],2);
        if sum(idx_ret)==1
            Data.RET(i,pac_flnms_ret) = Data.PAC_Ret(idx_ret,pac_flnms_ret);
        elseif sum(idx_ret)>1
            error('There should be just a single trial that matches, but there are more. Check what is happening')
        end
        % find the corresponding encoding trial
        idx_enc = all([Data.RET.Subject(i) Data.RET.Session(i) Data.RET.EncTrial(i)]==[Data.PAC_Enc.Subject Data.PAC_Enc.Session Data.PAC_Enc.EncTrial],2);
        if sum(idx_enc)==1
            Data.RET(i,pac_flnms_enc) = Data.PAC_Enc(idx_enc,pac_flnms_enc);
        elseif sum(idx_enc)>1
            error('There should be just a single trial that matches, but there are more. Check what is happening')
        end
    end
    % save the trial by trial data
    try
        writetable(Data.TF_byCond,[curexperiment.datafolder_output filesep 'Data_TbC.csv'])
    catch
    end
    try
        writetable(Data.RET,[curexperiment.datafolder_output filesep 'Data_TbT.csv'])
    catch
    end
    save(curexperiment.outputfile,'Data')
    outputfile2 = strcat(curexperiment.outputfile(1:end-4),'_',string(datetime('today')),'.mat');
    save(outputfile2,'Data')
    fprintf('\nANALYSIS DATA SAVED')
end

%%%%%%%%%%%%%%%%%%
%% EEG PLOTTING %%
%%%%%%%%%%%%%%%%%%
fprintf('\n###############\n')
fprintf('## PLOTTING ##\n')
fprintf('##############\n')
% Calculate the grand average and use them for plotting
if any(contains(curexperiment.analyses,'plt'))
    set(groot,'DefaultFigureRenderer','painters')
    % loop over the sessions
    for ses=1%:curexperiment.Nses
        % loop over the datasets
        for d=2%1:length(curexperiment.dataset_name)-1
            % loop over all levels of processing
            for l=1%1:curexperiment.(['data' num2str(d)]).levels-1
                TFfiles_lev = {dir([curexperiment.analysis_loc filesep '*_' num2str(ses)...
                    curexperiment.dataset_name{d} curexperiment.(['data', num2str(d)]).level_name{l} '_'...
                    '*TF.mat']).name};
                if ~isempty(TFfiles_lev)
                    %% TRIAL SELECTION
                    % Select trials based upon condition
                    curconname = curexperiment.(['data' num2str(d) 'l' num2str(l) '_name']);
                    for c=1:length(curconname)
                        % create an array to hold the inputfiles for the different conditions
                        curmatfiles = TFfiles_lev(contains(TFfiles_lev,['_',curconname{c},'_']));
                        % exclude participants
                        excl_ppn = {};
                        curmatfiles = curmatfiles(~contains(curmatfiles,excl_ppn));
                        for cf=1:length(curmatfiles)
                            inputfiles(c,cf) = fullfile(curexperiment.analysis_loc, curmatfiles(cf));
                        end
                        clear curmatfiles
                    end
                    clear c curdat curfiles cf
                    %% GRAND AVERAGE
                    cfg = [];
                    cfg.keepindividual = 'yes';
                    % calculate the grand average
                    % the timewindow is equal to that of the shortest of the inputs
                    for i=1:size(curconname,2)
                        cfg.inputfile = inputfiles(i, ~cellfun('isempty',inputfiles(i,:)));
                        fprintf('\nMaking grand average %s %s\n',curexperiment.dataset_name{d}(2:end),curconname{i});
                        try
                            GrandAverage_TF(i) = ft_freqgrandaverage(cfg);
                        catch
                            % if there are no subjects, create an empty struct
                            f = fieldnames(GrandAverage_TF(1))';
                            f{2,1}={[]};
                            GrandAverage_TF(i) = struct(f{:});
                        end
                    end
                    clear inputfiles i
                    %% DIFFERENCE PLOTS (SINGLEPLOTS & TOPOPLOTS)
                    try
                        % CALCULATE THE DIFFERENCE
                        if d==2 && l==1 
                           cur_diff = [3 6]; limt = [1, 1; .005, .04]; limt_t = [3 3; 3 3];% iHit-CR
                        else 
                        end
                        % store only the datasets that are currently needed
                        temp_GA = GrandAverage_TF(cur_diff);
                        % compare the participants that are in common between the GAs
                        inp.GA1 = temp_GA(1).cfg.inputfile;
                        inp.GA2 = temp_GA(2).cfg.inputfile;
                        % get only the numbers
                        inp.GA1 = regexp(inp.GA1,'\d+','match');
                        inp.GA2 = regexp(inp.GA2,'\d+','match');
                        inp.GA1 = str2double(vertcat(inp.GA1{:}));
                        inp.GA2 = str2double(vertcat(inp.GA2{:}));
                        % indicate which participants in the first dataset are not present in the second
                        idx = any(~ismember(inp.GA1,inp.GA2),2);
                        % and remove these from the first dataset
                        temp_GA(1).powspctrm(idx,:,:,:)=[];
                        clear idx
                        % indicate which participants in the second dataset are not present in the first
                        idx = any(~ismember(inp.GA2,inp.GA1),2);
                        % and remove these fromt the second dataset
                        temp_GA(2).powspctrm(idx,:,:,:)=[];
                        clear idx inp
                        % compute the difference
                        cfg = [];
                        cfg.parameter = 'powspctrm';
                        cfg.operation = 'subtract'; %'x1 - x2';
                        TF_diff = ft_math(cfg,temp_GA(1),temp_GA(2));
                        % compute the t-values
                        [h,p,ci,stats] = ttest(temp_GA(1).powspctrm,temp_GA(2).powspctrm);
                        TF_diff_t = temp_GA(1);
                        TF_diff_t.powspctrm = stats.tstat;
                        clear temp_GA
                        % SINGLEPLOT (difference plot)
                        % compute the average over trials
                        cfg = [];
                        GA_TF_plot = ft_freqdescriptives(cfg,TF_diff);
                        GA_TF_plot_t = ft_freqdescriptives(cfg,TF_diff_t);
                        % no baseline correction for the difference plots
                        for fr=1:size(curexperiment.pow.TbTFreqs,2)
                            cngrps = fieldnames(curexperiment.chngrp);
                            for cn=1:length(cngrps)
                                % plot both the raw values and the t-values
                                for ty=2%1:2
                                    if ty==1
                                        cur_data = GA_TF_plot;
                                    elseif ty==2
                                        cur_data = GA_TF_plot_t;
                                    end
                                    % compute the mean over the relevant channels
                                    pl.timewin = [-.2 1];
                                    if fr==1
                                        pl.freqwin = [1 13];
                                        cur_freq = 'theta';
                                        lim = limt(1,1);
                                        lim_t = limt_t(1,1);
                                        tcks = 2;
                                    elseif fr==2
                                        pl.freqwin = [30 55];
                                        cur_freq = 'gamma';
                                        lim = limt(2,1);
                                        lim_t = limt_t(2,1);
                                        tcks = 5;
                                    end
                                    pl.time = logical(cur_data.time >= pl.timewin(1) & cur_data.time <= pl.timewin(2));
                                    pl.freq = logical(cur_data.freq >= pl.freqwin(1) & cur_data.freq <= pl.freqwin(2));
                                    pl.chan = ismember(cur_data.label,curexperiment.chngrp.(cngrps{cn}));
                                    meanpow = squeeze(mean(cur_data.powspctrm(pl.chan,pl.freq,pl.time),1,'omitnan'));
                                    pl.time = cur_data.time(pl.time);
                                    pl.freq = cur_data.freq(pl.freq);
                                    % interpolate to make the plot smooth
                                    pl.time_inp = linspace(pl.timewin(1),pl.timewin(2),512);
                                    pl.freq_inp = linspace(pl.freqwin(1),pl.freqwin(2),512);
                                    [org_time,org_freq] = meshgrid(pl.time,pl.freq);
                                    [inp_time,inp_freq] = meshgrid(pl.time_inp,pl.freq_inp);
                                    meanpow = interp2(org_time,org_freq,meanpow,inp_time,inp_freq,'spline');
                                    % start the figure
                                    fig = figure();
                                    imagesc(pl.time,pl.freq,meanpow);
                                    xlim(pl.timewin);
                                    axis xy
                                    % add the labels
                                    xlabel('Time (s)')
                                    ylabel('Frequency (Hz)')
                                    % adjust the ticks
                                    yticks(pl.freq(1):tcks:pl.freq(end))
                                    % make the colorbar
                                    colormap(brewermap(256,'*RdBu'));
                                    h = colorbar();
                                    h.Location = 'southoutside';
                                    if ty==1
                                        clim([-lim lim+lim*.01]);
                                        ylabel(h,'Power difference (µV)')
                                    elseif ty==2
                                        clim([-lim_t lim_t+lim_t*.01]);
                                        ylabel(h,'t-values')
                                    end
                                    % line at t=0
                                    hold on;
                                    xline(0,'--k','LineWidth',2);
                                    % add lines to show important area
                                    if fr==1
                                        hold on;
                                        yline(curexperiment.pow.TbTFreqs(1,1),':k','LineWidth',2);
                                        hold on;
                                        yline(curexperiment.pow.TbTFreqs(1,2),':k','LineWidth',2);
                                    elseif fr==2
                                        hold on;
                                        yline(curexperiment.pow.TbTFreqs(2,1),':k','LineWidth',2);
                                        hold on;
                                        yline(curexperiment.pow.TbTFreqs(2,2),':k','LineWidth',2);
                                    end
                                    % change the font and plot size
                                    fontsize(fig, 25,"points")
                                    set(gcf,'Position', [1,1,500,1000])
                                    pause(.5)
                                    %set(gcf,'Position', get(0, 'Screensize')) % make full screen
                                    if ty==1
                                        exportgraphics(fig,[curexperiment.analysis_loc filesep 'Plots' filesep 'Session' num2str(ses) curexperiment.dataset_name{d} '_TF_DiffSing_GrandAverage' curexperiment.(['data' num2str(d)]).level_name{l} '_' cur_freq '_' cngrps{cn} '_' curconname{cur_diff(1)} '-' curconname{cur_diff(2)} '.png']);
                                        exportgraphics(fig,[curexperiment.analysis_loc filesep 'Plots' filesep 'Archive' filesep 'Session' num2str(ses) curexperiment.dataset_name{d} '_TF_DiffSing_GrandAverage' curexperiment.(['data' num2str(d)]).level_name{l} '_' cur_freq '_' cngrps{cn} '_' curconname{cur_diff(1)} '-' curconname{cur_diff(2)} '_' date '.png']);
                                    elseif ty==2
                                        exportgraphics(fig,[curexperiment.analysis_loc filesep 'Plots' filesep 'Session' num2str(ses) curexperiment.dataset_name{d} '_TF_DiffSing_GrandAverage_tvals' curexperiment.(['data' num2str(d)]).level_name{l} '_' cur_freq '_' cngrps{cn} '_' curconname{cur_diff(1)} '-' curconname{cur_diff(2)} '.png']);
                                        exportgraphics(fig,[curexperiment.analysis_loc filesep 'Plots' filesep 'Archive' filesep 'Session' num2str(ses) curexperiment.dataset_name{d} '_TF_DiffSing_GrandAverage_tvals' curexperiment.(['data' num2str(d)]).level_name{l} '_' cur_freq '_' cngrps{cn} '_' curconname{cur_diff(1)} '-' curconname{cur_diff(2)} '_' date '.png']);
                                    end
                                end
                            end
                            clear pl meanpow cur_freq
                        end
                        close all
                        % TOPOPLOT (difference plot)
                        cfg               = [];
                        cfg.xlim          = [.3 .8];
                        cfg.colorbar      = 'yes';
                        cfg.interactive   = 'no';
                        cfg.layout        = curexperiment.elec.lay;
                        cfg.renderer      = 'painters';                               
                        cfg.marker        = 'off';
                        cfg.highlight     = 'on';
                        cfg.highlightchannel = [curexperiment.chngrp.(cngrps{1}), curexperiment.chngrp.(cngrps{2})];
                        cfg.highlightsymbol = '*';
                        cfg.highlightsize = 15;
                        cfg.comment       = 'no';
                        cfg.gridscale     = 100;
                        for fr=1:size(curexperiment.pow.TbTFreqs,2)
                            % plot both the raw values and the t-values
                            for ty=2%1:2
                                if fr==1
                                    cfg.ylim = curexperiment.pow.TbTFreqs(1,:); % theta range
                                    cur_freq = 'theta';
                                    lim      = limt(1,2);
                                    lim_t    = limt_t(1,2);
                                elseif fr==2
                                    cfg.ylim = curexperiment.pow.TbTFreqs(2,:); % gamma range
                                    cur_freq = 'gamma';
                                    lim      = limt(2,2);
                                    lim_t    = limt_t(2,2);
                                end
                                if ty==1
                                    cur_data = GA_TF_plot;
                                    cfg.zlim = [-lim lim];
                                elseif ty==2
                                    cur_data = GA_TF_plot_t;
                                    cfg.zlim = [-lim_t lim_t];
                                end
                                figure;ft_topoplotTFR(cfg,cur_data);
                                colormap(brewermap(256,'*RdBu'));
                                h = colorbar();
                                ylabel(h,'t-values')
                                h.Location = 'southoutside';
                                set(gca,'fontsize',25)
                                set(gcf,'Position', [1,1,1000,1000])
                                pause(0.05) % in seconds %to save correctly
                                if ty==1
                                    exportgraphics(gcf,[curexperiment.analysis_loc filesep 'Plots' filesep 'Session' num2str(ses) curexperiment.dataset_name{d} '_TF_DiffTopo_GrandAverage' curexperiment.(['data' num2str(d)]).level_name{l} '_' cur_freq '_' curconname{cur_diff(1)} '-' curconname{cur_diff(2)} '.png']);
                                    exportgraphics(gcf,[curexperiment.analysis_loc filesep 'Plots' filesep 'Archive' filesep 'Session' num2str(ses) curexperiment.dataset_name{d} '_TF_DiffTopo_GrandAverage' curexperiment.(['data' num2str(d)]).level_name{l} '_' cur_freq '_' curconname{cur_diff(1)} '-' curconname{cur_diff(2)} '_' date '.png']);
                                elseif ty==2
                                    exportgraphics(gcf,[curexperiment.analysis_loc filesep 'Plots' filesep 'Session' num2str(ses) curexperiment.dataset_name{d} '_TF_DiffTopo_GrandAverage_tvals' curexperiment.(['data' num2str(d)]).level_name{l} '_' cur_freq '_' curconname{cur_diff(1)} '-' curconname{cur_diff(2)} '.png']);
                                    exportgraphics(gcf,[curexperiment.analysis_loc filesep 'Plots' filesep 'Archive' filesep 'Session' num2str(ses) curexperiment.dataset_name{d} '_TF_DiffTopo_GrandAverage_tvals' curexperiment.(['data' num2str(d)]).level_name{l} '_' cur_freq '_' curconname{cur_diff(1)} '-' curconname{cur_diff(2)} '_' date '.png']);
                                end
                                clear cur_freq
                            end
                        end
                        clear TF_diff cur_diff
                        close all
                    catch
                        fprintf('\nCould not make difference plot\n\n')
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%% BEHAVIORAL ANALYSES %%
%%%%%%%%%%%%%%%%%%%%%%%%%
if any(contains(curexperiment.analyses,'beh'))
    % get the data
    load(curexperiment.outputfile)
    %% ENCODING
    fprintf('\nENCODING SESSION 1:')
    % Imagining success rate
    fprintf('\nThe average imaging success rate was: %.4f (SD %.4f)', mean(Data.Kahn7Study.session1.EncAllImgSucc./3,'omitnan'), std(Data.Kahn7Study.session1.EncAllImgSucc./3,'omitnan'))
    fprintf('\nThe pleasant imaging success rate was: %.4f (SD %.4f)', mean(Data.Kahn7Study.session1.EncPleImgSucc./3,'omitnan'), std(Data.Kahn7Study.session1.EncPleImgSucc./3,'omitnan'))
    fprintf('\nThe place imaging success rate was: %.4f (SD %.4f)', mean(Data.Kahn7Study.session1.EncPlaImgSucc./3,'omitnan'), std(Data.Kahn7Study.session1.EncPlaImgSucc./3,'omitnan'))
    [~,p,~,stat] = ttest(Data.Kahn7Study.session1.EncPleImgSucc./3, Data.Kahn7Study.session1.EncPlaImgSucc./3);
    D = meanEffectSize(Data.Kahn7Study.session1.EncPleImgSucc./3,Data.Kahn7Study.session1.EncPlaImgSucc./3,Effect="cohen");
    fprintf('\nt-test pleasant vs. place imaging success: t(%d)=%.2f, p = %.4f, D = %.4f',stat.df, stat.tstat, p, D.Effect)
    % Reaction time
    fprintf('\n\nThe RT was: %.2f (SD %.2f)', mean(Data.Kahn7Study.session1.EncAllRT,'omitnan'), std(Data.Kahn7Study.session1.EncAllRT,'omitnan'))
    fprintf('\nThe pleasant RT was: %.2f (SD %.2f)', mean(Data.Kahn7Study.session1.EncPleRT,'omitnan'), std(Data.Kahn7Study.session1.EncPleRT,'omitnan'))
    fprintf('\nThe place RT was: %.2f (SD %.2f)', mean(Data.Kahn7Study.session1.EncPlaRT,'omitnan'), std(Data.Kahn7Study.session1.EncPlaRT,'omitnan'))
    [~,p,~,stat] = ttest(Data.Kahn7Study.session1.EncPleRT, Data.Kahn7Study.session1.EncPlaRT);
    D = meanEffectSize(Data.Kahn7Study.session1.EncPleRT,Data.Kahn7Study.session1.EncPlaRT,Effect="cohen");
    fprintf('\nt-test pleasant vs. place RT: t(%d)=%.2f, p = %.4f, D = %.4f',stat.df, stat.tstat, p, D.Effect)
    %% RETRIEVAL
    fprintf('\nRETRIEVAL SESSION 1:')
    % Item memory
    fprintf('\nThe average HR was: %.4f (SD %.4f)', mean(Data.Kahn7Test.session1.ihit_rate,'omitnan'), std(Data.Kahn7Test.session1.ihit_rate,'omitnan'))
    fprintf('\nThe average FA was: %.4f (SD %.4f)', mean(Data.Kahn7Test.session1.fa_rate,'omitnan'), std(Data.Kahn7Test.session1.fa_rate,'omitnan'))
    fprintf('\nThe average d-prime was: %.4f (SD %.4f)', mean(Data.Kahn7Test.session1.d_prime,'omitnan'), std(Data.Kahn7Test.session1.d_prime,'omitnan'))
    fprintf('\nThe pleasant d-prime was: %.4f (SD %.4f)', mean(Data.Kahn7Test.session1.d_prime_pleasant,'omitnan'), std(Data.Kahn7Test.session1.d_prime_pleasant,'omitnan'))
    fprintf('\nThe place d-prime was: %.4f (SD %.4f)', mean(Data.Kahn7Test.session1.d_prime_place,'omitnan'), std(Data.Kahn7Test.session1.d_prime_place,'omitnan'))
    [~,p,~,stat] = ttest(Data.Kahn7Test.session1.d_prime_pleasant, Data.Kahn7Test.session1.d_prime_place);
    D = meanEffectSize(Data.Kahn7Test.session1.d_prime_pleasant,Data.Kahn7Test.session1.d_prime_place,Effect="cohen");
    fprintf('\nt-test pleasant vs. place d-prime: t(%d)=%.2f, p = %.4f, D = %.4f',stat.df, stat.tstat, p, D.Effect)
    fprintf('\nThe average percentage of high-confident responses was: %.4f (SD %.4f)', mean(Data.Kahn7Test.session1.ON_HC_perc,'omitnan'), std(Data.Kahn7Test.session1.ON_HC_perc,'omitnan'))
    fprintf('\nThe average percentage of high-confident hits was: %.4f (SD %.4f)', mean(Data.Kahn7Test.session1.ON_HChit_perc,'omitnan'), std(Data.Kahn7Test.session1.ON_HChit_perc,'omitnan'))
    fprintf('\nThe average percentage of high-confident correct rejections was: %.4f (SD %.4f)', mean(Data.Kahn7Test.session1.ON_HCcr_perc,'omitnan'), std(Data.Kahn7Test.session1.ON_HCcr_perc,'omitnan'))
    % Item Reaction Time
    fprintf('\n\nThe RT was: %.2f (SD %.2f)', mean(Data.Kahn7Test.session1.ONAll_RT,'omitnan'), std(Data.Kahn7Test.session1.ONAll_RT,'omitnan'))
    fprintf('\nThe pleasant RT was: %.2f (SD %.2f)', mean(Data.Kahn7Test.session1.ONPleasant_RT,'omitnan'), std(Data.Kahn7Test.session1.ONPleasant_RT,'omitnan'))
    fprintf('\nThe place RT was: %.2f (SD %.2f)', mean(Data.Kahn7Test.session1.ONPlace_RT,'omitnan'), std(Data.Kahn7Test.session1.ONPlace_RT,'omitnan'))
    [~,p,~,stat] = ttest(Data.Kahn7Test.session1.ONPleasant_RT, Data.Kahn7Test.session1.ONPlace_RT);
    D = meanEffectSize(Data.Kahn7Test.session1.ONPleasant_RT,Data.Kahn7Test.session1.ONPlace_RT,Effect="cohen");
    fprintf('\nt-test pleasant vs. place RT: t(%d)=%.2f, p = %.4f, D = %.4f',stat.df, stat.tstat, p, D.Effect)
    % Source memory
    fprintf('\nThe average HR was: %.4f (SD %.4f)', mean(Data.Kahn7Test.session1.shit_rate,'omitnan'), std(Data.Kahn7Test.session1.shit_rate,'omitnan'))
    fprintf('\nThe average FA was: %.4f (SD %.4f)', mean(Data.Kahn7Test.session1.sfa_rate,'omitnan'), std(Data.Kahn7Test.session1.sfa_rate,'omitnan'))
    fprintf('\nThe average d-prime was: %.4f (SD %.4f)', mean(Data.Kahn7Test.session1.d_prime_source,'omitnan'), std(Data.Kahn7Test.session1.d_prime_source,'omitnan'))
    fprintf('\nThe average percentage of high-confident responses was: %.4f (SD %.4f)', mean(Data.Kahn7Test.session1.Source_HC_perc,'omitnan'), std(Data.Kahn7Test.session1.Source_HC_perc,'omitnan'))
    fprintf('\nThe average percentage of high-confident hits was: %.4f (SD %.4f)', mean(Data.Kahn7Test.session1.Source_HChit_perc,'omitnan'), std(Data.Kahn7Test.session1.Source_HChit_perc,'omitnan'))
    % Source Reaction Time
    fprintf('\n\nThe RT was: %.2f (SD %.2f)', mean(Data.Kahn7Test.session1.SourceAll_RT,'omitnan'), std(Data.Kahn7Test.session1.SourceAll_RT,'omitnan'))
    fprintf('\nThe pleasant RT was: %.2f (SD %.2f)', mean(Data.Kahn7Test.session1.SourcePleasant_RT,'omitnan'), std(Data.Kahn7Test.session1.SourcePleasant_RT,'omitnan'))
    fprintf('\nThe place RT was: %.2f (SD %.2f)', mean(Data.Kahn7Test.session1.SourcePlace_RT,'omitnan'), std(Data.Kahn7Test.session1.SourcePlace_RT,'omitnan'))
    [h,p,ci,stat] = ttest(Data.Kahn7Test.session1.SourcePleasant_RT, Data.Kahn7Test.session1.SourcePlace_RT);
    D = meanEffectSize(Data.Kahn7Test.session1.SourcePleasant_RT,Data.Kahn7Test.session1.SourcePlace_RT,Effect="cohen");
    fprintf('\nt-test pleasant vs. place RT: t(%d)=%.2f, p = %.4f, D = %.4f',stat.df, stat.tstat, p, D.Effect)


end

%%%%%%%%%%%%%%%%%%%%%%
%% CONTROL ANALYSES %%
%%%%%%%%%%%%%%%%%%%%%%
if any(contains(curexperiment.analyses,'ctl'))
    %% CHECK THE DIFFERENCES IN NUMBER OF BLINKS BETWEEN CONDITIONS
    blnks.data = Data.Blinks.Ret; % where the blink data is stored. Struct with for each participant a table with the conditions as rows
    blnks.ppns  = fieldnames(Data.Blinks.Ret); % participants
    blnks.conds = Data.Blinks.Ret.ppn_101.Properties.RowNames; % conditions
    blnks.blinks = {'BlinksMin', 'BlinksMinGood'}; % cell array of blink condition(s) to compare column labels of the table
    
    % make the table to store the output in
    for c=1:length(blnks.conds)
        blnks.outp.(blnks.conds{c}) = table('Size', [length(blnks.ppns),length(blnks.blinks)], ...
                                            'VariableNames',blnks.blinks, ...
                                            'VariableTypes',repmat({'double'},1,length(blnks.blinks)), ...
                                            'RowNames',blnks.ppns);
    end
    % loop over the participants
    for p=1:length(blnks.ppns)
        % loop over the conditions and store the blink values for that participant and condition
        for c=1:length(blnks.conds)
            % add the values to the output table
            blnks.outp.(blnks.conds{c})(blnks.ppns{p},:) = blnks.data.(blnks.ppns{p})(blnks.conds{c},blnks.blinks);
        end
    end
    % do the stats
    cond1 = 'iHitHC';    % condition to compare
    cond2 = 'iHitLC';    % condition to compare
    var   = 'BlinksMin'; % variable to research
    [h,p] = ttest(blnks.outp.(cond1).(var),blnks.outp.(cond2).(var))
    clear cond1 cond2 var h p
    %% COMPARE THE BLINK TFR BETWEEN CONDITIONS
    % loop over participants
    for f=101:curexperiment.Nsubs+100
        for s=1:curexperiment.Nses
            %% SUBJECT DATA
            % get the name of the current participant and session, and store it in the subjectdata
            subjectdata.nr   = num2str(f);
            subjectdata.ses  = num2str(s);
            fprintf('\n\nSubject: %s\n',subjectdata.nr)
            fprintf('Session: %s\n',subjectdata.ses)
            % locate the files of this participant
            subjectdata.dir_ses  = fullfile(curexperiment.datafolder_output, sprintf('%s',curexperiment.name,'_', subjectdata.nr));
            % load in the EOG data
            load([subjectdata.dir_ses filesep 'Session_' subjectdata.ses filesep subjectdata.nr '_Enc_ICA_EOG.mat'])
            % add EEGLAB to the path to use eeglab2fieldtrip
            addpath(fullfile(curexperiment.dirroot,'EEGLAB','eeglab2024_0','plugins','Fieldtrip-lite20201008','external','eeglab'))
            % save the averaged data
            eog_data = eeglab2fieldtrip(EEG,'timelock','none');
            eog_data.dimord = 'chan_time';
            % perform the time-frequency analysis
            cfg = [];
            cfg.method = 'mtmconvol';
            cfg.output = 'pow';
            eog_tfr = ft_freqanalysis(cfg,eog_data);
            % remove EEGLAB from the path
            rmpath(fullfile(curexperiment.dirroot,'EEGLAB','eeglab2024_0','plugins','Fieldtrip-lite20201008','external','eeglab'))
        end
    end
end