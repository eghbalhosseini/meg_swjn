%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 0: prepare the workspace 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
close all 
data_path='/Users/eghbalhosseini/MyData/meg_swjn';
subject_id='_017';
d= dir([data_path,strcat('/**/*',subject_id,'*_crunched.mat')]);
fprintf(' %d .mat files were found \n', length(d));
save_path='/Users/eghbalhosseini/MyData/meg_swjn/crunched/';
analysis_path='/Users/eghbalhosseini/MyData/meg_swjn/analysis/meg_swjn_Analysis_1_f_spect';
if ~exist(analysis_path)
    mkdir(analysis_path);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  STEP 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trial_condition=[];
% 
for k=1:length(d)
    %fprintf('adding %s from %s \n',d(k).name, strcat(d(k).folder,'/',d(k).name));
    subj=load(strcat(d(k).folder,'/',d(k).name));
    subj_id=fieldnames(subj);
    subj=subj.(subj_id{1});
    data=subj.data;
    info=subj.info;
    meg_chan=find(info.mag_chans + info.grad_chans);
    % step 1: compute mean across the word positions in each trial 
    %%%%%%%%%%%%%%%%%%% sentence 
    % convert the gamma into a tensor of shape : channel*trial*words
    sentence_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'S'),info.trial_type,'UniformOutput',false));    
    sentences=[data(sentence_trial_index)];
    stim_types=cellfun(@(x) x.stimuli_type,sentences,'uni',false);
    word_loc=cellfun(@(x) find(strcmp(x,'word')), stim_types, 'uni', false);
    signal_all_word=arrayfun(@(x) [sentences{x}.parsed{word_loc{x}}],1:length(word_loc),'uni',false);
    sentence_all_words=cellfun(@(x) x(meg_chan,:),signal_all_word','uni',false);
    % 
    nonwords_trial_index=~cellfun(@isempty,cellfun(@ (x) strfind(x,'N'),info.trial_type,'UniformOutput',false));    
    nonwords=[data(nonwords_trial_index)];
    stim_types=cellfun(@(x) x.stimuli_type,nonwords,'uni',false);
    word_loc=cellfun(@(x) find(strcmp(x,'nonword')), stim_types, 'uni', false);
    signal_all_word=arrayfun(@(x) [nonwords{x}.parsed{word_loc{x}}],1:length(word_loc),'uni',false);
    nonwords_all_words=cellfun(@(x) x(meg_chan,:),signal_all_word','uni',false);
end
%% spectral estimation 
Fs=info.sample_rate;
dt = 1/Fs  ; % sampling interval
fNQ = Fs/2 ; % Nyquist freq
N = length(x) ;
T = N*dt   ;
df = Fs/N  ;
n = 2^nextpow2(N);
f = Fs*(0:(n/2))/n;
fmax=300;
fmax_loc=find(f<=300);
sent_shapes=cell2mat(cellfun(@size,sentence_all_words,'uni',false));
nonw_shapes=cell2mat(cellfun(@size,nonwords_all_words,'uni',false));

num_rows=8;
num_columns=3;
total_plots=num_rows*num_columns;
p=0;
% 
close all 
for i=1:size(sentence_all_words{1},1)
    chan_sent=cell2mat(cellfun(@(x) x(i,1:min(sent_shapes)*[0;1]),sentence_all_words,'uni',false))';
    chan_nonw=cell2mat(cellfun(@(x) x(i,1:min(nonw_shapes)*[0;1]),nonwords_all_words,'uni',false))';
    pxx_sent=pmtm(chan_sent);
    pxx_nonw=pmtm(chan_nonw);
    fig=figure(fix((i-1)/total_plots)+1);
    set(gcf,'position',[1451,-285,957,1342]);
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    y=mean(pxx_sent,2);
    z=std(pxx_sent,[],2);
    bl = boundedline(f(fmax_loc(2:end)), y(fmax_loc(2:end)), z(fmax_loc(2:end)),'cmap',[.9,0,0],'alpha');
    bl.LineWidth=2;bl.DisplayName='sentences';
    hold on 
    y=mean(pxx_nonw,2);
    z=std(pxx_nonw,[],2);
    bl = boundedline(f(fmax_loc(2:end)), y(fmax_loc(2:end)), z(fmax_loc(2:end)),'cmap',[0,0,1],'alpha');
    bl.LineWidth=2;bl.DisplayName='nonwords';
    ax.XAxis.Scale='log';
    ax.XLim=[0,300];
    ax.XTick=[1,10,100];
    if ~mod(i,total_plots) | i==size(sentence_all_words{1},1)
        ax.XLabel.String='frequency (Hz)';
        ax.YLabel.String='power';
        ax.Children(2).Annotation.LegendInformation.IconDisplayStyle='off';
        ax.Children(4).Annotation.LegendInformation.IconDisplayStyle='off';
        legend('show')
        p=p+1;
        if ~exist(strcat(analysis_path,'/sub_',info.subject,'/'))
            mkdir(strcat(analysis_path,'/sub_',info.subject,'/'));
        end 
        print(fig, '-depsc', strcat(analysis_path,'/sub_',info.subject,'/','sub_',info.subject,'_spect_',num2str(p),'.eps')); 
        close(fig)
    end 
end 
%% compute spectrogram 
Fs=info.sample_rate;
dt = 1/Fs  ; % sampling interval
fNQ = Fs/2 ; % Nyquist freq
N = length(x) ;
T = N*dt   ;
df = Fs/N  ;
n = 2^nextpow2(N);
f = Fs*(0:(n/2))/n;
fmax=300;
fmax_loc=find(f<=300);
sent_shapes=cell2mat(cellfun(@size,sentence_all_words,'uni',false));
nonw_shapes=cell2mat(cellfun(@size,nonwords_all_words,'uni',false));

num_rows=8;
num_columns=3;
total_plots=num_rows*num_columns;
p=0;
% 
vir_color=plasma;
close all 
for i=1:size(sentence_all_words{1},1)
    chan_sent=(cellfun(@(x) x(i,1:min(sent_shapes)*[0;1]),sentence_all_words,'uni',false))';
    chan_nonw=(cellfun(@(x) x(i,1:min(nonw_shapes)*[0;1]),nonwords_all_words,'uni',false))';
    [sp,fp,tp] = cellfun(@(x) pspectrum(x,Fs,'spectrogram'),chan_sent,'uni',false);
    sp_all=cell2mat(permute(sp,[1,3,2]));
    sp_sent=mean(sp_all,3);
    % 
  
    fp=fp{1};
    fmax=300;
    tp=tp{1};
    fmax_loc=find((0<fp) .* (fp<=300));
    [sp,fp,tp] = cellfun(@(x) pspectrum(x,Fs,'spectrogram'),chan_nonw,'uni',false);
    
    sp_all=cell2mat(permute(sp,[1,3,2]));
    sp_nonw=mean(sp_all,3);
    fp=fp{1};
    fmax=300;
    tp=tp{1};
   
    fig=figure(fix((i-1)/total_plots)+1);
    set(gcf,'position',[1451,-285,957,1342]);
    ax=subplot(num_rows,num_columns,i-total_plots*fix((i-1)/total_plots));
    
    surf(tp,fp(fmax_loc),sp_nonw(fmax_loc,:),'LineStyle','none')
    text(max(tp),mean(fp(fmax_loc)),'Nonw','VerticalAlignment','bottom','FontSize',12,'HorizontalAlignment','left','Rotation',-90)
    hold on 
    gain_fac=1e3;
    surf(tp,gain_fac*(fp(fmax_loc)),sp_sent(fmax_loc,:),'LineStyle','none')
    text(max(tp),(mean(gain_fac*fp(fmax_loc))),'Sent','VerticalAlignment','bottom','FontSize',12,'HorizontalAlignment','left','Rotation',-90)
    view([0,90])
    ax.YAxis.Scale='log';
    ax.YTick=[1,10,100,gain_fac*[1,10,100]];
    ax.YTickLabel=arrayfun(@num2str,[1,10,100,[1,10,100]],'uni',false);
    axis tight
    grid off
    colormap(ax,vir_color)
    if ~mod(i,total_plots) | i==size(sentence_all_words{1},1)
    ax.XLabel.String='Time (s)';
    ax.YLabel.String='Frequency (Hz)';
        p=p+1;
        if ~exist(strcat(analysis_path,'/sub_',info.subject,'/'))
            mkdir(strcat(analysis_path,'/sub_',info.subject,'/'));
        end 
        print(fig, '-depsc', strcat(analysis_path,'/sub_',info.subject,'/','sub_',info.subject,'_spectgram_',num2str(p),'.eps')); 
        close(fig)
    end 
end 