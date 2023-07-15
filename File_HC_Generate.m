function [STIMS,stim_structure]=File_HC_Generate(stim_structure,HW)

%this code specific for TonePipSweep
% Example Stim Structure
% stim_structure
%         left_check: 0
%        right_check: 1
%         both_check: 0
%        frame_check: 1
%     linescan_check: 0
%               type: 'TonePipSweep'
%      stim_protocol: [1x1 struct]
%       protocol_str: [1x134 char]
%
% stim_structure.stim_protocol
%              frame_dur: 0.2000
%                pip_dur: 0.1250
%     total_frames_lines: 25
%             stim_start: 11
%                ipi_val: 0.2500
%                freq_lo: 10000
%                freq_hi: 40000
%         num_freq_steps: 9
%               level_lo: 30
%               level_hi: 10
%          num_att_steps: 3
%        npips_per_train: 4
%             total_reps: 5


Fs=HW.Fso;STIMS.Fs=HW.Fso;
HC_new=stim_structure.stim_protocol.play_3_tokenHCT;% added by ann -22/6/22
% following parameters for all stimuli %%%%%%%%%%%%%%%%%%%%%%%
% defining two conditions for HC only and HC-T combo as when
switch HC_new
    case 0
        nfs=stim_structure.stim_protocol.num_freq_steps;
        nas=stim_structure.stim_protocol.num_att_steps;
        numhc=stim_structure.stim_protocol.numhc;
        fl=stim_structure.stim_protocol.freq_lo;fh=stim_structure.stim_protocol.freq_hi;
        ll=stim_structure.stim_protocol.level_lo;lh=stim_structure.stim_protocol.level_hi;
        pipd=stim_structure.stim_protocol.pip_dur;
        ipi=stim_structure.stim_protocol.ipi_val;
        npips=stim_structure.stim_protocol.npips_per_train;
        STIMS.duration=npips*ipi;
        STIMS.lines_per_set=nfs*nas;
        STIMS.rep_sets=stim_structure.stim_protocol.total_reps;
        if stim_structure.left_check==1
            STIMS.ATTENS.Left=reshape(repmat(linspace(ll,lh,nas),nfs,1),1,nfs*nas);
            STIMS.ATTENS.Right=[];
        elseif stim_structure.right_check==1
            STIMS.ATTENS.Right=reshape(repmat(linspace(ll,lh,nas),nfs,1),1,nfs*nas);
            STIMS.ATTENS.Left=[];
        elseif stim_structure.both_check==1
            STIMS.ATTENS.Left=reshape(repmat(linspace(ll,lh,nas),nfs,1),1,nfs*nas);
            STIMS.ATTENS.Right=reshape(repmat(linspace(ll,lh,nas),nfs,1),1,nfs*nas);
        end
        STIMS.hard=1; % when 1 returns actual stimuli in STIMS.waveform.line1,2,3 etc
        % when 0 stimuli will be generated on the fly
        % when -1 stimuli will be loaded from designated directory on
        % the fly
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        STIMS.rise_fall_time=min([.01 (stim_structure.stim_protocol.pip_dur)/10]);
        STIMS.freqs=repmat(logspace(log10(fl),log10(fh),nfs),1,nas);
        tp=[0:(1/Fs):(ipi-1/Fs)];tpip=[0:(1/Fs):(pipd-1/Fs)];
        maxy=0;
        for line_num=1:nfs*nas
            %line_num
            om=2*pi*STIMS.freqs(line_num);
            y=zeros(size(tpip));
            for nhcc=1:numhc
                %nhcc
                y=y+sin(nhcc*om*tpip);
            end
            %figure(2011)
            %plot(y)
            %y=-sin(om*tpip);
            y=aud_gate_on(y,STIMS.rise_fall_time*1000,Fs);
            y=aud_gate_off(y,STIMS.rise_fall_time*1000,Fs);
            y=[y zeros(1,length(tp)-length(y))];y=repmat(y,1,npips);
            maxy=max([maxy max(abs(y))]);
            eval(sprintf('STIMS.waveform.line%i=y'';',line_num))
        end
        for line_num=1:nfs*nas
            %line_num
            eval(sprintf('STIMS.waveform.line%i=(1/maxy)*STIMS.waveform.line%i;',line_num,line_num))
        end
        STIMS.norm_factor=maxy;
        %save STIMS STIMS
        %save STIMS STIMS
    case 1 % for HC -T- combo
        %int_tokn_int_ms=(stim_structure.stim_protocol.inter_tok_int_ms)/1000;% added by ann -22/6/22
       int_tokn_int_ms= [60 90 150 280]/1000;
        nas=stim_structure.stim_protocol.num_att_steps;% should be always one
        numhc=stim_structure.stim_protocol.numhc;
        atten=stim_structure.stim_protocol.level_lo;% this case should have only one attenuation
        lh=stim_structure.stim_protocol.level_hi;
        pipd=stim_structure.stim_protocol.pip_dur;
        STIMS.rise_fall_time=min([.01 (pipd)/10]);

        STIMS.lines_per_set=16;
        STIMS.rep_sets=stim_structure.stim_protocol.total_reps;
        
        if stim_structure.left_check==1
            STIMS.ATTENS.Left=atten*ones(1,STIMS.lines_per_set);
            STIMS.ATTENS.Right=[];
        elseif stim_structure.right_check==1
            STIMS.ATTENS.Right=atten*ones(1,STIMS.lines_per_set);
            disp(STIMS.ATTENS.Right)
            STIMS.ATTENS.Left=[];
        elseif stim_structure.both_check==1
            STIMS.ATTENS.Left=atten*ones(1,STIMS.lines_per_set);
            STIMS.ATTENS.Right=atten*ones(1,STIMS.lines_per_set);
        end
        STIMS.hard=1; % when 1 returns actual stimuli in STIMS.waveform.line1,2,3 etc
        % when 0 stimuli will be generated on the fly
        % when -1 stimuli will be loaded from designated directory on
        % the fl
        STIMS.fre=stim_structure.stim_protocol.freq_f0; % single input
        tpip=0:(1/Fs):(pipd-1/Fs);
        om=2*pi*STIMS.fre;
        y1=sin(om*tpip);y2=sin(2*om*tpip); y1y2=y1+y2*.7071;
        t1=y1/(max(y1y2)/sqrt(2));
        t2= y2/(max(y1y2)/sqrt(2));
        t1t2= y1y2/sqrt(2);
        % for normalisation 
         STIMS.norm_factor=[(max(y1y2)/sqrt(2)) (max(y1y2)/sqrt(2)) sqrt(2)] ;

        t1=aud_gate_on(t1,STIMS.rise_fall_time*1000,Fs);
        t1=aud_gate_off(t1,STIMS.rise_fall_time*1000,Fs);
        % for t2
        t2=aud_gate_on(t2,STIMS.rise_fall_time*1000,Fs);
        t2=aud_gate_off(t2,STIMS.rise_fall_time*1000,Fs);
        % for t1t2
        t1t2=aud_gate_on(t1t2,STIMS.rise_fall_time*1000,Fs);
        t1t2=aud_gate_off(t1t2,STIMS.rise_fall_time*1000,Fs);
        % for t1, t2, t1t2
        STIMS.rmsvalues=[rms(t1) rms(t2) rms(t1t2)];

        %         %%%%%
                         HC_cell=cell(16,6);
                        %save HC_rand_mat HC_rand_mat;% store it in PARAMS
                        curdir=pwd;
                        cd C:\EXPERIMENTS\CODES\PARAMS
                        load HC_rand_mat;
                        cd (curdir)
        %% for saving the fre and gap values in the HC_cell for 16 stimulud with varying gap
        for i=1:16
            for j=1:6
                if HC_rand_mat(i,j)== 1
                    HC_cell{i,j}=t1;
                elseif HC_rand_mat(i,j)==2
                    HC_cell{i,j}=t2;
                elseif HC_rand_mat(i,j)==3
                    HC_cell{i,j}=t1t2;
                elseif HC_rand_mat(i,j)==4
                    HC_cell{i,j}=zeros(1,fix(int_tokn_int_ms(1)*Fs));
                elseif HC_rand_mat(i,j)==5
                    HC_cell{i,j}=zeros(1,fix(int_tokn_int_ms(2)*Fs));
                elseif HC_rand_mat(i,j)==6
                    HC_cell{i,j}=zeros(1,fix(int_tokn_int_ms(3)*Fs));
                elseif HC_rand_mat(i,j)==7
                    HC_cell{i,j}=zeros(1,fix(int_tokn_int_ms(4)*Fs));
                end
            end
        end
     
        curdir=pwd;
        cd C:\EXPERIMENTS\CODES
        save HC_cell HC_cell;
        cd (curdir)
%%
        STIMS.stimcode_HC=zeros(STIMS.rep_sets,STIMS.lines_per_set);
     
        for rep=1:STIMS.rep_sets
            for line_num=1:16
                rowvl=randi([1,16]);
                y=HC_cell(rowvl,:);
                STIMS.stimcode_HC(rep,line_num)=rowvl;
                y=cell2mat(y);
                eval(sprintf('STIMS.waveform.line%i=y'';',line_num))
            end
        end
     
        disp(STIMS)


end

