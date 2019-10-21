%% extract TMS peak-to-peak measurements 
%% Vibration blocks of current direction experiment 
%% Ainslie Johsntone
% N.B this is dependent on having a .txt file of the amplitude of EMG
% response from 1.15s - 1.30s with sampling 5000/s

%% Define some parameters
clear all
p1=5; %First ptp number
pn=30; %Last ptp number
samp=5000; %sampling per sec
prestart=1.172*samp; %Start of where we will look for precontractions 1.172s
preend=1.246*samp; %End precontraction window 1.247s
pulsestart=1.248*samp; %Start of where pulse artifact should be 1.2478s
pulseend=1.252*samp; %End 1.252s
MEPstart=1.266*samp; %Start of MEP 1.267s
MEPend=1.295*samp; %End 1.295s
MinPulse=0.1; %This should be the smallest possible size of pulse artifact
codebreak=[0,1,0,1,0,0,0,1,1,0,1,1,1,1,0,1,0,0,1,1,0,1,1,0,0,0]; %Is session 1 real stim?
ptpcount=0;
cd ~/../../Volumes/Ainslie_USB/VibData/; %Directory containing folder with extracted data

%% Loop around ptps, sessions, timepoints, states, and muscles

 for i=[p1:pn] %ptps
     ptpcount=ptpcount+1;
    for s=1:2 %sessions
        for t=1:4 %timepoints
            if t==1
              timept='Base';
            elseif t==2
              timept='During1';
            elseif t==3
              timept='During2'; 
            elseif t==4
              timept='Post';
            end          
        %open first file
        cd ~/../../Volumes/Ainslie_USB/VibData/;
        fileName=['P',num2str(i),'_S',num2str(s),'_',timept,'VIB.mat'];
        load(fileName);
           
            for state=1:6 %states  
                try
                for muscle=1:3 %musles, obvs
                  
        data=D.data(:,:,muscle);
        
        %find all the instances of the desired state and filter the data on
        %this basis
        thesecol=logical(D.state==state)';
        thisdata=data(:,thesecol);
        pulsedata=D.data(:,thesecol,:);
        [Srow, Scol]=size(thisdata);
        
        %get the mean root-mean-square for all the trials of a given
        %muscle, then calulate the value above which we will reject a trial
        RootMS=rms(data(prestart:preend,:));
        MaxPrecon=mean(RootMS)+2*std(RootMS);
        
        %go through each frame in the data file and get the rms of the
        %precontraction, pk2pk size of the largest pulse articact
        %and pk2pk MEP size, in the time windows defined above
        for frame=1:Scol
         preconsize=rms(thisdata(prestart:preend,frame));
         pulsesize=max(max(pulsedata(pulsestart:pulseend,frame,:))-min(pulsedata(pulsestart:pulseend,frame,:)));
         MEPsize=max(thisdata(MEPstart:MEPend,frame))-min(thisdata(MEPstart:MEPend,frame));
           
        %NaN any frames where there is precontraction over max or no pulse
        if preconsize<MaxPrecon && pulsesize>MinPulse 
            MEP=MEPsize;
        else 
            MEP=NaN;
        end
        
        %create a matrix of all the MEP sizes for this muslce in this state
        muscleMEP(frame,:)=MEP;
        end
        
        %for this state create a matrix of the MEP sizes for each muscle.
        %FDI, APB, ADM. This should be a 12x3 matrix as there should be 12
        %of each state in each protocol
        stateMEP(:,muscle)=muscleMEP(length(muscleMEP)-11:length(muscleMEP));
        end 
                
        %save the matrix with the MEP sizes for each state, at each timept, session and ptp            
        fileName=['P',num2str(i),'_S',num2str(s),'_VIB',num2str(t),'_st',num2str(state),'.txt'] ;
        cd ~/../../Volumes/Ainslie_USB/VibData/PreProcessedData/IndividualStates;
        dlmwrite(fileName, stateMEP ,'delimiter', ',', 'precision', 6);
        
        %calculate the mean for each muscle from the above.
        %filtering also done here
        initmeanMEP=mean(stateMEP,1,'omitnan');
        stdMEP=std(stateMEP, 1,'omitnan');
        noHighMEP=stateMEP;
        noHighMEP(stateMEP>initmeanMEP+2*stdMEP)=NaN; %remove MEP>2sd above mean
        noLowMEP=noHighMEP;
        noLowMEP(stateMEP<initmeanMEP-2*stdMEP)=NaN; %remove MEP<2sd below mean
        grubbsMEP=isoutlier(noLowMEP,1); %Grubbs test N.B. only on matlab_2018R
        noGrubbsMEP=noLowMEP;
        noGrubbsMEP(grubbsMEP==1)=NaN; 
        
        meanstateMEPs=mean(noGrubbsMEP,1,'omitnan');
        clear muscleMEP stateMEP 
        catch
              meanstateMEPs=[NaN,NaN,NaN]; 
                end 
        
        %Break the code, find out which stim condition
        if s==1
            stim=codebreak(ptpcount);
        else
            stim=codebreak(ptpcount)-1;
        end
                
        %now add these means into separate matricies for each state. These
        %are arranged with ptps in separate rows, then the columns
        %hierarchically stim type (sham, real), then timepoint then muscle. 
        startpos=abs(stim)*12+(t-1)*3+1;
        endpos=abs(stim)*12+(t-1)*3+3;        
        if state==1
            state1vals(ptpcount,startpos:endpos)=meanstateMEPs;
        elseif state==2
            state2vals(ptpcount,startpos:endpos)=meanstateMEPs;
        elseif state==3
            state3vals(ptpcount,startpos:endpos)=meanstateMEPs; 
        elseif state==4
            state4vals(ptpcount,startpos:endpos)=meanstateMEPs;
        elseif state==5
            state5vals(ptpcount,startpos:endpos)=meanstateMEPs;    
        elseif state==6
            state6vals(ptpcount,startpos:endpos)=meanstateMEPs;
        end     
            end
          
        end  
    end 
      disp(strcat('Analysis for ptp',num2str(i),' now complete. Total ptps analysed= ',num2str(ptpcount)))
 end 
 
 %% Now save these matrices with the information for each state 
 
cd ~/../../Volumes/Ainslie_USB/VibData/PreProcessedData;
 dlmwrite('MEPvibNO.txt', state1vals ,'delimiter', ',', 'precision', 6);      
 dlmwrite('MEPvibADM.txt', state3vals ,'delimiter', ',', 'precision', 6);                 
 dlmwrite('MEPvibFDI.txt', state5vals ,'delimiter', ',', 'precision', 6);  
 dlmwrite('state2valsVIB.txt', state2vals ,'delimiter', ',', 'precision', 6);      
 dlmwrite('state4valsVIB.txt', state4vals ,'delimiter', ',', 'precision', 6);      
 dlmwrite('state6valsVIB.txt', state6vals ,'delimiter', ',', 'precision', 6);
 
 %calculate the SICI measures and save them too
 SICIvibNO=state2vals./state1vals;
 SICIvibADM=state4vals./state3vals;
 SICIvibFDI=state6vals./state5vals;
 
 dlmwrite('SICIvibNO.txt', SICIvibNO ,'delimiter', ',', 'precision', 6);      
 dlmwrite('SICIvibADM.txt', SICIvibADM ,'delimiter', ',', 'precision', 6);                 
 dlmwrite('SICIvibFDI.txt', SICIvibFDI ,'delimiter', ',', 'precision', 6); 
 
 
 %% Save the change from BL values 
 
 BL_MEPvibNO=[state1vals(:,4:6)-state1vals(:,1:3),state1vals(:,7:9)-state1vals(:,1:3),...
     state1vals(:,10:12)-state1vals(:,1:3),state1vals(:,16:18)-state1vals(:,13:15),...
     state1vals(:,19:21)-state1vals(:,13:15),state1vals(:,22:24)-state1vals(:,13:15)];
 BL_MEPvibADM=[state3vals(:,4:6)-state3vals(:,1:3),state3vals(:,7:9)-state3vals(:,1:3),...
     state3vals(:,10:12)-state3vals(:,1:3),state3vals(:,16:18)-state3vals(:,13:15),...
     state3vals(:,19:21)-state3vals(:,13:15),state3vals(:,22:24)-state3vals(:,13:15)]; 
 BL_MEPvibFDI=[state5vals(:,4:6)-state5vals(:,1:3),state5vals(:,7:9)-state5vals(:,1:3),...
     state5vals(:,10:12)-state5vals(:,1:3),state5vals(:,16:18)-state5vals(:,13:15),...
     state5vals(:,19:21)-state5vals(:,13:15),state5vals(:,22:24)-state5vals(:,13:15)]; 
 
 BL_SICIvibNO=[SICIvibNO(:,4:6)-SICIvibNO(:,1:3),SICIvibNO(:,7:9)-SICIvibNO(:,1:3),...
     SICIvibNO(:,10:12)-SICIvibNO(:,1:3),SICIvibNO(:,16:18)-SICIvibNO(:,13:15),...
     SICIvibNO(:,19:21)-SICIvibNO(:,13:15),SICIvibNO(:,22:24)-SICIvibNO(:,13:15)]; 
 BL_SICIvibADM=[SICIvibADM(:,4:6)-SICIvibADM(:,1:3),SICIvibADM(:,7:9)-SICIvibADM(:,1:3),...
     SICIvibADM(:,10:12)-SICIvibADM(:,1:3),SICIvibADM(:,16:18)-SICIvibADM(:,13:15),...
     SICIvibADM(:,19:21)-SICIvibADM(:,13:15),SICIvibADM(:,22:24)-SICIvibADM(:,13:15)];  
 BL_SICIvibFDI=[SICIvibFDI(:,4:6)-SICIvibFDI(:,1:3),SICIvibFDI(:,7:9)-SICIvibFDI(:,1:3),...
     SICIvibFDI(:,10:12)-SICIvibFDI(:,1:3),SICIvibFDI(:,16:18)-SICIvibFDI(:,13:15),...
     SICIvibFDI(:,19:21)-SICIvibFDI(:,13:15),SICIvibFDI(:,22:24)-SICIvibFDI(:,13:15)];  
 
 dlmwrite('BL_MEPvibNO.txt',  BL_MEPvibNO ,'delimiter', ',', 'precision', 6);      
 dlmwrite('BL_MEPvibADM.txt',  BL_MEPvibADM ,'delimiter', ',', 'precision', 6);                 
 dlmwrite('BL_MEPvibFDI.txt',  BL_MEPvibFDI ,'delimiter', ',', 'precision', 6);  
  
 dlmwrite('BL_SICIvibNO.txt', BL_SICIvibNO ,'delimiter', ',', 'precision', 6);      
 dlmwrite('BL_SICIvibADM.txt', BL_SICIvibADM ,'delimiter', ',', 'precision', 6);                 
 dlmwrite('BL_SICIvibFDI.txt', BL_SICIvibFDI ,'delimiter', ',', 'precision', 6); 
 
 %% Convert the matricies to longform datasets for use in R
 
 %Data from T1 timepoint
 
T1_longform_SP=[reshape(state1vals(:,1:3),[],1);reshape(state1vals(:,13:15),[],1);...
     reshape(state3vals(:,1:3),[],1);reshape(state3vals(:,13:15),[],1);...
     reshape(state5vals(:,1:3),[],1);reshape(state5vals(:,13:15),[],1)];
 
T1_longform_SICI=[reshape(SICIvibNO(:,1:3),[],1);reshape(SICIvibNO(:,13:15),[],1);...
     reshape(SICIvibADM(:,1:3),[],1);reshape(SICIvibADM(:,13:15),[],1);...
     reshape(SICIvibFDI(:,1:3),[],1);reshape(SICIvibFDI(:,13:15),[],1)];
 
T1_ptp_VIB=num2cell(repmat([5:30]',18,1));
T1_tDCS_VIB=repmat([repelem({'sham'},78,1);repelem({'real'},78,1)],3,1);
T1_muscle_VIB=repmat([repelem({'FDI'},26,1);repelem({'APB'},26,1);repelem({'ADM'},26,1)],6,1);
T1_vib_VIB=repelem({'vibNO';'vibADM';'vibFDI'},156,1);

T1_SP_tableVIB=table(T1_ptp_VIB, T1_tDCS_VIB, T1_vib_VIB, T1_muscle_VIB,T1_longform_SP);
T1_SP_tableVIB.Properties.VariableNames = {'ptp','tDCS','vibCond','muscle','DATA'};
T1_SICI_tableVIB=table(T1_ptp_VIB, T1_tDCS_VIB, T1_vib_VIB, T1_muscle_VIB,T1_longform_SICI);
T1_SICI_tableVIB.Properties.VariableNames = {'ptp','tDCS','vibCond','muscle','DATA'};

writetable(T1_SP_tableVIB)
writetable(T1_SICI_tableVIB)

 %Column layout: Ptp, VibCondition, tDCS_stim, TimePT, Muscle, VibCondition, Data 
 
BL_longform_SP=[reshape(BL_MEPvibNO,[],1);reshape(BL_MEPvibADM,[],1);reshape(BL_MEPvibFDI,[],1)];
BL_longform_SICI=[reshape(BL_SICIvibNO,[],1);reshape(BL_SICIvibADM,[],1);reshape(BL_SICIvibFDI,[],1)];

BL_ptp_VIB=num2cell(repmat([5:30]',54,1));
BL_tDCS_VIB=repmat([repelem({'sham'},234,1);repelem({'real'},234,1)],3,1);
BL_muscle_VIB=repmat([repelem({'FDI'},26,1);repelem({'APB'},26,1);repelem({'ADM'},26,1)],18,1);
BL_vib_VIB=repelem({'vibNO';'vibADM';'vibFDI'},468,1);
 
BL_SP_tableVIB=table(BL_ptp_VIB, BL_tDCS_VIB, BL_vib_VIB, BL_muscle_VIB,BL_longform_SP);
BL_SP_tableVIB.Properties.VariableNames = {'ptp','tDCS','vibCond','muscle','DATA'};
BL_SICI_tableVIB=table(BL_ptp_VIB, BL_tDCS_VIB, BL_vib_VIB, BL_muscle_VIB,BL_longform_SICI);
BL_SICI_tableVIB.Properties.VariableNames = {'ptp','tDCS','vibCond','muscle','DATA'};

writetable(BL_SP_tableVIB)
writetable(BL_SICI_tableVIB)