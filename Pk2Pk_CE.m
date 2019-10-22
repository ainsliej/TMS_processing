%% extract TMS peak-to-peak measurements 
%% Corticospinal excitability blocks of current direction experiment 
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
        fileName=['P',num2str(i),'_S',num2str(s),'_',timept,'CE.mat'];
        load(fileName);
        try
        for muscle=1:3 %musles, obvs
                     
        data=D.data(:,:,muscle);
        
        %find all the instances of the desired state and filter the data on
        %this basis
        thesecol=logical(D.state==1)';
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
        
        %create a matrix of the MEP sizes for each muscle. FDI, APB, ADM. 
        %This should be a 15x3 matrix as there should be 15 CE pulses
        CE_MEP(:,muscle)=muscleMEP(length(muscleMEP)-14:length(muscleMEP));
        end
        
         %save the matrix with the MEP sizes for each timept, session and ptp            
        fileName=['P',num2str(i),'_S',num2str(s),'_CE',num2str(t),'.txt'] ;
        cd ~/../../Volumes/Ainslie_USB/VibData/PreProcessedData/IndividualStates;
        dlmwrite(fileName, CE_MEP ,'delimiter', ',', 'precision', 6);
        
        %filtering done here
        grubbsMEP=isoutlier(CE_MEP,1); %Grubbs test N.B. only on matlab_2018R
        noGrubbsMEP=CE_MEP;
        noGrubbsMEP(grubbsMEP==1)=NaN; 
        
        meanstateMEPs=mean(noGrubbsMEP,1,'omitnan');
        clear CE_MEP muscleMEP
        
        catch
              meanstateMEPs=[NaN,NaN,NaN]; 
                end 
        
        %Break the code, find out which stim condition
        if s==1
            stim=codebreak(ptpcount);
        else
            stim=codebreak(ptpcount)-1;
        end
        
        
       %now add these means into one matrix. These are arranged with ptps 
       %in separate rows, then the columns hierarchically stim type (sham, real),
       %then timepoint then muscle. 
       startpos=abs(stim)*12+(t-1)*3+1;
       endpos=abs(stim)*12+(t-1)*3+3;
        
       CEvals(ptpcount,startpos:endpos)=meanstateMEPs;
        end     
    end 
  disp(strcat('Analysis for ptp',num2str(i),' now complete. Total ptps analysed= ',num2str(ptpcount)))

    end 
 
  %% Now save this matrix
       
       cd ~/../../Volumes/Ainslie_USB/VibData/PreProcessedData;
       
     dlmwrite('MEP_CE.txt', CEvals ,'delimiter', ',', 'precision', 6);        
       
  BL_CEvals=[CEvals(:,4:6)-CEvals(:,1:3),CEvals(:,7:9)-CEvals(:,1:3),...
        CEvals(:,10:12)-CEvals(:,1:3),CEvals(:,16:18)-CEvals(:,13:15),...
     CEvals(:,19:21)-CEvals(:,13:15),CEvals(:,22:24)-CEvals(:,13:15)]; 
  
  dlmwrite('BL_CEvals.txt',  BL_CEvals ,'delimiter', ',', 'precision', 6); 
  
  %% Convert the matrix to longform for use in R
  
T1_longform_CE=[reshape(CEvals(:,1:3),[],1);reshape(CEvals(:,13:15),[],1)];
BL_longform_CE=reshape(BL_CEvals,[],1);

T1_ptp_CE=num2cell(repmat([5:30]',6,1));
T1_tDCS_CE=[repelem({'sham'},78,1);repelem({'real'},78,1)];
T1_muscle_CE=repmat([repelem({'FDI'},26,1);repelem({'APB'},26,1);repelem({'ADM'},26,1)],2,1);


BL_ptp_CE=num2cell(repmat([5:30]',18,1));
BL_tDCS_CE=repmat([repelem({'sham'},78,1);repelem({'real'},78,1)],3,1);
BL_muscle_CE=repmat([repelem({'FDI'},26,1);repelem({'APB'},26,1);repelem({'ADM'},26,1)],6,1);


T1_tableCE=table(T1_ptp_CE, T1_tDCS_CE, T1_muscle_CE,T1_longform_CE);
T1_tableCE.Properties.VariableNames = {'ptp','tDCS','muscle','DATA'};
BL_tableCE=table(BL_ptp_CE, BL_tDCS_CE, BL_muscle_CE,BL_longform_CE);
BL_tableCE.Properties.VariableNames = {'ptp','tDCS','muscle','DATA'};
 
writetable(T1_tableCE)
writetable(BL_tableCE)