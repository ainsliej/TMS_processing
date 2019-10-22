%% extract TMS peak-to-peak measurements 
%% Sensory afferent inhibition blocks of current direction experiment 
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
        fileName=['P',num2str(i),'_S',num2str(s),'_',timept,'SAI.mat'];
        load(fileName);
           
            for state=1:6 %states  
                try
                for muscle=1:3 %musles, obvs
                  
        data=D.data(:,:,muscle);
    
        %find all the instances of the desired state and filter the data on
        %this basis
        thesecol=logical(D.state==state)';
        pulsedata=D.data(:,thesecol,:);
        thisdata=data(:,thesecol);
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
        
        %NaN any frames where there is precontraction or no pulse   
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
        fileName=['P',num2str(i),'_S',num2str(s),'_SAI',num2str(t),'_st',num2str(state),'.txt'] ;
        cd ~/../../Volumes/Ainslie_USB/VibData/PreProcessedData/IndividualStates;
        dlmwrite(fileName, stateMEP ,'delimiter', ',', 'precision', 6);
        
        %filtering done here
        grubbsMEP=isoutlier(stateMEP,1); %Grubbs test N.B. only on matlab_2018R
        noGrubbsMEP=stateMEP;
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
        %hierarchically session, then timepoint then muscle.        
        startpos=abs(stim)*12+(t-1)*3+1;
        endpos=abs(stim)*12+(t-1)*3+3;
        
        if state==1
            state1vals(ptpcount,startpos:endpos)=meanstateMEPs;
        elseif state==2
            state2vals(ptpcount,startpos:endpos)=meanstateMEPs;
        end
           end
          
        end  
     
    end 
   disp(strcat('Analysis for ptp',num2str(i),' now complete. Total ptps analysed= ',num2str(ptpcount))) 
 end 
 
 %% Now save these matrices with the information
 
cd ~/../../Volumes/Ainslie_USB/VibData/PreProcessedData;
 dlmwrite('MEPnoSAI.txt', state1vals ,'delimiter', ',', 'precision', 6);      
 dlmwrite('MEPwithSAI.txt', state2vals ,'delimiter', ',', 'precision', 6);                 
 
 %calculate the SAI effect and save that too
 SAIeffect=state2vals./state1vals;
 dlmwrite('SAIeffect.txt', SAIeffect ,'delimiter', ',', 'precision', 6);      
 
 BL_SAIeffect=[SAIeffect(:,4:6)-SAIeffect(:,1:3),SAIeffect(:,7:9)-SAIeffect(:,1:3),...
        SAIeffect(:,10:12)-SAIeffect(:,1:3),SAIeffect(:,16:18)-SAIeffect(:,13:15),...
     SAIeffect(:,19:21)-SAIeffect(:,13:15),SAIeffect(:,22:24)-SAIeffect(:,13:15)]; 
 
 dlmwrite('BL_SAIeffect.txt',  BL_SAIeffect ,'delimiter', ',', 'precision', 6);  
 
%% Convert the matrix to longform datasets for use in R
T1_longform_SAI=[reshape(SAIeffect(:,1:3),[],1);reshape(SAIeffect(:,13:15),[],1)];
BL_longform_SAI=reshape(BL_SAIeffect,[],1);

T1_ptp_SAI=num2cell(repmat([5:30]',6,1));
T1_tDCS_SAI=[repelem({'sham'},78,1);repelem({'real'},78,1)];
T1_muscle_SAI=repmat([repelem({'FDI'},26,1);repelem({'APB'},26,1);repelem({'ADM'},26,1)],2,1);


BL_ptp_SAI=num2cell(repmat([5:30]',18,1));
BL_tDCS_SAI=repmat([repelem({'sham'},78,1);repelem({'real'},78,1)],3,1);
BL_muscle_SAI=repmat([repelem({'FDI'},26,1);repelem({'APB'},26,1);repelem({'ADM'},26,1)],6,1);


T1_tableSAI=table(T1_ptp_SAI, T1_tDCS_SAI, T1_muscle_SAI,T1_longform_SAI);
T1_tableSAI.Properties.VariableNames = {'ptp','tDCS','muscle','DATA'};
BL_tableSAI=table(BL_ptp_SAI, BL_tDCS_SAI, BL_muscle_SAI,BL_longform_SAI);
BL_tableSAI.Properties.VariableNames = {'ptp','tDCS','muscle','DATA'};

writetable(T1_tableSAI)
writetable(BL_tableSAI)