%% extract TMS peak-to-peak measurements 
%% Sensory afferent inhibition blocks of current direction experiment 
%% Ainslie Johsntone
% N.B this is dependent on having a .txt file of the amplitude of EMG
% response from 1.15s - 1.30s with sampling 5000/s

%% Define some parameters
clear all
p1=3; %First ptp number
pn=4; %Last ptp number
samp=5000; %sampling per sec
prestart=1.172*samp; %Start of where we will look for precontractions 1.172s
preend=1.246*samp; %End precontraction window 1.247s
pulsestart=1.248*samp; %Start of where pulse artifact should be 1.2478s
pulseend=1.252*samp; %End 1.252s
MEPstart=1.266*samp; %Start of MEP 1.267s
MEPend=1.295*samp; %End 1.295s
MinPulse=0.1; %This should be the smallest possible size of pulse artifact

cd ~/../../Volumes/Ainslie_USB/VibData/; %Directory containing folder with extracted data

%% Loop around ptps, sessions, timepoints, states, and muscles

 for i=[p1:pn] %ptps
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
        
        %now add these means into separate matricies for each state. These
        %are arranged with ptps in separate rows, then the columns
        %hierarchically session, then timepoint then muscle.        
        startpos=(s-1)*12+(t-1)*3+1;
        endpos=(s-1)*12+(t-1)*3+3;
        
        if state==1
            state1vals(i,startpos:endpos)=meanstateMEPs;
        elseif state==2
            state2vals(i,startpos:endpos)=meanstateMEPs;
        end
           end
          
        end  
     
    end 
 end 
 
 %% Now save these matrices with the information for each state 
 
 cd ../PreProcessedData;
 dlmwrite('MEPnoSAI.txt', state1vals ,'delimiter', ',', 'precision', 6);      
 dlmwrite('MEPwithSAI.txt', state2vals ,'delimiter', ',', 'precision', 6);                 
 
 %calculate the SAI effect and save that too
 SAIeffect=state2vals./state1vals;
 dlmwrite('SAIeffect.txt', SAIeffect ,'delimiter', ',', 'precision', 6);      
 
    