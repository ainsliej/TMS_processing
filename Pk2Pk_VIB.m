%% extract TMS peak-to-peak measurements 
%% Vibration blocks of current direction experiment 
%% Ainslie Johsntone
% N.B this is dependent on having a .txt file of the amplitude of EMG
% response from 1.15s - 1.30s with sampling 5000/s

%% Define some parameters
clear all
p1=3; %First ptp number
pn=4; %Last ptp number
prestart=111; %Start of where we will look for precontractions 111=1.172s
preend=486; %End precontraction window 486=1.247s
pulsestart=490; %Start of where pulse artifact should be 490=1.2478s
pulseend=511; %End 511=1.252s
MEPstart=586; %Start of MEP 586=1.267s
MEPend=726; %End 726=1.295s
MaxPrecon=0.4; %Exclude trials where precon over this amplitude
MinPulse=0.2; %This should be the smallest possible size of pulse artifact

cd ~/../../Volumes/Ainslie_USB/VibData/; %Directory containing folder with extracted data

%% Loop around ptps, sessions, timepoints, states, and muscles

 for i=[p1:pn] %ptps
    for s=1:2 %sessions
        for t=1:4 %timepoints
            for state=1:6 %states
               try 
                for muscle=1:3 %musles, obvs
   
        %open first file
        cd ~/../../Volumes/Ainslie_USB/VibData/ExtractedTxtData;
        fileName=['P',num2str(i),'_S',num2str(s),'_VIB',num2str(t),'_ch',num2str(muscle),'.txt'];
        data=load(fileName);
    
        %find all the instances of the desired state and filter the data on
        %this basis
        thesecol=logical(data(1,:)==state);
        thisdata=data(2:end,thesecol);
        [Srow, Scol]=size(thisdata);
        
        %go through each frame in the data file and get the pk2pk amplitude of
        %precontraction, pulse and MEP, in the time windows defined above
        for frame=1:Scol
           preconsize=max(thisdata(prestart:preend,frame))-min(thisdata(prestart:preend,frame));
           pulsesize=max(thisdata(pulsestart:pulseend,frame))-min(thisdata(pulsestart:pulseend,frame));
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
        
        %now add these means into separate matricies for each state. These
        %are arranged with ptps in separate rows, then the columns
        %hierarchically session, then timepoint then muscle. 
        startpos=(s-1)*12+(t-1)*3+1;
        endpos=(s-1)*12+(t-1)*3+3;        
        if state==1
            state1vals(i,startpos:endpos)=meanstateMEPs;
        elseif state==2
            state2vals(i,startpos:endpos)=meanstateMEPs;
        elseif state==3
            state3vals(i,startpos:endpos)=meanstateMEPs; 
        elseif state==4
            state4vals(i,startpos:endpos)=meanstateMEPs;
        elseif state==5
            state5vals(i,startpos:endpos)=meanstateMEPs;    
        elseif state==6
            state6vals(i,startpos:endpos)=meanstateMEPs;
        end     
            end
          
        end  
     
    end 
 end 
 
 %% Now save these matrices with the information for each state 
 
cd ../PreProcessedData;
 dlmwrite('MEPvibNO.txt', state1vals ,'delimiter', ',', 'precision', 6);      
 dlmwrite('MEPvibADM.txt', state2vals ,'delimiter', ',', 'precision', 6);                 
 dlmwrite('MEPvibFDI.txt', state3vals ,'delimiter', ',', 'precision', 6);  
 dlmwrite('state4vals.txt', state4vals ,'delimiter', ',', 'precision', 6);      
 dlmwrite('state5vals.txt', state5vals ,'delimiter', ',', 'precision', 6);      
 dlmwrite('state6vals.txt', state6vals ,'delimiter', ',', 'precision', 6);
 
 %calculate the SICI measures and save them too
 SICIvibNO=state4vals./state1vals;
 SICIvibADM=state5vals./state2vals;
 SICIvibFDI=state6vals./state3vals;
 
 dlmwrite('SICIvibNO.txt', SICIvibNO ,'delimiter', ',', 'precision', 6);      
 dlmwrite('SICIvibADM.txt', SICIvibADM ,'delimiter', ',', 'precision', 6);                 
 dlmwrite('SICIvibFDI.txt', SICIvibFDI ,'delimiter', ',', 'precision', 6); 
    