%% extract TMS peak-to-peak measurements 
%% Vibration blocks of current direction experiment 
%% Ainslie Johsntone


p1=3;
pn=4;
prestart=1;
preend=486;
pulsestart=490;
pulseend=506;
MEPstart=586;
MEPend=726;

cd ~/../../Volumes/Ainslie_USB/VibData/;

 for i=[p1:pn]
    for s=1:2
        for t=1:4
            for state=1:6
               try 
                for muscle=1:3
   
        cd ~/../../Volumes/Ainslie_USB/VibData/ExtractedTxtData;
        fileName=['P',num2str(i),'_S',num2str(s),'_VIB',num2str(t),'_ch',num2str(muscle),'.txt'];
        data=load(fileName);
    
        thesecol=logical(data(1,:)==state);
        thisdata=data(2:end,thesecol);
        [Srow, Scol]=size(thisdata);
        
        for frame=1:Scol
           preconsize=max(thisdata(prestart:preend,frame))-min(thisdata(prestart:preend,frame));
           pulsesize=max(thisdata(pulsestart:pulseend,frame))-min(thisdata(pulsestart:pulseend,frame));
           MEPsize=max(thisdata(MEPstart:MEPend,frame))-min(thisdata(MEPstart:MEPend,frame));
        if preconsize<0.4 && pulsesize>0.1
            MEP=MEPsize;
        else 
            MEP=NaN;
        end
            muscleMEP(frame,:)=MEP;
        end
            
        stateMEP(:,muscle)=muscleMEP(length(muscleMEP)-11:length(muscleMEP));
                end 
                
                    
        fileName=['P',num2str(i),'_S',num2str(s),'_VIB',num2str(t),'_st',num2str(state),'.txt'] ;
        cd ~/../../Volumes/Ainslie_USB/VibData/PreProcessedData/IndividualStates;
        dlmwrite(fileName, stateMEP ,'delimiter', ',', 'precision', 6);
                       
        meanstateMEPs=mean(stateMEP,1,'omitnan');
        clear muscleMEP stateMEP 
        catch
              meanstateMEPs=[NaN,NaN,NaN]; 
        end 
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
 
cd ../PreProcessedData;
 dlmwrite('MEPvibNO.txt', state1vals ,'delimiter', ',', 'precision', 6);      
 dlmwrite('MEPvibADM.txt', state2vals ,'delimiter', ',', 'precision', 6);                 
 dlmwrite('MEPvibFDI.txt', state3vals ,'delimiter', ',', 'precision', 6);  
 dlmwrite('state4vals.txt', state4vals ,'delimiter', ',', 'precision', 6);      
 dlmwrite('state5vals.txt', state5vals ,'delimiter', ',', 'precision', 6);      
 dlmwrite('state6vals.txt', state6vals ,'delimiter', ',', 'precision', 6);
 
 SICIvibNO=state4vals./state1vals;
 SICIvibADM=state5vals./state2vals;
 SICIvibFDI=state6vals./state3vals;
 
 dlmwrite('SICIvibNO.txt', SICIvibNO ,'delimiter', ',', 'precision', 6);      
 dlmwrite('SICIvibADM.txt', SICIvibADM ,'delimiter', ',', 'precision', 6);                 
 dlmwrite('SICIvibFDI.txt', SICIvibFDI ,'delimiter', ',', 'precision', 6); 
    