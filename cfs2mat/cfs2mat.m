% CFS2MAT
% =======
%
% Written by Gilad Jacobson, 9 May 2016
%
% Script for translating CFS files generated by CED's software to Matlab files
% 
% Prerequisites: Install the matcfs64c library available on the CED
% website (developed by J.G. Colebatch):
% http://ced.co.uk/downloads/contributed#MatlabCFS
% 
% Make sure to extract the files from the ZIP file into a folder and add
% that folder to the Matlab path.
% 
% To run:
% - Type cfs2mat in the Matlab prompt
% - Go to the wanted library and mark CFS files to be translated (multiple
%   files possible)
% 
% Results are automatically saved into .mat files with the same name as
% the .cfs files.
% Each .mat file contains a single structure variable named 'D'
% D.param: different parameters of the acquisition
% D.data:  the data, organised as a nPoints x nTrials x nChannels

% CFS parameters


LSTR = 7;
READ = 0;
INT4 = 4;  % **
DSVAR = 1; % **
nchan=1;


% choose files to convert

[fName, cfsDir] = uigetfile({'*.cfs','CFS files (*.cfs)'}, ...
    'Choose a file to load','multiselect','on');

if numel(fName) == 1 & all(fName == 0),
    error('No valid file selected');
end

if isstr(fName)
    fName = {fName};
end

% loop over all chosen files

for fCount = 1:length(fName)
    try
            fullfilename = [cfsDir fName{fCount}];
            fHandle = matcfs64c('cfsOpenFile',fullfilename,READ,0);
            if (fHandle < 0)
                error(['File opening error: ' int2str(fHandle)]);
            end
            D = struct; % initialise Matlab output structure
            
            % read file parameters
            [D.param.fTime,D.param.fDate,D.param.fComment] = matcfs64c('cfsGetGenInfo',fHandle);
            [D.param.channels,fileVars,DSVars,D.param.dataSections] = matcfs64c('cfsGetFileInfo',fHandle);
            stateDS = zeros(D.param.dataSections,1);  % ** allocate space
            % loop over all trials (data sections)
            for dsCount = 1:D.param.dataSections
%         dsCount
              %  flagSet = matcfs64c('cfsDSFlags', fHandle, dsCount, READ);  % setit = 0 to read    % ** not used
                stateDS(dsCount)=matcfs64c('cfsGetVarVal', fHandle, 1, DSVAR, dsCount, INT4);  % ** read state into vector
                dSbyteSize = matcfs64c('cfsGetDSSize',fHandle,dsCount);

                % loop over all channels

                for chCount = 1:nchan
%             chCount
                    % read trace parameters
                    [D.param.tOffset(chCount),points, ...
                        D.param.yScale(chCount), ...
                        D.param.yOffset(chCount), ...
                        D.param.xScale(chCount), ...
                        D.param.xOffset(chCount)] = ...
                        matcfs64c('cfsGetDSChan',fHandle,chCount-1,dsCount);
                    [D.param.channelName{chCount}, ...
                        D.param.yUnits{chCount},D.param.xUnits{chCount}, ...
                        dataType,dataKind,spacing,other] = ...
                        matcfs64c('cfsGetFileChan',fHandle,chCount-1);

                    % read actual data
                    D.data(:,dsCount,chCount) = ...
                        matcfs64c('cfsGetChanData',fHandle, ...
                        chCount-1,dsCount,0,points,dataType);
                    D.data(:,dsCount,chCount) = ...
                        (D.data(:,dsCount,chCount) * D.param.yScale(chCount)) + ...
                        D.param.yOffset(chCount);
                end % for chCount
            end  % for dsCount
           ret = matcfs64c('cfsCloseFile',fHandle); % close the CFS file
           D.state = stateDS;
           outName = [cfsDir fName{fCount}(1:end-4)];
           save(outName,'D') % save the Matlab structure
           disp(strcat('finished -',outName))
    catch
        disp(strcat('failed on -',fullfilename))
    end
end % for fCount