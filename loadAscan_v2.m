function [gain,AScan]= loadAScan_v2(sID,se,rID,re,mp,pathExperiment,varargin)
%
% Loads and returns Ascan for USCT v2
% If AScan could be loaded ascan is returned (without any preprocessing, i.e. unit16), but as (N,1)
%
% sID - Sender ID (TAS Nr)
% se - Sender Element
% rID - Receiver ID
% re - Receiver Element
% mp - aperture position id
% pathExperiment - path to experiment directory
% varargin - '': clear buffer,flags: update to current flags
% NOTE: do NOT use 'path' as identifier! It's an environment variable.
%
% last changed: $Date: 2010-07-29 $
% revision: $Rev: 2772 $


persistent fE AScanList TASMap receiverMap gainList timeStampList lastSID ...
	errors_occurred counter lastSE lastMP %List with ascans, only needed if new data format is used


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Active Clear needed between different Runs of Program!!!!
if isempty(varargin)  
    
    %deprecated CASE
    flags = getGlobalFlags;
else    
    %recreate flags from VARARGIN
    if isstruct(varargin{1})
        flags = varargin{1};
              
        if flags.useWaveletPulseDetection == 1
           disp('WARNING: use Wavelet pulse detectin in loadAScan not implemented vor USCT v2 yet!!!!!!!!!!!!');              
           return 
        end      
    % only clear
    else
       	AScanList=[]; TASMap = []; receiverMap = []; gainList =[]; timeStampList =[];lastSID=[]; lastSE=[]; lastMP =[];  fE=[];  AScan=[]; gain=[];
        disp('Cleared LoadAScan Buffers.');
        return
    end
end
%%%%%%%%%%%%%%%%%%%%%%END ACTIVE CLEAR


%checks for succesive of loadASCAN calls in ONE Run of Whole Program
AScan = []; gain = [];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%first enter of loadAScan
if isempty(fE)
    
    try eval('flags.FileFormat;'); catch
        flags.FileFormat= -1;  disp('Filled up missing flags!');end

    try eval('flags.OS_UNIX;'); catch
        flags.OS_UNIX = isunix;  disp('Filled up missing flags!');end

    try eval('flags.useWaveletPulseDetection;'); catch
       flags.useWaveletPulseDetection = 0;  disp('Filled up missing flags!');end

    try eval('flags.offsetElectronic;'); catch
        flags.offsetElectronic = 9; disp('Filled up missing flags!'); end

    try eval('flags.useWaveletPulseDetection;'); catch
        flags.useWaveletPulseDetection = 0; disp('Filled up missing flags!'); end

        if flags.useWaveletPulseDetection == 1       
             disp('WARNING: use Wavelet pulse detectin in loadAScan not implemented vor USCT v2 yet!!!!!!!!!!!!');              
             return
        end

    try eval('flags.expectedAScanLength;'); catch
        flags.expectedAScanLength = 3000; disp('Filled up missing flags!'); end

    try eval('flags.expectedAScanSampFreq;'); catch
        flags.expectedAScanSampFreq = 10000000; disp('Filled up missing flags!');end

    %assign('flags',caller);
    
    %msglen = 0; % msgbar does not delete previous output, because there wasn't any p.o.
    
    if ~exist(pathExperiment, 'dir')
        error(sprintf('target /source directory "%s" does not exist', pathExperiment));
    end
   
    fE=1; %first enter done!
    lastSID = -1; % needed for comparision that is not empty
    lastSE = -1;
    lastMP = -1;

	counter = 0;
	errors_occurred = 0;
   end
%%%%%%%%%%%%%%%%%%%END FIRST ENTER



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MAIN enter of LOADASCAN

if(isempty(AScanList) | sID~=lastSID | se ~= lastSE | mp~= lastMP) %new emitter needs to be read

	if isempty(errors_occurred)
		errors_occurred = 0;
		counter = 0;
	end
    counter = counter +1;
    subpath=sprintf('\\TAS%03d\\TASRotation%02d\\Emitter%02d.mat',sID,mp,se);
    if isunix subpath(strfind(subpath,'\'))='/';  end
    datapath=[pathExperiment,subpath];
	try load(datapath); catch
		errors_occurred = errors_occurred + 1;
		if (mod(errors_occurred, 10*ceil(sqrt(errors_occurred))) == 0) ...
			|| errors_occurred < 13
		%if (isprime(errors_occurred) == 0) ...
			%|| errors_occurred < 13
			disp(sprintf('Error loading a_scan (%3d times or %5.2f%% of a_scans / success for %4d). Details: sID %i mp %i se %i path:%s                          ', errors_occurred, floor(errors_occurred*10000/counter)/100, counter - errors_occurred, sID, mp, se, datapath));
		end
		return
    end
    if size(AScans, 1) == size(receiverIndices, 1)  % DIRTY WORKAROUND:
                                                % In case A scans are 
                                                % saved in wrong order,
                                                % transpose them. (Occurs
                                                % only in early version of
                                                % some simulations.
        AScans = AScans';
    end
    AScanList = AScans;
    TASMap = TASIndices;
    receiverMap = receiverIndices;
    if(exist('Amplification','var'))
        gainList = Amplification;
    end
    if(exist('TimeStamp','var'))
        timeStampList = TimeStamp;
    end
    lastSID = sID;
    lastSE = se;
    lastMP = mp;   
    
end
    
currentIndex = find((TASMap==rID)&(receiverMap==re));
if(isempty(currentIndex)|currentIndex>size(AScanList,2))
    %disp(sprintf('WARNING: No Ascan found for se %i sID %i re %i rID %i mp %i',se,sID,re,rID,mp));
    return
end

if(size(AScanList,1)==3001) %to be compatible to early USCT II measurements
   AScan = AScanList(2:end,currentIndex); 
   vgaControl = bitand(AScanList(1,currentIndex),2^12-1);
   TASDamping = 1; %?
   gain = CompleteGain(vgaControl,TASDamping);
    
else
    AScan = AScanList(:,currentIndex); 
    if(isempty(gainList))
        gain = 1;
    else
        gain = gainList(currentIndex,1);
    end
end



%%%%%%%%%%%%%%%%%%subfunction
function flags = getGlobalFlags()
global flags
%disp('Warning: Flags from global context, will be deprecated, use parameters!');
