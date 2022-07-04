function Data2=ReconstructBandpasssubsampling(Data)

%%%%%%%%%%%%%%test
%%%% figure; plot(-440:29559,interpft(ifft(fft(real(reconstructBandpasssubsampling(double(AScans(:,113))))).*conj(fft(flipud(CE(1:4:end)),3000))),30000),1:1:30000,interpft(ifft(fft(double(AScans2(:,113))).*conj(fft(CE(1:4:end),3000))),30000))

        downsamplingfactor=3;
        info.SampleRate=10e6/downsamplingfactor;
        flags.AScanReconstructionFreq=10e6;
        flags.expectedAScanLength=48;
        expectedAScanlengthDS=ceil(flags.expectedAScanLength/downsamplingfactor);
        
        %%%greater ascan than 3000!
        flags.expectedAScanLength=max(expectedAScanlengthDS*downsamplingfactor,size(Data,1)*downsamplingfactor);
        expectedAScanlengthDS=ceil(flags.expectedAScanLength/downsamplingfactor);
               
        s1=flags.AScanReconstructionFreq; %info.SampleRate*info.DownsamplingFactor/2; %20MHz to 10 MHz standard        
        f1=0:s1/(flags.expectedAScanLength-1):s1;
        f1_2=f1(1:1+length(f1)/2);
        s2=info.SampleRate;
        f2=0:s2/(size(Data,1)-1):s2; %if more precise assignment is required...
        f3=f2+info.SampleRate;
         f3_2=f3(1:1+length(f3)/2);
         
        Data2=zeros(flags.expectedAScanLength,size(Data,2));
        datafft=fftshift(fft(Data,expectedAScanlengthDS,1),1); %potential padding
        
        %% mirror-boundary freq
        datafft(expectedAScanlengthDS/2+1)=datafft(expectedAScanlengthDS/2)+1;
        datafft(expectedAScanlengthDS/2+1)=0        ;
        
        SR=info.SampleRate;
        SR_tol=info.SampleRate+diff(f1(1:2));
        SR_tol2=info.SampleRate-diff(f1(1:2));
         
        idx_l = find(f1>info.SampleRate/2 & f1 <= SR); %%% resolution of f1 as tolerance!1
        idx_r = find(f1>flags.AScanReconstructionFreq-SR & f1 <= flags.AScanReconstructionFreq-(SR/2));

        if isempty(idx_l) || isempty(idx_r) disp('downsampled data, you have to set flags.AScanReconstructionFreq higher'); return; end
        if size(datafft,1)/2 ~= length(idx_r) disp('error: Fourier resolution (->length) needs to be fixed by padding'); return; end

        % Data2([idx_l idx_r])=fftshift(fft(Data)); %([length(Data)/2:length(Data) 1 2:length(Data)/2+1])));
if 0
        for j=1:size(Data,2)
            Data2(idx_l(1:end-0)-1,j) = (datafft(1+end-(size(datafft,1)/2)-0:end,j));%(datafft(1:floor(size(datafft,1)/2)-0,j));
            Data2(idx_r(1:end-0)+1,j) = ((datafft(1:floor(size(datafft,1)/2)-0,j)));%(flipud(conj(datafft(1:(size(datafft,1)/2)-0,j))));
        end
else
        for j=1:size(Data,2)
            Data2(idx_l(1:end-0),j) = (datafft(1:floor(size(datafft,1)/2)-0,j));
            Data2(idx_r(1:end-0)+1,j) = (flipud(conj(datafft(1:(size(datafft,1)/2)-0,j))));
        end
end
        %datafft=fft(Data);
        %Data3=zeros(length(Data)*info.DownsamplingFactor/2,1);
        %Data3(find(f1<=s2 & f1>=s2/2))=datafft(find(f2<=s2 & f2>=s2/2));
        %Data3(find(s1/2+fliplr(f1)<=s1/2+s2 & s1/2+fliplr(f1)>=s1/2+s2/2))=datafft(find(f2>=0e6 & f2<=s2/2));

        %fixup & update
        Data2 = ifft(downsamplingfactor.*Data2,[],1);
end    