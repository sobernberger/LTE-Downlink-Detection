close all; clc; clear;

file = fopen('/mnt/local_data/sobernberger/bfconversion/ts_decim_full_2.sigmf-data');
in = fread(file, 'float32');

sr = 2e6;

waveform = in(1:2:end) + 1j*in(2:2:end);

%waveform = waveform';

if (isempty(waveform))
    fprintf('\nReceived signal must not be empty.\n');
    return;
end

if (~exist('channelFigure','var') || ~isvalid(channelFigure))
    channelFigure = figure('Visible','off');
end

scat = comm.ConstellationDiagram;
scat.ShowGrid = true;
scat.Position = [1350 100 400 400];
scat.XLimits = [-3 3];
scat.YLimits = [-3 3];

sA = dsp.SpectrumAnalyzer();
sA.Position = [0 500 320 200];
sA.SampleRate = sr;
corrplt = dsp.ArrayPlot();
corrplt.XLabel = 'Timing Offset';
corrplt.YLabel = 'Correlation';
corrplt.ShowGrid = true;
corrplt.Position = [0 1000-410 320 200];
corrplt.PlotType = 'Line';


WaveformSize = whos('waveform');
fprintf('Decoding of provided LTE Waveform starts. Waveform size is %4.2fGb\n', WaveformSize.bytes/1024^3)


enb = struct;           
enb.NDLRB = 6;     
enb.CyclicPrefix = 'Normal';
ofdmInfo = lteOFDMInfo(enb);
sep = repmat('*',1,50);
NumSamplesPerFrame = 10e-3*sr;
ApproxNum2Frames = floor(length(waveform)/NumSamplesPerFrame);
fprintf('\nPlotting received signal spectrum...\n');
sA(waveform);



cec.PilotAverage = 'UserDefined';    
cec.FreqWindow = 13;                 
cec.TimeWindow = 9;                
cec.InterpType = 'cubic';          
cec.InterpWindow = 'Centered';     
cec.InterpWinSize = 1;       
MIBDecoded = {};
CellIds = double.empty(0, 3);
MIBFLAG = 0;
RSRQdB = [];
RSRPdBm = [];
RSSIdBm = [];
x = 1:(floor(ApproxNum2Frames/2)-1);

for frameNum = 1:(floor(ApproxNum2Frames/2)-1)
    CellIdOk = 1;
    MIBok = 1;
    enb.NDLRB = 6;
    if (sr~=ofdmInfo.SamplingRate)
        if (sr < ofdmInfo.SamplingRate)
            warning('The received signal sampling rate (%0.3fMS/s) is lower than the desired sampling rate for cell search / MIB decoding (%0.3fMs/s); cell search / MIB decoding may fail.',sr/1e6,ofdmInfo.SamplingRate/1e6);
        end
        fprintf('\nResampling from %0.3fMs/s to %0.3fMS/s for cell search / MIB decoding for Sample Block nr. %d (2 Frames worth of Samples)...\n',sr/1e6,ofdmInfo.SamplingRate/1e6, frameNum);
    else
        fprintf('\nResampling not required; received signal is at desired sampling rate for cell search / MIB decoding (%0.3fMs/s).\n',sr/1e6);
    end
    sampIdxStart = (frameNum-1)*NumSamplesPerFrame*2 + 1;
    Samples = waveform(sampIdxStart:sampIdxStart+NumSamplesPerFrame*2);

    
   

    % Downsample Samples
    nSamples = ceil(ofdmInfo.SamplingRate/round(sr)*size(Samples,1));
    nRxAnts = size(Samples, 2);
    downsampled = zeros(nSamples, nRxAnts);
    for i=1:nRxAnts
         downsampled(:,i) = resample(Samples(:,i), ofdmInfo.SamplingRate, round(sr));
    end
    

%% CELL SEARCH
    fprintf('\nPerforming cell search for Sample Block nr. %d: \n', frameNum);
    duplexModes = {'TDD' 'FDD'};
    cyclicPrefixes = {'Normal' 'Extended'};
    
   
    %Cell Search for all CyclicPrefixes and all Duplex Modes
    searchalg.MaxCellCount = 1;
    searchalg.SSSDetection = 'PostFFT';
    peakMax = -Inf;
    %searchalg.CellIDs = 481;
    
    
    for duplexMode = duplexModes
        for cyclicPrefix = cyclicPrefixes
            enb.DuplexMode = duplexMode{1};
            enb.CyclicPrefix = cyclicPrefix{1};
            [enb.NCellID, offset, peak] = lteCellSearch(enb, downsampled, searchalg);
            enb.NCellID = enb.NCellID(1);
            offset = offset(1);
            peak = peak(1);
            if (peak>peakMax)
                enbMax = enb;
                offsetMax = offset;
                peakMax = peak;
            end
        end
    end

    enb = enbMax;
    offset = offsetMax;

    if ismember(enb.NCellID, CellIds(:,1))
        [~, idx] = ismember(enb.NCellID, CellIds(:,1));
        CellIds(idx, 2) = CellIds(idx, 2) + peakMax;
        CellIds(idx, 3) = CellIds(idx, 3) + 1;
        


    else 
       CellIds(end+1, 1) = enb.NCellID;
       CellIds(end, 2) = peakMax;
       CellIds(end, 3) = 1;
    end


    corr = cell(1,3);
    idGroup = floor(enbMax.NCellID/3);
    for i = 0:2
        enb.NCellID = idGroup*3 + mod(enbMax.NCellID + i,3);
        [~,corr{i+1}] = lteDLFrameOffset(enb, downsampled);
        corr{i+1} = sum(corr{i+1},2);
    end
    threshold = 1.3 * max([corr{2}; corr{3}]);
    if (max(corr{1})<threshold)    
        warning('Synchronization signal correlation was weak; detected CellID may be incorrect.');
        CellIdOk = 0;
    end
    
    % Return to originally detected cell identity
    enb.NCellID = enbMax.NCellID;

    corrplt.YLimits = [0 max([corr{1}; threshold])*1.1];
    corrplt([corr{1} threshold*ones(size(corr{1}))]);
     
    if offset < 0
        fprintf('Offset value of %d is invalid. Skipping to next dataset.', offset);
        continue
    end
    % Perform timing synchronization
    fprintf('Timing offset to frame start: %d samples\n',offset);
    downsampled = downsampled(1+offset:end,:); 
    enb.NSubframe = 0;
    
    fprintf('Cell-wide settings after cell search:\n');
    disp(enb);
    

    fprintf('\nPerforming frequency offset estimation...\n');

    if (strcmpi(enb.DuplexMode,'TDD'))
        enb.TDDConfig = 0;
        enb.SSC = 0;
    end
    delta_f = lteFrequencyOffset(enb, downsampled);
    fprintf('Frequency offset: %0.3fHz\n',delta_f);
    downsampled = lteFrequencyCorrect(enb, downsampled, delta_f); 

    % Assume 4 cell-specific reference signals for initial decoding attempt
    enb.CellRefP = 4;   
                        
    fprintf('Performing OFDM demodulation...\n\n');
    
    griddims = lteResourceGridSize(enb); % Resource grid dimensions
    L = griddims(2);                     % Number of OFDM symbols in a subframe 
    % OFDM demodulate signal 
    rxgrid = lteOFDMDemodulate(enb, downsampled);    
    
%     enbcp = enb;
%     enbcp.CellRefP = 1;
%     enbcp.NDLRB = 6;
%     meas1 = hRSMeasurements(enb, rxgrid(:,1:L,:));
%     RSRQdB(end+1) = meas1.RSRQdB;
%     RSRPdBm(end+1) = meas1.RSRPdBm;
%     RSSIdBm(end+1) = meas1.RSSIdBm;
    if (isempty(rxgrid))
        fprintf('After timing synchronization, signal is shorter than one subframe so no further demodulation will be performed.\n');
        return;
    end
    % Perform channel estimation
    [hest, nestsf] = lteDLChannelEstimate(enb, cec, rxgrid(:,1:L,:));

    fprintf('Performing MIB decoding...\n');
    pbchIndices = ltePBCHIndices(enb);
    [pbchRx, pbchHest] = lteExtractResources( ...
        pbchIndices, rxgrid(:,1:L,:), hest(:,1:L,:,:));
    
    % Decode PBCH
    [bchBits, pbchSymbols, nfmod4, mib, enb.CellRefP] = ltePBCHDecode(enb, pbchRx, pbchHest, nestsf); 
%     if exist('scat', 'var')
%         release(scat)
%     end   
    release(scat)
    step(scat, pbchSymbols);
%     drawnow;
    % Parse MIB bits
    enb = lteMIB(mib, enb); 
    

    enb.NFrame = enb.NFrame+nfmod4;
    
    fprintf('Cell-wide settings after MIB decoding:\n');
    disp(enb);
    
    if (enb.CellRefP==0)
        warning('MIB decoding failed (enb.CellRefP=0).\n\n');
        MIBok=0;
    end
    if (enb.NDLRB==0)
        warning('MIB decoding failed (enb.NDLRB=0).\n\n');
        MIBok=0;
        
    end
    
%     meas2 = hRSMeasurements(enbcp, rxgrid);
%     if MIBok
%         fprintf(sep)
%         fprintf('\n   MIB decoded')
%         CorrectMIB{frameNum} = enb;
%     end

    if CellIdOk && MIBok
        rxgridmib = lteOFDMDemodulate(enb, downsampled);
        %meas3 = hRSMeasurements(enb, rxgridmib(:,1:L,:));
        fprintf('\n   CellID has been found, and MIB has been successfully decoded. \n')
        fprintf(strcat('\n', sep, '\n'))
        disp(enb)
        MIBFLAG = MIBFLAG + 1;
        MIBDecoded{end+1} = enb;
    end
end
% 
% figure("Name",'RSRQdB')
% plot(x, RSRQdB);
% 
% figure("Name",'RSRPdBm');
% plot(x, RSRPdBm);
% 
% figure("Name",'RSSIdbm');
% plot(x, RSSIdBm);


Correlations = CellIds(:, 2)./CellIds(:, 3);

Result(:, 1) =CellIds(:, 1);

Result(:, 2) = Correlations;

Result(:, 3) = CellIds(:, 3);
Result = sortrows(Result, 3, 'descend');

if length(Result) >= 15
    disp(Result(1:15, :))
else
    disp(Result)
end



c = clock;
if ~exist('center_freq', 'var')
    center_freq = 0;
end
c = string(c);

for i = 2:length(c)-1
    if strlength(c(i)) < 2
        c(i) = strcat('0', c(i));
        
    end
end


if MIBFLAG > 0
    flag = '*';
else
    flag = '';
end

folder = strcat(c(1), c(2), c(3), '/');
strcat(c(1), c(2), c(3), '/');
mkdir(strcat('RecordedCellIdCorrelations/', folder));
filename = strcat('RecordedCellIdCorrelations/', folder, num2str(center_freq/1e6, ...
    '%.1f'),'__', c(1), '_', c(2), '_', c(3), '__', c(4), '_', c(5), flag, '.xlsx');
writematrix(Result, filename);


CellSize = size(MIBDecoded);
if MIBFLAG > 0
    writecell(fieldnames(MIBDecoded{1})', filename, WriteMode="append")
    for i = 1:CellSize(2)
        writetable(struct2table(MIBDecoded{i}), filename, WriteMode="append");
    end
    fprintf('   MIB was decoded %d time(s).\n', MIBFLAG);
end



% 
% 
% fprintf('Restarting reception now that bandwidth (NDLRB=%d) is known...\n',enb.NDLRB);
% 
% % Resample now we know the true bandwidth
% ofdmInfo = lteOFDMInfo(enb);
% if (sr~=ofdmInfo.SamplingRate)
%     if (sr < ofdmInfo.SamplingRate)
%         warning('The received signal sampling rate (%0.3fMs/s) is lower than the desired sampling rate for NDLRB=%d (%0.3fMs/s); PDCCH search / SIB1 decoding may fail.',sr/1e6,enb.NDLRB,ofdmInfo.SamplingRate/1e6);
%     end    
%     fprintf('\nResampling from %0.3fMs/s to %0.3fMs/s...\n',sr/1e6,ofdmInfo.SamplingRate/1e6);
% else
%     fprintf('\nResampling not required; received signal is at desired sampling rate for NDLRB=%d (%0.3fMs/s).\n',enb.NDLRB,sr/1e6);
% end
% nSamples = ceil(ofdmInfo.SamplingRate/round(sr)*size(waveform,1));
% resampled = zeros(nSamples, nRxAnts);
% for i = 1:nRxAnts
%     resampled(:,i) = resample(waveform(:,i), ofdmInfo.SamplingRate, round(sr));
% end
% 
% 
% fprintf('\nPerforming frequency offset estimation...\n');
% delta_f = lteFrequencyOffset(enb, resampled);
% fprintf('Frequency offset: %0.3fHz\n',delta_f);
% resampled = lteFrequencyCorrect(enb, resampled, delta_f);
% 
% % Find beginning of frame
% ThreeBeginningFrames = resampled(1:3*10e-3*ofdmInfo.SamplingRate);
% 
% 
% fprintf('\nPerforming timing offset estimation...\n');
% offset = lteDLFrameOffset(enb, ThreeBeginningFrames); 
% fprintf('Timing offset to frame start: %d samples\n',offset);
% % Aligning signal with the start of the frame
% resampled = resampled(1+offset:end,:);   
% 
% % OFDM demodulation
% fprintf('\nPerforming OFDM demodulation...\n\n');
% 
% 
% NumSamplesPerFrame = ofdmInfo.SamplingRate*10e-3;
% tailSamples = mod(length(resampled),NumSamplesPerFrame);
% resampled = resampled(1:end-tailSamples, :);
% numFullFrames = length(resampled(:,1))/NumSamplesPerFrame;
% 
% LFrame = 10*L; % OFDM symbols per frame
% 
% rxgrid = lteOFDMDemodulate(enb, resampled);
% 
% 
% if (mod(enb.NFrame,2)~=0)
%     if (size(rxgrid,2)>=(L*10))
%         rxgrid(:,1:(L*10),:) = [];
%         fprintf('Skipping frame %d (odd frame number does not contain SIB1).\n\n',enb.NFrame);
%     else
%         rxgrid = [];
%     end
%     enb.NFrame = enb.NFrame + 1;
% end
% 
% % Advance to subframe 5, or terminate if we have less than 5 subframes
% if (size(rxgrid,2)>=(L*5))
%     rxgrid(:,1:(L*5),:) = [];   % Remove subframes 0 to 4
% else
%     rxgrid = [];
% end
% enb.NSubframe = 5;
% 
% if (isempty(rxgrid))
%     fprintf('Received signal does not contain a subframe carrying SIB1.\n\n');
% end
% 
% % Reset the HARQ buffers
% decState = [];
% 
% % While we have more data left, attempt to decode SIB1
% while (size(rxgrid,2) > 0)
% 
%     fprintf('%s\n',sep);
%     fprintf('SIB1 decoding for frame %d\n',mod(enb.NFrame,1024));
%     fprintf('%s\n\n',sep);
% 
%     % Reset the HARQ buffer with each new set of 8 frames as the SIB1
%     % info may be different
%     if (mod(enb.NFrame,8)==0)
%         fprintf('Resetting HARQ buffers.\n\n');
%         decState = [];
%     end
% 
%     % Extract current subframe
%     rxsf = rxgrid(:,1:L,:);
% 
%     % Perform channel estimation
%     [hestsf, nestsf] = lteDLChannelEstimate(enb, cec, rxsf);
% 
%     fprintf('Decoding CFI...\n\n');
%     pcfichIndices = ltePCFICHIndices(enb);  
%     [pcfichRx, pcfichHest] = lteExtractResources(pcfichIndices, rxsf, hestsf);
%  
%     cfiBits = ltePCFICHDecode(enb, pcfichRx, pcfichHest, nestsf);
%     enb.CFI = lteCFIDecode(cfiBits); 
% 
%     
%     fprintf('Decoded CFI value: %d\n\n', enb.CFI);
% 
%     if (strcmpi(enb.DuplexMode,'TDD'))
%         tddConfigs = [1 6 0];
%     else
%         tddConfigs = 0;
%     end
%     dci = {};
%     while (isempty(dci) && ~isempty(tddConfigs))
%   
%         if (strcmpi(enb.DuplexMode,'TDD'))
%             enb.TDDConfig = tddConfigs(1);
%         end
%         tddConfigs(1) = [];
%         pdcchIndices = ltePDCCHIndices(enb); % Get PDCCH indices
%         [pdcchRx, pdcchHest] = lteExtractResources(pdcchIndices, rxsf, hestsf);
%         
%         [dciBits, pdcchSymbols] = ltePDCCHDecode(enb, pdcchRx, pdcchHest, nestsf);
% 
% 
%         fprintf('PDCCH search for SI-RNTI...\n\n');
%         pdcch = struct('RNTI', 65535);
%         dci = ltePDCCHSearch(enb, pdcch, dciBits); % Search PDCCH for DCI
%     end
% 
%     % If DCI was decoded, proceed with decoding PDSCH / DL-SCH
%     if ~isempty(dci)
% 
%         dci = dci{1};
%         fprintf('DCI message with SI-RNTI:\n');
%         disp(dci);
%         % Get the PDSCH configuration from the DCI
%         [pdsch, trblklen] = hPDSCHConfiguration(enb, dci, pdcch.RNTI);
%         pdsch.NTurboDecIts = 5;
%         fprintf('PDSCH settings after DCI decoding:\n');
%         disp(pdsch);
% 
%         % PDSCH demodulation and DL-SCH decoding to recover SIB bits.
%         % The DCI message is now parsed to give the configuration of the
%         % corresponding PDSCH carrying SIB1, the PDSCH is demodulated and
%         % finally the received bits are DL-SCH decoded to yield the SIB1
%         % bits.
% 
%         fprintf('Decoding SIB1...\n\n');
%         % Get PDSCH indices
%         [pdschIndices,pdschIndicesInfo] = ltePDSCHIndices(enb, pdsch, pdsch.PRBSet);
%         [pdschRx, pdschHest] = lteExtractResources(pdschIndices, rxsf, hest);
%         % Decode PDSCH
%         [dlschBits,pdschSymbols] = ltePDSCHDecode(enb, pdsch, pdschRx, pdschHest, nestsf);
%         % Decode DL-SCH with soft buffer input/output for HARQ combining
%         if ~isempty(decState)
%             fprintf('Recombining with previous transmission.\n\n');
%         end
%         [sib1, crc, decState] = lteDLSCHDecode(enb, pdsch, trblklen, dlschBits, decState);
% 
%         % Compute PDSCH EVM
%         recoded = lteDLSCH(enb, pdsch, pdschIndicesInfo.G, sib1);
%         remod = ltePDSCH(enb, pdsch, recoded);
%         [~,refSymbols] = ltePDSCHDecode(enb, pdsch, remod);
%         release(pdschEVM);
%         [rmsevm,peakevm] = step(pdschEVM, refSymbols{1}, pdschSymbols{1});
%         fprintf('PDSCH RMS EVM: %0.3f%%\n',rmsevm);
%         fprintf('PDSCH Peak EVM: %0.3f%%\n\n',peakevm);
% 
%         fprintf('SIB1 CRC: %d\n',crc);
%         if crc == 0
%             fprintf('Successful SIB1 recovery.\n\n');
%         else
%             fprintf('SIB1 decoding failed.\n\n');
%         end
% 
%     else
%         % indicate that DCI decoding failed
%         fprintf('DCI decoding failed.\n\n');
%     end
% 
%     % Update channel estimate plot
%     figure(channelFigure);
%     surf(abs(hest(:,:,1,1)));
%     hSIB1RecoveryExamplePlots(channelFigure);
%     channelFigure.CurrentAxes.XLim = [0 size(hest,2)+1];
%     channelFigure.CurrentAxes.YLim = [0 size(hest,1)+1];
% 
%     % Skip 2 frames and try SIB1 decoding again, or terminate if we
%     % have less than 2 frames left.
%     if (size(rxgrid,2)>=(L*20))
%         rxgrid(:,1:(L*20),:) = [];   % Remove 2 more frames
%     else
%         rxgrid = []; % Less than 2 frames left
%     end
%     enb.NFrame = enb.NFrame+2;
% 
% end