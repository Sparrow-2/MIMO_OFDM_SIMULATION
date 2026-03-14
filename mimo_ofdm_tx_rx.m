close all; clear; clc;
% System Parameters
Config.FFT_Size     = 512;
Config.Ovs_Factor   = 4;
Config.IFFT_Size    = Config.FFT_Size * Config.Ovs_Factor;
Config.CP_Ratio     = 1/4;
Config.Zeros_Left   = 43;
Config.Zeros_Right  = 43;
Config.Pilot_Step   = 8;
Config.Pilot_Value  = 3;
Config.ScramblerInit = [1 0 1 1 0 1 1];

T_u_us = 2000;
Config.T_sample = (T_u_us * 1e-6) / Config.IFFT_Size;
Ts = Config.T_sample;

fprintf('--- TX MIMO OFDM 2x2 (Oversampling x%d) ---\n', Config.Ovs_Factor);

% TX

textData = fileread("text.txt");
msgLen = length(textData);
fprintf('Text length %d chars \n', msgLen);

binStr = dec2bin(msgLen, 10); 
infoBits = (binStr - '0').';
Config.InfoBitsLen = length(infoBits); 

infoBitsAnt1 = infoBits(1 : end/2);
infoBitsAnt2 = infoBits(end/2 + 1 : end);

fprintf('System Information: %s (MIMO BPSK: 8 bits Ant1 + 8 bits Ant2)\n', binStr);

bits = ConvertTextToBits(textData);
bitsScrambled = ScrambleBits(bits, Config.ScramblerInit);

symbols = ModulateQAM16_Integer(bitsScrambled);
%MIMO Encoder
symbolsAnt1 = symbols(1:2:end);
symbolsAnt2 = symbols(2:2:end);

maxLength = max(length(symbolsAnt1), length(symbolsAnt2));
symbolsAnt1 = PadVector(symbolsAnt1, maxLength);
symbolsAnt2 = PadVector(symbolsAnt2, maxLength);

activeCarriers = Config.FFT_Size - Config.Zeros_Left - Config.Zeros_Right - 1;
pilotIndices = 1 : Config.Pilot_Step : activeCarriers;
pilotPowers = (-1).^(0:length(pilotIndices)-1)';
basePilots = Config.Pilot_Value * pilotPowers;

[txSignal1, numSyms, papr1_vec, gridMtx1] = GenerateOFDMFrame(symbolsAnt1, basePilots, 1, Config, infoBits);
[txSignal2, ~, papr2_vec, gridMtx2] = GenerateOFDMFrame(symbolsAnt2, basePilots, 2, Config, infoBits);

fprintf('\n PAPR [dB] (after oversampling)\n');
fprintf('Frame \t| Antenna 1\t| Antenna 2\n');
for k = 1:numSyms
    fprintf('Frame  %d\t| %6.2f\t| %6.2f\n', k, papr1_vec(k), papr2_vec(k));
end

% Radio Channel
fprintf('Radio Channel');

SNR_dB = 15;
CFO_Hz = 15;

% Adding zeros preamble 
noisePreambleLen = 1000;
rx1_Raw = [zeros(noisePreambleLen, 1); txSignal1];
rx2_Raw = [zeros(noisePreambleLen, 1); txSignal2];

% Adding AWGN 
sigPow1 = mean(abs(rx1_Raw).^2);
sigPow2 = mean(abs(rx2_Raw).^2);
noisePow1 = sigPow1 / (10^(SNR_dB/10));
noisePow2 = sigPow2 / (10^(SNR_dB/10));
noise1 = sqrt(noisePow1/2) * (randn(size(rx1_Raw)) + 1j*randn(size(rx1_Raw)));
noise2 = sqrt(noisePow2/2) * (randn(size(rx2_Raw)) + 1j*randn(size(rx2_Raw)));
rx1_Noisy = rx1_Raw + noise1;
rx2_Noisy = rx2_Raw + noise2;

% Adding Carrier Frequency Offset
tVecShift = (0:length(rx1_Noisy)-1).' * Ts;
freqShift = exp(1j * 2 * pi * CFO_Hz * tVecShift);
rx1_Final = rx1_Noisy .* freqShift;
rx2_Final = rx2_Noisy .* freqShift;

% RX
fprintf('RX ');

[startIndex, estCFO1] = SynchronizeAndEstimateFrequency(rx1_Final, Config);
[startIndex2, estCFO2] = SynchronizeAndEstimateFrequency(rx2_Final, Config);
estCFO = (estCFO1 + estCFO2) / 2;

fprintf('Frame start detected at sample %d\n', startIndex);
fprintf('CFO Est: %.2f Hz (Real: %.2f Hz)\n', estCFO, CFO_Hz);

tVector = (0:length(rx1_Final)-1).' * Ts;
correctionFactor = exp(-1j * 2 * pi * estCFO * tVector);
rx1_Corrected = rx1_Final .* correctionFactor;
rx2_Corrected = rx2_Final .* correctionFactor;

cpLen = round(Config.IFFT_Size * Config.CP_Ratio);
symLen = Config.IFFT_Size + cpLen;
totalFrameLen = numSyms * symLen;

if startIndex + totalFrameLen - 1 > length(rx1_Corrected)
    startIndex = length(rx1_Corrected) - totalFrameLen + 1;
end
rx1_Cut = rx1_Corrected(startIndex : startIndex + totalFrameLen - 1);
rx2_Cut = rx2_Corrected(startIndex : startIndex + totalFrameLen - 1);


[rxSymsAnt1, rxSymsAnt2, estInfoBits, rxSoft1, rxSoft2] = DOFDM_VBLAST(rx1_Cut, rx2_Cut, basePilots, Config);
rxMsgLen = bin2dec(num2str(estInfoBits'));
fprintf('System information');
fprintf('Received Length %d (Transmitted Lenght: %d)\n', rxMsgLen, msgLen);

rxSymsTotal = zeros(length(rxSymsAnt1) + length(rxSymsAnt2), 1);
rxSymsTotal(1:2:end) = rxSymsAnt1;
rxSymsTotal(2:2:end) = rxSymsAnt2;

rxSymsTotal = rxSymsTotal(1:length(symbols));

rxBitsRaw = DemodulateQAM16_Integer(rxSymsTotal);
rxBitsValid = rxBitsRaw(1:length(bitsScrambled));
rxBitsDescrambled = ScrambleBits(rxBitsValid, Config.ScramblerInit);
recoveredText = ConvertBitsToText(rxBitsDescrambled);

fprintf('Received text \n');
disp(recoveredText);

errorBits = sum(bits ~= rxBitsDescrambled);
if errorBits == 0
    fprintf('\nStatus: SUCCESS 0 ERRORS.\n');
else
    fprintf('\nStatus: ERRORS! %d bit errors.\n', errorBits);
end

bitsToShow = 50;
Fs_MHz = (1/Ts) / 1e6;

spec_xlim = [-Fs_MHz/2, Fs_MHz/2];

N_FFT_Plot = 8192;
f_axis_indices = (-N_FFT_Plot/2 : N_FFT_Plot/2-1);
f_axis = f_axis_indices * (Fs_MHz / N_FFT_Plot);
% Plots
% Input Data
figure('Name', '1. Input Bits', 'Color', 'w');
stairs(0:bitsToShow-1, bits(1:bitsToShow), 'b', 'LineWidth', 2);
grid on;
axis([-1 bitsToShow -0.2 1.2]);
title('Part of input bits');
xlabel('Index');
ylabel('Value');

% 2. Constellation TX
figure('Name', '2. Constellation Tx', 'Color', 'w');
plot(symbolsAnt1, 'b.', 'MarkerSize', 12);
hold on;
plot(symbolsAnt2, 'r.', 'MarkerSize', 12);
grid on;
axis equal;
xlim([-4 4]);
ylim([-4 4]);
title('Symbols QAM16 (Before OFDM)');
legend('Antenna 1', 'Antenna 2');

% 3. TX 1 - Time
figure('Name', '3. Tx1 - Time', 'Color', 'w');
tAxisTx = (0:length(txSignal1)-1) * Ts * 1e6;
plot(tAxisTx, abs(txSignal1), 'b');
grid on;
xlim([0 max(tAxisTx)]);
titleStr = sprintf('Signal Tx 1 (Oversampled x%d)', Config.Ovs_Factor);
title(titleStr);
xlabel('Time [\mus]');
ylabel('|Magnitude|');

% 4. TX 1 - Spectrum
figure('Name', '4. Tx1 - Spectrum', 'Color', 'w');
spec_tx1 = calc_spectrum_dB(txSignal1, N_FFT_Plot);
plot(f_axis, spec_tx1, 'b', 'LineWidth', 1.2);
grid on;
xlim(spec_xlim);
ylim([-100 0]);
title('Spectrum Tx 1 ');
xlabel('Frequency [MHz]');
ylabel('Magnitude [dB]');

% 5. TX 2 - Time
figure('Name', '5. Tx2 - Time', 'Color', 'w');
plot(tAxisTx, abs(txSignal2), 'r');
grid on;
xlim([0 max(tAxisTx)]);
titleStr = sprintf('Signal Tx 2 (Oversampled x%d)', Config.Ovs_Factor);
title(titleStr);
xlabel('Time [\mus]');
ylabel('|Magnitude|');

% 6. TX 2 - Spectrum
figure('Name', '6. Tx2 - Spectrum', 'Color', 'w');
spec_tx2 = calc_spectrum_dB(txSignal2, N_FFT_Plot);
plot(f_axis, spec_tx2, 'r', 'LineWidth', 1.2);
grid on;
xlim(spec_xlim);
ylim([-100 0]);
title('Spectrum Tx 2 ');
xlabel('Frequency [MHz]');
ylabel('Magnitude [dB]');

% 7. RX 1 - Time
figure('Name', '7. Rx1 - Time', 'Color', 'w');
tAxisRx = (0:length(rx1_Final)-1) * Ts * 1e6;
plot(tAxisRx, abs(rx1_Final), 'b');
grid on;
xlim([0 max(tAxisRx)]);
title('Signal Rx 1 (Noisy + Preamble)');
xlabel('Time [\mus]');
ylabel('|Magnitude|');

% 8. RX 1 - Spectrum
figure('Name', '8. Rx1 - Spectrum', 'Color', 'w');
spec_rx1 = calc_spectrum_dB(rx1_Final, N_FFT_Plot);
plot(f_axis, spec_rx1, 'b', 'LineWidth', 1);
grid on;
xlim(spec_xlim);
ylim([-90 10]);
title('Spectrum Rx 1 (with AWGN)');
xlabel('Frequency [MHz]');
ylabel('Magnitude [dB]');

% 9. RX 2 - Time
figure('Name', '9. Rx2 - Time', 'Color', 'w');
plot(tAxisRx, abs(rx2_Final), 'r');
grid on;
xlim([0 max(tAxisRx)]);
title('Signal Rx 2 (Noisy + Preamble)');
xlabel('Time [\mus]');
ylabel('|Magnitude|');

% 10. RX 2 - Spectrum
figure('Name', '10. Rx2 - Spectrum', 'Color', 'w');
spec_rx2 = calc_spectrum_dB(rx2_Final, N_FFT_Plot);
plot(f_axis, spec_rx2, 'r', 'LineWidth', 1);
grid on;
xlim(spec_xlim);
ylim([-90 10]);
title('Spectrum Rx 2 (with AWGN)');
xlabel('Frequency [MHz]');
ylabel('Magnitude [dB]');

% 11. Constellation RX 
figure('Name', '11. Constellation Rx', 'Color', 'w');
nSyms1 = ceil(length(symbols)/2);
nSyms2 = floor(length(symbols)/2);
rxSoftTotal = [rxSoft1(1:nSyms1); rxSoft2(1:nSyms2)];
plot(rxSoftTotal, 'g.', 'MarkerSize', 6);
hold on;

[GX, GY] = meshgrid([-3 -1 1 3], [-3 -1 1 3]);
plot(GX(:), GY(:), 'r+', 'MarkerSize', 12, 'LineWidth', 2);
grid on;
axis equal;
xlim([-5 5]);
ylim([-5 5]);
titleStr = sprintf('Received Symbols  SNR=%ddB', SNR_dB);
title(titleStr);
xlabel('I');
ylabel('Q');
legend('Received', 'Ideal');


%  Helper Functions

function spec_dB = calc_spectrum_dB(sig, N_FFT)
    fft_result = fft(sig, N_FFT);
    fft_shifted = ManualFFTShift(fft_result);
    magnitude = abs(fft_shifted);
    normalized = magnitude / length(sig);
    spec_dB = 20 * log10(normalized);
end

function v = PadVector(x, targetLen)
    v = x;
    if length(v) < targetLen
        v(end+1:targetLen) = 0;
    end
end

function [syms1, syms2, infoBitsDecoded, rawSyms1, rawSyms2] = DOFDM_VBLAST(rx1, rx2, pilots, cfg)
    cpLen = round(cfg.IFFT_Size * cfg.CP_Ratio);
    symbolLen = cfg.IFFT_Size + cpLen;
    numSyms = floor(length(rx1) / symbolLen);
    
    activeLen = cfg.FFT_Size - cfg.Zeros_Left - cfg.Zeros_Right - 1;
    dcIdx = cfg.FFT_Size / 2 + 1;
    
    allPilotIdx = 1 : cfg.Pilot_Step : activeLen;
    allDataIdx = find(~ismember(1:activeLen, allPilotIdx));
    
    nInfoBitsTotal = cfg.InfoBitsLen;
    nInfoBitsPerAnt = nInfoBitsTotal / 2;
    
    sysInfoIdx = allDataIdx(1:nInfoBitsPerAnt);
    payloadDataIdx = allDataIdx(nInfoBitsPerAnt+1:end);
    
    idxP1 = allPilotIdx(1:2:end);
    valP1 = pilots(1:2:end);
    idxP2 = allPilotIdx(2:2:end);
    valP2 = pilots(2:2:end);
    
    syms1 = []; syms2 = [];
    rawSyms1 = []; rawSyms2 = [];
    
    [~, qamMap] = ModulateQAM16_Integer([0 0 0 0]);
    infoSoftAcc1 = zeros(nInfoBitsPerAnt, 1); 
    infoSoftAcc2 = zeros(nInfoBitsPerAnt, 1); 

    for k = 1:numSyms
        idxStart = (k-1) * symbolLen + cpLen + 1;
        idxEnd = k * symbolLen;
        rawY1 = fft(rx1(idxStart : idxEnd));
        rawY2 = fft(rx2(idxStart : idxEnd));
        Y1_huge = ManualFFTShift(rawY1);
        Y2_huge = ManualFFTShift(rawY2);
        startIdx = (cfg.IFFT_Size - cfg.FFT_Size) / 2 + 1;
        Y1_full = Y1_huge(startIdx : startIdx + cfg.FFT_Size - 1);
        Y2_full = Y2_huge(startIdx : startIdx + cfg.FFT_Size - 1);
        Y1_part1 = Y1_full(cfg.Zeros_Left+1 : dcIdx-1);
        Y1_part2 = Y1_full(dcIdx+1 : cfg.FFT_Size-cfg.Zeros_Right);
        Y1 = [Y1_part1; Y1_part2];
        Y2_part1 = Y2_full(cfg.Zeros_Left+1 : dcIdx-1);
        Y2_part2 = Y2_full(dcIdx+1 : cfg.FFT_Size-cfg.Zeros_Right);
        Y2 = [Y2_part1; Y2_part2];
        
        H11_pilots = Y1(idxP1) ./ valP1;
        H21_pilots = Y2(idxP1) ./ valP1;
        H12_pilots = Y1(idxP2) ./ valP2;
        H22_pilots = Y2(idxP2) ./ valP2;
        queryPoints = 1:activeLen;
        H11 = ManualInterp1(idxP1, H11_pilots, queryPoints).';
        H21 = ManualInterp1(idxP1, H21_pilots, queryPoints).';
        H12 = ManualInterp1(idxP2, H12_pilots, queryPoints).';
        H22 = ManualInterp1(idxP2, H22_pilots, queryPoints).';

        for i = 1:length(sysInfoIdx)
             sc = sysInfoIdx(i);
             yVec = [Y1(sc); Y2(sc)];
             H_Mtx = [H11(sc), H12(sc); H21(sc), H22(sc)];
             ants = [1, 2];
             detected = [0; 0];
             G = pinv(H_Mtx);
             norms = sum(abs(G).^2, 2);
             [~, bestIdx] = min(norms);
             bestAnt = ants(bestIdx);
             z = G(bestIdx, :) * yVec;
             if real(z) > 0, est = 1; else, est = -1; end
             detected(bestAnt) = est;
             yRem = yVec - H_Mtx(:, bestIdx) * est;
             H_Mtx(:, bestIdx) = [];
             ants(bestIdx) = [];
             G2 = pinv(H_Mtx);
             z2 = G2 * yRem;
             if real(z2) > 0, est2 = 1; else, est2 = -1; end
             detected(ants(1)) = est2;
             infoSoftAcc1(i) = infoSoftAcc1(i) + detected(1);
             infoSoftAcc2(i) = infoSoftAcc2(i) + detected(2);
        end
        
        dec1 = zeros(length(payloadDataIdx), 1);
        dec2 = zeros(length(payloadDataIdx), 1);
        raw1 = zeros(length(payloadDataIdx), 1); 
        raw2 = zeros(length(payloadDataIdx), 1); 
        
        for i = 1:length(payloadDataIdx)
            sc = payloadDataIdx(i);
            yVec = [Y1(sc); Y2(sc)];
            H_Mtx = [H11(sc), H12(sc); H21(sc), H22(sc)];
            
            ants = [1, 2];
            detectedSyms = [0; 0];
            softVals = [0; 0]; 
            
            G = pinv(H_Mtx);
            [~, bestIdx] = min(sum(abs(G).^2, 2));
            bestAnt = ants(bestIdx);
            
            z = G(bestIdx, :) * yVec; 
            softVals(bestAnt) = z;    
            
            [~, minIdx] = min(abs(z - qamMap));
            symbolEst = qamMap(minIdx);
            detectedSyms(bestAnt) = symbolEst;
            
            yCancelled = yVec - H_Mtx(:, bestIdx) * symbolEst;
            H_Mtx(:, bestIdx) = [];
            ants(bestIdx) = [];
            
            G_rem = pinv(H_Mtx);
            z2 = G_rem * yCancelled; 
            softVals(ants(1)) = z2;  
            
            [~, minIdx2] = min(abs(z2 - qamMap));
            detectedSyms(ants(1)) = qamMap(minIdx2);
            
            dec1(i) = detectedSyms(1);
            dec2(i) = detectedSyms(2);
            raw1(i) = softVals(1);
            raw2(i) = softVals(2);
        end
        syms1 = [syms1; dec1];
        syms2 = [syms2; dec2];
        rawSyms1 = [rawSyms1; raw1];
        rawSyms2 = [rawSyms2; raw2];
    end
    
    bitsAnt1 = infoSoftAcc1 > 0;
    bitsAnt2 = infoSoftAcc2 > 0;
    infoBitsDecoded = [bitsAnt1; bitsAnt2];
end

function [txSig, numSyms, symPAPR, gridMatrix] = GenerateOFDMFrame(qamData, pilots, antID, cfg, infoBits)
    activeCarriers = cfg.FFT_Size - cfg.Zeros_Left - cfg.Zeros_Right - 1;
    pilotIdx = 1 : cfg.Pilot_Step : activeCarriers;
    
    isPilot = zeros(1, activeCarriers);
    isPilot(pilotIdx) = 1;
    allDataIdx = find(~isPilot);
   
    nInfoBitsTotal = length(infoBits); 
    nInfoBitsPerAnt = nInfoBitsTotal / 2; 
    
    sysInfoIdx = allDataIdx(1:nInfoBitsPerAnt);
    payloadDataIdx = allDataIdx(nInfoBitsPerAnt+1:end);
    
    numDataCarriers = length(payloadDataIdx);
    numSyms = ceil(length(qamData) / numDataCarriers);
    
    if antID == 1
        myInfoBits = infoBits(1 : nInfoBitsPerAnt);
    else
        myInfoBits = infoBits(nInfoBitsPerAnt+1 : end);
    end
    bpskInfo = 2 * double(myInfoBits) - 1; 

    paddingNeeded = numSyms * numDataCarriers - length(qamData);
    qamPadded = [qamData; zeros(paddingNeeded, 1)];
    dataBlocks = reshape(qamPadded, numDataCarriers, numSyms);
    
    dcIdx = cfg.FFT_Size / 2 + 1;
    symPAPR = zeros(numSyms, 1);
    cpLen = round(cfg.IFFT_Size * cfg.CP_Ratio);
    txSig = [];
    gridMatrix = zeros(cfg.FFT_Size, numSyms);
    
    for k = 1:numSyms
        symbolVec = zeros(activeCarriers, 1);
        symbolVec(payloadDataIdx) = dataBlocks(:, k);
        symbolVec(sysInfoIdx) = bpskInfo;
        
        currentPilots = zeros(length(pilotIdx), 1);
        if antID == 1
            activePilotIdx = 1:2:length(pilotIdx);
        else
            activePilotIdx = 2:2:length(pilotIdx);
        end
        currentPilots(activePilotIdx) = pilots(activePilotIdx);
        symbolVec(pilotIdx) = currentPilots;
        
        smallGrid = zeros(cfg.FFT_Size, 1);
        lenPart1 = (dcIdx - 1) - cfg.Zeros_Left;
        smallGrid(cfg.Zeros_Left+1 : dcIdx-1) = symbolVec(1:lenPart1);
        smallGrid(dcIdx+1 : cfg.FFT_Size-cfg.Zeros_Right) = symbolVec(lenPart1+1:end);
        gridMatrix(:, k) = smallGrid;
        
        bigGrid = zeros(cfg.IFFT_Size, 1);
        startIdx = (cfg.IFFT_Size - cfg.FFT_Size) / 2 + 1;
        bigGrid(startIdx : startIdx + cfg.FFT_Size - 1) = smallGrid;
        
        timeSymbol = ifft(ManualIFFTShift(bigGrid));
        peak = max(abs(timeSymbol).^2);
        avg = mean(abs(timeSymbol).^2);
        symPAPR(k) = 10 * log10(peak / avg);
        cp = timeSymbol(end-cpLen+1:end);
        txSig = [txSig; cp; timeSymbol];
    end
end

function [startIdx, estCFO] = SynchronizeAndEstimateFrequency(rxSignal, cfg)
    fftLen = cfg.IFFT_Size;
    cpLen = round(fftLen * cfg.CP_Ratio);
    
    searchWin = length(rxSignal) - fftLen - cpLen - 20;
    if searchWin < 10
        searchWin = 10;
    end
    
    metric = zeros(searchWin, 1);
    for d = 1:searchWin
        partCP = rxSignal(d : d + cpLen - 1);
        partEnd = rxSignal(d + fftLen : d + fftLen + cpLen - 1);
        metric(d) = sum(partCP .* conj(partEnd));
    end
    [~, startIdx] = max(abs(metric));
    phaseDiff = angle(metric(startIdx));
    estCFO = -phaseDiff / (2 * pi * cfg.T_sample * fftLen);
end

function [s, mapping] = ModulateQAM16_Integer(bits)
    I_levels = [-3, -1, 3, 1];
    Q_levels = [3, 1, -3, -1];
    mapping = zeros(16, 1);
    for i = 1:4
        for j = 1:4
            idx = (i-1) * 4 + j;
            mapping(idx) = I_levels(i) + 1j * Q_levels(j);
        end
    end
    remainder = mod(length(bits), 4);
    if remainder ~= 0
        bits = [bits, zeros(1, 4 - remainder)];
    end
    numSyms = length(bits) / 4;
    s = zeros(numSyms, 1);
    for k = 1:numSyms
        chunk = bits((k-1)*4 + 1 : k*4);
        idx = chunk(1)*8 + chunk(2)*4 + chunk(3)*2 + chunk(4)*1 + 1;
        s(k) = mapping(idx);
    end
end

function bits = DemodulateQAM16_Integer(symbols)
    bits = zeros(1, 4 * length(symbols));
    for k = 1:length(symbols)
        I = real(symbols(k));
        Q = imag(symbols(k));
        
        if I < -2,      bI = [0 0];
        elseif I < 0,   bI = [0 1];
        elseif I < 2,   bI = [1 1];
        else,           bI = [1 0];
        end
        
        if Q > 2,       bQ = [0 0];
        elseif Q > 0,   bQ = [0 1];
        elseif Q > -2,  bQ = [1 1];
        else,           bQ = [1 0];
        end
        
        bits((k-1)*4+1 : k*4) = [bI, bQ];
    end
end

function bits = ConvertTextToBits(strText)
    bytes = unicode2native(strText, 'UTF-8');
    binMatrix = dec2bin(bytes, 8).';
    bits = reshape(binMatrix - '0', 1, []);
end

function strText = ConvertBitsToText(bits)
    len = floor(length(bits) / 8) * 8;
    bits = bits(1:len);
    binMatrix = reshape(char(bits + '0'), 8, []).';
    bytes = uint8(bin2dec(binMatrix).');
    strText = native2unicode(bytes, 'UTF-8');
end

function outBits = ScrambleBits(inBits, initialState)
    reg = initialState;
    outBits = zeros(size(inBits));
    for k = 1:length(inBits)
        feedback = xor(reg(1), reg(5));
        outBits(k) = xor(inBits(k), feedback);
        reg = [feedback, reg(1:end-1)];
    end
end

function y = ManualFFTShift(x)
    mid = ceil(length(x) / 2);
    y = [x(mid+1:end); x(1:mid)];
end

function y = ManualIFFTShift(x)
    mid = floor(length(x) / 2);
    y = [x(mid+1:end); x(1:mid)];
end

function yi = ManualInterp1(known_x, known_y, query_x)
    yi = zeros(size(query_x));
    for k = 1:length(query_x)
        val = query_x(k);
        if val <= known_x(1)
            yi(k) = known_y(1);
            continue;
        end
        if val >= known_x(end)
            yi(k) = known_y(end);
            continue;
        end
        idx_left = find(known_x <= val, 1, 'last');
        x0 = known_x(idx_left);
        x1 = known_x(idx_left+1);
        y0 = known_y(idx_left);
        y1 = known_y(idx_left+1);
        slope = (y1 - y0) / (x1 - x0);
        yi(k) = y0 + slope * (val - x0);
    end
end