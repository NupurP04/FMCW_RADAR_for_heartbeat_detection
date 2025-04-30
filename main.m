% Multi-Person FMCW Radar for Heartbeat Detection
close all;
clc;
clear all;

% Importing Data
% load two_people_1.mat channel_data  % (taking one chirp per frame) frame sampling points

% Step 1: Reading the data from the .mat file
loadedData = load('two_people_1.mat');

% Step 2: Viewing the variable names in the loaded data
variableNames = fieldnames(loadedData);

% Step 3: Saveing the data as a new custom variable name
% Assuming there is only one variable in the .mat file
channel_data = variableNames{1};

% New custom variable name
channel_data = loadedData.(channel_data);

% Checking the data of the custom variable name
% disp(myCustomVarName);

% Radar Parameters
f0=600e8; % Start frequency
numADCSamples = 342; % Number of ADC samples
freqSlopeRate = 18e12; % Frequency slope rate, unit: MHz/us
adcSampleRate = 6e6; % Fast time ADC sampling rate, unit: Ksps
Ts = 50e-3; % Time between each frame
% Tc = numADCSamples/adcSampleRate;
Tc = 57e-6;
B = 1e9;
% B = Tc*freqSlopeRate;
useFramesNum=900;
rangeFFTNum=342;
deltaR = 3e8/(2*B);
d=0.025;
lambda=2*d;

% Data Processing
loop_cnt = floor(length(channel_data(:,1,1)) / useFramesNum);

for k = 1:loop_cnt

    use_channel_data = channel_data((k - 1) * useFramesNum + 1 : k * useFramesNum,:,:);

    % Mean Cancellation and Pulse Compression
    rangeProfile = MTI(use_channel_data,0,rangeFFTNum);
    sum_rangProfile = sum(abs(rangeProfile(:,:,1)),1);
    [~,targetIndex] = max(sum_rangProfile);
    
    % DOA Estimation
    searchAngleRange = 60; % Angle spectrum search range is +-searchAngleRange degrees

    [~,azimuSpectrogram,Rxv] = IWR6843ISK_DOA(rangeProfile,2,useFramesNum,searchAngleRange);
    azimuSpectrogram=flipud(azimuSpectrogram);
    % Finding Targets at Different Angles
    maxAzimu = max(azimuSpectrogram,[],2);
    %max(A,[],2) returns a column vector containing the maximum value of each row

    [values,peaks_index] = findpeaks(maxAzimu,'minpeakheight',mean(azimuSpectrogram,'all'));
    maxAzimu=sum(azimuSpectrogram,2);
   figure(1);  
   plot(1:length(maxAzimu),10*log10(maxAzimu));    
   hold on; grid on
   plot(peaks_index,10*log10(maxAzimu(peaks_index)),'bd');
   xlabel('Angle/°');ylabel('Gain (dB)');title('MVDR Angle Estimation');

   hold off

   figure(2);
   data_max=findLocalMaximaInIdx(azimuSpectrogram,45,2);
   imagesc(-searchAngleRange:searchAngleRange,(1 :rangeFFTNum)*deltaR,abs(azimuSpectrogram.'));
   ylabel('Distance (m)');xlabel('Angle (°)');title('Range-Angle Spectrum');
   axis xy

   % Range-Angle Fan Plot Drawing
    figure(3);
    R = (1 :rangeFFTNum)*deltaR;
    ang_ax =-searchAngleRange:searchAngleRange;
    X = R'*cosd(ang_ax); Y = R'*sind(ang_ax); %
    pcolor(Y,X,abs(azimuSpectrogram.'));
    axis equal tight  % x-axis unit scale and y-axis unit scale length are equal, best reflecting the actual curve;
    shading interp % Shading, making the colors transition smoothly
    axis off
    initialAz = -90; endAz = 90; % Label text
    text((max(R)+10)*cosd(initialAz),(max(R))*sind(initialAz),...
    [num2str(initialAz) '^o']);
    text((max(R)+10)*cosd(endAz),(max(R))*sind(endAz),...
    [num2str(endAz) '^o']);
    text((max(R)+10)*cosd(0),(max(R))*sind(0),[num2str(0) '^o']);

   % Selecting the Top Two Values with the Largest Magnitude in Angle
   [values_sort,index] =sort(values,'descend');
   peaks_index_sort =peaks_index(index);
   peaks_index_max =peaks_index_sort([1,5]);

   % Calculating Respiration and Heartbeat of Different Targets
   target_rangeProfile  = zeros(useFramesNum,length(peaks_index_max)); 
   % Storing the range data of the target after angle filtering, 
   % the subsequent respiration and heartbeat are processed based on this data, one column is the range data of one target  
    target_rangeProfile  = zeros(useFramesNum,size(data_max,1));
    for i=1:size(data_max,1)
        xt = squeeze(rangeProfile(:,data_max(i,2)-1,:));
        detAngle = -searchAngleRange + data_max(i,1) * (searchAngleRange * 2 / length(azimuSpectrogram(:,1)));

        fai = 2 * pi * sin(detAngle / 180 * pi) * d / lambda;

        aTheta = [1,exp(-1j*1*fai),exp(-1j*2*fai),exp(-1j*3*fai),exp(-1j*4*fai),exp(-1j*5*fai),exp(-1j*6*fai),exp(-1j*7*fai)].';%

        Wopt = (Rxv  * aTheta) / (aTheta' * Rxv  * aTheta);

        target_rangeProfile(:,i) = xt * Wopt; 
    end

    % Estimate the Respiration and Heartbeat of Multiple Targets Sequentially,
    % and displaying and plotting the results
    [breathRate,heartRate] = get_heartBreath_rate(target_rangeProfile,1/Ts);

end
 
