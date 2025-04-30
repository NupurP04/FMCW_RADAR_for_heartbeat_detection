function [pitchSpectrogram,azimuSpectrogram,Rxv] = IWR6843ISK_DOA(rangeProfile,doaFlag,useChirpNum,searchAngleRange,wopt_flag,wopt)
% Perform DOA estimation for azimuth and elevation directions on a frame of data from the IWR1642ODS radar,
% generating an angle spectrum.

FFT_DOA_MODE   = 0;
CBF_DOA_MODE   = 1;
MVDR_DOA_MODE  = 2;

% Calculate angle spectrum using CBF and MVDR
% Only one-dimensional CBF method is used to calculate azimuth and elevation angle spectra separately
rangNum = size(rangeProfile,2);
lambda = 5e-3; % Wavelength
d = lambda / 2; % Element spacing
if (doaFlag == CBF_DOA_MODE) || (doaFlag == MVDR_DOA_MODE)
    SEARCH_ANGLE_RANGE = searchAngleRange/180*pi; % Search angle range (in radians)
    SEARCH_ANGLE_SPACE = 1/180*pi; % Search angle interval (in radians)
    fai = 2 * pi * sin(-SEARCH_ANGLE_RANGE:SEARCH_ANGLE_SPACE:SEARCH_ANGLE_RANGE) * d / lambda; % fai = 2πd sin(θ)/λ 
    azimuSpectrogram = zeros(length(fai),rangNum);
    pitchSpectrogram = zeros(length(fai),rangNum);
    for r = 1:rangNum
        % Azimuth
         xt = squeeze(rangeProfile(1:useChirpNum,r,:))';
         Rx = xt * xt'; % Calculate covariance matrix 
        for an = 1:length(fai)
            aTheta = [1,exp(-j*1*fai(an)),exp(-j*2*fai(an)),exp(-j*3*fai(an)),exp(-j*4*fai(an)),exp(-j*5*fai(an)),exp(-j*6*fai(an)),exp(-j*7*fai(an))];%
            if (doaFlag == CBF_DOA_MODE) % CBF
                azimuSpectrogram(an,r) = abs(aTheta * Rx * aTheta'); 
            else % MVDR
                Rxv = pinv(Rx);  
                azimuSpectrogram(an,r) = 1 / abs(aTheta * Rxv * aTheta');
            end
        end
    end
end
AzimuSpectrogram = azimuSpectrogram / max(azimuSpectrogram);
end
