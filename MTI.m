function Cancel_PulseCompre_Data = MTI(framesData,cancelFlag,rangeFFTNum,pluseCancelSpace)
%  Perform clutter cancellation and pulse compression on a frame of radar echo data.
meanCancel_mode = 0;
numChirp = size(framesData,1);
% Mean Cancellation
if(cancelFlag == meanCancel_mode)
    meanCancelData = framesData - repmat(mean(framesData,1),[numChirp 1]); % Mean cancellation to remove DC component
    % Pulse compression X = ifft(Y,n,dim) returns the n-point inverse Fourier transform along dimension dim. 
    % For example, if Y is a matrix, then ifft(Y,n,2) returns the n-point inverse transform of each row.
    Cancel_PulseCompre_Data = fft(meanCancelData,rangeFFTNum,2); % Pulse compression X = ifft(Y,n,dim) returns the n-point inverse Fourier transform along dimension dim. For example, if Y is a matrix, then ifft(Y,n,2) returns the n-point inverse transform of each row.
end
