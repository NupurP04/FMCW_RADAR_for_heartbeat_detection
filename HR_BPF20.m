function Hd = HR_BPF20
%HR_BPF20 Returns a discrete-time filter object.
Fs = 20;  % Sampling Frequency
N   = 8;    % Order
Fc1 = 0.9;  % First Cutoff Frequency
Fc2 = 2;    % Second Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
Hd = design(h, 'butter');

