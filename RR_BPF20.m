function Hd = RR_BPF20
Fs = 20;  % Sampling Frequency
N   = 4;    % Order
Fc1 = 0.1;  % First Cutoff Frequency
Fc2 = 0.6;  % Second Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
Hd = design(h, 'butter');
