function Hd = design_filter(N,Fpass,Fstop,Fs)
%DESIGN_FIL Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 9.2 and the Signal Processing Toolbox 7.4.
% Generated on: 18-Jan-2018 16:52:25

% Equiripple Lowpass filter designed using the FIRPM function.

% All frequency values are in Hz.
Wpass = 1;      % Passband Weight
Wstop = 1;      % Stopband Weight
dens  = 20;     % Density Factor

% Calculate the coefficients using the FIRPM function.
%b  = firpm(N, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0], [Wpass Wstop], {dens});
b  = firls(N, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0], [Wpass Wstop]);

Hd = dfilt.dffir(b);

% [EOF]
end