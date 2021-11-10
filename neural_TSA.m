%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neural Time Series Analysis
% 
% Syntax:
% 
% 
% Inputs:
% 
% 
% Output:
% 
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2021  Sina Dabiri
% sdabiri@emory.edu
% 
% 
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A) Making and summing sine waves, Fourier analysis, power structure, compare nonzero freq with sine waves

fs = 1000;
T = 1/fs;
L = 1000;
t = (0:L-1)*T;

% Creating 5 sine waves and suming them to create signal 6
s1 = sin(2*pi*8*t);
s2 = 0.8*sin(2*pi*13*t);
s3 = 0.7*sin(2*pi*22*t);
s4 = 0.6*sin(2*pi*35*t);
s5 = 0.5*sin(2*pi*80*t);
s6 = s1+s2+s3+s4+s5;


% FFT on signal 6
data_f_two_sided = fft(s6);
d_two_sided = abs(data_f_two_sided/L);
half_L = round(L/2)+1;
d_one_side = d_two_sided(1:half_L);
d_one_side(2:end-1) = 2*d_one_side(2:end-1);
f = fs*(0:(L/2)) /L;

% Calculating the power of the signal
p = abs(f).^2/L;

%% 

% load dataset.mat eeg
% fs = 1000;
% [~, L] = size(eeg); %length of the signal
% 
% data_f_two_sided = fft(eeg);
% d_two_sided = abs(data_f_two_sided/L);
% half_L = round(L/2)+1;
% d_one_side = d_two_sided(1:half_L);
% d_one_side(2:end-1) = 2*d_one_side(2:end-1);
% f = fs*(0:(L/2)) /L;
%% Plots

% The separate sine plots
grid on;
figure(1)
xlabel('Time (ms)')
ylabel('Amplitude')
subplot(5,1,1);
plot(t,s1);
subplot(5,1,2);
plot(t,s2);
subplot(5,1,3);
plot(t,s3);
subplot(5,1,4);
plot(t,s4);
subplot(5,1,5);
plot(t,s5);

% Time
figure(2)
plot(t,s6);
xlabel('Time (ms)')
ylabel('Amplitude')

% Freq
figure(3)
% periodogram(data_base_cor)]
plot(f(1:100), d_one_side(1:100))
title('Single sided amplitute spectrum of the sum of the sins signal')
xlabel('f(Hz)')
ylabel('|s6|');

% Power
figure(4)
plot(f(1:100), p(1:100))
title('Single sided amplitute spectrum of the sum of the sins signal')
xlabel('f(Hz)')
ylabel('Power');