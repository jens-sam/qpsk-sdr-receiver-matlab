clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part I 
%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 1.0 load and analyze signal
% Figure 1: RF Spectrum

load('xRF1.mat');
% load('xRF2.mat');
% load('xRF3.mat');
% load('xRF4ans.mat');
% load('xRF5ans.mat');


figure;
spec_analysis(xRF, fs);
title('Spectrum of Received RF Signal xRF')
xlabel('Frequency (Hz)')
ylabel('Magnitude')
grid on

% 1.1

N = length(xRF); % Time vector
t = (0:N-1)'/fs;


xBB = xRF .* exp(-1j*2*pi*fc*t); % Downconvert


%Figure 2: Baseband Spectrum
figure;
spec_analysis(xBB, fs)
title('Spectrum after Downconversion (Baseband)')
grid on

pR = conj(flipud(pT));
xBB_filt = conv(xBB, pR);

% Figure 3: Filtered Spectrum
figure;
spec_analysis(xBB_filt, fs)
title('Spectrum after Matched Filtering')
grid on

%1.2 eye pattern

powers = zeros(L,1);

for k = 1:L
    samples = xBB_filt(k:L:end);
    powers(k) = sum(abs(samples).^2);
end

[~, best_k] = max(powers);

xBBd = xBB_filt(best_k:L:end);

% figure 4
figure;
plot(real(xBBd), imag(xBBd), '.')
title('Constellation after timing')
grid on


% 1.3 Find preamble location, then extract payload

cp_len = length(cp);
preamble = repmat(cp,4,1);   % known 128-symbol preamble


r = abs(conv(xBBd, flipud(conj(preamble)), 'valid')); % Correlate w/known preamble

[~, preamble_start] = max(r); % Find best match

payload_start = preamble_start + length(preamble); % Payload starts

payload = xBBd(payload_start:end); % Extract payload

% figure 5
figure;
plot(r)
title('Preamble correlation magnitude')
grid on


% 1.4

bits = QPSK2bits(payload);
bits(1:20)
length(payload)
length(bits)

% 1.5
bin2file(bits, 'recovered.txt')

