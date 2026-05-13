clear; clc; close all;

%%%%%%%%%%%%%%%%%%%
% Part II code with part I extended
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part I 
%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 1.0 load and analyze signal
% Figure 1: RF Spectrum

% load('xRF1.mat');
% load('xRF2.mat');
% load('xRF3.mat');
% load('xRF4ans.mat');
load('xRF5.mat');


% figure;
% spec_analysis(xRF, fs);
% title('Spectrum of Received RF Signal xRF')
% xlabel('Frequency (Hz)')
% ylabel('Magnitude')
% grid on

% 1.1

N = length(xRF); % Time vector
t = (0:N-1)'/fs;


xBB = xRF .* exp(-1j*2*pi*fc*t); % Downconvert


%Figure 2: Baseband Spectrum
% figure;
% spec_analysis(xBB, fs)
% title('Spectrum after Downconversion (Baseband)')
% grid on

pR = conj(flipud(pT));
xBB_filt = conv(xBB, pR);

% Figure 3: Filtered Spectrum
% figure;
% spec_analysis(xBB_filt, fs)
% title('Spectrum after Matched Filtering')
% grid on

%1.2 eye pattern

powers = zeros(L,1);

for k = 1:L
    samples = xBB_filt(k:L:end);
    powers(k) = sum(abs(samples).^2);
end

[~, best_k] = max(powers);

xBBd = xBB_filt(best_k:L:end);

% figure 4
% figure;
% plot(real(xBBd), imag(xBBd), '.')
% title('Constellation after timing')
% grid on


% 1.3 Find preamble location, then extract payload

cp_len = length(cp);
preamble = repmat(cp,4,1);   % known 128-symbol preamble


r = abs(conv(xBBd, flipud(conj(preamble)), 'valid')); % Correlate w/known preamble

[~, preamble_start] = max(r); % Find best match

payload_start = preamble_start + length(preamble); % Payload starts

payload = xBBd(payload_start:end); % Extract payload

% figure 5
% figure;
% plot(r)
% title('Preamble correlation magnitude')
% grid on



% 1.4

bits = QPSK2bits(payload);
bits(1:20)
length(payload)
length(bits)


%figure 6
% figure;
% stem(bits(1:50), 'filled')
% xlabel('Bit Index')
% ylabel('Bit Value')
% title('Recovered Bit Sequence (First 50 Bits)')

% 1.5
bin2file(bits, 'recovered.txt')


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part II 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.1 - detect preamble in xBBd
% using parts of code from chapter 11 problem 10

lag = length(cp);          % 32
ryy = zeros(length(xBBd),1);

for m = 2*lag:length(xBBd)
    acc = 0;
    for k = 0:lag-1
        acc = acc + xBBd(m-k) * conj(xBBd(m-lag-k));
    end
    ryy(m) = acc;
end

metric = abs(ryy);

beta = 0.6;
T = beta * max(metric);
aboveT = metric > T;

nstart = NaN;
for m = 1:length(aboveT)-lag+1
    if all(aboveT(m:m+lag-1))
        nstart = m;
        break;
    end
end

% figure 7
figure;
plot(metric,'LineWidth',1.2); hold on;
yline(T,'r--','LineWidth',1.2);
if ~isnan(nstart)
    xline(nstart,'k--','LineWidth',1.2);
end
grid on;
xlabel('Sample index n');
ylabel('|r_{yy}(n)|');
title('Part II Step 2.1: Running Autocorrelation for Preamble Detection');


% 2.2 - Cyclic equalizer adaptation 
% using code from chapter 11 problem 11

Neq = length(cp);          % equalizer / pilot length = 32
pilot = cp;                % known pilot sequence

% Extract one pilot cycle from detected preamble
if nstart + Neq - 1 <= length(xBBd)
    ypilot = xBBd(nstart : nstart + Neq - 1);
else
    error('Pilot extraction exceeds signal length');
end

% Build y0 = [y[n] y[n-1] ... y[n-N]]^T
y0 = flipud(ypilot(:));

% Initialize equalizer
w = zeros(Neq,1);
mu = 0.001;
Nit = 70000;

mse_hist = zeros(Nit,1);
yi = y0;

for ii = 1:Nit
    idx = mod(ii-1, Neq) + 1;      % desired pilot index

    e = pilot(idx) - w' * yi;      % error
    w = w + 2*mu*conj(e)*yi;       % LMS update

    mse_hist(ii) = abs(e)^2;

    % circular shift of input vector
    yi = [yi(end); yi(1:end-1)];
end

% Align largest tap to the middle
[~, idxMax] = max(abs(w));
center = ceil(Neq/2);
shift = center - idxMax;
w_aligned = circshift(w, shift);

% figure 8 and 9
figure;
semilogy(mse_hist,'LineWidth',1.2);
grid on;
xlabel('Iteration');
ylabel('|e[i]|^2');
title('Part II Step 2.2: Cyclic Equalizer Learning Curve');

figure;
stem(abs(w_aligned),'filled');
grid on;
xlabel('Tap index');
ylabel('|w[n]|');
title('Part II Step 2.2: Aligned Equalizer Tap Magnitudes');

% 2.3 - Apply equalizer to full signal 

xBBe = conv(xBBd, conj(w_aligned));


% figure 10
figure;
plot(real(xBBe), imag(xBBe), '.')
title('Part II Step 2.3: Equalized Constellation')
xlabel('In-phase')
ylabel('Quadrature')
grid on


% Eye diagram figure 11 and 12

% Reconstruct waveform (just like Chapter 11)
zUp = upsample(xBBe, L);
zWf = conv(zUp, pT);

Ns = 2*L;        % 2 symbols per trace
nTraces = 200;
start = 10*L;    % skip transient

% I (real part)
figure; hold on; grid on;
for k = 0:(nTraces-1)
    idx0 = start + k*L;
    if (idx0 + Ns - 1) <= length(zWf)
        plot(real(zWf(idx0:idx0+Ns-1)));
    end
end
xlabel('Samples (2 symbols)');
ylabel('Amplitude (I)');
title('Part II Step 2.3: Eye Diagram (I) After Equalization');

% Q (imag part)
figure; hold on; grid on;
for k = 0:(nTraces-1)
    idx0 = start + k*L;
    if (idx0 + Ns - 1) <= length(zWf)
        plot(imag(zWf(idx0:idx0+Ns-1)));
    end
end
xlabel('Samples (2 symbols)');
ylabel('Amplitude (Q)');
title('Part II Step 2.3: Eye Diagram (Q) After Equalization');


% 2.4 - Extract payload using equalized signal

Neq = length(cp);        % 32
Npre = 4*Neq;            % 128

[~, idxMax] = max(abs(w_aligned));
eqDelay = idxMax - 1;

payload_start = preamble_start + Npre + eqDelay;
payload_end   = payload_start + length(payload) - 1;

z_payload = xBBe(payload_start:payload_end);

bits = QPSK2bits(z_payload);

% bin2file(bits, 'recovered_eq_xRF1.txt');
% bin2file(bits, 'recovered_eq_xRF2.txt');
% bin2file(bits, 'recovered_eq_xRF3.txt');
% bin2file(bits, 'recovered_eq_xRF4.txt');
bin2file(bits, 'recovered_eq_xRF5.txt');
