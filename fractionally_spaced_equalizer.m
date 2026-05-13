clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part III: RXQPSK2
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose file
% load('xRF1.mat');
% load('xRF2.mat');
% load('xRF2ans.mat');
% load('xRF3.mat');
% load('xRF4.mat');
load('xRF5.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Downconvert and matched filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(xRF);
t = (0:N-1)'/fs;

xBB = xRF .* exp(-1j*2*pi*fc*t);     % Downconvert to baseband

pR = conj(flipud(pT));
xBB_filt = conv(xBB, pR);            % Matched filter output

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Decimate to 2 samples/symbol
%%%%%%%%%%%%%%%%%%%%%%%%%%%

k0 = 1;              % timing phase
M2 = L/2;            % decimation factor 

xBBd = xBB_filt(k0:M2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Preamble identifier
%%%%%%%%%%%%%%%%%%%%%%%%%%%

cp_len = length(cp);             % 32 symbols
preamble = repmat(cp,4,1);       % 128-symbol preamble

% Create half-symbol-spaced version of preamble
preamble2 = zeros(2*length(preamble),1);
preamble2(1:2:end) = preamble;

% Correlate with known half-symbol-spaced preamble
r = abs(conv(xBBd, flipud(conj(preamble2)), 'valid'));

[~, preamble_start] = max(r);

figure;
plot(r,'LineWidth',1.2);
grid on;
xlabel('Sample index');
ylabel('Correlation magnitude');
title('Part III: Half-Symbol-Spaced Preamble Correlation');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Half-symbol-spaced cyclic equalizer
%%%%%%%%%%%%%%%%%%%%%%%%%%%

pilot = cp;                    % desired symbols stay symbol-spaced
Ns = length(pilot);           % 32 symbols
Neq = 2*length(cp);           % 64 taps, T/2 spacing
mu = 0.01;
Nit = 7000;

% one 2-sps pilot cycle
if preamble_start + Neq - 1 <= length(xBBd)
    ypilot = xBBd(preamble_start : preamble_start + Neq - 1);
else
    error('Pilot extraction exceeds signal length');
end

y0 = flipud(ypilot(:));

w = zeros(Neq,1);
yi = y0;
mse_hist = zeros(Nit,1);

for ii = 1:Nit
    idx = mod(ii-1, Ns) + 1;      % desired symbol index: 1...32

    e = pilot(idx) - w' * yi;     
    w = w + 2*mu*conj(e)*yi;
    mse_hist(ii) = abs(e)^2;

    yi = circshift(yi, 2);        % advance by 2 samples = 1 symbol
end

[~, idxMax] = max(abs(w));
center = ceil(Neq/2);             % 32
shift = center - idxMax;
w_aligned = circshift(w, shift);


figure;
semilogy(mse_hist,'LineWidth',1.2);
grid on;
xlabel('Iteration');
ylabel('|e[i]|^2');
title('Part III: Half-Symbol-Spaced Equalizer Learning Curve');

figure;
stem(abs(w_aligned),'filled');
grid on;
xlabel('Tap index');
ylabel('|w[n]|');
title('Part III: Half-Symbol-Spaced Equalizer Tap Magnitudes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Apply equalizer to full signal at 2 sps
%%%%%%%%%%%%%%%%%%%%%%%%%%%

xBBe = conv(xBBd, conj(w_aligned), 'same');

figure;
plot(real(xBBe), imag(xBBe), '.');
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('Part III: Equalized Constellation (All Samples)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.5 Eye diagrams at equalizer output (2 sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ns_eye = 4;              
nTraces = 150;
start_eye = 20;

figure; hold on; grid on;
for k = 0:(nTraces-1)
    idx0 = start_eye + 2*k;
    if idx0 + Ns_eye - 1 <= length(xBBe)
        plot(real(xBBe(idx0:idx0+Ns_eye-1)));
    end
end
xlabel('Samples (2 symbols at 2 sps)');
ylabel('Amplitude (I)');
title('Part III: Eye Diagram (I) at Equalizer Output');

figure; hold on; grid on;
for k = 0:(nTraces-1)
    idx0 = start_eye + 2*k;
    if idx0 + Ns_eye - 1 <= length(xBBe)
        plot(imag(xBBe(idx0:idx0+Ns_eye-1)));
    end
end
xlabel('Samples (2 symbols at 2 sps)');
ylabel('Amplitude (Q)');
title('Part III: Eye Diagram (Q) at Equalizer Output');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Choose symbol-spaced samples from equalized output
%%%%%%%%%%%%%%%%%%%%%%%%%%%

z1 = xBBe(1:2:end);
z2 = xBBe(2:2:end);

E1 = mean(abs(z1).^2);
E2 = mean(abs(z2).^2);

if E2 > E1
    zSym = z2;
    preamble_start_sym = ceil((preamble_start-1)/2);
else
    zSym = z1;
    preamble_start_sym = ceil(preamble_start/2);
end


figure;
plot(real(zSym), imag(zSym), '.');
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('Part III: Symbol Decisions After Half-Symbol Equalizer');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. Eye diagrams after symbol selection
%%%%%%%%%%%%%%%%%%%%%%%%%%%

zUp = upsample(zSym, L);
zWf = conv(zUp, pT);

Ns_eye = 2*L;
nTraces = 200;
start_eye = 10*L;

figure; hold on; grid on;
for k = 0:(nTraces-1)
    idx0 = start_eye + k*L;
    if (idx0 + Ns_eye - 1) <= length(zWf)
        plot(real(zWf(idx0:idx0+Ns_eye-1)));
    end
end
xlabel('Samples (2 symbols)');
ylabel('Amplitude (I)');
title('Part III: Eye Diagram (I) After Symbol Selection');

figure; hold on; grid on;
for k = 0:(nTraces-1)
    idx0 = start_eye + k*L;
    if (idx0 + Ns_eye - 1) <= length(zWf)
        plot(imag(zWf(idx0:idx0+Ns_eye-1)));
    end
end
xlabel('Samples (2 symbols)');
ylabel('Amplitude (Q)');
title('Part III: Eye Diagram (Q) After Symbol Selection');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. Extract payload and recover bits
%%%%%%%%%%%%%%%%%%%%%%%%%%%

payload_start_sym = preamble_start_sym + length(preamble);
payload = zSym(payload_start_sym:end);

figure;
plot(real(payload), imag(payload), '.');
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('Part III: Payload Constellation');

bits = QPSK2bits(payload);
% bin2file(bits, 'recovered_fse_part3_xRF1.txt');
% bin2file(bits, 'recovered_fse_part3_xRF2.txt');
% bin2file(bits, 'recovered_fse_part3_xRF3.txt');
%bin2file(bits, 'recovered_fse_part3_xRF4.txt');
bin2file(bits, 'recovered_fse_part3_xRF5.txt');



