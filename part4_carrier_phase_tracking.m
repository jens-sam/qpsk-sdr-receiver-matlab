% final for part 4, 1-3
clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part I 
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('xRF1.mat');
% load('xRF2.mat');
% load('xRF3.mat');
% load('xRF4.mat');
% load('xRF5ans.mat');
% load('xRF6.mat');
load('xRF7.mat');
% load('xRF8.mat');


N = length(xRF);
t = (0:N-1)'/fs;

% Downconvert
xBB = xRF .* exp(-1j*2*pi*fc*t);

% Matched filter
pR = conj(flipud(pT));
xBB_filt = conv(xBB, pR);

% Timing recovery
powers = zeros(L,1);
for k = 1:L
    samples = xBB_filt(k:L:end);
    powers(k) = sum(abs(samples).^2);
end

[~, best_k] = max(powers);
xBBd = xBB_filt(best_k:L:end);

figure;
plot(real(xBBd), imag(xBBd), '.');
title('Part IV 4.1: Constellation BEFORE Carrier Correction');
xlabel('In-phase'); ylabel('Quadrature');
grid on;
axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part IV Step 4.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tb_sym = Tb;                 % symbol period 
Npilot = length(cp);

% Known repeated preamble
preamble = repmat(cp,4,1);

% Detect preamble
r = abs(conv(xBBd, flipud(conj(preamble)), 'valid'));
[~, preamble_start] = max(r);

% Extract 4 pilot repetitions
pilot_rx = xBBd(preamble_start : preamble_start + 4*Npilot - 1);

% Compute J using adjacent periods
J1 = sum(pilot_rx(Npilot+1:2*Npilot) .* conj(pilot_rx(1:Npilot)));
J2 = sum(pilot_rx(2*Npilot+1:3*Npilot) .* conj(pilot_rx(Npilot+1:2*Npilot)));
J3 = sum(pilot_rx(3*Npilot+1:4*Npilot) .* conj(pilot_rx(2*Npilot+1:3*Npilot)));

J = J1 + J2 + J3;

Dfc_est = angle(J) / (2*pi*Npilot*Tb_sym);
fprintf('Estimated carrier offset = %.6f Hz\n', Dfc_est);

% Remove carrier offset
n = (0:length(xBBd)-1).';
xBBd_ca = xBBd .* exp(-1j*2*pi*Dfc_est*n*Tb_sym);


figure;
plot(real(xBBd_ca), imag(xBBd_ca), '.');
title('Part IV 4.1: Constellation AFTER Carrier Correction');
xlabel('In-phase'); ylabel('Quadrature');
grid on;
axis equal;




%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Continue with corrected signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Redo preamble detection with corrected signal
r = abs(conv(xBBd_ca, flipud(conj(preamble)), 'valid'));
[~, preamble_start] = max(r);

payload_start = preamble_start + length(preamble);
payload = xBBd_ca(payload_start:end);

bits = QPSK2bits(payload);
bin2file(bits, 'recovered.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part II
%%%%%%%%%%%%%%%%%%%%%%%%%%%

lag = length(cp);
ryy = zeros(length(xBBd_ca),1);

for m = 2*lag:length(xBBd_ca)
    acc = 0;
    for k = 0:lag-1
        acc = acc + xBBd_ca(m-k) * conj(xBBd_ca(m-lag-k));
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




% Equalizer
Neq = length(cp);
pilot = cp;

if nstart + Neq - 1 <= length(xBBd_ca)
    ypilot = xBBd_ca(nstart : nstart + Neq - 1);
else
    error('Pilot extraction exceeds signal length');
end

y0 = flipud(ypilot(:));

w = zeros(Neq,1);
mu = 0.001;
Nit = 100000;
mse_hist = zeros(Nit,1);
yi = y0;

for ii = 1:Nit
    idx = mod(ii-1, Neq) + 1;
    e = pilot(idx) - w' * yi;
    w = w + 2*mu*conj(e)*yi;
    mse_hist(ii) = abs(e)^2;
    yi = [yi(end); yi(1:end-1)];
end

[~, idxMax] = max(abs(w));
center = ceil(Neq/2);
shift = center - idxMax;
w_aligned = circshift(w, shift);

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


% Apply equalizer
xBBe = conv(xBBd_ca, conj(w_aligned));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part IV Step 4.2
% Decision-directed phase tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = zeros(size(xBBe));
xBBd_dd = zeros(size(xBBe));   

mu_dd = 0.001;

for k = 1:length(xBBe)-1
    xBBd_dd(k) = xBBe(k) * exp(-1j*phi(k));

    I = real(xBBd_dd(k));
    Q = imag(xBBd_dd(k));
    e = sign(I)*Q - sign(Q)*I;

    phi(k+1) = phi(k) + mu_dd * e;
end

xBBd_dd(end) = xBBe(end) * exp(-1j*phi(end));

figure;
plot(real(xBBe), imag(xBBe), '.');
title('Before Decision-Directed Phase Tracking');
xlabel('In-phase'); ylabel('Quadrature');
grid on; axis equal;

figure;
plot(real(xBBd_dd), imag(xBBd_dd), '.');
title('After Decision-Directed Phase Tracking');
xlabel('In-phase'); ylabel('Quadrature');
grid on; axis equal;

figure;
plot(phi);
title('Estimated Phase Correction');
xlabel('Sample index'); ylabel('\phi[n]');
grid on;


% figure 10
figure;
plot(real(xBBe), imag(xBBe), '.')
title('Part II Step 2.3: Equalized Constellation')
xlabel('In-phase')
ylabel('Quadrature')
grid on

% Extract payload

Neq = length(cp);        % 32
Npre = 4*Neq;            % 128

% Get payload length from corrected unequalized signal
payload_start_raw = preamble_start + length(preamble);
payload_raw = xBBd_ca(payload_start_raw:end);
payload_len = length(payload_raw);

% Original equalizer-delay logic
[~, idxMax] = max(abs(w_aligned));
eqDelay = idxMax - 1;

payload_start_eq = preamble_start + Npre + eqDelay;

% Make sure the end index does not exceed signal length
payload_len_final = min(payload_len, length(xBBd_dd) - payload_start_eq + 1);
payload_end_eq = payload_start_eq + payload_len_final - 1;

fprintf('length(xBBd_dd) = %d\n', length(xBBd_dd));
fprintf('payload_start_eq = %d\n', payload_start_eq);
fprintf('payload_end_eq = %d\n', payload_end_eq);

z_payload = xBBd_dd(payload_start_eq:payload_end_eq);

bits = QPSK2bits(z_payload);

% bin2file(bits, 'recovered_eq_p4_xRF1.txt');
% bin2file(bits, 'recovered_eq_p4_xRF2.txt');
% bin2file(bits, 'recovered_eq_p4_xRF3.txt');
% bin2file(bits, 'recovered_eq_p4_xRF4.txt');
% bin2file(bits, 'recovered_eq_p4_xRF5.txt');
% bin2file(bits, 'recovered_eq_p4_xRF6.txt');
bin2file(bits, 'recovered_eq_p4_xRF7.txt');
% bin2file(bits, 'recovered_eq_p4_xRF8.txt');


clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part IV: From Part III code RXQPSK2
%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part III: RXQPSK2
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose file
% load('xRF1.mat');
% load('xRF2.mat');
% load('xRF2ans.mat');
% load('xRF3.mat');
% load('xRF4.mat');
 load('xRF8.mat');

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

k0 = 1;              
M2 = L/2;            % decimation factor to get 2 samples/symbol

xBBd = xBB_filt(k0:M2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part IV Step 4.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tb_sym = Tb;
Npilot = length(cp);
preamble = repmat(cp,4,1);

% Build 2-sps preamble for carrier acquisition
preamble2 = zeros(2*length(preamble),1);
preamble2(1:2:end) = preamble;

% Detect preamble on 2-sps stream
r = abs(conv(xBBd, flipud(conj(preamble2)), 'valid'));
[~, preamble_start] = max(r);

% Extract repeated pilot region at 2 sps
pilot_rx = xBBd(preamble_start : preamble_start + 2*4*Npilot - 1);

% Correlate one pilot period apart (2N samples apart at 2 sps)
J1 = sum(pilot_rx(2*Npilot+1:4*Npilot) .* conj(pilot_rx(1:2*Npilot)));
J2 = sum(pilot_rx(4*Npilot+1:6*Npilot) .* conj(pilot_rx(2*Npilot+1:4*Npilot)));
J3 = sum(pilot_rx(6*Npilot+1:8*Npilot) .* conj(pilot_rx(4*Npilot+1:6*Npilot)));

J = J1 + J2 + J3;

Dfc_est = angle(J) / (2*pi*Npilot*Tb_sym);
fprintf('Estimated carrier offset = %.6f Hz\n', Dfc_est);

n = (0:length(xBBd)-1).';
xBBd_ca = xBBd .* exp(-1j*2*pi*Dfc_est*n*(Tb_sym/2));

figure;
plot(real(xBBd), imag(xBBd), '.');
title('Part IV 4.1: Constellation BEFORE Carrier Correction');
xlabel('In-phase'); ylabel('Quadrature');
grid on; axis equal;

figure;
plot(real(xBBd_ca), imag(xBBd_ca), '.');
title('Part IV 4.1: Constellation AFTER Carrier Correction');
xlabel('In-phase'); ylabel('Quadrature');
grid on; axis equal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Preamble identifier
%%%%%%%%%%%%%%%%%%%%%%%%%%%

cp_len = length(cp);             % 32 symbols
preamble = repmat(cp,4,1);       % 128-symbol preamble

% Create half-symbol-spaced version of preamble
preamble2 = zeros(2*length(preamble),1);
preamble2(1:2:end) = preamble;

% Correlate with known half-symbol-spaced preamble
r = abs(conv(xBBd_ca, flipud(conj(preamble2)), 'valid'));

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

pilot = cp;                   
Ns = length(pilot);           % 32 symbols
Neq = 2*length(cp);           % 64 taps, T/2 spacing
mu = 0.001;
Nit = 50000;

% one 2 sps pilot cycle
if preamble_start + Neq - 1 <= length(xBBd_ca)
    ypilot = xBBd_ca(preamble_start : preamble_start + Neq - 1);
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

xBBe = conv(xBBd_ca, conj(w_aligned), 'same');

figure;
plot(real(xBBe), imag(xBBe), '.');
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('Part III: Equalized Constellation (All Samples)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.5 Eye diagrams at equalizer output (2 sps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ns_eye = 4;              % 2 symbols at 2 sps
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
% Part IV Step 4.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = zeros(size(zSym));
zSym_dd = zeros(size(zSym));

mu_dd = 0.003;

for k = 1:length(zSym)-1
    zSym_dd(k) = zSym(k) * exp(-1j*phi(k));

    I = real(zSym_dd(k));
    Q = imag(zSym_dd(k));
    e = sign(I)*Q - sign(Q)*I;

    phi(k+1) = phi(k) + mu_dd * e;
end

zSym_dd(end) = zSym(end) * exp(-1j*phi(end));

figure;
plot(real(zSym), imag(zSym), '.');
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('Part IV 4.2: Before Decision-Directed Phase Tracking');

figure;
plot(real(zSym_dd), imag(zSym_dd), '.');
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('Part IV 4.2: After Decision-Directed Phase Tracking');

figure;
plot(phi);
grid on;
xlabel('Sample index');
ylabel('\phi[n]');
title('Part IV 4.2: Estimated Phase Correction');

figure;
plot(phi(end-500:end));
grid on;
title('Last 500 samples of phase estimate');


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
payload = zSym_dd(payload_start_sym:end);

figure;
plot(real(payload), imag(payload), '.');
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('Part III: Payload Constellation');

bits = QPSK2bits(payload);

% bin2file(bits, 'recovered_fse_p4_xRF1.txt');
% bin2file(bits, 'recovered_fse_p4_xRF2.txt');
% bin2file(bits, 'recovered_fse_p4_xRF3.txt');
% bin2file(bits, 'recovered_fse_p4_xRF4.txt');
% bin2file(bits, 'recovered_fse_p4_xRF5.txt');
% bin2file(bits, 'recovered_fse_p4_xRF6.txt');
% bin2file(bits, 'recovered_fse_p4_xRF7.txt');
bin2file(bits, 'recovered_fse_p4_xRF8.txt');




