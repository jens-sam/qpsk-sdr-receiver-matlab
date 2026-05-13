clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part VI: RXQPSK_part6.m
% Use for xRF7-xRF10
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% Load file and settings
% 
file_tag = 'xRF7';
load('xRF7.mat');

mu_train = 0.001;
mu_dd    = 0.0001;
mu_phase = 0.002;
mu_f     = 5e-5;

fprintf('Running %s\n', file_tag);

% 
% Downconvert
% 
Nrf = length(xRF);
t = (0:Nrf-1)'/fs;
xBB = xRF .* exp(-1j*2*pi*fc*t);

% 
% Matched filter
% 
pR = conj(flipud(pT));
xBB_filt = conv(xBB, pR);

% 
% Symbol-rate timing selection
% 
powers = zeros(L,1);

for k = 1:L
    samples = xBB_filt(k:L:end);
    powers(k) = sum(abs(samples).^2);
end

[~, best_k] = max(powers);
xBBd = xBB_filt(best_k:L:end);

fprintf('best_k = %d\n', best_k);

% 
% Preamble detection
% 
cp_len = length(cp);
Neq = cp_len;
Npre = 4*Neq;

preamble = repmat(cp,4,1);

r = abs(conv(xBBd, flipud(conj(preamble)), 'valid'));
[~, preamble_start] = max(r);

lag = cp_len;
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

if isnan(nstart)
    error('Preamble not detected.');
end

fprintf('preamble_start (corr) = %d\n', preamble_start);
fprintf('nstart (autocorr)     = %d\n', nstart);

figure;
plot(r,'LineWidth',1.2);
grid on;
xlabel('Symbol index');
ylabel('|r[n]|');
title([file_tag ' Preamble Correlation']);

figure;
plot(metric,'LineWidth',1.2); hold on;
yline(T,'r--','LineWidth',1.2);
xline(nstart,'k--','LineWidth',1.2);
grid on;
xlabel('Symbol index');
ylabel('|r_{yy}(n)|');
title([file_tag ' Running Autocorrelation']);

% 
% Coarse CFO estimate from repeated pilot
% 
if nstart + Npre - 1 > length(xBBd)
    error('Preamble exceeds signal length.');
end

zpre = xBBd(nstart:nstart+Npre-1);

J = 0;
for m = 1:3
    blk1 = zpre((m-1)*Neq + (1:Neq));
    blk2 = zpre(m*Neq     + (1:Neq));
    J = J + sum(blk2 .* conj(blk1));
end

Delta_f = angle(J) / (2*pi*Neq*Tb);
fprintf('Estimated coarse CFO = %.6f Hz\n', Delta_f);

n_sym = (0:length(xBBd)-1).';
xBBd_cfo = xBBd .* exp(-1j*2*pi*Delta_f*Tb*n_sym);

% 
% Initial cyclic equalizer training on one pilot cycle
% 
pilot = cp(:);

if nstart + Neq - 1 > length(xBBd_cfo)
    error('Pilot extraction exceeds signal length.');
end

ypilot = xBBd_cfo(nstart:nstart+Neq-1);
y0 = flipud(ypilot(:));

w = zeros(Neq,1);
Nit = 70000;

mse_hist = zeros(Nit,1);
yi = y0;

for ii = 1:Nit
    idx = mod(ii-1, Neq) + 1;

    z = w' * yi;
    e = pilot(idx) - z;
    w = w + 2*mu_train*conj(e)*yi;

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
title([file_tag ' Equalizer Learning Curve']);

figure;
stem(abs(w_aligned),'filled');
grid on;
xlabel('Tap index');
ylabel('|w[n]|');
title([file_tag ' Initial Aligned Equalizer Taps']);

% 
% Equalize full stream
% 
N = length(xBBd_cfo);
xBBe = zeros(N,1);

yi = zeros(Neq,1);
err_hist = zeros(N,1);
tap_idx_hist = zeros(N,1);

for n = 1:N

    yi = [xBBd_cfo(n); yi(1:end-1)];

    z = w_aligned' * yi;
    xBBe(n) = z;

    if n >= nstart && n < nstart + Npre
        idx = mod(n - nstart, Neq) + 1;
        d = cp(idx);
    else
        d = (2*(real(z)>=0)-1) + 1j*(2*(imag(z)>=0)-1);
    end

    e = d - z;
    err_hist(n) = abs(e)^2;

    % tiny adaptation only during preamble
    if n < nstart + Npre
        w_aligned = w_aligned + 2*mu_dd*conj(e)*yi;
    end

    [~, idxMax] = max(abs(w_aligned));
    tap_idx_hist(n) = idxMax;

    shift = center - idxMax;
    w_aligned = circshift(w_aligned, shift);
end

figure;
semilogy(err_hist,'LineWidth',1.2);
grid on;
xlabel('Symbol index');
ylabel('|e[n]|^2');
title([file_tag ' Equalizer Error']);

figure;
plot(tap_idx_hist,'LineWidth',1.2);
grid on;
xlabel('Symbol index');
ylabel('Largest tap index');
title([file_tag ' Largest Tap Position']);

% 
% DD phase/frequency tracking
% 
phi = 0;
freq = 0;

z_corr = zeros(size(xBBe));
phi_hist = zeros(length(xBBe),1);
freq_hist = zeros(length(xBBe),1);

for n = 1:length(xBBe)

    z_rot = xBBe(n) * exp(-1j*phi);
    z_corr(n) = z_rot;

    %a_hat = (2*(real(z_rot)>=0)-1) + 1j*(2*(imag(z_rot)>=0)-1);
    a_hat = sign(real(z_rot)) + 1j*sign(imag(z_rot));

    %e_phase = sign(real(z_rot))*imag(z_rot) - sign(imag(z_rot))*real(z_rot);
    e_phase = imag(z_rot * conj(a_hat));

    freq = freq + mu_f * e_phase;
    phi  = phi + freq + mu_phase * e_phase;

    phi_hist(n) = phi;
    freq_hist(n) = freq;
end

figure;
plot(phi_hist,'LineWidth',1.2);
grid on;
xlabel('Symbol index');
ylabel('\phi[n]');
title([file_tag ' Phase Tracking']);

figure;
plot(freq_hist,'LineWidth',1.2);
grid on;
xlabel('Symbol index');
ylabel('freq[n]');
title([file_tag ' Frequency Tracking']);

% 
% Extract payload
% 
payload_start = nstart + Npre;

if payload_start < 1 || payload_start > length(z_corr)
    error('Payload start exceeds equalized stream length.');
end

payload_end = length(z_corr);
z_payload = z_corr(payload_start:payload_end);

figure;
plot(real(z_payload), imag(z_payload), '.');
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title([file_tag ' Payload Constellation']);

%
% Final
%
bits = QPSK2bits(z_payload);

outname = ['recovered_part6_' file_tag '.txt'];
bin2file(bits, outname);

disp(bits(1:min(20,length(bits))));
fprintf('Recovered %d bits\n', length(bits));
fprintf('Saved %s\n', outname);


clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part VI: RXQPSK2
% Half-symbol-spaced adaptive equalizer
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose file
load('xRF7.mat');
% load('xRF8.mat');
% load('xRF9.mat');
% load('xRF10.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Downconvert and matched filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(xRF);
t = (0:N-1)'/fs;

xBB = xRF .* exp(-1j*2*pi*fc*t);

pR = conj(flipud(pT));
xBB_filt = conv(xBB, pR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Decimate to 2 samples/symbol
%%%%%%%%%%%%%%%%%%%%%%%%%%%

k0 = 1;                  % fixed timing phase
M2 = L/2;                % decimate to 2 sps
xBBd = xBB_filt(k0:M2:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Detect cyclic preamble
%%%%%%%%%%%%%%%%%%%%%%%%%%%

cp_len = length(cp);             % 32
preamble = repmat(cp,4,1);       % 128 symbols

% zero-stuffed preamble at 2 sps
preamble2 = zeros(2*length(preamble),1);
preamble2(1:2:end) = preamble;

r = abs(conv(xBBd, flipud(conj(preamble2)), 'valid'));
[~, preamble_start] = max(r);

figure;
plot(r,'LineWidth',1.2);
grid on;
xlabel('Sample index');
ylabel('Correlation magnitude');
title('Part VI: Half-Symbol-Spaced Preamble Correlation');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. Initial cyclic equalizer training on preamble
%%%%%%%%%%%%%%%%%%%%%%%%%%%

pilot = cp;                    % symbol-spaced desired pilot
Ns   = length(pilot);          % 32 symbols
Neq  = 2*length(cp);           % 64 taps, T/2 spaced
mu0  = 0.0001;                  % initial training step size
Nit0 = 100000;                  % initial training iterations

if preamble_start + Neq - 1 <= length(xBBd)
    ypilot = xBBd(preamble_start : preamble_start + Neq - 1);
else
    error('Pilot extraction exceeds signal length.');
end

y0 = flipud(ypilot(:));
w = zeros(Neq,1);
yi = y0;
mse_hist0 = zeros(Nit0,1);

for ii = 1:Nit0
    idx = mod(ii-1, Ns) + 1;
    d = pilot(idx);

    yhat = w' * yi;
    e = d - yhat;

    w = w + 2*mu0*conj(e)*yi;
    mse_hist0(ii) = abs(e)^2;

    yi = circshift(yi, 2);     % advance by 1 symbol = 2 samples
end

% align largest tap to center
[~, idxMax] = max(abs(w));
center = ceil(Neq/2);
shift = center - idxMax;
w = circshift(w, shift);

figure;
semilogy(mse_hist0,'LineWidth',1.2);
grid on;
xlabel('Iteration');
ylabel('|e[i]|^2');
title('Part VI: Initial Half-Symbol Equalizer Training Curve');

figure;
stem(abs(w),'filled');
grid on;
xlabel('Tap index');
ylabel('|w[n]|');
title('Part VI: Initial Equalizer Tap Magnitudes');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Adaptive equalization over full packet
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Equalizer adaptation is used to track slow carrier/timing drift

mu_dd = 0.00001;               
recenter_thresh = 2;         
freeze_len = 2*length(preamble);  

xBBe = zeros(size(xBBd));
w_hist = zeros(Neq,length(xBBd));   
err_hist = zeros(length(xBBd),1);

for n = Neq:2:length(xBBd)

    yi = flipud(xBBd(n-Neq+1:n));
    yhat = w' * yi;
    xBBe(n) = yhat;

    % Save current taps
    w_hist(:,n) = w;

    % Decide desired symbol
    sym_idx = (n+1)/2;   

    if n >= preamble_start && n < preamble_start + 2*length(preamble)
        rel = n - preamble_start;
        ksym = floor(rel/2) + 1;

        if ksym >= 1 && ksym <= length(preamble)
            d = preamble(ksym);
        else
            d = sign(real(yhat)) + 1j*sign(imag(yhat));
        end
    else
        
        d = sign(real(yhat)) + 1j*sign(imag(yhat));
        if real(d) == 0
            d = d + 1;
        end
        if imag(d) == 0
            d = real(d) + 1j;
        end
    end

    e = d - yhat;
    err_hist(n) = abs(e)^2;

 
    if n > preamble_start + freeze_len
        w = w + 2*mu_dd*conj(e)*yi;
    end

    % recenter if largest tap drifts too far
    [~, idxMax] = max(abs(w));
    if abs(idxMax - center) > recenter_thresh
        shift = center - idxMax;
        w = circshift(w, shift);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. Fill odd samples for plotting only
%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n = Neq+1:2:length(xBBd)-1
    yi = flipud(xBBd(n-Neq+1:n));
    xBBe(n) = w' * yi;
end

figure;
plot(real(xBBe), imag(xBBe), '.');
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('Part VI: Adaptive Equalized Constellation (All 2-sps Samples)');

figure;
semilogy(err_hist,'LineWidth',1.2);
grid on;
xlabel('2-sps sample index');
ylabel('|e[n]|^2');
title('Part VI: Equalizer Error During Packet');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 7. Eye diagrams at equalizer output
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Ns_eye = 4;          % 2 symbols at 2 sps
% nTraces = 150;
% start_eye = max(Neq+20, preamble_start);
% 
% figure; hold on; grid on;
% for k = 0:(nTraces-1)
%     idx0 = start_eye + 2*k;
%     if idx0 + Ns_eye - 1 <= length(xBBe)
%         plot(real(xBBe(idx0:idx0+Ns_eye-1)));
%     end
% end
% xlabel('Samples (2 symbols at 2 sps)');
% ylabel('Amplitude (I)');
% title('Part VI: Eye Diagram (I) at Adaptive Equalizer Output');
% 
% figure; hold on; grid on;
% for k = 0:(nTraces-1)
%     idx0 = start_eye + 2*k;
%     if idx0 + Ns_eye - 1 <= length(xBBe)
%         plot(imag(xBBe(idx0:idx0+Ns_eye-1)));
%     end
% end
% xlabel('Samples (2 symbols at 2 sps)');
% ylabel('Amplitude (Q)');
% title('Part VI: Eye Diagram (Q) at Adaptive Equalizer Output');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. Pick symbol-spaced branch
%%%%%%%%%%%%%%%%%%%%%%%%%%%

z1 = xBBe(1:2:end);
z2 = xBBe(2:2:end);

E1 = mean(abs(z1).^2);
E2 = mean(abs(z2).^2);

if E2 > E1
    zSym = z2;
    preamble_start_sym = ceil((preamble_start-1)/2);
    chosen_branch = 2;
else
    zSym = z1;
    preamble_start_sym = ceil(preamble_start/2);
    chosen_branch = 1;
end

disp(['Chosen symbol branch = ', num2str(chosen_branch)]);

figure;
plot(real(zSym), imag(zSym), '.');
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('Part VI: Symbol-Spaced Samples After Adaptive Half-Symbol Equalizer');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 9. Re-detect preamble on zSym
%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_sym = abs(conv(zSym, flipud(conj(preamble)), 'valid'));
[~, pstart_check] = max(r_sym);

figure;
plot(r_sym,'LineWidth',1.2);
grid on;
xlabel('Symbol index');
ylabel('Correlation magnitude');
title('Part VI: Symbol-Spaced Preamble Correlation After Equalization');

disp(['Initial preamble_start_sym = ', num2str(preamble_start_sym)]);
disp(['Re-detected preamble_start_sym = ', num2str(pstart_check)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10. Extract payload
%%%%%%%%%%%%%%%%%%%%%%%%%%%

payload_start_sym = pstart_check + length(preamble);

if payload_start_sym > length(zSym)
    error('Payload start exceeds symbol sequence length.');
end

payload = zSym(payload_start_sym:end);

figure;
plot(real(payload), imag(payload), '.');
grid on;
xlabel('In-phase');
ylabel('Quadrature');
title('Part VI: Payload Constellation');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 11. Recover bits
%%%%%%%%%%%%%%%%%%%%%%%%%%%

bits = QPSK2bits(payload);
disp(['Recovered ', num2str(length(bits)), ' bits']);

% choose one output file at a time

bin2file(bits, 'recovered_part6_fse_xRF7.txt');
% bin2file(bits, 'recovered_part6_fse_xRF8.txt');
% bin2file(bits, 'recovered_part6_fse_xRF9.txt');
% bin2file(bits, 'recovered_part6_fse_xRF10.txt');



