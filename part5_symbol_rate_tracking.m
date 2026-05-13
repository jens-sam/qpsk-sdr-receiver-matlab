
clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part V: Symbol Rate Tracking
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load('xRF1.mat');
load('xRF9.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Downconvert and match filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = length(xRF);
t = (0:N-1)'/fs;

xBB = xRF .* exp(-1j*2*pi*fc*t);
pR = conj(flipud(pT));
xBB = conv(xBB, pR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coarse timing and payload start
%%%%%%%%%%%%%%%%%%%%%%%%%%%

powers = zeros(L,1);
for k = 1:L
    powers(k) = sum(abs(xBB(k:L:end)).^2);
end
[~, best_k] = max(powers);

xBBd = xBB(best_k:L:end);

cp = cp(:);
preamble = repmat(cp,4,1);

r = abs(conv(xBBd(:), flipud(conj(preamble)), 'valid'));
[~, preamble_start] = max(r);

payload_start_sym = preamble_start + length(preamble);
Npayload_sym = length(xBBd(payload_start_sym:end));

payload_start_sample = best_k + (payload_start_sym - 1)*L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract payload block at sample rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%

front_margin_sym = 4;
front_margin_samp = front_margin_sym * L;

payload_len_samp = Npayload_sym * L;

block_start = max(1, payload_start_sample - front_margin_samp);
block_end = min(length(xBB), payload_start_sample + payload_len_samp - 1);

y = xBB(block_start:block_end);
y = [y; zeros(4*L,1)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decision-directed timing loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%

mu = 0.05;
dtau = 6;

start = 1 + front_margin_samp;
Ly = length(y);

tau = zeros(Npayload_sym+1,1);
z = zeros(Npayload_sym,1);

kk = 1;
for k = start:L:(start + (Npayload_sym-1)*L)

    tauSamp = tau(kk)*L;

    n_now   = k + tauSamp;
    n_late  = n_now + dtau;
    n_early = n_now - dtau;

    i_now   = floor(n_now);
    i_late  = floor(n_late);
    i_early = floor(n_early);

    a_now   = n_now   - i_now;
    a_late  = n_late  - i_late;
    a_early = n_early - i_early;

    if i_now < 1 || i_now+1 > Ly || ... % bounds check
       i_late < 1 || i_late+1 > Ly || ...
       i_early < 1 || i_early+1 > Ly
        break;
    end

    y_now = (1-a_now)*y(i_now) + a_now*y(i_now+1);
    y_late = (1-a_late)*y(i_late) + a_late*y(i_late+1);
    y_early = (1-a_early)*y(i_early) + a_early*y(i_early+1);

    sk = sign(real(y_now)) + 1j*sign(imag(y_now));

    z(kk) = y_now;

    tau(kk+1) = tau(kk) + mu * real(conj(sk - y_now) * (y_late - y_early));

    kk = kk + 1;
end

z = z(1:kk-1);
tau_hist = tau(1:kk-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decode
%%%%%%%%%%%%%%%%%%%%%%%%%%%

bits = QPSK2bits(z);
bin2file(bits, 'recovered_part5_xRF9.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
plot(tau_hist)
title('Tracked Timing Phase')
xlabel('Symbol Index')
ylabel('\tau[n]')
grid on

figure;
stem(bits(1:min(50,length(bits))), 'filled')
xlabel('Bit Index')
ylabel('Bit Value')
title('Recovered Bit Sequence (First 50 Bits)')




