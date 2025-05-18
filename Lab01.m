close all
clear all
% ????? ??? ???? ???
%msg = [1 0 1];
msg = round(rand(1,1000));
disp(msg)
% 1/2 rated convolutional Encoder
trellis = poly2trellis(3, [6 7]); % trellis is a structure 

user = convenc(msg, trellis);
disp('Encoded output: ')
disp(user)

% Convert binary to bipolar NRZ format
length_user = length(user);
disp('length_user')
disp(length_user)
for i = 1:length_user
    if user(i) == 0
        user(i) = -1;
    end
end
disp('user polar:')
disp(user)

fc = 5000;         % Carrier frequency in Hz (??? ????????????? ???? ??????? sine wave ?? ????????????)
eb = 0.5;          % Energy per bit (???? ??? ?????????? ???? ???? ????????? ?????, ??????? ????????????? purposes ?? ???? ??????? ??)
bitrate = 1000;    % Bit rate in bps (bits per second) ? ????? ???????? ???? ??? ?????? ?????
tb = 1/bitrate;    % Bit duration (seconds per bit) ? ???? ????? ?????? (bit period), ?????? ????? ??? ?? ??? ??? ?????? ???

chiprate = 10000;  % Chip rate in chips per second ? ??? spread spectrum ??????? ??, ????? ??? ????? ???????? ???? chip ?????? ?? ?? ?????
tc = 1/chiprate;   % Chip duration (seconds per chip) ? ???? chip ?? ??? ??? ?????? ??

% Time vector
t = tc:tc:tb*length_user;

% Baseband signal generation
basebandsig = []; %baseband vector
for i = 1:length_user
    for j = tc:tc:tb
        if user(i) == 1 
            basebandsig = [basebandsig 1];
        else 
            basebandsig = [basebandsig -1];
        end
    end
end

disp(basebandsig)

% === Figure 1: Baseband Signal ===
figure(1)
stairs(t(1:800), basebandsig(1:800), 'r') % red
xlabel('Time (sec)')
ylabel('Binary value')
set(gca, 'ytick', [-1 1])
axis([0 max(t(1:800)) -1.5 1.5]);
title('A segment of original binary sequence for a single user')
grid on

% BPSK Modulation
bpskmod = []; % vector
for i = 1:length_user
    for j = tc:tc:tb
        bpskmod = [bpskmod sqrt(2*eb)*user(i)*cos(2*pi*fc*j)];
    end
end

disp('bpsk output: ')
disp(bpskmod)

% Frequency Domain Analysis
number = length(t);
spectrum = abs(fft(bpskmod)); % convert time domain signal to frequency domain
sampling_frequency = 2 * fc;
sampling_interval = 1.0 / sampling_frequency;
nyquest_frequency = 1.0 / (2.0 * sampling_interval);

for i = 1:number
    frequency(i) = (1.0 / (number * sampling_interval)) * i;
end

% === Figure 2: Frequency Domain ===
figure(2)
plot(frequency, spectrum, 'b') % blue
title('Frequency Domain analysis of BPSK modulated signal for a single user')

xlabel('Frequency (Hz)')
ylabel('Magnitude')
grid on

% PN sequence generation for spreading
seed = [1 -1 1 -1];
pn = [];
for i = 1:length_user
    for j = 1:10
        pn = [pn seed(4)];
        if seed(4) == seed(3)
            temp = -1;
        else
            temp = 1;
        end
        seed(4) = seed(3);
        seed(3) = seed(2);
        seed(2) = seed(1);
        seed(1) = temp;
    end
end
% same the size
% Upsample PN sequence
pnupsampled = [];
len_pn = length(pn);
for i = 1:len_pn
    for j = 10*tc:10*tc:tb
        if pn(i) == 1 
            pnupsampled = [pnupsampled 1];
        else 
            pnupsampled = [pnupsampled -1];
        end
    end
end

% Multiply BPSK with PN sequence (Spreading)
sigtx = bpskmod .* pnupsampled;

% === Figure 3: Transmitted Signal ===
figure(3)
plot(t(1:200), sigtx(1:200), 'g') % green
title('A segment of Transmitted DS CDMA signal')
axis([0 max(t(1:200)) -1.5 1.5]);
xlabel('Time (sec)')
ylabel('Amplitude')
grid on

% BER Simulation under AWGN
snr_in_dBs = 0:1.0:10;
for m = 1:length(snr_in_dBs)
    ber(m) = 0.0;

    % Add AWGN noise
    composite_signal = awgn(sigtx, snr_in_dBs(m), 'measured');

    % Receiver multiplies with same PN
    rx = composite_signal .* pnupsampled;

    % BPSK Demodulation
    demodcar = [];
    for i = 1:length_user
        for j = tc:tc:tb
            demodcar = [demodcar sqrt(2*eb)*cos(2*pi*fc*j)];
        end
    end

    bpskdemod = rx .* demodcar;

    % Integrate and dump
    len_dmod = length(bpskdemod);
    sum = zeros(1, len_dmod/10);
    for i = 1:len_dmod/10
        for j = (i-1)*10 + 1:i*10
            sum(i) = sum(i) + bpskdemod(j);
        end
    end

    % Hard decision
    rxbits = [];
    for i = 1:length_user
        if sum(i) > 0
            rxbits = [rxbits 1];
        else
            rxbits = [rxbits 0];
        end
    end

    % Viterbi decoding
    tblen = 3;
    delay = tblen;
    decoded = vitdec(rxbits, trellis, tblen, 'cont', 'hard');

    % BER Calculation
    [~, rat] = biterr(decoded(delay+1:end), msg(1:end-delay));
    ber(m) = rat;
end

% === Figure 4: BER vs SNR ===
figure(4)
plot(snr_in_dBs, ber, 'm', 'LineWidth', 2) % magenta
xlabel('Signal to noise ratio (dB)')  
ylabel('BER')
legend('BER simulation for a single user')
title('Coded BER simulation under AWGN channel')
grid on




