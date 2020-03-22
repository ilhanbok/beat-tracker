%{
format long;

sig_start = 1 * 1000;
sig_end = 20 * 1000;
length = sig_end - sig_start + 1;
fs = 44100;
period = 1/fs;

[data, fs] = audioread('MapleLeafRag.ogg');
snippet = data(sig_start:sig_end);
data_fft = fft(snippet);

P2 = angle(data_fft/length);
P1 = P2(1:length/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = fs*(0:(length/2))/length;
plot(f,P1) 
title('Phase Spectrum of FFT (Maple Leaf Rag)')
xlabel('f (Hz)')
ylabel('Phase')
%}
format long;

frame_start = 2001;
frame_end = 4001;
length = frame_end - frame_start;

[data, fs] = audioread('rhythm_2.wav');

beat_start = 0;
beat_end = 1000;

b = zeros(1,beat_end - beat_start + 1);

for i = beat_start:beat_end
    data_fft = fft(data(frame_start:frame_end));
    P2 = angle(data_fft/length);
    P1 = P2(1:length/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = fs*(0:(length/2))/length;
    %plot(f,P1);
    xpad = ones(1,length/2 + 1);
    lsr = ([xpad; f]')\(P1');
    
    b(i - beat_start + 1) = lsr(2);
    
    frame_start = frame_start + 250;
    frame_end = frame_end + 250;
end

plot((beat_start:beat_end), b);