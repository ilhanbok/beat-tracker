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

% Read the audio file
[data, Fs] = audioread('4321.wav');

% ///// ENERGY MEASURE ///// %

% Start and end markers for energy derivation frame
env_start = 0;
env_end = 2000;

len = env_end - env_start;

% How much to shift the frames
env_gap = 250;

% How many sample frames to take
frames = 1000;

% Previous energy value to determine slope
prev_energy = 0;

% Place to store the calculated slopes
slopes = zeros(1, frames);

% Keep shifting window and finding slope of energy
for i = 1:frames
    % Take the sum of squares of audio magnitude
    energy = sum(data(env_start + 1:env_end + 1).^2);
    slopes(i) = energy - prev_energy;
    prev_energy = energy;
    % Increment frame indicies to shift window
    env_start = env_start + env_gap;
    env_end = env_end + env_gap;
end

% Plot the results
figure(1);

plot(1:frames, slopes);

% Plot logistics
title('*** Energy Measure ***');
xlabel('Frame Number');
ylabel('Energy Derivative');

% ///// GROUP DELAY ///// %

% Start and end markers for FFT window (moving)
frame_start = 0;
frame_end = 2000;

len = frame_end - frame_start;

% How much to move the frame per iteration
frame_gap = 250;

% How many iterations during which the frame shifts
frames = 1000;

% Place to store the calculated slopes
slopes = zeros(1, frames);

% Keep shifting window and finding slope
for i = 1:frames
    % Take the DFT, extract angle, and store slope
    data_fft = fft(data(frame_start + 1:frame_end + 1));
    p = polyfit(0:len, angle(data_fft), 1);
    slopes(i) = p(1);
    % Increment frame indicies to shift window
    frame_start = frame_start + frame_gap;
    frame_end = frame_end + frame_gap;
end

% Plot the results
figure(2);

plot(1:frames, slopes);

% Plot logistics
title('*** Group Delay ***');
xlabel('Frame Number');
ylabel('Group Delay');

% ///// SPECTRAL CENTER ///// %

% Start and end markers for spectral center derivation frame
spec_start = 0;
spec_end = 2000;

len = spec_end - spec_start;

% How much to shift the frames
spec_gap = 250;

% How many sample frames to take
frames = 1000;

% Previous center value to determine slope
prev_center = 0;

% Place to store the calculated slopes
slopes = zeros(1, frames);

% Keep shifting window and finding slope of spectral mean
for i = 1:frames
    % Take the DFT to find the magnitude coeffs
    mag_fft = abs(fft(data(spec_start + 1:spec_end + 1)));
    % Find the sum, and iterate from the start to find halfway point
    mag_sum = sum(mag_fft);
    accum_mag = 0;
    center = 0;
    for n = 1:length(mag_fft)
        accum_mag = accum_mag + mag_fft(n);
        if accum_mag >= mag_sum/2
            center = n;
            break;
        end
    end
    % Take the derivative of the center
    slopes(i) = center - prev_center;
    prev_center = center;
    % Increment frame indicies to shift window
    spec_start = spec_start + spec_gap;
    spec_end = spec_end + spec_gap;
end

% Plot the results
figure(3);

plot(1:frames, slopes);

% Plot logistics
title('*** Spectral Center ***');
xlabel('Frame Number');
ylabel('Center Derivative');

%{
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
%}