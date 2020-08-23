% Attempt to track the beat of a song using an adaptive wavetable
% oscillator indexing into a cosine wavetable

% Ilhan Bok, 2020
% CC By Attribution

format long;

% Read the audio file
[data, Fs] = audioread('MapleLeafRag.ogg');

% How long to sample
num_chunks = 10000;

% How much to shift the frames
env_gap = 250;

% Play audio to hear sample
player = audioplayer(data(1:num_chunks*env_gap), Fs);
%play(player);

% Apply SPECTRAL CENTER

% Start and end markers for energy derivation frame
env_start = 0;
env_end = 2000;

len = env_end - env_start;

% How many sample frames to take
frames = num_chunks;

% Previous center value to determine slope
prev_center = 0;

% Previous dispersion value to determine slope
prev_disp = 0;

% Previous energy value to determine slope
prev_energy = 0;

% Place to store the calculated slopes
slopes = zeros(1, frames);

% Keep shifting window and finding slope of spectral mean
for i = 1:frames
    % Take the DFT to find the magnitude coeffs
    mag_fft = abs(fft(data(env_start + 1:env_end + 1)));
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
    env_start = env_start + env_gap;
    env_end = env_end + env_gap;
end

data = slopes;

% normalize values!
data = data ./ max(abs(data));

figure(1);
plot(data);

% Done applying Energy Measure

% Set constants
N = 200; % Starting frequency
windowSize = 70;
k = 1 + windowSize;
mu = 1; % change as needed
mu_alpha = 1;

% Generate wavetable (technically not a table)
x = linspace(0,2*pi,N);
figure(2);
plot(x);
w = cos(x);
dw = -sin(x);

% Apply LPF
b = (1/windowSize)*ones(1,windowSize);
a = 1;

% Set holding arrays
beta_1 = [zeros(1,windowSize), 2];
alpha = [zeros(1,windowSize), 2];

s = [zeros(1,windowSize), 1];
o = [zeros(1,windowSize), 1];

% Loop through whole song
for index = 1:num_chunks-k

    % Find max fit
    beta_possible = linspace(0,5,1000);
    dw_dbeta = dw(mod(k, N)+1) ./ (beta_possible - beta_1(k-1));
    beta_1(k+1) = beta_1(k) + mu * max(mean(filter(b,a,data(k-windowSize:k))) * dw_dbeta);
    
    % Adapt frequency
    alpha_possible = linspace(0,5,1000);
    dw_ds = dw(mod(k,N)+1) / (s(k) - s(k-1));
    ds_dalpha = (s(k) - s(k-1)) ./ (alpha_possible - alpha(k-1));
    
    alpha(k+1) = alpha(k) + mu_alpha .* max(mean(filter(b,a,data(k-windowSize:k))) * dw_ds * ds_dalpha);
    
    % Index into wavetable
    s(k+1) = mod(s(k-1) + alpha(k), N)+1;
    o(k+1) = w(floor(mod(s(k) + beta_1(k), N))+1);
    
    k = k + 1;

end

% Read the audio file
[data, Fs] = audioread('MapleLeafRag.ogg');

figure(3);
plot(1:num_chunks, o);

figure(4);
% Only preserve a click on the tracked beat
o_click = repelem(islocalmax(o), env_gap);
plot(1:num_chunks*env_gap, o_click);

beat_added = o_click*3 + data(1:num_chunks*env_gap);

beat_tracked_audio = audioplayer(beat_added, Fs);
play(beat_tracked_audio);

audiowrite('BeatTracked_Maple_Ilhan.wav', beat_added, Fs);