% Attempt to track the beat of a song using an adaptive wavetable
% oscillator indexing into a cosine wavetable

% Ilhan Bok, 2020
% CC By Attribution

format long;

% Read the audio file
[data, Fs] = audioread('MapleLeafRag.ogg');

% Play audio to hear sample
player = audioplayer(data(1:5000*250), Fs);
%play(player);

% Apply Energy Measure

% Start and end markers for energy derivation frame
env_start = 0;
env_end = 2000;

len = env_end - env_start;

% How much to shift the frames
env_gap = 250;

% How many sample frames to take
frames = 5000;

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

data = slopes;

% Done applying Energy Measure

% Set constants
N = 3000; % Starting frequency
windowSize = 50;
k = 1 + windowSize;
mu = 1; % change as needed
mu_alpha = 1;

% Generate wavetable (technically not a table)
x = linspace(0,2*pi,N);
figure(1);
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
for index = 1:5000-k

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

figure(2);
plot(1:5000, o);

figure(3);
% Only preserve a click on the tracked beat
o_click =  repelem(o > 0.99, env_gap);
plot(1:5000*250, o_click);

beat_tracked_audio = audioplayer(o_click + data(1:5000*250), Fs);
play(beat_tracked_audio);