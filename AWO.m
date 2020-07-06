% Attempt to track the beat of a song using an adaptive wavetable
% oscillator indexing into a cosine wavetable

% Ilhan Bok, 2020
% CC By Attribution

format long;

% Read the audio file
[data, Fs] = audioread('pure_beat.wav');

% Set constants
N = 100;
windowSize = 10; 
k = 1 + windowSize;
mu = 1; % change as needed
mu_alpha = 1;

% Generate wavetable (technically not a table)
x = linspace(0,2*pi,N);
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
for index = 1:200-k

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

plot(1:200, o);