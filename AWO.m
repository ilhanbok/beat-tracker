% Attempt to track the beat of a song using an adaptive wavetable
% oscillator indexing into a cosine wavetable

% Ilhan Bok, 2020
% CC By Attribution

% Set constants
N = 100;
windowSize = 10; 
k = 1 + windowSize;
mu = 1; % change as needed

% Generate wavetable (technically not a table)
x = linspace(0,2*pi,N);
w = cos(x);
dw = -sin(x);

% Index into wavetable
s(k) = mod(s(k-1) + alpha, N);
o(k) = w(mod(s(k) + beta_1), N);

% Apply LPF
b = (1/windowSize)*ones(1,windowSize);
a = 1;

% Find max fit
beta_possible = linspace(0,5,1000);
dw_dbeta = (w(k) - w(k-1)) / (beta_1(k) - beta_1(k-1));
beta_1(k+1) = beta_1(k) + mu * max(filter(b,a,input * dw_dbeta));
