T0 = 4;                                    % Period of the signal
dt = 0.001;                                % Sampling period 
t = -T0/2:dt:T0/2;                         % Time vector spanning one period [-T0/2, T0/2]
omega0 = 2*pi/T0;                          % Fundamental angular frequency

% Number of harmonics to compute 
N = 20;  

% Define original signal x(t) over one period
x = zeros(size(t));                        % Initialize x(t)
x(t >= -2 & t < -1) = 1;                   % x(t) = 1 for t in [-2, -1)
x(t >= -1 & t < 0) = t(t >= -1 & t < 0);   % x(t) = t for t in [-1, 0)
x(t >= 0 & t < 1) = t(t >= 0 & t < 1);     % x(t) = t for t in [0, 1)
x(t >= 1 & t <= 2) = 1;                    % x(t) = 1 for t in [1, 2]

cn_original = zeros(1, 2*N+1);  % Fourier coefficients (for n = -N to N)

for n = -N:N
    cn_original(n+N+1) = (1/T0) * trapz(t, x .* exp(-1i * n * omega0*t));
end

% Time-shifted signal (by T0/4 to the right)
shift = T0/4;               % Time shift amount
x_shifted = x;              % The signal remains the same, only time shifts

% Using the Fourier shift property
cn_shifted = cn_original .* exp(-1i * (-N:N) * omega0 * shift);  

% Magnitude and Phase of Fourier coefficients for original signal
magnitude_original = abs(cn_original);
phase_original = angle(cn_original);

% Magnitude and Phase of Fourier coefficients for shifted signal
magnitude_shifted = abs(cn_shifted);
phase_shifted = angle(cn_shifted);

% Frequency axis 
frequencies = (-N:N);

figure;

% Magnitude spectrum (original and shifted) 
subplot(2, 2, 1);
stem(frequencies, magnitude_original, 'b', 'LineWidth', 1.5);
title('Magnitude Spectrum of Original Signal', 'Interpreter', 'latex');
xlabel('$f_n$ (rad/s)', 'Interpreter', 'latex');
ylabel('$|C_n|$', 'Interpreter', 'latex');
ylim("padded")
grid on;

subplot(2, 2, 2);
stem(frequencies, magnitude_shifted, 'r', 'LineWidth', 1.5);
title('Magnitude Spectrum of Shifted Signal', 'Interpreter', 'latex');
xlabel('$f_n$ (rad/s)', 'interpreter', 'latex');
ylabel('$|C_n|$', 'interpreter', 'latex');
ylim("padded")
grid on;

% Phase spectrum (original and shifted) 
subplot(2, 2, 3);
stem(frequencies, phase_original, 'b', 'LineWidth', 1.5);
title('Phase Spectrum of Original Signal', 'Interpreter', 'latex');
xlabel('$f_n$ (rad/s)', 'interpreter', 'latex');
ylabel('$\arg(C_n)$', 'Interpreter', 'latex');
grid on;

subplot(2, 2, 4);
stem(frequencies, phase_shifted, 'r', 'LineWidth', 1.5);
title('Phase Spectrum of Shifted Signal', 'Interpreter', 'latex');
xlabel('$f_n$ (rad/s)', 'Interpreter', 'latex');
ylabel('$\arg(C_n)$', 'Interpreter', 'latex');
grid on;

sgtitle('Magnitude and Phase Spectra', 'Interpreter', 'latex');
