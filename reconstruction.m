T0 = 4;                                   % Period of the original signal
N = 1000;                                 % Number of harmonics for reconstruction
fs = 2*N / T0;                            % Sampling frequency
dt = 1 / fs;                              % Sampling period (1/fs)
t = -T0/2:dt:T0/2;                        % Time vector spanning one period [-T0/2, T0/2]
omega0 = 2*pi/T0;                         % Fundamental angular frequency

% Define original signal x(t) over one period
x = zeros(size(t));                       % Initialize x(t)
x(t >= -2 & t < -1) = 1;                  % x(t) = 1 for t in [-2, -1)
x(t >= -1 & t < 0) = t(t >= -1 & t < 0);  % x(t) = t for t in [-1, 0)
x(t >= 0 & t < 1) = t(t >= 0 & t < 1);    % x(t) = t for t in [0, 1)
x(t >= 1 & t <= 2) = 1;                   % x(t) = 1 for t in [1, 2]

x_reconstructed = zeros(size(t));

for n = -N:N
    c_n = (1 / T0) * trapz(t, x .* exp(-1i * n * omega0 * t));
    x_reconstructed = x_reconstructed + c_n * exp(1i * n * omega0 * t);
end

x_reconstructed = real(x_reconstructed);

figure;
subplot(2, 1, 1);
plot(t, x, 'LineWidth', 1.5);
title('Original Signal x(t)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
xlim([-2, 2]);
ylim([-1.5, 1.5]);

subplot(2, 1, 2);
plot(t, x_reconstructed, 'r', 'LineWidth', 1.5);
title(['Reconstructed Signal with ', num2str(N), ' Harmonics']);
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
xlim([-2, 2]);
ylim([-1.5, 1.5]);

harmonics = [1, 5, 10, 20, 100, 1000];  
figure;
for i = 1:length(harmonics)
    N_current = harmonics(i);
    x_reconstructed = zeros(size(t));
    for n = -N_current:N_current
        c_n = (1 / T0) * trapz(t, x .* exp(-1i * n * omega0 * t));
        x_reconstructed = x_reconstructed + c_n * exp(1i * n * omega0 * t);
    end

    x_reconstructed = real(x_reconstructed);
    
    subplot(3, 2, i);
    plot(t, x, 'b', t, x_reconstructed, 'r', 'LineWidth', 1.2);
    title(['N = ', num2str(N_current)]);
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Original', 'Reconstructed');
    grid on;
    xlim([-2, 2]);
    ylim([-1.5, 1.5]);
end
