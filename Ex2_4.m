% Define the input signal
n = 0:100;
x = sin(0.1*pi*n); % Example input signal (sine wave)

% Compute the output of the filter
y = zeros(size(x));
y(1) = x(1);
for i = 2:length(x)
    y(i) = x(i) + x(i-1);
end

% Plot the input signal
subplot(3, 1, 1);
stem(n, x);
xlabel('n');
ylabel('x(n)');
title('Input Signal');

% Plot the output of the filter
subplot(3, 1, 2);
stem(n, y);
xlabel('n');
ylabel('y(n)');
title('Output Signal');

% Plot the block diagram
subplot(3, 1, 3);
blockDiagram = imread('block_diagram.png'); % Replace 'block_diagram.png' with the actual filename of your block diagram image
imshow(blockDiagram);
title('Implementation Scheme');
