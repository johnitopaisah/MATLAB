% Image Registration using Cross-Correlation
%
%   STUDENT_1 -> John Itopa ISAH - 63362
%   STUDENT_2 -> Blesso Danny Jeba Stephen - 63291
%
% Description:
% This MATLAB script performs image registration using the cross-correlation
% technique. It reads an input image from a BMP file, calculates the
% cross-correlation between adjacent rows of the image to determine the
% shift required for alignment, and applies the calculated shifts to
% register the image. The script displays both the original and registered
% images for comparison.
%
% Usage:
% 1. Ensure the input image file is saved in BMP format and placed in the
%    same directory as this script.
% 2. Modify the filename in the script to match the name of your input image
%    file, if different from the default ('fichier2.bmp').
% 3. Run the script in MATLAB.
% 4. View the original and registered images displayed in separate figures
%    to observe the effect of image registration.
%
% Implementation Details:
% - The script reads the input image from a BMP file and converts it to a
%   grayscale image.
% - Cross-correlation is performed between adjacent rows of the image to
%   calculate the shift required for alignment.
% - Cumulative shift values are computed and applied to register the image.
%
% Files:
% - `ImageRegistration.m`: MATLAB script for image registration.
% - Sample image file ('fichier2.bmp') for demonstration.
%
% Dependencies:
% - MATLAB Image Processing Toolbox (for image read/write and manipulation).
%
% Note:
% - Ensure that the input image file is present in the same directory as
%   the script and is in BMP format.


% Clear Workspace and Close Figures
clc;
clearvars;
close all;

% Read Input Image
input_image = imread('fichier2.bmp', 'bmp');

% Convert Image to Grayscale
input_image_gray = 255 * input_image;

% Display Original Image
figure('Name', 'Original Image');
image(input_image_gray);
colormap(gray);
title('Original Image');

% Initialize Variables
[num_rows, num_cols] = size(input_image_gray);
shifts = zeros(num_rows, 1);

% Calculate Shifts using Cross-Correlation
for row = 1:num_rows - 1
    cross_corr = xcorr(input_image_gray(row + 1, :), input_image_gray(row, :));
    [~, max_index] = max(cross_corr);
    shifts(row + 1) = num_cols - max_index;
end

% Compute Cumulative Shifts
cumulative_shifts = cumsum(shifts);

% Display Cumulative Shifts
disp(cumulative_shifts);

% Register Image using Shifts
registered_image = zeros(num_rows, num_cols);
registered_image(1, :) = input_image_gray(1, :);
for k = 2:num_rows
    registered_image(k, :) = circshift(input_image_gray(k, :), cumulative_shifts(k));
end

% Display Registered Image
figure('Name', 'Registered Image');
image(registered_image);
colormap(gray);
title('Registered Image');
