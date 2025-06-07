clear;
clc;
close all;

parameter = [
    10, 72, 640;
    ];

data_folder = "data/hmf/";
img_folder = "img/hmf/";
size_block = 10;
size_pitch = 5;
sigma = 0;
flag_delta_xy = false;
flag_rotation = false;
test_count = 1000;
config = {data_folder, img_folder, size_block, size_pitch, sigma, flag_delta_xy, flag_rotation, test_count};
ref_g = readmatrix(data_folder + string(10) + '+' + string(10) + '_G.txt');
FUNC_EMBEDDED_SAMPLE_V2(parameter, config, ref_g);


function [image, x, y, r] = FUNC_EMBEDDED_SAMPLE_V2(parameter, config, ref_g)
size_window = parameter(1);
size_field = parameter(2);
size_acquisition = parameter(3);
data_folder = config{1};
img_folder = config{2};
size_block = config{3};
size_pitch = config{4};
sigma = config{5};
ref_x = ref_g(:, 1);
ref_y = ref_g(1, :);

k = size_acquisition / size_block / 2;
safe_margin = ceil(ceil(k / 2 * sqrt(2)) / size_pitch) + 1;

x = rand() * (size_field - 2 * safe_margin - 1) + safe_margin + 1;
y = rand() * (size_field - 2 * safe_margin - 1) + safe_margin + 1;

x = 0;
y = 0;

r = rand() * 2 * pi;

mesh_x = (0:1:size_acquisition-1)-size_acquisition/2+0.5;
mesh_y = ((size_acquisition-1:-1:0)-size_acquisition/2+0.5)';
[mesh_x, mesh_y] = meshgrid(mesh_x, mesh_y);
phi_x = x * pi * size_pitch;
phi_y = y * pi * size_pitch;
detail_level = 10;
f = k / size_acquisition;

%rotate
mesh_x_rot = mesh_x*cos(r)+mesh_y*sin(r);
mesh_y_rot = -mesh_x*sin(r)+mesh_y*cos(r);

ratio = 0.4;
square_mesh_x_rot = 2 * (2 * pi * f * mesh_x_rot + phi_x);
square_mesh_y_rot = 2 * (2 * pi * f * mesh_y_rot + phi_y);
square_x = - ratio;
square_y = ratio;
for i = 1:detail_level
    square_x = square_x - 2 / i / pi * sin(i * pi * ratio) * cos(i * square_mesh_x_rot);
    square_y = square_y + 2 / i / pi * sin(i * pi * ratio) * cos(i * square_mesh_y_rot);
end
square = square_x .* square_y;
mesh_pitch_index_x = int64(floor((square_mesh_x_rot / 2 / pi + 0.5) / size_pitch));
mesh_pitch_index_y = int64(floor((square_mesh_y_rot / 2 / pi + 0.5) / size_pitch));
mesh_block_index_x = int64(mod(floor(square_mesh_x_rot / 2 / pi + 0.5), size_pitch));
mesh_block_index_y = int64(mod(floor(square_mesh_y_rot / 2 / pi + 0.5), size_pitch));
square(mesh_block_index_x ~= 0) = 0;
square(mesh_block_index_y ~= 0) = 0;
square(xor(ref_x(mesh_pitch_index_x), ref_y(mesh_pitch_index_y)) == 0) = 0;
if mod(size_pitch, 2) ~= 0
    square(mod(mesh_pitch_index_x + mesh_pitch_index_y, 2) == 1) = - square(mod(mesh_pitch_index_x + mesh_pitch_index_y, 2) == 1);
end

checkerboard = 0.5;
mesh_checkerboard_x_rot = (2*pi*f*mesh_x_rot + phi_x + pi / 2);
mesh_checkerboard_y_rot = (2*pi*f*mesh_y_rot + phi_y + pi / 2);
for i = 1:detail_level
    for j = 1:detail_level
        checkerboard = checkerboard + 4/pi^2/(2*i-1)/(2*j-1)*(cos((2*i-1)*mesh_checkerboard_x_rot-(2*j-1)*mesh_checkerboard_y_rot)-cos((2*i-1)*mesh_checkerboard_x_rot+(2*j-1)*mesh_checkerboard_y_rot));
    end
end

image = checkerboard + square;
gray_min = min(image, [], 'all');
gray_max = max(image, [], 'all');
image = (image - gray_min) ./ (gray_max - gray_min);
imshow(image)
end