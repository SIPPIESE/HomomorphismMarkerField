function [image, x, y, r] = FUNC_EMBEDDED_SAMPLE_V2(ref_g, size_field, size_acquisition, size_block, size_pitch, size_orient_pitch_count, sigma, detail_level)
k = size_acquisition / size_block / 2;
safe_margin = ceil(k * sqrt(2) / size_pitch);

x = rand() * (size_field - 2 * safe_margin - 1) + safe_margin + 1;
y = rand() * (size_field - 2 * safe_margin - 1) + safe_margin + 1;
r = rand() * 2 * pi;
% 
% x = 31.04;
% y = 62.72;
% r = 0;
% x = floor(floor(rand() * (size_field - 2 * safe_margin - 1) + safe_margin + 1) / 2) * 2;
% y = floor(floor(rand() * (size_field - 2 * safe_margin - 1) + safe_margin + 1) / 2) * 2;
% r = 0;

mesh_x = ((0:1:size_acquisition-1)-size_acquisition/2+0.5)';
mesh_y = (0:1:size_acquisition-1)-size_acquisition/2+0.5;
[mesh_y, mesh_x] = meshgrid(mesh_y, mesh_x);

phi_x = x * pi * size_pitch;
phi_y = y * pi * size_pitch;
f = k / size_acquisition;

%rotate
mesh_x_rot = -mesh_y*sin(-r)+mesh_x*cos(-r);
mesh_y_rot = mesh_y*cos(-r)+mesh_x*sin(-r);

ratio = 0.4;
square_mesh_x_rot = 2 * (2 * pi * f * mesh_x_rot + phi_x);
square_mesh_y_rot = 2 * (2 * pi * f * mesh_y_rot + phi_y);
square_x = ratio;
square_y = - ratio;
for i = 1:detail_level
    square_x = square_x + 2 / i / pi * sin(i * pi * ratio) * cos(i * square_mesh_x_rot);
    square_y = square_y - 2 / i / pi * sin(i * pi * ratio) * cos(i * square_mesh_y_rot);
end
square = square_y .* square_x;

square_mask = ones(size(square));
mesh_phase_x = (square_mesh_x_rot / 2 / pi + 0.5);
mesh_phase_y = (square_mesh_y_rot / 2 / pi + 0.5);
mesh_pitch_index_x = int64(floor(mesh_phase_x / size_pitch));
mesh_pitch_index_y = int64(floor(mesh_phase_y / size_pitch));
mesh_block_index_x = int64(mod(floor(mesh_phase_x), size_pitch));
mesh_block_index_y = int64(mod(floor(mesh_phase_y), size_pitch));
square_mask(mesh_block_index_x ~= 0) = 0;
square_mask(mesh_block_index_y ~= 0) = 0;
square_mask(ref_g((mesh_pitch_index_y - 1) .* size(ref_g, 1) + mesh_pitch_index_x) == 0) = 0;

if size_orient_pitch_count > 0
    mesh_orientation_block_index_x = int64(mod(floor(mesh_phase_x), size_orient_pitch_count * size_pitch));
    mesh_orientation_block_index_y = int64(mod(floor(mesh_phase_y), size_orient_pitch_count * size_pitch));
    square_mask(and(mesh_orientation_block_index_x == 1, mesh_orientation_block_index_y == 1)) = 1;
    square_mask(and(mesh_orientation_block_index_x == 2, mesh_orientation_block_index_y == 1)) = -1;
end

if mod(size_pitch, 2) ~= 0
    square_invert = mod(mesh_pitch_index_y + mesh_pitch_index_x, 2) == 1;
    square_mask(square_invert) = - square_mask(square_invert);
end

square = square .* square_mask;

checkerboard = 0.5;
mesh_checkerboard_x_rot = (2*pi*f*mesh_x_rot + phi_x + pi / 2);
mesh_checkerboard_y_rot = (2*pi*f*mesh_y_rot + phi_y + pi / 2);

for i = 1:detail_level
    for j = 1:detail_level
        checkerboard = checkerboard + 4/pi^2/(2*i-1)/(2*j-1)*(cos((2*i-1)*mesh_checkerboard_y_rot-(2*j-1)*mesh_checkerboard_x_rot)-cos((2*i-1)*mesh_checkerboard_y_rot+(2*j-1)*mesh_checkerboard_x_rot));
    end
end

image = checkerboard + square;
gray_min = min(image, [], 'all');
gray_max = max(image, [], 'all');
image = (image - gray_min) ./ (gray_max - gray_min);
end