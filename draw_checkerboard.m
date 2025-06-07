clear;
clc;

% This script file embeds the coding into a common grid pattern with noise
% and blur, which is used for measuring phase estimation precision loss.
% Here, select the available window length, recommended 4 - 7

hmf_data_folder = "data/hmf/";
umf_data_folder = "data/umf/";
dbt_data_folder = "data/dbt/";
hmf_img_folder = "img/hmf/";
umf_img_folder = "img/umf/";
dbt_img_folder = "img/dbt/";
empty_img_folder = "img/empty/";

field_size_range = 16:8:256;
pitch = 4;
ratio = 2/5;
block_size = 10;

pitch = pitch + 1;
for field_size = field_size_range
    DRAW_hmf(field_size, pitch, ratio, block_size, hmf_data_folder, hmf_img_folder);
    DRAW_UMF(field_size, pitch, ratio, block_size, umf_data_folder, umf_img_folder);
    DRAW_DBT(field_size, pitch, ratio, block_size, dbt_data_folder, dbt_img_folder)
    DRAW_EMPTY(field_size, pitch, ratio, block_size, empty_img_folder);
end

function [img] = EMBED(G, L, pitch, ratio, block, rotation)
img = zeros(L, L, 'uint8');
x = (block + 1) / 2;
y = (block + 1) / 2;
padding = (block - 1) / 2;
code = floor(block * ratio);
if (mod(block, 2) ~= mod(code, 2))
    code = code - 1;
end
code_padding = (code - 1) / 2;
x_coord = 1;
y_coord = 1;
while (x < L)
    while (y < L)
        if (mod(x_coord + y_coord, 2) == 0)
            img(x - padding:x + padding, y - padding:y + padding) = 255;
        end
        if (mod(x_coord, pitch) == 0 && mod(y_coord, pitch) == 0)
            x_index = x_coord / pitch;
            y_index = y_coord / pitch;
            if (G(x_index, y_index) == 1)
                if (mod(x_coord + y_coord, 2) == 0)
                    img(x - code_padding:x + code_padding, y - code_padding:y + code_padding) = 0;
                else
                    img(x - code_padding:x + code_padding, y - code_padding:y + code_padding) = 255;
                end
            end
        end
        if (length(rotation) == 2)
            if (mod(x_coord, rotation(1) * pitch) == 1 && mod(y_coord, rotation(2) * pitch) == 1)
                if (mod(x_coord + y_coord, 2) == 0)
                    img(x - code_padding:x + code_padding, y - code_padding:y + code_padding) = 0;
                else
                    img(x - code_padding:x + code_padding, y - code_padding:y + code_padding) = 255;
                end
            end
            if (mod(x_coord, rotation(1) * pitch) == 2 && mod(y_coord, rotation(2) * pitch) == 1)
                if (mod(x_coord + y_coord, 2) == 0)
                    img(x - code_padding:x + code_padding, y - code_padding:y + code_padding) = 0;
                else
                    img(x - code_padding:x + code_padding, y - code_padding:y + code_padding) = 255;
                end
            end
        end
        y = y + block;
        y_coord = y_coord + 1;
    end
    x = x + block;
    x_coord = x_coord + 1;
    y = (block + 1) / 2;
    y_coord = 1;
end
img = im2uint8(img);
end

function [field_size] = DRAW_hmf(field_size, pitch, ratio, block_size, hmf_data_folder, hmf_img_folder)
field_size_list = [4, 6, 10, 16, 28, 39, 91, 210, 385, 838, 1637, 3600, 7069];
window_size_list = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16];
window_size = -1;
file_name = '';
flag_find = false;
for j = 1:size(window_size_list, 2)
    if (field_size_list(j) >= field_size)
        window_size = window_size_list(j);
        file_name = string(window_size) + '+' + string(window_size) + '_' + string(field_size_list(j)) + '_G.txt';
        flag_find = true;
        break;
    end
end
if (~flag_find)
    return
end
G = readmatrix(hmf_data_folder + file_name);
G = G(1:field_size, 1:field_size);
L = field_size * pitch * block_size + (pitch - 1) * block_size;
img = EMBED(G, L, pitch, ratio, block_size, []);
disp("hmf >> Writing " + string(window_size) + '+' + string(window_size) + '_' + string(field_size));
imwrite(img, hmf_img_folder + string(window_size) + '+' + string(window_size) + '_' + string(field_size) + '.png');
disp("hmf >> Complete");
end

function [] = DRAW_UMF(field_size, pitch, ratio, block_size, umf_data_folder, umf_img_folder)
field_size_list = [77];
window_size_list = [4];
window_size = -1;
file_name = '';
flag_find = false;
for j = 1:size(window_size_list, 2)
    if (field_size_list(j) >= field_size)
        window_size = window_size_list(j);
        file_name = string(window_size) + 'x' + string(window_size) + '_' + string(field_size_list(j)) + '_G.txt';
        flag_find = true;
        break;
    end
end
if (~flag_find)
    return
end
G = readmatrix(umf_data_folder + file_name);
G = G(1:field_size, 1:field_size);
L = field_size * pitch * block_size + (pitch - 1) * block_size;
img = EMBED(G, L, pitch, ratio, block_size, []);
disp("UMF >> Writing " + string(window_size) + 'x' + string(window_size) + '_' + string(field_size));
imwrite(img, umf_img_folder + string(window_size) + 'x' + string(window_size) + '_' + string(field_size) + '.png');
disp("UMF >> Complete");
end

function [] = DRAW_DBT(field_size, pitch, ratio, block_size, dbt_data_folder, dbt_img_folder)
field_size_list = [16, 256, 4096, 8192];
window_size_list = [3, 3; 4, 4; 5, 5; 5, 6];
window_size = [];
file_name = '';
flag_find =  false;
for j = 1:size(window_size_list, 1)
    if (field_size_list(j) >= field_size)
        window_size = window_size_list(j, :);
        file_name = string(window_size(1)) + 'x' + string(window_size(2)) + '_' + string(field_size_list(j)) + '_G.txt';
        flag_find = true;
        break;
    end
end
if (~flag_find)
    return
end
G = readmatrix(dbt_data_folder + file_name);
G = G(1:field_size, 1:field_size);
L = field_size * pitch * block_size + (pitch - 1) * block_size;
img = EMBED(G, L, pitch, ratio, block_size, window_size);
disp("DBT >> Writing " + string(window_size(1)) + "x" + string(window_size(2)) + '_' + string(field_size));
imwrite(img, dbt_img_folder + string(window_size(1)) + "x" + string(window_size(2)) + '_' + string(field_size) + '.png');
disp("DBT >> Complete");
end

function [] = DRAW_EMPTY(field_size, pitch, ratio, block_size, empty_img_folder)
G = zeros(field_size, field_size, 'uint8');
L = field_size * pitch * block_size + (pitch - 1) * block_size;
img = EMBED(G, L, pitch, ratio, block_size, []);
disp("EMPTY >> Writing " + string(field_size));
imwrite(img, empty_img_folder + string(field_size) + '.png');
disp("EMPTY >> Complete");
end