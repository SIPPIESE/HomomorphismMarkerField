%% parameters for comparing decoding time
clear;
clc;
close all;

parameter = [
    4, 4, 32, 640;
    4, 4, 40, 640;
    4, 4, 48, 640;
    4, 4, 56, 640;
    4, 4, 64, 640;
    4, 4, 72, 640;
    4, 4, 80, 640;
    4, 4, 88, 640;
    4, 4, 96, 640;
    4, 4, 104, 640;
    4, 4, 112, 640;
    4, 4, 120, 640;
    4, 4, 128, 640;
    ];

data_folder = "data/dbt/";
size_block = 10;
size_pitch = 5;
sigma = 0;
test_count = 100;
config = {data_folder, size_block, size_pitch, sigma, test_count};

%% parameters for comparing precision impact
clear;
clc;
close all;

parameter = [
    4, 4, 72, 640;
    4, 4, 72, 720;
    4, 4, 72, 800;
    4, 4, 72, 880;
    4, 4, 72, 960;
    4, 4, 72, 1040;
    4, 4, 72, 1120;
    4, 4, 72, 1200;
    4, 4, 72, 1280;
    4, 4, 72, 1360;
    4, 4, 72, 1440;
    ];

data_folder = "data/dbt/";
size_block = 10;
size_pitch = 5;
sigma = 0;
test_count = 100;
config = {data_folder, size_block, size_pitch, sigma, test_count};

%% simulations
time = zeros(size(parameter, 1), 1);
rmse = zeros(size(parameter, 1), 3);
for i = 1:size(parameter, 1)
    [rmse(i, :), time(i)] = SIMULATE_DBT(parameter(i, :), config);
end
disp(time);
disp(rmse);

%% do table
times = ' \times ';
for i = 1:size(parameter, 1)
    a_s = "$" + string(parameter(i, 4)) + times + string(parameter(i, 4)) + "$";
    w_s = "$" + string(parameter(i, 1)) + times + string(parameter(i, 2)) + "$";
    p_s = "$" + string(parameter(i, 3)) + times + string(parameter(i, 3)) + "$";
    rmse_x_s = "$" + strrep(num2str(rmse(i, 1) .* 40, '%.3e'), 'e-0', '\times 10^{-') + '}' + "$";
    rmse_y_s = "$" + strrep(num2str(rmse(i, 2) .* 40, '%.3e'), 'e-0', '\times 10^{-') + '}' + "$";
    rmse_r_s = "$" + strrep(num2str(rmse(i, 3), '%.3e'), 'e-0', '\times 10^{-') + '}' + "$";
    t_s = "$" + num2str(time(i) * 1000, '%.4f') + "$";
    row_s = p_s + " & " + a_s + " & " + w_s + " & " + rmse_x_s + " & " + rmse_y_s + " & " + rmse_r_s + " & " + t_s + "\\";
    disp(row_s)
end


function [rmse, time] = SIMULATE_DBT(parameter, config)
size_window = parameter(1:2);
size_field = parameter(3);
size_acquisition = parameter(4);
data_folder = config{1};
size_block = config{2};
size_pitch = config{3};
sigma = config{4};
test_count = config{5};

ref_g = readmatrix(data_folder + string(size_window(1)) + 'x' + string(size_window(2)) + '_G.txt');
ref_g = ref_g(1:size_field, 1:size_field);
size_g = size(ref_g);

real = zeros(test_count, 3);
esti = zeros(test_count, 3);
error = zeros(test_count, 3);
time = zeros(test_count, 1);

for i = 1:test_count
    [image, x, y, r] = FUNC_EMBEDDED_SAMPLE_V2(ref_g, size_field, size_acquisition, size_block, size_pitch, size_window(1), sigma, 1);
    [vec_m, vec_n, bias_m, bias_n, theta] = FUNC_EMBEDDED_ESTIMATE(image);
    try
        [index_m, index_n, seq_p, indent_m, indent_n, flag_m_is_x, flip_m, flip_n] = EXTRACT_DBT(image, size(image), size_window, size_pitch, bias_m, bias_n, vec_m, vec_n);
        tic;
        [est_x, est_y, est_r] = LOCATE_DBT(ref_g, size_g, size_window, size_pitch, bias_m, bias_n, theta, index_m, index_n, seq_p, indent_m, indent_n, flag_m_is_x, flip_m, flip_n);
        time(i) = toc;
        if (est_r - r > pi)
            est_r = est_r - 2 * pi;
        end
        if (est_r - r < -pi)
            est_r = est_r + 2 * pi;
        end
        real(i, :) = [x * size_pitch / 2, y * size_pitch / 2, r];
        esti(i, :) = [est_x, est_y, est_r];
        error(i, :) = esti(i, :) - real(i, :);
        if mod(i, 2) == 1
            temp = error(i, 1);
            error(i, 1) = error(i, 2);
            error(i, 2) = temp;
        end
        if (error(i, 3) > 0.1)
            disp("Angle error")
            pause
        end
        fprintf('Test %i completes with error: %d, %d, %d\n', i, error(i, 1), error(i, 2), error(i, 3));
    catch
        disp('Error occured during simulation')
        pause;
    end
end

rmse = sqrt(1/test_count*sum(error.^2, 1));
time = sum(time) / test_count;
disp('Test complete, RMSE:');
disp(rmse);
disp('Average extraction + decoding time:');
disp(time);
end

function [index_m, index_n, ext_p, indent_m, indent_n, flag_m_is_x, flip_m, flip_n] = EXTRACT_DBT(image, size, window, pitch, bias_m, bias_n, vec_m, vec_n)
% gray_mean = mean(image, 'all');
gray_mean = 0.5;
coord_origin = [(size(1) + 1) / 2, (size(2) + 1) / 2] - bias_m .* vec_m - bias_n .* vec_n;
[index_m, index_n] = FIND_CODE_DBT(image, size, window, pitch, vec_m, vec_n, coord_origin, gray_mean);

ext_p = -ones(window(1) * 2 - 1, window(2) * 2 - 1);
indent_m = window(1);
indent_n = window(2);
flag_find_orientation_code = false;
flag_m_is_x = false;
flip_m = false;
flip_n = false;
% debug = image;

for i = -(window(1)-1):window(1)-1
    for j = -(window(2)-1):window(2)-1
        [bound, value] = GET_CODE_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m + pitch * i, index_n + pitch * j);
        if (~flag_find_orientation_code && bound)
            [bound, corner] = GET_CORNER_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m + pitch * i, index_n + pitch * j);
            if (bound && norm(corner) ~= 0)
                if (corner(1) == 1)
                    flip_m = false;
                else
                    flip_m = true;
                end
                if (corner(2) == 1)
                    flip_n = false;
                else
                    flip_n = true;
                end
                [bound, nearby] = GET_NEARBY_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m + pitch * i + corner(1), index_n + pitch * j + corner(2));
                if (bound && norm(nearby) ~= 0)
                    flag_find_orientation_code = true;
                end
                if (abs(nearby(1)) > 0)
                    flag_m_is_x = true;
                else
                    flag_m_is_x = false;
                end
            end
        end
        if (bound)
            if (value)
                ext_p(i + window(1), j + window(2)) = 1;
            else
                ext_p(i + window(1), j + window(2)) = 0;
            end
        end
    end
end
if (~flag_find_orientation_code)
    disp("No orientation code found!")
end
for i = 1:window(1)
    for j = 1:window(2)
        flag = true;
        for x = 1:window(1)
            for y = 1:window(2)
                if (ext_p(i + x - 1, j + y - 1) == -1)
                    flag = false;
                end
            end
        end
        if (flag)
            ext_p = ext_p(i:i + window(1) - 1, j:j + window(2) - 1);
            indent_m = indent_m - (i - 1);
            indent_n = indent_n - (j - 1);
            return
        end
    end
end
end

function [index_m, index_n] = FIND_CODE_DBT(image, size, window, pitch, vec_m, vec_n, coord_origin, gray_mean)
% mathematically, for every n + 1 window there is a coded block.
index_m = 0;
index_n = 0;
for i = 0:pitch * window(1) % n + 1 points
    for j = 0:pitch * window(2)  % n + 1 points
        [bound_a, value] = GET_CODE_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, i, j);
        [bound_b, nearby] = GET_NEARBY_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, i, j);
        if (bound_a && bound_b && value && norm(nearby) == 0)
            index_m = mod(i, pitch);
            index_n = mod(j, pitch);
            return
        end
        [bound_a, value] = GET_CODE_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, i, -j);
        [bound_b, nearby] = GET_NEARBY_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, i, -j);
        if (bound_a && bound_b && value && norm(nearby) == 0)
            index_m = mod(i, pitch);
            index_n = mod(-j, pitch);
            return
        end
        [bound_a, value] = GET_CODE_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, -i, -j);
        [bound_b, nearby] = GET_NEARBY_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, -i, -j);
        if (bound_a && bound_b && value && norm(nearby) == 0)
            index_m = mod(-i, pitch);
            index_n = mod(-j, pitch);
            return
        end
        [bound_a, value] = GET_CODE_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, -i, j);
        [bound_b, nearby] = GET_NEARBY_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, -i, j);
        if (bound_a && bound_b && value && norm(nearby) == 0)
            index_m = mod(-i, pitch);
            index_n = mod(j, pitch);
            return
        end
    end
end
disp("Search failure, please increase sample area to include window + 1 blocks")
end

function [bound, corner] = GET_CORNER_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m, index_n)
corner = [0, 0];
bound = true;
[bound_a, value] = GET_CODE_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m + 1, index_n + 1);
if (bound_a && value)
    corner = [1, 1];
    return
end
[bound_b, value] = GET_CODE_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m - 1, index_n + 1);
if (bound_b && value)
    corner = [-1, 1];
    return
end
[bound_c, value] = GET_CODE_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m - 1, index_n - 1);
if (bound_c && value)
    corner = [-1, -1];
    return
end
[bound_d, value] = GET_CODE_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m + 1, index_n - 1);
if (bound_d && value)
    corner = [1, -1];
    return
end
bound = bound_a && bound_b && bound_c && bound_d;
end

function [bound, nearby] = GET_NEARBY_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m, index_n)
nearby = [0, 0];
bound = true;
[bound_a, value] = GET_CODE_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m + 1, index_n);
if (bound_a && value)
    nearby = [1, 0];
    return
end
[bound_b, value] = GET_CODE_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m, index_n + 1);
if (bound_b && value)
    nearby = [0, 1];
    return
end
[bound_c, value] = GET_CODE_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m - 1, index_n);
if (bound_c && value)
    nearby = [-1, 0];
    return
end
[bound_d, value] = GET_CODE_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m, index_n - 1);
if (bound_d && value)
    nearby = [0, -1];
    return
end
bound = bound_a && bound_b && bound_c && bound_d;
end

function [bound, value] = GET_CODE_DBT(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m, index_n)
bound = false;
value = false;
coord = round(coord_origin + vec_m ./2 .* index_m + vec_n ./2 .* index_n);
if (coord(1) >= 1 && coord(1) <= size(1) && coord(2) >= 1 && coord(2) <= size(2))
    gray = image(coord(1), coord(2));
    if ((mod(index_m + index_n, 2) == 0 && gray < gray_mean) || (mod(index_m + index_n, 2) ~= 0 && gray > gray_mean))
        bound = true;
        value = true;
        return
    end
    bound = true;
    value = false;
end
end

% Search function
function [index_x, index_y]= DECODE_DBT(ref_g, size_g, window, seq_p)
index_x = 0;
index_y = 0;
for i = 1:size_g - window(1) + 1
    for j = 1:size_g - window(2) + 1
        ref_p = ref_g(i:i + window(1) - 1, j:j + window(2) - 1);
        if (all(ref_p == seq_p))
            index_x = i;
            index_y = j;
        end
    end
end
end

% Decode sequence
function [period_x, period_y, est_r] = LOCATE_DBT(ref_g, size_g, window, pitch, bias_m, bias_n, theta, index_m, index_n, seq_p, indent_m, indent_n, flag_m_is_x, flip_m, flip_n)
if (flip_m)
    seq_p = flip(seq_p, 1);
end
if (flip_n)
    seq_p = flip(seq_p, 2);
end
flip_x = flip_m;
flip_y = flip_n;
if (~flag_m_is_x)
    seq_p = seq_p';
    flip_x = flip_n;
    flip_y = flip_m;
end
[index_x, index_y] = DECODE_DBT(ref_g, size_g, window, seq_p);
if (flag_m_is_x)
    if (flip_x)
        period_x = ((index_x + (window(1) - indent_m)) * pitch + index_m) / 2;
        period_x = period_x - bias_m;
    else
        period_x = ((index_x + (indent_m - 1)) * pitch - index_m) / 2;
        period_x = period_x + bias_m;
    end
    if (flip_y)
        period_y = ((index_y + (window(2) - indent_n)) * pitch + index_n) / 2;
        period_y = period_y - bias_n;
    else
        period_y = ((index_y + (indent_n - 1)) * pitch - index_n) / 2;
        period_y = period_y + bias_n;
    end
else
    if (flip_x)
        period_x = ((index_x + (window(2) - indent_n)) * pitch + index_n) / 2;
        period_x = period_x - bias_n;
    else
        period_x = ((index_x + (indent_n - 1)) * pitch - index_n) / 2;
        period_x = period_x + bias_n;
    end
    if (flip_y)
        period_y = ((index_y + (window(1) - indent_m)) * pitch + index_m) / 2;
        period_y = period_y - bias_m;
    else
        period_y = ((index_y + (indent_m - 1)) * pitch - index_m) / 2;
        period_y = period_y + bias_m;
    end
end
% derive orientation from the direction of axes
est_r = 0 + theta * 180 / pi;
if isequal([0, 1], [flip_x, flip_y])
    est_r = 90 + theta * 180 / pi;
elseif isequal([1, 0], [flip_x, flip_y])
    est_r = 270 + theta * 180 / pi;
elseif isequal([1, 1], [flip_x, flip_y])
    est_r = 180 + theta * 180 / pi;
end
est_r = est_r / 360 * (2 * pi);
end