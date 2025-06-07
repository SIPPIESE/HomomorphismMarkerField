%% parameters for comparing decoding time
clear;
clc;
close all;

parameter = [
    9, 32, 640;
    10, 40, 640;
    10, 48, 640;
    10, 56, 640;
    10, 64, 640;
    10, 72, 640;
    10, 80, 640;
    10, 88, 640;
    11, 96, 640;
    11, 104, 640;
    11, 112, 640;
    11, 120, 640;
    11, 128, 640;
    ];

data_folder = "data/hmf/";
size_block = 10;
size_pitch = 5;
sigma = 0;
test_count = 1000;
config = {data_folder, size_block, size_pitch, sigma, test_count};

%% parameters for comparing precision impact
clear;
clc;
close all;

parameter = [
    10, 72, 640;
    10, 72, 720;
    10, 72, 800;
    10, 72, 880;
    10, 72, 960;
    10, 72, 1040;
    10, 72, 1120;
    10, 72, 1200;
    10, 72, 1280;
    10, 72, 1360;
    10, 72, 1440;
    ];

data_folder = "data/hmf/";
size_block = 10;
size_pitch = 5;
sigma = 0;
test_count = 1000;
config = {data_folder, size_block, size_pitch, sigma, test_count};

%% simulations
time = zeros(size(parameter, 1), 1);
rmse = zeros(size(parameter, 1), 3);
for i = 1:size(parameter, 1)
    [rmse(i, :), time(i)] = SIMULATE_HMF(parameter(i, :), config);
end
disp(time);
disp(rmse);

%% do table
times = ' \times ';
for i = 1:size(parameter, 1)
    a_s = "$" + string(parameter(i, 3)) + times + string(parameter(i, 3)) + "$";
    w_s = "$" + string(parameter(i, 1)) + "+" + string(parameter(i, 1)) + "$";
    p_s = "$" + string(parameter(i, 2)) + times + string(parameter(i, 2)) + "$";
    rmse_x_s = "$" + strrep(num2str(rmse(i, 1) .* 40, '%.3e'), 'e-0', '\times 10^{-') + '}' + "$";
    rmse_y_s = "$" + strrep(num2str(rmse(i, 2) .* 40, '%.3e'), 'e-0', '\times 10^{-') + '}' + "$";
    rmse_r_s = "$" + strrep(num2str(rmse(i, 3), '%.3e'), 'e-0', '\times 10^{-') + '}' + "$";
    t_s = "$" + num2str(time(i) * 1000, '%.4f') + "$";
    row_s = p_s + " & " + a_s + " & " + w_s + " & " + rmse_x_s + " & " + rmse_y_s + " & " + rmse_r_s + " & " + t_s + "\\";
    disp(row_s)
end

function [rmse, time] = SIMULATE_HMF(parameter, config)
size_window = parameter(1);
size_field = parameter(2);
size_acquisition = parameter(3);
data_folder = config{1};
size_block = config{2};
size_pitch = config{3};
sigma = config{4};
test_count = config{5};

ref_x = readmatrix(data_folder + string(size_window) + '+' + string(size_window) + '_X.txt');
ref_y = readmatrix(data_folder + string(size_window) + '+' + string(size_window) + '_Y.txt');
ref_g = readmatrix(data_folder + string(size_window) + '+' + string(size_window) + '_G.txt');
len_x = length(ref_x);
len_y = length(ref_y);

real = zeros(test_count, 3);
esti = zeros(test_count, 3);
error = zeros(test_count, 3);
time = zeros(test_count, 1);

for i = 1:test_count
    [image, x, y, r] = FUNC_EMBEDDED_SAMPLE_V2(ref_g, size_field, size_acquisition, size_block, size_pitch, 0, sigma, 7);
    
    image = imfilter(image, fspecial('disk', 2));
    image = image + randn(size(image)) * 0.025;
    imshow(image)
    imwrite(image, "sample_l.png")
    [vec_m, vec_n, bias_m, bias_n, theta] = FUNC_EMBEDDED_ESTIMATE(image);
    try
        [seq_m, seq_n, index_m, index_n, indent_m, indent_n] = EXTRACT_HMF(image, size(image), size_window, size_pitch, bias_m, bias_n, vec_m, vec_n);
        tic;
        [est_x, est_y, est_r] = LOCATE_HMF(ref_x, ref_y, len_x, len_y, size_window, size_pitch, bias_m, bias_n, seq_m, seq_n, index_m, index_n, indent_m, indent_n, theta);
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

function [seq_m, seq_n, index_m, index_n, indent_m, indent_n] = EXTRACT_HMF(image, size, window, pitch, bias_m, bias_n, vec_m, vec_n)
% gray_mean = mean(image, 'all');
gray_mean = 0.5;
coord_origin = [(size(1) + 1) / 2, (size(2) + 1) / 2] - bias_m .* vec_m - bias_n .* vec_n;
[index_m, index_n] = FIND_CODE_HMF(image, size, window, pitch, vec_m, vec_n, coord_origin, gray_mean);
ext_m = zeros(1, window);
ext_n = zeros(1, window);
indent_m = 0;
indent_n = 0;
cur_m = 1;
cur_n = 1;
% debug = image;
for i = -(window - 1):(window - 1)
    if (i == 0)
        indent_m = cur_m;
        indent_n = cur_n;
    end
    [flag, value] = GET_CODE_HMF(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m + i * pitch, index_n);
    if (flag && cur_m <= window)
        ext_m(cur_m) = value;
        cur_m = cur_m + 1;
    end
    [flag, value] = GET_CODE_HMF(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m, index_n + i * pitch);
    if (flag && cur_n <= window)
        ext_n(cur_n) = value;
        cur_n = cur_n + 1;
    end
end
seq_m = D(ext_m);
seq_n = D(ext_n);
end

function [index_m, index_n] = FIND_CODE_HMF(image, size, window, pitch, vec_m, vec_n, coord_origin, gray_mean)
% mathematically, every window array contains a coded block.
index_m = 0;
index_n = 0;
for i = 0:pitch * window - 1 % window points
    for j = 0:pitch - 1  % pitch points
        [bound, value] = GET_CODE_HMF(image, size, vec_m, vec_n, coord_origin, gray_mean, i, j);
        if (bound && value)
            index_m = mod(i, pitch);
            index_n = mod(j, pitch);
            return
        end
        [bound, value] = GET_CODE_HMF(image, size, vec_m, vec_n, coord_origin, gray_mean, -i, j);
        if (bound && value)
            index_m = mod(-i, pitch);
            index_n = mod(j, pitch);
            return
        end
        [bound, value] = GET_CODE_HMF(image, size, vec_m, vec_n, coord_origin, gray_mean, i, -j);
        if (bound && value)
            index_m = mod(i, pitch);
            index_n = mod(-j, pitch);
            return
        end
        [bound, value] = GET_CODE_HMF(image, size, vec_m, vec_n, coord_origin, gray_mean, i, -j);
        if (bound && value)
            index_m = mod(-i, pitch);
            index_n = mod(-j, pitch);
            return
        end
    end
end
disp("Search failure, please increase sample area to include window * (pitch + 1) blocks")
end

function [bound, value] = GET_CODE_HMF(image, size, vec_m, vec_n, coord_origin, gray_mean, index_m, index_n)
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

% D morphism
function [S] = D(X)
S = [];
for i = 1:length(X) - 1
    if X(i) == X(i + 1)
        S = [S, 0];
    else
        S = [S, 1];
    end
end
end

% Search function
function [flag_m_is_x, index_x, flip_x, index_y, flip_y] = DECODE_HMF(ref_x, ref_y, len_x, len_y, n, seq_m, seq_m_f, seq_n, seq_n_f)
flag_m_is_x = false;
for i = 1:len_x - n + 1
    if isequal(seq_m, ref_x(i:i + n - 1))
        flip_x = false;
        index_x = i;
        flag_m_is_x = true;
        break;
    end
    if isequal(seq_m_f, ref_x(i:i + n - 1))
        flip_x = true;
        index_x = i;
        flag_m_is_x = true;
        break;
    end
end
if flag_m_is_x == false
    for i = 1:len_x - n + 1
        if isequal(seq_n, ref_x(i:i + n - 1))
            flip_x = false;
            index_x = i;
            break;
        end
        if isequal(seq_n_f, ref_x(i:i + n - 1))
            flip_x = true;
            index_x = i;
            break;
        end
    end
    for i = 1:len_y - n + 1
        if isequal(seq_m, ref_y(i:i + n - 1))
            flip_y = false;
            index_y = i;
            return;
        end
        if isequal(seq_m_f, ref_y(i:i + n - 1))
            flip_y = true;
            index_y = i;
            return;
        end
    end
else
    for i = 1:len_y - n + 1
        if isequal(seq_n, ref_y(i:i + n - 1))
            flip_y = false;
            index_y = i;
            return;
        end
        if isequal(seq_n_f, ref_y(i:i + n - 1))
            flip_y = true;
            index_y = i;
            return;
        end
    end
end
end

% Decode sequence
function [period_x, period_y, est_r] = LOCATE_HMF(ref_x, ref_y, len_x, len_y, window, pitch, bias_m, bias_n, seq_m, seq_n, index_m, index_n, indent_m, indent_n, theta)
seq_m_f = flip(seq_m);
seq_n_f = flip(seq_n);
% find the axes and directions A and B represent
[flag_m_is_x, index_x, flip_x, index_y, flip_y] = DECODE_HMF(ref_x, ref_y, len_x, len_y, window - 1, seq_m, seq_m_f, seq_n, seq_n_f);
if (flag_m_is_x)
    if (flip_x)
        period_x = ((index_x + (window - indent_m)) * pitch + index_m) / 2;
        period_x = period_x - bias_m;
    else
        period_x = ((index_x + (indent_m - 1)) * pitch - index_m) / 2;
        period_x = period_x + bias_m;
    end
    if (flip_y)
        period_y = ((index_y + (window - indent_n)) * pitch + index_n) / 2;
        period_y = period_y - bias_n;
    else
        period_y = ((index_y + (indent_n - 1)) * pitch - index_n) / 2;
        period_y = period_y + bias_n;
    end
else
    if (flip_x)
        period_x = ((index_x + (window - indent_n)) * pitch + index_n) / 2;
        period_x = period_x - bias_n;
    else
        period_x = ((index_x + (indent_n - 1)) * pitch - index_n) / 2;
        period_x = period_x + bias_n;
    end
    if (flip_y)
        period_y = ((index_y + (window - indent_m)) * pitch + index_m) / 2;
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

function [] = PREVIEW(image, coord_m, coord_n, ext_m, ext_n, bias_m, bias_n, vec_m, vec_n, index_m, index_n)
coord_origin = [(size(image, 1) + 1) / 2, (size(image, 2) + 1) / 2];
fig = figure(1);
fig.Position = [200, 200, 800, 800];
axs = subplot(1, 1, 1);
axs.Position = [0, 0, 1, 1];
set(0, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultLegendInterpreter', 'latex');
set(0,'defaultAxesFontSize',20)
hold on
box on
colormap gray
image(1, 1) = 5;
imagesc(abs(image))
axis([0, size(image, 1) + 1, 0, size(image, 2) + 1])

xticks([]);
xticklabels([]);
yticks([]);
yticklabels([])
vec_b =  round(- bias_m .* vec_m - bias_n .* vec_n);
vec_i = round(vec_m ./2 .* index_m + vec_n ./2 .* index_n);
p = nsidedpoly(32, 'Center', [coord_origin(2), coord_origin(1)], 'SideLength', 1);
plot(p, "EdgeColor", "#888", "FaceColor", "none", "LineWidth", 2);
quiver(coord_origin(2), coord_origin(1), vec_b(2), vec_b(1), 1, 'Color', "#fff", "LineWidth", 2, 'MaxHeadSize',10);
quiver(coord_origin(2) + vec_b(2), coord_origin(1) + vec_b(1), vec_i(2), vec_i(1), 1, 'Color', "#fff", "LineWidth", 2, 'MaxHeadSize',0.5);
quiver(coord_m(1, 2), coord_m(1, 1), coord_m(end, 2) - coord_m(1, 2), coord_m(end, 1) - coord_m(1, 1), 1, 'Color', "#fff", "LineWidth", 2);
quiver(coord_n(1, 2), coord_n(1, 1), coord_n(end, 2) - coord_n(1, 2), coord_n(end, 1) - coord_n(1, 1), 1, 'Color', "#fff", "LineWidth", 2);
for i=1:size(coord_m, 1)
    if (ext_m(i) == 1)
        p = nsidedpoly(3, 'Center', [coord_m(i, 2), coord_m(i, 1)], 'SideLength', 15);
        p = rotate(p, 180, [coord_m(i, 2), coord_m(i, 1)]);
        text(coord_m(i, 2) + 25, coord_m(i, 1), "$\textbf{1}$", 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Units', 'data', 'Color', 'white');
        plot(p, "EdgeColor", "#fff", "FaceColor", "none", "LineWidth", 2);
    else
        p = nsidedpoly(3, 'Center', [coord_m(i, 2), coord_m(i, 1)], 'SideLength', 15);
        text(coord_m(i, 2) + 25, coord_m(i, 1), "$\textbf{0}$", 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Units', 'data', 'Color', "#aaa");
        plot(p, "EdgeColor", "#888", "FaceColor", "none", "LineWidth", 2);
    end
end
for i=1:size(coord_n, 1)
    if (ext_n(i) == 1)
        p = nsidedpoly(3, 'Center', [coord_n(i, 2), coord_n(i, 1)], 'SideLength', 15);
        p = rotate(p, 180, [coord_n(i, 2), coord_n(i, 1)]);
        text(coord_n(i, 2) + 25, coord_n(i, 1), "$\textbf{1}$", 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Units', 'data', 'Color', 'white');
        plot(p, "EdgeColor", "#fff", "FaceColor", "none", "LineWidth", 2);
    else
        p = nsidedpoly(3, 'Center', [coord_n(i, 2), coord_n(i, 1)], 'SideLength', 15);
        text(coord_n(i, 2) + 25, coord_n(i, 1), "$\textbf{0}$", 'Interpreter', 'latex', 'FontSize', 20, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Units', 'data', 'Color', "#aaa");
        plot(p, "EdgeColor", "#888", "FaceColor", "none", "LineWidth", 2);
    end
end
end

