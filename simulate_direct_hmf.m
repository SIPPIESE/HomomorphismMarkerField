clear;
clc;
close all;
hmf_data_folder = "data/hmf/";
%% parameters for direct decoding hmf
parameter = [
    6, -1;
    7, -1;
    8, -1;
    9, -1;
    10, -1;
    11, -1;
    12, -1;
    13, -1;
    14, -1;
    15, -1;
    16, -1;
];
test_count = 10000;
%% parameters for comparing decoding time
parameter = [
    7, 16;
    9, 32;
    10, 64;
    11, 128;
    12, 256;
    13, 512;
    14, 1024;
    15, 2048;
    16, 4096;
];
test_count = 10000;
%% simulations
% warm up
TEST_hmf([6, -1], hmf_data_folder, 1000);

time = zeros(size(parameter, 1), 1);
for i = 1:size(parameter, 1)
    [parameter(i, 2), time(i)] = TEST_hmf(parameter(i, :), hmf_data_folder, test_count);
end


function [S, T] = TEST_hmf(parameter, hmf_data_folder, test_count)
n = parameter(1);
S = parameter(2);
X = readmatrix(hmf_data_folder + string(n) + "+" + string(n) + '_X.txt') == 1;
Y = readmatrix(hmf_data_folder + string(n) + "+" + string(n) + '_Y.txt') == 1;
G = readmatrix(hmf_data_folder + string(n) + "+" + string(n) + '_G.txt') == 1;
if (S == -1)
    S = size(G, 1);
else
    G = G(1:S, 1:S);
end
L = length(X);

disp("hmf >> " + string(S));
success_count = 0;
tot_time = 0;
for i = 1:test_count
    % sample A and B
    [x, y, r, A, B] = SAMPLE(S, G, n);
    % start decoding process
    tic;
    [x_est, y_est, r_est] = DECODE(A, B, X, Y, L, n);
    % add trial time to tot_time
    tot_time = tot_time + toc;
    if (i == 1)
        disp("hmf >> hmf order - " + string(n) + " sequence length - " + string(length(A)))
    end
    % decide whether the measurement is accurate, add success count
    if isequal([x, y, r], [x_est, y_est, r_est])
        success_count = success_count + 1;
    end
end
disp("hmf >> the average decoding time is:                " + string(tot_time / test_count) + "s");
disp("hmf >> the success count versus total trials is:    " + string(success_count) + "/" + string(test_count));
if success_count == test_count
    disp("hmf >> success")
else
    disp("hmf >> failure")
end
T = tot_time / test_count;
end

function [S] = D(X, n)
S = xor(X(1:n - 1), X(2:n));
end

function [x_est, y_est, r_est] = DECODE(A, B, X, Y, L, n)
    A = D(A, n);
    B = D(B, n);
    A_F = flip(A);
    B_F = flip(B);
    % find the axes and directions A and B represent
    [I_X, f_X, I_Y, f_Y]= FIND(A, A_F, B, B_F, n - 1, X, Y, L);
    % derive orientation from the direction of axes
    r_est = 0;
    if isequal([0, 1], [f_X, f_Y])
        r_est = 90;
    elseif isequal([1, 0], [f_X, f_Y])
        r_est = 270;
    elseif isequal([1, 1], [f_X, f_Y])
        r_est = 180;
    end
    x_est = I_X;
    y_est = I_Y;
end

function [x, y, r, A, B] = SAMPLE(S, G, n)
    % select a point
    % x_R, y_R is the index of the point after rotation
    x = randi([1, S + 1 - n]);
    y = randi([1, S + 1 - n]);
    G_S = G(x:x + n - 1, y:y + n - 1);
    r = randi([0, 3]);
    for j = 1:r
        G_S = rot90(G_S);
    end
    r = r * 90;
    A = G_S(1, :);
    B = G_S(:, 1)';
end

%the function of D-mapping
function [I_X, f_X, I_Y, f_Y] = FIND(A, A_F, B, B_F, n, X, Y, L)
A_flag = false;
for i = 1:L - n + 1
    if isequal(A, X(i:i + n - 1))
        f_X = false;
        I_X = i;
        A_flag = true;
        break;
    end
    if isequal(A_F, X(i:i + n - 1))
        f_X = true;
        I_X = i;
        A_flag = true;
        break;
    end
end
if A_flag == false
    for i = 1:L - n + 1
        if isequal(B, X(i:i + n - 1))
            f_X = false;
            I_X = i;
            break;
        end
        if isequal(B_F, X(i:i + n - 1))
            f_X = true;
            I_X = i;
            break;
        end
    end
    for i = 1:L - n + 1
        if isequal(A, Y(i:i + n - 1))
            f_Y = false;
            I_Y = i;
            return;
        end
        if isequal(A_F, Y(i:i + n - 1))
            f_Y = true;
            I_Y = i;
            return;
        end
    end
else
    for i = 1:L - n + 1
        if isequal(B, Y(i:i + n - 1))
            f_Y = false;
            I_Y = i;
            return;
        end
        if isequal(B_F, Y(i:i + n - 1))
            f_Y = true;
            I_Y = i;
            return;
        end
    end
end
end

