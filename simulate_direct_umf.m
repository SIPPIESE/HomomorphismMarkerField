clear;
clc;
close all;
umf_data_folder = "data/umf/";
%% parameters for comparing decoding time
parameter = [
    4, 4, 16;
    4, 4, 32;
    4, 4, 64;
];
test_count = 10000;
%% simulations
% warm up
TEST_UMF([4, 4, 16], umf_data_folder, 1000);

time = zeros(size(parameter, 1), 1);
for i = 1:size(parameter, 1)
    [time(i)] = TEST_UMF(parameter(i, :), umf_data_folder, test_count);
end


function [T] = TEST_UMF(parameter, umf_data_folder, test_count)
w = parameter(1:2);
S = parameter(3);
G = readmatrix(umf_data_folder + string(w(1)) + "x" + string(w(2)) + '_G.txt') == 1;
if (S == -1)
    S = size(G, 1);
else
    G = G(1:S, 1:S);
end

disp("umf >> " + string(S));
success_count = 0;
tot_time = 0;
for i = 1:test_count
    % sample A and B
    [x, y, r, G_S] = SAMPLE(S, G, w);
    % start decoding process
    tic;
    [x_est, y_est, r_est] = DECODE(S, G_S, G, w);
    % add trial time to tot_time
    tot_time = tot_time + toc;
    if (i == 1)
        disp("umf >> umf order - " + string(w(1)) + "x" + string(w(2)) + " array size - " + string(size(G_S, 1)) + "x" + string(size(G_S, 2)))
    end
    % decide whether the measurement is accurate, add success count
    if isequal([x, y, r], [x_est, y_est, r_est])
        success_count = success_count + 1;
    end
end
disp("umf >> the average decoding time is:                " + string(tot_time / test_count) + "s");
disp("umf >> the success count versus total trials is:    " + string(success_count) + "/" + string(test_count));
if success_count == test_count
    disp("umf >> success")
else
    disp("umf >> failure")
end
T = tot_time / test_count;
end

function [I_X, I_Y, R] = DECODE(S, G_S, G, w)
G_S_90 = rot90(G_S, -1);
G_S_180 = rot90(G_S, -2);
G_S_270 = rot90(G_S, -3);
% find the axes and directions A and B represent
[I_X, I_Y, R] = FIND(S, G_S, G_S_90, G_S_180, G_S_270, G, w);
end

function [x, y, r, G_S] = SAMPLE(S, G, w)
% select a point
% x_R, y_R is the index of the point after rotation
x = randi([1, S + 1 - w(1)]);
y = randi([1, S + 1 - w(2)]);
G_S = G(x:x + w(1) - 1, y:y + w(2) - 1);
r = randi([0, 3]);
for j = 1:r
    G_S = rot90(G_S);
end
r = r * 90;
end

%the function of D-mapping
function [I_X, I_Y, R] = FIND(S, G_S, G_S_90, G_S_180, G_S_270, G, w)
for i = 1:S - w(1) + 1
    for j = 1:S - w(2) + 1
        G_P = G(i:i + w(1) - 1, j:j + w(2) - 1);
        if isequal(G_S, G_P)
            I_X = i;
            I_Y = j;
            R = 0;
            return;
        end
        if isequal(G_S_90, G_P)
            I_X = i;
            I_Y = j;
            R = 90;
            return;
        end
        if isequal(G_S_180, G_P)
            I_X = i;
            I_Y = j;
            R = 180;
            return;
        end
        if isequal(G_S_270, G_P)
            I_X = i;
            I_Y = j;
            R = 270;
            return;
        end
    end
end
end

