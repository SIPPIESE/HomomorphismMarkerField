clear;
clc;
close all;
dbt_data_folder = "data/dbt/";
%% parameters for comparing decoding time
parameter = [
    3, 3, 16, 10000;
    4, 4, 32, 10000;
    4, 4, 64, 10000;
    4, 4, 128, 1000;
    4, 4, 256, 500;
    5, 5, 512, 250;
    5, 5, 1024, 100;
    5, 5, 2048, 50;
    5, 5, 4096, 25;
];
test_count = 1000;
%% simulations
% warm up
TEST_DBT([3, 3, -1, 10], dbt_data_folder, 1000);

time = zeros(size(parameter, 1), 1);
for i = 1:size(parameter, 1)
    [time(i)] = TEST_DBT(parameter(i, :), dbt_data_folder, test_count);
end


function [T] = TEST_DBT(parameter, dbt_data_folder, test_count)
w = parameter(1:2);
S = parameter(3);
test_count = parameter(4);
G = readmatrix(dbt_data_folder + string(w(1)) + "x" + string(w(2)) + '_G.txt') == 1;
if (S == -1)
    S = size(G, 1);
else
    G = G(1:S, 1:S);
end

disp("dbt >> " + string(S));
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
        disp("dbt >> dbt order - " + string(w(1)) + "x" + string(w(2)) + " array size - " + string(size(G_S, 1)) + "x" + string(size(G_S, 2)))
    end
    % decide whether the measurement is accurate, add success count
    if isequal([x, y, r], [x_est, y_est, r_est])
        success_count = success_count + 1;
    end
end
disp("dbt >> the average decoding time is:                " + string(tot_time / test_count) + "s");
disp("dbt >> the success count versus total trials is:    " + string(success_count) + "/" + string(test_count));
if success_count == test_count
    disp("dbt >> success")
else
    disp("dbt >> failure")
end
T = tot_time / test_count;
end

function [I_X, I_Y, R] = DECODE(S, G_S, G, w)
% find the axes and directions A and B represent
[I_X, I_Y, R] = FIND(S, G_S, G, w);
end

function [x, y, r, G_S] = SAMPLE(S, G, w)
% select a point
% x_R, y_R is the index of the point after rotation
x = randi([1, S + 1 - w(1)]);
y = randi([1, S + 1 - w(2)]);
G_S = G(x:x + w(1) - 1, y:y + w(2) - 1);
r = 0;
end

%the function of D-mapping
function [I_X, I_Y, R] = FIND(S, G_S, G, w)
for i = 1:S - w(1) + 1
    for j = 1:S - w(2) + 1
        G_P = G(i:i + w(1) - 1, j:j + w(2) - 1);
        if isequal(G_S, G_P)
            I_X = i;
            I_Y = j;
            R = 0;
            return;
        end
    end
end
end

