clear;
clc;
seq_data_folder = "data/aos/";
hmf_data_folder = "data/hmf/";
% This function constructs 2D self-location & self-orientation coding
% pattern from the oriented half debruijn sequence, saving the X, Y and
% coding matrix G to /data.
% window in this file refers to the window array size after inverse D
% morphism, which is sequence window array size + 1.

for window_size = 4:16
    CONSTRUCT(window_size, seq_data_folder, hmf_data_folder);
end

function [] = CONSTRUCT(window_size, seq_data_folder, hmf_data_folder)
n = window_size - 1;
S = readmatrix(seq_data_folder + string(n) + '.txt');
a = floor(length(S) / 2) + floor((n - 1) / 2);

X = S(1:a);
Y = S(length(S) - a + 1:length(S));

X_0 = INVERSE_D(X, 0);
X_1 = INVERSE_D(X, 1);
Y_0 = INVERSE_D(Y, 0);

G = [];

for y = Y_0
    if y == 0
        G = [G, X_0];
    else
        G = [G, X_1];
    end
end

G = reshape(G, a + 1, a + 1);

writematrix(X, hmf_data_folder + string(window_size) + '+' + string(window_size) + '_X.txt');
writematrix(Y, hmf_data_folder + string(window_size) + '+' + string(window_size) + '_Y.txt');
writematrix(G, hmf_data_folder + string(window_size) + '+' + string(window_size) + '_G.txt');
end

function [S] = INVERSE_D(X, START)
S = START;
for x = X
    if x == 1
        if S(end) == 1
            S = [S, 0];
        else
            S = [S, 1];
        end
    else
        if S(end) == 1
            S = [S, 1];
        else
            S = [S, 0];
        end
    end
end
end