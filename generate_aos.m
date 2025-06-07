clear;
clc;
seq_data_folder = "data/aos/";
% This function generates oriented half debruijn sequence, saving the
% result to data/aos/n.txt
% In this file, n refers to the window array size in the generated 1-D
% orientable de Bruijn sequence.
n_max = 8;

for n = 3:n_max
    tic;
    S = INITIATE(n);
    toc
    disp(join(string(S)));
    writematrix(S, seq_data_folder + string(n) + '.txt');
end

function [S] = INITIATE(n)
global H_MAX;
global T_MAX;
T_MAX = 0;
X = [zeros(1, n-1), 1];
H = X;
H = SEARCH(X, H, 1);
S = [H_MAX(1, 1:n-1), H_MAX(:, size(H_MAX, 2))'];
end

function [H, L] = SEARCH(X, H, L)
global T_MAX;
T_MAX = toc;
global H_MAX;
Y = X(2:length(X));
X_0 = [Y, 0];
X_1 = [Y, 1];
H_0 = H;
H_1 = H;
L_0 = L;
L_1 = L;
if T_MAX > 240
    return;
end
if isequal(X_0, flip(X_0)) == false
    if FIND(X_0, H, L) == false
        if FIND(flip(X_0), H, L) == false
            [H_0, L_0] = SEARCH(X_0, [H; X_0], L + 1);
        end
    end
end
if T_MAX > 240
    return;
end
if isequal(X_1, flip(X_1)) == false
    if FIND(X_1, H, L) == false
        if FIND(flip(X_1), H, L) == false
            [H_1, L_1] = SEARCH(X_1, [H; X_1], L + 1);
        end
    end
end

if T_MAX > 240
    return;
end
if L_0 > L_1
    H = H_0;
    L = L_0;
    H_MAX = H;
    disp(L);
else
    H = H_1;
    L = L_1;
    H_MAX = H;
    disp(L);
end
end

function [T] = FIND(X, H, L)
    T = false;
    for i = 1:L
        if isequal(X, H(i, :))
            T = true;
            return
        end
    end
end