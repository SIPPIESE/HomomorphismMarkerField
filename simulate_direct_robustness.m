clear;
clc;

window_range = 6:10;
error_range = 0.0:0.02:0.4;
R_matrix = [];

for i = window_range
    [R] = TEST(i, error_range);
    R_matrix = [R_matrix; R];
end

disp(R_matrix);

function [success_rate] = TEST(n, error_range)
X = readmatrix('data/hmf/' + string(n) + "+" + string(n) + '_X.txt');
Y = readmatrix('data/hmf/' + string(n) + "+" + string(n) + '_Y.txt');
G = readmatrix('data/hmf/' + string(n) + "+" + string(n) + '_G.txt');

G = G(1:n, 1:n); % choose the first window as window locations are irrelavant

test = 2000;
success_rate = [];
for e = error_range
    success_count = 0;
    total_count = 0;
    for i = 1:test
        p = rand([n, n]) <= e;
        g = abs(G - p);
        result_matrix = [];
        for x = 1:n
            for y = 1:n
                A = D(g(x, :));
                B = D(g(:, y)');
                try
                    [x_est, y_est, r_est] = DECODE(A, B, X, Y, n - 1);
                    if (isempty(result_matrix))
                        result_matrix = [x_est, y_est, r_est, 1];
                    else
                        flag_find = false;
                        for j = 1:size(result_matrix, 1)
                            temp = result_matrix(j, :);
                            if (temp(1:3) == [x_est, y_est, r_est])
                                result_matrix(j, 4) = result_matrix(j, 4) + 1;
                                flag_find = true;
                                break;
                            end
                        end
                        if (flag_find == false)
                            result_matrix = [result_matrix; [x_est, y_est, r_est, 1]];
                        end
                    end
                catch
                end
            end
        end
        count_max = 0;
        result_max = [];
        for j = 1:size(result_matrix, 1)
            if (count_max == 0)
                count_max = result_matrix(j, 4);
                result_max = result_matrix(j, 1:3);
            else
                if (result_matrix(j, 4) > count_max)
                    count_max = result_matrix(j, 4);
                    result_max = result_matrix(j, 1:3);
                end
            end
        end
        try
            if (result_max == [1, 1, 0])
                success_count = success_count + 1;
            end
        catch
        end
        total_count = total_count + 1;
    end
    success_rate = [success_rate, success_count / total_count];
    disp(success_rate)
end
end

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

function [x_est, y_est, r_est] = DECODE(A, B, X, Y, n)
    A_F = flip(A);
    B_F = flip(B);
    [I_X, f_X, I_Y, f_Y]= FIND(A, A_F, B, B_F, n, X, Y);
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

%the function of D-mapping
function [I_X, f_X, I_Y, f_Y] = FIND(A, A_F, B, B_F, n, X, Y)
A_flag = false;
L = length(X);
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

