clear;
clc;

load("mat/embedded_hmf_32_128_10000.mat");
time_hmf = time;
window_hmf = [parameter(:, 1), parameter(:, 1)];
field_size = parameter(:, 2);

load("mat/embedded_umf_32_72_10000.mat");
time_umf = time;
window_umf = [parameter(:, 1), parameter(:, 2)];

load("mat/embedded_dbt_32_128_10000.mat");
time_dbt = time;
window_dbt = [parameter(:, 1), parameter(:, 2)];

times = ' \times ';
for i = 1:size(field_size)
    f_s = "$" + string(field_size(i)) + times + string(field_size(i)) + "$";
    w_hmf = "$" + string(window_hmf(i, 1)) + "+" + string(window_hmf(i, 2)) + "$";
    t_hmf = "$" + num2str(time_hmf(i) * 1000, '%.3f') + "$";
    if (i <= size(window_umf, 1))
        w_umf = "$" + string(window_umf(i, 1)) + times + string(window_umf(i, 2)) + "$";
        t_umf = "$" + num2str(time_umf(i) * 1000, '%.3f') + "$";
    else
        w_umf = "/";
        t_umf = "/";
    end
    w_dbt = "$" + string(window_dbt(i, 1)) + times + string(window_dbt(i, 2)) + "$";
    t_dbt = "$" + num2str(time_dbt(i) * 1000, '%.3f') + "$";
    row_s = f_s + " & " + w_umf + " & " + t_umf + " & " + w_dbt + " & " + t_dbt + " & "  + w_hmf + " & " + t_hmf + "\\";
    disp(row_s)
end