function [vec_m, vec_n, bias_m, bias_n, theta] = FUNC_EMBEDDED_ESTIMATE(image)
freq = abs(fftshift(fft2(image)));
[k_o, l_o] = find(freq == max(max(freq)));
freq(k_o(1) - 1:k_o(1) + 1, l_o(1) - 1:l_o(1) + 1) = 0;
[k_a, l_a] = find(freq == max(max(freq)));
freq(k_a(1) - 1:k_a(1) + 1, l_a(1) - 1:l_a(1) + 1) = 0;
freq(k_a(2) - 1:k_a(2) + 1, l_a(2) - 1:l_a(2) + 1) = 0;
[k_b, l_b] = find(freq == max(max(freq)));
k_a = k_a - k_o;
l_a = l_a - l_o;
k_b = k_b - k_o;
l_b = l_b - l_o;
if (k_a(1) * l_a(1) >= 0)
    k_1 = abs(k_a(1));
    l_1 = abs(l_a(1));
    k_4 = abs(k_b(1));
    l_4 = - abs(l_b(1));
else
    k_1 = abs(k_b(1));
    l_1 = abs(l_b(1));
    k_4 = abs(k_a(1));
    l_4 = - abs(l_a(1));
end
[k_1, l_1, p_1] = ipmdw_2d(image, k_1, l_1);
[k_4, l_4, p_4] = ipmdw_2d(image, k_4, l_4);
r = k_1 * l_4 - k_4 * l_1;
v_1 = [l_4 * size(image, 1), - k_4 * size(image, 2)] ./ r;
v_4 = [- l_1 * size(image, 1), k_1 * size(image, 2)] ./ r;
vec_m = (v_1 + v_4);
vec_n = (v_1 - v_4);
p_m = (p_4 - p_1) / 2;
p_n = (p_4 + p_1) / 2;
bias_m = p_m / (2 * pi);
bias_n = p_n / (2 * pi);
theta = angle(vec_m(1) + 1i * vec_m(2));
end