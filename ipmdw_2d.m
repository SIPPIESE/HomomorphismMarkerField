function [k_est, l_est, phase] = ipmdw_2d(image, k, l)
[N, M] = size(image);
H = 2;
m = (0:1:M-1)/M;
n = (N-1:-1:0)'/N;
wx = zeros(1, M);
wy = zeros(N, 1);
%创建窗
for h = 0:1:H-1
    if h==0
        alpha = nchoosek(2*H-2, H-1)/2^(2*H-2);
    else
        alpha = nchoosek(2*H-2, H-h-1)/2^(2*H-3);
    end
    wx = wx+(-1)^h*alpha*cos(2*pi*h*m);
    wy = wy+(-1)^h*alpha*cos(2*pi*h*n);
end
w = wy*wx;
img_win = w.*image;
x = double((0:1:M-1)');
y = double(N-1:-1:0);
S_M = exp(-1j*2*pi/M);
S_x = S_M .^ x;
S_x_l = (S_x .^ (k - 1));
S_x_o = (S_x .^ (k));
S_x_u = (S_x .^ (k + 1));
S_N = exp(-1j*2*pi/N);
S_y = S_N .^ y;
S_y_l = (S_y .^ (l - 1));
S_y_o = (S_y .^ (l));
S_y_u = (S_y .^ (l + 1));
T_x = img_win*S_x_o;
T_y = S_y_o*img_win;
F_o_o = S_y_o*T_x;
F_l_o = T_y*S_x_l;
F_u_o = T_y*S_x_u;
F_o_l = S_y_l*T_x;
F_o_u = S_y_u*T_x;
Mat1 = [(2*H-1)*H, 2*H-1, F_l_o-F_o_o;
    -k^2-H^2, 2*k, F_o_o; 
    (2*H-1)*H, -(2*H-1), F_u_o-F_o_o];
Mat2 = [(2*H-1), F_l_o-F_o_o;
    -(2*H-1), F_u_o-F_o_o];
k_est = real(sqrt(det(Mat1)/det(Mat2)));
if (k < 0)
    k_est = - k_est;
end
if (k == 0)
    if abs(sum(F_u_o, 'all')) < abs(sum(F_l_o, 'all'))
        k_est = - k_est;
    end
end
Mat1 = [(2*H-1)*H, 2*H-1, F_o_l-F_o_o;
    -l^2-H^2, 2*l, F_o_o;
    (2*H-1)*H, -(2*H-1), F_o_u-F_o_o];
Mat2 = [(2*H-1), F_o_l-F_o_o;
    -(2*H-1), F_o_u-F_o_o];
l_est = real(sqrt(det(Mat1)/det(Mat2)));
if (l < 0)
    l_est = - l_est;
end
if (l == 0)
    if abs(sum(F_o_u, 'all')) < abs(sum(F_o_l, 'all'))
        l_est = - l_est;
    end
end
img_hann = (hann(N) * hann(M)') .* image;
x_m = double((0:1:M-1)' - (M-1)/2);
y_m = double((N-1:-1:0) - (N-1)/2);
phase = angle((S_N .^ y_m .^ l_est)*img_hann*(S_M .^ x_m .^ k_est));
end