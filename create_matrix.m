function [H, R, u0] = create_matrix(x1, x2, y1, y2, M, N, a, b, q, f, phi_down, phi_up, phi_left, phi_right, rdown, rup, rleft, rright, sdown, sup, sleft, sright)
H = zeros(5, M, N);
R = zeros(M, N);
hx = (x2 - x1) / (M - 1);
hy = (y2 - y1) / (N - 1);
hy_hx = hy / hx;
px = zeros(M - 1, N);
py = zeros(M, N - 1);
for i = 1 : M - 1
    for j = 1 : N
        px(i, j) = a(x1 + (i - 0.5) * hx, y1 + (j - 1) * hy);
    end
end
for j = 1 : N - 1
    for i = 1 : M
        py(i, j) = b(x1 + (i - 1) * hx, y1 + (j - 0.5) * hy);
    end
end
u0 = zeros(M, N);
%down
for i = 1 : M
    u0(i, 1) = phi_down(x1 + (i - 1) * hx) / sdown;
end
%up
for i = 1 : M
    u0(i, N) = phi_up(x1 + (i - 1) * hx) / sup;
end
%left
for j = 2 : N - 1
    u0(1, j) = phi_left(y1 + (j - 1) * hy) / sleft;
end
%right
for j = 2 : N - 1
    u0(M, j) = phi_right(y1 + (j - 1) * hy) / sright;
end
for i = 2 : M - 1
    for j = 2 : N - 1
        if i > 2
            H(1, i, j) = -px(i - 1, j) * hy_hx;
        else
            R(i, j) = R(i, j) + px(i - 1, j) * hy_hx * u0(1, j);
        end
        if j > 2
            H(2, i, j) = -py(i, j - 1) / hy_hx;
        else
            R(i, j) = R(i, j) + py(i, j - 1) / hy_hx * u0(i, 1);
        end
        if j < N - 1
            H(4, i, j) = -py(i, j) / hy_hx;
        else
            R(i, j) = R(i, j) + py(i, j) / hy_hx * u0(i, N);
        end
        if i < M - 1
            H(5, i, j) = -px(i, j) * hy_hx;
        else
            R(i, j) = R(i, j) + px(i, j) * hy_hx * u0(M, j);
        end
        H(3, i, j) = q(x1 + (i - 1) * hx, y1 + (j - 1) * hy) * hx * hy + px(i - 1, j) * hy_hx + px(i, j) * hy_hx + py(i, j - 1) / hy_hx + py(i, j) / hy_hx;
        R(i, j) = R(i, j) + f(x1 + (i - 1) * hx, y1 + (j - 1) * hy) * hx * hy;
    end
end