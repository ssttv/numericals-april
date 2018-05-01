function [H, R] = create_matrix(x1, x2, y1, y2, M, N, a, b, q, f, phi_down, phi_up, phi_left, phi_right, rdown, rup, rleft, rright, sdown, sup, sleft, sright)
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
for i = 2 : M - 1
    for j = 2 : N - 1
        H(1, i, j) = -px(i - 1, j) * hy_hx;
        H(2, i, j) = -py(i, j - 1) / hy_hx;
        H(4, i, j) = -py(i, j) / hy_hx;
        H(5, i, j) = -px(i, j) * hy_hx;
        H(3, i, j) = q(x1 + (i - 1) * hx, y1 + (j - 1) * hy) * hx * hy + px(i - 1, j) * hy_hx + px(i, j) * hy_hx + py(i, j - 1) / hy_hx + py(i, j) / hy_hx;
        R(i, j) = f(x1 + (i - 1) * hx, y1 + (j - 1) * hy) * hx * hy;
    end
end
%down
for i = 2 : M - 1
    if (rdown == 0)
        H(3, i, 1) = sdown;
        R(i, 1) = phi_down(x1 + (i - 1) * hx);
    else
        pborder = b(x1 + (i - 1) * hx, y1) * hx;
        H(1, i, 1) = -px(i - 1, 1);
        H(4, i, 1) = q(x1 + (i - 1) * hx, y1 + hy) / 4 * hx * hy - 2 * py(i, 1);
        H(5, i, 1) = -px(i, 1);
        H(3, i, 1) = 3 * q(x1 + (i - 1) * hx, y1) / 4 * hx * hy + 2 * pborder * sdown / rdown + px(i - 1, 1) + px(i, 1) + 2 * py(i, 1);
        R(i, 1) = f(x1 + (i - 1) * hx, y1 + hy / 4) * hx * hy + 2 * pborder * phi_down(x1 + (i - 1) * hx) / rdown;
    end
end
%up
for i = 2 : M - 1
    if (rup == 0)
        H(3, i, N) = sup;
        R(i, N) = phi_up(x1 + (i - 1) * hx);
    else
        pborder = b(x1 + (i - 1) * hx, y2) * hx;
        H(1, i, N) = -px(i - 1, N);
        H(2, i, N) = q(x1 + (i - 1) * hx, y2 - hy) / 4 * hx * hy - 2 * py(i, N - 1);
        H(5, i, N) = -px(i, N);
        H(3, i, N) = 3 * q(x1 + (i - 1) * hx, y2) / 4 * hx * hy + 2 * pborder * sup / rup + px(i - 1, N) + px(i, N) + 2 * py(i, N - 1);
        R(i, N) = f(x1 + (i - 1) * hx, y2 - hy / 4) * hx * hy + 2 * pborder * phi_up(x1 + (i - 1) * hx) / rup;
    end
end
%left
for j = 2 : N - 1
    if (rleft == 0)
        H(3, 1, j) = sleft;
        R(1, j) = phi_left(y1 + (j - 1) * hy);
    else
        pborder = a(x1, y1 + (j - 1) * hy) * hy;
        H(2, 1, j) = -py(1, j - 1);
        H(4, 1, j) = -py(1, j);
        H(5, 1, j) = q(x1 + hx, y1 + (j - 1) * hy) / 4 * hx * hy - 2 * px(1, j);
        H(3, 1, j) = 3 * q(x1, y1 + (j - 1) * hy) / 4 * hx * hy + 2 * pborder * sleft / rleft + 2 * px(1, j) + py(1, j) + py(1, j - 1);
        R(1, j) = f(x1 + hx / 4, y1 + (j - 1) * hy) * hx * hy + 2 * pborder * phi_left(y1 + (j - 1) * hy) / rleft;
    end
end
%right
for j = 2 : N - 1
    if (rright == 0)
        H(3, M, j) = sright;
        R(M, j) = phi_right(y1 + (j - 1) * hy);
    else
        pborder = a(x2, y1 + (j - 1) * hy) * hy;
        H(2, M, j) = -py(M, j - 1);
        H(4, M, j) = -py(M, j);
        H(1, M, j) = q(x2 - hx, y1 + (j - 1) * hy) / 4 * hx * hy - 2 * px(M - 1, j);
        H(3, M, j) = 3 * q(x2, y1 + (j - 1) * hy) / 4 * hx * hy + 2 * pborder * sright / rright + 2 * px(M - 1, j) + py(M, j) + py(M, j - 1);
        R(M, j) = f(x2 - hx / 4, y1 + (j - 1) * hy) * hx * hy + 2 * pborder * phi_right(y1 + (j - 1) * hy) / rright;
    end
end

if (rleft == 0)
    H(3, 1, 1) = sleft;
    R(1, 1) = phi_left(y1);
elseif (rdown == 0)
    H(3, 1, 1) = sdown;
    R(1, 1) = phi_down(x1);
else
    pl = a(x1, y1 + hy / 4) * hy;
    pd = b(x1 + hx / 4, y1) * hx;
    H(4, 1, 1) = q(x1, y1 + hy) / 4 * hx * hy - 2 * py(1, 1) + 2 * pl * sleft / rleft / 4;
    H(5, 1, 1) = q(x1 + hx, y1) / 4 * hx * hy - 2 * px(1, 1) + 2 * pd * sdown / rdown / 4;
    H(3, 1, 1) = q(x1, y1) / 2 * hx * hy + 2 * pl * sleft / rleft * 3 / 4 + 2 * pd * sdown / rdown * 3 / 4 + 2 * py(1, 1) + 2 * px(1, 1);
    R(1, 1) = f(x1 + hx / 4, y1 + hy / 4) * hx * hy + 2 * pl * phi_left(y1 + hy / 4) / rleft + 2 * pd * phi_down(x1 + hx / 4) / rdown;
end

if (rright == 0)
    H(3, M, 1) = sright;
    R(M, 1) = phi_right(y1);
elseif (rdown == 0)
    H(3, M, 1) = sdown;
    R(M, 1) = phi_down(x2);
else
    pr = a(x2, y1 + hy / 4) * hy;
    pd = b(x2 - hx / 4, y1) * hx;
    H(4, M, 1) = q(x2, y1 + hy) / 4 * hx * hy - 2 * py(M, 1) + 2 * pr * sright / rright / 4;
    H(1, M, 1) = q(x2 - hx, y1) / 4 * hx * hy - 2 * px(M - 1, 1) + 2 * pd * sdown / rdown / 4;
    H(3, M, 1) = q(x2, y1) / 2 * hx * hy + 2 * pr * sright / rright * 3 / 4 + 2 * pd * sdown / rdown * 3 / 4 + 2 * py(M, 1) + 2 * px(M - 1, 1);
    R(M, 1) = f(x2 - hx / 4, y1 + hy / 4) * hx * hy + 2 * pr * phi_right(y1 + hy / 4) / rright + 2 * pd * phi_down(x2 - hx / 4) / rdown ;
end

if (rleft == 0)
    H(3, 1, N) = sleft;
    R(1, N) = phi_left(y2);
elseif (rup == 0)
    H(3, 1, N) = sup;
    R(1, N) = phi_up(x1);
else
    pl = a(x1, y2 - hy / 4) * hy;
    pu = b(x1 + hx / 4, y2) * hx;
    H(2, 1, N) = q(x1, y2 - hy) / 4 * hx * hy - 2 * py(1, N - 1) + 2 * pl * sleft / rleft / 4;
    H(5, 1, N) = q(x1 + hx, y2) / 4 * hx * hy - 2 * px(1, N) + 2 * pu * sup / rup / 4;
    H(3, 1, N) = q(x1, y2) / 2 * hx * hy + 2 * pl * sleft / rleft * 3 / 4 + 2 * pu * sup / rup * 3 / 4 + 2 * py(1, N - 1) + 2 * px(1, N);
    R(1, N) = f(x1 + hx / 4, y2 - hy / 4) * hx * hy + 2 * pl * phi_left(y2 - hy / 4) / rleft + 2 * pu * phi_up(x1 + hx / 4) / rup;
end

if (rright == 0)
    H(3, M, N) = sright;
    R(M, N) = phi_right(y2);
elseif (rup == 0)
    H(3, M, N) = sup;
    R(M, N) = phi_up(x2);
else
    pr = a(x2, y2 - hy / 4) * hx;
    pu = b(x2 - hx / 4, y2) * hy;
    H(2, M, N) = q(x2, y2 - hy) / 4 * hx * hy - 2 * py(M, N - 1) + 2 * pr * sright / rright / 4;
    H(1, M, N) = q(x2 - hx, y2) / 4 * hx * hy - 2 * px(M - 1, N) + 2 * pu * sup / rup / 4;
    H(3, M, N) = q(x2, y2) / 2 * hx * hy + 2 * pr * sright / rright * 3 / 4 + 2 * pu * sup / rup * 3 / 4 + 2 * py(M, N - 1) + 2 * px(M - 1, N);
    R(M, N) = f(x2 - hx / 4, y2 - hy / 4) * hx * hy + 2 * pr * phi_right(y2 - hy / 4) / rright + 2 * pd * phi_up(x1 + hx / 4) / rup;
end