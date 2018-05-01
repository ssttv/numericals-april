clearvars
clc
rdown = 0;
rup = 0;
rleft = 0;
rright = 0;
sdown = 0.25;
sup = 1 / sqrt(7);
sleft = 1 / cos(1);
sright = 1 / -cos(3);
x1 = 1;
x2 = 3;
y1 = 0;
y2 = 3;
M = 65;
N = 129;
hx = (x2 - x1) / (M - 1);
hy = (y2 - y1) / (N - 1);
X = x1 : hx : x2;
Y = y1 : hy : y2;
K = 4;
nu = 50;
H_array = cell(K + 1, 1);
[H_array{1, 1}, R] = create_matrix(x1, x2, y1, y2, M, N, @a, @b, @q, @f, @phi_down, @phi_up, @phi_left, @phi_right, rdown, rup, rleft, rright, sdown, sup, sleft, sright);
for k = 1 : K
    [H_array{k + 1, 1}, ~] = create_matrix(x1, x2, y1, y2, M, N, @a, @b, @q, @f, @phi_down, @phi_up, @phi_left, @phi_right, rdown, rup, rleft, rright, sdown, sup, sleft, sright);
end
v = zeros(M, N);
for i = 1 : M
    for j = 1 : N
        v(i, j) = cos(1 + hx * (i - 1)) * sqrt(16 - (hy * (j - 1)) ^ 2);
    end
end

% for k = 1 : 30
%     iter(1, k) = 0;
%     nu(1, k) = 5 * k;
%     e = 1;
%     u = zeros(M, N);
%     while e > 3e-5
%         u = multigrid(M, N, H_array, R, u, nu(1, k), 0, K, @jacobi, @reduc, @interp);
%         iter(1, k) = iter(1, k) + 1;
%         e = 0;
%         for i = 1 : M
%             for j = 1 : N
%                 e = e + (v(i, j) - u(i, j)) ^ 2;
%             end
%         end
%         e = sqrt(e * hx * hy);
%     end
% end
% for k = 1 : 30
%     iter(2, k) = 0;
%     nu(2, k) = 5 * k;
%     e = 1;
%     u = zeros(M, N);
%     while e > 3e-5
%         u = multigrid(M, N, H_array, R, u, nu(2, k), 0, K, @seidel, @reduc, @interp);
%         iter(2, k) = iter(2, k) + 1;
%         e = 0;
%         for i = 1 : M
%             for j = 1 : N
%                 e = e + (v(i, j) - u(i, j)) ^ 2;
%             end
%         end
%         e = sqrt(e * hx * hy);
%     end
% end
% for k = 3 : 20
%     iter(3, k) = 0;
%     nu(3, k) = 5 * k;
%     e = 1;
%     u = zeros(M, N);
%     while e > 3e-5
%         u = multigrid(M, N, H_array, R, u, nu(3, k), 0, K, @relax, @reduc, @interp);
%         iter(3, k) = iter(3, k) + 1;
%         e = 0;
%         for i = 1 : M
%             for j = 1 : N
%                 e = e + (v(i, j) - u(i, j)) ^ 2;
%             end
%         end
%         e = sqrt(e * hx * hy);
%     end
% end
% figure
% hold on
% grid on
% plot(nu(1, :), iter(1, :), nu(2, :), iter(2, :), nu(3, 3 : 20), iter(3, 3 : 20))
% xlabel('число сглаживающих итераций')
% ylabel('число многосеточных итераций')
% legend('jacobi', 'seidel', 'relax')

% for k = 1 : 4
%     M(k) = 2 ^ (k + 3) + 1;
%     N(k) = 2 ^ (k + 4) + 1;
%     H_array = cell(K + 1, 1);
%     [H_array{1, 1}, R] = create_matrix(x1, x2, y1, y2, M(k), N(k), @a, @b, @q, @f, @phi_down, @phi_up, @phi_left, @phi_right, rdown, rup, rleft, rright, sdown, sup, sleft, sright);
%     for l = 1 : K
%         [H_array{l + 1, 1}, ~] = create_matrix(x1, x2, y1, y2, M(k), N(k), @a, @b, @q, @f, @phi_down, @phi_up, @phi_left, @phi_right, rdown, rup, rleft, rright, sdown, sup, sleft, sright);
%     end
%     iter(k) = 0;
%     e = 1;
%     u = zeros(M(k), N(k));
%     v = zeros(M(k), N(k));
%     hx = (x2 - x1) / (M(k) - 1);
%     hy = (y2 - y1) / (N(k) - 1);
%     for i = 1 : M(k)
%         for j = 1 : N(k)
%             v(i, j) = cos(1 + hx * (i - 1)) * sqrt(16 - (hy * (j - 1)) ^ 2);
%         end
%     end
%     while e > 1e-3
%         u = multigrid(M(k), N(k), H_array, R, u, nu, 0, K, @seidel, @reduc, @interp);
%         iter(k) = iter(k) + 1;
%         e = 0;
%         for i = 1 : M(k)
%             for j = 1 : N(k)
%                 e = e + (v(i, j) - u(i, j)) ^ 2;
%             end
%         end
%         e = sqrt(e * hx * hy);
%     end
% end
% figure
% hold on
% grid on
% plot(M .* N, iter);
% xlabel('число узлов')
% ylabel('число многосеточных итераций')

iter = 0;
e = 1;
u = zeros(M, N);
while e > 3e-5
    u = multigrid(M, N, H_array, R, u, nu, 0, K, @seidel, @reduc, @interp);
    iter = iter + 1;
    e(iter) = 0;
    for i = 1 : M
        for j = 1 : N
            e(iter) = e(iter) + (v(i, j) - u(i, j)) ^ 2;
        end
    end
    e(iter) = sqrt(e(iter) * hx * hy);
    if iter > 1
        ro(iter) = e(iter) / e(iter - 1);
    end
end
% figure
% semilogy(1 : iter, e);
% hold on
% grid on
% xlabel('номер итерации')
% ylabel('погрешность')
% figure
% hold on
% grid on
% plot(2 : iter, ro(2 : iter));
% xlabel('номер итерации')
% ylabel('ro')
%u = conj_grad(M, N, H_array{1, 1}, R, v, hx, hy);
% figure
% mesh(Y, X, u)
% figure
% mesh(Y, X, v)