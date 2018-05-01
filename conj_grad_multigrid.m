function [u, iter] = conj_grad_multigrid(M, N, H_array, R, u0, v, hx, hy, seidel, reduc, interp)
u = u0;
r = R;
q = multigrid(M, N, H_array, r, zeros(M, N), 20, 0, 3, seidel, reduc, interp);
p = q;
e = 1;
r_q = 0;
for i = 2 : M - 1
    for j = 2 : N - 1
        r_q = r_q + r(i, j) * q(i, j);
    end
end
iter = 0;
while e > 3e-5
    iter = iter + 1;
    Hp_p = 0;
    Hp = zeros(M, N);
    for i = 2 : M - 1
        for j = 2 : N - 1
            if (i > 2)
                Hp(i, j) = Hp(i, j) + H_array{1, 1}(1, i, j) * p(i - 1, j);
            end
            if (j > 2)
                Hp(i, j) = Hp(i, j) + H_array{1, 1}(2, i, j) * p(i, j - 1);
            end
            if (j < N - 1)
                Hp(i, j) = Hp(i, j) + H_array{1, 1}(4, i, j) * p(i, j + 1);
            end
            if (i < M - 1)
                Hp(i, j) = Hp(i, j) + H_array{1, 1}(5, i, j) * p(i + 1, j);
            end
            Hp(i, j) = Hp(i, j) + H_array{1, 1}(3, i, j) * p(i, j);
            Hp_p = Hp_p + Hp(i, j) * p(i, j);
        end
    end
    alpha = r_q / Hp_p;
    u = u + alpha * p;
    r = r - alpha * Hp;
    q = multigrid(M, N, H_array, r, zeros(M, N), 20, 0, 3, seidel, reduc, interp);
    r_qold = r_q;
    r_q = 0;
    for i = 2 : M - 1
        for j = 2 : N - 1
            r_q = r_q + r(i, j) * q(i, j);
        end
    end
    beta = r_q / r_qold;
    p = q + beta * p;
    e = 0;
    for i = 1 : M
        for j = 1 : N
            e = e + (v(i, j) - u(i, j)) ^ 2;
        end
    end
    e = sqrt(e * hx * hy);
end