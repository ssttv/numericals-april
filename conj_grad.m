function [u, iter] = conj_grad(M, N, H, R, u0, v, hx, hy)
u = u0;
r = R;
p = r;
e = 1;
norm2r = 0;
for i = 2 : M - 1
    for j = 2 : N - 1
        norm2r = norm2r + r(i, j) ^ 2;
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
                Hp(i, j) = Hp(i, j) + H(1, i, j) * p(i - 1, j);
            end
            if (j > 2)
                Hp(i, j) = Hp(i, j) + H(2, i, j) * p(i, j - 1);
            end
            if (j < N - 1)
                Hp(i, j) = Hp(i, j) + H(4, i, j) * p(i, j + 1);
            end
            if (i < M - 1)
                Hp(i, j) = Hp(i, j) + H(5, i, j) * p(i + 1, j);
            end
            Hp(i, j) = Hp(i, j) + H(3, i, j) * p(i, j);
            Hp_p = Hp_p + Hp(i, j) * p(i, j);
        end
    end
    alpha = norm2r / Hp_p;
    u = u + alpha * p;
    r = r - alpha * Hp;
    norm2rold = norm2r;
    norm2r = 0;
    for i = 2 : M - 1
        for j = 2 : N - 1
            norm2r = norm2r + r(i, j) ^ 2;
        end
    end
    beta = norm2r / norm2rold;
    p = r + beta * p;
    e = 0;
    for i = 1 : M
        for j = 1 : N
            e = e + (v(i, j) - u(i, j)) ^ 2;
        end
    end
    e = sqrt(e * hx * hy);
end