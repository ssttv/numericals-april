function u = conj_grad(M, N, H, R, v, hx, hy)
u = zeros(M, N);
r = R;
p = r;
e = 1;
norm2r = 0;
for i = 1 : M
    for j = 1 : N
        norm2r = norm2r + r(i, j) ^ 2;
    end
end
iter = 0;
while e > 3e-5
    iter = iter + 1;
    Hp_p = 0;
    Hp = zeros(M, N);
    for i = 1 : M
        for j = 1 : N
            if (i > 1)
                Hp(i, j) = Hp(i, j) + H(1, i, j) * p(i - 1, j);
            end
            if (j > 1)
                Hp(i, j) = Hp(i, j) + H(2, i, j) * p(i, j - 1);
            end
            if (j < N)
                Hp(i, j) = Hp(i, j) + H(4, i, j) * p(i, j + 1);
            end
            if (i < M)
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
    for i = 1 : M
        for j = 1 : N
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