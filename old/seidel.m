function u = seidel(M, N, H, R, uold, nu)
u = uold;
for n = 1 : nu
    uold = u;
    for i = 1 : M
        for j = 1 : N
            u(i, j) = R(i, j);
            if (i > 1)
                u(i, j) = u(i, j) - H(1, i, j) * u(i - 1, j);
            end
            if (j > 1)
                u(i, j) = u(i, j) - H(2, i, j) * u(i, j - 1);
            end
            if (j < N)
                u(i, j) = u(i, j) - H(4, i, j) * uold(i, j + 1);
            end
            if (i < M)
                u(i, j) = u(i, j) - H(5, i, j) * uold(i + 1, j);
            end
            u(i, j) = u(i, j) / H(3, i, j);
        end
    end
end