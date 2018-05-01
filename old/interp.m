function u = interp(M, N, x)
u = zeros(2 * M - 1, 2 * N - 1);
for i = 1 : M
    for j = 1 : N
        u(i * 2 - 1, j * 2 - 1) = x(i, j);
        if (i > 1)
            u(i * 2 - 2, j * 2 - 1) = u(i * 2 - 2, j * 2 - 1) + 0.5 * x(i, j);
            if (j > 1)
                u(i * 2 - 2, j * 2 - 2) =  u(i * 2 - 2, j * 2 - 2) + 0.25 * x(i, j);
            end
            if (j < N)
                u(i * 2 - 2, j * 2) = u(i * 2 - 2, j * 2) + 0.25 * x(i, j);
            end
        end
        if (j > 1)
            u(i * 2 - 1, j * 2 - 2) = u(i * 2 - 1, j * 2 - 2) + 0.5 * x(i, j);
        end
        if (j < N)
            u(i * 2 - 1, j * 2) = u(i * 2 - 1, j * 2) + 0.5 * x(i, j);
        end
        if (i < M)
            u(i * 2, j * 2 - 1) = u(i * 2, j * 2 - 1) + 0.5 * x(i, j);
            if (j > 1)
                u(i * 2, j * 2 - 2) = u(i * 2, j * 2 - 2) + 0.25 * x(i, j);
            end
            if (j < N)
                u(i * 2, j * 2) = u(i * 2, j * 2) + 0.25 * x(i, j);
            end
        end
    end
end