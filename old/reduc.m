function x = reduc(M, N, u)
x = zeros((M + 1) / 2, (N + 1) / 2);
for i = 1 : (M + 1) / 2
    for j = 1 : (N + 1) / 2
        if (i == 1 || j == 1 || i == (M + 1) / 2 || j == (N + 1) / 2)
            x(i, j) = u(i * 2 - 1, j * 2 - 1);
            continue;
        end
        x(i, j) = 0.25 * u(i * 2 - 1, j * 2 - 1);
        x(i, j) = x(i, j) + 0.125 * u(i * 2 - 2, j * 2 - 1);
        x(i, j) = x(i, j) + 0.0625 * u(i * 2 - 2, j * 2 - 2);
        x(i, j) = x(i, j) + 0.0625 * u(i * 2 - 2, j * 2);
        x(i, j) = x(i, j) + 0.125 * u(i * 2 - 1, j * 2 - 2);
        x(i, j) = x(i, j) + 0.125 * u(i * 2 - 1, j * 2);
        x(i, j) = x(i, j) + 0.125 * u(i * 2, j * 2 - 1);
        x(i, j) = x(i, j) + 0.0625 * u(i * 2, j * 2 - 2);
        x(i, j) = x(i, j) + 0.0625 * u(i * 2, j * 2);
    end
end