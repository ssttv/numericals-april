function x = reduc(M, N, u)
x = zeros((M + 1) / 2, (N + 1) / 2);
for i = 2 : (M + 1) / 2 - 1
    for j = 2 : (N + 1) / 2 - 1
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