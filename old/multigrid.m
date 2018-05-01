function v = multigrid(M, N, H_array, R, u, nu, k, K, smoother, reduc, interp)
if k < K
    v = smoother(M, N, H_array{k + 1, 1}, R, u, nu);
    r = R;
    for i = 1 : M
        for j = 1 : N
            if (i > 1)
                r(i, j) = r(i, j) - H_array{k + 1, 1}(1, i, j) * v(i - 1, j);
            end
            if (j > 1)
                r(i, j) = r(i, j) - H_array{k + 1, 1}(2, i, j) * v(i, j - 1);
            end
            if (j < N)
                r(i, j) = r(i, j) - H_array{k + 1, 1}(4, i, j) * v(i, j + 1);
            end
            if (i < M)
                r(i, j) = r(i, j) - H_array{k + 1, 1}(5, i, j) * v(i + 1, j);
            end
            r(i, j) = r(i, j) - H_array{k + 1, 1}(3, i, j) * v(i, j);
        end
    end
    d = reduc(M, N, r);
    for i = 2 : (M + 1) / 2 - 1
        for j = 2 : (M + 1) / 2 - 1
            d(i, j) = d(i, j) * 4;
        end
    end
    e = multigrid((M + 1) / 2, (N + 1) / 2, H_array, d, zeros((M + 1) / 2, (N + 1) / 2), nu, k + 1, K, smoother, reduc, interp);
    v = v + interp((M + 1) / 2, (N + 1) / 2, e);
else
    l = 1;
    r = zeros(M * N, 1);
    for i = 1 : M
        for j = 1 : N
            idx1(l) = (i - 1) * N + j;
            idx2(l) = (i - 1) * N + j;
            val(l) = H_array{k + 1, 1}(3, i, j);
            l = l + 1;
            if (i > 1)
                idx1(l) = (i - 1) * N + j;
                idx2(l) = (i - 2) * N + j;
                val(l) = H_array{k + 1, 1}(1, i, j);
                l = l + 1;
            end
            if (j > 1)
                idx1(l) = (i - 1) * N + j;
                idx2(l) = (i - 1) * N + j - 1;
                val(l) = H_array{k + 1, 1}(2, i, j);
                l = l + 1;
            end
            if (j < N)
                idx1(l) = (i - 1) * N + j;
                idx2(l) = (i - 1) * N + j + 1;
                val(l) = H_array{k + 1, 1}(4, i, j);
                l = l + 1;
            end
            if (i < M)
                idx1(l) = (i - 1) * N + j;
                idx2(l) = i * N + j;
                val(l) = H_array{k + 1, 1}(5, i, j);
                l = l + 1;
            end
            r((i - 1) * N + j) = R(i, j);
        end
    end
    H = sparse(idx1, idx2, val);
    V = H \ r;
    v = zeros(M, N);
    for i = 1 : M
        for j = 1 : N
            v(i, j) = V((i - 1) * N + j);
        end
    end
end