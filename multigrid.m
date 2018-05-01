function v = multigrid(M, N, H_array, R, u, nu, k, K, smoother, reduc, interp)
if k < K
    v = smoother(M, N, H_array{k + 1, 1}, R, u, nu);
    r = R;
    for i = 2 : M - 1
        for j = 2 : N - 1
            if (i > 2)
                r(i, j) = r(i, j) - H_array{k + 1, 1}(1, i, j) * v(i - 1, j);
            end
            if (j > 2)
                r(i, j) = r(i, j) - H_array{k + 1, 1}(2, i, j) * v(i, j - 1);
            end
            if (j < N - 1)
                r(i, j) = r(i, j) - H_array{k + 1, 1}(4, i, j) * v(i, j + 1);
            end
            if (i < M - 1)
                r(i, j) = r(i, j) - H_array{k + 1, 1}(5, i, j) * v(i + 1, j);
            end
            r(i, j) = r(i, j) - H_array{k + 1, 1}(3, i, j) * v(i, j);
        end
    end
    d = reduc(M, N, r);
    d = 4 * d;
    e = multigrid((M + 1) / 2, (N + 1) / 2, H_array, d, zeros((M + 1) / 2, (N + 1) / 2), nu, k + 1, K, smoother, reduc, interp);
    v = v + interp((M + 1) / 2, (N + 1) / 2, e);
    v = smoother(M, N, H_array{k + 1, 1}, R, v, nu);
else
    l = 1;
    r = zeros((M - 2) * (N - 2), 1);
    for i = 2 : M - 1
        for j = 2 : N - 1
            idx1(l) = (i - 2) * (N - 2) + j - 1;
            idx2(l) = (i - 2) * (N - 2) + j - 1;
            val(l) = H_array{k + 1, 1}(3, i, j);
            l = l + 1;
            if (i > 2)
                idx1(l) = (i - 2) * (N - 2) + j - 1;
                idx2(l) = (i - 3) * (N - 2) + j - 1;
                val(l) = H_array{k + 1, 1}(1, i, j);
                l = l + 1;
            end
            if (j > 2)
                idx1(l) = (i - 2) * (N - 2) + j - 1;
                idx2(l) = (i - 2) * (N - 2) + j - 2;
                val(l) = H_array{k + 1, 1}(2, i, j);
                l = l + 1;
            end
            if (j < N - 1)
                idx1(l) = (i - 2) * (N - 2) + j - 1;
                idx2(l) = (i - 2) * (N - 2) + j;
                val(l) = H_array{k + 1, 1}(4, i, j);
                l = l + 1;
            end
            if (i < M - 1)
                idx1(l) = (i - 2) * (N - 2) + j - 1;
                idx2(l) = (i - 1) * (N - 2) + j - 1;
                val(l) = H_array{k + 1, 1}(5, i, j);
                l = l + 1;
            end
            r((i - 2) * (N - 2) + j - 1) = R(i, j);
        end
    end
    H = sparse(idx1, idx2, val);
    V = H \ r;
    v = u;
    for i = 2 : M - 1
        for j = 2 : N - 1
            v(i, j) = V((i - 2) * (N - 2) + j - 1);
        end
    end
end