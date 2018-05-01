function z = f(x, y)
z = -x .* cos(x) .* y .^ 2 ./ sqrt(16 - y .^ 2) + (sin(x) + 3 .* x .* cos(x)) .* sqrt(16 - y .^ 2);