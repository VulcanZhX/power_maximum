function H = approx_hessian_fd(gradfun, x, epsilon)
    % 近似计算 Hessian 矩阵
    % gradfun: 计算梯度的函数句柄，输入 x，输出梯度向量
    % x: 计算点（列向量）
    % epsilon: 差分步长，默认 1e-6

    if nargin < 3
        epsilon = 1e-3;
    end

    n = length(x);
    g0 = gradfun(x);
    H = zeros(n, n);

    for i = 1:n
        x_eps = x;
        x_eps(i) = x_eps(i) + epsilon;
        g1 = gradfun(x_eps);
        H(:, i) = (g1 - g0) / epsilon;
    end

    % 保证 Hessian 对称
    H = (H + H') / 2;
end
