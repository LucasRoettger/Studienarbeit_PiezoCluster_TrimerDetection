function v = getGaussianLike(val, left, right)
    step = 1000;
    x = linspace(1, left+right, step);
    [s, m] = std(x);
    y = normpdf(x, m, s);
    y = y/max(y);
    [~, idx] = min(abs(x-val));

    v = y(idx);

end