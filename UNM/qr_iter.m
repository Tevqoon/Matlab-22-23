function X = qr_iter(X, tol, N)
    [x,y] = size(X);
    if x == 1
        return
    end
    k = 0;
    while k < N && abs(X(x, y - 1)) >= (abs(X(x - 1, y - 1)) + abs(X(x,y))) * tol
        s = X(x, y);
        [Q, R] = qr(X - s * eye(x,y));
        X = R * Q + s * eye(x, y);
        k = k + 1;
    end
    X(1:x - 1, 1:y - 1) = qr_iter(X(1:x - 1, 1:y - 1), tol, N);
end

