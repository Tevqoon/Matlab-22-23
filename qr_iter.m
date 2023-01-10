function X = qr_iter(X, tol)
    [x,y] = size(X);
    if x == 1
        return
    end
    while abs(X(x, y - 1)) >= (abs(X(x - 1, y - 1)) + abs(X(x,y))) * tol
        s = X(x, y);
        [Q, R] = qr(X - s * eye(x,y));
        X = R * Q + s * eye(x, y);
    end
    X(1:x - 1, 1:y - 1) = qr_iter(X(1:x - 1, 1:y - 1), tol);
end

