function Yd1Fd2 = NDiffFd2(X,Y)
    X = squeeze(X);
    Y = squeeze(Y);
    if size(X,1) ~= 1
        X = X.';
    end
    if size(Y,1) ~= 1
        Y = Y.';
    end
    dX     = X(2)-X(1);
    Yd1Fd2 = (circshift(Y,-1)-circshift(Y,1))/(2*dX);
end