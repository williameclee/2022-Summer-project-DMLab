function Yd1Fd4 = NDiffFd4(X,Y)
    X = squeeze(X);
    Y = squeeze(Y);
    if size(X,1) ~= 1
        X = X.';
    end
    if size(Y,1) ~= 1
        Y = Y.';
    end
    dX     = X(2)-X(1);
    Yd1Fd4 = (+8*[Y(2:end),Y(1)]-8*[Y(end),Y(1:end-1)] ...
              -1*[Y(3:end),Y(1:2)]+1*[Y(end-1:end),Y(1:end-2)])/(12*dX);
end