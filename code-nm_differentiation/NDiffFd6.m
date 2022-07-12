function Yd1Fd6 = NDiffFd6(X,Y)
    X = squeeze(X);
    Y = squeeze(Y);
    if size(X,1) ~= 1
        X = X.';
    end
    if size(Y,1) ~= 1
        Y = Y.';
    end
    dX     = X(2)-X(1);
    Yd1Fd6 = (+45*[Y(2:end),Y(1:1)]-45*[Y(end-0:end),Y(1:end-1)] ...
              -9*[Y(3:end),Y(1:2)]+9*[Y(end-1:end),Y(1:end-2)] ...
              +1*[Y(4:end),Y(1:3)]-1*[Y(end-2:end),Y(1:end-3)])/(60*dX);
end