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
    Yd1Fd4 = (+8*circshift(Y,-1)-8*circshift(Y,1) ...
              -1*circshift(Y,-2)+1*circshift(Y,2))/(12*dX);
end