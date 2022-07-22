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
    Yd1Fd6 = (+45*circshift(Y,-1)-45*circshift(Y,1) ...
               -9*circshift(Y,-2) +9*circshift(Y,2) ...
               +1*circshift(Y,-3) -1*circshift(Y,3))/(60*dX);
end