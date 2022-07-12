function Yd1Fft = NDiffFft(X,Y,od)
    if ~exist('order','var')
         od = 1;
    end
    X = squeeze(X);
    Y = squeeze(Y);
    if size(X,1) ~= 1
        X = X.';
    end
    if size(Y,1) ~= 1
        Y = Y.';
    end
    npts = length(X);
    L    = X(end)-X(1)+X(2)-X(1);
    if mod(length(Y),2) == 0
        omega = fftshift(1i*(2*pi/L)*(-(npts/2):(npts/2)-1));
    else
        omega = fftshift(1i*(2*pi/L)*(-floor(npts/2):floor(npts/2)));
    end
    Yhat   = fft(Y);
    Yd1hat = omega.^od.*Yhat;
    Yd1Fft = real(ifft(Yd1hat));
end