%% Setup default parameters
Omg     = 0.1;
L       = 10;
N       = 2^7;
X       = linspace(-L,L,N+1);
dX      = X(2)-X(1);
X       = X(1:end-1);
Y       = X.';
[Xg,Yg] = meshgrid(X,Y);
U       = -Omg*Yg;
V       = Omg*Xg;

%% Define initial field
dT = 50*Omg/N;
% dT = 0.05;
T  = 0:dT:2*(pi/Omg);
% Q  = zeros(length(T),Npts,Npts);
Q0 = zeros(N);
s  = 0.6;
for i = 1:N
    for j = 1:N
        r = sqrt((X(i)-L/3)^2+(Y(j))^2);
        if r/s < pi
            Q0(j,i) = (cos(r/s)+1);
        end
    end
end
% Q(1,:,:) = Q0;
Q = Q0;
% Q(1,:,:) = exp(-5*(sqrt((mod(Xg-L/2+L,2*L)-L).^2+(mod(Yg-L/2+L,2*L)-L).^2)));

%% Intergrating
for t = 1:length(T)-1
%     Q(t+1,:,:) = NIntgRk4(T(t),squeeze(Q(t,:,:)),dT,@adv,U,V,X,Y);
    Q = NIntgRk4(T(t),Q,dT,@adv,U,V,X,Y);
end

%% Plot
circulation_plot;

%% Function
function Qt = adv(~,Q,U,V,X,Y)
    Nx = length(X);
    Ny = length(Y);
    Qx = zeros(size(Q));
    Qy = zeros(size(Q));
    for i = 1:Ny
        Qx(i,:) = NDiffFft(X,Q(i,:));
    end
    for j = 1:Nx
        Qy(:,j) = NDiffFft(Y,Q(:,j)).'; 
    end
    Qt = -(U.*Qx+V.*Qy);
end