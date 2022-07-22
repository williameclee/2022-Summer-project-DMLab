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
T  = 0:dT:2*(2*pi/Omg);
Q0   = CalQAna(X,Y,0);
QNum = Q0;
QFd4 = Q0;
E = zeros(size(T));
Ea = zeros(size(T));
Ef = zeros(size(T));
M = VideoWriter('advection_fftfd1.mp4','MPEG-4');
M.Quality    = 100;
M.FrameRate  = 30;
open(M)

customcolour
customcolourmap
% figure(1)
% clf;
% subplot(1,2,1)
% contourf(X,Y,QFft,-2:0.02:2,'linecolor','none'); 
% colorbar; caxis([0 1]); axis equal
% xlim([-L,L]); ylim([-L,L]);
% subplot(1,2,2)
% contourf(X,Y,QFft-CalQAna(X,Y,0),-2:0.02:2,'linecolor','none'); 
% colorbar; caxis([0 1]); axis equal
% xlim([-L,L]); ylim([-L,L]);
% writeVideo(M,getframe(gcf));

%% Intergrating
for t = 2:length(T)
    QAna = CalQAna(X,Y,Omg*T(t));
    QNum = NIntgFd1(T(t-1),QNum,dT,@Adv,U,V,X,Y);
    E(t) = rms(QNum-QAna,"all")^2;
    sn  = std(QNum,[],"all");
    sa  = std(QAna,[],"all");
    rho = corrcoef(QNum,QAna);
    rho = rho(2);
    Ea(t) = (mean(QNum,"all")-mean(QAna,"all"))^2+(sn-sa)^2;
    Ef(t) = 2*(1-rho)*sn*sa;
    if mod(t,20) == 1
        EX_advection_2dn_plot;
        writeVideo(M,getframe(gcf));
    end
end

EX_advection_2dn_plot;
writeVideo(M,getframe(gcf));
M.close;

%% Function
function Qt = Adv(~,Q,U,V,X,Y)
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

function QAna = CalQAna(X,Y,phi)
    QAna = zeros(length(Y),length(X));
    R    = 10/3;
    Rx   = R*cos(phi);
    Ry   = R*sin(phi);
    s  = 0.6;
    for i = 1:length(X)
        for j = 1:length(Y)
            r = sqrt((X(i)-Rx)^2+(Y(j)-Ry)^2);
            if r/s < pi
                QAna(j,i) = (cos(r/s)+1)/2;
            end
        end
    end
end