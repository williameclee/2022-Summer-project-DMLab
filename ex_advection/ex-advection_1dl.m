L = 10;
Npts = 400;
Tlen = 100;
X = linspace(-L,L,Npts+1);
X = X(1:end-1);
dX = X(2)-X(1);
dT = 0.1*dX;
T = 0:dT:Tlen;
c = 0.2;

U      = zeros(length(T),length(X));
U(1,:) = exp(-5*X.^2);
% UFft      = zeros(length(T),length(X));
% UFft(1,:) = exp(-5*X.^2);
% UFd6      = zeros(length(T),length(X));
% UFd6(1,:) = exp(-5*X.^2);

for i = 1:length(T)-1
    U(i+1,:) = NIntgRk4(0,U(i,:),dT,@fFft,X,c);
%    UFft(i+1,:) = NIntgRk4(0,UFft(i,:),dT,@fFft,X,c);
%    UFd6(i+1,:) = NIntgRk4(0,UFd6(i,:),dT,@fFd6,X,c);
end

clf;
contourf(X,T,U,'linestyle','none'); colorbar
% plot(X,U(2,:))
% hold on
% plot(T,max(UFft,[],2),'DisplayName','FFT')
% plot(T,max(UFd6,[],2),'DisplayName','FD6')
% hold off
% legend('show','box','off')

function Y = fFft(~,U,X,c,~,~,~,~)
    if ~exist('c','var')
         c = 1;
    end
    Y = -c*NDiffFft(X,U);
end
% function Y = fFd6(~,U,X,c,~,~,~,~)
%     if ~exist('c','var')
%          c = 1;
%     end
%     Y = -c*NDiffFd6(X,U);
% end
