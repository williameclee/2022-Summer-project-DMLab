figure(1); clf
subplot(1,2,1)
% contourf(X,Y,squeeze(Q(end,:,:)),-2:0.02:2,'linecolor','none'); 
contourf(X,Y,Q,-2:0.02:2,'linecolor','none'); 
colorbar; caxis([0 1]); axis equal
xlim([-L,L]); ylim([-L,L]);
subplot(1,2,2)
% contourf(X,Y,squeeze(Q(1,:,:)),-2:0.02:2,'linecolor','none'); 
contourf(X,Y,Q-Q0,-2:0.02:2,'linecolor','none'); 
colorbar; caxis([-1 1]); axis equal
xlim([-L,L]); ylim([-L,L]);
% subplot(2,2,3)
% contourf(X,Y,squeeze(Qx),-2:0.02:2,'linecolor','none'); 
% colorbar; caxis([-2 2]); axis equal
% xlim([-L,L]); ylim([-L,L]);
% subplot(2,2,4)
% contourf(X,Y,squeeze(Qy),-2:0.02:2,'linecolor','none'); 
% colorbar; caxis([-2 2]); axis equal
% % contourf(X,Y,adv([],[],U,V,Qx,Qy),'linecolor','none'); colorbar; caxis([-2 2]); axis equal
% xlim([-L,L]); ylim([-L,L]);