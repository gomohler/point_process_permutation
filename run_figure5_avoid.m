rng(123)

[times x]=generate_avoiding_data(.2,100);

subplot(1,5,1);
plot(x(:,1),times,'.');

kpts=100;
kcut=.3;
[Kf]=K_Fun_Diff(x,times,kpts,kcut,1);
subplot(1,5,2);
plot(Kf);

Npts=max(size(x));

tperm=randsample(times,Npts,'false');
subplot(1,5,3);
plot(x(:,1),tperm,'.');

tperm=randsample(times,Npts,'false');
[tperm1,Kftemp1,errs1]=run_2_opt(x,tperm,Npts,kpts,kcut,Kf,.01,1,6000);
tperm=randsample(times,Npts,'false');
[tperm2,Kftemp2,errs2]=run_2_opt(x,tperm,Npts,kpts,kcut,Kf,.01,1,6000);
tperm=randsample(times,Npts,'false');
[tperm3,Kftemp3,errs3]=run_2_opt(x,tperm,Npts,kpts,kcut,Kf,.01,1,6000);

subplot(1,5,4);
plot(x(:,1),tperm1,'r.',x(:,1),tperm2,'b.',x(:,1),tperm3,'k.');

subplot(1,5,2);
[Kf1]=K_Fun_Diff(x,tperm1,kpts,kcut,1);
[Kf2]=K_Fun_Diff(x,tperm2,kpts,kcut,1);
[Kf3]=K_Fun_Diff(x,tperm3,kpts,kcut,1);
[Kfr]=K_Fun_Diff(x,tperm,kpts,kcut,1);
plot([1:kpts]*kcut/kpts,Kf,'-','LineWidth',1);
hold on
plot([1:kpts]*kcut/kpts,Kf1,'-r','LineWidth',1);
plot([1:kpts]*kcut/kpts,Kf2,'-b','LineWidth',1);
plot([1:kpts]*kcut/kpts,Kf3,'-k','LineWidth',1);
plot([1:kpts]*kcut/kpts,Kfr,'-','LineWidth',1);
hold off

subplot(2,2,1)
plot(x(:,1),times,'k.','MarkerSize',10); axis([0 1 0 1]);
title('Data (Regular)')
xlabel('x') 
ylabel('t') 
subplot(2,2,2)
plot(x(:,1),tperm,'b.','MarkerSize',10);axis([0 1 0 1]);
title('Random permutation')
xlabel('x') 
ylabel('t')
subplot(2,2,3)
plot(x(:,1),tperm1,'.','color',[0.3010 0.7450 0.9330],'MarkerSize',10);axis([0 1 0 1]);
hold on
plot(x(:,1),tperm2,'.','color',[0.4660 0.6740 0.1880],'MarkerSize',10);axis([0 1 0 1]);
plot(x(:,1),tperm3,'.','color',[0.8500 0.3250 0.0980],'MarkerSize',10); axis([0 1 0 1]);
title('Three L-function preserving permutations')
xlabel('x') 
ylabel('t')
hold off
subplot(2,2,4)
plot([1:kpts]*kcut/kpts,Kf,'-k','LineWidth',2);
hold on
plot([1:kpts]*kcut/kpts,Kf1,'-.','LineWidth',2,'color',[0.4660 0.6740 0.1880]);
[Kfr]=K_Fun_Diff(x,tperm,kpts,kcut,1);
plot([1:kpts]*kcut/kpts,Kfr,'-.b','LineWidth',2);
title('L-function')
xlabel('r') 
ylabel('L(r)')
legend('Data','Optimized','Random','Location','northwest')
hold off

print(gcf,'regular.png','-dpng','-r1000')

