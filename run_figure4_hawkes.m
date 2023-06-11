
rng(123)
mu=40;
T=1;
sig=.01;
w=100;
k0=.75;
%back=50;
back=-1;

[tn x]=generate_hawkes_data(mu,k0,w,T,sig,back);


xn=x(:,1);
yn=x(:,2);

[tf x2]=generate_hawkes_data(mu,k0,w,T,sig,back);

xf=x2(:,1);
yf=x2(:,2);


plot(xf,tf,'.',xn,tn,'r.')

tcut=.1;
dcut=.1;
tol=.08;
maxiter=6000;

[KS]=knox_statistic_euc(tcut,dcut,yn,xn,tn,yf,xf,tf);

Lf=1;
kpts=100;
kcut=.3;
Npts=max(size(xn));



Tmax=max(tn);
xns=(xn-min(xn))/(max(xn)-min(xn));
yns=(yn-min(yn))/(max(yn)-min(yn));
[Kf]=K_Fun_Diff([yns xns],tn/Tmax,kpts,kcut,Lf);



Miter=200;
KSvec=zeros(Miter,1); 
KSvecopt=zeros(Miter,1); 

Kfvec=zeros(Miter,kpts);
Kfvecopt=zeros(Miter,kpts);

time_perm_starts=[];

for i=1:Miter
    i
    
idx = randperm(max(size(tn)));
ttmp=tn(idx);

time_perm_starts{i}=ttmp;

[Kfrand]=K_Fun_Diff([yns xns],ttmp/Tmax,kpts,kcut,Lf);
[KSt]=knox_statistic_euc(tcut,dcut,yn,xn,ttmp,yf,xf,tf);

KSvec(i)=KSt;
Kfvec(i,:)=Kfrand';

subplot(2,2,1);
h1 = histogram(KSvec(1:i),10,'FaceColor','r');
title('Knox test random permutation')
xlabel('\kappa') 
hold on
xline(KS,'--','LineWidth',2);
hold off

subplot(2,2,3);
hold on
plot(Kfvec(1:i,:)','r');
plot(Kf,'k--','LineWidth',2);
hold off
title('L function rand permutation')
xlabel('r') 
ylabel('L(r)') 

drawnow

end

Kfmean=mean(Kfvec,1);

for i=1:Miter
    i
    
    ttmp=time_perm_starts{i};
    
    errK=Kfvec(i,:)-Kfmean;

[tperm1,Kftemp1,errs1]=run_2_opt([yns xns],ttmp/Tmax,Npts,kpts,kcut,Kf+errK',tol,Lf,maxiter);




[KSt2]=knox_statistic_euc(tcut,dcut,yn,xn,tperm1*Tmax,yf,xf,tf);


KSvecopt(i)=KSt2;


Kfvecopt(i,:)=Kftemp1';

save('hawkes_vars.mat','KS','Kf','KSvec','KSvecopt','Kfvec','Kfvecopt');


subplot(2,2,2);
h2 = histogram(KSvecopt(1:i),10,'FaceColor','b');
title('Knox test SOP permutation')
xlabel('\kappa') 
hold on
xline(KS,'--','LineWidth',2);
hold off

subplot(2,2,4);
hold on
plot(Kfvecopt(1:i,:)','b');
plot(Kf,'k--','LineWidth',2);
hold off
title('L function SOP permutation')
xlabel('r') 
ylabel('L(r)') 

drawnow

end


prand=2*(1-sum(KS<KSvec(1:Miter))/Miter)
popt=2*(1-sum(KS<KSvecopt(1:Miter))/Miter)

subplot(2,2,1);
h1 = histogram(KSvec(1:Miter),10,'FaceColor','r');
title('Knox test random permutation')
xlabel('\kappa') 
text(200,40,strcat("p=",num2str(prand)))
hold on
xline(KS,'--','LineWidth',2);
hold off
subplot(2,2,2);
h2 = histogram(KSvecopt(1:Miter),10,'FaceColor','b');
title('Knox test SOP permutation')
text(200,40,strcat("p=",num2str(popt)))
xlabel('\kappa') 
hold on
xline(KS,'--','LineWidth',2);
hold off
subplot(2,2,3);
plot([1:kpts]*kcut/kpts,Kf,'k--','LineWidth',2);
hold on
plot([1:kpts]*kcut/kpts,Kfvec(1:Miter,:)','r');
hold off
title('L function random permutation')
xlabel('r') 
ylabel('L(r)') 
legend('Data','Rand Permutation','Location','southeast')
subplot(2,2,4);
plot([1:kpts]*kcut/kpts,Kf,'k--','LineWidth',2);
hold on
plot([1:kpts]*kcut/kpts,Kfvecopt(1:Miter,:)','b');
plot([1:kpts]*kcut/kpts,Kf,'k--','LineWidth',2);
hold off
title('L function SOP permutation')
xlabel('r') 
ylabel('L(r)') 
legend('Data','SOP Permutation','Location','southeast')


print(gcf,'hawkes_knox_fig.png','-dpng','-r1000')