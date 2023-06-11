function [tperm,Kfnew,errs]=run_2_opt(x,tperm,Npts,kpts,kcut,Kf,tol,L,maxiter)
if nargin < 8
  L=0;
end
if nargin <9
     maxiter=10000;
end

[Kfnew]=K_Fun_Diff(x,tperm,kpts,kcut,L);
errs=zeros(2,1);
errs(1)=norm(Kfnew-Kf);

iters=1;



while errs(iters)>tol && iters<maxiter
    iters=iters+1;
    ttemp=tperm;
    [u]=randsample(Npts,2,'false');
    t1=tperm(u(1));
    t2=tperm(u(2));
    ttemp(u(1))=t2;
    ttemp(u(2))=t1;
    [Kftemp]=K_Fun_Diff(x,ttemp,kpts,kcut,L);
    errs(iters)=norm(Kftemp-Kf);
  
    if(errs(iters)<errs(iters-1))
        tperm=ttemp;
        Kfnew=Kftemp;
        
        
%        subplot(1,5,4);
%        plot(x(:,1),tperm,'.');
%        subplot(1,5,5);
%        plot(errs(1:iters));
%        subplot(1,5,2);
%       plot([1:kpts],Kf,[1:kpts],Kftemp,'or');
%        drawnow;


    else
        errs(iters)=errs(iters-1);
    end
    
end