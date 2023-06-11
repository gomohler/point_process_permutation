function [KS]=knox_statistic_euc(tcut,dcut,lat1,long1,t1,lat2,long2,t2)

KS=0;
N=max(size(lat1));
M=max(size(lat2));
for i=1:N
    
    for j=1:M
    dt=t2(j)-t1(i);
    
    if(abs(dt)<=tcut)
    
        
   
    dx=((lat1(i)-lat2(j))^2+(long1(i)-long2(j))^2)^.5;
    
    
    
    if(dx<=dcut)
        KS=KS+1;
    end
    
    end
    end
end

end