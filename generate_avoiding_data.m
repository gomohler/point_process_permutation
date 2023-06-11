function [times x]=generate_avoiding_data(r,N)
x=zeros(N,2);
times=zeros(N,1);
x(1,1)=rand();
x(1,2)=rand();
times(1)=rand();
for i=2:N
    
    check=0;
    while check==0
    tp=rand();
    xp1=rand();
    xp2=rand();
    dist=((tp-times(1:i)).^2+(xp1-x(1:i,1)).^2+(xp2-x(1:i,2)).^2).^.5;
    if min(dist>r)
        times(i)=tp;
        x(i,1)=xp1;
        x(i,2)=xp2;
        check=1;
    end
    end
    
end


end

