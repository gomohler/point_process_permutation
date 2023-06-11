function [times x]=generate_hawkes_data(mu,k0,w,T,sig,back)
%this program simulates event times, called "times", according
%to a self-exciting point process with exponential triggering kernel and
%parameters mu (const. background rate)
%k0 (branching ratio) and w (exp parameter) on the time interval [0,T]
if nargin < 6
  back=-1;
end


times=zeros(5000,1);
x=zeros(5000,2);

%first simulate "background" events
%this is done by picking p points where p is Poisson with parameter mu*T
%and then distributing the points uniformly in the interval [0,T]
if(back<0)
p=pois(mu*T);
else
p=back;
end
times(1:p,1)=rand(p,1)*T;
x(1:p,1)=rand(p,1);
x(1:p,2)=rand(p,1);


counts=1;
countf=p;

%Next loop through every event and simulate the "offspring"
%even the offspring events can generate their own offspring

while((countf-counts)>-1)
p=pois(k0); %each event generates p offspring according to a Poisson r.v. with parameter k0
for j=1:p
    temp=times(counts)-log(rand())/w; % this generates an exponential r.v. on [t_counts,infty]
    temp2=x(counts,1)+sig*randn(); % inter-point distances are gaussian
    temp3=x(counts,2)+sig*randn();
    if(temp<T)                        % we only keep this time if it is in [t_counts,T]
        countf=countf+1;
        times(countf)=temp;
        x(countf,1)=temp2;
        x(countf,2)=temp3;
     
    else
    end
end
counts=counts+1;
end
data=[times(1:countf) x(1:countf,1) x(1:countf,2)];
data=sortrows(data,1);
times=data(:,1);
x=data(:,2:3);





end


function p=pois(S)

if(S<=100)
temp=-S;
L=exp(temp);
k=0;
p=1; 
while(p > L)
k=k+1;
p=p*rand();
end
p=k-1;
else
p=floor(S+S^.5*randn());
end
end