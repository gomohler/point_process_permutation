function [K]=K_Fun_Diff(x,t,npnts,ctoff,L)
if nargin < 5
  L = 0;
end

N=max(size(x));
K=zeros(npnts,1);
v=[ctoff/npnts:ctoff/npnts:ctoff];
z=[x(:,1)';x(:,2)'; t';];
IP = z' * z;
dmat = sqrt(diag(IP)+diag(IP)' - 2 * IP);
for i=1:npnts
   K(i)=sum(sum(logit(100*(v(i)-dmat)))); 
end

% for j=1:N
%   l=1; %for l=1:N
%        %if(j~=l)
%           dist=((x(j)-x(l))^2+(t(j)-t(l))^2)^.5;
%           for i=1:npnts
%               %if(dist<v(i))
%               %K(i)=K(i)+1;
%               %end
%               K(i)=K(i)+sigmoid(100*(v(i)-dist));
%           end
%        %end
%    %end
% end
K=K/N^2;
if(L==1)
    K=K.^.5;
end
%dKdt=dlgradient(K,t);
end


function y=logit(x)
y=1./(1+exp(-x));
end

             
            
