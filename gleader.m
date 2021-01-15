function  [gu,guprime]=future_calc(x0,u,A,B, T, dt)

C= [0 0 1 0;
    0 0 0 1];

M=T/deltat;

x=zeros(4,M);

x(:,1)=x0;

w=zeros(4,2,M);

for i=2:M
    x(:,i)=x(:,i-1)+(A*x(:,i-1)+B*u)*dt;    
    w(:,:,i)=w(:,:,i-1)+(A*w(:,:,i-1)+B)*dt;

end

gu=C*x(:,M);
guprime=C*w(:,:,M);

end





    


