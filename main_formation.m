close all
%function [xleaderl,xleader2,xfollower1,xfollower2,uleader,ufollower,rleader,rfollower]=reg(tf,dt,alpha)
%function reg(tf,dt,alpha)
tf=27; dt=0.01; alpha=30; alpha_tail=30;
% follower trajectory

speed_up=20;

a=0.2;
T=0.2;
deltat=0.01*T;

L=1.5;

A=[ 0 0 1 0;
    0 0 0 1;
    0 0 -a 0;
    0 0 0 -a];

B= [0 0;
    0 0;
    1 0;
    0 1];

sat=5;

% C= [0 0 1 0;
%     0 0 0 1];

C=[1 0 0 0;
   0 1 0 0];

t=0:dt:tf;

N=size(t,2);
Nbar=(N-1)*23/27;


x_1_1=zeros(4,N);
x_2_1=zeros(4,N);
x_2_2=zeros(4,N);
x_3_1=zeros(4,N);
x_3_2=zeros(4,N);
x_3_3=zeros(4,N);


X0_1_1=[0;0;0;6*speed_up/tf];
X0_2_1=[-L/2;-sqrt(3)*L/2;0;6*speed_up/tf];
X0_2_2=[L/2;-sqrt(3)*L/2;0;6*speed_up/tf];
X0_3_1=[-2*L/2;-2*sqrt(3)*L/2;0;6*speed_up/tf];
X0_3_2=[0;-2*sqrt(3)*L/2;0;6*speed_up/tf];
X0_3_3=[2*L/2;-2*sqrt(3)*L/2;0;6*speed_up/tf];


x_1_1(:,1)=X0_1_1;
x_2_1(:,1)=X0_2_1;
x_2_2(:,1)=X0_2_2;
x_3_1(:,1)=X0_3_1;
x_3_2(:,1)=X0_3_2;
x_3_3(:,1)=X0_3_3;

u_1_1=zeros(2,N);
u_2_1=zeros(2,N);
u_2_2=zeros(2,N);
u_3_1=zeros(2,N);
u_3_2=zeros(2,N);
u_3_3=zeros(2,N);

rvel1_1=zeros(2,N);
r1_1=zeros(2,N);
predict1_1=zeros(2,N);

rvel2_1=zeros(2,N);
rfollower2_1=zeros(2,N);
predictfollower2_1=zeros(2,N);

rvel2_1=zeros(2,N);
rfollower2=zeros(2,N);
predictfollower2=zeros(2,N);

nleader=zeros(1,N);
nfollower1=zeros(1,N);
nfollower2=zeros(1,N);
nfollower21=zeros(1,N);
nfollower22=zeros(1,N);
nfollower23=zeros(1,N);


for i=1:(N-1)/27*10
    rvel1_1(:,i)=[2*pi*sin(2*pi*i/(N*20/27))/tf,2*pi*cos(2*pi*i/(N*20/27))/tf]';
end

for i=(N-1)/27*10+1:(N-1)*20/27
    rvel1_1(:,i)=[-2*pi*sin(2*pi*i/(N*20/27))/tf,2*pi*cos(2*pi*i/(N*20/27))/tf]';
end

for i=(N-1)*20/27+1:N
    rvel1_1(:,i)=[0,2*pi/tf]';
end    

rvel1_1=speed_up*rvel1_1;

for i=2:N
    r1_1(:,i)=r1_1(:,i-1)+rvel1_1(:,i-1)*dt;
end    

for i=1:N-T/dt
   predict1_1(:,i)=r1_1(:,i+T/dt); 
end

for i=N-T/dt+1:N
   predict1_1(:,i)=predict1_1(:,N-T/dt); 
end    

reffollower1(:,1)=C*(x_1_1(:,1)-x_2_1(:,1));
reffollower2(:,1)=C*(x_1_1(:,1)-x_2_2(:,1));

for i=2:N 
    xdot=dxdt(x_1_1(:,i-1),u_1_1(:,i-1),A,B);
    x_1_1(:,i)=x_1_1(:,i-1)+xdot*dt;
    real_vel_11=x_1_1(3:4,i-1);
    [gu(:,i),guprime]=future_calc(x_1_1(:,i-1),u_1_1(:,i-1),A,B,C, T, deltat);
    
    deltau=guprime\(predict1_1(:,i)-gu(:,i))*alpha*dt;
    u_1_1(:,i)=u_1_1(:,i-1)+deltau;
    
     % Begin  saturation u 
   % if norm(uleader(:,i))>1
   %     uleader(:,i)=1*uleader(:,i)/norm(uleader(:,i));
    % end
     %  end saturation
    
     
    
     xdot=dxdt(x_2_1(:,i-1),u_2_1(:,i-1),A,B);
     x_2_1(:,i)=x_2_1(:,i-1)+xdot*dt;     
     real_vel_21=x_2_1(3:4,i-1);
     [gu1(:,i),guprime1]=future_calc(x_2_1(:,i-1),u_2_1(:,i-1),A,B,C, T, deltat);
     
     xdot=dxdt(x_2_2(:,i-1),u_2_2(:,i-1),A,B);
     x_2_2(:,i)=x_2_2(:,i-1)+xdot*dt;     
     [gu2(:,i),guprime2]=future_calc(x_2_2(:,i-1),u_2_2(:,i-1),A,B,C, T, deltat);
     real_vel_22=x_2_2(3:4,i-1);
     
     xdot=dxdt(x_3_1(:,i-1),u_3_1(:,i-1),A,B);
     x_3_1(:,i)=x_3_1(:,i-1)+xdot*dt;     
     [gu21(:,i),guprime21]=future_calc(x_3_1(:,i-1),u_3_1(:,i-1),A,B,C, T, deltat);
     real_vel_31=x_3_1(3:4,i-1);
     
     
     xdot=dxdt(x_3_2(:,i-1),u_3_2(:,i-1),A,B);
     x_3_2(:,i)=x_3_2(:,i-1)+xdot*dt;     
     [gu22(:,i),guprime22]=future_calc(x_3_2(:,i-1),u_3_2(:,i-1),A,B,C, T, deltat);
     real_vel_32=x_3_2(3:4,i-1);
     
     xdot=dxdt(x_3_3(:,i-1),u_3_3(:,i-1),A,B);
     x_3_3(:,i)=x_3_3(:,i-1)+xdot*dt;     
     [gu23(:,i),guprime23]=future_calc(x_3_3(:,i-1),u_3_3(:,i-1),A,B,C, T, deltat);
     real_vel_33=x_3_3(3:4,i-1);
    
    if i<2*T/dt+1
        ki=i-1;
    else
        ki=2*T/dt;
    end
    
    
%     if i<2*T/dt
%         p=xleader(1:2,i)-xleader(1:2,1);
%     else    
%         p=xleader(1:2,i)-xleader(1:2,i-ki+T/dt);
%     end    
        
        p=x_1_1(1:2,i)-x_1_1(1:2,i-1);
        p=real_vel_11*dt;

        theta1=5*pi/6;
        theta2=-5*pi/6;
        rot1=[cos(theta1),-sin(theta1);sin(theta1),cos(theta1)];
        rot2=[cos(theta2),-sin(theta2);sin(theta2),cos(theta2)];
        
        if norm(p)==0
            q1=0;
            q2=0;
        else
            q1=rot1*p/norm(p);
            q2=rot2*p/norm(p);
%               q1=rot1*p;
%               q2=rot2*p;

        end
        
        reffollower1(:,i)=C*x_1_1(:,i)+L*q1+p*(T-dt)/dt;
        pos21(:,i)=C*x_1_1(:,i)+L*q1;
        reffollower2(:,i)=C*x_1_1(:,i)+L*q2+p*(T-dt)/dt;
        pos22(:,i)=C*x_1_1(:,i)+L*q2;
        
        pos31(:,i)=C*x_1_1(:,i)+2*L*q1;
        pos33(:,i)=C*x_1_1(:,i)+2*L*q2;
        pos32(:,i)=C*x_1_1(:,i)+L*q1+L*q2;
    
    deltau1=guprime1\(reffollower1(:,i)-gu1(:,i))*alpha*dt;
    u_2_1(:,i)=u_2_1(:,i-1)+deltau1;
    
    deltau2=guprime2\(reffollower2(:,i)-gu2(:,i))*alpha*dt;
    u_2_2(:,i)=u_2_2(:,i-1)+deltau2;
    
    
    
%%   for 3rd layer

%% follower 21 and 22_1
    if i<2*T/dt+1
        ki=i-1;
    else
        ki=2*T/dt;
    end
    
    
%     if i<2*T/dt
%         p=xfollower1(1:2,i)-xfollower1(1:2,1);
%     else    
%         p=xfollower1(1:2,i)-xfollower1(1:2,i-ki+T/dt);
%     end    
        
        p=x_2_1(1:2,i)-x_2_1(1:2,i-1);
        p=real_vel_21*dt;
        
        theta1=5*pi/6;
        theta2=-5*pi/6;
        rot1=[cos(theta1),-sin(theta1);sin(theta1),cos(theta1)];
        rot2=[cos(theta2),-sin(theta2);sin(theta2),cos(theta2)];
        
        if norm(p)==0
            q1=0;
            q2=0;
        else
            q1=rot1*p/norm(p);
            q2=rot2*p/norm(p);
        end
        
        reffollower21(:,i)=C*x_2_1(:,i)+L*q1+p*(T-dt)/dt;
        reffollower22_1(:,i)=C*x_2_1(:,i)+L*q2+p*(T-dt)/dt;
  %      pos31(:,i)=C*xfollower1(:,i)+L*q1;
  %      pos32_1(:,i)=C*xfollower1(:,i)+L*q2;
        
%% follower 22_2 and 23

    if i<2*T/dt+1
        ki=i-1;
    else
        ki=2*T/dt;
    end
    
    
%     if i<2*T/dt
%         p=xfollower2(1:2,i)-xfollower2(1:2,1);
%     else    
%         p=xfollower2(1:2,i)-xfollower2(1:2,i-ki+T/dt);
%     end    

        
        p=x_2_2(1:2,i)-x_2_2(1:2,i-1);
        p=real_vel_22*dt;
        
        theta1=5*pi/6;
        theta2=-5*pi/6;
        rot1=[cos(theta1),-sin(theta1);sin(theta1),cos(theta1)];
        rot2=[cos(theta2),-sin(theta2);sin(theta2),cos(theta2)];
        
        if norm(p)==0
            q1=0;
            q2=0;
        else
            q1=rot1*p/norm(p);
            q2=rot2*p/norm(p);
        end
        
        reffollower22_2(:,i)=C*x_2_2(:,i)+L*q1+p*(T-dt)/dt;
        reffollower23(:,i)=C*x_2_2(:,i)+L*q2+p*(T-dt)/dt;       
        
   %     pos32_2(:,i)=C*xfollower2(:,i)+L*q1;
    %    pos33(:,i)=C*xfollower2(:,i)+L*q2;
        
    %    pos32(:,i)=0.5*(pos32_2(:,i)+pos32_1(:,i));
        
    reffollower22(:,i)=(reffollower22_1(:,i)+reffollower22_2(:,i))/2;    
        
    deltau21=guprime21\(reffollower21(:,i)-gu21(:,i))*alpha_tail*dt;
    u_3_1(:,i)=u_3_1(:,i-1)+deltau21;
    
    deltau22=guprime22\(reffollower22(:,i)-gu22(:,i))*alpha_tail*dt;
    u_3_2(:,i)=u_3_2(:,i-1)+deltau22;
    
    deltau23=guprime23\(reffollower23(:,i)-gu23(:,i))*alpha_tail*dt;
    u_3_3(:,i)=u_3_3(:,i-1)+deltau23;
        
     % Begin  saturation u 
   % if norm(ufollower1(:,i))>3
   %     ufollower1(:,i)=3*ufollower1(:,i)/norm(ufollower1(:,i));
    % end
     %  end saturation
   
   
     % Begin  saturation u 
   % if norm(ufollower(:,i))>3
   %     ufollower(:,i)=3*ufollower(:,i)/norm(ufollower(:,i));
    % end
     %  end saturation
     
     %% add saturation
     
     if norm(u_1_1(:,i))>sat
       u_1_1(:,i)=sat*u_1_1(:,i)/norm(u_1_1(:,i));
     end
     
     if norm(u_2_1(:,i))>sat
       u_2_1(:,i)=sat*u_2_1(:,i)/norm(u_2_1(:,i));
     end
  
     if norm(u_2_2(:,i))>sat
       u_2_2(:,i)=sat*u_2_2(:,i)/norm(u_2_2(:,i));
     end
     
     if norm(u_3_1(:,i))>sat
       u_3_1(:,i)=sat*u_3_1(:,i)/norm(u_3_1(:,i));
     end
     
     if norm(u_3_2(:,i))>sat
       u_3_2(:,i)=sat*u_3_2(:,i)/norm(u_3_2(:,i));
     end
     
     if norm(u_3_3(:,i))>sat
       u_3_3(:,i)=sat*u_3_3(:,i)/norm(u_3_3(:,i));
     end
     
   
     nleader(i)=norm(u_1_1(:,i));
     nfollower1(i)=norm(u_2_1(:,i));
     nfollower2(i)=norm(u_2_2(:,i));
     nfollower21(i)=norm(u_3_1(:,i));
     nfollower22(i)=norm(u_3_2(:,i));
     nfollower23(i)=norm(u_3_3(:,i));
end



%% input norm graph
    figure(4)

    plot(t(1:Nbar),nleader(1:Nbar), 'LineWidth',1.5);
    hold on
    plot(t(1:Nbar),nfollower1(1:Nbar), 'LineWidth',1.5);
    plot(t(1:Nbar),nfollower1(1:Nbar), 'LineWidth',1.5);
    plot(t(1:Nbar),nfollower21(1:Nbar), 'LineWidth',1.5);
    plot(t(1:Nbar),nfollower22(1:Nbar), 'LineWidth',1.5);
    plot(t(1:Nbar),nfollower23(1:Nbar), 'LineWidth',1.5);
    
    x1=xlabel('Time$~[s]$');
 y1=ylabel('Input norm$~[N]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
leg1=legend('$A_{1,1}$','$A_{2,1}$','$A_{2,2}$','$A_{3,1}$','$A_{3,2}$','$A_{3,3}$');
 set(leg1,'Interpreter','latex')
 
%  set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('Input_norm','-dsvg','-r0')


%% r-g

for i=1:N
r_g_leader(i)=norm(predict1_1(:,i)-gu(:,i));
r_g_follower1(i)=norm(reffollower1(:,i)-gu1(:,i));
r_g_follower2(i)=norm(reffollower2(:,i)-gu2(:,i));
r_g_follower21(i)=norm(reffollower21(:,i)-gu21(:,i));
r_g_follower22(i)=norm(reffollower22(:,i)-gu22(:,i));
r_g_follower23(i)=norm(reffollower23(:,i)-gu23(:,i));
end

figure(42)

    plot(t(2:Nbar),r_g_leader(2:Nbar), 'LineWidth',1.5);
    hold on
    plot(t(2:Nbar),r_g_follower1(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),r_g_follower2(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),r_g_follower21(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),r_g_follower22(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),r_g_follower23(2:Nbar), 'LineWidth',1.5);


    
        x1=xlabel('Time$~[s]$');
 y1=ylabel('Control Error$~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
leg1=legend('$A_{1,1}$','$A_{2,1}$','$A_{2,2}$','$A_{3,1}$','$A_{3,2}$','$A_{3,3}$');
 set(leg1,'Interpreter','latex')
 
%  set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('r_g','-dsvg','-r0')

%% total error

%{

for i=1:N
error_1_1(i)=norm(r1_1(:,i)-C*x_1_1(:,i));
error_2_1(i)=norm(pos21(:,i)-C*x_2_1(:,i));
error_2_2(i)=norm(pos22(:,i)-C*x_2_2(:,i));
error_3_1(i)=norm(pos31(:,i)-C*x_3_1(:,i));
error_3_2(i)=norm(pos32(:,i)-C*x_3_2(:,i));
error_3_3(i)=norm(pos33(:,i)-C*x_3_3(:,i));
end

figure(21)

    plot(t(2:Nbar),error_1_1(2:Nbar), 'LineWidth',1.5);
    hold on
    plot(t(2:Nbar),error_2_1(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),error_2_2(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),error_3_1(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),error_3_2(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),error_3_3(2:Nbar), 'LineWidth',1.5);


    
        x1=xlabel('Time$~[s]$');
 y1=ylabel('Total Error$~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
leg1=legend('$A_{1,1}$','$A_{2,1}$','$A_{2,2}$','$A_{3,1}$','$A_{3,2}$','$A_{3,3}$');
 set(leg1,'Interpreter','latex')
 
%  set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
title('Total Error vs time')
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('total_error','-dsvg','-r0')
%}

%% total error new

for i=1:N-T/dt
error_1_1(i)=norm(predict1_1(:,i)-C*x_1_1(:,i+T/dt));
error_2_1(i)=norm(reffollower1(:,i)-C*x_2_1(:,i+T/dt));
error_2_2(i)=norm(reffollower2(:,i)-C*x_2_2(:,i+T/dt));
error_3_1(i)=norm(reffollower21(:,i)-C*x_3_1(:,i+T/dt));
error_3_2(i)=norm(reffollower22(:,i)-C*x_3_2(:,i+T/dt));
error_3_3(i)=norm(reffollower23(:,i)-C*x_3_3(:,i+T/dt));
end

figure(21)

    plot(t(2:Nbar),error_1_1(2:Nbar), 'LineWidth',1.5);
    hold on
    plot(t(2:Nbar),error_2_1(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),error_2_2(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),error_3_1(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),error_3_2(2:Nbar), 'LineWidth',1.5);
    plot(t(2:Nbar),error_3_3(2:Nbar), 'LineWidth',1.5);


    
        x1=xlabel('Time$~[s]$');
 y1=ylabel('Total Error$~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
leg1=legend('$A_{1,1}$','$A_{2,1}$','$A_{2,2}$','$A_{3,1}$','$A_{3,2}$','$A_{3,3}$');
 set(leg1,'Interpreter','latex')
 
%  set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('total_error','-dsvg','-r0')



%% inter_agent distance

figure(97)
for i=1:N
   dL_11(i)= norm(x_1_1(1:2,i)-x_2_1(1:2,i));
   dL_12(i)= norm(x_1_1(1:2,i)-x_2_2(1:2,i));
   
   d11_21(i)= norm(x_2_1(1:2,i)-x_3_1(1:2,i));
   d11_22(i)= norm(x_2_1(1:2,i)-x_3_2(1:2,i));
   d12_22(i)= norm(x_2_2(1:2,i)-x_3_2(1:2,i));
   d12_23(i)= norm(x_2_2(1:2,i)-x_3_3(1:2,i));
end    

plot (t(1:Nbar),dL_11(1:Nbar), 'LineWidth',1.5);
hold on
plot (t(1:Nbar),dL_12(1:Nbar), 'LineWidth',1.5);
plot (t(1:Nbar),d11_21(1:Nbar), 'LineWidth',1.5);
plot (t(1:Nbar),d11_22(1:Nbar), 'LineWidth',1.5);
plot (t(1:Nbar),d12_22(1:Nbar), 'LineWidth',1.5);
plot (t(1:Nbar),d12_23(1:Nbar), 'LineWidth',1.5);

       x1=xlabel('Time$~[s]$');
 y1=ylabel('Inter-agent distance$~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
leg1=legend('$d_{A_{1,1}-A_{2,1}}$','$d_{A_{1,1}-A_{2,2}}$','$d_{A_{2,1}-A_{3,1}}$','$d_{A_{2,1}-A_{3,2}}$','$d_{A_{2,2}-A_{3,2}}$','$d_{A_{2,2}-A_{3,3}}$');
 set(leg1,'Interpreter','latex')
 
%  set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
hold off

ylim([1.4 1.7])

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('d_ij','-dsvg','-r0')

%% lateral error
%{
figure (77)
for i=1:100
    r1(:,i)=[predict1_1(1,1)-1+i/100;predict1_1(2,1)];
    r2(:,i)=[reffollower1(1,1)-1+i/100;reffollower1(2,1)];
    r3(:,i)=[reffollower2(1,1)-1+i/100;reffollower2(2,1)];
    r4(:,i)=[reffollower21(1,1)-1+i/100;reffollower21(2,1)];
    r5(:,i)=[reffollower22(1,1)-1+i/100;reffollower22(2,1)];
    r6(:,i)=[reffollower23(1,1)-1+i/100;reffollower23(2,1)];
    
end    
    
for i=1:N
    r1(:,i+100)=predict1_1(:,i);  
    r2(:,i+100)=reffollower1(:,i);
    r3(:,i+100)=reffollower2(:,i);
    r4(:,i+100)=reffollower21(:,i);
    r5(:,i+100)=reffollower22(:,i);
    r6(:,i+100)=reffollower23(:,i);
end



    
    [close_point,lateral_error_norm(1,:),arc]=distance2curve(r1',x_1_1(1:2,:)');
    [close_point,lateral_error_norm(2,:),arc]=distance2curve(r2',x_2_1(1:2,:)');
    [close_point,lateral_error_norm(3,:),arc]=distance2curve(r3',x_2_2(1:2,:)');
    [close_point,lateral_error_norm(4,:),arc]=distance2curve(r4',x_3_1(1:2,:)');
    [close_point,lateral_error_norm(5,:),arc]=distance2curve(r5',x_3_2(1:2,:)');
    [close_point,lateral_error_norm(6,:),arc]=distance2curve(r6',x_3_3(1:2,:)');
    
    

 plot(t(1:Nbar),lateral_error_norm(1,1:Nbar),'LineWidth',1.5)
 hold on 
 plot(t(1:Nbar),lateral_error_norm(2,1:Nbar),'LineWidth',1.5)
 plot(t(1:Nbar),lateral_error_norm(3,1:Nbar),'LineWidth',1.5)
 plot(t(1:Nbar),lateral_error_norm(4,1:Nbar),'LineWidth',1.5)
 plot(t(1:Nbar),lateral_error_norm(5,1:Nbar),'LineWidth',1.5)
 plot(t(1:Nbar),lateral_error_norm(6,1:Nbar),'LineWidth',1.5)

 
  x1=xlabel('Time$~[s]$');
 y1=ylabel('Lateral Error$~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
 
leg1=legend('$A_{1,1}$','$A_{2,1}$','$A_{2,2}$','$A_{3,1}$','$A_{3,2}$','$A_{3,3}$');
 set(leg1,'Interpreter','latex')
 
%  set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
hold off

pbaspect([2.5 1 1])
fig.PaperUnits = 'inches';
print('lat_error','-dsvg','-r0')

%}

 %% trajectories
    figure(13)
        
    plot(x_1_1(1,1:Nbar),x_1_1(2,1:Nbar),'LineWidth',1.5)
    hold on
    plot(x_2_1(1,1:Nbar),x_2_1(2,1:Nbar),'LineWidth',1.5)
    plot(x_2_2(1,1:Nbar),x_2_2(2,1:Nbar),'LineWidth',1.5)
    plot(x_3_1(1,1:Nbar),x_3_1(2,1:Nbar),'LineWidth',1.5)
    plot(x_3_2(1,1:Nbar),x_3_2(2,1:Nbar),'LineWidth',1.5)
    plot(x_3_3(1,1:Nbar),x_3_3(2,1:Nbar),'LineWidth',1.5)
    pbaspect([2.5 1 1])
    
    
    xlim([-5 95]);
    ylim([-20 20]);
    
    
        x1=xlabel('$Z_1~[m]$');
 y1=ylabel('$Z_2~[m]$');
  set(x1,'Interpreter','latex')
 set(y1,'Interpreter','latex')
%     
s1 = plot(x_1_1(1,1),x_1_1(2,1),'o','MarkerSize', 4,'MarkerFaceColor','red');      % bot 1
%q = plot(r(1,1),r(2,1),'o','MarkerFaceColor','blue');                          % ref
s2 = plot(x_2_1(1,1),x_2_1(2,1),'o','MarkerSize', 4,'MarkerFaceColor','green');
s3 = plot(x_2_2(1,1),x_2_2(2,1),'o','MarkerSize', 4,'MarkerFaceColor','green');
s4 = plot(x_3_1(1,1),x_3_1(2,1),'o','MarkerSize', 4,'MarkerFaceColor','blue');
s5 = plot(x_3_2(1,1),x_3_2(2,1),'o','MarkerSize', 4,'MarkerFaceColor','blue');
s6 = plot(x_3_3(1,1),x_3_3(2,1),'o','MarkerSize', 4,'MarkerFaceColor','blue');
 
for k = 2:Nbar
    s1.XData = x_1_1(1,k);
    s1.YData = x_1_1(2,k);
      
%     q.XData = r(1,k);
%     q.YData = r(2,k);
    
    s2.XData = x_2_1(1,k);
    s2.YData = x_2_1(2,k);
    
    s3.XData = x_2_2(1,k);
    s3.YData = x_2_2(2,k);

    s4.XData = x_3_1(1,k);
    s4.YData = x_3_1(2,k);
    
    
    s5.XData = x_3_2(1,k);
    s5.YData = x_3_2(2,k);
    
    
    s6.XData = x_3_3(1,k);
    s6.YData = x_3_3(2,k);
    
    
    drawnow
end


 
leg1=legend('$A_{1,1}$','$A_{2,1}$','$A_{2,2}$','$A_{3,1}$','$A_{3,2}$','$A_{3,3}$');
 set(leg1,'Interpreter','latex')
  
%  set(gcf, 'color', 'none');
% set(gca, 'color', 'none');
hold off


fig.PaperUnits = 'inches';
print('Path','-dsvg','-r0')










