%% 连杆刚体1-7 head的全部信息 包括坐标和转动角
all_data=xlsread('straight_line_wire.xlsx');
framenum=1662;
link1=all_data(1:framenum,1:8);
link2=all_data(1*framenum+1:2*framenum,1:8);
link3=all_data(2*framenum+1:3*framenum,1:8);
link4=all_data(3*framenum+1:4*framenum,1:8);
link5=all_data(4*framenum+1:5*framenum,1:8);
link6=all_data(5*framenum+1:6*framenum,1:8);
link7=all_data(6*framenum+1:7*framenum,1:8);
head=all_data(7*framenum+1:8*framenum,1:8);
link1=link1(3:end,:);link2=link2(3:end,:);link3=link3(3:end,:);link4=link4(3:end,:);
link5=link5(3:end,:);link6=link6(3:end,:);link7=link7(3:end,:);head=head(3:end,:);



%% 数据连续化
for i=1:size(link4,1)
    for j=1:size(link4,2)
    if (link1(i,j)>5000)||(link1(i,j)==0)
        link1(i,j)=link1(i-1,j);
    end
    if (link2(i,j)>5000)||(link2(i,j)==0)
        link2(i,j)=link2(i-1,j);
    end
    if (link3(i,j)>5000)||(link3(i,j)==0)
        link3(i,j)=link3(i-1,j);
    end
    if (link4(i,j)>5000)||(link4(i,j)==0)
        link4(i,j)=link4(i-1,j);
    end
    if (link5(i,j)>5000)||(link5(i,j)==0)
        link5(i,j)=link5(i-1,j);
    end
    if (link6(i,j)>5000)||(link6(i,j)==0)
        link6(i,j)=link6(i-1,j);
    end
    if (link7(i,j)>5000)||(link7(i,j)==0)
        link7(i,j)=link7(i-1,j);
    end
    if (head(i,j)>5000)||(head(i,j)==0)
        head(i,j)=head(i-1,j);
    end
    end
end



%% 连杆位置坐标
x1=link1(:,2);y1=link1(:,3);
x2=link2(:,2);y2=link2(:,3);
x3=link3(:,2);y3=link3(:,3);
x4=link4(:,2);y4=link4(:,3);
x5=link5(:,2);y5=link5(:,3);
x6=link6(:,2);y6=link6(:,3);
x7=link7(:,2);y7=link7(:,3);
x8=head(:,2);y8=head(:,3);
x=[x1,x2,x3,x4,x5,x6,x7,x8];
y=[y1,y2,y3,y4,y5,y6,y7,y8];
N=8;
e=ones(N,1);

for i=1:length(x1)

    xrobot(i)=(1/N)*e'*x(i,:)';
    yrobot(i)=(1/N)*e'*y(i,:)';

end
plot(xrobot,yrobot)

dt=1/60;
%% 连杆角度
theta1=link1(:,7);
theta2=link2(:,7);
theta3=link3(:,7);
theta4=link4(:,7);
theta5=link5(:,7);
theta6=link6(:,7);
theta7=link7(:,7);
theta8=head(:,7);
theta=[theta1,theta2,theta3,theta4,theta5,theta6,theta7,theta8];
%% 计算速度
for i=1:size(x,1)-1
%     dx1(i)=(x1(i+1)-x1(i))/dt;
    dx(i,:)=(x(i+1,:)-x(i,:))./dt;
    dy(i,:)=(y(i+1,:)-y(i,:))./dt;
end
dx=[zeros(1,N);dx];dy=[zeros(1,N);dy];
%% 计算加速度
for i=1:size(dx,1)-1
%     dx1(i)=(x1(i+1)-x1(i))/dt;
    ddx(i,:)=(dx(i+1,:)-dx(i,:))./dt;
    ddy(i,:)=(dy(i+1,:)-dy(i,:))./dt;
end
ddx=[zeros(2,N);ddx];ddy=[zeros(2,N);ddy];
for i=1:length(x1)


    ddxrobot(i)=(1/N)*e'*ddx(i,:)';
    ddyrobot(i)=(1/N)*e'*ddy(i,:)';
end
ddp=[ddxrobot;ddyrobot];

N=8;m=0.08;
E=[ones(N,1),zeros(N,1);zeros(N,1),ones(N,1)];
fRR=N*m*ddp;
dt=1/60;

save('stra_dx','dx');
save('stra_dy','dy');
save('stra_ddx','ddx');
save('stra_ddy','ddy');
save('stra_theta','theta');
load('ff_smc_stra');
the1=ff(1,:);
plot(the1)
hold on
plot(theta1*pi/180)