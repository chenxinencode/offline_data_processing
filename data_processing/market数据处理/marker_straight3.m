%%这是导入marker点的数据，可以知道位置、速度、加速度，不需要计算。与刚体数据的结果吻合
%%需要将导出的marker数据删减前面3行
all_data=xlsread('sin_straight3');
t=all_data(4:end,2);
num=11;
l1=all_data(4:end,2+0*num+1:2+1*num);
l2=all_data(4:end,2+1*num+1:2+2*num);
l3=all_data(4:end,2+2*num+1:2+3*num);
l4=all_data(4:end,2+3*num+1:2+4*num);
l5=all_data(4:end,2+4*num+1:2+5*num);
l6=all_data(4:end,2+5*num+1:2+6*num);
l7=all_data(4:end,2+6*num+1:2+7*num);
l8=all_data(4:end,2+7*num+1:2+8*num);
head_marker=all_data(4:end,2+8*num+1:2+9*num);
%% marker位置
l1_x=l1(:,1);l1_y=l1(:,2);
l2_x=l2(:,1);l2_y=l2(:,2);
l3_x=l3(:,1);l3_y=l3(:,2);
l4_x=l4(:,1);l4_y=l4(:,2);
l5_x=l5(:,1);l5_y=l5(:,2);
l6_x=l6(:,1);l6_y=l6(:,2);
l7_x=l7(:,1);l7_y=l7(:,2);
l8_x=l8(:,1);l8_y=l8(:,2);
head_x=head_marker(:,1);head_y=head_marker(:,2);

x=[l1_x,l2_x,l3_x,l4_x,l5_x,l6_x,l7_x,head_x];
y=[l1_y,l2_y,l3_y,l4_y,l5_y,l6_y,l7_y,head_y];
N=8;
%% 求角度

angle1=atan2((head_y-l8_y),(head_x-l8_x));

%% 质心坐标
e=ones(N,1);

for i=1:length(l1_x)

    xrobot(i)=(1/N)*e'*x(i,:)';
    yrobot(i)=(1/N)*e'*y(i,:)';
    
end
theta=-90;
xy=[xrobot;yrobot];

% figure(1)
% plot(xrobot,yrobot);
% hold on
% scatter(xrobot(1),yrobot(1),'r*')
R=[cos(theta),-sin(theta);sin(theta),cos(theta)];
for i=1:length(l1_x)
    xy_rot(:,i)=R*xy(:,i);
end
xy_rot=xy_rot*0.001;
figure(2)
scatter(xy_rot(1,1),xy_rot(2,1),'r*','linewidth',2)
hold on
plot(xy_rot(1,:),xy_rot(2,:),'linewidth',2,'color','b');


legend('起点','质心')
xlabel('px/m');ylabel('py/m');
title('转弯90度实验')
%% marker点速度
l1_dx=l1(:,4);l1_dy=l1(:,5);l1_dxy=l1(:,7);
l2_dx=l2(:,4);l2_dy=l2(:,5);l2_dxy=l2(:,7);
l3_dx=l3(:,4);l3_dy=l3(:,5);l3_dxy=l3(:,7);
l4_dx=l4(:,4);l4_dy=l4(:,5);l4_dxy=l4(:,7);
l5_dx=l5(:,4);l5_dy=l5(:,5);l5_dxy=l5(:,7);
l6_dx=l6(:,4);l6_dy=l6(:,5);l6_dxy=l6(:,7);
l7_dx=l7(:,4);l7_dy=l7(:,5);l7_dxy=l7(:,7);
%% marker点加速度
l1_ddx=l1(:,8);l1_ddy=l1(:,9);l1_ddxy=l1(:,11);
l2_ddx=l2(:,8);l2_ddy=l2(:,9);l2_ddxy=l2(:,11);
l3_ddx=l3(:,8);l3_ddy=l3(:,9);l3_ddxy=l3(:,11);
l4_ddx=l4(:,8);l4_ddy=l4(:,9);l4_ddxy=l4(:,11);
l5_ddx=l5(:,8);l5_ddy=l5(:,9);l5_ddxy=l5(:,11);
l6_ddx=l6(:,8);l6_ddy=l6(:,9);l6_ddxy=l6(:,11);
l7_ddx=l7(:,8);l7_ddy=l7(:,9);l7_ddxy=l7(:,11);

head_marker_dx=head_marker(:,4);head_marker_dx=head_marker(:,5);head_marker_dx=head_marker(:,7);
% plot(l2_dx,'r')
% plot(l2_ddx,'r');
% plot()