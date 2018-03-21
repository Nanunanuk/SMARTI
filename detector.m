function [angle_stats]=detector(T_bot)

figure(3)
hold on
polar(0,90,'-k')
hold on


v=-[T_bot.dir_vector];
x=v(1,:);
y=v(2,:);
z=v(3,:);

theta=acosd(z);

% 0<alpha<90
a=x>0;
b=y>=0;
x1=a.*b.*x;
y1=a.*b.*y;
alpha1=asin(y1./sqrt(x1.^2+y1.^2));
alpha1(isnan(alpha1))=0;
% 90<alpha<180
a=x<=0;
b=y>0;
x2=a.*b.*x;
y2=a.*b.*y;
alpha2=pi()*ones(1,length(x))-asin(y2./sqrt(x2.^2+y2.^2));
alpha2(isnan(alpha2))=0;
% 180<alpha<270
a=x<0;
b=y<=0;
x3=a.*b.*x;
y3=a.*b.*y;
alpha3=pi()*ones(1,length(x))-asin(y3./sqrt(x3.^2+y3.^2));
alpha3(isnan(alpha3))=0;
% 270<alpha<360
a=x>=0;
b=y<0;
x4=a.*b.*x;
y4=a.*b.*y;
alpha4=2*pi()*ones(1,length(x))+asin(y4./sqrt(x4.^2+y4.^2));
alpha4(isnan(alpha4))=0;

alpha=alpha1+alpha2+alpha3+alpha4;

polar(alpha,theta,'.')

angles=[1:1:90];
[angle_stats] = histc([theta],angles);
figure(4)
hold on
bar(angle_stats/sum(angle_stats))

