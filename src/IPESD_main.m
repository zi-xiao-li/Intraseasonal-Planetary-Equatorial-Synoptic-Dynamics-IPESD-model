%新增经向模态
%% 也就是说在固定位置(warm pool)天气尺度对流建立-盛期-衰退过程中MJO的响应；
%% 这种过程主要表现在最大加热中心的位置(z)变化上，也即通过改变层云占比控制周期；
%% build up phase of the organized synoptic activity:a(t)从0-1,在20天达到最大值
%% decay: a(t) from 1 to 0
%% symetric about EQ，对于第一斜压模的经向谱系数，所有奇数项都为0

clear all; close all; clc


%*******************************************************

a        = 1/3;       %层云占比,后面可以根据研究内容进行t、X上的改变 
e        = 0.125;     %Froude number
L        = 5000;      %unit:km；暖池半径
y0       = 0;         %也就是距离赤道多少km，这里先选对称的，即y0=0
HT       = 16;        %对流层高度；
ro       = 0.65;      %波包典型length scale
x0       = 0.5;       %层云降水滞后距离
x_scale  = 1500;      %unit:km
y_scale  = 1500;      %unit:km
z_scale  = HT/pi;     %unit:km
X_scale  = 12000;     %unit:km
t_scale  = 3;         %unit:days；所以后面时间所有都要先化成days，然后再无量纲化；
u_scale  = 6.25;      %unit:m/s；行星尺度天气尺度一样的
v_scale  = 6.25;      %注意：因为行星尺度v太小，所以这个是和天气尺度扰动尺度一样的
w_scale  = 2.5;       %unit:cm/s；同上，是天气尺度
p_scale  = 312;       %unit:(m/s)^2；行星尺度天气尺度一样的
d        = 0.18;      %有争议，选取的是cooling相关，耗散速度要比drag慢
L        = L/x_scale; %注意这里仍旧选取天气尺度特征值，比较方便
c        = 5/u_scale; %envelope移动速度

%*******************************************************

N          = 31;                        %如：l经向展开0-4，则v=5,r=6；矩阵需要6+1=7;展开20阶需要20-1+3=22
dx         = 100/x_scale;               %空间：150km分辨率
dz         = 0.1/z_scale;               %高度：0.5km分辨率
dt         = 8.3*0.000694/t_scale;      %时间：16.6min或8.3min、4.15
xgrids1    = -10000/x_scale; 
xgrids2    = 10000/x_scale;
ygrids1    = -3000/y_scale; 
ygrids2    = 3000/y_scale;
zgrids1    = 0./z_scale; 
zgrids2    = HT/z_scale;
ngrids_x   = floor((xgrids2-xgrids1)/dx)+1;
ngrids_y   = floor((ygrids2-ygrids1)/dx)+1;
ngrids_z   = floor((zgrids2-zgrids1)/dz)+1;
ngrids_t   = 7200;                      %16.6min:积分3600步，41.5天;8.3min,7200步，41.5天；4.15min,3600*4
n_time     = 3;                         %自定积分时间是几个41.5d

for i=1:ngrids_z
    zcords(i)=zgrids1+(i-1)*dz;
end
for i=1:ngrids_y
    ycords(i)=ygrids1+(i-1)*dx;
end
for i=1:ngrids_x
    xcords(i)=xgrids1+(i-1)*dx;
end

%********************************************************





%                     1.正交基底
%=======================================================

syms y m

B1 = sqrt(1./(2.^m.*factorial(m).*sqrt(pi))).*hermiteH(m,y).*exp(-y.^2./2); 
B3 = sqrt(sqrt(3)./(2.^m*factorial(m).*sqrt(pi))).*hermiteH(m,sqrt(3).*y).*exp(-3.*y.^2./2);
B1 = matlabFunction(B1);
B3 = matlabFunction(B3);
%test the orthogonality
%这是标准正交向量
%注意：函数相乘时必须加上自变量，比如B1(m,y)*B3(m,y)，然后再matlabFunction(ans)
test(1)=double(int(B1(1,y)*B1(2,y),y,-inf,inf));
test(2)=double(int(B1(1,y)*B1(1,y),y,-inf,inf));
test(3)=double(int(B3(1,y)*B3(2,y),y,-inf,inf));
test(4)=double(int(B3(1,y)*B3(1,y),y,-inf,inf));




%                  2.fr、fl的经向谱系数
%========================================================

syms x z t

F   = cos(pi.*(x-c.*t)./(2.*L)).*(x>=-L&x<=L) + 0.*(x>L|x<-L); %envelope function
H   = sqrt(10).*exp(-(y+y0).^2);
Hy  = diff(H);

k   = 3/4.*F.^2.*a.*ro.*sin(x0./ro);
FU1 = k.*(2.*H.^2+y.*H.*Hy);
FT1 = k.*(5.*y.^2.*H.^2+4.*y.^3.*Hy.*H);
FU3 = -k.*(2.*H.^2+y.*H.*Hy);
FT3 = (k./9).*(15.*y.^2.*H.^2+4.*y.^3.*Hy.*H);

fr1 = matlabFunction(0.5.*(FT1+FU1));
fl1 = matlabFunction(0.5.*(FT1-FU1));
fr3 = matlabFunction(0.5.*(FT3+FU3));
fl3 = matlabFunction(0.5.*(FT3-FU3));

%% 经向分解，m
%  第一斜压模fr系数：fr1_m = sqrt(1./(2.^m.*factorial(m).*sqrt(pi))).*...
%			  				  int(fr1(t,x,y).*hermiteH(m,y).*exp(-y.^2./2),y,-inf,inf);
%  第一斜压模fl系数：fl1_m = sqrt(1./(2.^m.*factorial(m).*sqrt(pi))).*...
%						      int(fl1(t,x,y).*hermiteH(m,y).*exp(-y.^2./2),y,-inf,inf);
%  第三斜压模fr系数：fr3_m = sqrt(sqrt(3)./(2.^m.*factorial(m).*sqrt(pi))).*...
%							  int(fr3(t,x,y).*hermiteH(m,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);		
%  第三斜压模fl系数：fl3_m = sqrt(sqrt(3)./(2.^m.*factorial(m).*sqrt(pi))).*...
%						      int(fl3(t,x,y).*hermiteH(m,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);

q = 0;
fr1_0 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_0 = matlabFunction(fr1_0);
q = 1;
fr1_1 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_1 = matlabFunction(fr1_1);
q = 2;
fr1_2 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_2 = matlabFunction(fr1_2);
q = 3;
fr1_3 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_3 = matlabFunction(fr1_3);
q = 4;
fr1_4 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_4 = matlabFunction(fr1_4);
q = 5;
fr1_5 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_5 = matlabFunction(fr1_5);
q = 6;
fr1_6 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_6 = matlabFunction(fr1_6);
q = 7;
fr1_7 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_7 = matlabFunction(fr1_7);
q = 8;
fr1_8 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_8 = matlabFunction(fr1_8);
q = 9;
fr1_9 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_9 = matlabFunction(fr1_9);
q = 10;
fr1_10 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_10 = matlabFunction(fr1_10);
q = 11;
fr1_11 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_11 = matlabFunction(fr1_11);
q = 12;
fr1_12 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_12 = matlabFunction(fr1_12);
q = 13;
fr1_13 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_13 = matlabFunction(fr1_13);
q = 14;
fr1_14 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_14 = matlabFunction(fr1_14);
q = 15;
fr1_15 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_15 = matlabFunction(fr1_15);
q = 16;
fr1_16 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_16 = matlabFunction(fr1_16);
q = 17;
fr1_17 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_17 = matlabFunction(fr1_17);
q = 18;
fr1_18 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_18 = matlabFunction(fr1_18);
q = 19;
fr1_19 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_19 = matlabFunction(fr1_19);
q = 20;
fr1_20 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_20 = matlabFunction(fr1_20);
q = 21;
fr1_21 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_21 = matlabFunction(fr1_21);
q = 22;
fr1_22 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_22 = matlabFunction(fr1_22);
q = 23;
fr1_23 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_23 = matlabFunction(fr1_23);
q = 24;
fr1_24 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_24 = matlabFunction(fr1_24);
q = 25;
fr1_25 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_25 = matlabFunction(fr1_25);
q = 26;
fr1_26 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_26 = matlabFunction(fr1_26);
q = 27;
fr1_27 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_27 = matlabFunction(fr1_27);
q = 28;
fr1_28 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_28 = matlabFunction(fr1_28);
q = 29;
fr1_29 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_29 = matlabFunction(fr1_29);
q = 30;
fr1_30 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_30 = matlabFunction(fr1_30);
q = 31;
fr1_31 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr1_31 = matlabFunction(fr1_31);



q = 0;
fl1_0 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_0 = matlabFunction(fl1_0);
q = 1;
fl1_1 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_1 = matlabFunction(fl1_1);
q = 2;
fl1_2 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_2 = matlabFunction(fl1_2);
q = 3;
fl1_3 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_3 = matlabFunction(fl1_3);
q = 4;
fl1_4 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_4 = matlabFunction(fl1_4);
q = 5;
fl1_5 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_5 = matlabFunction(fl1_5);
q = 6;
fl1_6 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_6 = matlabFunction(fl1_6);
q = 7;
fl1_7 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_7 = matlabFunction(fl1_7);
q = 8;
fl1_8 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_8 = matlabFunction(fl1_8);
q = 9;
fl1_9 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_9 = matlabFunction(fl1_9);
q = 10;
fl1_10 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_10 = matlabFunction(fl1_10);
q = 11;
fl1_11 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_11 = matlabFunction(fl1_11);
q = 12;
fl1_12 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_12 = matlabFunction(fl1_12);
q = 13;
fl1_13 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_13 = matlabFunction(fl1_13);
q = 14;
fl1_14 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_14 = matlabFunction(fl1_14);
q = 15;
fl1_15 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_15 = matlabFunction(fl1_15);
q = 16;
fl1_16 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_16 = matlabFunction(fl1_16);
q = 17;
fl1_17 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_17 = matlabFunction(fl1_17);
q = 18;
fl1_18 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_18 = matlabFunction(fl1_18);
q = 19;
fl1_19 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_19 = matlabFunction(fl1_19);
q = 20;
fl1_20 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_20 = matlabFunction(fl1_20);
q = 21;
fl1_21 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_21 = matlabFunction(fl1_21);
q = 22;
fl1_22 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_22 = matlabFunction(fl1_22);
q = 23;
fl1_23 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_23 = matlabFunction(fl1_23);
q = 24;
fl1_24 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_24 = matlabFunction(fl1_24);
q = 25;
fl1_25 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_25 = matlabFunction(fl1_25);
q = 26;
fl1_26 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_26 = matlabFunction(fl1_26);
q = 27;
fl1_27 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_27 = matlabFunction(fl1_27);
q = 28;
fl1_28 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_28 = matlabFunction(fl1_28);
q = 29;
fl1_29 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_29 = matlabFunction(fl1_29);
q = 30;
fl1_30 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_30 = matlabFunction(fl1_30);
q = 31;
fl1_31 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl1(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl1_31 = matlabFunction(fl1_31);



q = 0;
fr3_0 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fr3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);	
fr3_0 = matlabFunction(fr3_0);
q = 1;
fr3_1 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fr3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);	
fr3_1 = matlabFunction(fr3_1);
q = 2;
fr3_2 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fr3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);	
fr3_2 = matlabFunction(fr3_2);
q = 3;
fr3_3 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fr3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);	
fr3_3 = matlabFunction(fr3_3);
q = 4;
fr3_4 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fr3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);	
fr3_4 = matlabFunction(fr3_4);
q = 5;
fr3_5 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fr3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);	
fr3_5 = matlabFunction(fr3_5);
q = 6;
fr3_6 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fr3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);	
fr3_6 = matlabFunction(fr3_6);
q = 7;
fr3_7 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fr3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);	
fr3_7 = matlabFunction(fr3_7);
q = 8;
fr3_8 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fr3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);	
fr3_8 = matlabFunction(fr3_8);
q = 9;
fr3_9 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fr3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);	
fr3_9 = matlabFunction(fr3_9);
q = 10;
fr3_10 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fr3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);	
fr3_10 = matlabFunction(fr3_10);
q = 11;
fr3_11 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fr3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);	
fr3_11 = matlabFunction(fr3_11);
q = 12;
fr3_12 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fr3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);	
fr3_12 = matlabFunction(fr3_12);
q = 13;
fr3_13 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fr3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);	
fr3_13 = matlabFunction(fr3_13);
q = 14;
fr3_14 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_14 = matlabFunction(fr3_14);
q = 15;
fr3_15 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_15 = matlabFunction(fr3_15);
q = 16;
fr3_16 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_16 = matlabFunction(fr3_16);
q = 17;
fr3_17 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_17 = matlabFunction(fr3_17);
q = 18;
fr3_18 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_18 = matlabFunction(fr3_18);
q = 19;
fr3_19 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_19 = matlabFunction(fr3_19);
q = 20;
fr3_20 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_20 = matlabFunction(fr3_20);
q = 21;
fr3_21 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_21 = matlabFunction(fr3_21);
q = 22;
fr3_22 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_22 = matlabFunction(fr3_22);
q = 23;
fr3_23 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_23 = matlabFunction(fr3_23);
q = 24;
fr3_24 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_24 = matlabFunction(fr3_24);
q = 25;
fr3_25 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_25 = matlabFunction(fr3_25);
q = 26;
fr3_26 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_26 = matlabFunction(fr3_26);
q = 27;
fr3_27 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_27 = matlabFunction(fr3_27);
q = 28;
fr3_28 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_28 = matlabFunction(fr3_28);
q = 29;
fr3_29 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_29 = matlabFunction(fr3_29);
q = 30;
fr3_30 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_30 = matlabFunction(fr3_30);
q = 31;
fr3_31 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fr3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fr3_31 = matlabFunction(fr3_31);



q = 0;
fl3_0 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);
fl3_0 = matlabFunction(fl3_0);
q = 1;
fl3_1 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);
fl3_1 = matlabFunction(fl3_1);
q = 2;
fl3_2 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);
fl3_2 = matlabFunction(fl3_2);
q = 3;
fl3_3 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);
fl3_3 = matlabFunction(fl3_3);
q = 4;
fl3_4 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);
fl3_4 = matlabFunction(fl3_4);
q = 5;
fl3_5 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);
fl3_5 = matlabFunction(fl3_5);
q = 6;
fl3_6 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);
fl3_6 = matlabFunction(fl3_6);
q = 7;
fl3_7 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);
fl3_7 = matlabFunction(fl3_7);
q = 8;
fl3_8 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);
fl3_8 = matlabFunction(fl3_8);
q = 9;
fl3_9 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);
fl3_9 = matlabFunction(fl3_9);
q = 10;
fl3_10 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);
fl3_10 = matlabFunction(fl3_10);
q = 11;
fl3_11 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);
fl3_11 = matlabFunction(fl3_11);
q = 12;
fl3_12 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);
fl3_12 = matlabFunction(fl3_12);
q = 13;
fl3_13 = sqrt(sqrt(3)./(2.^q.*factorial(q).*sqrt(pi))).*...
		int(fl3(t,x,y).*hermiteH(q,sqrt(3).*y).*exp(-3.*y.^2./2),y,-inf,inf);
fl3_13 = matlabFunction(fl3_13);
q = 14;
fl3_14 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_14 = matlabFunction(fl3_14);
q = 15;
fl3_15 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_15 = matlabFunction(fl3_15);
q = 16;
fl3_16 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_16 = matlabFunction(fl3_16);
q = 17;
fl3_17 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_17 = matlabFunction(fl3_17);
q = 18;
fl3_18 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_18 = matlabFunction(fl3_18);
q = 19;
fl3_19 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_19 = matlabFunction(fl3_19);
q = 20;
fl3_20 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_20 = matlabFunction(fl3_20);
q = 21;
fl3_21 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_21 = matlabFunction(fl3_21);
q = 22;
fl3_22 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_22 = matlabFunction(fl3_22);
q = 23;
fl3_23 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_23 = matlabFunction(fl3_23);
q = 24;
fl3_24 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_24 = matlabFunction(fl3_24);
q = 25;
fl3_25 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_25 = matlabFunction(fl3_25);
q = 26;
fl3_26 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_26 = matlabFunction(fl3_26);
q = 27;
fl3_27 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_27 = matlabFunction(fl3_27);
q = 28;
fl3_28 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_28 = matlabFunction(fl3_28);
q = 29;
fl3_29 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_29 = matlabFunction(fl3_29);
q = 30;
fl3_30 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_30 = matlabFunction(fl3_30);
q = 31;
fl3_31 = sqrt(1./(2.^q.*factorial(q).*sqrt(pi))).*...
			int(fl3(t,x,y).*hermiteH(q,y).*exp(-y.^2./2),y,-inf,inf);
fl3_31 = matlabFunction(fl3_31);



%=============================================================================
%=============================================================================
%
%
%                           1st Baroclinic Mode
%
%
%==============================================================================
%==============================================================================


%                       数值积分
%========================================================

% 1) lm,vm+1,rm+2  !!!!!!!!!
c1 = matlabFunction((1./(2.*m+3))./e);
c2 = matlabFunction(sqrt((m+1).*(m+2))./(2.*m+3));
c3 = matlabFunction((m+2)./(2.*m+3)); 
d0 = (1./sqrt(2)).*(sqrt(m+2)+(m+1)./sqrt(m+2));
d1 = (2.*sqrt((m+1)./(m+2)))./d0;
d1 = matlabFunction(d1./e);
d2 = matlabFunction((sqrt((m+1)./(m+2)))./d0);
d0 = matlabFunction(1./d0);

% 1.1) l0,v1,r2:q=0  !!!!!!!!!
q    = 0;
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_2(n.*dt,xcords(ii))+...
					c3(q).*fl1_0(n.*dt,xcords(ii));
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); %16.6min则87;8.6min则87*2;
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_0(n.*dt,xcords(ii))-...
							d0(q).*fr1_2(n.*dt,xcords(ii));
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l01.txt','a'); 
		fid2 = fopen('C:\Users\pc\Desktop\data\v11.txt','a');
		fid3 = fopen('C:\Users\pc\Desktop\data\r21.txt','a');
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l01.txt');
vvv = load('C:\Users\pc\Desktop\data\v11.txt');
rrr = load('C:\Users\pc\Desktop\data\r21.txt');
%存储l、r以画p、u时间-经度图，记住这个只需在第一次的时候设置
l        = zeros(ngrids_x,length(lll)/ngrids_x,N);
v        = zeros(ngrids_x,length(lll)/ngrids_x,N);
r        = zeros(ngrids_x,length(lll)/ngrids_x,N);
l(:,:,1) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);
v(:,:,2) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);
r(:,:,3) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);



% 1.2) l1,v2,r3:q=1   !!!!!!!!!
q    = 1;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1);
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx));
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l11.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v21.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r31.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l11.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v21.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r31.txt');%改
l(:,:,2) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,3) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,4) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改



% 1.3) l2,v3,r4:q=2   !!!!!!!!!
q    = 2;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_4(n.*dt,xcords(ii))+... %改
					c3(q).*fl1_2(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_2(n.*dt,xcords(ii))-...%改
							d0(q).*fr1_4(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l21.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v31.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r41.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l21.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v31.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r41.txt');%改
l(:,:,3) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,4) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,5) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改



% 1.4) l3,v4,r5:q=3   !!!!!!!!!
q    = 3;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1);
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx));
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l31.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v41.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r51.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l31.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v41.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r51.txt');%改
l(:,:,4) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,5) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,6) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改



% 1.5) l4,v5,r6:q=4   !!!!!!!!!
q    = 4;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_6(n.*dt,xcords(ii))+... %改
					c3(q).*fl1_4(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_4(n.*dt,xcords(ii))-...%改
							d0(q).*fr1_6(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l41.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v51.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r61.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l41.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v51.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r61.txt');%改
l(:,:,5) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,6) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,7) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改

% 1.6) l5,v6,r7:q=5   !!!!!!!!!
q    = 5;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1);
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx));
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l51.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v61.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r71.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l51.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v61.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r71.txt');%改
l(:,:,6) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,7) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,8) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改

% 1.7) l6,v7,r8:q=6   !!!!!!!!!
q    = 6;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_8(n.*dt,xcords(ii))+... %改
					c3(q).*fl1_6(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_6(n.*dt,xcords(ii))-...%改
							d0(q).*fr1_8(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l61.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v71.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r81.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l61.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v71.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r81.txt');%改
l(:,:,7) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,8) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,9) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.8) l7,v8,r9:q=7   !!!!!!!!!
q    = 7;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1);
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx));
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l71.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v81.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r91.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l71.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v81.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r91.txt');%改
l(:,:,8) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,9) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,10) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.9) l8,v9,r10:q=8   !!!!!!!!!
q    = 8;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_10(n.*dt,xcords(ii))+... %改
					c3(q).*fl1_8(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_8(n.*dt,xcords(ii))-...%改
							d0(q).*fr1_10(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l81.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v91.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r101.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l81.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v91.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r101.txt');%改
l(:,:,9) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,10) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,11) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.11) l10,v11,r12:q=10   !!!!!!!!!
q    = 10;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_12(n.*dt,xcords(ii))+... %改
					c3(q).*fl1_10(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_10(n.*dt,xcords(ii))-...%改
							d0(q).*fr1_12(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l101.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v111.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r121.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l101.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v111.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r121.txt');%改
l(:,:,11) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,12) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,13) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.13) l12,v13,r14:q=12   !!!!!!!!!
q    = 12;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_14(n.*dt,xcords(ii))+... %改
					c3(q).*fl1_12(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_12(n.*dt,xcords(ii))-...%改
							d0(q).*fr1_14(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l121.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v131.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r141.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l121.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v131.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r141.txt');%改
l(:,:,13) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,14) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,15) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.15) l14,v15,r16:q=14   !!!!!!!!!
q    = 14;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_16(n.*dt,xcords(ii))+... %改
					c3(q).*fl1_14(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_14(n.*dt,xcords(ii))-...%改
							d0(q).*fr1_16(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l141.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v151.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r161.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l141.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v151.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r161.txt');%改
l(:,:,15) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,16) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,17) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.17) l16,v17,r18:q=16   !!!!!!!!!
q    = 16;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_18(n.*dt,xcords(ii))+... %改
					c3(q).*fl1_16(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_16(n.*dt,xcords(ii))-...%改
							d0(q).*fr1_18(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l161.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v171.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r181.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l161.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v171.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r181.txt');%改
l(:,:,17) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,18) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,19) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.19) l18,v19,r20:q=18   !!!!!!!!!
q    = 18;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_20(n.*dt,xcords(ii))+... %改
					c3(q).*fl1_18(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_18(n.*dt,xcords(ii))-...%改
							d0(q).*fr1_20(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l181.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v191.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r201.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l181.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v191.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r201.txt');%改
l(:,:,19) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,20) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,21) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.21) l20,v21,r22:q=20   !!!!!!!!!
q    = 20;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_22(n.*dt,xcords(ii))+... %改
					c3(q).*fl1_20(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_20(n.*dt,xcords(ii))-...%改
							d0(q).*fr1_22(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l201.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v211.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r221.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l201.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v211.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r221.txt');%改
l(:,:,21) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,22) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,23) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.23) l22,v23,r24:q=22   !!!!!!!!!
q    = 22;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_24(n.*dt,xcords(ii))+... %改
					c3(q).*fl1_22(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_22(n.*dt,xcords(ii))-...%改
							d0(q).*fr1_24(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l221.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v231.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r241.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l221.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v231.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r241.txt');%改
l(:,:,23) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,24) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,25) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.25) l24,v25,r26:q=24   !!!!!!!!!
q    = 24;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_26(n.*dt,xcords(ii))+... %改
					c3(q).*fl1_24(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_24(n.*dt,xcords(ii))-...%改
							d0(q).*fr1_26(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l241.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v251.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r261.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l241.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v251.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r261.txt');%改
l(:,:,25) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,26) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,27) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.27) 
q    = 26;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_28(n.*dt,xcords(ii))+... %改
					c3(q).*fl1_26(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_26(n.*dt,xcords(ii))-...%改
							d0(q).*fr1_28(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l261.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v271.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r281.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l261.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v271.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r281.txt');%改
l(:,:,27) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,28) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,29) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.29) 
q    = 28;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr1_30(n.*dt,xcords(ii))+... %改
					c3(q).*fl1_28(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl1_28(n.*dt,xcords(ii))-...%改
							d0(q).*fr1_30(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l281.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v291.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r301.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l281.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v291.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r301.txt');%改
l(:,:,29) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,30) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,31) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 2) v0  !!!!!!!!!
n  = 0;
nn = 1;
while n<=ngrids_t*n_time
	n_mod = mod(n,87*2); 
	if(n_mod==0)
		for ii=2:ngrids_x-1
			%asymmetric:v(ii,nn,1)=-sqrt(2).*fr1_1(n.*dt,xcords(ii));
			v(ii,nn,1)=0;
		end
		nn = nn+1;
	end
	n = n+1;
end


% 3) r0  !!!!!!!!!
n    = 0;
h1   = 1/e;
rr   = zeros(ngrids_x,2);
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = -h1.*((rr(ii+1,1)-rr(ii-1,1))./(2.*dx))-d.*rr(ii,1)+...
				fr1_0(n.*dt,xcords(ii)); %改
		rr(ii,2) = rr(ii,1)+c4.*dt;
		rr(ii,1) = rr(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(rr(ii,1));
		ll(ii,2) = smoothdata(rr(ii,2));
	end
	if(n_mod==0)
		fid3 = fopen('C:\Users\pc\Desktop\data\r01.txt','a');%改
		fprintf(fid3,'%e\n',rr(:,1));
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
rrr = load('C:\Users\pc\Desktop\data\r01.txt');%改
r(:,:,1) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改



%                  Sum of Meridinal Modes
%========================================================

s1 = [];
s2 = [];
s3 = [];
for ii=1:size(l,2)
	s3 = 0;
	for jj=1:N-2  % ！！
		s1 = B1(jj-1,ycords);
		s2 = l(:,ii,jj);
		s3 = s3+s1.*s2;
	end
	l_N(:,:,ii) = s3;
end
s1 = [];
s2 = [];
s3 = [];
for ii=1:size(l,2)
	s3 = 0;
	for jj=1:N-2 
		s1 = B1(jj-1,ycords);
		s2 = r(:,ii,jj);
		s3 = s3+s1.*s2;
	end
	r_N(:,:,ii) = s3;
end
s1 = [];
s2 = [];
s3 = [];
for ii=1:size(l,2)
	s3 = 0;
	for jj=1:N-2
		s1 = B1(jj-1,ycords);
		s2 = v(:,ii,jj);
		s3 = s3+s1.*s2;
	end
	v_N(:,:,ii) = s3;
end
u_N1 = r_N-l_N;
p_N1 = r_N+l_N;
v_N1 = v_N;






%=============================================================================
%=============================================================================
%
%
%                           3th Baroclinic Mode
%
%
%==============================================================================
%==============================================================================

%             数值积分:B1改B3;fr1、fl1改fl1、fl3
%========================================================


% 1) lm,vm+1,rm+2  !!!!!!!!!
c1 = matlabFunction((1./3.*(1./(2.*m+3)))./e);
c2 = matlabFunction(sqrt((m+1).*(m+2))./(2.*m+3));
c3 = matlabFunction((m+2)./(2.*m+3)); 
d0 = (sqrt(6)./18).*(sqrt(m+2)+(m+1)./sqrt(m+2));
d1 = ((2./3).*sqrt((m+1)./(m+2)))./d0;
d1 = matlabFunction(d1./e);
d2 = matlabFunction((sqrt((m+1)./(m+2)))./d0);
d0 = matlabFunction(1./d0);

% 1.1) l0,v1,r2:q=0  !!!!!!!!!
q    = 0;
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_2(n.*dt,xcords(ii))+...
					c3(q).*fl3_0(n.*dt,xcords(ii));
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_0(n.*dt,xcords(ii))-...
							d0(q).*fr3_2(n.*dt,xcords(ii));
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l03.txt','a'); 
		fid2 = fopen('C:\Users\pc\Desktop\data\v13.txt','a');
		fid3 = fopen('C:\Users\pc\Desktop\data\r23.txt','a');
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l03.txt');
vvv = load('C:\Users\pc\Desktop\data\v13.txt');
rrr = load('C:\Users\pc\Desktop\data\r23.txt');
%存储l、r以画p、u时间-经度图，记住这个只需在第一次的时候设置
l        = zeros(ngrids_x,length(lll)/ngrids_x,N);
v        = zeros(ngrids_x,length(lll)/ngrids_x,N);
r        = zeros(ngrids_x,length(lll)/ngrids_x,N);
l(:,:,1) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);
v(:,:,2) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);
r(:,:,3) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);



% 1.2) l1,v2,r3:q=1   !!!!!!!!!
q    = 1;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1);
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx));
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l13.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v23.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r33.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l13.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v23.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r33.txt');%改
l(:,:,2) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,3) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,4) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改



% 1.3) l2,v3,r4:q=2   !!!!!!!!!
q    = 2;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_4(n.*dt,xcords(ii))+... %改
					c3(q).*fl3_2(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_2(n.*dt,xcords(ii))-...%改
							d0(q).*fr3_4(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l23.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v33.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r43.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l23.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v33.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r43.txt');%改
l(:,:,3) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,4) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,5) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改



% 1.4) l3,v4,r5:q=3   !!!!!!!!!
q    = 3;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1);
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx));
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l33.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v43.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r53.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l33.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v43.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r53.txt');%改
l(:,:,4) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,5) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,6) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改



% 1.5) l4,v5,r6:q=4   !!!!!!!!!
q    = 4;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_6(n.*dt,xcords(ii))+... %改
					c3(q).*fl3_4(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_4(n.*dt,xcords(ii))-...%改
							d0(q).*fr3_6(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l43.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v53.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r63.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l43.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v53.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r63.txt');%改
l(:,:,5) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,6) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,7) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.6) l5,v6,r7:q=5   !!!!!!!!!
q    = 5;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1);
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx));
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l53.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v63.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r73.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l53.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v63.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r73.txt');%改
l(:,:,6) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,7) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,8) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.7) l6,v7,r8:q=6   !!!!!!!!!
q    = 6;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_8(n.*dt,xcords(ii))+... %改
					c3(q).*fl3_6(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_6(n.*dt,xcords(ii))-...%改
							d0(q).*fr3_8(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l63.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v73.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r83.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l63.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v73.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r83.txt');%改
l(:,:,7) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,8) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,9) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.8) l7,v8,r9:q=7   !!!!!!!!!
q    = 7;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1);
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx));
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l73.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v83.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r93.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l73.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v83.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r93.txt');%改
l(:,:,8) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,9) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,10) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.9) l8,v9,r10:q=8   !!!!!!!!!
q    = 8;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_10(n.*dt,xcords(ii))+... %改
					c3(q).*fl3_8(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_8(n.*dt,xcords(ii))-...%改
							d0(q).*fr3_10(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l83.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v93.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r103.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l83.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v93.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r103.txt');%改
l(:,:,9) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,10) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,11) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.11) l10,v11,r12:q=10   !!!!!!!!!
q    = 10;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_12(n.*dt,xcords(ii))+... %改
					c3(q).*fl3_10(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_10(n.*dt,xcords(ii))-...%改
							d0(q).*fr3_12(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l103.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v113.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r123.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l103.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v113.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r123.txt');%改
l(:,:,11) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,12) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,13) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.13) l12,v13,r14:q=12   !!!!!!!!!
q    = 12;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_14(n.*dt,xcords(ii))+... %改
					c3(q).*fl3_12(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_12(n.*dt,xcords(ii))-...%改
							d0(q).*fr3_14(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l123.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v133.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r143.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l123.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v133.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r143.txt');%改
l(:,:,13) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,14) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,15) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.15) l14,v15,r16:q=14   !!!!!!!!!
q    = 14;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_16(n.*dt,xcords(ii))+... %改
					c3(q).*fl3_14(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_14(n.*dt,xcords(ii))-...%改
							d0(q).*fr3_16(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l143.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v153.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r163.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l143.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v153.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r163.txt');%改
l(:,:,15) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,16) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,17) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.17) l16,v17,r18:q=16   !!!!!!!!!
q    = 16;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_18(n.*dt,xcords(ii))+... %改
					c3(q).*fl3_16(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_16(n.*dt,xcords(ii))-...%改
							d0(q).*fr3_18(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l163.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v173.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r183.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l163.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v173.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r183.txt');%改
l(:,:,17) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,18) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,19) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.19) l18,v19,r20:q=18   !!!!!!!!!
q    = 18;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_20(n.*dt,xcords(ii))+... %改
					c3(q).*fl3_18(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_18(n.*dt,xcords(ii))-...%改
							d0(q).*fr3_20(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l183.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v193.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r203.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l183.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v193.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r203.txt');%改
l(:,:,19) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,20) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,21) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.21) l20,v21,r22:q=20   !!!!!!!!!
q    = 20;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_22(n.*dt,xcords(ii))+... %改
					c3(q).*fl3_20(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_20(n.*dt,xcords(ii))-...%改
							d0(q).*fr3_22(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l203.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v213.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r223.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l203.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v213.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r223.txt');%改
l(:,:,21) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,22) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,23) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.23) l22,v23,r24:q=22   !!!!!!!!!
q    = 22;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_24(n.*dt,xcords(ii))+... %改
					c3(q).*fl3_22(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_22(n.*dt,xcords(ii))-...%改
							d0(q).*fr3_24(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l223.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v233.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r243.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l223.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v233.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r243.txt');%改
l(:,:,23) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,24) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,25) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.25) l24,v25,r26:q=24   !!!!!!!!!
q    = 24;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_26(n.*dt,xcords(ii))+... %改
					c3(q).*fl3_24(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_24(n.*dt,xcords(ii))-...%改
							d0(q).*fr3_26(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l243.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v253.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r263.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l243.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v253.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r263.txt');%改
l(:,:,25) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,26) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,27) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.27) 
q    = 26;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_28(n.*dt,xcords(ii))+... %改
					c3(q).*fl3_26(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_26(n.*dt,xcords(ii))-...%改
							d0(q).*fr3_28(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l263.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v273.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r283.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l263.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v273.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r283.txt');%改
l(:,:,27) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,28) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,29) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改


% 1.29) 
q    = 28;%改
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = c1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))-d.*ll(ii,1)+...
				c2(q).*fr3_30(n.*dt,xcords(ii))+... %改
					c3(q).*fl3_28(n.*dt,xcords(ii));%改
		ll(ii,2) = ll(ii,1)+c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if(n_mod==0)
		for ii=2:ngrids_x-1
			vv(ii,1) = d1(q).*((ll(ii+1,1)-ll(ii-1,1))./(2.*dx))+...
						d2(q).*fl3_28(n.*dt,xcords(ii))-...%改
							d0(q).*fr3_30(n.*dt,xcords(ii));%改
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l283.txt','a');%改 
		fid2 = fopen('C:\Users\pc\Desktop\data\v293.txt','a');%改
		fid3 = fopen('C:\Users\pc\Desktop\data\r303.txt','a');%改
		fprintf(fid1,'%e\n',ll(:,1));  
		fprintf(fid2,'%e\n',vv(:,1));
		fprintf(fid3,'%e\n',sqrt((q+1)./(q+2)).*ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
lll = load('C:\Users\pc\Desktop\data\l283.txt');%改
vvv = load('C:\Users\pc\Desktop\data\v293.txt');%改
rrr = load('C:\Users\pc\Desktop\data\r303.txt');%改
l(:,:,29) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);%改
v(:,:,30) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);%改
r(:,:,31) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改




% 2) v0  !!!!!!!!!
n  = 0;
nn = 1;
while n<=ngrids_t*n_time
	n_mod = mod(n,87*2); 
	if(n_mod==0)
		for ii=2:ngrids_x-1
			%symetric:v(ii,nn,1)=-sqrt(2).*fr1_1(n.*dt,xcords(ii));
			v(ii,nn,1)=0;
		end
		nn = nn+1;
	end
	n = n+1;
end


% 3) r0  !!!!!!!!!
n    = 0;
h1   = 1/(3.*e);
rr   = zeros(ngrids_x,2);
rrr  = [];
while n<=ngrids_t*n_time
	for ii=2:ngrids_x-1
		c4 = -h1.*((rr(ii+1,1)-rr(ii-1,1))./(2.*dx))-d.*rr(ii,1)+...
				fr3_0(n.*dt,xcords(ii)); %改
		rr(ii,2) = rr(ii,1)+c4.*dt;
		rr(ii,1) = rr(ii,2);
	end
	%每天输出为一个文件
	n_mod = mod(n,87*2); 
	k_mod = mod(n,87);
	if(k_mod==0)
		ll(ii,1) = smoothdata(rr(ii,1));
		ll(ii,2) = smoothdata(rr(ii,2));
	end
	if(n_mod==0)
		fid3 = fopen('C:\Users\pc\Desktop\data\r03.txt','a');%改
		fprintf(fid3,'%e\n',rr(:,1));
		fclose(fid3);
	end
	n = n+1;
end
%重新读入成矩阵形式
rrr = load('C:\Users\pc\Desktop\data\r03.txt');%改
r(:,:,1) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);%改



%                  Sum of Meridinal Modes
%========================================================

s1 = [];
s2 = [];
s3 = [];
for ii=1:size(l,2)
	s3 = 0;
	for jj=1:N-2
		s1 = B3(jj-1,ycords);
		s2 = l(:,ii,jj);
		s3 = s3+s1.*s2;  
	end
	l_N(:,:,ii) = s3;
end
s1 = [];
s2 = [];
s3 = [];
for ii=1:size(l,2)
	s3 = 0;
	for jj=1:N-2
		s1 = B3(jj-1,ycords);
		s2 = r(:,ii,jj);
		s3 = s3+s1.*s2;
	end
	r_N(:,:,ii) = s3;
end
s1 = [];
s2 = [];
s3 = [];
for ii=1:size(l,2)
	s3 = 0;
	for jj=1:N-2
		s1 = B3(jj-1,ycords);
		s2 = v(:,ii,jj);
		s3 = s3+s1.*s2;
	end
	v_N(:,:,ii) = s3;
end
u_N3 = r_N-l_N; 
p_N3 = r_N+l_N;
v_N3 = v_N;




%=============================================================================
%=============================================================================
%
%
%                            mode 1  +mode 3
%
%
%==============================================================================
%==============================================================================

U         = [];
V         = [];
P         = [];

for ii=1:size(v_N3,1)
	for jj=1:size(v_N3,2)
		for kk=1:size(v_N3,3)
			U(ii,jj,:,kk) = u_N1(ii,jj,kk).*cos(zcords)+u_N3(ii,jj,kk).*cos(3.*zcords);
			V(ii,jj,:,kk) = v_N1(ii,jj,kk).*cos(zcords)+v_N3(ii,jj,kk).*cos(3.*zcords);
			P(ii,jj,:,kk) = p_N1(ii,jj,kk).*cos(zcords)+p_N3(ii,jj,kk).*cos(3.*zcords);
		end
	end
end































%                       PLOT: fig.3.
%========================================================


[~,I0]=min(abs(ycords(:)-0./y_scale));
[~,I1]=min(abs(zcords(:)-0./z_scale));
[~,I2]=min(abs(zcords(:)-2./z_scale));
[~,I3]=min(abs(zcords(:)-4./z_scale));
[~,I4]=min(abs(zcords(:)-12./z_scale));



%*******************************************************************

days      = 80;                             %选择要查看的天数
H         = 10;                             %五点平滑迭代次数
S         = 1./2;                           %五点平滑系数
xmin      = -6;                             %x轴取最小值
xmax      = 6;                              %x轴取最大值
ymin      = -3000/x_scale;                  %y轴取最小值
ymax      = 3000/x_scale;                   %y轴取最大值
xtick1    = num2str(xmin.*x_scale./1000);   %x轴最小值地理标号
xtick2    = num2str(xmax.*x_scale./1000);   %x轴最大值地理标号
ytick1    = num2str(ymin.*y_scale./1000);   %y轴最小值地理标号
ytick2    = num2str(ymax.*y_scale./1000);   %y轴最大值地理标号
lowC      = -1;                             %colorbar取值
highC     = 1;                              %colorbar取值
LevelStep = 1;                              %等值线间距
LineWidth = 0.8;                            %等值线宽
stride1   = 5;
stride2   = 2;

%*******************************************************************


%                       PLOT: fig.5.
%========================================================
figure

[xx1 yy]  = meshgrid(linspace(xgrids1,xgrids2,ngrids_x),...
	linspace(ygrids1,ygrids2,ngrids_y)); 
s         = [];
s(:,:)    = P(:,:,I1,days);
s_sm      = s;
for kk = 1:H
	for ii=2:size(s,1)-1
		for jj=2:size(s,2)-1
			s_sm(ii,jj) = s(ii,jj)+0.5.*S.*(1-S).*(s(ii+1,jj)+s(ii-1,jj)+s(ii,jj+1)+s(ii,jj-1)-4.*s(ii,jj))...
				+0.25.*S.^2.*(s(ii+1,jj+1)+s(ii+1,jj-1)+s(ii-1,jj+1)+s(ii-1,jj-1)-4.*s(ii,jj));
		end
	end
	s = s_sm;
end
subplot(4,1,1)
aa1 = s_sm;
aa2 = s_sm;
aa1(aa1<0) = 0;
aa2(aa2>0) = 0;
contour(xx1,yy,aa1','-k','LineWidth',LineWidth);
hold on
contour(xx1,yy,aa2','--k','LineWidth',LineWidth);
title('(a) 0km')
set(gca,'XTick',[xmin 0 xmax])
set(gca,'XTicklabel',{xtick1,'0',xtick2}) 
set(gca,'YTick',[ymin 0 ymax])
set(gca,'YTicklabel',{ytick1,'0',ytick2})
set(gca,'clim',[lowC highC])
axis([xmin xmax ymin ymax])
hold on
s1      = [];
u_p     = U(:,:,I1,days);
v_p     = V(:,:,I1,days);
u_pp    = u_p([1:stride1:size(u_p,1)],[1:stride2:size(u_p,2)]);
v_pp    = v_p([1:stride1:size(v_p,1)],[1:stride2:size(v_p,2)]);
uv      = sqrt(u_p.^2+v_p.^2);
s1(1)   = max(max(uv));

[x_hrz_wnd y_hrz_wnd] = meshgrid(linspace(xgrids1,xgrids2,size(u_pp,2)),...
	linspace(ygrids1,ygrids2,size(u_pp,1)));
quiver(x_hrz_wnd,y_hrz_wnd,u_pp,v_pp,'-k','MaxHeadSize',3,'AutoScale','on')



s      = [];
S      = 1/2;
s(:,:) = P(:,:,I2,days);
s_sm   = s;
for kk = 1:H
	for ii=2:size(s,1)-1
		for jj=2:size(s,2)-1
			s_sm(ii,jj) = s(ii,jj)+0.5.*S.*(1-S).*(s(ii+1,jj)+s(ii-1,jj)+s(ii,jj+1)+s(ii,jj-1)-4.*s(ii,jj))...
				+0.25.*S.^2.*(s(ii+1,jj+1)+s(ii+1,jj-1)+s(ii-1,jj+1)+s(ii-1,jj-1)-4.*s(ii,jj));
		end
	end
	s = s_sm;
end
hold on
subplot(4,1,2)
aa1 = s_sm;
aa2 = s_sm;
aa1(aa1<0) = 0;
aa2(aa2>0) = 0;
contour(xx1,yy,aa1','-k','LineWidth',LineWidth);
hold on
contour(xx1,yy,aa2','--k','LineWidth',LineWidth);
title('(b) 2km')
set(gca,'XTick',[xmin 0 xmax])
set(gca,'XTicklabel',{xtick1,'0',xtick2}) 
set(gca,'YTick',[ymin 0 ymax])
set(gca,'YTicklabel',{ytick1,'0',ytick2})
set(gca,'clim',[lowC highC])
axis([xmin xmax ymin ymax])
hold on
u_p     = U(:,:,I2,days);
v_p     = V(:,:,I2,days);
u_pp    = u_p([1:stride1:size(u_p,1)],[1:stride2:size(u_p,2)]);
v_pp    = v_p([1:stride1:size(v_p,1)],[1:stride2:size(v_p,2)]);
uv      = sqrt(u_p.^2+v_p.^2);
s1(2)   = max(max(uv));

[x_hrz_wnd y_hrz_wnd] = meshgrid(linspace(xgrids1,xgrids2,size(u_pp,2)),...
	linspace(ygrids1,ygrids2,size(u_pp,1)));
quiver(x_hrz_wnd,y_hrz_wnd,u_pp,v_pp,'-k','MaxHeadSize',3,'AutoScale','on')



s      = [];
S      = 1/2;
s(:,:) = P(:,:,I3,days);
s_sm   = s;
for kk = 1:H
	for ii=2:size(s,1)-1
		for jj=2:size(s,2)-1
			s_sm(ii,jj) = s(ii,jj)+0.5.*S.*(1-S).*(s(ii+1,jj)+s(ii-1,jj)+s(ii,jj+1)+s(ii,jj-1)-4.*s(ii,jj))...
				+0.25.*S.^2.*(s(ii+1,jj+1)+s(ii+1,jj-1)+s(ii-1,jj+1)+s(ii-1,jj-1)-4.*s(ii,jj));
		end
	end
	s = s_sm;
end
hold on
subplot(4,1,3)
aa1 = s_sm;
aa2 = s_sm;
aa1(aa1<0) = 0;
aa2(aa2>0) = 0;
contour(xx1,yy,aa1','-k','LineWidth',LineWidth);
hold on
contour(xx1,yy,aa2','--k','LineWidth',LineWidth);
title('(c) 4km')
set(gca,'XTick',[xmin 0 xmax])
set(gca,'XTicklabel',{xtick1,'0',xtick2}) 
set(gca,'YTick',[ymin 0 ymax])
set(gca,'YTicklabel',{ytick1,'0',ytick2})
set(gca,'clim',[lowC highC])
axis([xmin xmax ymin ymax])
hold on
u_p     = U(:,:,I3,days);
v_p     = V(:,:,I3,days);
u_pp    = u_p([1:stride1:size(u_p,1)],[1:stride2:size(u_p,2)]);
v_pp    = v_p([1:stride1:size(v_p,1)],[1:stride2:size(v_p,2)]);
uv      = sqrt(u_p.^2+v_p.^2);
s1(3)   = max(max(uv));

[x_hrz_wnd y_hrz_wnd] = meshgrid(linspace(xgrids1,xgrids2,size(u_pp,2)),...
	linspace(ygrids1,ygrids2,size(u_pp,1)));

quiver(x_hrz_wnd,y_hrz_wnd,u_pp,v_pp,'-k','MaxHeadSize',3,'AutoScale','on')


s      = [];
S      = 1/2;
s(:,:) = P(:,:,I4,days);
s_sm   = s;
for kk = 1:H
	for ii=2:size(s,1)-1
		for jj=2:size(s,2)-1
			s_sm(ii,jj) = s(ii,jj)+0.5.*S.*(1-S).*(s(ii+1,jj)+s(ii-1,jj)+s(ii,jj+1)+s(ii,jj-1)-4.*s(ii,jj))...
				+0.25.*S.^2.*(s(ii+1,jj+1)+s(ii+1,jj-1)+s(ii-1,jj+1)+s(ii-1,jj-1)-4.*s(ii,jj));
		end
	end
	s = s_sm;
end
hold on
subplot(4,1,4)
aa1 = s_sm;
aa2 = s_sm;
aa1(aa1<0) = 0;
aa2(aa2>0) = 0;
contour(xx1,yy,aa1','-k','LineWidth',LineWidth);
hold on
contour(xx1,yy,aa2','--k','LineWidth',LineWidth);
title('(d) 12km')
xlabel('Zonal(x 1000km)')
ylabel('Horizonal(x 1000km)')
%set(gca,'position',[0.3 0.1 0.4 0.23])
set(gca,'XTick',[xmin 0 xmax])
set(gca,'XTicklabel',{xtick1,'0',xtick2}) 
set(gca,'YTick',[ymin 0 ymax])
set(gca,'YTicklabel',{ytick1,'0',ytick2})
set(gca,'clim',[lowC highC])
axis([xmin xmax ymin ymax])
hold on
u_p     = U(:,:,I4,days);
v_p     = V(:,:,I4,days);
u_pp    = u_p([1:stride1:size(u_p,1)],[1:stride2:size(u_p,2)]);
v_pp    = v_p([1:stride1:size(v_p,1)],[1:stride2:size(v_p,2)]);
uv      = sqrt(u_p.^2+v_p.^2);
s1(3)   = max(max(uv));

[x_hrz_wnd y_hrz_wnd] = meshgrid(linspace(xgrids1,xgrids2,size(u_pp,2)),...
	linspace(ygrids1,ygrids2,size(u_pp,1)));

quiver(x_hrz_wnd,y_hrz_wnd,u_pp,v_pp,'-k','MaxHeadSize',3,'AutoScale','on')




%                       PLOT: fig.6.
%========================================================
figure 

[xx2 yy2]  = meshgrid(linspace(xgrids1,xgrids2,ngrids_x),...
	linspace(zgrids1,zgrids2,ngrids_z)); 
s         = [];
s(:,:)    = U(:,I0,:,days);
s_sm      = s;
H         = 100;                             %五点平滑迭代次数
kk        = 1;
for kk = 1:H
	for ii=2:size(s,1)-1
		for jj=2:size(s,2)-1
			s_sm(ii,jj) = s(ii,jj)+0.5.*S.*(1-S).*(s(ii+1,jj)+s(ii-1,jj)+s(ii,jj+1)+s(ii,jj-1)-4.*s(ii,jj))...
				+0.25.*S.^2.*(s(ii+1,jj+1)+s(ii+1,jj-1)+s(ii-1,jj+1)+s(ii-1,jj-1)-4.*s(ii,jj));
		end
	end
	s = s_sm;
end

subplot(2,1,1)
contourf(xx2,yy2,s_sm');
colormap gray
title('(a)')
xlabel('x(1000km)')
ylabel('z(km)')
%set(gca,'clim',[lowC highC])
set(gca,'XTick',[-6 -3 0 3 6])
set(gca,'XTicklabel',{'-9','-4.5','0','4.5','9'})
set(gca,'YTick',[0 max(zcords)])
set(gca,'YTicklabel',{'0','16'})
set(gca,'clim',[-10 10])
%绘制垂直于x轴的垂线
line([-3 -3],[0 max(zcords)],'Color','k');
line([-2 -2],[0 max(zcords)],'Color','k');
line([-1 -1],[0 max(zcords)],'Color','k');
line([0 0],[0 max(zcords)],'Color','k');
line([1 1],[0 max(zcords)],'Color','k');
line([2 2],[0 max(zcords)],'Color','k');
line([3 3],[0 max(zcords)],'Color','k');

subplot(2,1,2)
[~,II1]=min(abs(xcords(:)-(-3)));
[~,II2]=min(abs(xcords(:)-(-2)));
[~,II3]=min(abs(xcords(:)-(-1)));
[~,II4]=min(abs(xcords(:)-(0)));
[~,II5]=min(abs(xcords(:)-(1)));
[~,II6]=min(abs(xcords(:)-(2)));
[~,II7]=min(abs(xcords(:)-(3)));

uu(1,:) = s_sm(II1,:);
uu(2,:) = s_sm(II2,:)+10;
uu(3,:) = s_sm(II3,:)+20;
uu(4,:) = s_sm(II4,:)+30;
uu(5,:) = s_sm(II5,:)+40;
uu(6,:) = s_sm(II6,:)+50;
uu(7,:) = s_sm(II7,:)+60;
plot(uu','k')
title('(b)')
xlabel('z(km)')
ylabel('U(m/s)')
view(90,-90)
line([0 200],[0 0],'Color','k');
line([0 200],[10 10],'Color','k');
line([0 200],[20 20],'Color','k');
line([0 200],[30 30],'Color','k');
line([0 200],[40 40],'Color','k');
line([0 200],[50 50],'Color','k');
line([0 200],[60 60],'Color','k');
set(gca,'YTick',[])
set(gca,'XTick',[0 200])
set(gca,'XTicklabel',{'0','16'})























