% Add meridional modes
%% That is, the MJO response during the build-up, mature, and decay phases 
%% of synoptic-scale convection at a fixed location (warm pool);
%% This process mainly appears as the variation of the maximum heating center (z),
%% i.e., controlling the period by changing the fractional coverage of layer clouds;
%% Build-up phase of organized synoptic activity: a(t) goes from 0 to 1, 
%% reaching maximum in ~20 days
%% Decay phase: a(t) goes from 1 to 0
%% Symmetric about the equator; for the first baroclinic meridional mode, 
%% all odd terms of the meridional spectral coefficients are zero

clear all; close all; clc

%*******************************************************

a        = 1/3;       % Fractional coverage of layer clouds, can later vary in t and X
e        = 0.125;     % Froude number
L        = 5000;      % Unit: km; warm pool radius
y0       = 0;         % Distance from equator; choose symmetric here (y0 = 0)
HT       = 16;        % Troposphere height
ro       = 0.65;      % Typical wave packet length scale
x0       = 0.5;       % Layer cloud precipitation lag distance
x_scale  = 1500;      % Unit: km
y_scale  = 1500;      % Unit: km
z_scale  = HT/pi;     % Unit: km
X_scale  = 12000;     % Unit: km
t_scale  = 3;         % Unit: days; all times later should be converted to days, then nondimensionalized
u_scale  = 6.25;      % Unit: m/s; planetary-scale and synoptic-scale same
v_scale  = 6.25;      % Unit: m/s; note: planetary-scale v is too small, so use synoptic-scale magnitude
w_scale  = 2.5;       % Unit: cm/s; same, synoptic-scale
p_scale  = 312;       % Unit: (m/s)^2; planetary-scale and synoptic-scale same
d        = 0.18;      % Controversial; chosen for cooling-related effect, dissipation slower than drag
L        = L/x_scale; % Still using synoptic-scale characteristic value for convenience
c        = 5/u_scale; % Envelope propagation speed

%*******************************************************

N          = 31;                        % e.g., for meridional expansion 0-4, v=5, r=6; matrix needs 6+1=7
                                        % Expanding to 20th order requires 20-1+3=22
dx         = 100/x_scale;               % Spatial resolution: 150 km
dz         = 0.1/z_scale;               % Vertical resolution: 0.5 km
dt         = 8.3*0.000694/t_scale;      % Time step: 16.6 min (or 8.3 min, 4.15 min)
xgrids1    = -10000/x_scale; 
xgrids2    =  10000/x_scale;
ygrids1    = -3000/y_scale; 
ygrids2    =  3000/y_scale;
zgrids1    = 0./z_scale; 
zgrids2    = HT/z_scale;

ngrids_x   = floor((xgrids2-xgrids1)/dx)+1;
ngrids_y   = floor((ygrids2-ygrids1)/dx)+1;
ngrids_z   = floor((zgrids2-zgrids1)/dz)+1;

ngrids_t   = 7200;                      % Time steps: for 16.6 min: integrate 3600 steps (~41.5 days)
                                        % 8.3 min: 7200 steps (~41.5 days); 4.15 min: 3600*4 steps
n_time     = 3;                         % Number of custom integration periods, each ~41.5 days

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
%                     1. Orthogonal Basis
%=======================================================

syms y m

% Define the standard Hermite function (orthonormal basis)
B1 = sqrt(1./(2.^m .* factorial(m) .* sqrt(pi))) .* hermiteH(m, y) .* exp(-y.^2./2);

% Define a scaled Hermite function for a different meridional mode
B3 = sqrt(sqrt(3)./(2.^m .* factorial(m) .* sqrt(pi))) .* hermiteH(m, sqrt(3).*y) .* exp(-3.*y.^2./2);

% Convert symbolic expressions to MATLAB functions
B1 = matlabFunction(B1);
B3 = matlabFunction(B3);

% Test the orthogonality
% These are standard orthonormal vectors.
% Note: When multiplying functions, always include the variable, e.g., B1(m,y)*B3(m,y), 
% then convert using matlabFunction(ans) if needed.

test(1) = double(int(B1(1,y) * B1(2,y), y, -inf, inf));  % Should be ~0
test(2) = double(int(B1(1,y) * B1(1,y), y, -inf, inf));  % Should be ~1
test(3) = double(int(B3(1,y) * B3(2,y), y, -inf, inf));  % Should be ~0
test(4) = double(int(B3(1,y) * B3(1,y), y, -inf, inf));  % Should be ~1



%========================================================
%               2. Meridional Spectral Coefficients for fr and fl
%========================================================

syms x z t

% Envelope function along x (zonal direction)
F = cos(pi.*(x - c.*t) ./ (2.*L)) .* (x >= -L & x <= L) + 0.*(x > L | x < -L); 

% Meridional structure (Gaussian)
H  = sqrt(10) .* exp(-(y + y0).^2);  
Hy = diff(H);  % derivative of H w.r.t y

% Scaling factor for convective heating
k = 3/4 .* F.^2 .* a .* ro .* sin(x0 ./ ro);

% Zonal & temporal components of heating for first and third baroclinic modes
FU1 = k .* (2 .* H.^2 + y .* H .* Hy);         % zonal flux for first mode
FT1 = k .* (5 .* y.^2 .* H.^2 + 4 .* y.^3 .* Hy .* H);  % temporal flux for first mode

FU3 = -k .* (2 .* H.^2 + y .* H .* Hy);        % zonal flux for third mode
FT3 = (k ./ 9) .* (15 .* y.^2 .* H.^2 + 4 .* y.^3 .* Hy .* H); % temporal flux for third mode

% Construct fr and fl (right- and left-moving components)
fr1 = matlabFunction(0.5 .* (FT1 + FU1));  % first baroclinic right-moving
fl1 = matlabFunction(0.5 .* (FT1 - FU1));  % first baroclinic left-moving
fr3 = matlabFunction(0.5 .* (FT3 + FU3));  % third baroclinic right-moving
fl3 = matlabFunction(0.5 .* (FT3 - FU3));  % third baroclinic left-moving

%% Meridional decomposition in terms of Hermite modes (m)
%  First baroclinic fr coefficient:
% fr1_m = sqrt(1./(2.^m .* factorial(m) .* sqrt(pi))) .* ...
%         int(fr1(t,x,y) .* hermiteH(m, y) .* exp(-y.^2./2), y, -inf, inf);
%
%  First baroclinic fl coefficient:
% fl1_m = sqrt(1./(2.^m .* factorial(m) .* sqrt(pi))) .* ...
%         int(fl1(t,x,y) .* hermiteH(m, y) .* exp(-y.^2./2), y, -inf, inf);
%
%  Third baroclinic fr coefficient:
% fr3_m = sqrt(sqrt(3)./(2.^m .* factorial(m) .* sqrt(pi))) .* ...
%         int(fr3(t,x,y) .* hermiteH(m, sqrt(3).*y) .* exp(-3.*y.^2./2), y, -inf, inf);
%
%  Third baroclinic fl coefficient:
% fl3_m = sqrt(sqrt(3)./(2.^m .* factorial(m) .* sqrt(pi))) .* ...
%         int(fl3(t,x,y) .* hermiteH(m, sqrt(3).*y) .* exp(-3.*y.^2./2), y, -inf, inf);

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
%                           1st Baroclinic Mode
%
%=============================================================================
%=============================================================================


%==============================
%        Numerical Integration
%==============================

% Define coefficients for finite difference scheme along zonal direction
% Note: indices: lm, vm+1, rm+2
c1 = matlabFunction((1./(2.*m + 3)) ./ e);
c2 = matlabFunction(sqrt((m + 1) .* (m + 2)) ./ (2.*m + 3));
c3 = matlabFunction((m + 2) ./ (2.*m + 3)); 

d0 = (1./sqrt(2)) .* (sqrt(m + 2) + (m + 1) ./ sqrt(m + 2));
d1 = (2 .* sqrt((m + 1) ./ (m + 2))) ./ d0;
d1 = matlabFunction(d1 ./ e);
d2 = matlabFunction((sqrt((m + 1) ./ (m + 2))) ./ d0);
d0 = matlabFunction(1 ./ d0);

%-------------------------------
% Initial mode: l0, v1, r2 for q=0
%-------------------------------
q    = 0;
n    = 0;

ll   = zeros(ngrids_x, 2); % left-moving component
vv   = zeros(ngrids_x, 1); % meridional velocity
lll  = [];
vvv  = [];
rrr  = [];

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        % Finite difference update for left-moving component (ll)
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) ...
             - d .* ll(ii,1) ...
             + c2(q) .* fr1_2(n .* dt, xcords(ii)) ...
             + c3(q) .* fl1_0(n .* dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);
    end

    % Smooth daily output (adjust depending on time step)
    n_mod = mod(n, 87*2);  % For 16.6 min timestep: 87 steps per day; adjust for other timesteps
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(:,1) = smoothdata(ll(:,1));
        ll(:,2) = smoothdata(ll(:,2));
    end

    % Save data every day (or fixed interval)
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) ...
                       + d2(q) .* fl1_0(n .* dt, xcords(ii)) ...
                       - d0(q) .* fr1_2(n .* dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l01.txt', 'a'); 
        fid2 = fopen('C:\Users\pc\Desktop\data\v11.txt', 'a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r21.txt', 'a');
        fprintf(fid1, '%e\n', ll(:,1));  
        fprintf(fid2, '%e\n', vv(:,1));
        fprintf(fid3, '%e\n', sqrt((q + 1) ./ (q + 2)) .* ll(:,1));
        fclose(fid1);
        fclose(fid2);
        fclose(fid3);
    end

    n = n + 1;
end

%-------------------------------
% Reload saved files into matrices
%-------------------------------
lll = load('C:\Users\pc\Desktop\data\l01.txt');
vvv = load('C:\Users\pc\Desktop\data\v11.txt');
rrr = load('C:\Users\pc\Desktop\data\r21.txt');

% Store left, right, and meridional components for time–longitude plots
% This is done only the first time
l = zeros(ngrids_x, length(lll)/ngrids_x, N);
v = zeros(ngrids_x, length(lll)/ngrids_x, N);
r = zeros(ngrids_x, length(lll)/ngrids_x, N);

l(:,:,1) = reshape(lll, [ngrids_x, length(lll)/ngrids_x]);
v(:,:,2) = reshape(vvv, [ngrids_x, length(lll)/ngrids_x]);
r(:,:,3) = reshape(rrr, [ngrids_x, length(lll)/ngrids_x]);



%===============================
% 1.2) Left, meridional, right components: q = 1
%===============================
q    = 1; % mode index
n    = 0;

ll   = zeros(ngrids_x, 2); % left-moving component
vv   = zeros(ngrids_x, 1); % meridional velocity
lll  = [];
vvv  = [];
rrr  = [];

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        % Central difference update for the left-moving component
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) - d .* ll(ii,1);
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);
    end

    % Smooth daily output
    n_mod = mod(n, 87*2);  % adjust for time step (16.6 min: 87*2)
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(:,1) = smoothdata(ll(:,1));
        ll(:,2) = smoothdata(ll(:,2));
    end

    % Save output to text files every day (or interval)
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l11.txt', 'a'); 
        fid2 = fopen('C:\Users\pc\Desktop\data\v21.txt', 'a'); 
        fid3 = fopen('C:\Users\pc\Desktop\data\r31.txt', 'a'); 
        fprintf(fid1, '%e\n', ll(:,1));  
        fprintf(fid2, '%e\n', vv(:,1));
        fprintf(fid3, '%e\n', sqrt((q + 1) ./ (q + 2)) .* ll(:,1));
        fclose(fid1);
        fclose(fid2);
        fclose(fid3);
    end

    n = n + 1;
end

%-------------------------------
% Reload saved data into matrices
%-------------------------------
lll = load('C:\Users\pc\Desktop\data\l11.txt');
vvv = load('C:\Users\pc\Desktop\data\v21.txt');
rrr = load('C:\Users\pc\Desktop\data\r31.txt');

% Reshape into 3D matrices: [x, time, mode]
l(:,:,2) = reshape(lll, [ngrids_x, length(lll)/ngrids_x]);
v(:,:,3) = reshape(vvv, [ngrids_x, length(lll)/ngrids_x]);
r(:,:,4) = reshape(rrr, [ngrids_x, length(lll)/ngrids_x]);



%===============================================================================
% 1.3) Left, meridional, right components: q = 2 (3rd baroclinic mode)
%===============================================================================
q    = 2; % mode index
n    = 0;

ll   = zeros(ngrids_x, 2); % left-moving component
vv   = zeros(ngrids_x, 1); % meridional velocity
lll  = [];
vvv  = [];
rrr  = [];

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        % Central difference update with forcing terms
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) ...
             - d .* ll(ii,1) ...
             + c2(q) .* fr1_4(n*dt, xcords(ii)) ...
             + c3(q) .* fl1_2(n*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);
    end

    % Smooth daily output
    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(:,1) = smoothdata(ll(:,1));
        ll(:,2) = smoothdata(ll(:,2));
    end

    % Save output to files
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) ...
                        + d2(q) .* fl1_2(n*dt, xcords(ii)) ...
                        - d0(q) .* fr1_4(n*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l21.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v31.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r41.txt','a');
        fprintf(fid1, '%e\n', ll(:,1));
        fprintf(fid2, '%e\n', vv(:,1));
        fprintf(fid3, '%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end

    n = n + 1;
end

% Load saved data and reshape to [x, time, mode]
lll = load('C:\Users\pc\Desktop\data\l21.txt');
vvv = load('C:\Users\pc\Desktop\data\v31.txt');
rrr = load('C:\Users\pc\Desktop\data\r41.txt');

l(:,:,3) = reshape(lll, [ngrids_x, length(lll)/ngrids_x]);
v(:,:,4) = reshape(vvv, [ngrids_x, length(lll)/ngrids_x]);
r(:,:,5) = reshape(rrr, [ngrids_x, length(lll)/ngrids_x]);


%===============================================================================
% 1.4) Left, meridional, right components: q = 3 (4th baroclinic mode)
%===============================================================================
q    = 3; % mode index
n    = 0;

ll   = zeros(ngrids_x, 2);
vv   = zeros(ngrids_x, 1);
lll  = [];
vvv  = [];
rrr  = [];

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        % Central difference update, no forcing terms for this mode
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) - d .* ll(ii,1);
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);
    end

    % Smooth daily output
    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(:,1) = smoothdata(ll(:,1));
        ll(:,2) = smoothdata(ll(:,2));
    end

    % Save output to files
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l31.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v41.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r51.txt','a');
        fprintf(fid1, '%e\n', ll(:,1));
        fprintf(fid2, '%e\n', vv(:,1));
        fprintf(fid3, '%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end

    n = n + 1;
end

% Load saved data and reshape to [x, time, mode]
lll = load('C:\Users\pc\Desktop\data\l31.txt');
vvv = load('C:\Users\pc\Desktop\data\v41.txt');
rrr = load('C:\Users\pc\Desktop\data\r51.txt');

l(:,:,4) = reshape(lll, [ngrids_x, length(lll)/ngrids_x]);
v(:,:,5) = reshape(vvv, [ngrids_x, length(lll)/ngrids_x]);
r(:,:,6) = reshape(rrr, [ngrids_x, length(lll)/ngrids_x]);



%===============================================================================
% 1.5) Left, meridional, right components: q = 4 (5th baroclinic mode)
%===============================================================================
q = 4; % mode index
n = 0;
ll = zeros(ngrids_x, 2);
vv = zeros(ngrids_x, 1);
lll = [];
vvv = [];
rrr = [];

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) ...
             - d .* ll(ii,1) ...
             + c2(q) .* fr1_6(n*dt, xcords(ii)) ...
             + c3(q) .* fl1_4(n*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);
    end

    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(:,1) = smoothdata(ll(:,1));
        ll(:,2) = smoothdata(ll(:,2));
    end

    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) ...
                        + d2(q) .* fl1_4(n*dt, xcords(ii)) ...
                        - d0(q) .* fr1_6(n*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l41.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v51.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r61.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end

    n = n + 1;
end

lll = load('C:\Users\pc\Desktop\data\l41.txt');
vvv = load('C:\Users\pc\Desktop\data\v51.txt');
rrr = load('C:\Users\pc\Desktop\data\r61.txt');

l(:,:,5) = reshape(lll, [ngrids_x, length(lll)/ngrids_x]);
v(:,:,6) = reshape(vvv, [ngrids_x, length(lll)/ngrids_x]);
r(:,:,7) = reshape(rrr, [ngrids_x, length(lll)/ngrids_x]);


%===============================================================================
% 1.6) Left, meridional, right components: q = 5 (6th baroclinic mode)
%===============================================================================
q = 5;
n = 0;
ll = zeros(ngrids_x,2);
vv = zeros(ngrids_x,1);
lll = [];
vvv = [];
rrr = [];

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        % Advection + damping only, no forcing terms
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) - d .* ll(ii,1);
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);
    end

    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(:,1) = smoothdata(ll(:,1));
        ll(:,2) = smoothdata(ll(:,2));
    end

    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l51.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v61.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r71.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

lll = load('C:\Users\pc\Desktop\data\l51.txt');
vvv = load('C:\Users\pc\Desktop\data\v61.txt');
rrr = load('C:\Users\pc\Desktop\data\r71.txt');

l(:,:,6) = reshape(lll, [ngrids_x, length(lll)/ngrids_x]);
v(:,:,7) = reshape(vvv, [ngrids_x, length(lll)/ngrids_x]);
r(:,:,8) = reshape(rrr, [ngrids_x, length(lll)/ngrids_x]);


%===============================================================================
% 1.7) Left, meridional, right components: q = 6 (7th baroclinic mode)
%===============================================================================
q = 6;
n = 0;
ll = zeros(ngrids_x,2);
vv = zeros(ngrids_x,1);
lll = [];
vvv = [];
rrr = [];

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) ...
             - d .* ll(ii,1) ...
             + c2(q) .* fr1_8(n*dt, xcords(ii)) ...
             + c3(q) .* fl1_6(n*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);
    end

    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(:,1) = smoothdata(ll(:,1));
        ll(:,2) = smoothdata(ll(:,2));
    end

    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) ...
                        + d2(q) .* fl1_6(n*dt, xcords(ii)) ...
                        - d0(q) .* fr1_8(n*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l61.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v71.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r81.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

lll = load('C:\Users\pc\Desktop\data\l61.txt');
vvv = load('C:\Users\pc\Desktop\data\v71.txt');
rrr = load('C:\Users\pc\Desktop\data\r81.txt');

l(:,:,7) = reshape(lll, [ngrids_x, length(lll)/ngrids_x]);
v(:,:,8) = reshape(vvv, [ngrids_x, length(lll)/ngrids_x]);
r(:,:,9) = reshape(rrr, [ngrids_x, length(lll)/ngrids_x]);


%===============================================================================
% 1.8) Left, meridional, right components: q = 7 (8th baroclinic mode)
%===============================================================================
q = 7;
n = 0;
ll = zeros(ngrids_x,2);  % left component
vv = zeros(ngrids_x,1);  % meridional component
lll = [];
vvv = [];
rrr = [];

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        % Advection + damping only (no convective forcing)
        c4 = c1(q) .* ((ll(ii+1,1)-ll(ii-1,1)) ./ (2 .* dx)) - d .* ll(ii,1);
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);
    end
    
    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(:,1) = smoothdata(ll(:,1));
        ll(:,2) = smoothdata(ll(:,2));
    end
    
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1)-ll(ii-1,1)) ./ (2 .* dx));
        end
        % Append to output files
        fid1 = fopen('C:\Users\pc\Desktop\data\l71.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v81.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r91.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    
    n = n + 1;
end

% Load files as matrices
lll = load('C:\Users\pc\Desktop\data\l71.txt');
vvv = load('C:\Users\pc\Desktop\data\v81.txt');
rrr = load('C:\Users\pc\Desktop\data\r91.txt');

l(:,:,8) = reshape(lll,[ngrids_x, length(lll)/ngrids_x]);
v(:,:,9) = reshape(vvv,[ngrids_x, length(lll)/ngrids_x]);
r(:,:,10) = reshape(rrr,[ngrids_x, length(lll)/ngrids_x]);


%===============================================================================
% 1.9) Left, meridional, right components: q = 8 (9th baroclinic mode)
% Convective forcing included
%===============================================================================
q = 8;
n = 0;
ll = zeros(ngrids_x,2); vv = zeros(ngrids_x,1); lll = []; vvv = []; rrr = [];

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        c4 = c1(q) .* ((ll(ii+1,1)-ll(ii-1,1))./(2 .* dx)) - d .* ll(ii,1) ...
             + c2(q) .* fr1_10(n*dt, xcords(ii)) ...
             + c3(q) .* fl1_8(n*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4*dt;
        ll(ii,1) = ll(ii,2);
    end
    
    n_mod = mod(n, 87*2); k_mod = mod(n,87);
    if k_mod==0, ll(:,1)=smoothdata(ll(:,1)); ll(:,2)=smoothdata(ll(:,2)); end
    
    if n_mod==0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q)*((ll(ii+1,1)-ll(ii-1,1))/(2*dx)) ...
                        + d2(q)*fl1_8(n*dt,xcords(ii)) ...
                        - d0(q)*fr1_10(n*dt,xcords(ii));
        end
        fid1=fopen('C:\Users\pc\Desktop\data\l81.txt','a');
        fid2=fopen('C:\Users\pc\Desktop\data\v91.txt','a');
        fid3=fopen('C:\Users\pc\Desktop\data\r101.txt','a');
        fprintf(fid1,'%e\n',ll(:,1)); fprintf(fid2,'%e\n',vv(:,1));
        fprintf(fid3,'%e\n',sqrt((q+1)/(q+2))*ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    
    n = n+1;
end

lll = load('C:\Users\pc\Desktop\data\l81.txt');
vvv = load('C:\Users\pc\Desktop\data\v91.txt');
rrr = load('C:\Users\pc\Desktop\data\r101.txt');

l(:,:,9) = reshape(lll,[ngrids_x, length(lll)/ngrids_x]);
v(:,:,10)= reshape(vvv,[ngrids_x, length(lll)/ngrids_x]);
r(:,:,11)= reshape(rrr,[ngrids_x, length(lll)/ngrids_x]);


%===============================================================================
% 1.11) Left, meridional, right components: q = 10 (11th baroclinic mode)
% Convective forcing included
%===============================================================================
q = 10;
n = 0;
ll = zeros(ngrids_x,2); vv = zeros(ngrids_x,1); lll = []; vvv = []; rrr = [];

while n <= ngrids_t*n_time
    for ii = 2:ngrids_x-1
        c4 = c1(q) * ((ll(ii+1,1)-ll(ii-1,1))/(2*dx)) - d*ll(ii,1) ...
             + c2(q)*fr1_12(n*dt,xcords(ii)) ...
             + c3(q)*fl1_10(n*dt,xcords(ii));
        ll(ii,2) = ll(ii,1) + c4*dt; ll(ii,1) = ll(ii,2);
    end
    
    n_mod = mod(n,87*2); k_mod=mod(n,87);
    if k_mod==0, ll(:,1)=smoothdata(ll(:,1)); ll(:,2)=smoothdata(ll(:,2)); end
    
    if n_mod==0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q)*((ll(ii+1,1)-ll(ii-1,1))/(2*dx)) ...
                        + d2(q)*fl1_10(n*dt,xcords(ii)) ...
                        - d0(q)*fr1_12(n*dt,xcords(ii));
        end
        fid1=fopen('C:\Users\pc\Desktop\data\l101.txt','a');
        fid2=fopen('C:\Users\pc\Desktop\data\v111.txt','a');
        fid3=fopen('C:\Users\pc\Desktop\data\r121.txt','a');
        fprintf(fid1,'%e\n',ll(:,1)); fprintf(fid2,'%e\n',vv(:,1));
        fprintf(fid3,'%e\n',sqrt((q+1)/(q+2))*ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    
    n = n+1;
end

lll = load('C:\Users\pc\Desktop\data\l101.txt');
vvv = load('C:\Users\pc\Desktop\data\v111.txt');
rrr = load('C:\Users\pc\Desktop\data\r121.txt');

l(:,:,11) = reshape(lll,[ngrids_x, length(lll)/ngrids_x]);
v(:,:,12)= reshape(vvv,[ngrids_x, length(lll)/ngrids_x]);
r(:,:,13)= reshape(rrr,[ngrids_x, length(lll)/ngrids_x]);


% 1.13) l12,v13,r14:q=12   !!!!!!!!!
q    = 12; % change
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n <= ngrids_t * n_time
	for ii = 2:ngrids_x-1
		c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) + ...
				c2(q) .* fr1_14(n.*dt, xcords(ii)) + ... % change
					c3(q) .* fl1_12(n.*dt, xcords(ii)); % change
		ll(ii,2) = ll(ii,1) + c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	% Output one file per day
	n_mod = mod(n, 87*2);
	k_mod = mod(n, 87);
	if (k_mod == 0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if (n_mod == 0)
		for ii = 2:ngrids_x-1
			vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) + ...
						d2(q) .* fl1_12(n.*dt, xcords(ii)) - ... % change
							d0(q) .* fr1_14(n.*dt, xcords(ii)); % change
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l121.txt','a'); % change
		fid2 = fopen('C:\Users\pc\Desktop\data\v131.txt','a'); % change
		fid3 = fopen('C:\Users\pc\Desktop\data\r141.txt','a'); % change
		fprintf(fid1,'%e\n', ll(:,1));  
		fprintf(fid2,'%e\n', vv(:,1));
		fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n + 1;
end
% Reload into matrix format
lll = load('C:\Users\pc\Desktop\data\l121.txt'); % change
vvv = load('C:\Users\pc\Desktop\data\v131.txt'); % change
rrr = load('C:\Users\pc\Desktop\data\r141.txt'); % change
l(:,:,13) = reshape(lll, [ngrids_x length(lll)/ngrids_x]); % change
v(:,:,14) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]); % change
r(:,:,15) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]); % change

% 1.15) l14,v15,r16:q=14   !!!!!!!!!
q    = 14; % change
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n <= ngrids_t * n_time
	for ii = 2:ngrids_x-1
		c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) + ...
				c2(q) .* fr1_16(n.*dt, xcords(ii)) + ... % change
					c3(q) .* fl1_14(n.*dt, xcords(ii)); % change
		ll(ii,2) = ll(ii,1) + c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	% Output one file per day
	n_mod = mod(n, 87*2);
	k_mod = mod(n, 87);
	if (k_mod == 0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if (n_mod == 0)
		for ii = 2:ngrids_x-1
			vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) + ...
						d2(q) .* fl1_14(n.*dt, xcords(ii)) - ... % change
							d0(q) .* fr1_16(n.*dt, xcords(ii)); % change
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l141.txt','a'); % change
		fid2 = fopen('C:\Users\pc\Desktop\data\v151.txt','a'); % change
		fid3 = fopen('C:\Users\pc\Desktop\data\r161.txt','a'); % change
		fprintf(fid1,'%e\n', ll(:,1));  
		fprintf(fid2,'%e\n', vv(:,1));
		fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n + 1;
end
% Reload into matrix format
lll = load('C:\Users\pc\Desktop\data\l141.txt'); % change
vvv = load('C:\Users\pc\Desktop\data\v151.txt'); % change
rrr = load('C:\Users\pc\Desktop\data\r161.txt'); % change
l(:,:,15) = reshape(lll, [ngrids_x length(lll)/ngrids_x]); % change
v(:,:,16) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]); % change
r(:,:,17) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]); % change

% 1.17) l16,v17,r18:q=16   !!!!!!!!!
q    = 16; % change
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n <= ngrids_t * n_time
	for ii = 2:ngrids_x-1
		c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) + ...
				c2(q) .* fr1_18(n.*dt, xcords(ii)) + ... % change
					c3(q) .* fl1_16(n.*dt, xcords(ii)); % change
		ll(ii,2) = ll(ii,1) + c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	% Output one file per day
	n_mod = mod(n, 87*2);
	k_mod = mod(n, 87);
	if (k_mod == 0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if (n_mod == 0)
		for ii = 2:ngrids_x-1
			vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) + ...
						d2(q) .* fl1_16(n.*dt, xcords(ii)) - ... % change
							d0(q) .* fr1_18(n.*dt, xcords(ii)); % change
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l161.txt','a'); % change
		fid2 = fopen('C:\Users\pc\Desktop\data\v171.txt','a'); % change
		fid3 = fopen('C:\Users\pc\Desktop\data\r181.txt','a'); % change
		fprintf(fid1,'%e\n', ll(:,1));  
		fprintf(fid2,'%e\n', vv(:,1));
		fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n + 1;
end
% Reload into matrix format
lll = load('C:\Users\pc\Desktop\data\l161.txt'); % change
vvv = load('C:\Users\pc\Desktop\data\v171.txt'); % change
rrr = load('C:\Users\pc\Desktop\data\r181.txt'); % change
l(:,:,17) = reshape(lll, [ngrids_x length(lll)/ngrids_x]); % change
v(:,:,18) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]); % change
r(:,:,19) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]); % change

% 1.19) l18,v19,r20:q=18   !!!!!!!!!
q    = 18; % change
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n <= ngrids_t * n_time
	for ii = 2:ngrids_x-1
		c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) + ...
				c2(q) .* fr1_20(n.*dt, xcords(ii)) + ... % change
					c3(q) .* fl1_18(n.*dt, xcords(ii)); % change
		ll(ii,2) = ll(ii,1) + c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	% Output one file per day
	n_mod = mod(n, 87*2);
	k_mod = mod(n, 87);
	if (k_mod == 0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if (n_mod == 0)
		for ii = 2:ngrids_x-1
			vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) + ...
						d2(q) .* fl1_18(n.*dt, xcords(ii)) - ... % change
							d0(q) .* fr1_20(n.*dt, xcords(ii)); % change
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l181.txt','a'); % change
		fid2 = fopen('C:\Users\pc\Desktop\data\v191.txt','a'); % change
		fid3 = fopen('C:\Users\pc\Desktop\data\r201.txt','a'); % change
		fprintf(fid1,'%e\n', ll(:,1));  
		fprintf(fid2,'%e\n', vv(:,1));
		fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n + 1;
end
% Reload into matrix format
lll = load('C:\Users\pc\Desktop\data\l181.txt'); % change
vvv = load('C:\Users\pc\Desktop\data\v191.txt'); % change
rrr = load('C:\Users\pc\Desktop\data\r201.txt'); % change
l(:,:,19) = reshape(lll, [ngrids_x length(lll)/ngrids_x]); % change
v(:,:,20) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]); % change
r(:,:,21) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]); % change

% 1.21) l20,v21,r22:q=20   !!!!!!!!!
q    = 20; % change
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n <= ngrids_t * n_time
	for ii = 2:ngrids_x-1
		c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) + ...
				c2(q) .* fr1_22(n.*dt, xcords(ii)) + ... % change
					c3(q) .* fl1_20(n.*dt, xcords(ii)); % change
		ll(ii,2) = ll(ii,1) + c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	% Output one file per day
	n_mod = mod(n, 87*2); 
	k_mod = mod(n, 87);
	if (k_mod == 0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if (n_mod == 0)
		for ii = 2:ngrids_x-1
			vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) + ...
						d2(q) .* fl1_20(n.*dt, xcords(ii)) - ... % change
							d0(q) .* fr1_22(n.*dt, xcords(ii)); % change
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l201.txt','a'); % change
		fid2 = fopen('C:\Users\pc\Desktop\data\v211.txt','a'); % change
		fid3 = fopen('C:\Users\pc\Desktop\data\r221.txt','a'); % change
		fprintf(fid1,'%e\n', ll(:,1));  
		fprintf(fid2,'%e\n', vv(:,1));
		fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n + 1;
end
% Reload into matrix format
lll = load('C:\Users\pc\Desktop\data\l201.txt'); % change
vvv = load('C:\Users\pc\Desktop\data\v211.txt'); % change
rrr = load('C:\Users\pc\Desktop\data\r221.txt'); % change
l(:,:,21) = reshape(lll, [ngrids_x length(lll)/ngrids_x]); % change
v(:,:,22) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]); % change
r(:,:,23) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]); % change

% 1.23) l22,v23,r24:q=22   !!!!!!!!!
q    = 22; % change
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n <= ngrids_t * n_time
	for ii = 2:ngrids_x-1
		c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) + ...
				c2(q) .* fr1_24(n.*dt, xcords(ii)) + ... % change
					c3(q) .* fl1_22(n.*dt, xcords(ii)); % change
		ll(ii,2) = ll(ii,1) + c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	% Output one file per day
	n_mod = mod(n, 87*2); 
	k_mod = mod(n, 87);
	if (k_mod == 0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if (n_mod == 0)
		for ii = 2:ngrids_x-1
			vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) + ...
						d2(q) .* fl1_22(n.*dt, xcords(ii)) - ... % change
							d0(q) .* fr1_24(n.*dt, xcords(ii)); % change
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l221.txt','a'); % change
		fid2 = fopen('C:\Users\pc\Desktop\data\v231.txt','a'); % change
		fid3 = fopen('C:\Users\pc\Desktop\data\r241.txt','a'); % change
		fprintf(fid1,'%e\n', ll(:,1));  
		fprintf(fid2,'%e\n', vv(:,1));
		fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n + 1;
end
% Reload into matrix format
lll = load('C:\Users\pc\Desktop\data\l221.txt'); % change
vvv = load('C:\Users\pc\Desktop\data\v231.txt'); % change
rrr = load('C:\Users\pc\Desktop\data\r241.txt'); % change
l(:,:,23) = reshape(lll, [ngrids_x length(lll)/ngrids_x]); % change
v(:,:,24) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]); % change
r(:,:,25) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]); % change

% 1.25) l24,v25,r26:q=24   !!!!!!!!!
q    = 24; % change
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n <= ngrids_t * n_time
	for ii = 2:ngrids_x-1
		c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) + ...
				c2(q) .* fr1_26(n.*dt, xcords(ii)) + ... % change
					c3(q) .* fl1_24(n.*dt, xcords(ii)); % change
		ll(ii,2) = ll(ii,1) + c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	% Output one file per day
	n_mod = mod(n, 87*2); 
	k_mod = mod(n, 87);
	if (k_mod == 0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if (n_mod == 0)
		for ii = 2:ngrids_x-1
			vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) + ...
						d2(q) .* fl1_24(n.*dt, xcords(ii)) - ... % change
							d0(q) .* fr1_26(n.*dt, xcords(ii)); % change
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l241.txt','a'); % change
		fid2 = fopen('C:\Users\pc\Desktop\data\v251.txt','a'); % change
		fid3 = fopen('C:\Users\pc\Desktop\data\r261.txt','a'); % change
		fprintf(fid1,'%e\n', ll(:,1));  
		fprintf(fid2,'%e\n', vv(:,1));
		fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n + 1;
end
% Reload into matrix format
lll = load('C:\Users\pc\Desktop\data\l241.txt'); % change
vvv = load('C:\Users\pc\Desktop\data\v251.txt'); % change
rrr = load('C:\Users\pc\Desktop\data\r261.txt'); % change
l(:,:,25) = reshape(lll, [ngrids_x length(lll)/ngrids_x]); % change
v(:,:,26) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]); % change
r(:,:,27) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]); % change

% 1.27)
q    = 26; % change
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n <= ngrids_t * n_time
	for ii = 2:ngrids_x-1
		c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) + ...
				c2(q) .* fr1_28(n.*dt, xcords(ii)) + ... % change
					c3(q) .* fl1_26(n.*dt, xcords(ii)); % change
		ll(ii,2) = ll(ii,1) + c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	% Output one file per day
	n_mod = mod(n, 87*2); 
	k_mod = mod(n, 87);
	if (k_mod == 0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if (n_mod == 0)
		for ii = 2:ngrids_x-1
			vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) + ...
						d2(q) .* fl1_26(n.*dt, xcords(ii)) - ... % change
							d0(q) .* fr1_28(n.*dt, xcords(ii)); % change
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l261.txt','a'); % change
		fid2 = fopen('C:\Users\pc\Desktop\data\v271.txt','a'); % change
		fid3 = fopen('C:\Users\pc\Desktop\data\r281.txt','a'); % change
		fprintf(fid1,'%e\n', ll(:,1));  
		fprintf(fid2,'%e\n', vv(:,1));
		fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n + 1;
end
% Reload into matrix format
lll = load('C:\Users\pc\Desktop\data\l261.txt'); % change
vvv = load('C:\Users\pc\Desktop\data\v271.txt'); % change
rrr = load('C:\Users\pc\Desktop\data\r281.txt'); % change
l(:,:,27) = reshape(lll, [ngrids_x length(lll)/ngrids_x]); % change
v(:,:,28) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]); % change
r(:,:,29) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]); % change


% 1.29) 
q    = 28; % change
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];
while n <= ngrids_t * n_time
	for ii = 2:ngrids_x-1
		c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) + ...
				c2(q) .* fr1_30(n.*dt, xcords(ii)) + ... % change
					c3(q) .* fl1_28(n.*dt, xcords(ii)); % change
		ll(ii,2) = ll(ii,1) + c4.*dt;
		ll(ii,1) = ll(ii,2);
	end
	% Output one file per day
	n_mod = mod(n, 87*2); 
	k_mod = mod(n, 87);
	if (k_mod == 0)
		ll(ii,1) = smoothdata(ll(ii,1));
		ll(ii,2) = smoothdata(ll(ii,2));
	end
	if (n_mod == 0)
		for ii = 2:ngrids_x-1
			vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) + ...
						d2(q) .* fl1_28(n.*dt, xcords(ii)) - ... % change
							d0(q) .* fr1_30(n.*dt, xcords(ii)); % change
		end
		fid1 = fopen('C:\Users\pc\Desktop\data\l281.txt','a'); % change
		fid2 = fopen('C:\Users\pc\Desktop\data\v291.txt','a'); % change
		fid3 = fopen('C:\Users\pc\Desktop\data\r301.txt','a'); % change
		fprintf(fid1,'%e\n', ll(:,1));  
		fprintf(fid2,'%e\n', vv(:,1));
		fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
		fclose(fid1);
		fclose(fid2);
		fclose(fid3);
	end
	n = n + 1;
end
% Reload into matrix format
lll = load('C:\Users\pc\Desktop\data\l281.txt'); % change
vvv = load('C:\Users\pc\Desktop\data\v291.txt'); % change
rrr = load('C:\Users\pc\Desktop\data\r301.txt'); % change
l(:,:,29) = reshape(lll, [ngrids_x length(lll)/ngrids_x]); % change
v(:,:,30) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]); % change
r(:,:,31) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]); % change


% 2) v0  !!!!!!!!!
n  = 0;
nn = 1;
while n <= ngrids_t * n_time
	n_mod = mod(n, 87*2); 
	if (n_mod == 0)
		for ii = 2:ngrids_x-1
			% asymmetric: v(ii, nn, 1) = -sqrt(2) .* fr1_1(n.*dt, xcords(ii));
			v(ii, nn, 1) = 0;
		end
		nn = nn + 1;
	end
	n = n + 1;
end


% 3) r0  !!!!!!!!!
n    = 0;
h1   = 1 / e;
rr   = zeros(ngrids_x,2);
rrr  = [];
while n <= ngrids_t * n_time
	for ii = 2:ngrids_x-1
		c4 = -h1 .* ((rr(ii+1,1) - rr(ii-1,1)) ./ (2.*dx)) - d.*rr(ii,1) + ...
				fr1_0(n.*dt, xcords(ii)); % change
		rr(ii,2) = rr(ii,1) + c4.*dt;
		rr(ii,1) = rr(ii,2);
	end
	% Output one file per day
	n_mod = mod(n, 87*2); 
	k_mod = mod(n, 87);
	if (k_mod == 0)
		ll(ii,1) = smoothdata(rr(ii,1));
		ll(ii,2) = smoothdata(rr(ii,2));
	end
	if (n_mod == 0)
		fid3 = fopen('C:\Users\pc\Desktop\data\r01.txt','a'); % change
		fprintf(fid3,'%e\n', rr(:,1));
		fclose(fid3);
	end
	n = n + 1;
end
% Reload into matrix format
rrr = load('C:\Users\pc\Desktop\data\r01.txt'); % change
r(:,:,1) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]); % change



%                  Sum of Meridional Modes
%========================================================

s1 = [];
s2 = [];
s3 = [];
for ii = 1:size(l,2)
    s3 = 0;
    for jj = 1:N-2  % iterate over meridional modes
        s1 = B1(jj-1, ycords);  % evaluate the meridional basis function
        s2 = l(:, ii, jj);      % select the corresponding mode coefficient
        s3 = s3 + s1 .* s2;     % sum the contribution of each mode
    end
    l_N(:,:,ii) = s3;           % store the reconstructed l-field
end

s1 = [];
s2 = [];
s3 = [];
for ii = 1:size(l,2)
    s3 = 0;
    for jj = 1:N-2
        s1 = B1(jj-1, ycords);
        s2 = r(:, ii, jj);
        s3 = s3 + s1 .* s2;
    end
    r_N(:,:,ii) = s3;           % store the reconstructed r-field
end

s1 = [];
s2 = [];
s3 = [];
for ii = 1:size(l,2)
    s3 = 0;
    for jj = 1:N-2
        s1 = B1(jj-1, ycords);
        s2 = v(:, ii, jj);
        s3 = s3 + s1 .* s2;
    end
    v_N(:,:,ii) = s3;           % store the reconstructed v-field
end

% Compute combined fields
u_N1 = r_N - l_N;  % zonal wind component
p_N1 = r_N + l_N;  % pressure-like field
v_N1 = v_N;        % meridional wind component





%=============================================================================
%=============================================================================
%
%                           3rd Baroclinic Mode
%
%=============================================================================
%=============================================================================

%             Numerical integration: Replace B1 with B3; fr1, fl1 -> fr3, fl3
%========================================================

% 1) lm, vm+1, rm+2  !!!!!!!!!
c1 = matlabFunction((1./3 .* (1./(2.*m+3)))./e);        % coefficient c1 for q
c2 = matlabFunction(sqrt((m+1).*(m+2))./(2.*m+3));      % coefficient c2 for q
c3 = matlabFunction((m+2)./(2.*m+3));                  % coefficient c3 for q
d0 = (sqrt(6)./18).*(sqrt(m+2) + (m+1)./sqrt(m+2));    % auxiliary coefficient d0
d1 = ((2./3).*sqrt((m+1)./(m+2)))./d0;
d1 = matlabFunction(d1./e);                            % d1 as a function of m
d2 = matlabFunction((sqrt((m+1)./(m+2)))./d0);         % d2 as a function of m
d0 = matlabFunction(1./d0);                            % invert d0 for later use

% 1.1) l0, v1, r2 : q = 0  !!!!!!!!!
q    = 0;
n    = 0;
ll   = zeros(ngrids_x,2);   % l-field integration array
vv   = zeros(ngrids_x,1);   % v-field integration array
lll  = [];                   % to store l-field for output
vvv  = [];                   % to store v-field for output
rrr  = [];                   % to store r-field for output

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        % finite difference in x plus damping and forcing
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) ...
             - d .* ll(ii,1) ...
             + c2(q) .* fr3_2(n.*dt, xcords(ii)) ...   % forcing from fr3
             + c3(q) .* fl3_0(n.*dt, xcords(ii));      % forcing from fl3
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);   % update for next timestep
    end

    % Output to file every "day" (time interval)
    n_mod = mod(n, 87*2); 
    k_mod = mod(n, 87);
    if (k_mod == 0)
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if (n_mod == 0)
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) ...
                       + d2(q) .* fl3_0(n.*dt, xcords(ii)) ...
                       - d0(q) .* fr3_2(n.*dt, xcords(ii));
        end
        % Write daily output
        fid1 = fopen('C:\Users\pc\Desktop\data\l03.txt','a'); 
        fid2 = fopen('C:\Users\pc\Desktop\data\v13.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r23.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));  
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1);
        fclose(fid2);
        fclose(fid3);
    end
    n = n + 1;
end

% Reload as matrices
lll = load('C:\Users\pc\Desktop\data\l03.txt');
vvv = load('C:\Users\pc\Desktop\data\v13.txt');
rrr = load('C:\Users\pc\Desktop\data\r23.txt');

% Store l, r fields for plotting p and u as time-longitude plots
% This initialization is only needed the first time
l = zeros(ngrids_x, length(lll)/ngrids_x, N);
v = zeros(ngrids_x, length(lll)/ngrids_x, N);
r = zeros(ngrids_x, length(lll)/ngrids_x, N);
l(:,:,1) = reshape(lll, [ngrids_x length(lll)/ngrids_x]);
v(:,:,2) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]);
r(:,:,3) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]);



% 1.2) l1, v2, r3: q = 1   !!!!!!!!!
q    = 1; % change q
n    = 0;
ll   = zeros(ngrids_x,2);   % l-field temporary array
vv   = zeros(ngrids_x,1);   % v-field temporary array
lll  = [];                   % store l-field for output
vvv  = [];                   % store v-field for output
rrr  = [];                   % store r-field for output

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        % finite difference along x plus damping
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) - d .* ll(ii,1);
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);   % update for next timestep
    end

    % Write daily output to files
    n_mod = mod(n, 87*2); 
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l13.txt','a'); 
        fid2 = fopen('C:\Users\pc\Desktop\data\v23.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r33.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

% Reload as matrices
lll = load('C:\Users\pc\Desktop\data\l13.txt');
vvv = load('C:\Users\pc\Desktop\data\v23.txt');
rrr = load('C:\Users\pc\Desktop\data\r33.txt');
l(:,:,2) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);
v(:,:,3) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);
r(:,:,4) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);


% 1.3) l2, v3, r4: q = 2   !!!!!!!!!
q    = 2; % change q
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        % finite difference along x plus damping and forcing
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) - d .* ll(ii,1) ...
             + c2(q) .* fr3_4(n.*dt, xcords(ii)) ...  % forcing from fr3
             + c3(q) .* fl3_2(n.*dt, xcords(ii));     % forcing from fl3
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);
    end

    % Write daily output to files
    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) ...
                       + d2(q) .* fl3_2(n.*dt, xcords(ii)) ...
                       - d0(q) .* fr3_4(n.*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l23.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v33.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r43.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

% Reload as matrices
lll = load('C:\Users\pc\Desktop\data\l23.txt');
vvv = load('C:\Users\pc\Desktop\data\v33.txt');
rrr = load('C:\Users\pc\Desktop\data\r43.txt');
l(:,:,3) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);
v(:,:,4) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);
r(:,:,5) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);


% 1.4) l3, v4, r5: q = 3   !!!!!!!!!
q    = 3; % change q
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        % finite difference along x plus damping
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) - d .* ll(ii,1);
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);
    end

    % Write daily output to files
    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l33.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v43.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r53.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

% Reload as matrices
lll = load('C:\Users\pc\Desktop\data\l33.txt');
vvv = load('C:\Users\pc\Desktop\data\v43.txt');
rrr = load('C:\Users\pc\Desktop\data\r53.txt');
l(:,:,4) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);
v(:,:,5) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);
r(:,:,6) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);



% 1.5) l4, v5, r6: q = 4   !!!!!!!!!
q    = 4; % change q
n    = 0;
ll   = zeros(ngrids_x,2);   % temporary array for l-field
vv   = zeros(ngrids_x,1);   % temporary array for v-field
lll  = [];                   % store l-field for output
vvv  = [];                   % store v-field for output
rrr  = [];                   % store r-field for output

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        % finite difference along x + damping + 3rd baroclinic mode forcing
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) - d .* ll(ii,1) ...
             + c2(q) .* fr3_6(n.*dt, xcords(ii)) ...
             + c3(q) .* fl3_4(n.*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);   % update for next timestep
    end

    % Daily output to files
    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) ...
                       + d2(q) .* fl3_4(n.*dt, xcords(ii)) ...
                       - d0(q) .* fr3_6(n.*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l43.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v53.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r63.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

% Reload as matrices
lll = load('C:\Users\pc\Desktop\data\l43.txt');
vvv = load('C:\Users\pc\Desktop\data\v53.txt');
rrr = load('C:\Users\pc\Desktop\data\r63.txt');
l(:,:,5) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);
v(:,:,6) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);
r(:,:,7) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);


% 1.6) l5, v6, r7: q = 5   !!!!!!!!!
q    = 5; % change q
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        % finite difference along x + damping
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) - d .* ll(ii,1);
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);
    end

    % Daily output to files
    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l53.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v63.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r73.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

% Reload as matrices
lll = load('C:\Users\pc\Desktop\data\l53.txt');
vvv = load('C:\Users\pc\Desktop\data\v63.txt');
rrr = load('C:\Users\pc\Desktop\data\r73.txt');
l(:,:,6) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);
v(:,:,7) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);
r(:,:,8) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);


% 1.7) l6, v7, r8: q = 6   !!!!!!!!!
q    = 6; % change q
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];

while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        % finite difference along x + damping + forcing
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) - d .* ll(ii,1) ...
             + c2(q) .* fr3_8(n.*dt, xcords(ii)) ...
             + c3(q) .* fl3_6(n.*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4 .* dt;
        ll(ii,1) = ll(ii,2);
    end

    % Daily output to files
    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2 .* dx)) ...
                       + d2(q) .* fl3_6(n.*dt, xcords(ii)) ...
                       - d0(q) .* fr3_8(n.*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l63.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v73.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r83.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

% Reload as matrices
lll = load('C:\Users\pc\Desktop\data\l63.txt');
vvv = load('C:\Users\pc\Desktop\data\v73.txt');
rrr = load('C:\Users\pc\Desktop\data\r83.txt');
l(:,:,7) = reshape(lll,[ngrids_x length(lll)/ngrids_x]);
v(:,:,8) = reshape(vvv,[ngrids_x length(lll)/ngrids_x]);
r(:,:,9) = reshape(rrr,[ngrids_x length(lll)/ngrids_x]);


% 1.8) l7, v8, r9: q = 7   !!!!!!!!!
q    = 7; % change q
n    = 0;
ll   = zeros(ngrids_x,2);  % temporary array for l-field
vv   = zeros(ngrids_x,1);  % temporary array for v-field
lll  = [];                  % store l-field for output
vvv  = [];                  % store v-field for output
rrr  = [];                  % store r-field for output

while n <= ngrids_t*n_time
    for ii = 2:ngrids_x-1
        % finite difference along x + damping
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1);
        ll(ii,2) = ll(ii,1) + c4.*dt;
        ll(ii,1) = ll(ii,2);
    end

    % daily output to files
    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l73.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v83.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r93.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

% reload as matrices
lll = load('C:\Users\pc\Desktop\data\l73.txt');
vvv = load('C:\Users\pc\Desktop\data\v83.txt');
rrr = load('C:\Users\pc\Desktop\data\r93.txt');
l(:,:,8)  = reshape(lll, [ngrids_x length(lll)/ngrids_x]);
v(:,:,9)  = reshape(vvv, [ngrids_x length(lll)/ngrids_x]);
r(:,:,10) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]);


% 1.9) l8, v9, r10: q = 8   !!!!!!!!!
q    = 8; % change q
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];

while n <= ngrids_t*n_time
    for ii = 2:ngrids_x-1
        % finite difference + damping + 3rd baroclinic mode forcing
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) ...
             + c2(q) .* fr3_10(n.*dt, xcords(ii)) ...
             + c3(q) .* fl3_8(n.*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4.*dt;
        ll(ii,1) = ll(ii,2);
    end

    % daily output to files
    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) ...
                       + d2(q) .* fl3_8(n.*dt, xcords(ii)) ...
                       - d0(q) .* fr3_10(n.*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l83.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v93.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r103.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

% reload as matrices
lll = load('C:\Users\pc\Desktop\data\l83.txt');
vvv = load('C:\Users\pc\Desktop\data\v93.txt');
rrr = load('C:\Users\pc\Desktop\data\r103.txt');
l(:,:,9)  = reshape(lll, [ngrids_x length(lll)/ngrids_x]);
v(:,:,10) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]);
r(:,:,11) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]);


% 1.11) l10, v11, r12: q = 10   !!!!!!!!!
q    = 10; % change q
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];

while n <= ngrids_t*n_time
    for ii = 2:ngrids_x-1
        % finite difference + damping + 3rd baroclinic mode forcing
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) ...
             + c2(q) .* fr3_12(n.*dt, xcords(ii)) ...
             + c3(q) .* fl3_10(n.*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4.*dt;
        ll(ii,1) = ll(ii,2);
    end

    % daily output to files
    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) ...
                       + d2(q) .* fl3_10(n.*dt, xcords(ii)) ...
                       - d0(q) .* fr3_12(n.*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l103.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v113.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r123.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

% reload as matrices
lll = load('C:\Users\pc\Desktop\data\l103.txt');
vvv = load('C:\Users\pc\Desktop\data\v113.txt');
rrr = load('C:\Users\pc\Desktop\data\r123.txt');
l(:,:,11) = reshape(lll, [ngrids_x length(lll)/ngrids_x]);
v(:,:,12) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]);
r(:,:,13) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]);

% 1.13) l12, v13, r14: q = 12   !!!!!!!!!
q    = 12; % change q
n    = 0;
ll   = zeros(ngrids_x,2);  % temporary array for l-field
vv   = zeros(ngrids_x,1);  % temporary array for v-field
lll  = [];                  % store l-field for output
vvv  = [];                  % store v-field for output
rrr  = [];                  % store r-field for output

while n <= ngrids_t*n_time
    for ii = 2:ngrids_x-1
        % finite difference along x + damping + 3rd baroclinic mode forcing
        c4 = c1(q).* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) ...
             + c2(q).* fr3_14(n.*dt, xcords(ii)) ...
             + c3(q).* fl3_12(n.*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4.*dt;
        ll(ii,1) = ll(ii,2);
    end

    % daily output to files
    n_mod = mod(n,87*2);
    k_mod = mod(n,87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q).* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) ...
                       + d2(q).* fl3_12(n.*dt, xcords(ii)) ...
                       - d0(q).* fr3_14(n.*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l123.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v133.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r143.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)).* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

% reload as matrices
lll = load('C:\Users\pc\Desktop\data\l123.txt');
vvv = load('C:\Users\pc\Desktop\data\v133.txt');
rrr = load('C:\Users\pc\Desktop\data\r143.txt');
l(:,:,13) = reshape(lll, [ngrids_x length(lll)/ngrids_x]);
v(:,:,14) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]);
r(:,:,15) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]);


% 1.15) l14, v15, r16: q = 14   !!!!!!!!!
q    = 14; % change q
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];

while n <= ngrids_t*n_time
    for ii = 2:ngrids_x-1
        c4 = c1(q).* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) ...
             + c2(q).* fr3_16(n.*dt, xcords(ii)) ...
             + c3(q).* fl3_14(n.*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4.*dt;
        ll(ii,1) = ll(ii,2);
    end

    n_mod = mod(n,87*2);
    k_mod = mod(n,87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q).* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) ...
                       + d2(q).* fl3_14(n.*dt, xcords(ii)) ...
                       - d0(q).* fr3_16(n.*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l143.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v153.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r163.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)).* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

lll = load('C:\Users\pc\Desktop\data\l143.txt');
vvv = load('C:\Users\pc\Desktop\data\v153.txt');
rrr = load('C:\Users\pc\Desktop\data\r163.txt');
l(:,:,15) = reshape(lll, [ngrids_x length(lll)/ngrids_x]);
v(:,:,16) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]);
r(:,:,17) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]);


% 1.17) l16, v17, r18: q = 16   !!!!!!!!!
q    = 16; % change q
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];

while n <= ngrids_t*n_time
    for ii = 2:ngrids_x-1
        c4 = c1(q).* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) ...
             + c2(q).* fr3_18(n.*dt, xcords(ii)) ...
             + c3(q).* fl3_16(n.*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4.*dt;
        ll(ii,1) = ll(ii,2);
    end

    n_mod = mod(n,87*2);
    k_mod = mod(n,87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q).* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) ...
                       + d2(q).* fl3_16(n.*dt, xcords(ii)) ...
                       - d0(q).* fr3_18(n.*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l163.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v173.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r183.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

lll = load('C:\Users\pc\Desktop\data\l163.txt');
vvv = load('C:\Users\pc\Desktop\data\v173.txt');
rrr = load('C:\Users\pc\Desktop\data\r183.txt');
l(:,:,17) = reshape(lll, [ngrids_x length(lll)/ngrids_x]);
v(:,:,18) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]);
r(:,:,19) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]);


% 1.19) l18, v19, r20: q = 18   !!!!!!!!!
q    = 18; % change q
n    = 0;
ll   = zeros(ngrids_x,2);  % temporary array for l-field
vv   = zeros(ngrids_x,1);  % temporary array for v-field
lll  = [];                  % store l-field for output
vvv  = [];                  % store v-field for output
rrr  = [];                  % store r-field for output

while n <= ngrids_t*n_time
    for ii = 2:ngrids_x-1
        % finite difference along x + damping + 3rd baroclinic mode forcing
        c4 = c1(q).* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) ...
             + c2(q).* fr3_20(n.*dt, xcords(ii)) ...
             + c3(q).* fl3_18(n.*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4.*dt;
        ll(ii,1) = ll(ii,2);
    end

    % daily output to files
    n_mod = mod(n,87*2);
    k_mod = mod(n,87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q).* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) ...
                       + d2(q).* fl3_18(n.*dt, xcords(ii)) ...
                       - d0(q).* fr3_20(n.*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l183.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v193.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r203.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

% reload as matrices
lll = load('C:\Users\pc\Desktop\data\l183.txt');
vvv = load('C:\Users\pc\Desktop\data\v193.txt');
rrr = load('C:\Users\pc\Desktop\data\r203.txt');
l(:,:,19) = reshape(lll, [ngrids_x length(lll)/ngrids_x]);
v(:,:,20) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]);
r(:,:,21) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]);


% 1.21) l20, v21, r22: q = 20   !!!!!!!!!
q    = 20; % change q
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];

while n <= ngrids_t*n_time
    for ii = 2:ngrids_x-1
        c4 = c1(q).* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) ...
             + c2(q).* fr3_22(n.*dt, xcords(ii)) ...
             + c3(q).* fl3_20(n.*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4.*dt;
        ll(ii,1) = ll(ii,2);
    end

    n_mod = mod(n,87*2);
    k_mod = mod(n,87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q).* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) ...
                       + d2(q).* fl3_20(n.*dt, xcords(ii)) ...
                       - d0(q).* fr3_22(n.*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l203.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v213.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r223.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

lll = load('C:\Users\pc\Desktop\data\l203.txt');
vvv = load('C:\Users\pc\Desktop\data\v213.txt');
rrr = load('C:\Users\pc\Desktop\data\r223.txt');
l(:,:,21) = reshape(lll, [ngrids_x length(lll)/ngrids_x]);
v(:,:,22) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]);
r(:,:,23) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]);


% 1.23) l22, v23, r24: q = 22   !!!!!!!!!
q    = 22; % change q
n    = 0;
ll   = zeros(ngrids_x,2);
vv   = zeros(ngrids_x,1);
lll  = [];
vvv  = [];
rrr  = [];

while n <= ngrids_t*n_time
    for ii = 2:ngrids_x-1
        c4 = c1(q).* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) - d.*ll(ii,1) ...
             + c2(q).* fr3_24(n.*dt, xcords(ii)) ...
             + c3(q).* fl3_22(n.*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4.*dt;
        ll(ii,1) = ll(ii,2);
    end

    n_mod = mod(n,87*2);
    k_mod = mod(n,87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q).* ((ll(ii+1,1) - ll(ii-1,1)) ./ (2.*dx)) ...
                       + d2(q).* fl3_22(n.*dt, xcords(ii)) ...
                       - d0(q).* fr3_24(n.*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l223.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v233.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r243.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)./(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

lll = load('C:\Users\pc\Desktop\data\l223.txt');
vvv = load('C:\Users\pc\Desktop\data\v233.txt');
rrr = load('C:\Users\pc\Desktop\data\r243.txt');
l(:,:,23) = reshape(lll, [ngrids_x length(lll)/ngrids_x]);
v(:,:,24) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]);
r(:,:,25) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]);


%=======================================================================
% 3rd Baroclinic Mode: q = 24, 26, 28
%=======================================================================

%% 1.25) l24, v25, r26: q = 24
q = 24; 
n = 0;
ll = zeros(ngrids_x,2);
vv = zeros(ngrids_x,1);
lll = []; vvv = []; rrr = [];

while n <= ngrids_t*n_time
    for ii = 2:ngrids_x-1
        % finite difference + damping + 3rd baroclinic forcing
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) / (2*dx)) - d.*ll(ii,1) ...
             + c2(q).*fr3_26(n*dt, xcords(ii)) ...
             + c3(q).*fl3_24(n*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4*dt;
        ll(ii,1) = ll(ii,2);
    end

    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) / (2*dx)) ...
                       + d2(q).*fl3_24(n*dt, xcords(ii)) ...
                       - d0(q).*fr3_26(n*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l243.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v253.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r263.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)/(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

lll = load('C:\Users\pc\Desktop\data\l243.txt');
vvv = load('C:\Users\pc\Desktop\data\v253.txt');
rrr = load('C:\Users\pc\Desktop\data\r263.txt');
l(:,:,25) = reshape(lll, [ngrids_x length(lll)/ngrids_x]);
v(:,:,26) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]);
r(:,:,27) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]);


%% 1.27) l26, v27, r28: q = 26
q = 26; 
n = 0;
ll = zeros(ngrids_x,2);
vv = zeros(ngrids_x,1);
lll = []; vvv = []; rrr = [];

while n <= ngrids_t*n_time
    for ii = 2:ngrids_x-1
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) / (2*dx)) - d.*ll(ii,1) ...
             + c2(q).*fr3_28(n*dt, xcords(ii)) ...
             + c3(q).*fl3_26(n*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4*dt;
        ll(ii,1) = ll(ii,2);
    end

    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) / (2*dx)) ...
                       + d2(q).*fl3_26(n*dt, xcords(ii)) ...
                       - d0(q).*fr3_28(n*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l263.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v273.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r283.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)/(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

lll = load('C:\Users\pc\Desktop\data\l263.txt');
vvv = load('C:\Users\pc\Desktop\data\v273.txt');
rrr = load('C:\Users\pc\Desktop\data\r283.txt');
l(:,:,27) = reshape(lll, [ngrids_x length(lll)/ngrids_x]);
v(:,:,28) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]);
r(:,:,29) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]);


%% 1.29) l28, v29, r30: q = 28
q = 28; 
n = 0;
ll = zeros(ngrids_x,2);
vv = zeros(ngrids_x,1);
lll = []; vvv = []; rrr = [];

while n <= ngrids_t*n_time
    for ii = 2:ngrids_x-1
        c4 = c1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) / (2*dx)) - d.*ll(ii,1) ...
             + c2(q).*fr3_30(n*dt, xcords(ii)) ...
             + c3(q).*fl3_28(n*dt, xcords(ii));
        ll(ii,2) = ll(ii,1) + c4*dt;
        ll(ii,1) = ll(ii,2);
    end

    n_mod = mod(n, 87*2);
    k_mod = mod(n, 87);
    if k_mod == 0
        ll(ii,1) = smoothdata(ll(ii,1));
        ll(ii,2) = smoothdata(ll(ii,2));
    end
    if n_mod == 0
        for ii = 2:ngrids_x-1
            vv(ii,1) = d1(q) .* ((ll(ii+1,1) - ll(ii-1,1)) / (2*dx)) ...
                       + d2(q).*fl3_28(n*dt, xcords(ii)) ...
                       - d0(q).*fr3_30(n*dt, xcords(ii));
        end
        fid1 = fopen('C:\Users\pc\Desktop\data\l283.txt','a');
        fid2 = fopen('C:\Users\pc\Desktop\data\v293.txt','a');
        fid3 = fopen('C:\Users\pc\Desktop\data\r303.txt','a');
        fprintf(fid1,'%e\n', ll(:,1));
        fprintf(fid2,'%e\n', vv(:,1));
        fprintf(fid3,'%e\n', sqrt((q+1)/(q+2)) .* ll(:,1));
        fclose(fid1); fclose(fid2); fclose(fid3);
    end
    n = n + 1;
end

lll = load('C:\Users\pc\Desktop\data\l283.txt');
vvv = load('C:\Users\pc\Desktop\data\v293.txt');
rrr = load('C:\Users\pc\Desktop\data\r303.txt');
l(:,:,29) = reshape(lll, [ngrids_x length(lll)/ngrids_x]);
v(:,:,30) = reshape(vvv, [ngrids_x length(lll)/ngrids_x]);
r(:,:,31) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]);




% 2) v0  !!!!!!!!!
n  = 0;
nn = 1;
while n <= ngrids_t * n_time
    n_mod = mod(n, 87*2); 
    if (n_mod == 0)
        for ii = 2:ngrids_x-1
            % symmetric: v(ii,nn,1) = -sqrt(2) .* fr1_1(n.*dt, xcords(ii));
            v(ii, nn, 1) = 0;
        end
        nn = nn + 1;
    end
    n = n + 1;
end

% 3) r0  !!!!!!!!!
n    = 0;
h1   = 1 / (3 * e);
rr   = zeros(ngrids_x, 2);
rrr  = [];
while n <= ngrids_t * n_time
    for ii = 2:ngrids_x-1
        c4 = -h1 .* ((rr(ii+1,1) - rr(ii-1,1)) ./ (2*dx)) - d .* rr(ii,1) + ...
             fr3_0(n*dt, xcords(ii)); % modified
        rr(ii,2) = rr(ii,1) + c4 .* dt;
        rr(ii,1) = rr(ii,2);
    end
    % output to a file every day
    n_mod = mod(n, 87*2); 
    k_mod = mod(n, 87);
    if (k_mod == 0)
        ll(ii,1) = smoothdata(rr(ii,1));
        ll(ii,2) = smoothdata(rr(ii,2));
    end
    if (n_mod == 0)
        fid3 = fopen('C:\Users\pc\Desktop\data\r03.txt', 'a'); % modified
        fprintf(fid3, '%e\n', rr(:,1));
        fclose(fid3);
    end
    n = n + 1;
end

% reload as a matrix
rrr = load('C:\Users\pc\Desktop\data\r03.txt'); % modified
r(:,:,1) = reshape(rrr, [ngrids_x length(lll)/ngrids_x]); % modified

%========================================================
% Sum of Meridional Modes
%========================================================

s1 = [];
s2 = [];
s3 = [];

for ii = 1:size(l,2)
    s3 = 0;
    for jj = 1:N-2
        s1 = B3(jj-1, ycords);
        s2 = l(:, ii, jj);
        s3 = s3 + s1 .* s2;  
    end
    l_N(:,:,ii) = s3;
end

s1 = [];
s2 = [];
s3 = [];

for ii = 1:size(l,2)
    s3 = 0;
    for jj = 1:N-2
        s1 = B3(jj-1, ycords);
        s2 = r(:, ii, jj);
        s3 = s3 + s1 .* s2;
    end
    r_N(:,:,ii) = s3;
end

s1 = [];
s2 = [];
s3 = [];

for ii = 1:size(l,2)
    s3 = 0;
    for jj = 1:N-2
        s1 = B3(jj-1, ycords);
        s2 = v(:, ii, jj);
        s3 = s3 + s1 .* s2;
    end
    v_N(:,:,ii) = s3;
end

u_N3 = r_N - l_N;
p_N3 = r_N + l_N;
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

[~, I0] = min(abs(ycords(:) - 0 ./ y_scale));
[~, I1] = min(abs(zcords(:) - 0 ./ z_scale));
[~, I2] = min(abs(zcords(:) - 2 ./ z_scale));
[~, I3] = min(abs(zcords(:) - 4 ./ z_scale));
[~, I4] = min(abs(zcords(:) - 12 ./ z_scale));

%*******************************************************************

days      = 80;                             % Number of days to view
H         = 10;                             % Number of five-point smoothing iterations
S         = 1./2;                           % Five-point smoothing coefficient
xmin      = -6;                             % Minimum x-axis value
xmax      = 6;                              % Maximum x-axis value
ymin      = -3000 / y_scale;                % Minimum y-axis value
ymax      = 3000 / y_scale;                 % Maximum y-axis value
xtick1    = num2str(xmin .* x_scale ./ 1000); % Geographic label for x-axis min
xtick2    = num2str(xmax .* x_scale ./ 1000); % Geographic label for x-axis max
ytick1    = num2str(ymin .* y_scale ./ 1000); % Geographic label for y-axis min
ytick2    = num2str(ymax .* y_scale ./ 1000); % Geographic label for y-axis max
lowC      = -1;                             % Colorbar minimum value
highC     = 1;                              % Colorbar maximum value
LevelStep = 1;                              % Contour interval
LineWidth = 0.8;                            % Contour line width
stride1   = 5;                              % Stride for quiver plot x-direction
stride2   = 2;                              % Stride for quiver plot y-direction

%*******************************************************************

%                       PLOT: fig.5.
%========================================================
figure

[xx1, yy] = meshgrid(linspace(xgrids1, xgrids2, ngrids_x), ...
                     linspace(ygrids1, ygrids2, ngrids_y));

s      = [];
s(:,:) = P(:,:,I1,days);
s_sm   = s;

for kk = 1:H
    for ii = 2:size(s,1)-1
        for jj = 2:size(s,2)-1
            s_sm(ii,jj) = s(ii,jj) + 0.5 .* S .* (1 - S) .* ...
                (s(ii+1,jj) + s(ii-1,jj) + s(ii,jj+1) + s(ii,jj-1) - 4 .* s(ii,jj)) + ...
                0.25 .* S.^2 .* ...
                (s(ii+1,jj+1) + s(ii+1,jj-1) + s(ii-1,jj+1) + s(ii-1,jj-1) - 4 .* s(ii,jj));
        end
    end
    s = s_sm;
end

subplot(4,1,1)
aa1 = s_sm;
aa2 = s_sm;
aa1(aa1 < 0) = 0;
aa2(aa2 > 0) = 0;
contour(xx1, yy, aa1', '-k', 'LineWidth', LineWidth);
hold on
contour(xx1, yy, aa2', '--k', 'LineWidth', LineWidth);
title('(a) 0km')
set(gca, 'XTick', [xmin 0 xmax])
set(gca, 'XTickLabel', {xtick1, '0', xtick2}) 
set(gca, 'YTick', [ymin 0 ymax])
set(gca, 'YTickLabel', {ytick1, '0', ytick2})
set(gca, 'clim', [lowC highC])
axis([xmin xmax ymin ymax])
hold on

s1    = [];
u_p   = U(:,:,I1,days);
v_p   = V(:,:,I1,days);
u_pp  = u_p(1:stride1:end, 1:stride2:end);
v_pp  = v_p(1:stride1:end, 1:stride2:end);
uv    = sqrt(u_p.^2 + v_p.^2);
s1(1) = max(max(uv));

[x_hrz_wnd, y_hrz_wnd] = meshgrid(linspace(xgrids1, xgrids2, size(u_pp,2)), ...
                                  linspace(ygrids1, ygrids2, size(u_pp,1)));

quiver(x_hrz_wnd, y_hrz_wnd, u_pp, v_pp, '-k', 'MaxHeadSize', 3, 'AutoScale', 'on')



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

[xx2, yy2] = meshgrid(linspace(xgrids1, xgrids2, ngrids_x), ...
                       linspace(zgrids1, zgrids2, ngrids_z)); 

s         = [];
s(:,:)    = U(:, I0, :, days);
s_sm      = s;
H         = 100;       % Number of five-point smoothing iterations
kk        = 1;

for kk = 1:H
    for ii = 2:size(s,1)-1
        for jj = 2:size(s,2)-1
            s_sm(ii,jj) = s(ii,jj) + 0.5 .* S .* (1 - S) .* ...
                (s(ii+1,jj) + s(ii-1,jj) + s(ii,jj+1) + s(ii,jj-1) - 4 .* s(ii,jj)) + ...
                0.25 .* S.^2 .* ...
                (s(ii+1,jj+1) + s(ii+1,jj-1) + s(ii-1,jj+1) + s(ii-1,jj-1) - 4 .* s(ii,jj));
        end
    end
    s = s_sm;
end

subplot(2,1,1)
contourf(xx2, yy2, s_sm');
colormap gray
title('(a)')
xlabel('x (1000 km)')
ylabel('z (km)')
% set(gca,'clim',[lowC highC])
set(gca, 'XTick', [-6 -3 0 3 6])
set(gca, 'XTickLabel', {'-9','-4.5','0','4.5','9'})
set(gca, 'YTick', [0 max(zcords)])
set(gca, 'YTickLabel', {'0','16'})
set(gca, 'clim', [-10 10])

% Draw vertical lines perpendicular to x-axis
line([-3 -3], [0 max(zcords)], 'Color', 'k');
line([-2 -2], [0 max(zcords)], 'Color', 'k');
line([-1 -1], [0 max(zcords)], 'Color', 'k');
line([0 0], [0 max(zcords)], 'Color', 'k');
line([1 1], [0 max(zcords)], 'Color', 'k');
line([2 2], [0 max(zcords)], 'Color', 'k');
line([3 3], [0 max(zcords)], 'Color', 'k');

subplot(2,1,2)
[~, II1] = min(abs(xcords(:) - (-3)));
[~, II2] = min(abs(xcords(:) - (-2)));
[~, II3] = min(abs(xcords(:) - (-1)));
[~, II4] = min(abs(xcords(:) - 0));
[~, II5] = min(abs(xcords(:) - 1));
[~, II6] = min(abs(xcords(:) - 2));
[~, II7] = min(abs(xcords(:) - 3));

uu(1,:) = s_sm(II1,:);
uu(2,:) = s_sm(II2,:) + 10;
uu(3,:) = s_sm(II3,:) + 20;
uu(4,:) = s_sm(II4,:) + 30;
uu(5,:) = s_sm(II5,:) + 40;
uu(6,:) = s_sm(II6,:) + 50;
uu(7,:) = s_sm(II7,:) + 60;

plot(uu','k')
title('(b)')
xlabel('z (km)')
ylabel('U (m/s)')
view(90, -90)

% Add horizontal reference lines
line([0 200], [0 0], 'Color', 'k');
line([0 200], [10 10], 'Color', 'k');
line([0 200], [20 20], 'Color', 'k');
line([0 200], [30 30], 'Color', 'k');
line([0 200], [40 40], 'Color', 'k');
line([0 200], [50 50], 'Color', 'k');
line([0 200], [60 60], 'Color', 'k');

set(gca, 'YTick', [])
set(gca, 'XTick', [0 200])
set(gca, 'XTickLabel', {'0','16'})






















