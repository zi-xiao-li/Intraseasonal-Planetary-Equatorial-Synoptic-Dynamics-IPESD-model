% ----------------------------------------------------------------
% Set the 3D synoptic-scale heat source over the equatorial region
% Specify upward momentum and heat transport induced by the heat source
% ----------------------------------------------------------------


% ============================================== 
%                synoptic Q                  
% ==============================================
clear all; clc

phi0    = pi/4;     % stratiform lag pi/4
a       = 1/3;      % ratio of stratiform to convective component
phit    = 0;        % arbitrary wave speed
e       = 0.125;
y0      = 0;
x_scale = 1000;     % km
X_scale = 12500;
y_scale = 1000;
t_scale = 3.3;      % days
z_scale = 16/pi;
u_scale = 6.25;     % m/s

% ===== Warm Pool =====

Lx      = 5000;     % Width of warm pool (dimensionless, synoptic scale)
Ly      = 2000; 
xmin    = -Lx/x_scale;
xmax    = Lx/x_scale;
ymin    = -Ly/y_scale;
ymax    = Ly/y_scale;
zmin    = 0;
zmax    = pi;

% [x, y, z] = meshgrid(linspace(xmin, xmax, 101), linspace(ymin, ymax, 81), ...
% linspace(zmin, zmax, 30));  % Note: in meshgrid, input1 corresponds to columns (longitude), input2 to rows (latitude)

syms x y z 

H      = exp(-2.*y.^2);
Hy     = diff(H);
F      = cos(pi.*x./(2.*Lx./x_scale));
G1     = H.*sin(x);
G2 	   = -a.*H.*sin(x+phi0);
G1y    = diff(G1,y);
G2y    = diff(G2,y);
G1x    = diff(G1,x);
G2x    = diff(G2,x);
Q_snp  = H.*(cos(x).*sin(z)-a.*cos(x+phi0).*sin(2.*z));
S_snp  = Q_snp.*F;
u_snp  = -((2.*G1+y.*G1y).*cos(z)+2.*(2.*G2+y.*G2y).*cos(2.*z));
v_snp  = y.*(G1x.*cos(z)+G2x.*cos(2.*z));
w_snp  = G1x.*sin(z)+G2x.*sin(2.*z);
p_snp  = y.^2.*(G1.*cos(z)+2.*G2.*cos(2.*z));
t_snp  = -y.^2.*(G1.*sin(z)+4.*G2.*sin(2*z));


% ============================================== 
%                      EQ (X-Z)                
% ==============================================

%--matlabFunction
Q_snp = matlabFunction(Q_snp);
S_snp = matlabFunction(S_snp);
u_snp = matlabFunction(u_snp);
v_snp = matlabFunction(v_snp);
w_snp = matlabFunction(w_snp);
p_snp = matlabFunction(p_snp);
t_snp = matlabFunction(t_snp);
%--eq
Q_snp_eq_func = matlabFunction(Q_snp(x,0,z));
S_snp_eq_func = matlabFunction(S_snp(x,0,z));
u_snp_eq_func = matlabFunction(u_snp(x,0,z));
v_snp_eq_func = matlabFunction(v_snp(x,0,z));
w_snp_eq_func = matlabFunction(w_snp(x,0,z));
p_snp_eq_func = matlabFunction(p_snp(x,0,z));
t_snp_eq_func = matlabFunction(t_snp(x,0,z));
%--plot
Q_snp_eq = [];
S_snp_eq = [];
u_snp_eq = [];
v_snp_eq = [];
w_snp_eq = [];
p_snp_eq = [];
t_snp_eq = [];
xx = linspace(xmin,xmax,101);
zz = linspace(zmin,zmax,30);
for ii=1:length(xx)
    for jj=1:length(zz)
        Q_snp_eq(ii,jj) = Q_snp_eq_func(xx(ii),zz(jj));
        S_snp_eq(ii,jj) = S_snp_eq_func(xx(ii),zz(jj));
        u_snp_eq(ii,jj) = u_snp_eq_func(xx(ii),zz(jj));
        w_snp_eq(ii,jj) = w_snp_eq_func(xx(ii),zz(jj));
        %v_snp_eq(ii,jj) = v_snp_eq_func(xx(ii),zz(jj)); = 0
        %p_snp_eq(ii,jj) = p_snp_eq_func(xx(ii),zz(jj)); = 0
        %t_snp_eq(ii,jj) = t_snp_eq_func(xx(ii),zz(jj));
     end
end
Q_snp_eq = Q_snp_eq';
S_snp_eq = S_snp_eq';
u_snp_eq = u_snp_eq';
w_snp_eq = w_snp_eq';
%v_snp_eq = v_snp_eq';
%p_snp_eq = p_snp_eq';
%t_snp_eq = t_snp_eq';

% -- Realistic plot

figure

% -- Color table
colortable = textread('WhiteBlueGreenYellowRed.txt');  % Read colormap from file
colormap(colortable);

% -- Coordinates for contour plot
[x_vtc, z_vtc] = meshgrid(linspace(xmin, xmax, 101), linspace(zmin, zmax, 30));

% -- Heat source (contour fill of synoptic variable)
contourf(x_vtc, z_vtc, S_snp_eq, 'LineStyle', 'none')
hold on

% -- Wind field (subsampled for clarity)
stride1 = 3;  % subsample in vertical
stride2 = 3;  % subsample in horizontal
u_snpp_eq = u_snp_eq(3:stride1:end-3, 3:stride2:end-3);
w_snpp_eq = w_snp_eq(3:stride1:end-3, 3:stride2:end-3);

% -- Coordinates for quiver plot
[x_vtc_wnd, z_vtc_wnd] = meshgrid(linspace(xmin, xmax, size(u_snpp_eq, 2)), ...
                                  linspace(zmin, zmax, size(u_snpp_eq, 1)));

% -- Overlay wind vectors
quiver(x_vtc_wnd, z_vtc_wnd, u_snpp_eq, w_snpp_eq, '-k', 'MaxHeadSize', 0.05)

% -- Axis labels
xlabel('Zonal (x 1000 km)')
ylabel('Vertical (z km)')
yticklabels({'0','2.5','5','7.5','10','12.5','15'})

% -- Colorbar and axes formatting
ax = gca;
c = colorbar;
ax.FontSize = 12;
c.LineWidth = 1.0;
c.FontName = 'Arial';
c.FontSize = 12;






% ============================================== 
%                FU FT                  
% ==============================================

k  = 3/4.*a.*sin(phi0);
FU = k.*(cos(z)-cos(3.*z)).*(2.*H.^2+y.*H.*Hy);
FU = matlabFunction(FU);
FT = -k.*(sin(z).*(5.*y.^2.*H.^2+4.*y.^3.*Hy.*H)+...
            (sin(3.*z)./3).*(15.*y.^2.*H.^2+4.*y.^3.*Hy.*H));
FT = matlabFunction(FT);



figure


colorbar_start_1 = -0.7;
colorbar_end_1   = 0.7;
colorbar_int_1   = 0.1;
colorbar_start_2 = -0.1;
colorbar_end_2   = 0.1;
colorbar_int_2   = 0.025;
colortable = textread('WhiteBlueGreenYellowRed.txt');
colormap(colortable);
% -- Coordinates
[y_vtc, z_vtc] = meshgrid(linspace(ymin, ymax, 81), linspace(zmin, zmax, 30));
% -- Variables / Data arrays
yy     = [];
zz     = [];
FU_vtc = [];
FT_vtc = [];
yy = linspace(ymin,ymax,81);
zz = linspace(zmin,zmax,30);
for ii=1:length(yy)
    for jj=1:length(zz)
        FU_vtc(ii,jj) = FU(yy(ii),zz(jj));
        FT_vtc(ii,jj) = FT(yy(ii),zz(jj));
    end
end

subplot(2,2,1)
contourf(y_vtc,z_vtc,FU_vtc','LineStyle','none')
xlabel('Meridional(x 1000km)')
ylabel('z(km)')
yticklabels({'0','5','10','15'})
title('a.                                                                 ','FontSize',16)
ax             = gca;
c              = colorbar;
ax.FontSize    = 12;
c.LineWidth    = 1.;
c.FontName     = 'Arial';
c.FontSize     = 12;
caxis([colorbar_start_1,colorbar_end_1])
set(c,'YTick',colorbar_start_1:colorbar_int_1:colorbar_end_1);


subplot(2,2,3)
contourf(y_vtc,z_vtc,FT_vtc','LineStyle','none')
xlabel('Meridional(x 1000km)')
ylabel('z(km)')
yticklabels({'0','5','10','15'})
title('b.                                                                 ','FontSize',16)
ax             = gca;
c              = colorbar;
ax.FontSize    = 12;
c.LineWidth    = 1.;
c.FontName     = 'Arial';
c.FontSize     = 12;
caxis([colorbar_start_2,colorbar_end_2])
set(c,'YTick',colorbar_start_2:colorbar_int_2:colorbar_end_2);

