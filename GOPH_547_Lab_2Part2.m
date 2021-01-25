clear all; close all;
G = 6.674e-11;
load('goph_547_lab_2_data_w2020.mat');
base_point=grav_survey_data(:,1);
day=grav_survey_data(:,2);
time=grav_survey_data(:,3);
total_time=grav_survey_data(:,4);
x=grav_survey_data(:,5);
y=grav_survey_data(:,6);
z=grav_survey_data(:,7);
g=grav_survey_data(:,8);
dg_tide=grav_survey_data(:,9);

%Data Preparation
%Q1
figure
g(1)=0;
p = polyfit(total_time,g,1);
f = polyval(p,total_time);
plot(total_time,g,'v',total_time,f,'-',...
                                        'MarkerSize',0.5,...
                                        'LineWidth',1.5);
hold on
plot(total_time,dg_tide,'LineWidth',1.5);
xlabel('total time[h]');
ylabel('gravity effect');
legend('gravity mearsurements','tidal viration','best fit line');
title('time vs gravity effect');
%Q2
figure
contourf(Xt,Yt,Zt);
colorbar;
xlabel('Eastings[km]');ylabel('Northings[km]');
title('elevation data[m]');
%display elevation information
disp(string(strcat('Mean elevation:',{' '},mat2str(min(min(Zt))))));
disp(string(strcat('Range of elevations:',{' '},mat2str(min(min(Zt))),{' '},'to',{' '},mat2str(max(max(Zt))))));
z_datum=min(min(Zt));% set the datum

%Q3
x_sort=zeros(size(x,1),3);
x_sort =[x y z];
[xsorted,ind_sort]=sortrows(x_sort,[1 2]);
[uniquex,ind_uniq,~]=unique(xsorted,'rows');
Xg=uniquex(:,1);
Yg=uniquex(:,2);
Zg=uniquex(:,3);
Nx=sqrt(length(Xg));
Ny=sqrt(length(Yg));

% creat the grids with reshape function
Xg1=reshape(Xg,[Nx,Ny]);
Yg1=reshape(Yg,[Nx,Ny]);
Zg1=reshape(Zg,[Nx,Ny]);
%creat the grids with meshgrid function 
Xgrid=0:50;
Ygrid=0:50;
[Xg2,Yg2]=meshgrid(Xgrid,Ygrid);
%verifing step
xdiff=Xg1-Xg2;
ydiff=Yg1-Yg2;
verifyx=find(xdiff);
verifyy=find(ydiff);

g(1)=grav_survey_data(1,8);
g_corr=g;
g_raw=g_corr;
g_raw(2:end)=g_raw(2:end)+g_raw(1);
g_raw=g_raw(ind_sort);
g_raw=g_raw(ind_uniq);
g_raw=reshape(g_raw,[Nx,Ny]);

figure
contourf(Xg1,Yg1,g_raw);
title('Measured Raw Gravity Effect');
xlabel('Easting/km');
ylabel('Northing/km');
colorbar;
%normal gravity correction
ge=9.7803253359;
k=0.00193185265241;
e2=0.0069437999013;
gt=ge*(1+k*(sin(49.1286*pi/180))^2)/sqrt(1-e2*(sin(49.1286*pi/180))^2);%using IGF compute gt
g_corr(1)=g_corr(1)-gt*10^8;
g_norm=g_corr;
g_norm(2:end)=g_norm(2:end)+g_norm(1);
g_norm=g_norm(ind_sort);
g_norm=g_norm(ind_uniq);
g_norm=reshape(g_norm,[Nx,Ny]);
figure
contourf(Xg1,Yg1,g_norm);
title('Normal Gravity Correction');
xlabel('Easting/km');
ylabel('Northing/km');
colorbar;
%%drift corrections
%tidal drift correction
g_corr=g_corr-dg_tide;
g_tide=g_corr;
g_tide(2:end)=g_tide(2:end)+g_tide(1);
g_tide=g_tide(ind_sort);
g_tide=g_tide(ind_uniq);
g_tide=reshape(g_tide,[Nx,Ny]);
figure
contourf(Xg1,Yg1,g_tide);
title('Normal Gravity and Tide Drift Correction');
xlabel('Easting/km');
ylabel('Northing/km');
colorbar;
%instrument drift correction
st=1;
for i = 2:length(base_point)  
    if base_point(i) == 1
        st = i;
        drift = g_corr(i); 
        g_corr(st:end) = g_corr(st:end) - drift; %subtract overnight drift
    end 
    
    if base_point(i) == 2
        drift = g_corr(i);       
        drift_rate = drift /(total_time(i) - total_time(st));
        g_corr(st+1:i-1) = g_corr(st+1:i-1) - drift_rate * ...
            (total_time(st+1:i-1) - total_time(st)); % subtract the time drift 
        g_corr(i:end) = g_corr(i:end) - drift;  % subtract the drift from points that are not the start nor end of loop
    end
end

% verify g_corr=0
for i=2:length(g_corr)
if (base_point(i)==1||base_point(i)==2)
        if(g_corr(i)~=0)
            disp('error');
        end
    end
end
g_corr(2:end)=g_corr(2:end)+g_corr(1);
%plot normal gravity and drift corrected data
figure
plot(total_time,g_corr);
xlabel('total time/h');
ylabel('gravity');
title('Normal gravity, tide and instrument drift corrected gravity')

g_corr=g_corr(ind_sort);
g_corr=g_corr(ind_uniq);
g_corr=reshape(g_corr,[Nx,Ny]);
g_drft=g_corr;
figure
contourf(Xg1,Yg1,g_drft);
title('normal gravity, tide and instrument drift corrected gravity');
xlabel('Easting/km');
ylabel('Northing/km');
colorbar;

%%elevation and terrain corrections
% free air correction
dg_FA=zeros(size(Zg));  
dz=Zg-z_datum;  
dg_FA=(-0.3086*dz)*10^3;  % convert mGal to uGal
dg_FA=reshape(dg_FA,[Nx,Ny]);
g_corr=g_corr-dg_FA;
g_FA=g_corr;


figure
contourf(Xg1,Yg1,g_FA);
title('free air correction');
xlabel('Easting/km');ylabel('Northing/km');
colorbar;
%Bourger plate correction
rho=2.65;
dg_BP=(0.04193*rho*dz)*10^3; %convert density from g/cc to kg/m^3
dg_BP((Zg-z_datum)>0)=0; 
dg_BP=reshape(dg_BP,[Nx,Ny]);
g_corr=g_corr-dg_BP;
g_elev=g_corr;


figure
contourf(Xg1,Yg1,g_elev);
title('Bourger plate correction');
xlabel('Easting/km');
ylabel('Northing/km');
colorbar;

%terrain correction
dg_terr=zeros(size(g_corr));
dx=1000;
dy=1000;
for i=1:numel(dg_terr)
    xi=[Xg1(i)*1000, Yg1(i)*1000,z_datum];
    for j=1:numel(Xt)
           xm=[Xt(j)*1000, Yt(j)*1000, 0.5*(Zt(j)+z_datum)];
           dm=rho*10^3*(Zt(j)-z_datum)*dx*dy;
           dg_terr(i)=dg_terr(i)+abs(grav_eff_point(xi,xm,dm,G));
    end
end

dg_terr=dg_terr*10^8; %convert unit from m/s^2 to uGal
g_corr=g_corr+dg_terr;
g_terr=g_corr;
%plot g_terr

figure
contourf(Xg1,Yg1,g_terr);
title('terrain correction');
xlabel('Easting/km');
ylabel('Northing/km');
colorbar;


dg_rgn1=mean(mean(g_corr));
g_corr=g_corr-dg_rgn1;
g_anom=g_corr;


figure
contourf(Xg1,Yg1,g_anom);
title('Bouguer anomaly after regional correction');
xlabel('Easting/km');
ylabel('Northing/km');
colorbar
%range of gravity anomaly
disp('Range of the gravity anomaly in uGal:');
fprintf('%.4e : %.4e\n\n',min(min(g_anom)),max(max(g_anom)));

%calculate the derivatives of gravity anomaly
dgdx=zeros(size(g_anom));
dgdy=zeros(size(g_anom));
d2gdz2=zeros(size(g_anom));
dgdx(:,2:end)=(g_anom(:,2:end)-g_anom(:,1:end-1))/dx;
dgdy(2:end,:)=(g_anom(2:end,:)-g_anom(1:end-1,:))/dy;
d2gdz2(2:end-1,2:end-1) = (2*(dx^2+dy^2)*g_anom(2:end-1,2:end-1)-( dx^2*(g_anom(3:end,2:end-1)+g_anom(1:end-2,2:end-1))+dy^2*(g_anom(2:end-1,3:end)+g_anom(2:end-1,1:end-2))))/(dx^2*dy^2);
%plot the first derivatives

figure
subplot(1,2,1);
contourf(Xg1,Yg1,dgdx);colorbar;axis equal;
title('dg/dx');xlabel('Easting/km');ylabel('Northing/km');
subplot(1,2,2);
contourf(Xg1,Yg1,dgdy);colorbar;axis equal;
title('dg/dy');xlabel('Easting/km');ylabel('Northing/km');

%plot the second derivative

figure
contourf(Xg1,Yg1,d2gdz2);colorbar;
title('second derivative of gravity effect');xlabel('Easting/km');ylabel('Northing/km');

%find the max g_anom 
[xi_anom,yi_anom]=find(g_anom==max(max(g_anom)));
maxanom=g_anom(xi_anom,yi_anom);
x_anom=Xg1(xi_anom,yi_anom);
y_anom=Yg1(xi_anom,yi_anom);
masslocation=[x_anom,y_anom,z_datum-500];
slocation=[x_anom,y_anom,z_datum];
mass=maxanom*10^(-8)*norm(masslocation-slocation)^3/(G*500);
disp('The total mass of anomalies is around:');
fprintf('%.4e kg\n',mass);


