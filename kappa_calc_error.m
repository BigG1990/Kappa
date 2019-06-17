function [kappa,kappa_geomean,kappa_geostd_down,kappa_geostd_up,time_CCN,Dc_mean_old,Dc_error_old,CCN,T_inlet,SS_kali,kappa_p] = kappa_calc_error(time_CCN,CCN,time_SD,SD,d,SS_kali,T_inlet)

% This function is used to calculate the particle hygroscopicity properties, a sigle parameter kappa.
% The Monte Carlo simulation is included to get measurement uncertainties.

% Call as [kappa,kappa_geomean,kappa_geostd_down,kappa_geostd_up,time_CCN,Dc_mean_old,Dc_error_old,CCN,T_inlet,SS_kali,kappa_p] = kappa_calc_error(time_CCN,CCN,time_SD,SD,d,SS_kali,T_inlet)

% where kappa               include Monte Carlo simulation results
%       kappa_geomean       is the geomean value of the kappa
%       kappa_geostd_down	is one lower geometric standard deviation of kappa  
%       kappa_geostd_up     is one upper geometric standard deviation of kappa
%       time_CCN            is the time of CCN data
%       Dc_mean_old         is the median value of Dc
%       Dc_error_old        is one standard deviation of Dc
%       CCN                 is CCN number concentration
%       T_inlet             is inlet temperature
%       SS_kali             is calibrated SS
%       kappa_p             is 1st 2.5th 5th 10th 50-34.1th 25th 50th 75th 50+34.1th 90th 95th 97.5th 99th percentile of kappa
%       time_SD             is time of particle number size distribution (PNSD)
%       SD                  is PNSD
%       d                   is particle diameter


%
% Copyright (C) Feb 2017 - present Xianda Gong, Leibniz Institute for Tropospheric Research (TROPOS)



SD=fliplr(SD); %flip SD and diameter
d=fliplr(d);                   

Dc1=NaN*ones(length(CCN),1);
Dc2=Dc1;
for j=1:length(CCN)
    [rm cm]=min(abs(time_SD-time_CCN(j))); % find closest time periods
    if rm<0.5/24                                                              % calculate d_crit only within a certain time intervall
        for i=1:1000                                                        % Monte Carlo 1000 times
            gauss=randn(size(SD,2),1);                                      % randn gives a normal distribution with std=1 and mean=0
            SD_error_1=SD(cm,:)'*0.05;                                      % 5% error in the measured number concentration
            SD_error_2=gauss.*SD_error_1;                                   % multiplying the gaussian with the error of N
            SD_error=cumsum((SD(cm,:)+SD_error_2'));                        % adding the error of N to the signal of N
            [r1 s1]=find(SD_error>=CCN(j));                                 % find closest cummulative CN concentration to CCN concentration for the N error
            if size(s1,2)>0 & s1(1)>1                                    % do the calculations only if CCN concentrations exist, that are larger than the CN concentration. make sure subscribe induices must to larger than 0
                Dc1(j,i)=d(s1(1));                                          % closest discrete critical diameter for the N error
                if Dc1(j,i)>15 & Dc1(j,i)<600                               % here you can set the size range for which the interpolation is performed
                    Dc2(j,i)=interp1([SD_error(s1(1)),SD_error(s1(1)-1)],[d(s1(1)),d(s1(1)-1)],CCN(j));                        % lineare Interpolation: Formel mit Hand nachgerechnet --> passt    
                else
                    Dc2(j,:)=NaN;                                           % if conditions above are not fulfill then write NaN 
                end             
            else
                Dc2(j,:)=NaN;
            end
        end
    else
        Dc2(j,:)=NaN;
    end
    fprintf('%s %f %s\n','Completed', j/(length(CCN))*100,'%');
end

Dc_mean=nanmean(Dc2,2);
Dc_std=nanstd(Dc2,0,2);
for p=1:length(Dc_mean)
    Dc_error(p,:)=sqrt((0.03)^2+(Dc_std(p)/Dc_mean(p))^2)*Dc_mean(p);         % error propagation of the Dc_error due to N and the 3% (Ali-->Heike-->me) in d itselfe
end

% dele NaN rows:
[row, col] = find(isnan(Dc_mean));
Dc_mean(row,:)=[];
Dc_std(row,:)=[];
Dc_error(row,:)=[];
SS_kali(row,:)=[];
T_inlet(row,:)=[];
time_CCN(row,:)=[];
time_SD(row,:)=[];
CCN(row,:)=[];

N = 10000;              % Number of data points
R=8.31451;              %J/mol*K
Mw=0.018;               %kg/mol
rho_w=997.242;          %kg/m?
sigma=0.07274;          %N/m
Dc_mean_old=Dc_mean;
Dc_error_old=Dc_error;
Dc_mean=Dc_mean*1e-9;  %critical diameter converted to m
Dc_error=Dc_error*1e-9; %error on the critical diameter converted to m
kappa=NaN*ones(length(CCN),N);
for k=1:length(CCN)
    gauss2=randn(1,N);
    A(k,:)=4*sigma*Mw/R/rho_w/T_inlet(k);
    Dc3=Dc_mean(k)+Dc_error(k)*gauss2;
    %Dc3=Dc_mean(k)*ones(1,N);                      % d_crit ohne Monte Carlo um den Einfluss von SS zu testen
    if SS_kali<0.2
        SS_kali_error=SS_kali(k)+0.007*gauss2;      %die 0.007 sind die 1sigma (0.035*0.2) weil unterhalb von SS=0.2% ein absoluter Fehler auftritt
        %SS_kali_error=SS_kali(k)*ones(1,N);        % konst. SS (also ohne Monte Carlo) um den Einfluss von d_crit zu testen 
        kappa(k,:)=4*(A(k)^3)./(27*(Dc3.^3).*(log(SS_kali_error*10^-2+1)).^2);
    else
        SS_kali_error=SS_kali(k)+0.035*SS_kali(k)*gauss2;   %die 3,5% sind 1sigma, Heike hat die Unsicherheit f?r 2sigma als 7% bestimmt.
        kappa(k,:)=4*(A(k)^3)./(27*(Dc3.^3).*(log(SS_kali_error*10^-2+1)).^2);
    end
    [m n]=find(kappa(k,:)>1.3);
    kappa(k,n)=NaN;
    SS_test(k,:)=SS_kali_error;
end
id=find(kappa<0);
kappa(id)=NaN;
kappa_p=prctile(kappa,[1 2.5 5 10 50-34.1 25 50 75 50+34.1 90 95 97.5 99],2);  %percentiles for investigation of kappa distribution
kappa_std=nanstd(kappa,0,2);
kappa_mean=nanmean(kappa,2);
kappa_median=nanmedian(kappa,2);
kappa_geomean=geomean(kappa,2);
kappa_geomean=exp(nanmean(log(kappa),2));
kappa_geostd=geostd(kappa,0,2);     %geometrische Std

% The geometric / n std / mean is taken in the case of a log normal distribution to\
% take into account the distortion to large values (similar to the median).
% The geostd is calculated by logarithmic, then std / mean / calculates and exposes.
% For the unsymmetrical error bars, however, first calculate mean +/- std and then exp.
kappa_geostd_up=exp(log(kappa_geomean)+log(kappa_geostd));          %Mittlerer geomean der Monte Carlo Simulationen ( oben und unten)
kappa_geostd_down=exp(log(kappa_geomean)-log(kappa_geostd));

% kappa_geomean(isnan(kappa_geomean(:,1)),:)=[];
kappa_geostd_end_up=exp(log(geomean(kappa_geomean))+log(geostd(kappa_geomean)));
kappa_geostd_end_down=exp(log(geomean(kappa_geomean))-log(geostd(kappa_geomean)));