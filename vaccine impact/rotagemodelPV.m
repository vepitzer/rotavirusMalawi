clear global all
global B wm wi1 wi2 u um beta b1 phi1 d1 d2 rr1 rr2 ri2 ri3 al v1 v2 v3 sc1 sc2 sc3 sc2n sc3n wv reintro wA; % immig deaths b2 phi2 

%%% FIXED PARAMETERS %%%
age=[0:1/12:23/12 2:4 5:5:75];
%age=[0:1/12:11/12 1:4 5:5:75];
al=length(age);
au2=24; %13; %
au5=27; %16; %
avgage=[1:24 30 42 54];
%avgage=[1:12 18 30 42 54];
agep=[(1/12)*malpop08(1)*ones(1,12) (1/48)*malpop08(2)*ones(1,12) (1/4)*malpop08(2)*ones(1,3) malpop08(3:17).*ones(1,15)];
%agep=[(1/12)*malpop08(1)*ones(1,12) (1/4)*malpop08(2)*ones(1,4) malpop08(3:17).*ones(1,15)];
agep=agep/sum(agep);
u=[1/4.3*ones(1,24) 1/52*ones(1,3) 1/(52*5)*ones(1,14) 1/(52*25)];
%u=[1/4.3*ones(1,12) 1/52*ones(1,4) 1/(52*5)*ones(1,14) 1/(52*25)];

tmax=round(52.18*55); 
t0=round(52.18*30); 
tvacc=757;
tvf=t0+tvacc+547;
datesim=(datenum([1967 7 2 0 0 0]):7:datenum([2022 7 1 0 0 0]))';

malawi_cbr=[47.6; 48.3; 49.9; 51.5; 53.4; 55; 53.9; 52.2; 48.2; 47.1; 43.4; 41.8; 39.9; 38.3; 36.6; 35; 33.2; 31.3]; % crude birth rate from 1965-70 to 2020-25
malawi_cdr=[27.8; 26.5; 25.5; 24.9; 23.2; 22.5; 20.9; 18.7; 18.3; 18.2; 16.6; 13.5; 11.5; 10.5; 9.6; 8.8; 8.2; 7.7]; % crude death rate from 1965-70 to 2020-25
datepop=[1950 1957:5:2037; 1 7*ones(1,17); ones(1,18)]';
N1=870000;  % Blantyre population in 2000
N=(N1/exp(.04*30))*agep; % Adjust initial population for rate of exponential growth and multiple by age distribution
B=interp1(datenum(datepop),malawi_cbr/1000,datesim); % interpolate weekly birth rate
if t0>0 && length(B)<tmax
    B=[B; B(end)*ones(tmax-length(B),1)]; % pad end with last available value
end
B=[B zeros(tmax,al-1)]; % birth rate by age (only into first age group)
cdr=interp1(datenum(datepop),malawi_cdr/1000,datesim); % interpolate weekly death rate
cdr=log(1+cdr)/52.18; % convert from annual rate to weekly rate
immig=.012*ones(length(cdr),1); % annual immigration rate (estimated)
immig=log(1+immig)/52.18; % weekly immigration rate
um=.010*ones(1,al); % annual mortality rate from all age groups (estimated)
um=log(1+um)/52.18; % weekly mortality rate

dur=1; %duration of infectiousness
d1=1/dur; %rate of recovery from primary infection (per week)
d2=2*d1; %rate of recovery from subsequent infection (per week)
rr1=0.62;  %relative risk of second infection 
rr2=0.35; %relative risk of third infection 
ri2=0.5;   %relative infectiousness of secondary infection
%ri3=0.1; %relative infectiousness of asymptomatic infection
wm=1/4.3; %rate of waning maternal immunity (avg duration of immunity = 1mo) 
wi1=1/13; %rate of waning immunity following primary infection
wi2=1/13; %rate of waning immunity following 2nd infection
wA=0; %rate of waning immunity against symptomatic infection in adults


%%% ESTIMATED PARAMETERS %%%
pfit_mal=[78.7719 47.3943 34.0096 26.5937 21.8076;...
    0.1740 0.1450 0.1345 0.1301 0.1275;...
    6.8559 7.7299 8.0841 8.2030 8.3139;...
    0.0171 0.0171 0.0171 0.0170 0.0171];   

for relinf=1%:5

pars=pfit_mal(:,relinf); 
ri3=relinf/10; 

ptrans=pars(1); %~R0
b1=pars(2); %seasonal forcing 
phi1=pars(3); %seasonal offset
ar1=1; %relative risk <1y
ar2=1; %relative risk 1-<2y

%Variation in reporting effort over time (by age)
repeff=[repeff_prevacc;...
    ones(26,al);... 
    repeff_prevacc08;...
    ones(104,al);...
    repeffPV ones(547,al-au5);...
    ones(tmax-tvf,1)*repeffPV(end,:) ones(tmax-tvf,al-au5)];
    
%c1=100*ones(al); 
c2=100*[ar1*ones(al,12) ar2*ones(al,au2-12) ones(al,al-au2)];

beta=(ptrans/100)*c2; 

h=pars(4); %proportion of severe diarrhea cases hospitalized at QECH
hosp1=.13*h*repeff; %proportion of first infections resulting in severe disease and reported
hosp2=.063*h*repeff; %proportion of second infections resulting in severe disease and reported
hosp3=zeros(tmax-t0,al); %proportion of subsequent infections resulting in severe disease and reported
delta1=0.41*h*ones(1,al); %proportion of first infections that are symptomatic
delta2=0.35*h*ones(1,al); %proportion of second infections that are symptomatic
delta3=0.21*h*ones(1,al); %rate of detection of subsequent infection

reintro=0;

R0=max(eig(dur*beta.*(ones(al,1)*N)));


%%% VACCINATION PARAMETERS

for targ=3 %1=all-or-nothing protection following first dose; 2=all-or-nothing protection following second dose; 3=incremental protection following each dose; 4=booster dose
for v=1:2 

if v==1
    vcov=zeros(length(vcovPV_mavg1),2);
else
    vcov=vcovPV_mavg1;
end
avacc=[1 3 3 3; %1=birth, 2=1mo, 3=2mo, etc
       0 0 4 4; 
       0 0 0 10]; 
sc1=.687; %probability of responding to first dose
sc2=.687; %probability of responding to second dose
sc3=.687; %probability of responding to third dose
%res=sc1/sc2; %probability of being a "responder"
sc2n=sc2; %(1-sc2)*sc1/(1-sc1); %probability of responding to second dose given did not respond to first dose
sc3n=sc3; %(1-sc2)^2*sc1/(1-res+res*(1-sc2)^2); %probability of responding to third dose given did not respond to first or second dose

v1=zeros(tmax,al); v2=zeros(tmax,al); v3=zeros(tmax,al);
v1(:,avacc(1,targ))=[zeros(t0+tvacc,1); vcov(:,1); mean(vcov(end-52:end,1))*ones(tmax-t0-tvacc-length(vcov),1)];
if avacc(2,targ)>0
v2(:,avacc(2,targ))=[zeros(t0+tvacc,1); vcov(:,2); mean(vcov(end-52:end,2))*ones(tmax-t0-tvacc-length(vcov),1)];
end
if avacc(3,targ)>0
v3(:,avacc(3,targ))=[zeros(t0+tvacc+315,1); .81*ones(tmax-t0-tvacc-315,1)];
end

wv=0; %rate of waning of vaccine-induced immunity

% Model-predicted vaccine efficacy
VE_predict=1-((1-sc1)*(1-sc2n)*1 + (sc1*(1-sc2)+(1-sc1)*sc2n)*.03/.13 + sc1*sc2*0);
VE_predboost=1-((1-sc1)*(1-sc2n)*(1-sc3n)*1 + (sc1*(1-sc2)*(1-sc3)+(1-sc1)*sc2n*(1-sc3)+(1-sc1)*(1-sc2n)*sc3n)*.03/.13 + (sc1*sc2 + sc1*(1-sc2)*sc3 + (1-sc1)*sc2n*sc3)*0);

%Initialize vector to keep track of the number of people in each state
St0=zeros(30*al,1);
St0(1:al,1)=[N(1) zeros(1,al-1)]; %Maternal immunity
St0(al+1:2*al,1)=[0 N(2:al)-[ones(1,al-11) zeros(1,10)]]; %Susceptible_0
St0(2*al+1:3*al,1)=[0 ones(1,al-11) zeros(1,10)]; %Infectious_1 (primary) 
St0(3*al+1:4*al,1)=zeros(1,al); %Recovered_1
St0(4*al+1:5*al,1)=zeros(1,al); %Susceptible_1
St0(5*al+1:6*al,1)=zeros(1,al); %Infectious_2 (2nd time)
St0(6*al+1:7*al,1)=zeros(1,al); %Recovered_2
St0(7*al+1:8*al,1)=zeros(1,al); %Susceptible-Resistant
St0(8*al+1:9*al,1)=zeros(1,al); %Asymptomatic Infectious_3 (subsequent)
St0(9*al+1:10*al,1)=zeros(1,al); %Temp Resistant
St0(10*al+1:11*al,1)=zeros(1,al); %M2
St0(11*al+1:12*al,1)=zeros(1,al); %M3
St0(12*al+1:13*al,1)=zeros(1,al); %M4
St0(13*al+1:14*al,1)=zeros(1,al); %M5
St0(14*al+1:15*al,1)=zeros(1,al); %M6
St0(15*al+1:16*al,1)=zeros(1,al); %SV1
St0(16*al+1:17*al,1)=zeros(1,al); %IV0
St0(17*al+1:18*al,1)=zeros(1,al); %RV0
St0(18*al+1:19*al,1)=zeros(1,al); %SV1
St0(19*al+1:20*al,1)=zeros(1,al); %IV1
St0(20*al+1:21*al,1)=zeros(1,al); %RV1
St0(21*al+1:22*al,1)=zeros(1,al); %SV2
St0(22*al+1:23*al,1)=zeros(1,al); %IV2
St0(23*al+1:24*al,1)=zeros(1,al); %RV2
St0(24*al+1:25*al,1)=zeros(1,al); %MV1
St0(25*al+1:26*al,1)=zeros(1,al); %MV2
St0(26*al+1:27*al,1)=zeros(1,al); %MV3
St0(27*al+1:28*al,1)=zeros(1,al); %MV4
St0(28*al+1:29*al,1)=zeros(1,al); %MV5
St0(29*al+1:30*al,1)=zeros(1,al); %MV6

clear St lambda H

options=odeset('NonNegative',1:length(St0));
[time St]=ode45('rasisV1',1:tmax,St0,options);

if t0>0
time(1:t0,:)=[];
St(1:t0,:)=[];
end

%Force of infection
lambda=zeros(tmax-t0,al);
for t=1:size(St,1)
    lambda(t,:)=(1+b1*cos(2*pi*(time(t)-phi1)/52.18))*((St(t,2*al+1:3*al)+ri2*St(t,5*al+1:6*al)+ri3*St(t,8*al+1:9*al)+St(t,16*al+1:17*al)+ri2*St(t,19*al+1:20*al)+ri3*St(t,22*al+1:23*al))*beta)./sum(St(t,:)); %+b2*cos(2*pi*(time(t)-phi2)/26)
end

%Incident RVGE cases (C) and hospitalizations (H)
Cu=zeros(tmax-t0,al); Hu=zeros(tmax-t0,al); Cv=zeros(tmax-t0,al); Hv=zeros(tmax-t0,al); Incid=zeros(tmax-t0,al); Prev=zeros(tmax-t0,al);
for i=1:al
    Cu(:,i)=max(0,delta1(i)*St(:,al+i).*lambda(:,i)+delta2(i)*rr1*St(:,4*al+i).*lambda(:,i)+delta3(i)*rr2*St(:,7*al+i).*lambda(:,i));
    Cv(:,i)=max(0,delta1(i)*St(:,15*al+i).*lambda(:,i)+delta2(i)*rr1*St(:,18*al+i).*lambda(:,i)+delta3(i)*rr2*St(:,21*al+i).*lambda(:,i));
    Hu(:,i)=max(0,hosp1(:,i).*St(:,al+i).*lambda(:,i)+hosp2(:,i).*rr1.*St(:,4*al+i).*lambda(:,i)+hosp3(:,i).*rr2.*St(:,7*al+i).*lambda(:,i));
    Hv(:,i)=max(0,hosp1(:,i).*St(:,15*al+i).*lambda(:,i)+hosp2(:,i).*rr1.*St(:,18*al+i).*lambda(:,i)+hosp3(:,i).*rr2.*St(:,21*al+i).*lambda(:,i));
end
C=Cu+Cv;
H=Hu+Hv;

%Prevalent RV infections
Pu=St(:,2*al+1:3*al)+St(:,5*al+1:6*al)+St(:,8*al+1:9*al);
Pv=St(:,16*al+1:17*al)+St(:,19*al+1:20*al)+St(:,22*al+1:23*al);

%Total population in each age group
pop=St(:,1:al);
for i=1:29 %32 %
    pop=pop+St(:,i*al+1:(i+1)*al);
end
popagedist=sum(pop)/sum(sum(pop,2));
popagedist=[sum(popagedist(1:au5)) popagedist(au5+1:al)];

%Individuals with maternal immunity
Mu=St(:,1:al);
for i=10:14
    Mu=Mu+St(:,i*al+1:(i+1)*al);
end
Mv=St(:,24*al+1:25*al);
for i=25:29
    Mv=Mv+St(:,i*al+1:(i+1)*al);
end

%Total unvaccinated (U) and vaccinated (V) populations (by age)
U=St(:,1:al);
for i=1:14
    U=U+St(:,i*al+1:(i+1)*al);
end
V=St(:,15*al+1:16*al);
for i=16:29 %32 %
    V=V+St(:,i*al+1:(i+1)*al);
end

if v==1
    Hnovacc=H;
    Unovacc=U;
end
end

A=(H(:,1:au5)./(sum(H(:,1:au5),2)*ones(1,au5)))*avgage';
A_prevacc=A(tvacc:tvacc+48,1)'*sum(H(tvacc:tvacc+48,1:au5),2)/sum(sum(H(tvacc:tvacc+48,1:au5)));
A_postvacc=A(tvacc+49:tvf-t0,1)'*sum(H(tvacc+49:tvf-t0,1:au5),2)/sum(sum(H(tvacc+49:tvf-t0,1:au5)));
agedist_prevacc=sum(H(1:tvacc,1:au5)./sum(sum(H(1:tvacc,1:au5),2)));
agedist_postvacc=sum(H(tvacc+49:tvf-t0,1:au5))./sum(sum(H(tvacc+49:tvf-t0,1:au5),2));
agedist_novacc=sum(Hnovacc(tvacc+49:tvf-t0,1:au5))./sum(sum(Hnovacc(tvacc+49:tvf-t0,1:au5),2));
agedistV=sum(Hv(tvacc+49:tvf-t0,1:au5))./sum(sum(Hv(tvacc+49:tvf-t0,1:au5),2));
agedistU=sum(Hu(tvacc+49:tvf-t0,1:au5))./sum(sum(Hu(tvacc+49:tvf-t0,1:au5),2));

if targ==3
Hrelinf0_nobooster(:,relinf)=sum(H(:,1:12),2);
Hrelinf1_nobooster(:,relinf)=sum(H(:,13:au2),2);
Hrelinf2_nobooster(:,relinf)=sum(H(:,au2+1:au5),2);
elseif targ==4
Hrelinf0_booster(:,relinf)=sum(H(:,1:12),2);
Hrelinf1_booster(:,relinf)=sum(H(:,13:au2),2);
Hrelinf2_booster(:,relinf)=sum(H(:,au2+1:au5),2);
end

VEtot_homog(relinf,1)=1-((sum(Hv(tvacc+52:tvf-t0,4:au5))./sum(V(tvacc+52:tvf-t0,4:au5)))*agedistV(4:au5)')/((sum(Hu(tvacc+52:tvf-t0,4:au5))./sum(U(tvacc+52:tvf-t0,4:au5)))*agedistU(4:au5)');
VEage_homog(relinf,:)=1-[(sum(sum(Hv(tvacc+52:tvf-t0,4:12),2))/sum(sum(V(tvacc+52:tvf-t0,4:12),2)))/(sum(sum(Hu(tvacc+52:tvf-t0,4:12),2))/sum(sum(U(tvacc+52:tvf-t0,4:12),2)))...
    (sum(sum(Hv(tvacc+52:tvf-t0,13:au2),2))/sum(sum(V(tvacc+52:tvf-t0,13:au2),2)))/(sum(sum(Hu(tvacc+52:tvf-t0,13:au2),2))/sum(sum(U(tvacc+52:tvf-t0,13:au2),2)))...
    (sum(Hv(tvacc+52:tvf-t0,au2+1))/sum(V(tvacc+52:tvf-t0,au2+1)))/(sum(Hu(tvacc+52:tvf-t0,au2+1))/sum(U(tvacc+52:tvf-t0,au2+1)))];

IEtot_homog(relinf,1)=sum((1-(sum(Hu(tvacc+52:tvf-t0,1:au5).*Hu(tvacc+52:tvf-t0,1:au5)./U(tvacc+52:tvf-t0,1:au5),2)./sum(Hu(tvacc+52:tvf-t0,1:au5),2))./(sum(Hnovacc(tvacc+52:tvf-t0,1:au5).*Hnovacc(tvacc+52:tvf-t0,1:au5)./Unovacc(tvacc+52:tvf-t0,1:au5),2)./sum(Hnovacc(tvacc+52:tvf-t0,1:au5),2))).*sum(Hu(tvacc+52:tvf-t0,1:au5),2))/sum(sum(Hu(tvacc+52:tvf-t0,1:au5),2));
IEage_homog(relinf,:)=1-[(sum(sum(Hu(tvacc+52:tvf-t0,1:12),2))/sum(sum(U(tvacc+52:tvf-t0,1:12),2)))/(sum(sum(Hnovacc(tvacc+52:tvf-t0,1:12),2))/sum(sum(Unovacc(tvacc+52:tvf-t0,1:12),2)))...
   (sum(sum(Hu(tvacc+52:tvf-t0,13:au2),2))/sum(sum(U(tvacc+52:tvf-t0,13:au2),2)))/(sum(sum(Hnovacc(tvacc+52:tvf-t0,13:au2),2))/sum(sum(Unovacc(tvacc+52:tvf-t0,13:au2),2)))...
    (sum(sum(Hu(tvacc+52:tvf-t0,au2+1:au5),2))/sum(sum(U(tvacc+52:tvf-t0,au2+1:au5),2)))/(sum(sum(Hnovacc(tvacc+52:tvf-t0,au2+1:au5),2))/sum(sum(Unovacc(tvacc+52:tvf-t0,au2+1:au5),2)))];

OEtot_homog(relinf,1)=1-sum(sum(H(tvacc+1:tvf-t0,1:au5)))/sum(sum(Hnovacc(tvacc+1:tvf-t0,1:au5)));
OEage_homog(relinf,:)=1-[sum(sum(H(tvacc+1:tvf-t0,1:12)))/sum(sum(Hnovacc(tvacc+1:tvf-t0,1:12)))...
    sum(sum(H(tvacc+1:tvf-t0,13:au2)))/sum(sum(Hnovacc(tvacc+1:tvf-t0,13:au2)))...
    sum(sum(H(tvacc+1:tvf-t0,au2+1:au5)))/sum(sum(Hnovacc(tvacc+1:tvf-t0,au2+1:au5)))];

OEtot_obs=1-sum(sum(rotamalPV(45:end,:)))/sum(sum(Hnovacc(tvacc+1:tvf-t0,1:au5)));
OEage_obs=1-[sum(sum(rotamalPV(45:end,1:12)))/sum(sum(Hnovacc(tvacc+1:tvf-t0,1:12)))...
    sum(sum(rotamalPV(45:end,13:24)))/sum(sum(Hnovacc(tvacc+1:tvf-t0,13:au2)))...
    sum(sum(rotamalPV(45:end,25:27)))/sum(sum(Hnovacc(tvacc+1:tvf-t0,au2+1:au5)))];

if targ==3
VEtot_nobooster(relinf,1)=1-((sum(Hv(tvacc+315:tvacc+419,4:au5))./sum(V(tvacc+315:tvacc+419,4:au5)))*agedistV(4:au5)')/((sum(Hu(tvacc+315:tvacc+419,4:au5))./sum(U(tvacc+315:tvacc+419,4:au5)))*agedistU(4:au5)');
VEage_nobooster(relinf,:)=1-[(sum(sum(Hv(tvacc+315:tvacc+419,4:12),2))/sum(sum(V(tvacc+315:tvacc+419,4:12),2)))/(sum(sum(Hu(tvacc+315:tvacc+419,4:12),2))/sum(sum(U(tvacc+315:tvacc+419,4:12),2)))...
    (sum(sum(Hv(tvacc+315:tvacc+419,13:au2),2))/sum(sum(V(tvacc+315:tvacc+419,13:au2),2)))/(sum(sum(Hu(tvacc+315:tvacc+419,13:au2),2))/sum(sum(U(tvacc+315:tvacc+419,13:au2),2)))...
    (sum(Hv(tvacc+315:tvacc+419,au2+1))/sum(V(tvacc+315:tvacc+419,au2+1)))/(sum(Hu(tvacc+315:tvacc+419,au2+1))/sum(U(tvacc+315:tvacc+419,au2+1)))];
elseif targ==4
VEtot_booster(relinf,1)=1-((sum(Hv(tvacc+315:tvacc+419,4:au5))./sum(V(tvacc+315:tvacc+419,4:au5)))*agedistV(4:au5)')/((sum(Hu(tvacc+315:tvacc+419,4:au5))./sum(U(tvacc+315:tvacc+419,4:au5)))*agedistU(4:au5)');
VEage_booster(relinf,:)=1-[(sum(sum(Hv(tvacc+315:tvacc+419,4:12),2))/sum(sum(V(tvacc+315:tvacc+419,4:12),2)))/(sum(sum(Hu(tvacc+315:tvacc+419,4:12),2))/sum(sum(U(tvacc+315:tvacc+419,4:12),2)))...
    (sum(sum(Hv(tvacc+315:tvacc+419,13:au2),2))/sum(sum(V(tvacc+315:tvacc+419,13:au2),2)))/(sum(sum(Hu(tvacc+315:tvacc+419,13:au2),2))/sum(sum(U(tvacc+315:tvacc+419,13:au2),2)))...
    (sum(Hv(tvacc+315:tvacc+419,au2+1))/sum(V(tvacc+315:tvacc+419,au2+1)))/(sum(Hu(tvacc+315:tvacc+419,au2+1))/sum(U(tvacc+315:tvacc+419,au2+1)))];
end

%%
Hyr_novacc=zeros(11,4);
Hyr_nobooster=zeros(11,4);
Hyr_booster=zeros(11,4);
for y=1:11
    for j=1:length(rotamalPV)
        if malPVdate(j,1)==2011+y
            Hyr_novacc(y,:)=Hyr_novacc(y,:)+[sum(Hnovacc(tvacc+j,1:12)) sum(Hnovacc(tvacc+j,13:au2)) sum(Hnovacc(tvacc+j,au2+1:au5)) sum(Hnovacc(tvacc+j,1:au5))];
            if targ==3
                Hyr_nobooster(y,:)=Hyr_nobooster(y,:)+[Hrelinf0_nobooster(tvacc+j,1) Hrelinf1_nobooster(tvacc+j,1) Hrelinf2_nobooster(tvacc+j,1) sum(H(tvacc+j,1:au5),2)];
            elseif targ==4
                Hyr_booster(y,:)=Hyr_booster(y,:)+[Hrelinf0_booster(tvacc+j,1) Hrelinf1_booster(tvacc+j,1) Hrelinf2_booster(tvacc+j,1) sum(H(tvacc+j,1:au5),2)];
            end
        end
    end
end
   
OEyr_homog=1-Hyr_nobooster./Hyr_novacc;
OEyr_obs=1-Hyr_obsPV./Hyr_novacc;

for yr=1:9
    agedistV_byyr(yr,:)=sum(Hv(tvacc+52*yr+1:tvacc+52*(yr+1),1:au5))./sum(sum(Hv(tvacc+52*yr+1:tvacc+52*(yr+1),1:au5),2));
    agedistU_byyr(yr,:)=sum(Hu(tvacc+52*yr+1:tvacc+52*(yr+1),1:au5))./sum(sum(Hu(tvacc+52*yr+1:tvacc+52*(yr+1),1:au5),2));
    
    VEtot_byyr(yr,1)=1-((sum(Hv(tvacc+52*yr+1:tvacc+52*(yr+1),4:au5))./sum(V(tvacc+52*yr+1:tvacc+52*(yr+1),4:au5)))*agedistV_byyr(yr,4:au5)')/((sum(Hu(tvacc+52*yr+1:tvacc+52*(yr+1),4:au5))./sum(U(tvacc+52*yr+1:tvacc+52*(yr+1),4:au5)))*agedistU_byyr(yr,4:au5)');
    IEtot_byyr(yr,1)=sum((1-(sum(Hu(tvacc+52*yr+1:tvacc+52*(yr+1),1:au5).*Hu(tvacc+52*yr+1:tvacc+52*(yr+1),1:au5)./U(tvacc+52*yr+1:tvacc+52*(yr+1),1:au5),2)./sum(Hu(tvacc+52*yr+1:tvacc+52*(yr+1),1:au5),2))./(sum(Hnovacc(tvacc+52*yr+1:tvacc+52*(yr+1),1:au5).*Hnovacc(tvacc+52*yr+1:tvacc+52*(yr+1),1:au5)./Unovacc(tvacc+52*yr+1:tvacc+52*(yr+1),1:au5),2)./sum(Hnovacc(tvacc+52*yr+1:tvacc+52*(yr+1),1:au5),2))).*sum(Hu(tvacc+52*yr+1:tvacc+52*(yr+1),1:au5),2))/sum(sum(Hu(tvacc+52*yr+1:tvacc+52*(yr+1),1:au5),2));
    
    VEage_mal1_byyr(yr,:)=1-[(sum(sum(Hv(tvacc+52*yr+1:tvacc+52*(yr+1),4:12),2))/sum(sum(V(tvacc+52*yr+1:tvacc+52*(yr+1),4:12),2)))/(sum(sum(Hu(tvacc+52*yr+1:tvacc+52*(yr+1),4:12),2))/sum(sum(U(tvacc+52*yr+1:tvacc+52*(yr+1),4:12),2)))...
        (sum(sum(Hv(tvacc+52*yr+1:tvacc+52*(yr+1),13:au2),2))/sum(sum(V(tvacc+52*yr+1:tvacc+52*(yr+1),13:au2),2)))/(sum(sum(Hu(tvacc+52*yr+1:tvacc+52*(yr+1),13:au2),2))/sum(sum(U(tvacc+52*yr+1:tvacc+52*(yr+1),13:au2),2)))...
        (sum(sum(Hv(tvacc+52*yr+1:tvacc+52*(yr+1),au2+1:au5),2))/sum(sum(V(tvacc+52*yr+1:tvacc+52*(yr+1),au2+1:au5),2)))/(sum(sum(Hu(tvacc+52*yr+1:tvacc+52*(yr+1),au2+1:au5),2))/sum(sum(U(tvacc+52*yr+1:tvacc+52*(yr+1),au2+1:au5),2)))];
    IEage_mal1_byyr(yr,:)=1-[(sum(sum(Hu(tvacc+52*yr+1:tvacc+52*(yr+1),1:12),2))/sum(sum(U(tvacc+52*yr+1:tvacc+52*(yr+1),1:12),2)))/(sum(sum(Hnovacc(tvacc+52*yr+1:tvacc+52*(yr+1),1:12),2))/sum(sum(Unovacc(tvacc+52*yr+1:tvacc+52*(yr+1),1:12),2)))...
       (sum(sum(Hu(tvacc+52*yr+1:tvacc+52*(yr+1),13:au2),2))/sum(sum(U(tvacc+52*yr+1:tvacc+52*(yr+1),13:au2),2)))/(sum(sum(Hnovacc(tvacc+52*yr+1:tvacc+52*(yr+1),13:au2),2))/sum(sum(Unovacc(tvacc+52*yr+1:tvacc+52*(yr+1),13:au2),2)))...
       (sum(sum(Hu(tvacc+52*yr+1:tvacc+52*(yr+1),au2+1:au5),2))/sum(sum(U(tvacc+52*yr+1:tvacc+52*(yr+1),au2+1:au5),2)))/(sum(sum(Hnovacc(tvacc+52*yr+1:tvacc+52*(yr+1),au2+1:au5),2))/sum(sum(Unovacc(tvacc+52*yr+1:tvacc+52*(yr+1),au2+1:au5),2)))];
end

end
end

%%
datesim(1:t0,:)=[];

%%
figure
subplot(2,3,1:3)
hold on
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 10],'--k')
plot(malweekPV,sum(rotamalPV,2),'Color',[.7 .7 .7])
plot(datesim(tvacc:end),sum(Hnovacc(tvacc:end,1:au5),2),'b','LineWidth',2)
plot(datesim(tvacc:end),sum(H(tvacc:end,1:au5),2),'r','LineWidth',2)
datetick('x','mmm-yy')
ylabel('Number of RVGE cases (per week)')

subplot(2,3,4)
hold on
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 10],'--k')
plot(malweekPV,sum(rotamalPV(:,1:12),2),'Color',[.7 .7 .7])
plot(datesim(tvacc:end),sum(Hnovacc(tvacc:end,1:12),2),'b','LineWidth',2)
plot(datesim(tvacc:end),sum(H(tvacc:end,1:12),2),'r','LineWidth',2)
datetick('x','mmm-yy')
ylabel('Number of RVGE cases (per week)')
title('<1 yr olds')

subplot(2,3,5)
hold on
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 6],'--k')
plot(malweekPV,sum(rotamalPV(:,13:au2),2),'Color',[.7 .7 .7])
plot(datesim(tvacc:end),sum(Hnovacc(tvacc:end,13:au2),2),'b','LineWidth',2)
plot(datesim(tvacc:end),sum(H(tvacc:end,13:au2),2),'r','LineWidth',2)
datetick('x','mmm-yy')
ylabel('Number of RVGE cases (per week)')
title('1-<2 yr olds')

subplot(2,3,6)
hold on
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 4],'--k')
plot(malweekPV,sum(rotamalPV(:,au2+1:au5),2),'Color',[.7 .7 .7])
plot(datesim(tvacc:end),sum(Hnovacc(tvacc:end,au2+1:au5),2),'b','LineWidth',2)
plot(datesim(tvacc:end),sum(H(tvacc:end,au2+1:au5),2),'r','LineWidth',2)
datetick('x','mmm-yy')
ylabel('Number of RVGE cases (per week)')
title('2-<5 yr olds')

%%
%figure 
%hold on
%plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 10],'--k')
%plot(datenum(malPVdate),[sum(rotamalPV(:,1:12),2) sum(rotamalPV(:,13:24),2) sum(rotamalPV(:,25:end),2)]); %,'k','LineWidth',1)
%plot(datesim(tvacc:end),[sum(H(tvacc:end,1:12),2) sum(H(tvacc:end,13:au2),2) sum(H(tvacc:end,au2+1:au5),2)],'LineWidth',2)
%plot(datesim(tvacc:end),[Hrelinf0_homog(tvacc:end,relinf) Hrelinf1_homog(tvacc:end,relinf) Hrelinf2_homog(tvacc:end,relinf)],'--','LineWidth',1.5)
%plot(datesim(tvacc:tvf-t0),[Hrelinf0_nobooster(tvacc:tvf-t0,relinf) Hrelinf1_nobooster(tvacc:tvf-t0,relinf) Hrelinf2_nobooster(tvacc:tvf-t0,relinf)],'LineWidth',1.5)
%plot(datesim(tvacc:tvacc+233),[Hrelinf0_nobooster_noimm(tvacc:tvacc+233,relinf) Hrelinf1_nobooster_noimm(tvacc:tvacc+233,relinf) Hrelinf2_nobooster_noimm(tvacc:tvacc+233,relinf)],'LineWidth',1.5)
%plot(datesim(tvacc:end),[Hrelinf0_booster(tvacc:end,relinf) Hrelinf1_booster(tvacc:end,relinf) Hrelinf2_booster(tvacc:end,relinf)],'LineWidth',1.5)
%datetick('x','mmm-yy')
%ylabel('Number of RVGE cases (per week)')

%%
%figure; subplot(2,1,1); hold on;
%plot(datenum([blantyre_pop(:,1:2) ones(length(blantyre_pop),1) zeros(length(blantyre_pop),3)]),blantyre_pop(:,5),'xb')
%plot(datenum([(1998:2014)' 7*ones(length(blantyre_pop2),1) ones(length(blantyre_pop2),1) zeros(length(blantyre_pop2),3)]),sum(blantyre_pop2,2),'ok')
%plot(datesim,sum(pop,2),'r'); 
%datetick('x','yyyy')
%ylabel('Population of Blantyre district'); legend('Data','Model')
%subplot(2,1,2); y1=bar([malawi_agedist(1:16)' popagedist']);
%set(y1(1),'FaceColor','b'); set(y1(2),'FaceColor','r'); ylim([0 .31])
%ylabel('Proportion of population'); xlabel('Age group (years)'); legend('Data','Model')