%%% MODEL 4 ASSUMING A NEW VACCINE WITH IMPROVED IMMUNOGENICITY IS INTRODUCED IN JANUARY 2018 

% NOTE: You need to run rotagemodelMVw_nr first

clear global all
global B wm wi1 wi2 u um beta b1 phi1 d1 d2 rr1 rr2 ri2 ri3 al v1 v2 v3 sc1 sc2 sc3 sc2n sc3n wv reintro wA; 

%%% FIXED POPULATION PARAMETERS %%%
age=[0:1/12:23/12 2:4 5:5:75];
al=length(age); %Number of age groups
au2=24; %Number of age groups <2 yrs old
au5=27; %Number of age groups <5 yrs old
agep=[(1/12)*malpop(1)*ones(1,12) (1/48)*malpop(2)*ones(1,12) (1/4)*malpop(2)*ones(1,3) malpop(3:17).*ones(1,15)];
agep=agep/sum(agep); %Proportion of pop'n in each age group 
u=[1/4.3*ones(1,24) 1/52*ones(1,3) 1/(52*5)*ones(1,14) 1/(52*25)]; %Rate of aging out of each age group

t0=round(52.18*30); %burn-in period
tvacc=757; %start of post-vaccination data
tvf=t0+tvacc+length(rotamalVS3); %end of post-vaccination data
tmax=round(52.18*55); %length of simulation
datepop=[1950 1 1 0 0 0; (1957:5:2037)' 7*ones(17,1) ones(17,1) zeros(17,3)]; %Date corresponding to census estimates
datesim=(datenum([1967 7 2 0 0 0]):7:datenum([2022 7 1 0 0 0]))'; %Date of simulation (weekly time step)

N1=870000; %Blantyre population in 2000
N=(N1/exp(.04*30))*agep; %Adjust initial population size for exponential growth occurring over the burn-in period
B=interp1(datenum(datepop),malawi_cbr/1000,datesim); %Interpolate weekly birth rate
if t0>0 && length(B)<tmax
    B=[B; B(end)*ones(tmax-length(B),1)]; %Assume birth rate is equal to last available data if simulating for longer
end
B=[B zeros(tmax,al-1)]; %Birth rate for all age groups (only non-zero for youngest age group)
um=.010*ones(1,al); %Crude death/immigration rate (chosen to maintain population size consistent with Blantyre population growth)
um=log(1+um)/52.18; %Convert from annual to weekly rate

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
pfit_mal=[78.77 47.39 34.01 26.59 21.81; %R0 Best-fit parameters for different values of ri3 (=0.1 to 0.5)
    0.1740 0.1450 0.1345 0.1301 0.1275; %Seasonal amplitude
    6.856 7.730 8.084 8.203 8.314; %Seasonal offset
    0.0171 0.0171 0.0171 0.0171 0.0171]; %Mean reporting fraction

pvacc_mal=[0.4957; 0.7882; 0.0222]; %Best-fit vaccination parameters (did not vary by ri3)
 
for relinf=1
    
ri3=relinf/10; %relative infectiousness of asymptomatic infection
pars=pfit_mal(:,relinf); 
pvacc=pvacc_mal;

ptrans=pars(1); %beta0~R0
b1=pars(2); %seasonal forcing 
phi1=pars(3); %seasonal offset
ar1=1; %relative risk of infection for <1 yr olds
ar2=1; %relative risk of infection for 1-2 yr olds 

%Variation in reporting effort over time (by age)
repeff=[(negtestage_movavg(:,1)/mean(negtest_age(:,1)))*ones(1,12) (negtestage_movavg(:,2)/mean(negtest_age(:,2)))*ones(1,au2-12) (negtestage_movavg(:,3)/mean(negtest_age(:,3)))*ones(1,au5-au2) ones(522,al-au5);...
    ones(26,al); (negtestage08_movavg(:,1)/mean(negtest_age(:,1)))*ones(1,12) (negtestage08_movavg(:,2)/mean(negtest_age(:,2)))*ones(1,au2-12) (negtestage08_movavg(:,3)/mean(negtest_age(:,3)))*ones(1,au5-au2) ones(105,al-au5);...
    ones(104,al); (negtestVS3_movavg(:,1)/mean(negtest_age(:,1)))*ones(1,12) (negtestVS3_movavg(:,2)/mean(negtest_age(:,2)))*ones(1,au2-12) (negtestVS3_movavg(:,3)/mean(negtest_age(:,3)))*ones(1,au5-au2) ones(295,al-au5);...
    (negtestVS3_movavg(end,1)/mean(negtest_age(:,1)))*ones(tmax-tvf,12) (negtestVS3_movavg(end,2)/mean(negtest_age(:,2)))*ones(tmax-tvf,au2-12) (negtestVS3_movavg(end,3)/mean(negtest_age(:,3)))*ones(tmax-tvf,au5-au2) ones(tmax-tvf,al-au5)];

%c1=100*ones(al); %Homogeneous mixing
c2=100*[ar1*ones(al,12) ar2*ones(al,au2-12) ones(al,al-au2)]; %Allows for different risk of infection for <1 and 1-2 yr olds (if ar1, ar2 =/=1)

beta=(ptrans/100)*c2; %transmission rate matrix

immunity=[0.13 0.063; 0.03 0.077]; %Probability of moderate-to-severe RVGE upon first and second infection according to Velazquez et al (NEJM 1996) and Gladstone et al (NEJM 2011)
im=1; %Velazquez et al immunity parameters provide better fit to data

h=pars(4); %proportion of moderate-to-severe diarrhea cases hospitalized (at QECH) and reported
hosp1=immunity(1,im)*h*repeff; %reporting rate over time for first infections
hosp2=immunity(2,im)*h*repeff; %reporting rate over time for second infections
hosp3=zeros(tmax-t0,al); %reporting rate over time for subsequent infections (~0 when estimated)
delta1=0.41*h*ones(1,al); %proportion of primary infections that are symptomatic (any RVGE)
delta2=0.35*h*ones(1,al); %proportion of secondary infections that are symptomatic (any RVGE)
delta3=0.21*h*ones(1,al); %rate of detection of subsequent infection (any RVGE)

reintro=0; %Allows for possible constant low background risk of infection if =/=0

R0=max(eig(dur*beta.*(ones(al,1)*N))); %R0 = max eigenvalue of the next-generation matrix


%%% VACCINATION PARAMETERS

for targ=3 %1=all-or-nothing protection following first dose; 2=all-or-nothing protection following second dose; 3=incremental protection following each dose
for v=2 %[1 10]
%c(v)=(v-1)/10; %[.5 .7 .9]; %(v-11+90)/100; %
if v==1
    vcov=zeros(length(vcov_mavg),2);
else
    vcov=vcov_mavg;
end
avacc=[1 3 3 3; %1=birth, 2=1mo, 3=2mo, etc
       0 0 4 4; 
       0 0 0 10]; %[1 5 9]

v1=zeros(tmax,al); v2=zeros(tmax,al); v3=zeros(tmax,al);
%v1(:,avacc(1,targ))=[zeros(t0+tvacc,1); ones(tmax-t0-tvacc,1)]*c(v);
v1(:,avacc(1,targ))=[zeros(t0+tvacc,1); vcov(:,1); mean(vcov(end-52:end,1))*ones(tmax-t0-tvacc-length(vcov),1)];
if avacc(2,targ)>0
%v2(:,avacc(2,targ))=[zeros(t0+tvacc+4,1); ones(tmax-t0-tvacc-4,1)*c(v)];
v2(:,avacc(2,targ))=[zeros(t0+tvacc+4,1); vcov(:,2); mean(vcov(end-52:end,2))*ones(tmax-t0-tvacc-4-length(vcov),1)];
end
if avacc(3,targ)>0
%v3(:,avacc(3,targ))=[zeros(t0+tvacc+267,1); .87*ones(tmax-t0-tvacc-267,1)];
v3(:,avacc(3,targ))=[zeros(t0+tvacc+315,1); .81*ones(tmax-t0-tvacc-315,1)];
end

wv=pvacc(3); %vpfit3r_mal; %1/(52*prctile(keep_pars_Vnr(:,3,1),50)); %vpfit2_mal(3); %0; %1/78; %rate of waning of vaccine-induced immunity

for im=1:26
sc1=[pvacc(1)*ones(t0+tvacc+315,1); (pvacc(1)+(im-1)/100)*ones(tmax-t0-tvacc-315,1)];
sc2=pvacc(2)*ones(tmax,1); %[pvacc(2)*ones(t0+tvacc+315,1); (pvacc(2)+(im-1)/100)*ones(tmax-t0-tvacc-315,1)];
sc2n=(1-sc2).*sc1./(1-sc1).*ones(tmax,1); %*[ones(t0+tvacc+315,1); (1+im/100)*ones(tmax-t0-tvacc-315,1)]; %sc2; %probability of responding to second dose given did not respond to first dose
res=sc1(1)/sc2(1); %probability of being a "responder"
sc3=pvacc(2); 
sc3n=(1-sc2(1))^2*sc1(1)/(1-res+res*(1-sc2(1))^2); %sc3; %probability of responding to third dose given did not respond to first or second dose


%Initialize vector to keep track of the number of people in each state
St0=zeros(33*al,1);
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
St0(10*al+1:11*al,1)=zeros(1,al); %Maternal immunity
St0(11*al+1:12*al,1)=zeros(1,al); %
St0(12*al+1:13*al,1)=zeros(1,al); %
St0(13*al+1:14*al,1)=zeros(1,al); %
St0(14*al+1:15*al,1)=zeros(1,al); %
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
[time St]=ode45('rasisV1wim',1:tmax,St0,options);

if t0>0
time(1:t0,:)=[];
St(1:t0,:)=[];
end

lambda=zeros(tmax-t0,al);
for t=1:size(St,1)
    lambda(t,:)=(1+b1*cos(2*pi*(time(t)-phi1)/52.18))*((St(t,2*al+1:3*al)+ri2*St(t,5*al+1:6*al)+ri3*St(t,8*al+1:9*al)+St(t,16*al+1:17*al)+ri2*St(t,19*al+1:20*al)+ri3*St(t,22*al+1:23*al))*beta)./sum(St(t,:)); %+b2*cos(2*pi*(time(t)-phi2)/26)
end

%Incident RVGE cases (C) and hospitalizations (H)
Hu=zeros(tmax-t0,al); Hv=zeros(tmax-t0,al); 
for i=1:al
    Hu(:,i)=max(0,hosp1(:,i).*St(:,al+i).*lambda(:,i)+hosp2(:,i).*rr1.*St(:,4*al+i).*lambda(:,i)+hosp3(:,i).*rr2.*St(:,7*al+i).*lambda(:,i));
    Hv(:,i)=max(0,hosp1(:,i).*St(:,15*al+i).*lambda(:,i)+hosp2(:,i).*rr1.*St(:,18*al+i).*lambda(:,i)+hosp3(:,i).*rr2.*St(:,21*al+i).*lambda(:,i));
end
H=Hu+Hv;

%Prevalent RV infections
Pu=St(:,2*al+1:3*al)+St(:,5*al+1:6*al)+St(:,8*al+1:9*al);
Pv=St(:,16*al+1:17*al)+St(:,19*al+1:20*al)+St(:,22*al+1:23*al);

%Total population in each age group
pop=St(:,1:al);
for i=1:32 %29 %
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
for i=16:32 
    V=V+St(:,i*al+1:(i+1)*al);
end

if v==1
    Hnovacc_wane=H;
    Unovacc_wane=U;
end

%% Calculate summary measures of vaccine impact

%Calculate the number of annual RVGE cases in each age group for the 10
%years following vaccine introduction
for y=1:10
    Hyr0_imprsc4(y,im)=sum(sum(H(tvacc+round(52.18*(y-1)):tvacc+round(52.18*(y-1))+51,1:4)));
    Hyr1_imprsc4(y,im)=sum(sum(H(tvacc+round(52.18*(y-1)):tvacc+round(52.18*(y-1))+51,5:12)));
    Hyr2_imprsc4(y,im)=sum(sum(H(tvacc+round(52.18*(y-1)):tvacc+round(52.18*(y-1))+51,13:au2)));
    Hyr3_imprsc4(y,im)=sum(sum(H(tvacc+round(52.18*(y-1)):tvacc+round(52.18*(y-1))+51,au2+1:au5)));
    HyrU5_imprsc4(y,im)=sum(sum(H(tvacc+round(52.18*(y-1)):tvacc+round(52.18*(y-1))+51,1:au5)));
end

%Calculate the overall effectiveness for all children <5yo and by age group
OE_imprsc(im,4)=1-sum(HyrU5_imprsc4(7:9,im))./sum(sum(Hyr_wane_nr(7:9,:,1),2));
OEage_imprsc4(im,:)=1-[(sum(Hyr0_imprsc4(7:9,im)+sum(Hyr1_imprsc4(7:9,im))))./sum(sum(Hyr_wane_nr(7:9,1:2,1),2)) sum(Hyr2_imprsc4(7:9,im))./sum(Hyr_wane_nr(7:9,3,1)) sum(Hyr3_imprsc4(7:9,im))./sum(Hyr_wane_nr(7:9,4,1))];

end
end
end
end

