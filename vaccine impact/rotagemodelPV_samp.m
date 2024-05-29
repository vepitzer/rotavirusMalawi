% Simulates base-case rotavirus model and estimates vaccine effects while
% sampling from the posterior distribution of estimated parameters
% NOTE THAT THIS CODE TAKES ~1 DAY TO RUN

clear global all
global B wm wi1 wi2 u um beta b1 phi1 d1 d2 rr1 rr2 ri2 ri3 al v1 v2 v3 sc1 sc2 sc3 sc2n sc3n wv reintro wA; 

%load('malawidata.mat')
load('malawi_mcmc_output1.mat','keep_pars_1')
load('malawi_mcmc_outputV_bc.mat','keep_pars_V_bc')

%%% FIXED POPULATION PARAMETERS %%%
age=[0:1/12:23/12 2:4 5:5:75];
al=length(age); %Number of age groups
au2=24; %13; %
au5=27; %16; %
agep=[(1/12)*malpop(1)*ones(1,12) (1/48)*malpop(2)*ones(1,12) (1/4)*malpop(2)*ones(1,3) malpop(3:17).*ones(1,15)];
agep=agep/sum(agep); %Proportion of pop'n in each age group 
u=[1/4.3*ones(1,24) 1/52*ones(1,3) 1/(52*5)*ones(1,14) 1/(52*25)]; %Rate of aging out of each age group

tmax=round(52.18*55); 
t0=round(52.18*30); 
tvacc=757;
tvf=t0+tvacc+length(rotamalPV);
datesim=(datenum([1967 7 2 0 0 0]):7:datenum([2022 7 1 0 0 0]))';
datepop=[1950 1 1 0 0 0; (1957:5:2037)' 7*ones(17,1) ones(17,1) zeros(17,3)]; %Date corresponding to census estimates

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
for relinf=1%:5
for samp=1:100%00
    
pars=keep_pars_1(10*samp,:,1); 
ri3=relinf/10; 

ptrans=pars(1); %beta0~R0
b1=pars(2); %seasonal forcing 
phi1=pars(3); %seasonal offset
ar1=1; %relative risk of infection for <1 yr olds
ar2=1; %relative risk of infection for 1-2 yr olds 

%Variation in reporting effort over time (by age)
repeff=[(negtestage_movavg(:,1)/mean(negtest_age(:,1)))*ones(1,12) (negtestage_movavg(:,2)/mean(negtest_age(:,2)))*ones(1,au2-12) (negtestage_movavg(:,3)/mean(negtest_age(:,3)))*ones(1,au5-au2) ones(522,al-au5);...
    ones(26,al);... 
    (negtestage08_movavg(:,1)/mean(negtest_age(:,1)))*ones(1,12) (negtestage08_movavg(:,2)/mean(negtest_age(:,2)))*ones(1,au2-12) (negtestage08_movavg(:,3)/mean(negtest_age(:,3)))*ones(1,au5-au2) ones(105,al-au5);...
    ones(104,al);...
    %(negtestVS3_movavg(:,1)/mean(negtest_age(:,1)))*ones(1,12) (negtestVS3_movavg(:,2)/mean(negtest_age(:,2)))*ones(1,au2-12) (negtestVS3_movavg(:,3)/mean(negtest_age(:,3)))*ones(1,au5-au2) ones(295,al-au5);...
    %(negtestVS3_movavg(end,1)/mean(negtest_age(:,1)))*ones(tmax-tvf,12) (negtestVS3_movavg(end,2)/mean(negtest_age(:,2)))*ones(tmax-tvf,au2-12) (negtestVS3_movavg(end,3)/mean(negtest_age(:,3)))*ones(tmax-tvf,au5-au2) ones(tmax-tvf,al-au5)];
    repeffPV ones(547,al-au5);...
    ones(tmax-tvf,1)*repeffPV(end,:) ones(tmax-tvf,al-au5)];

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

avacc=[1 3 3 3; %1=birth, 2=1mo, 3=2mo, etc
       0 0 4 4; 
       0 0 0 10]; %[1 5 9]
sc1=betarnd(24.7,11.3); 
sc2=betarnd(24.7,11.3); 
sc3=betarnd(24.7,11.3); 
%sc1=keep_pars_V_bc(10*samp-5000*floor((samp-.5)/500),1,ceil((samp-.5)/500)); %.528; %.687; %vpfit1_mal(1); %.568; %betarnd(24.7,11.3); %.63; %betarnd(192,117); %
%sc2=keep_pars_V_bc(10*samp-5000*floor((samp-.5)/500),2,ceil((samp-.5)/500)); %.895; %.687; %vpfit1_mal(2); %.736; %betarnd(24.7,11.3); %%.63; %betarnd(192,117); %
%sc3=keep_pars_V_bc(10*samp-5000*floor((samp-.5)/500),2,ceil((samp-.5)/500)); %.895; %.687; %vpfit1_mal(2); %.736; %betarnd(55,25); %
res=sc1/sc2; %probability of being a "responder"
sc2n=sc2; %(1-sc2)*sc1/(1-sc1); %probability of responding to second dose given did not respond to first dose
sc3n=sc3; %(1-sc2)^2*sc1/(1-res+res*(1-sc2)^2); %probability of responding to third dose given did not respond to first or second dose

wv=0; %rate of waning of vaccine-induced immunity

for targ=3:4 %1=all-or-nothing protection following first dose; 2=all-or-nothing protection following second dose; 3=incremental protection following each dose
for v=1:2 
if v==1
    vcov=zeros(length(vcovPV_mavg1),2);
else
    vcov=vcovPV_mavg1;
end

v1=zeros(tmax,al); v2=zeros(tmax,al); v3=zeros(tmax,al);
v1(:,avacc(1,targ))=[zeros(t0+tvacc,1); vcov(:,1); mean(vcov(end-52:end,1))*ones(tmax-t0-tvacc-length(vcov),1)];
if avacc(2,targ)>0
v2(:,avacc(2,targ))=[zeros(t0+tvacc,1); vcov(:,2); mean(vcov(end-52:end,2))*ones(tmax-t0-tvacc-length(vcov),1)];
end
if avacc(3,targ)>0
%v3(:,avacc(3,targ))=[zeros(t0+tvacc+267,1); .87*ones(tmax-t0-tvacc-267,1)];
v3(:,avacc(3,targ))=[zeros(t0+tvacc+315,1); .81*ones(tmax-t0-tvacc-315,1)];
end

%Initialize vector to keep track of the number of people in each state
St0=zeros(30*al,1);
St0(1:al,1)=[N(1) zeros(1,al-1)]; %Maternal immunity
St0(al+1:2*al,1)=[0 N(2:al)-[ones(1,al-11) zeros(1,10)]]; %Susceptible_0
St0(2*al+1:3*al,1)=[0 ones(1,al-11) zeros(1,10)]; %Infectious_1 (primary) 


clear St lambda H

options=odeset('NonNegative',1:length(St0));
[time St]=ode45('rasisV1',1:tmax,St0,options);

if t0>0
time(1:t0,:)=[];
St(1:t0,:)=[];
end

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
for i=16:29 
    V=V+St(:,i*al+1:(i+1)*al);
end

agedist_prevacc=sum(H(1:tvacc,1:au5)./sum(sum(H(1:tvacc,1:au5),2)));
agedist_postvacc=sum(H(tvacc+52:tvf-t0,1:au5))./sum(sum(H(tvacc+52:tvf-t0,1:au5),2));
agedistV=sum(Hv(tvacc+52:tvf-t0,1:au5))./sum(sum(Hv(tvacc+52:tvf-t0,1:au5),2));
agedistU=sum(Hu(tvacc+52:tvf-t0,1:au5))./sum(sum(Hu(tvacc+52:tvf-t0,1:au5),2));

if v==1
    Hnovacc=H;
    Unovacc=U;
end
end

Hsamp0_novacc(:,samp)=sum(Hnovacc(:,1:12),2);
Hsamp1_novacc(:,samp)=sum(Hnovacc(:,13:au2),2);
Hsamp2_novacc(:,samp)=sum(Hnovacc(:,au2+1:au5),2);
if targ==3
Hsamp0(:,samp)=sum(H(:,1:12),2);
Hsamp1(:,samp)=sum(H(:,13:au2),2);
Hsamp2(:,samp)=sum(H(:,au2+1:au5),2);
elseif targ==4
Hsamp0_boost(:,samp)=sum(H(:,1:12),2);
Hsamp1_boost(:,samp)=sum(H(:,13:au2),2);
Hsamp2_boost(:,samp)=sum(H(:,au2+1:au5),2);
end    


VEtot_samp(samp,1)=1-((sum(Hv(tvacc+52:tvf-t0,4:au5))./sum(V(tvacc+52:tvf-t0,4:au5)))*agedistV(4:au5)')/((sum(Hu(tvacc+52:tvf-t0,4:au5))./sum(U(tvacc+52:tvf-t0,4:au5)))*agedistU(4:au5)');
VEage_samp(samp,:)=1-[(sum(sum(Hv(tvacc+52:tvf-t0,4:12),2))/sum(sum(V(tvacc+52:tvf-t0,4:12),2)))/(sum(sum(Hu(tvacc+52:tvf-t0,4:12),2))/sum(sum(U(tvacc+52:tvf-t0,4:12),2)))...
    (sum(sum(Hv(tvacc+52:tvf-t0,13:au2),2))/sum(sum(V(tvacc+52:tvf-t0,13:au2),2)))/(sum(sum(Hu(tvacc+52:tvf-t0,13:au2),2))/sum(sum(U(tvacc+52:tvf-t0,13:au2),2)))...
    (sum(Hv(tvacc+52:tvf-t0,au2+1))/sum(V(tvacc+52:tvf-t0,au2+1)))/(sum(Hu(tvacc+52:tvf-t0,au2+1))/sum(U(tvacc+52:tvf-t0,au2+1)))];

%IEtot_samp(samp,1)=1-((sum(Hu(tvacc+52:tvf-t0,1:au5))./sum(U(tvacc+52:tvf-t0,1:au5)))*agedistU(1:au5)')/((sum(Hnovacc(tvacc+52:tvf-t0,1:au5))./sum(Unovacc(tvacc+52:tvf-t0,1:au5)))*agedist_prevacc(1:au5)');
IEtot_samp(samp,1)=sum((1-(sum(Hu(tvacc+52:tvf-t0,1:au5).*Hu(tvacc+52:tvf-t0,1:au5)./U(tvacc+52:tvf-t0,1:au5),2)./sum(Hu(tvacc+52:tvf-t0,1:au5),2))./(sum(Hnovacc(tvacc+52:tvf-t0,1:au5).*Hnovacc(tvacc+52:tvf-t0,1:au5)./Unovacc(tvacc+52:tvf-t0,1:au5),2)./sum(Hnovacc(tvacc+52:tvf-t0,1:au5),2))).*sum(Hu(tvacc+52:tvf-t0,1:au5),2))/sum(sum(Hu(tvacc+52:tvf-t0,1:au5),2));
IEage_samp(samp,:)=1-[(sum(sum(Hu(tvacc+52:tvf-t0,1:12),2))/sum(sum(U(tvacc+52:tvf-t0,1:12),2)))/(sum(sum(Hnovacc(tvacc+52:tvf-t0,1:12),2))/sum(sum(Unovacc(tvacc+52:tvf-t0,1:12),2)))...
    (sum(sum(Hu(tvacc+52:tvf-t0,13:au2),2))/sum(sum(U(tvacc+52:tvf-t0,13:au2),2)))/(sum(sum(Hnovacc(tvacc+52:tvf-t0,13:au2),2))/sum(sum(Unovacc(tvacc+52:tvf-t0,13:au2),2)))...
   (sum(sum(Hu(tvacc+52:tvf-t0,au2+1:au5),2))/sum(sum(U(tvacc+52:tvf-t0,au2+1:au5),2)))/(sum(sum(Hnovacc(tvacc+52:tvf-t0,au2+1:au5),2))/sum(sum(Unovacc(tvacc+52:tvf-t0,au2+1:au5),2)))];

VEtot_boostsamp(samp,1)=1-((sum(Hv(tvacc+315:tvacc+419,4:au5))./sum(V(tvacc+315:tvacc+419,4:au5)))*agedistV(4:au5)')/((sum(Hu(tvacc+315:tvacc+419,4:au5))./sum(U(tvacc+315:tvacc+419,4:au5)))*agedistU(4:au5)');
VEage_boostsamp(samp,:)=1-[(sum(sum(Hv(tvacc+315:tvacc+419,4:12),2))/sum(sum(V(tvacc+315:tvacc+419,4:12),2)))/(sum(sum(Hu(tvacc+315:tvacc+419,4:12),2))/sum(sum(U(tvacc+315:tvacc+419,4:12),2)))...
    (sum(sum(Hv(tvacc+315:tvacc+419,13:au2),2))/sum(sum(V(tvacc+315:tvacc+419,13:au2),2)))/(sum(sum(Hu(tvacc+315:tvacc+419,13:au2),2))/sum(sum(U(tvacc+315:tvacc+419,13:au2),2)))...
    (sum(Hv(tvacc+315:tvacc+419,au2+1))/sum(V(tvacc+315:tvacc+419,au2+1)))/(sum(Hu(tvacc+315:tvacc+419,au2+1))/sum(U(tvacc+315:tvacc+419,au2+1)))];

%pars=[ptrans; b1; phi1; al; 1; 1; h];
%sampled_logLL_Vnofit(samp,1)=rotamalV_LL(St(tvacc+1:tvf-t0,:),pars,rotamalVS3_vacc,rotamalVS3_unvacc,repeffV3,time(tvacc+1:tvf-t0,1));

end
end
end

%%
Hsamp_novacc=Hsamp0_novacc+Hsamp1_novacc+Hsamp2_wane;
Hsamp=Hsamp0+Hsamp1+Hsamp2;

Hobserr=zeros(tmax-t0,10000);
Hobserr0=zeros(tmax-t0,10000);
Hobserr1=zeros(tmax-t0,10000);
Hobserr2=zeros(tmax-t0,10000);
for i=1:100
    Hobserr(:,100*(i-1)+1:100*i)=poissrnd(Hsamp);
    Hobserr0(:,100*(i-1)+1:100*i)=poissrnd(Hsamp0);
    Hobserr1(:,100*(i-1)+1:100*i)=poissrnd(Hsamp1);
    Hobserr2(:,100*(i-1)+1:100*i)=poissrnd(Hsamp2);
end
Hpi=[prctile(Hobserr,2.5,2) prctile(Hobserr,97.5,2)];
Hpi0=[prctile(Hobserr0,2.5,2) prctile(Hobserr0,97.5,2)];
Hpi1=[prctile(Hobserr1,2.5,2) prctile(Hobserr1,97.5,2)];
Hpi2=[prctile(Hobserr2,2.5,2) prctile(Hobserr2,97.5,2)];

%%
outside_pi=zeros(520,4);
for i=1:431 
    if sum(rotamalPV(i,:),2)<Hpi(tvacc+i,1) || sum(rotamalPV(i,:),2)>Hpi(tvacc+i,2)
        outside_pi(i,1)=1;
    end
    if sum(rotamalPV(i,1:12),2)<Hpi0(tvacc+i,1) || sum(rotamalPV(i,1:12),2)>Hpi0(tvacc+i,2)
        outside_pi(i,2)=1;
    end
    if sum(rotamalPV(i,13:au2),2)<Hpi1(tvacc+i,1) || sum(rotamalPV(i,13:au2),2)>Hpi1(tvacc+i,2)
        outside_pi(i,3)=1;
    end
    if sum(rotamalPV(i,au2+1:au5),2)<Hpi2(tvacc+i,1) || sum(rotamalPV(i,au2+1:au5),2)>Hpi2(tvacc+i,2)
        outside_pi(i,4)=1;
    end
end
for i=459:547
    if sum(rotamalPV(i,:),2)<Hpi(tvacc+i,1) || sum(rotamalPV(i,:),2)>Hpi(tvacc+i,2)
        outside_pi(i-27,1)=1;
    end
    if sum(rotamalPV(i,1:12),2)<Hpi0(tvacc+i,1) || sum(rotamalPV(i,1:12),2)>Hpi0(tvacc+i,2)
        outside_pi(i-27,2)=1;
    end
    if sum(rotamalPV(i,13:au2),2)<Hpi1(tvacc+i,1) || sum(rotamalPV(i,13:au2),2)>Hpi1(tvacc+i,2)
        outside_pi(i-27,3)=1;
    end
    if sum(rotamalPV(i,au2+1:au5),2)<Hpi2(tvacc+i,1) || sum(rotamalPV(i,au2+1:au5),2)>Hpi2(tvacc+i,2)
        outside_pi(i-27,4)=1;
    end
end
mean(outside_pi)
[mean(outside_pi(1:295,:)); mean(outside_pi(296:end,:))]

%%
%r_samp=zeros(1,100);
%r_pval=zeros(1,100);
for i=1:100
    [r_samp(i),r_pval(i)]=corr(sum(rotamalPV([1:431 459:end],:),2),Hsamp([tvacc+1:tvacc+431 tvacc+459:(end-1)],i),'Type','Spearman');
    r_samp_early(i)=corr(sum(rotamalPV(1:295,:),2),Hsamp(tvacc+1:tvacc+295,i),'Type','Spearman');
    r_samp_late(i)=corr(sum(rotamalPV([296:431 459:end],:),2),Hsamp([tvacc+296:tvacc+431 tvacc+459:(end-1)],i),'Type','Spearman');
    [r_samp0(i),r_pval0(i)]=corr(sum(rotamalPV([1:431 459:end],:),2),Hsamp0([tvacc+1:tvacc+431 tvacc+459:(end-1)],i),'Type','Spearman');
    [r_samp1(i),r_pval1(i)]=corr(sum(rotamalPV([1:431 459:end],:),2),Hsamp1([tvacc+1:tvacc+431 tvacc+459:(end-1)],i),'Type','Spearman');
    [r_samp2(i),r_pval2(i)]=corr(sum(rotamalPV([1:431 459:end],:),2),Hsamp2([tvacc+1:tvacc+431 tvacc+459:(end-1)],i),'Type','Spearman');
end
mean([r_samp' r_pval'])
mean([r_samp_early' r_samp_late'])
%prctile(r_samp',[2.5 97.5])
%mean([r_samp0' r_samp1' r_samp2'])

%%
save('malawi_pvsamp_output.mat','Hsamp0','Hsamp1','Hsamp2','Hsamp0_boost','Hsamp1_boost','Hsamp2_boost',...
    'Hsamp0_novacc','Hsamp1_novacc','Hsamp2_novacc',...
    'VEtot_samp','VEage_samp','IEtot_samp','IEage_samp','VEtot_boostsamp','VEage_boostsamp') %,'sampled_logLL_Vnofit')

%%
datesim(1:t0,:)=[];

%%
figure
subplot(2,3,1:3)
hold on
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi(tvacc:end,1); Hpi(end:-1:tvacc,2); Hpi(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 15],'--k')
plot(malweekPV,sum(rotamalPV,2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp0(tvacc:end,:)+Hsamp1(tvacc:end,:)+Hsamp2(tvacc:end,:),'r')
%plot(datesim(tvacc:end),Hpi(tvacc:end,:),'--r')
plot(datesim(tvacc:end),Hsamp0_novacc(tvacc:end,:)+Hsamp1_novacc(tvacc:end,:)+Hsamp2_novacc(tvacc:end,:),'b')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','mmm-yy')
xlim([datesim(tvacc-1) datesim(end)+7])
ylabel('Number of RVGE cases (per week)')
title('Homogeneity in vaccine response (no waning)')

subplot(2,3,4)
hold on
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi0(tvacc:end,1); Hpi0(end:-1:tvacc,2); Hpi0(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 10],'--k')
plot(malweekPV,sum(rotamalPV(:,1:12),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp0_novacc(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp0(tvacc:end,:),'r')
%plot(datesim(tvacc:end),Hpi0(tvacc:end,:),'--r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','mmm-yy')
xlim([datesim(tvacc-1) datesim(end)+7])
ylabel('Number of RVGE cases (per week)')
title('<1 yr olds')

subplot(2,3,5)
hold on
plot(malweekPV,sum(rotamalPV(:,13:au2),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp1(tvacc:end,1),'r')
plot(datesim(tvacc:end),Hsamp1_novacc(tvacc:end,1),'b')
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi1(tvacc:end,1); Hpi1(end:-1:tvacc,2); Hpi1(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 6],'--k')
plot(malweekPV,sum(rotamalPV(:,13:au2),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp1_novacc(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp1(tvacc:end,:),'r')
%plot(datesim(tvacc:end),Hpi1(tvacc:end,:),'--r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','mmm-yy')
xlim([datesim(tvacc-1) datesim(end)+7])
%ylabel('Number of RVGE cases (per week)')
title('1-<2 yr olds')
legend('Observed','Model predicted','Model predicted (no vaccination)','Orientation','Horizontal','Location','SO')

subplot(2,3,6)
hold on
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi2(tvacc:end,1); Hpi2(end:-1:tvacc,2); Hpi2(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 4],'--k')
plot(malweekPV,sum(rotamalPV(:,au2+1:au5),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp2_novacc(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp2(tvacc:end,:),'r')
%plot(datesim(tvacc:end),Hpi2(tvacc:end,:),'--r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','mmm-yy')
xlim([datesim(tvacc-1) datesim(end)+7])
%ylabel('Number of RVGE cases (per week)')
title('2-<5 yr olds')

%%
%samp_posterior_nofit=zeros(1000,1);
%for samp=1:1000
    %samp_posterior_bc(samp,1)=sampled_logLL_V_bc(500+10*samp-5000*floor((samp-.5)/500),ceil((samp-.5)/550))+log(betapdf(keep_pars_V_bc(10*samp-5000*floor((samp-.5)/500),1,ceil((samp-.5)/500)),24.7,11.3))+log(betapdf(keep_pars_V_bc(10*samp-5000*floor((samp-.5)/500),2,ceil((samp-.5)/500)),24.7,11.3));
%    samp_posterior_nofit(samp,1)=sampled_logLL_Vnofit(samp,1)+log(betapdf(betarnd(24.7,11.3),24.7,11.3))+log(betapdf(betarnd(24.7,11.3),24.7,11.3));
%end
%model_bic_nofit=log(295*6)*0-2*mean(samp_posterior_nofit)
