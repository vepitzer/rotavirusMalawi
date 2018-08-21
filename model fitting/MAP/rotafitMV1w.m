function [LL,Dv]=rotafitMV1w(unkp,pars,dataV,dataU,malpop,Bmal,repeffV,vcov,varargin)
% This is for the model assuming waning of vaccine-induced immunity and
% homogeneity in vaccine response. 

% Calculates the negative log-likelihood of the post-vaccination data
% ("dataV" and "dataU", representing cases among vaccinated and
% unvaccinated individuals, respectively) given the model parameters
% ("unkp", estimated, and "pars", estimated from the pre-vaccination data),
% Malawi population age distribution ("malpop"), birth rate through time 
% ("Bmal"), variation in reporting effort through time ("repeffV"), vaccine
% coverage ("vcov"), and other possible inputs ("varargin", see below).

% Call using "fminsearch" to obtain the MAP estimate. 

global B wm wi1 wi2 u um beta b1 phi1 d1 d2 rr1 rr2 ri2 ri3 al v1 v2 v3 sc1 sc2 sc3 sc2n sc3n wv reintro wA; 

priorLL=log(betapdf(unkp(1),24.7,11.3))+log(betapdf(unkp(2),24.7,11.3)); %Prior distribution for the probability of responding to each vaccine dose

if priorLL<-1000000 %Do not bother simulating model if prior is unreasonable
    LL=Inf;    
else

if isempty(varargin)
    ri3=0.1; %Default assumption for relative infectiousness of subsequent infections
else
    ri3=varargin{1};
end
    
%%% FIXED POPULATION PARAMETERS %%%
age=[0:1/12:23/12 2:4 5:5:75];
al=length(age); %Number of age groups
agep=[(1/12)*malpop(1)*ones(1,12) (1/48)*malpop(2)*ones(1,12) (1/4)*malpop(2)*ones(1,3) malpop(3:17).*ones(1,15)];
agep=agep/sum(agep); %Proportion of pop'n in each age group 
u=[1/4.3*ones(1,24) 1/52*ones(1,3) 1/(52*5)*ones(1,14) 1/(52*25)]; %Rate of aging out of each age group

t0=round(52.18*30); %burn-in period
tvacc=757; %start of post-vaccination data
tmax=t0+tvacc+length(dataV); %length of simulation
datepop=[1950 1 1 0 0 0; (1957:5:2037)' 7*ones(17,1) ones(17,1) zeros(17,3)]; %Date corresponding to census estimates
datesim=(datenum([1967 7 2 0 0 0]):7:datenum([2017 8 20 0 0 0]))'; %Date of simulation (weekly time step)

N1=870000; %Blantyre population in 2000
N=(N1/exp(.04*30))*agep; %Adjust initial population size for exponential growth occurring over the burn-in period
B=interp1(datenum(datepop),Bmal/1000,datesim);
B=[B zeros(tmax,al-1)];
um=.010*ones(1,al);
um=log(1+um)/52.18;


%%% INFECTION PARAMETERS %%%
dur=1; %duration of infectiousness
d1=1/dur; %rate of recovery from primary infection (per week)
d2=2*d1; %rate of recovery from subsequent infection (per week)
rr1=0.62;  %relative risk of second infection 
rr2=0.35; %relative risk of third infection 
ri2=0.5;   %relative infectiousness of secondary infection
%ri3=0.1;   %relative infectiousness of asymptomatic infection 
wm=1/4.3; %rate of waning maternal immunity (distribution gamma(6,1))
wi1=1/13; %unkp(7); %rate of waning immunity following primary infection
wi2=1/13; %unkp(7); %rate of waning immunity following 2nd infection
wA=0; %rate of waning of immunity to symptomatic infections (S2->S1)

%%% INPUT PARAMETERS %%%
ptrans=pars(1); %probability of transmission given contact
b1=pars(2); %seasonal forcing 
phi1=pars(3); %seasonal offset
ar1=1; %unkp(4); %RR for <1 yr olds
ar2=1; %unkp(5); %RR for <2 yr olds

%c1=100*ones(al); %Homogeneous mixing
c2=100*[ar1*ones(al,12) ar2*ones(al,12) ones(al,al-24)]; %Age-related acquisition for <2 yr olds

beta=(ptrans/100)*c2; %/sum(N) transmission matrix -- can change type of mixing

h=pars(4); %proportion of severe diarrhea cases hospitalized
hosp1=0.13*h*ones(1,al); %proportion of primary infections with severe diarrhea who are hospitalized
hosp2=0.03*h*ones(1,al); %proportion of secondary infections with severe diarrhea who are hospitalized
hosp3=zeros(1,al); %proportion of subsequent infections with severe diarrhea who are hospitalized

reintro=0; %background risk of infection

%%% VACCINATION PARAMETERS %%%
v1=zeros(tmax,al); v2=zeros(tmax,al); v3=zeros(tmax,al); %initialize vaccination rate across all ages

avacc=[3 4 0]; 
sc1=unkp(1);  
sc2=unkp(2);  
sc3=unkp(2);  
sc2n=sc2; %probability of responding to second dose given did not respond to first dose
sc3n=sc3; %probability of responding to third dose given did not respond to first or second dose

v1(:,avacc(1))=[zeros(t0+tvacc,1); vcov(:,1)]; 
v2(:,avacc(2))=[zeros(t0+tvacc,1); vcov(:,2)]; 

wv=unkp(3); %rate of waning of vaccine-induced immunity

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
St0(10*al+1:11*al,1)=zeros(1,al); %Maternal immunity 2
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

clear St lambda H %clear outcome variables which may be in memory

options=odeset('NonNegative',1:length(St0)); %force solutions to differential equations to be non-negative
[time St]=ode45('rasisV1w',1:tmax,St0,options);

time(1:t0+tvacc,:)=[]; %delete output from from burn-in period
St(1:t0+tvacc,:)=[]; 


lambda=zeros(tmax-t0-tvacc,al);
for t=1:size(St,1)
    lambda(t,:)=(1+b1*cos(2*pi*(time(t)-phi1)/52.18))*((St(t,2*al+1:3*al)+ri2*St(t,5*al+1:6*al)+ri3*St(t,8*al+1:9*al)+St(t,16*al+1:17*al)+ri2*St(t,19*al+1:20*al)+ri3*St(t,22*al+1:23*al))*beta)./sum(St(t,:)); 
end

Hu=zeros(tmax-t0-tvacc,al); Hv=zeros(tmax-t0-tvacc,al);
for i=1:al %calculate number of rotavirus hospitalizations in each age group across time
    Hu(:,i)=max(0,hosp1(i)*St(:,al+i).*lambda(:,i)+hosp2(i)*rr1*St(:,4*al+i).*lambda(:,i)+hosp3(i)*rr2*St(:,7*al+i).*lambda(:,i));
    Hv(:,i)=max(0,hosp1(i)*St(:,15*al+i).*lambda(:,i)+hosp2(i)*rr1*St(:,18*al+i).*lambda(:,i)+hosp3(i)*rr2*St(:,21*al+i).*lambda(:,i));
end
H=Hu+Hv;

obsV=dataV; %observed # of hospitalizations among vaccinated children
obsU=dataU; %observed # of hospitalizations among unvaccinated children
Dv=repeffV(:,[1 13 25]).*[sum(Hv(:,1:12),2) sum(Hv(:,13:24),2) sum(Hv(:,25:end),2)]; %predicted # of cases by year of age among vaccinated children
Du=repeffV(:,[1 13 25]).*[sum(Hu(:,1:12),2) sum(Hu(:,13:24),2) sum(Hu(:,25:end),2)]; %predicted # of cases by year of age among unvaccinated children

%%% CALCULATE MODEL FIT (Negative log-likelihood) %%%

llikl=0;
%fit to post-vaccination data
for i=1:size(obsV,1)
    for j=1:size(obsV,2)
        if Dv(i,j)>0
        llikl=llikl+obsV(i,j).*log(Dv(i,j)) - Dv(i,j) - sum(log(1:obsV(i,j))); %Log-likelihood of data assuming number of cases at time i, age j is Poisson-distributed
        end
        llikl=llikl+obsU(i,j).*log(Du(i,j)) - Du(i,j) - sum(log(1:obsU(i,j))); %Log-likelihood of data assuming number of cases at time i, age j is Poisson-distributed
    end
end

LL=-llikl-priorLL; %Negative log-likelihood
    
end