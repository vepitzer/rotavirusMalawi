function [LL,D,D2]=rotafitM1(unkp,data,data2,malpop,Bmal,varargin)
% Calculates the negative log-likelihood of the pre-vaccination data
% ("data" and "data2", representing cases from 1997-2007 and 2008-2009, 
% respectively) given the model parameters ("unkp", estimated),
% Malawi population age distribution ("malpop"), birth rate through time 
% ("Bmal"), and other possible inputs ("varargin", see below).

% Call using "fminsearch" to obtain the MAP estimate (e.g. see "malawimodelfit.m") 
 
global B wm wi1 wi2 u um d1 d2 rr1 rr2 ri2 ri3 v1 v2 v3 reintro wA; 

if nargin>=5
    negtest=varargin{1};
    negtest2=varargin{2};
else
    negtest=ones(length(data),1); %reporting effort
    negtest2=ones(length(data2),1); %reporting effort
end
if nargin>=7
    ri3=varargin{3};
else ri3=0.1; %default relative infectiousness of asymptomatic infections
end

if sum(unkp<0)>0 %Constrain all parameters to be positive
    LL=Inf;    
else
 
%%% FIXED POPULATION PARAMETERS %%%
age=[0:1/12:23/12 2:4 5:5:75];
al=length(age); %number of age groups
agep=[(1/12)*malpop(1)*ones(1,12) (1/48)*malpop(2)*ones(1,12) (1/4)*malpop(2)*ones(1,3) malpop(3:17).*ones(1,15)];
agep=agep/sum(agep); %Proportion of pop'n in each age group (square age distribution)
u=[1/4.3*ones(1,24) 1/52*ones(1,3) 1/(52*5)*ones(1,14) 1/(52*25)]; %rate of aging out of each age group

t0=round(52.18*30); %burn-in period
tmax=t0+length(data)+length(data2)+26; %length of simulation
datepop=[1950 1 1 0 0 0; (1957:5:2037)' 7*ones(17,1) ones(17,1) zeros(17,3)]; %date corresponding to census estimates
datesim=(datenum([1967 7 2 0 0 0]):7:datenum([2009 12 27 0 0 0]))'; %date of simulation (weekly time step)

N1=870000; %Blantyre population in 2000
N=(N1/exp(.04*30))*agep; %Adjust initial population size for exponential growth occurring over the burn-in period
B=interp1(datenum(datepop),Bmal/1000,datesim); %Interpolate weekly birth rate from census data
B=[B zeros(tmax,al-1)]; %Birth rate for each age group
um=.010*ones(1,al); %Annual crude death rate (to adjust population size to approximate Blantyre population growth) 
um=log(1+um)/52.18; %Weekly crude death rate


%%% INFECTION PARAMETERS %%%
dur=1; %duration of infectiousness
d1=1/dur; %rate of recovery from primary infection (per week)
d2=2*d1; %rate of recovery from subsequent infection (per week)
rr1=0.62;  %relative risk of second infection 
rr2=0.35; %relative risk of third infection 
ri2=0.5;   %relative infectiousness of secondary infection
%ri3=0.1;   %relative infectiousness of asymptomatic infection 
wm=1/4.3; %rate of waning maternal immunity (distribution gamma(6,1))
wi1=1/13; %1; %rate of waning immunity following primary infection
wi2=1/13; %1; %rate of waning immunity following 2nd infection
wA=0; %rate of waning of immunity to symptomatic infections (S2->S1)

%%% INPUT PARAMETERS %%%
ptrans=unkp(1); %probability of transmission given contact
b1=unkp(2); %0; %seasonal forcing 
phi1=unkp(3); %0; %seasonal offset
ar1=1; %unkp(4); %RR for <1 yr olds
ar2=1; %unkp(5); %RR for <2 yr olds

%c1=100*ones(al); %Homogeneous mixing
c2=100*[ar1*ones(al,12) ar2*ones(al,12) ones(al,al-24)]; %Age-related acquisition for <2 yr olds

beta=(ptrans/100)*c2; %/sum(N) transmission matrix -- can change type of mixing

h=unkp(4); %unkp(2); %proportion of severe diarrhea cases hospitalized
hosp1=0.13*h*ones(1,al); %proportion of primary infections with severe diarrhea who are hospitalized
hosp2=0.03*h*ones(1,al); %proportion of secondary infections with severe diarrhea who are hospitalized
hosp3=zeros(1,al); %proportion of subsequent infections with severe diarrhea who are hospitalized

pars=[ptrans; b1; phi1; al; ar1; ar2];

reintro=0;

%%% VACCINATION PARAMETERS %%%
v1=zeros(tmax,al); v2=zeros(tmax,al); v3=zeros(tmax,al); %initialize vaccination rate across all ages

%Initialize vector to keep track of the number of people in each state
St0=[];
St0(1:al,1)=[N(1) zeros(1,al-1)]; %Maternal immunity
St0(al+1:2*al,1)=[0 N(2:al)-ones(1,al-1)]; %Susceptible_0
St0(2*al+1:3*al,1)=[0 ones(1,al-1)]; %Infectious_1 (primary) 
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

clear St lambda H %clear outcome variables which may be in memory

options=odeset('NonNegative',1:length(St0)); %force solutions to differential equations to be non-negative
[time St]=ode45(@rasisM1,1:tmax,St0,options,pars); %solve differential equations defined in 'rasisM1'

time(1:t0,:)=[]; %delete output from from burn-in period
St(1:t0,:)=[]; 


lambda=zeros(tmax-t0,al);
for t=1:size(St,1)
    lambda(t,:)=(1+b1*cos(2*pi*(time(t)-phi1)/52.18))*((St(t,2*al+1:3*al)+ri2*St(t,5*al+1:6*al)+ri3*St(t,8*al+1:9*al))*beta)./sum(St(t,:));
end

H=zeros(tmax-t0,al);
for i=1:al %calculate number of rotavirus hospitalizations in each age group across time
    H(:,i)=max(0,hosp1(i)*St(:,al+i).*lambda(:,i)+hosp2(i)*rr1*St(:,4*al+i).*lambda(:,i)+hosp3(i)*rr2*St(:,7*al+i).*lambda(:,i));
end

obs=data; %observed # of hospitalizations in children <5 yrs old (1997-2007)
obs2=data2; %observed # of hospitalizations in children <5 yrs old (2008-2009)
repeff=(negtest/mean([negtest; negtest2]))*ones(1,size(obs,2)); %relative reporting effort for 1997-2007
repeff2=(negtest2/mean([negtest; negtest2]))*ones(1,size(obs,2)); %relative reporting effort for 2008-2009
D=repeff.*H(1:size(obs,1),1:size(obs,2)); %predicted # of hospitalizations in children <5 yrs old (1997-2007)
D2=repeff2.*H(size(obs,1)+27:size(obs,1)+size(obs2,1)+26,1:size(obs2,2)); %predicted # of hospitalizations in children <5 yrs old (2008-2009)

%%% CALCULATE MODEL FIT (Negative log-likelihood) %%%

%nbr=unkp(5); %shape parameter for negative binomial distribution
llikl=0;
for i=1:size(obs,1)
    for j=1:size(obs,2)
        llikl=llikl+obs(i,j).*log(D(i,j)) - D(i,j) - sum(log(1:obs(i,j))); %Log-likelihood of data assuming number of cases at time i, age j is Poisson-distributed
        %llikl=llikl+log(nbinpdf(obs(i,j),nbr,nbr/(nbr+D(i,j))));    
    end
end
for i=1:size(obs2,1)
    for j=1:size(obs2,2)
        llikl=llikl+obs2(i,j).*log(D2(i,j)) - D2(i,j) - sum(log(1:obs2(i,j))); %Log-likelihood of data assuming number of cases at time i, age j is Poisson-distributed
        %llikl=llikl+log(nbinpdf(obs2(i,j),nbr,nbr/(nbr+D2(i,j))));    
    end
end

LL=-llikl; %Negative log-likelihood
    
end