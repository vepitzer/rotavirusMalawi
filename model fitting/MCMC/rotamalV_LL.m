function [logLL,H]=rotamalV_LL(St,p,dataV,dataU,repeffV,time)
% This function is called from "rotamalawi_mcmcV...". It calculates the
% log-likelihood of the post-vaccination data ("dataV" and "dataU") given the
% model-predicted state variables ("St") and model parameters ("p"), while
% controlling for trends in reporting through time ("repeffV").

global rr1 rr2 ri2 ri3;  

b=p(1);
b1=p(2);
phi1=p(3);
al=p(4);
ar1=p(5);
ar2=p(6);
h=p(7);

%c1=100*ones(al); %Homogeneous mixing
c2=100*[ar1*ones(al,12) ar2*ones(al,12) ones(al,al-24)]; %Age-related acquisition for <2 yr olds
beta=(b/100)*c2; %/sum(N) type of mixing

immunity=[0.13; 0.03];
im=1;

hosp1=immunity(1,im)*h*ones(1,al); %proportion of first infections that are severe AND hospitalized
hosp2=immunity(2,im)*h*ones(1,al); %proportion of second infections that are severe AND hospitalized
hosp3=zeros(1,al); %proportion of subsequent infections that are severe AND hospitalized

%Initialize matrices to keep track of the outcomes of interest
lambda=zeros(size(St,1),al);
for t=1:size(St,1)
    lambda(t,:)=(1+b1*cos(2*pi*(time(t)-phi1)/52.18))*((St(t,2*al+1:3*al)+ri2*St(t,5*al+1:6*al)+ri3*St(t,8*al+1:9*al)+St(t,16*al+1:17*al)+ri2*St(t,19*al+1:20*al)+ri3*St(t,22*al+1:23*al))*beta)./sum(St(t,:)); 
end

Hu=zeros(size(St,1),al); Hv=zeros(size(St,1),al);
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
for i=1:size(obsV,1)
    for j=1:size(obsV,2)
        if Dv(i,j)>0
        llikl=llikl+obsV(i,j).*log(Dv(i,j)) - Dv(i,j) - sum(log(1:obsV(i,j))); %Log-likelihood of data assuming number of cases at time i, age j is Poisson-distributed
        end
        llikl=llikl+obsU(i,j).*log(Du(i,j)) - Du(i,j) - sum(log(1:obsU(i,j))); %Log-likelihood of data assuming number of cases at time i, age j is Poisson-distributed
    end
end

logLL=llikl; %Log-likelihood

