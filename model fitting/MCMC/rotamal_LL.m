function [logLL,H]=rotamal_LL(St,p,data,data2,negtest,negtest2,time)
% This function is called from "rotamalawi_mcmc1". It calculates the
% log-likelihood of the prevaccination data ("data" and "data2") given the
% model-predicted state variables ("St") and model parameters ("p"), while
% controlling for trends in reporting through time (as indicated by the
% number of rotavirus-negative cases "negtest" and "negtest2").

global rr1 rr2 ri2 ri3; %al beta b1 phi1 

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

immunity=[0.13 0.063; 0.03 0.077];
im=1;

hosp1=immunity(1,im)*h*ones(1,al); %proportion of first infections that are severe AND hospitalized
hosp2=immunity(2,im)*h*ones(1,al); %proportion of second infections that are severe AND hospitalized
hosp3=zeros(1,al); %proportion of subsequent infections that are severe AND hospitalized

%Initialize matrices to keep track of the outcomes of interest
lambda=zeros(length(St),al); %Force of infection
H=zeros(length(St),al); %Number of hospitalizations

for t=1:length(St)
    lambda(t,:)=(1+b1*cos(2*pi*(time(t)-phi1)/52.18))*((St(t,2*al+1:3*al)+ri2*St(t,5*al+1:6*al)+ri3*St(t,8*al+1:9*al))*beta)./sum(St(t,:));
end

for i=1:al %calculate number of rotavirus hospitalizations in each age group across time
    H(:,i)=max(0,hosp1(i)*St(:,al+i).*lambda(:,i)+hosp2(i)*rr1*St(:,4*al+i).*lambda(:,i)+hosp3(i)*rr2*St(:,7*al+i).*lambda(:,i));
end

%observed # of hospitalizations in children <5 yrs old
obs=data; 
obs2=data2;
%Adjust for trend in 105-week moving average of the number of rotavirus-negative cases
repeff=(negtest/mean([negtest; negtest2]))*ones(1,size(obs,2)); 
repeff2=(negtest2/mean([negtest; negtest2]))*ones(1,size(obs,2));
%model-predicted number of hospitalizations in children <5 yrs old
D=repeff.*H(1:size(obs,1),1:size(obs,2));
D2=repeff2.*H(size(obs,1)+27:size(obs,1)+size(obs2,1)+26,1:size(obs2,2));

%%% CALCULATE MODEL FIT %%%

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

logLL=llikl; %Log-likelihood

