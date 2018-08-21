load('malawidata.mat')

clear global all
global B wm wi1 wi2 u um d1 d2 rr1 rr2 ri2 ri3 v1 v2 v3 reintro wA;  

%%% FIXED POPULATION PARAMETERS %%%
age=[0:1/12:23/12 2:4 5:5:75];
al=length(age); %number of age groups
agep=[(1/12)*malpop(1)*ones(1,12) (1/48)*malpop(2)*ones(1,12) (1/4)*malpop(2)*ones(1,3) malpop(3:17).*ones(1,15)];
agep=agep/sum(agep); %Proportion of pop'n in each age group (square age distribution)
u=[1/4.3*ones(1,24) 1/52*ones(1,3) 1/(52*5)*ones(1,14) 1/(52*25)]; %rate of aging out of each age group

t0=round(52.18*30); %burn-in period
tmax=t0+length(rotamalT2)+length(rotamal08T2)+26; %length of simulation
datepop=[1950 1 1 0 0 0; (1957:5:2037)' 7*ones(17,1) ones(17,1) zeros(17,3)]; %date corresponding to census estimates
datesim=(datenum([1967 7 2 0 0 0]):7:datenum([2009 12 27 0 0 0]))'; %date of simulation (weekly time step)

N1=870000; %Blantyre population in 2000
N=(N1/exp(.04*30))*agep; %Adjust initial population size for exponential growth occurring over the burn-in period
B=interp1(datenum(datepop),malawi_cbr/1000,datesim); %Interpolate weekly birth rate
if t0>0 && length(B)<tmax
    B=[B; B(end)*ones(tmax-length(B),1)]; %Assume birth rate is equal to last available data if simulating for longer
end
B=[B zeros(tmax,al-1)]; %Birth rate for each age group
um=.010*ones(1,al); %Annual crude death rate (to adjust population size to approximate Blantyre population growth) 
um=log(1+um)/52.18; %Weekly crude death rate


%%% FIXED INFECTION PARAMETERS %%%
c1=100*ones(al); %homogeneous mixing
dur=1; %duration of infectiousness
d1=1/dur; %rate of recovery from primary infection (per week)
d2=2*d1; %rate of recovery from subsequent infection (per week)
rr1=0.62;  %relative risk of second infection 
rr2=0.35; %relative risk of third infection 
ri2=0.5;   %relative infectiousness of secondary infection
ri3=0.3;   %relative infectiousness of asymptomatic infection 
wm=1/4.3; %rate of waning maternal immunity (distribution gamma(6,1))
wi1=1/13; %1; %rate of waning immunity following primary infection
wi2=1/13; %1; %rate of waning immunity following 2nd infection
wA=0; %rate of waning of immunity to symptomatic infections (S2->S1)
reintro=0;

%%% VACCINATION PARAMETERS %%%
v1=zeros(tmax,al); v2=zeros(tmax,al); v3=zeros(tmax,al); %initialize vaccination rate across all ages


%%% DEFINE CHARACTERISTICS OF THE MCMC %%%
npar=4; %Number of parameters to be estimated
nchains=2;
burnin_steps=1000;
run_steps=10000;
total_steps=burnin_steps+run_steps;
wr=0.1;
 
% Initialize a vector to store the sampled parameters
sampled_pars_3=zeros(total_steps,npar,nchains);
sampled_logLL_3=zeros(total_steps,nchains);
accept=zeros(total_steps,nchains);
accept_rate_3=zeros(total_steps,nchains);
keep_pars_3=zeros(total_steps-burnin_steps,npar,nchains);

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
St0(10*al+1:11*al,1)=zeros(1,al); %Maternal immunity
St0(11*al+1:12*al,1)=zeros(1,al); %
St0(12*al+1:13*al,1)=zeros(1,al); %
St0(13*al+1:14*al,1)=zeros(1,al); %
St0(14*al+1:15*al,1)=zeros(1,al); %

%%% ESTIMATED PARAMETERS %%%

for chain=1:nchains
if chain==1
old_ptrans=30; 
old_b1=.1;
old_phi=9;
old_h=.02;
else
old_ptrans=40; 
old_b1=.2;
old_phi=7;
old_h=.015;
end

%%% INPUT PARAMETERS %%%
ptrans=old_ptrans; %~R0
beta=(ptrans/100)*c1; %/sum(N) %type of mixing

b1=old_b1; %seasonal forcing 
phi1=old_phi; %seasonal offset
h=old_h; %proportion of severe diarrhea cases hospitalized

pars=[ptrans; b1; phi1; al; 1; 1; h];

%clear St lambda H %clear outcome variables which may be in memory

options=odeset('NonNegative',1:length(St0)); %force solutions to differential equations to be non-negative
[time St]=ode45(@(t,St) rasisM1(t,St,pars),1:tmax,St0,options); %solve differential equations defined in 'rasisM1'

time(1:t0,:)=[]; %delete output from from burn-in period
St(1:t0,:)=[]; 

old_logLL=rotamal_LL(St,pars,rotamalT2,rotamal08T2,negtest_movavg,negtest08_movavg,time);
init_logLL=old_logLL;


for i=1:total_steps
    
    log_proposal_ratio=0;
    log_prior_ratio=0;
    
    % Propose new value for ptrans
    current_ptrans=old_ptrans*exp(wr*randn);
    log_proposal_ratio=log_proposal_ratio+log(current_ptrans)-log(old_ptrans);
    log_prior_ratio=log_prior_ratio+log(unifpdf(current_ptrans,1,100))-log(unifpdf(old_ptrans,1,100)); 
    ptrans=current_ptrans;
    beta=(ptrans/100)*c1; %/sum(N)type of mixing

    % Propose new value for b1
    current_b1=old_b1*exp(wr*randn);
    log_proposal_ratio=log_proposal_ratio+log(current_b1)-log(old_b1);
    log_prior_ratio=log_prior_ratio+log(unifpdf(current_b1,0,1))-log(unifpdf(old_b1,0,1)); %Uniform(0,1) prior
    b1=current_b1;
       
    % Propose new value for phi
    current_phi=old_phi*exp(wr*randn);
    log_proposal_ratio=log_proposal_ratio+log(current_phi)-log(old_phi);
    log_prior_ratio=log_prior_ratio+log(unifpdf(current_phi,0,52.18))-log(unifpdf(old_phi,0,52.18)); 
    phi1=current_phi;
 
    % Propose new value for h
    current_h=old_h*exp(wr*randn);
    log_proposal_ratio=log_proposal_ratio+log(current_h)-log(old_h);
    log_prior_ratio=log_prior_ratio+log(unifpdf(current_h,0,1))-log(unifpdf(old_h,0,1)); %Uniform(0,1) prior (i.e. all values between 0 and 1 are equally likely)
    h=current_h;
    
    pars=[ptrans; b1; phi1; al; 1; 1; h];

    options=odeset('NonNegative',1:length(St0)); %force solutions to differential equations to be non-negative
    [time St]=ode45(@(t,St) rasisM1(t,St,pars),1:tmax,St0,options); %solve differential equations defined in 'rasisM1'

    time(1:t0,:)=[]; %delete output from from burn-in period
    St(1:t0,:)=[]; 

    current_logLL=rotamal_LL(St,pars,rotamalT2,rotamal08T2,negtest_movavg,negtest08_movavg,time);
    log_likelihood_ratio=current_logLL-old_logLL;
    acceptance_prob=exp(log_proposal_ratio+log_likelihood_ratio+log_prior_ratio);
    if rand < acceptance_prob
        old_ptrans=current_ptrans;
        old_b1=current_b1;
        old_phi=current_phi;
        old_h=current_h;
        old_logLL=current_logLL;
        accept(i,chain)=1;
    else
        ptrans = old_ptrans;
        b1 = old_b1;
        phi1 = old_phi;
        h = old_h;
    end
    sampled_pars_3(i,:,chain)=[old_ptrans old_b1 old_phi old_h];  
    sampled_logLL_3(i,chain)=old_logLL;
end
end

% Calculate acceptance rate
for i=2:total_steps
    accept_rate_3(i,:)=mean(accept(1:i,:));
end

% Only keep values after the burn-in period
for chain=1:2
    keep_pars_3(:,:,chain)=sampled_pars_3(burnin_steps+1:total_steps,:,chain);
end

save('malawi_mcmc_output3.mat','sampled_pars_3','sampled_logLL_3','accept_rate_3','keep_pars_3')

%%
figure
subplot(2,4,1); plot([sampled_pars_3(:,1,1) sampled_pars_3(:,1,2)])
subplot(2,4,2); plot([sampled_pars_3(:,2,1) sampled_pars_3(:,2,2)])
subplot(2,4,3); plot([sampled_pars_3(:,3,1) sampled_pars_3(:,3,2)])
subplot(2,4,4); plot([sampled_pars_3(:,4,1) sampled_pars_3(:,4,2)])
subplot(2,4,5); hist([keep_pars_3(:,1,1) keep_pars_3(:,1,2)])
subplot(2,4,6); hist([keep_pars_3(:,2,1) keep_pars_3(:,2,2)])
subplot(2,4,7); hist([keep_pars_3(:,3,1) keep_pars_3(:,3,2)])
subplot(2,4,8); hist([keep_pars_3(:,4,1) keep_pars_3(:,4,2)])

%%
figure
subplot(2,1,1); plot(sampled_logLL_3)
subplot(2,1,2); plot(accept_rate_3)

%% Gelman-Rubin diagnostic
run_steps=10000;
x=var(keep_pars_3(:,:,1));
x(2,:)=var(keep_pars_3(:,:,2))
y=mean(keep_pars_3(:,:,1));
y(2,:)=mean(keep_pars_3(:,:,2))

grdiag_3=sqrt(((run_steps-1)/run_steps*mean(x) + var(y))./mean(x))
