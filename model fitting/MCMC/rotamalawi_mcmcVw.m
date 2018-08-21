% NOTE: This code will take up to a week to run

load('malawidata.mat')

global B wm wi1 wi2 u um beta b1 phi1 d1 d2 rr1 rr2 ri2 ri3 al v1 v2 v3 sc1 sc2 sc3 sc2n sc3n wv reintro wA; % immig deaths b2 phi2 

%%% FIXED POPULATION PARAMETERS %%%
age=[0:1/12:23/12 2:4 5:5:75];
al=length(age); %Number of age groups
agep=[(1/12)*malpop(1)*ones(1,12) (1/48)*malpop(2)*ones(1,12) (1/4)*malpop(2)*ones(1,3) malpop(3:17).*ones(1,15)];
agep=agep/sum(agep); %Proportion of pop'n in each age group 
u=[1/4.3*ones(1,24) 1/52*ones(1,3) 1/(52*5)*ones(1,14) 1/(52*25)]; %Rate of aging out of each age group

t0=round(52.18*30); %burn-in period
tvacc=757; %start of post-vaccination data
tmax=t0+tvacc+length(rotamalVS3); %length of simulation
datepop=[1950 1 1 0 0 0; (1957:5:2037)' 7*ones(17,1) ones(17,1) zeros(17,3)]; %Date corresponding to census estimates
datesim=(datenum([1967 7 2 0 0 0]):7:datenum([2017 8 20 0 0 0]))'; %Date of simulation (weekly time step)

N1=870000; %Blantyre population in 2000
N=(N1/exp(.04*30))*agep; %Adjust initial population size for exponential growth occurring over the burn-in period
B=interp1(datenum(datepop),malawi_cbr/1000,datesim); %Interpolate weekly birth rate
if t0>0 && length(B)<tmax
    B=[B; B(end)*ones(tmax-length(B),1)]; %Assume birth rate is equal to last available data if simulating for longer
end
B=[B zeros(tmax,al-1)];
um=.010*ones(1,al);
um=log(1+um)/52.18;

%FIXED PARAMETERS

c1=100*ones(al); %homogeneous mixing
dur=1; %duration of infectiousness
d1=1/dur; %rate of recovery from primary infection (per week)
d2=2*d1; %rate of recovery from subsequent infection (per week)
rr1=0.62;  %relative risk of second infection 
rr2=0.35; %relative risk of third infection 
ri2=0.5;   %relative infectiousness of secondary infection
ri3=0.1; %relative infectiousness of asymptomatic infection
wm=1/4.3; %rate of waning maternal immunity (avg duration of immunity = 1mo) 
wi1=1/13; %rate of waning immunity following primary infection
wi2=1/13; %rate of waning immunity following 2nd infection
wA=0; %rate of waning immunity against symptomatic infection in adults

reintro=0;

%%% PREVIOUSLY ESTIMAED PARAMETERS %%%
fixed_pars=pfit_mal(:,1);

ptrans=fixed_pars(1); %probability of transmission given contact
b1=fixed_pars(2); %seasonal forcing 
phi1=fixed_pars(3); %seasonal offset

beta=(ptrans/100)*c1; %/sum(N) transmission matrix -- can change type of mixing

h=fixed_pars(4); %proportion of severe diarrhea cases hospitalized
hosp1=0.13*h*ones(1,al); %0.063*h*ones(1,al); %proportion of primary infections with severe diarrhea who are hospitalized
hosp2=0.03*h*ones(1,al); %0.077*h*ones(1,al); %proportion of secondary infections with severe diarrhea who are hospitalized
hosp3=zeros(1,al); %unkp(6)*h*ones(1,al); %proportion of subsequent infections with severe diarrhea who are hospitalized

pars=[ptrans; b1; phi1; al; 1; 1; h];

%%% VACCINATION PARAMETERS %%%
v1=zeros(tmax,al); v2=zeros(tmax,al); v3=zeros(tmax,al); %initialize vaccination rate across all ages
sc3=0; sc3n=0;

vcov=vcov_mavg;
avacc=[3 4 0]; 
v1(:,avacc(1))=[zeros(t0+tvacc,1); vcov(:,1)]; 
v2(:,avacc(2))=[zeros(t0+tvacc,1); vcov(:,2)]; 


%%% DEFINE CHARACTERISTICS OF THE MCMC %%%
npar=3; %Number of parameters to be estimated
nchains=2;
burnin_steps=500;
run_steps=5000;
total_steps=burnin_steps+run_steps;
wr=0.5;
 
% Initialize a vector to store the sampled parameters
sampled_pars_V=zeros(total_steps,npar,nchains); 
sampled_logLL_V=zeros(total_steps,nchains); 
accept=zeros(total_steps,nchains);
accept_rate_V=zeros(total_steps,nchains); 
keep_pars_V=zeros(total_steps-burnin_steps,npar,nchains);

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


%%% ESTIMATED PARAMETERS %%%

for chain=1:nchains
if chain==1
old_sc1=.6; 
old_sc2=.8;
old_wv=1;
else
old_sc1=.7; 
old_sc2=.7;
old_wv=3;
end

%%% INPUT PARAMETERS %%%
sc1=old_sc1;
sc2=old_sc2;
sc2n=sc2;
wv=1/(52*old_wv);


options=odeset('NonNegative',1:length(St0)); %force solutions to differential equations to be non-negative
[time,St]=ode45('rasisV1w',1:tmax,St0,options);

time(1:t0+tvacc,:)=[]; %delete output from from burn-in period
St(1:t0+tvacc,:)=[]; 

old_logLL=rotamalV_LL(St,pars,rotamalVS3_vacc,rotamalVS3_unvacc,repeffV3,time);
init_logLL=old_logLL;


for i=1:total_steps
    
    log_proposal_ratio=0;
    log_prior_ratio=0;
    
    % Propose new value for sc1
    current_sc1=old_sc1+.05*randn; %*exp(wr*randn);
    log_prior_ratio=log_prior_ratio+log(betapdf(current_sc1,24.7,11.3))-log(betapdf(old_sc1,24.7,11.3)); 
    sc1=current_sc1;
    
    % Propose new value for sc2
    current_sc2=old_sc2+.05*randn; %*exp(wr*randn);
    log_prior_ratio=log_prior_ratio+log(betapdf(current_sc2,24.7,11.3))-log(betapdf(old_sc2,24.7,11.3)); 
    sc2=current_sc2; 
    sc2n=sc2;
       
    % Propose new value for wv
    current_wv=old_wv*exp(wr*randn);
    log_proposal_ratio=log_proposal_ratio+log(current_wv)-log(old_wv);
    log_prior_ratio=log_prior_ratio+log(unifpdf(current_wv,0,100))-log(unifpdf(old_wv,0,100)); 
    wv=1/(52*current_wv);
 

    options=odeset('NonNegative',1:length(St0)); %force solutions to differential equations to be non-negative
    [time,St]=ode45('rasisV1w',1:tmax,St0,options);

    time(1:t0+tvacc,:)=[]; %delete output from from burn-in period
    St(1:t0+tvacc,:)=[]; 
    
    current_logLL=rotamalV_LL(St,pars,rotamalVS3_vacc,rotamalVS3_unvacc,repeffV3,time);
    
    log_likelihood_ratio=current_logLL-old_logLL;
    acceptance_prob=exp(log_proposal_ratio+log_likelihood_ratio+log_prior_ratio);
    if rand < acceptance_prob
        old_sc1=current_sc1;
        old_sc2=current_sc2;
        old_wv=current_wv;
        old_logLL=current_logLL;
        accept(i,chain)=1;
    else
        sc1 = old_sc1;
        sc2 = old_sc2;
        wv = 1/(52*old_wv);
    end
    sampled_pars_V(i,:,chain)=[old_sc1 old_sc2 old_wv];  
    sampled_logLL_V(i,chain)=old_logLL;
end
end

% Calculate acceptance rate
for i=2:total_steps
    accept_rate_V(i,:)=mean(accept(1:i,:));
end

% Only keep values after the burn-in period
for chain=1:2
    keep_pars_V(:,:,chain)=sampled_pars_V(burnin_steps+1:total_steps,:,chain);
end

%%
save('malawi_mcmc_outputV.mat','sampled_pars_V','sampled_logLL_V','accept_rate_V','keep_pars_V')

%%
figure
subplot(2,3,1); plot([sampled_pars_V(:,1,1) sampled_pars_V(:,1,2)]); hold on; plot([burnin_steps burnin_steps],[0 1],'--k'); xlabel('Iteration'); ylabel('\its_1'); xlim([0 6000])
subplot(2,3,2); plot([sampled_pars_V(:,2,1) sampled_pars_V(:,2,2)]); hold on; plot([burnin_steps burnin_steps],[0 1],'--k'); xlabel('Iteration'); ylabel('\its_2'); xlim([0 6000])
subplot(2,3,3); plot([sampled_pars_V(:,3,1) sampled_pars_V(:,3,2)]); hold on; plot([burnin_steps burnin_steps],[0 2],'--k'); xlabel('Iteration'); ylabel({'Duration of immunity';'(years)'}); xlim([0 6000])
subplot(2,3,4); hist([keep_pars_V(:,1,1) keep_pars_V(:,1,2)]); xlabel('\its_1'); ylabel('No. of posterior samples')
subplot(2,3,5); hist([keep_pars_V(:,2,1) keep_pars_V(:,2,2)]); xlabel('\its_2'); ylabel('No. of posterior samples')
subplot(2,3,6); hist([keep_pars_V(:,3,1) keep_pars_V(:,3,2)]); xlabel('Duration of immunity (years)'); ylabel('No. of posterior samples')

%%
figure
subplot(2,1,1); plot(sampled_logLL_V)
subplot(2,1,2); plot(accept_rate_V)

%% Gelman-Rubin diagnostic(?)
x=var(keep_pars_V(:,:,1));
x(2,:)=var(keep_pars_V(:,:,2))
y=mean(keep_pars_V(:,:,1));
y(2,:)=mean(keep_pars_V(:,:,2))

grdiag_V=sqrt(((run_steps-1)/run_steps*mean(x) + var(y))./mean(x))

