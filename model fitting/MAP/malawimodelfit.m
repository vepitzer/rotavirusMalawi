load('malawidata.mat')

%% Including 2008-09 data, Homogeneous mixing, 1-month age categories up to 2 yrs of age, using negative tests to control for reporting trend
for i=1:5
    pfit_mal(:,i)=fminsearch(@(p) rotafitM1(p,rotamalT2,rotamal08T2,malpop,malawi_cbr,negtest_movavg,negtest08_movavg,i/10),[50; .1; 10; .02])
end

%%
for i=1:5
LLmalfit(i)=rotafitM1(pfit_mal(:,i),rotamalT2,rotamal08T2,malpop,malawi_cbr,negtest_movavg,negtest08_movavg,i/10)
end

%% Same as above, but using negative binomial distribution instead of Poisson
for i=1:5
    pfitNB_mal(:,i)=fminsearch(@(p) rotafitM1nb(p,rotamalT2,rotamal08T2,malpop,malawi_cbr,negtest_movavg,negtest08_movavg,i/10),[pfit_mal(:,i); 10])
end

%%
for i=1:5
LLmalfitNB(i)=rotafitM1nb(pfitNB_mal(:,i),rotamalT2,rotamal08T2,malpop,malawi_cbr,negtest_movavg,negtest08_movavg,i/10)
end


%% Log-likelihood of base-case model predictions validated against post-vaccination data
[LLvacpred_bc,Hvacpred_bc]=rotafitMV1([.687; .687],pfit_mal(:,1),rotamalVS3_vacc,rotamalVS3_unvacc,malpop,malawi_cbr,repeffV3,vcov_mavg);
[LLvacpred_nr,Hvacpred_nr]=rotafitMV1nr([.528; .895],pfit_mal(:,1),rotamalVS3_vacc,rotamalVS3_unvacc,malpop,malawi_cbr,repeffV3,vcov_mavg);

%% Fit model with waning vaccine-induced immunity to the post-vaccination data
pvacc1w_mal=fminsearch(@(p) rotafitMV1w(p,pfit_mal(:,1),rotamalVS3_vacc,rotamalVS3_unvacc,malpop,malawi_cbr,repeffV3,vcov_mavg),[.687; .687; 1/52])
pvacc1wnr_mal=fminsearch(@(p) rotafitMV1wnr(p,pfit_mal(:,1),rotamalVS3_vacc,rotamalVS3_unvacc,malpop,malawi_cbr,repeffV3,vcov_mavg),[.687; .687; 1/52])

[LLvacfit_w,Hvacfit_w]=rotafitMV1w(pvacc1w_mal,pfit_mal(:,1),rotamalVS3_vacc,rotamalVS3_unvacc,malpop,malawi_cbr,repeffV3,vcov_mavg);
[LLvacfit_wnr,Hvacfit_wnr]=rotafitMV1wnr(pvacc1wnr_mal,pfit_mal(:,1),rotamalVS3_vacc,rotamalVS3_unvacc,malpop08,malawi_cbr,repeffV3,vcov_mavg);

%% Fit base-case model to the post-vaccination data
pvacc1_mal=fminsearch(@(p) rotafitMV1(p,pfit_mal(:,1),rotamalVS3_vacc,rotamalVS3_unvacc,malpop,malawi_cbr,repeffV3,vcov_mavg),[.687; .687])
pvacc1nr_mal=fminsearch(@(p) rotafitMV1nr(p,pfit_mal(:,1),rotamalVS3_vacc,rotamalVS3_unvacc,malpop,malawi_cbr,repeffV3,vcov_mavg),[.687; .687])

[LLvacfit_bc,Hvacfit_bc]=rotafitMV1(pvacc1_mal,pfit_mal(:,1),rotamalVS3_vacc,rotamalVS3_unvacc,malpop,malawi_cbr,repeffV3,vcov_mavg);
[LLvacfit_nr,Hvacfit_nr]=rotafitMV1nr(pvacc1nr_mal,pfit_mal(:,1),rotamalVS3_vacc,rotamalVS3_unvacc,malpop,malawi_cbr,repeffV3,vcov_mavg);
