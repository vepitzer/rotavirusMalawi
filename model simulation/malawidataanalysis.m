load('malawidata.mat')

%% Fit "constant incidence" Poisson model to the pre-vaccination data
blantyre_wkpop=interp1([malweeknum(1); datenum([blantyre_pop(:,1:2) ones(length(blantyre_pop),1) zeros(length(blantyre_pop),3)])],[800000; blantyre_pop(:,5)],[malweeknum; malweek08]);
agep=[(1/12)*malpop(1)*ones(1,12) (1/48)*malpop(2)*ones(1,12) (1/4)*malpop(2)*ones(1,3) malpop(3:17).*ones(1,15)];
agep=agep/sum(agep); %Proportion of pop'n in each age group 

for j=1:27
simple_coeff(j,:)=glmfit([malweeknum-malweeknum(1); malweek08-malweeknum(1)],[rotamalT2(:,j); rotamal08T2(:,j)],'Poisson');
simple_predict(:,j)=glmval(simple_coeff(j,:)',[malweeknum-malweeknum(1); malweek08-malweeknum(1)],'log');
simple_rep(:,j)=simple_predict(:,j).*[negtest_movavg/mean([negtest_movavg; negtest08_movavg]); negtest08_movavg/mean([negtest_movavg; negtest08_movavg])]; 
end
simple_coeff_all=glmfit([malweeknum-malweeknum(1); malweek08-malweeknum(1)],sum([rotamalT2; rotamal08T2],2),'Poisson');
simple_predict_all=glmval(simple_coeff_all,[malweeknum-malweeknum(1); malweek08-malweeknum(1)],'log');
simple_rep_all=(simple_predict_all.*[negtest_movavg/mean([negtest_movavg; negtest08_movavg]); negtest08_movavg/mean([negtest_movavg; negtest08_movavg])])*(agep(1:27)/sum(agep(1:27))); 

%%
LL_simple=0;
for i=1:size(rotamalT2,1)
    for j=1:size(rotamalT2,2)
        LL_simple=LL_simple+log(poisspdf(rotamalT2(i,j),simple_rep(i,j)));
    end
end
for i=1:size(rotamal08T2,1)
    for j=1:size(rotamal08T2,2)
        LL_simple=LL_simple+log(poisspdf(rotamal08T2(i,j),simple_rep(length(rotamalT2)+i,j)));
    end
end
LL_simple=-LL_simple;
%%
LL_simple_all=0;
for i=1:size(rotamalT2,1)
    for j=1:size(rotamalT2,2)
        LL_simple_all=LL_simple_all+log(poisspdf(rotamalT2(i,j),simple_rep_all(i,j)));
    end
end
for i=1:size(rotamal08T2,1)
    for j=1:size(rotamal08T2,2)
        LL_simple_all=LL_simple_all+log(poisspdf(rotamal08T2(i,j),simple_rep_all(length(rotamalT2)+i,j)));
    end
end
LL_simple_all=-LL_simple_all;

%%
mdl_GLMall=GeneralizedLinearModel.fit([malweeknum; malweek08],[sum(rotamalT2,2); sum(rotamal08T2,2)],'linear','Distribution','poisson','Offset',blantyre_wkpop)

%% THE CODE BELOW WAS USED TO CALCULATE THE OBSERVED VACCINE EFFECTIVENESS FROM THE CASE-CONTROL STUDY DATA
% HOWEVER, THE INDIVIDUAL-LEVEL DATA NECESSARY TO RUN THE POISSON
% REGRESSION MODELS CANNOT BE SHARED DUE TO PRIVACY CONCERNS

% Vaccine effectiveness by age <1 or 1 yo
x1=~malposVS3(:,3).*malVS3_veligible.*malVS3_vaccard.*(malVS3_vacc(:,1)==malVS3_vacc(:,2)).*(malageVS3(:,1)<12);
x2=~malposVS3(:,3).*malVS3_veligible.*malVS3_vaccard.*(malVS3_vacc(:,1)==malVS3_vacc(:,2)).*(malageVS3(:,1)>=12).*(malageVS3(:,1)<24);

glmfit([malVS3_vacc(find(x1),2) malageVS3(find(x1),1) malVS3_admdate(find(x1),1:2)],malposVS3(find(x1),1),'binomial','link','logit');
exp(ans(2:4));
1-ans(1)

glmfit([malVS3_vacc(find(x2),2) malageVS3(find(x2),1) malVS3_admdate(find(x2),1:2)],malposVS3(find(x2),1),'binomial','link','logit');
exp(ans(2:4));
1-ans(1)

mdl_age1=GeneralizedLinearModel.fit([malVS3_vacc(find(x1),2) malageVS3(find(x1),1) malVS3_admdate(find(x1),1:2)],malposVS3(find(x1),1),'linear','Distribution','binomial','Link','logit')
mdl_age2=GeneralizedLinearModel.fit([malVS3_vacc(find(x2),2) malageVS3(find(x2),1) malVS3_admdate(find(x2),1:2)],malposVS3(find(x2),1),'linear','Distribution','binomial','Link','logit')

mdl1ci=coefCI(mdl_age1);
mdl2ci=coefCI(mdl_age2);
1-exp(mdl1ci(2,:))
1-exp(mdl2ci(2,:))


%% Vaccine effectiveness against any RVGE, all ages 
x1=~malposVS3(:,3).*malVS3_veligible.*malVS3_vaccard.*(malVS3_vacc(:,1)==malVS3_vacc(:,2));

glmfit([malVS3_vacc(find(x1),2) malageVS3(find(x1),1) malVS3_admdate(find(x1),1:2)],malposVS3(find(x1),1),'binomial','link','logit')
exp(ans(2:4))
1-ans(1)

mdl_all1=GeneralizedLinearModel.fit([malVS3_vacc(find(x1),2) malageVS3(find(x1),1) malVS3_admdate(find(x1),1:2)],malposVS3(find(x1),1),'linear','Distribution','binomial','Link','logit')
mdl1ci=coefCI(mdl_all1);
1-exp(mdl1ci(2,:))


