%% Generate model predictions + observation process error (Poisson)

for i=1:1000
    
    Hobserror0_nofit(:,i)=poissrnd(Hsamp0_nofit(:,i));
    Hobserror1_nofit(:,i)=poissrnd(Hsamp1_nofit(:,i));
    Hobserror2_nofit(:,i)=poissrnd(Hsamp2_nofit(:,i));
    
    Hobserror0_nfnr(:,i)=poissrnd(Hsamp0_nfnr(:,i));
    Hobserror1_nfnr(:,i)=poissrnd(Hsamp1_nfnr(:,i));
    Hobserror2_nfnr(:,i)=poissrnd(Hsamp2_nfnr(:,i));
    
    Hobserror0_bc(:,i)=poissrnd(Hsamp0_bc(:,i));
    Hobserror1_bc(:,i)=poissrnd(Hsamp1_bc(:,i));
    Hobserror2_bc(:,i)=poissrnd(Hsamp2_bc(:,i));
    
    Hobserror0_nr(:,i)=poissrnd(Hsamp0_nr(:,i));
    Hobserror1_nr(:,i)=poissrnd(Hsamp1_nr(:,i));
    Hobserror2_nr(:,i)=poissrnd(Hsamp2_nr(:,i));
    
    Hobserror0_wane(:,i)=poissrnd(Hsamp0_wane(:,i));
    Hobserror1_wane(:,i)=poissrnd(Hsamp1_wane(:,i));
    Hobserror2_wane(:,i)=poissrnd(Hsamp2_wane(:,i));
    
    Hobserror0_wane_nr(:,i)=poissrnd(Hsamp0_wane_nr(:,i));
    Hobserror1_wane_nr(:,i)=poissrnd(Hsamp1_wane_nr(:,i));
    Hobserror2_wane_nr(:,i)=poissrnd(Hsamp2_wane_nr(:,i));
    
    %HobserrorA(:,i)=nbinrnd(nbr,nbr./(nbr+(Hsamp0_bc(:,i).*repeff(:,1)+Hsamp1b(:,i).*repeff(:,13)+Hsamp2b(:,i)).*repeff(:,25)));
    %Hobserror0a(:,i)=nbinrnd(nbr,nbr./(nbr+Hsamp0_bc(:,i)));
    %Hobserror1a(:,i)=nbinrnd(nbr,nbr./(nbr+Hsamp1_bc(:,i)));
    %Hobserror2a(:,i)=nbinrnd(nbr,nbr./(nbr+Hsamp2_bc(:,i)));
end

%% Overall effectiveness for entire post-vaccination period

OE_samp=zeros(400,1);
OE_samp_age=zeros(400,3);
for i=1:100
    OE_samp(i,1)=1-sum(sum(rotamalPV(45:end,:),2))/sum(Hsamp_novacc(tvacc+1:end,i));
    OE_samp_age(i,1)=1-sum(sum(rotamalPV(45:end,1:12),2))/sum(Hsamp0_novacc(tvacc+1:end,i));
    OE_samp_age(i,2)=1-sum(sum(rotamalPV(45:end,13:24),2))/sum(Hsamp1_novacc(tvacc+1:end,i));
    OE_samp_age(i,3)=1-sum(sum(rotamalPV(45:end,25:27),2))/sum(Hsamp2_novacc(tvacc+1:end,i));

    OE_samp(i+100,1)=1-sum(sum(rotamalPV(45:end,:),2))/sum(Hsamp_novacc_nr(tvacc+1:end,i));
    OE_samp_age(i+100,1)=1-sum(sum(rotamalPV(45:end,1:12),2))/sum(Hsamp0_novacc_nr(tvacc+1:end,i));
    OE_samp_age(i+100,2)=1-sum(sum(rotamalPV(45:end,13:24),2))/sum(Hsamp1_novacc_nr(tvacc+1:end,i));
    OE_samp_age(i+100,3)=1-sum(sum(rotamalPV(45:end,25:27),2))/sum(Hsamp2_novacc_nr(tvacc+1:end,i));
    
    OE_samp(i+200,1)=1-sum(sum(rotamalPV(45:end,:),2))/sum(Hsamp_novacc_wane(tvacc+1:end,i));
    OE_samp_age(i+200,1)=1-sum(sum(rotamalPV(45:end,1:12),2))/sum(Hsamp0_novacc_wane(tvacc+1:end,i));
    OE_samp_age(i+200,2)=1-sum(sum(rotamalPV(45:end,13:24),2))/sum(Hsamp1_novacc_wane(tvacc+1:end,i));
    OE_samp_age(i+200,3)=1-sum(sum(rotamalPV(45:end,25:27),2))/sum(Hsamp2_novacc_wane(tvacc+1:end,i));
    
    OE_samp(i+300,1)=1-sum(sum(rotamalPV(45:end,:),2))/sum(Hsamp_novacc_wane_nr(tvacc+1:end,i));
    OE_samp_age(i+300,1)=1-sum(sum(rotamalPV(45:end,1:12),2))/sum(Hsamp0_novacc_wane_nr(tvacc+1:end,i));
    OE_samp_age(i+300,2)=1-sum(sum(rotamalPV(45:end,13:24),2))/sum(Hsamp1_novacc_wane_nr(tvacc+1:end,i));
    OE_samp_age(i+300,3)=1-sum(sum(rotamalPV(45:end,25:27),2))/sum(Hsamp2_novacc_wane_nr(tvacc+1:end,i));
end

%% Model-predicted cases (no vaccination) by year post-vaccination

Hyr_samp_novacc=zeros(11,400);
Hyr0_samp_novacc=zeros(11,400);
Hyr1_samp_novacc=zeros(11,400);
Hyr2_samp_novacc=zeros(11,400);
for y=1:11
    for j=1:length(rotamalPV)
        if malPVdate(j,1)==2011+y
            for i=1:100
                Hyr_samp_novacc(y,i)=Hyr_samp_novacc(y,i)+sum(Hsamp_novacc(tvacc+j,i));
                Hyr0_samp_novacc(y,i)=Hyr0_samp_novacc(y,i)+sum(Hsamp0_novacc(tvacc+j,i));
                Hyr1_samp_novacc(y,i)=Hyr1_samp_novacc(y,i)+sum(Hsamp1_novacc(tvacc+j,i));
                Hyr2_samp_novacc(y,i)=Hyr2_samp_novacc(y,i)+sum(Hsamp2_novacc(tvacc+j,i));
                            
                Hyr_samp_novacc(y,i+100)=Hyr_samp_novacc(y,i+100)+sum(Hsamp_novacc_nr(tvacc+j,i));
                Hyr0_samp_novacc(y,i+100)=Hyr0_samp_novacc(y,i+100)+sum(Hsamp0_novacc_nr(tvacc+j,i));
                Hyr1_samp_novacc(y,i+100)=Hyr1_samp_novacc(y,i+100)+sum(Hsamp1_novacc_nr(tvacc+j,i));
                Hyr2_samp_novacc(y,i+100)=Hyr2_samp_novacc(y,i+100)+sum(Hsamp2_novacc_nr(tvacc+j,i));
                                
                Hyr_samp_novacc(y,i+200)=Hyr_samp_novacc(y,i+200)+sum(Hsamp_novacc_wane(tvacc+j,i));
                Hyr0_samp_novacc(y,i+200)=Hyr0_samp_novacc(y,i+200)+sum(Hsamp0_novacc_wane(tvacc+j,i));
                Hyr1_samp_novacc(y,i+200)=Hyr1_samp_novacc(y,i+200)+sum(Hsamp1_novacc_wane(tvacc+j,i));
                Hyr2_samp_novacc(y,i+200)=Hyr2_samp_novacc(y,i+200)+sum(Hsamp2_novacc_wane(tvacc+j,i));
                                
                Hyr_samp_novacc(y,i+300)=Hyr_samp_novacc(y,i+300)+sum(Hsamp_novacc_wane_nr(tvacc+j,i));
                Hyr0_samp_novacc(y,i+300)=Hyr0_samp_novacc(y,i+300)+sum(Hsamp0_novacc_wane_nr(tvacc+j,i));
                Hyr1_samp_novacc(y,i+300)=Hyr1_samp_novacc(y,i+300)+sum(Hsamp1_novacc_wane_nr(tvacc+j,i));
                Hyr2_samp_novacc(y,i+300)=Hyr2_samp_novacc(y,i+300)+sum(Hsamp2_novacc_wane_nr(tvacc+j,i));
            end
        end
    end
end

%% Overall effectiveness by year post-vaccination

OE_samp_yr=zeros(400,11);
OE_samp0_yr=zeros(400,11);
OE_samp1_yr=zeros(400,11);
OE_samp2_yr=zeros(400,11);
for y=1:11
    for i=1:400
        OE_samp_yr(i,y)=1-Hyr_obsPV(y,4)/Hyr_samp_novacc(y,i);
        OE_samp0_yr(i,y)=1-Hyr_obsPV(y,1)/Hyr0_samp_novacc(y,i);
        OE_samp1_yr(i,y)=1-Hyr_obsPV(y,2)/Hyr1_samp_novacc(y,i);
        OE_samp2_yr(i,y)=1-Hyr_obsPV(y,3)/Hyr2_samp_novacc(y,i);
    end
end

OEyr_ci=prctile(OE_samp_yr,[2.5 97.5])';
OEyr_ci0=prctile(OE_samp0_yr,[2.5 97.5])';
OEyr_ci1=prctile(OE_samp1_yr,[2.5 97.5])';
OEyr_ci2=prctile(OE_samp2_yr,[2.5 97.5])';

%% Overall effectiveness for 2017-2022

for i=1:400
    OE_samp_post2017(i,1)=1-sum(Hyr_obsPV(6:end,4))/sum(Hyr_samp_novacc(6:end,i));
    OE_samp0_post2017(i,1)=1-sum(Hyr_obsPV(6:end,1))/sum(Hyr0_samp_novacc(6:end,i));
    OE_samp1_post2017(i,1)=1-sum(Hyr_obsPV(6:end,2))/sum(Hyr1_samp_novacc(6:end,i));
    OE_samp2_post2017(i,1)=1-sum(Hyr_obsPV(6:end,3))/sum(Hyr2_samp_novacc(6:end,i));
end
OE_post2017=prctile(OE_samp_post2017,[50 2.5 97.5]);
OEage_post2017=[prctile(OE_samp0_post2017,[50 2.5 97.5]); prctile(OE_samp1_post2017,[50 2.5 97.5]); prctile(OE_samp2_post2017,[50 2.5 97.5])];

%% Model-predicted overall effectiveness 

OE_pred=zeros(100,4);
OE_pred_age=zeros(100,3,4);
for i=1:100
    OE_pred(i,1)=1-sum(Hsamp(tvacc+1:end,i))/sum(Hsamp_novacc(tvacc+1:end,i));
    OE_pred_age(i,1)=1-sum(Hsamp0(tvacc+1:end,i))/sum(Hsamp0_novacc(tvacc+1:end,i));
    OE_pred_age(i,2)=1-sum(Hsamp1(tvacc+1:end,i))/sum(Hsamp1_novacc(tvacc+1:end,i));
    OE_pred_age(i,3)=1-sum(Hsamp2(tvacc+1:end,i))/sum(Hsamp2_novacc(tvacc+1:end,i));

    OE_pred(i,2)=1-sum(Hsamp_nr(tvacc+1:end,i))/sum(Hsamp_novacc_nr(tvacc+1:end,i));
    OE_pred_age(i,1,2)=1-sum(Hsamp0_nr(tvacc+1:end,i))/sum(Hsamp0_novacc_nr(tvacc+1:end,i));
    OE_pred_age(i,2,2)=1-sum(Hsamp1_nr(tvacc+1:end,i))/sum(Hsamp1_novacc_nr(tvacc+1:end,i));
    OE_pred_age(i,3,2)=1-sum(Hsamp2_nr(tvacc+1:end,i))/sum(Hsamp2_novacc_nr(tvacc+1:end,i));
    
    OE_pred(i,3)=1-sum(Hsamp_wane(tvacc+1:end,i))/sum(Hsamp_novacc_wane(tvacc+1:end,i));
    OE_pred_age(i,1,3)=1-sum(Hsamp0_wane(tvacc+1:end,i))/sum(Hsamp0_novacc_wane(tvacc+1:end,i));
    OE_pred_age(i,2,3)=1-sum(Hsamp1_wane(tvacc+1:end,i))/sum(Hsamp1_novacc_wane(tvacc+1:end,i));
    OE_pred_age(i,3,3)=1-sum(Hsamp2_wane(tvacc+1:end,i))/sum(Hsamp2_novacc_wane(tvacc+1:end,i));
    
    OE_pred(i,4)=1-sum(Hsamp_wane_nr(tvacc+1:end,i))/sum(Hsamp_novacc_wane_nr(tvacc+1:end,i));
    OE_pred_age(i,1,4)=1-sum(Hsamp0_wane_nr(tvacc+1:end,i))/sum(Hsamp0_novacc_wane_nr(tvacc+1:end,i));
    OE_pred_age(i,2,4)=1-sum(Hsamp1_wane_nr(tvacc+1:end,i))/sum(Hsamp1_novacc_wane_nr(tvacc+1:end,i));
    OE_pred_age(i,3,4)=1-sum(Hsamp2_wane_nr(tvacc+1:end,i))/sum(Hsamp2_novacc_wane_nr(tvacc+1:end,i));
end

mean(OE_pred)
prctile(OE_pred,[2.5 97.5])

mean(OE_pred_age)
prctile(OE_pred_age,[2.5 97.5])

%% Spearman's rank correlation between observed and predicted cases

pvi=[1:431 459:547]; % index for Jan 2012-Jun 2022 excluding Apr-Oct 2020
insi=1:295; % index for in-sample validation period
outsi=[296:431 459:547]; % index for in-sample validation period

for i=1:100
    %r_samp(i)=corr(Hsamp(tvacc+pvi,i),sum(rotamalPV(pvi,:),2),'Type','Spearman');
    %r_samp_nr(i)=corr(Hsamp_nr(tvacc+pvi,i),sum(rotamalPV(pvi,:),2),'Type','Spearman');
    %r_samp_wane(i)=corr(Hsamp_wane(tvacc+pvi,i),sum(rotamalPV(pvi,:),2),'Type','Spearman');
    %r_samp_wane_nr(i)=corr(Hsamp_wane_nr(tvacc+pvi,i),sum(rotamalPV(pvi,:),2),'Type','Spearman');

    %r_insamp(i)=corr(Hsamp(tvacc+insi,i),sum(rotamalPV(insi,:),2),'Type','Spearman');
    %r_insamp_nr(i)=corr(Hsamp_nr(tvacc+insi,i),sum(rotamalPV(insi,:),2),'Type','Spearman');
    %r_insamp_wane(i)=corr(Hsamp_wane(tvacc+insi,i),sum(rotamalPV(insi,:),2),'Type','Spearman');
    %r_insamp_wane_nr(i)=corr(Hsamp_wane_nr(tvacc+insi,i),sum(rotamalPV(insi,:),2),'Type','Spearman');
    
    %r_outsamp(i)=corr(Hsamp(tvacc+outsi,i),sum(rotamalPV(outsi,:),2),'Type','Spearman');
    %r_outsamp_nr(i)=corr(Hsamp_nr(tvacc+outsi,i),sum(rotamalPV(outsi,:),2),'Type','Spearman');
    %r_outsamp_wane(i)=corr(Hsamp_wane(tvacc+outsi,i),sum(rotamalPV(outsi,:),2),'Type','Spearman');
    %r_outsamp_wane_nr(i)=corr(Hsamp_wane_nr(tvacc+outsi,i),sum(rotamalPV(outsi,:),2),'Type','Spearman');
    
    r_samp(i)=corr([Hsamp0(tvacc+pvi,i); Hsamp1(tvacc+pvi,i); Hsamp2(tvacc+pvi,i)],[sum(rotamalPV(pvi,1:12),2); sum(rotamalPV(pvi,13:24),2); sum(rotamalPV(pvi,25:27),2)],'Type','Spearman');
    r_samp_nr(i)=corr([Hsamp0_nr(tvacc+pvi,i); Hsamp1_nr(tvacc+pvi,i); Hsamp2_nr(tvacc+pvi,i)],[sum(rotamalPV(pvi,1:12),2); sum(rotamalPV(pvi,13:24),2); sum(rotamalPV(pvi,25:27),2)],'Type','Spearman');
    r_samp_wane(i)=corr([Hsamp0_wane(tvacc+pvi,i); Hsamp1_wane(tvacc+pvi,i); Hsamp2_wane(tvacc+pvi,i)],[sum(rotamalPV(pvi,1:12),2); sum(rotamalPV(pvi,13:24),2); sum(rotamalPV(pvi,25:27),2)],'Type','Spearman');
    r_samp_wane_nr(i)=corr([Hsamp0_wane_nr(tvacc+pvi,i); Hsamp1_wane_nr(tvacc+pvi,i); Hsamp2_wane_nr(tvacc+pvi,i)],[sum(rotamalPV(pvi,1:12),2); sum(rotamalPV(pvi,13:24),2); sum(rotamalPV(pvi,25:27),2)],'Type','Spearman');

    r_insamp(i)=corr([Hsamp0(tvacc+insi,i); Hsamp1(tvacc+insi,i); Hsamp2(tvacc+insi,i)],[sum(rotamalPV(insi,1:12),2); sum(rotamalPV(insi,13:24),2); sum(rotamalPV(insi,25:27),2)],'Type','Spearman');
    r_insamp_nr(i)=corr([Hsamp0_nr(tvacc+insi,i); Hsamp1_nr(tvacc+insi,i); Hsamp2_nr(tvacc+insi,i)],[sum(rotamalPV(insi,1:12),2); sum(rotamalPV(insi,13:24),2); sum(rotamalPV(insi,25:27),2)],'Type','Spearman');
    r_insamp_wane(i)=corr([Hsamp0_wane(tvacc+insi,i); Hsamp1_wane(tvacc+insi,i); Hsamp2_wane(tvacc+insi,i)],[sum(rotamalPV(insi,1:12),2); sum(rotamalPV(insi,13:24),2); sum(rotamalPV(insi,25:27),2)],'Type','Spearman');
    r_insamp_wane_nr(i)=corr([Hsamp0_wane_nr(tvacc+insi,i); Hsamp1_wane_nr(tvacc+insi,i); Hsamp2_wane_nr(tvacc+insi,i)],[sum(rotamalPV(insi,1:12),2); sum(rotamalPV(insi,13:24),2); sum(rotamalPV(insi,25:27),2)],'Type','Spearman');

    r_outsamp(i)=corr([Hsamp0(tvacc+outsi,i); Hsamp1(tvacc+outsi,i); Hsamp2(tvacc+outsi,i)],[sum(rotamalPV(outsi,1:12),2); sum(rotamalPV(outsi,13:24),2); sum(rotamalPV(outsi,25:27),2)],'Type','Spearman');
    r_outsamp_nr(i)=corr([Hsamp0_nr(tvacc+outsi,i); Hsamp1_nr(tvacc+outsi,i); Hsamp2_nr(tvacc+outsi,i)],[sum(rotamalPV(outsi,1:12),2); sum(rotamalPV(outsi,13:24),2); sum(rotamalPV(outsi,25:27),2)],'Type','Spearman');
    r_outsamp_wane(i)=corr([Hsamp0_wane(tvacc+outsi,i); Hsamp1_wane(tvacc+outsi,i); Hsamp2_wane(tvacc+outsi,i)],[sum(rotamalPV(outsi,1:12),2); sum(rotamalPV(outsi,13:24),2); sum(rotamalPV(outsi,25:27),2)],'Type','Spearman');
    r_outsamp_wane_nr(i)=corr([Hsamp0_wane_nr(tvacc+outsi,i); Hsamp1_wane_nr(tvacc+outsi,i); Hsamp2_wane_nr(tvacc+outsi,i)],[sum(rotamalPV(outsi,1:12),2); sum(rotamalPV(outsi,13:24),2); sum(rotamalPV(outsi,25:27),2)],'Type','Spearman');

end
mean([r_samp; r_samp_nr; r_samp_wane; r_samp_wane_nr],2)
mean([r_insamp; r_insamp_nr; r_insamp_wane; r_insamp_wane_nr],2)
mean([r_outsamp; r_outsamp_nr; r_outsamp_wane; r_outsamp_wane_nr],2)

%% Mean squared error for in-sample and out-of-sample periods and overall

for i=1:100
    %rmse_samp(i)=rmse(Hsamp(tvacc+pvi,i),sum(rotamalPV(pvi,:),2),);
    %rmse_samp_nr(i)=rmse(Hsamp_nr(tvacc+pvi,i),sum(rotamalPV(pvi,:),2));
    %rmse_samp_wane(i)=rmse(Hsamp_wane(tvacc+pvi,i),sum(rotamalPV(pvi,:),2));
    %rmse_samp_wane_nr(i)=rmse(Hsamp_wane_nr(tvacc+pvi,i),sum(rotamalPV(pvi,:),2));

    %rmse_insamp(i)=rmse(Hsamp(tvacc+insi,i),sum(rotamalPV(insi,:),2));
    %rmse_insamp_nr(i)=rmse(Hsamp_nr(tvacc+insi,i),sum(rotamalPV(insi,:),2));
    %rmse_insamp_wane(i)=rmse(Hsamp_wane(tvacc+insi,i),sum(rotamalPV(insi,:),2));
    %rmse_insamp_wane_nr(i)=rmse(Hsamp_wane_nr(tvacc+insi,i),sum(rotamalPV(insi,:),2));
    
    %rmse_outsamp(i)=corr(Hsamp(tvacc+outsi,i),sum(rotamalPV(outsi,:),2));
    %rmse_outsamp_nr(i)=corr(Hsamp_nr(tvacc+outsi,i),sum(rotamalPV(outsi,:),2));
    %rmse_outsamp_wane(i)=corr(Hsamp_wane(tvacc+outsi,i),sum(rotamalPV(outsi,:),2));
    %rmse_outsamp_wane_nr(i)=corr(Hsamp_wane_nr(tvacc+outsi,i),sum(rotamalPV(outsi,:),2));
    
    rmse_samp(i)=rmse([Hsamp0(tvacc+pvi,i); Hsamp1(tvacc+pvi,i); Hsamp2(tvacc+pvi,i)],[sum(rotamalPV(pvi,1:12),2); sum(rotamalPV(pvi,13:24),2); sum(rotamalPV(pvi,25:27),2)]);
    rmse_samp_nr(i)=rmse([Hsamp0_nr(tvacc+pvi,i); Hsamp1_nr(tvacc+pvi,i); Hsamp2_nr(tvacc+pvi,i)],[sum(rotamalPV(pvi,1:12),2); sum(rotamalPV(pvi,13:24),2); sum(rotamalPV(pvi,25:27),2)]);
    rmse_samp_wane(i)=rmse([Hsamp0_wane(tvacc+pvi,i); Hsamp1_wane(tvacc+pvi,i); Hsamp2_wane(tvacc+pvi,i)],[sum(rotamalPV(pvi,1:12),2); sum(rotamalPV(pvi,13:24),2); sum(rotamalPV(pvi,25:27),2)]);
    rmse_samp_wane_nr(i)=rmse([Hsamp0_wane_nr(tvacc+pvi,i); Hsamp1_wane_nr(tvacc+pvi,i); Hsamp2_wane_nr(tvacc+pvi,i)],[sum(rotamalPV(pvi,1:12),2); sum(rotamalPV(pvi,13:24),2); sum(rotamalPV(pvi,25:27),2)]);

    rmse_insamp(i)=rmse([Hsamp0(tvacc+insi,i); Hsamp1(tvacc+insi,i); Hsamp2(tvacc+insi,i)],[sum(rotamalPV(insi,1:12),2); sum(rotamalPV(insi,13:24),2); sum(rotamalPV(insi,25:27),2)]);
    rmse_insamp_nr(i)=rmse([Hsamp0_nr(tvacc+insi,i); Hsamp1_nr(tvacc+insi,i); Hsamp2_nr(tvacc+insi,i)],[sum(rotamalPV(insi,1:12),2); sum(rotamalPV(insi,13:24),2); sum(rotamalPV(insi,25:27),2)]);
    rmse_insamp_wane(i)=rmse([Hsamp0_wane(tvacc+insi,i); Hsamp1_wane(tvacc+insi,i); Hsamp2_wane(tvacc+insi,i)],[sum(rotamalPV(insi,1:12),2); sum(rotamalPV(insi,13:24),2); sum(rotamalPV(insi,25:27),2)]);
    rmse_insamp_wane_nr(i)=rmse([Hsamp0_wane_nr(tvacc+insi,i); Hsamp1_wane_nr(tvacc+insi,i); Hsamp2_wane_nr(tvacc+insi,i)],[sum(rotamalPV(insi,1:12),2); sum(rotamalPV(insi,13:24),2); sum(rotamalPV(insi,25:27),2)]);

    rmse_outsamp(i)=rmse([Hsamp0(tvacc+outsi,i); Hsamp1(tvacc+outsi,i); Hsamp2(tvacc+outsi,i)],[sum(rotamalPV(outsi,1:12),2); sum(rotamalPV(outsi,13:24),2); sum(rotamalPV(outsi,25:27),2)]);
    rmse_outsamp_nr(i)=rmse([Hsamp0_nr(tvacc+outsi,i); Hsamp1_nr(tvacc+outsi,i); Hsamp2_nr(tvacc+outsi,i)],[sum(rotamalPV(outsi,1:12),2); sum(rotamalPV(outsi,13:24),2); sum(rotamalPV(outsi,25:27),2)]);
    rmse_outsamp_wane(i)=rmse([Hsamp0_wane(tvacc+outsi,i); Hsamp1_wane(tvacc+outsi,i); Hsamp2_wane(tvacc+outsi,i)],[sum(rotamalPV(outsi,1:12),2); sum(rotamalPV(outsi,13:24),2); sum(rotamalPV(outsi,25:27),2)]);
    rmse_outsamp_wane_nr(i)=rmse([Hsamp0_wane_nr(tvacc+outsi,i); Hsamp1_wane_nr(tvacc+outsi,i); Hsamp2_wane_nr(tvacc+outsi,i)],[sum(rotamalPV(outsi,1:12),2); sum(rotamalPV(outsi,13:24),2); sum(rotamalPV(outsi,25:27),2)]);

end
mean([rmse_samp; rmse_samp_nr; rmse_samp_wane; rmse_samp_wane_nr],2)
mean([rmse_insamp; rmse_insamp_nr; rmse_insamp_wane; rmse_insamp_wane_nr],2)
mean([rmse_outsamp; rmse_outsamp_nr; rmse_outsamp_wane; rmse_outsamp_wane_nr],2)
