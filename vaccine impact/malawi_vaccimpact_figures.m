%% Vaccine coverage with 95% CI (Figure 1)

figure; hold on 
plot(datenum(malPVdate),vcovPV_mavg1(:,1),'b')
plot(datenum(malPVdate),vcovPV_mavg1(:,2),'r')
fill([datenum(malPVdate); flipud(datenum(malPVdate)); datenum(malPVdate(1,:))],[vcovPV_dose2ci(:,1); flipud(vcovPV_dose2ci(:,2)); vcovPV_dose2ci(1,1)],[1 .7 .7],'EdgeColor',[1 .7 .7])
fill([datenum(malPVdate); flipud(datenum(malPVdate)); datenum(malPVdate(1,:))],[vcovPV_dose1ci(:,1); flipud(vcovPV_dose1ci(:,2)); vcovPV_dose1ci(1,1)],[.7 .7 1],'EdgeColor',[.7 .7 1])
plot(datenum(malPVdate),vcovPV_mavg1(:,1),'b')
plot(datenum(malPVdate),vcovPV_mavg1(:,2),'r')
legend('1st dose','2nd dose')
ylabel('Proportion vaccinated')
datetick('x','mmm-yy')

%% Model 1-4 predictions (Figure 2)

figure
subplot(2,2,1)
hold on
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 15; 15; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi(tvacc:end,1); Hpi(end:-1:tvacc,2); Hpi(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 15],'--k')
plot(malweekPV,sum(rotamalPV,2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp0(tvacc:end,:)+Hsamp1(tvacc:end,:)+Hsamp2(tvacc:end,:),'r')
plot(datesim(tvacc:end),Hsamp0_novacc(tvacc:end,:)+Hsamp1_novacc(tvacc:end,:)+Hsamp2_novacc(tvacc:end,:),'b')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
ylabel('Number of RVGE cases (per week)')
title('Homogeneity in vaccine response (no waning)')

subplot(2,2,2)
hold on
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 15; 15; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi_nr(tvacc:end,1); Hpi_nr(end:-1:tvacc,2); Hpi_nr(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 15],'--k')
plot(malweekPV,sum(rotamalPV,2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp0_novacc_nr(tvacc:end,:)+Hsamp1_novacc_nr(tvacc:end,:)+Hsamp2_novacc_nr(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp0_nr(tvacc:end,:)+Hsamp1_nr(tvacc:end,:)+Hsamp2_nr(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
ylabel('Number of RVGE cases (per week)')
title('Heterogeneity in vaccine response (no waning)')

subplot(2,2,3)
hold on
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 15; 15; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi_wane(tvacc:end,1); Hpi_wane(end:-1:tvacc,2); Hpi_wane(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 15],'--k')
plot(malweekPV,sum(rotamalPV,2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp0_novacc_wane(tvacc:end,:)+Hsamp1_novacc_wane(tvacc:end,:)+Hsamp2_novacc_wane(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp0_wane(tvacc:end,:)+Hsamp1_wane(tvacc:end,:)+Hsamp2_wane(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
ylabel('Number of RVGE cases (per week)')
title({'Homogeneity in vaccine response';'and waning of vaccine-induced immunity'})

subplot(2,2,4)
hold on
plot(malweekPV,sum(rotamalPV,2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp0(tvacc:end,1)+Hsamp1(tvacc:end,1)+Hsamp2(tvacc:end,1),'r')
plot(datesim(tvacc:end),Hsamp0_novacc(tvacc:end,1)+Hsamp1_novacc(tvacc:end,1)+Hsamp2_novacc(tvacc:end,1),'b')
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 15; 15; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi_wane_nr(tvacc:end,1); Hpi_wane_nr(end:-1:tvacc,2); Hpi_wane_nr(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 15],'--k')
plot(malweekPV,sum(rotamalPV,2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp0_novacc_wane_nr(tvacc:end,:)+Hsamp1_novacc_wane_nr(tvacc:end,:)+Hsamp2_novacc_wane_nr(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp0_wane_nr(tvacc:end,:)+Hsamp1_wane_nr(tvacc:end,:)+Hsamp2_wane_nr(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
ylabel('Number of RVGE cases (per week)')
title({'Heterogeneity in vaccine response';'and waning of vaccine-induced immunity'})
legend('Observed','Model-predicted','Model-predicted (no vaccination)','Orientation','Horizontal','Location','SO')

%% Overall effectiveness by year post-vaccination (Figure 3)

figure
subplot(2,2,1)
hold on
fill([2012:2022 2022:-1:2012 2012],100*[OEyr_ci(:,1); flipud(OEyr_ci(:,2)); OEyr_ci(1,1)],[.7 .7 .7],'EdgeColor',[.7 .7 .7])
plot(2012:2022,zeros(11,1),'--k')
plot(2012:2022,100*OEyr_obs(:,4),'k','LineWidth',2)
ylim([-20 100])
ylabel('Overall effectiveness (%)')
title('<5 years old')

subplot(2,2,2)
hold on
fill([2012:2022 2022:-1:2012 2012],100*[OEyr_ci0(:,1); flipud(OEyr_ci0(:,2)); OEyr_ci0(1,1)],[.7 .7 1],'EdgeColor',[.7 .7 1])
plot(2012:2022,zeros(11,1),'--k')
plot(2012:2022,100*OEyr_obs(:,1),'b','LineWidth',2)
ylim([-20 100])
ylabel('Overall effectiveness (%)')
title('<1 year old')

subplot(2,2,3)
hold on
fill([2012:2022 2022:-1:2012 2012],100*[OEyr_ci1(:,1); flipud(OEyr_ci1(:,2)); OEyr_ci1(1,1)],[.7 1 .7],'EdgeColor',[.7 1 .7])
plot(2012:2022,zeros(11,1),'--k')
plot(2012:2022,100*OEyr_obs(:,2),'Color',[0 .7 0],'LineWidth',2)
ylabel('Overall effectiveness (%)')
title('1-<2 years old')

subplot(2,2,4)
hold on
fill([2012:2022 2022:-1:2012 2012],100*[OEyr_ci2(:,1); flipud(OEyr_ci2(:,2)); OEyr_ci2(1,1)],[1 .7 .7],'EdgeColor',[1 .7 .7])
plot(2012:2022,zeros(11,1),'--k')
plot(2012:2022,100*OEyr_obs(:,3),'r','LineWidth',2)
ylabel('Overall effectiveness (%)')
title('2-<5 years old')
%%
gtext('(a)')
gtext('(b)')
gtext('(c)')
gtext('(d)')

%% Model fit to pre-vaccination data (Figure S3)
fig=figure; 
set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]);

subplot(3,1,1)
hold on
fill([datesim(1:653); flipud(datesim(1:653)); datesim(1)],[prctile(Hobserr(1:653,:)',2.5)'; flipud(prctile(Hobserr(1:653,:)',97.5)'); prctile(Hobserr(1,:),2.5)],[.7 .7 1],'EdgeColor',[.7 .7 1])
plot(malweeknum,sum(rotamalT2,2),'k') %'Color',[.5 .5 .5])
plot(malweek08,sum(rotamal08T2,2),'k') %'Color',[.5 .5 .5])
plot(datesim(1:653),Hsamp_novacc(1:653,:),'b')
set(gca,'FontSize',8)
datetick('x','yyyy')
xlim([min(malweeknum)-7 max(malweek08)+7])
ylim([-.1 20])
xlabel('Year','FontSize',8)
ylabel({'Number of RVGE cases'; '(per week)'},'FontSize',8)

%agedistmal2=(sum(rotamalT2)+sum(rotamal08T2))/sum(sum(rotamalT2)+sum(rotamal08T2));
%agedistHmalfit9=(sum(Hmalfit9)+sum(Hmalfit9_08))/sum(sum(Hmalfit9)+sum(Hmalfit9_08));
subplot(3,1,2)
y=bar([agedistmal2' agedistHmalfit9']); 
set(y(1),'FaceColor','k')
set(y(2),'FaceColor',[.7 .7 1])
set(gca,'FontSize',8,'XLim',[0 28],'XTick',1:2:27,'XTickLabel',agecatlabel2(1:2:27))
xlabel('Age group','FontSize',8)
legend('Observed','Predicted')
ylabel('Proportion of cases','FontSize',8)

subplot(3,1,3)
hold on
yyaxis left
set(gca,'FontSize',8)
plot(malweeknum,sum(negtest_age,2),'Color',[.7 .7 .7])
plot(malweek08,sum(negtest08_age,2),'-','Color',[.7 .7 .7])
xlim([min(malweeknum)-7 max(malweek08)+7])
ylim([-.1 35])
ylabel({'Number of rotavirus-negative';'cases (per week)'},'FontSize',8)
yyaxis right
set(gca,'FontSize',8)
plot(malweeknum,negtest_movavg/mean(negtest_movavg),'k')
plot(malweek08,negtest08_movavg/mean(negtest_movavg),'k')
plot([malweeknum(1) malweek08(end)],[1 1],'--k')
ylim([0 2])
datetick('x','yyyy')
xlim([min(malweeknum)-7 max(malweek08)+7])
ylabel('Relative reporting effort','FontSize',8)
%%
gtext('(a)')
gtext('(b)')
gtext('(c)')

%% Rotavirus-positive and negative cases and reporting effort (Figure S3)
fig=figure; 
set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]);

yyaxis left
hold on 
plot(malweekPV,sum(negtestPV,2),'Color',[.7 .7 .7])
plot(malweekPV,sum(rotamalPV,2),'-k')
ylabel({'No. of rotavirus-positive and negative cases';'(per week)'})

yyaxis right
plot(malweekPV,mean(repeffPV,2),'Color',[.5 0 .8],'LineWidth',1)
plot(malweekPV,ones(length(malweekPV),1),'--k')
datetick('x','mmm-yy')
xlim([malweekPV(1)-7 malweekPV(end)+7])
legend('RV-negative','RV-positive','reporting effort')
ylabel('Relative reporting effort')

%% Model 1-4 predictions by age (Figure S4)

figure
subplot(4,3,1)
hold on
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 10; 10; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi0(tvacc:end,1); Hpi0(end:-1:tvacc,2); Hpi0(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 10],'--k')
plot(malweekPV,sum(rotamalPV(:,1:12),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp0_novacc(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp0(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
ylabel({'Number of RVGE cases';'(per week)'})
title('<1 yr olds')

subplot(4,3,2)
hold on
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 6; 6; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi1(tvacc:end,1); Hpi1(end:-1:tvacc,2); Hpi1(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 6],'--k')
plot(malweekPV,sum(rotamalPV(:,13:au2),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp1_novacc(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp1(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
title('1-<2 yr olds')

subplot(4,3,3)
hold on
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 4; 4; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi2(tvacc:end,1); Hpi2(end:-1:tvacc,2); Hpi2(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 4],'--k')
plot(malweekPV,sum(rotamalPV(:,au2+1:au5),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp2_novacc(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp2(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
title('2-<5 yr olds')

subplot(4,3,4)
hold on
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 10; 10; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi0_nr(tvacc:end,1); Hpi0_nr(end:-1:tvacc,2); Hpi0_nr(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 10],'--k')
plot(malweekPV,sum(rotamalPV(:,1:12),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp0_novacc_nr(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp0_nr(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
ylabel({'Number of RVGE cases';'(per week)'})
title('<1 yr olds')

subplot(4,3,5)
hold on
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 6; 6; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi1_nr(tvacc:end,1); Hpi1_nr(end:-1:tvacc,2); Hpi1_nr(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 6],'--k')
plot(malweekPV,sum(rotamalPV(:,13:au2),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp1_novacc_nr(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp1_nr(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
title('1-<2 yr olds')

subplot(4,3,6)
hold on
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 4; 4; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi2_nr(tvacc:end,1); Hpi2_nr(end:-1:tvacc,2); Hpi2_nr(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 4],'--k')
plot(malweekPV,sum(rotamalPV(:,au2+1:au5),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp2_novacc_nr(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp2_nr(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
title('2-<5 yr olds')

subplot(4,3,7)
hold on
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 10; 10; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi0_wane(tvacc:end,1); Hpi0_wane(end:-1:tvacc,2); Hpi0_wane(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 10],'--k')
plot(malweekPV,sum(rotamalPV(:,1:12),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp0_novacc_wane(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp0_wane(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
ylabel({'Number of RVGE cases';'(per week)'})
title('<1 yr olds')

subplot(4,3,8)
hold on
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 6; 6; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi1_wane(tvacc:end,1); Hpi1_wane(end:-1:tvacc,2); Hpi1_wane(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 6],'--k')
plot(malweekPV,sum(rotamalPV(:,13:au2),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp1_novacc_wane(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp1_wane(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
title('1-<2 yr olds')

subplot(4,3,9)
hold on
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 4; 4; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi2_wane(tvacc:end,1); Hpi2_wane(end:-1:tvacc,2); Hpi2_wane(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 4],'--k')
plot(malweekPV,sum(rotamalPV(:,au2+1:au5),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp2_novacc_wane(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp2_wane(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
title('2-<5 yr olds')

subplot(4,3,10)
hold on
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 10; 10; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi0_wane_nr(tvacc:end,1); Hpi0_wane_nr(end:-1:tvacc,2); Hpi0_wane_nr(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 10],'--k')
plot(malweekPV,sum(rotamalPV(:,1:12),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp0_novacc_wane_nr(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp0_wane_nr(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
ylabel({'Number of RVGE cases';'(per week)'})
title('<1 yr olds')

subplot(4,3,11)
hold on
plot(malweekPV,sum(rotamalPV(:,13:au2),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp1_wane_nr(tvacc:end,1),'r')
plot(datesim(tvacc:end),Hsamp1_novacc_wane_nr(tvacc:end,1),'b')
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 6; 6; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi1_wane_nr(tvacc:end,1); Hpi1_wane_nr(end:-1:tvacc,2); Hpi1_wane_nr(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 6],'--k')
plot(malweekPV,sum(rotamalPV(:,13:au2),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp1_novacc_wane_nr(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp1_wane_nr(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
title('1-<2 yr olds')
legend('Observed','Model predicted','Model predicted (no vaccination)','Orientation','Horizontal','Location','SO')

subplot(4,3,12)
hold on
fill([datesim(tvacc+295); datesim(end); datesim(end); datesim(tvacc+295); datesim(tvacc+295)],[.01; .01; 4; 4; .01],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
fill([datesim(tvacc:end); datesim(end:-1:tvacc); datesim(tvacc)],[Hpi2_wane_nr(tvacc:end,1); Hpi2_wane_nr(end:-1:tvacc,2); Hpi2_wane_nr(tvacc,1)],[1 .8 .8],'EdgeColor',[1 .7 .7])
plot([datenum([2012 10 29]) datenum([2012 10 29])],[0 4],'--k')
plot(malweekPV,sum(rotamalPV(:,au2+1:au5),2),'Color',[.5 .5 .5])
plot(datesim(tvacc:end),Hsamp2_novacc_wane_nr(tvacc:end,:),'b')
plot(datesim(tvacc:end),Hsamp2_wane_nr(tvacc:end,:),'r')
plot([datesim(tvacc-1) datesim(end)+7],[0 0],'k')
datetick('x','yyyy')
xlim([datesim(tvacc-1) datesim(end)+7])
title('2-<5 yr olds')

%%
gtext('(a)')
gtext('(b)')
gtext('(c)')
gtext('(d)')

%% Annual rotavirus-positive and negative cases by age

figure 
subplot(1,2,1)
bar([mean(Hyr_obsT2([2:10 12:13],1:3)); Hyr_obsPV(:,1:3)],'stacked')
set(gca,'XTickLabel',{'Pre-vaccination',2012:2022})
ylabel('Rotavirus-positive cases')
%ylim([0 700])
%legend('<1 yr old','1-<2 yrs old','2-<5 yrs old')
subplot(1,2,2)
bar([mean(negtestT2_yr([2:10 12:13],1:3)); negtestPV_yr(:,1:3)],'stacked')
set(gca,'XTickLabel',{'Pre-vaccination',2012:2022})
ylabel('Rotavirus-negative cases')
legend('<1 yr old','1-<2 yrs old','2-<5 yrs old')

%%
gtext('(a)')
gtext('(b)')
