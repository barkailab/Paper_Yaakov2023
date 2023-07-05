% combined
clear all
load('TFdata.mat')
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Rotem/promoterLengthsORF.mat','promoterLengthsORF')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/promoterIDXvecFimp.mat','promoterIDXvecF')
metaPro=metaProfilePromLenDivya(normProfile,'promEnd','position','afterTss',0,'promLen',promoterLengthsORF);
sumProm=squeeze(sum(metaPro,2,'omitnan'));
sumProm(isnan(promoterLengthsORF),:)=NaN;
metaPro=metaProfilePromLenDivya(promoterIDXvecF,'promEnd','position','afterTss',0,'promLen',promoterLengthsORF);
promType=mode(metaPro,2)

cmbSmp=unique(sTable(sTable.spike==0&ismember(sTable.tp,[0,30,60,90]),{'tf','gt','tp'}))
[~,sTable.csId]=ismember(sTable(:,{'tf','gt','tp'}),cmbSmp);
cmbSmp.n=accumarray(sTable.csId(sTable.csId>0),sTable.spike(sTable.csId>0)==0);
cmbSmp.nSpike=accumarray(sTable.csId(sTable.csId>0),sTable.spike(sTable.csId>0)==1);


for i=1:height(cmbSmp)
    selSmp=find(sTable.csId==i&sTable.spike==0)
    if numel(selSmp)>1
        cmbSmp.cr(i)=1-mean(squareform(1-corr(sumProm(promType<3,selSmp))));
    else
        cmbSmp.cr(i)=nan;
    end
    cmbProfile(:,i)=mean(normProfile(:,selSmp),2);
end
writetable(cmbSmp,'checSmp.xlsx')
metaPro=metaProfilePromLenDivya(cmbProfile,'promEnd','position','afterTss',0,'promLen',promoterLengthsORF);
cmbProm=squeeze(sum(metaPro,2,'omitnan'));
cmbProm(isnan(promoterLengthsORF),:)=NaN;
%% Figure 1B
figure('Color',[1 1 1],'Renderer','painters','Position',[2127 288 560 361])
selSmp=ismember(cmbSmp.tf,{'msn2','reb1','sok2','tbp','yap1'})%true(height(cmbSmp),1)
crAll=corr(cmbProm(promType<3,selSmp));
hold off
imagesc(crAll,[0 1])
set(gca,'ytick',1:numel(selSmp),'yticklabel',strcat(cmbSmp.tf(selSmp),'-',cmbSmp.gt(selSmp),'-',num2str(cmbSmp.tp(selSmp))));
hold on
plot([1:numel(selSmp)]+[0.5;0.5],ylim','-','Color',[1 1 1].*0.85)
plot(xlim',[1:numel(selSmp)]+[0.5;0.5],'-','Color',[1 1 1].*0.85)
bL=find(cmbSmp.tp(selSmp)==90)

plot(bL'+[0.5;0.5],ylim','-','Color',[1 1 1].*0.15,'Linewidth',1.5)
plot(xlim',bL'+[0.5;0.5],'-','Color',[1 1 1].*0.15,'Linewidth',1.5)

%text([1:sum(selSmp)],[1:sum(selSmp)],num2str(cmbSmp.cr(selSmp),'%.2f'),'HorizontalAlignment','center','FontSize',6)
colorbar()
title('sumProm correlation noSpike in means')
saveas(gcf,sprintf('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/svg/Fig6SumProm.fig'))
%% Figure 5C
load('dataFS-single.mat', 'nucMeta')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/bsAll.mat');
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/promoterIDXvecFimp.mat','promoterIDXvecF')
[~,nucIdx]=min(abs(bsAll.pos-nucMeta.pos'),[],2);
bsAll.gene=nucMeta.gene(nucIdx);

bsAll.pos=round(bsAll.pos);
bsAll.type=promoterIDXvecF(bsAll.pos);
bsAll.dir=ismember(bsAll.dir,'+')*2-1;

gtOrder={'aaino80','aaspt6','aaspt16','aachd1','aaspt6spt16'};
tfOrder={'msn2','reb1'};
gtColor=lines(numel(gtOrder))

g=3
t=1
selSmp=find(cmbSmp.tp<=90 &ismember(cmbSmp.gt,gtOrder(g))&ismember(cmbSmp.tf,tfOrder(t)))
subplot(1,1,1)
currBs=bsAll(contains(lower(bsAll.tf),tfOrder(t))&ismember(bsAll.type,[1,2])&(bsAll.gene>0),:);
[~,idx]=sortrows(table(currBs.type,-mean(reshape(cmbProfile(currBs.pos+[-50:50],selSmp(end)),height(currBs),101),2)))

bsSur=150;
figure('Color',[1 1 1],'Renderer','painters','Position',[2127 288 560 361])
subplot(1,1,1)
c=0;
for s=selSmp'
    c=c+1;
    subplot(4,numel(selSmp),c+[4,8,12])
    imageMat=reshape(cmbProfile(currBs.pos+[-bsSur:bsSur].*currBs.dir,s),height(currBs),2*bsSur+1);
    bL=sum(currBs.type==1)
    imagesc(movmean(imageMat(idx,:),5,2),'XData',[-1 1].*bsSur,[0 30])
    hold on
    plot(xlim,[bL]+[.5;.5],'k-','linewidth',1)
    %title(cmbSmp.tp(s))
    subplot(4,numel(selSmp),c)
    hold off
    for p=1:2
        plot([-bsSur:bsSur]',movmean(mean(imageMat(currBs.type==p,:,:),1),5,2),'-','linewidth',2)
        hold on
    end
    ylim([0 max(movmean(mean(imageMat(currBs.type==1,:,:),1),5,2)).*1.1])
    title(cmbSmp.tp(s))    
end
sgtitle([tfOrder{t} ,'-',gtOrder{g}])
saveas(gcf,sprintf('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/Fig6MotifRaw.fig'))

%% Figure S6A
subplot(1,1,1)
bsSur=150;
gtOrder={'aaino80','aaspt6','aaspt16','aachd1','aaspt6spt16'};
tfOrder={'msn2','reb1'};
for t=1:2
    subplot(1,1,1)
    hold off
    currBs=bsAll(contains(lower(bsAll.tf),tfOrder(t))&ismember(bsAll.type,[1,2])&(bsAll.gene>0),:);    
    for g=1:numel(gtOrder)
        selSmp=find(cmbSmp.tp<=90 &ismember(cmbSmp.gt,gtOrder(g))&ismember(cmbSmp.tf,tfOrder(t)))
        
        c=0;
        for s=selSmp'
            imageMat=reshape(cmbProfile(currBs.pos+[-bsSur:bsSur].*currBs.dir,s),height(currBs),2*bsSur+1);
            c=c+1;
            subplot(numel(gtOrder),numel(selSmp),c+(g-1)*numel(selSmp))
            hold off
            for p=1:2
                plot([-bsSur:bsSur]',movmean(mean(imageMat(currBs.type==p,:,:),1),5,2),'-','linewidth',2)
                hold on
            end
            ylim([0 max(movmean(mean(imageMat(currBs.type==1,:,:),1),5,2)).*1.1])
            title(sprintf('%s-%d',cmbSmp.gt{s},cmbSmp.tp(s)))
            axis tight
            ylim([0 max(ylim)])
        end
    end
    sgtitle(tfOrder{t})
    saveas(gcf,sprintf('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/svg/FigS6MotifMean-%s.fig',tfOrder{t}))
    pause(0.5)
end



%% Figure 5D
surLen=100;
figure('Color',[1 1 1],'Renderer','painters','Position',[2127 288 560 361])
subplot(1,1,1)
c=0;
typeCol=brighten(lines(2),0.5);
typeCol(1,:)=[1 1 1].*0.65;
for t=1:numel(tfOrder)
    c=c+1;
    yL=[];
    currBs=bsAll(contains(lower(bsAll.tf),tfOrder(t))&bsAll.type<3.,:);
    for g=1:numel(gtOrder)
        i=find(cmbSmp.tp==0 &ismember(cmbSmp.gt,gtOrder(g))&ismember(cmbSmp.tf,tfOrder(t)))
        tp90=find(cmbSmp.tp==90&ismember(cmbSmp(:,{'tf','gt'}),cmbSmp(i,{'tf','gt'})))
        
        c=c+1;
        
        if numel(tp90)==1
            bsSur0=reshape(cmbProfile(currBs.pos+[-surLen:surLen].*currBs.dir,i),height(currBs),2*surLen+1,[]);
            bsSur90=reshape(cmbProfile(currBs.pos+[-surLen:surLen].*currBs.dir,tp90),height(currBs),2*surLen+1,[]);
            for p=[1,2]
                subplot(numel(tfOrder),numel(gtOrder)+1,(t-1)*(numel(gtOrder)+1)+1)
                plot([-surLen:surLen]',movmean(mean(bsSur0(currBs.type==p,:,:),1),5,2),'-','Color',typeCol(p,:));%brighten([1 1 1].*0.25,(2-p)*0.8))
                hold on
                subplot(numel(tfOrder),numel(gtOrder)+1,c)
                plot([-surLen:surLen]',movmean(mean(bsSur90(currBs.type==p,:,:),1),5,2),'-','Color',typeCol(p,:),'Linewidth',2)
                hold on
                title(strcat(cmbSmp.tf{i},'-',cmbSmp.gt{i},'-90'))
                ylim([0 max(movmean(mean(bsSur90(currBs.type==1,:,:),1),5,2))*1.1])
                yL=[yL, max(movmean(mean(bsSur0(currBs.type==1,:,:),1),5,2))];
                xlabel('distance from motif')
                ylabel(sprintf('mean checSignal (%d,%d)',sum(currBs.type==1),sum(currBs.type==2)))
            end
        end
    end
    subplot(numel(tfOrder),numel(gtOrder)+1,(t-1)*(numel(gtOrder)+1)+1)
    title([tfOrder{t},'-TP 0'])
    ylim([0 max(yL).*1.1])
    xlabel('distance from motif')
    ylabel('mean checSignal')
end
saveas(gcf,sprintf('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/Fig6Motif.fig'))



%% Figure 5E
cScales={'greens','blues','Reds','oranges','greys'}
figure('Color',[1 1 1],'Renderer','painters','Position',[2127 288 560 361])
bsSur=25;
subplot(1,1,1)
gtCol=lines(numel(gtOrder))
for t=1:numel(tfOrder)
    currBs=bsAll(contains(lower(bsAll.tf),tfOrder(t))&ismember(bsAll.type,[1,2])&(bsAll.gene>0),:);
    motScore=squeeze(mean(reshape(cmbProfile(currBs.pos+[-bsSur:-6,6:bsSur],:),height(currBs),2*bsSur-10,width(cmbProfile)),2));
    motProm=mean(motScore(currBs.type==1,:))
    motGB=mean(motScore(currBs.type==2,:))
    subplot(1,2,t)
    for g=1:numel(gtOrder)
        
        selSmp=find(cmbSmp.tp<=90 &ismember(cmbSmp.gt,gtOrder(g))&ismember(cmbSmp.tf,tfOrder(t)))
        cMap=brewermap(numel(selSmp)+2,cScales{g})
        blankVal=-log2(motProm(selSmp(end)))+log2(motGB(selSmp(end)))+log2(motProm(selSmp(1)))-log2(motGB(selSmp(1)))
        scatter(cmbSmp.tp(selSmp),-log2(motProm(selSmp))+log2(motGB(selSmp)),100,gtCol(g,:),'filled','MarkerEdgeColor',[1 1 1].*0.35,'MarkerFaceAlpha',rescale(blankVal,0.25,1,'InputMin',0,'InputMax',2))
        hold on
        plot(cmbSmp.tp(selSmp),-log2(motProm(selSmp))+log2(motGB(selSmp)),'-','Color',gtCol(g,:))
       
    end
            xlabel('time(min)')
        ylabel('ratio (gb vs. prom)')
    %xlim([min([xlim,ylim]) max([xlim,ylim])])
    %ylim(xlim)
    %plot(xlim,xlim,'k-')
     title([tfOrder{t}])
end
saveas(gcf,sprintf('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/Fig6Ratio.fig'))

%% Figure S6B
gtOrder={'aaino80','aaspt6','aaspt16','aachd1','aaspt6spt16'};
figure('Color',[1 1 1],'Renderer','painters','Position',[2127 288 560 361])
bsSur=25;
subplot(1,1,1)
tfOrder={'sok2','yap1','tbp'};
gtCol=lines(numel(gtOrder));
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/promoterIDXvecFimp.mat','promoterIDXvecF')
motProm=mean(cmbProfile(promoterIDXvecF==1,:));
motGB=mean(cmbProfile(promoterIDXvecF==2,:));
for t=1:numel(tfOrder)
    subplot(3,3,t)
    hold off
    clear plotLine
    c=0
    for g=1:numel(gtOrder)        
        selSmp=find(cmbSmp.tp<=90 &ismember(cmbSmp.gt,gtOrder(g))&ismember(cmbSmp.tf,tfOrder(t)))
        if numel(selSmp)>0
            c=c+1
        blankVal=log2(motGB(selSmp(end)))-log2(motProm(selSmp(end)))+log2(motProm(selSmp(1)))-log2(motGB(selSmp(1)))
        scatter(cmbSmp.tp(selSmp),log2(motGB(selSmp))-log2(motProm(selSmp)),100,gtCol(g,:),'filled','MarkerEdgeColor',[1 1 1].*0.35,'MarkerFaceAlpha',rescale(blankVal,0.25,1,'InputMin',0,'InputMax',2))
        hold on
        plotLine(c)=plot(cmbSmp.tp(selSmp),-log2(motProm(selSmp))+log2(motGB(selSmp)),'-','Color',gtCol(g,:),'DisplayName',gtOrder{g})
        end
    end
            xlabel('time(min)')
        ylabel('ratio (gb vs. prom (not motif))')
    %xlim([min([xlim,ylim]) max([xlim,ylim])])
    %ylim(xlim)
    %plot(xlim,xlim,'k-')
     title([tfOrder{t}])
     ylim(mean(ylim)+[-2 2])
     legend(plotLine,'Location','best')
end
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/Fig6PromGBAll.fig')

%% Figure 5D
%define Arou80 targets
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/promoterIDXvecFimp.mat','promoterIDXvecF')
tfTrace=mean(normProfile(:,ismember(sTable.tf,'aro80')),2);
tfTrace(promoterIDXvecF~=1,:)=0;
[pVal,pPos]=findpeaks(movmean(tfTrace,101),'NPeaks',100,'MinPeakDistance',100,'SortStr','descend')
keepIdx=1:find([diff(log2(pVal));-2]<-1,1)
aro80Pos=pPos(keepIdx)
aro80Sig(1,:)=sum(normProfile(acol(aro80Pos+[-100:100]),:));

load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/bsAll.mat');
bsAll.pos=round(bsAll.pos);
bsAll.type=promoterIDXvecF(bsAll.pos);
bsAll.dir=ismember(bsAll.dir,'+')*2-1;
selRows=sTable.spike==1 & contains(sTable.rpt,'x');
[spikeTcs,~,sTable.stc(selRows)]=unique(sTable(selRows,{'tf','gt','rpt'}));
tfOrder=unique(spikeTcs.tf)
gtOrder=unique(spikeTcs.gt)
c=0;
bsSur=25;
subplot(1,1,1)
motCols=[([1 1 1].*.65);0.8500    0.3250    0.0980];
for t=1:numel(tfOrder)
    currBs=bsAll(contains(lower(bsAll.tf),tfOrder(t))&ismember(bsAll.type,[1,2])&min(abs(bsAll.pos-aro80Pos')>200,[],2),:);
    motScore=squeeze(mean(reshape(normProfile(currBs.pos+[-bsSur:-6,6:bsSur],:),height(currBs),2*bsSur-10,width(normProfile)),2));
    motProm=mean(motScore(currBs.type==1,:));
    motGB=mean(motScore(currBs.type==2,:));
    for g=1:numel(gtOrder)
        c=c+1;
        subplot(numel(tfOrder),numel(gtOrder),c)
        hold off
        for tc=find(ismember(spikeTcs.tf,tfOrder(t))&ismember(spikeTcs.gt,gtOrder(g)))'
            selSmp=find(sTable.stc==tc)'
            scatter(sTable.tp(selSmp),log2(motProm(selSmp))-log2(aro80Sig(selSmp)),100,motCols(1,:),'filled')
            hold on
            plot(sTable.tp(selSmp),log2(motProm(selSmp))-log2(aro80Sig(selSmp)),'--','Color',motCols(1,:))
            
            scatter(sTable.tp(selSmp),log2(motGB(selSmp))-log2(aro80Sig(selSmp)),100,motCols(2,:),'filled')
            hold on
            plot(sTable.tp(selSmp),log2(motGB(selSmp))-log2(aro80Sig(selSmp)),'--','Color',motCols(2,:))
        end
        title(sprintf('%s -%s',tfOrder{t},gtOrder{g}))
    end
end
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/Fig6Abs.fig')

%% Figure 5G , Figure S6C
tfOrder={'msn2','reb1'}
gtOrder={'aaspt6','aaspt16','aaspt6spt16'};
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/promoterIDXvecFimp.mat','promoterIDXvecF')
h3=load('dataFScomb.mat')
figure('Color',[1 1 1],'Renderer','painters','Position',[2127 288 560 361])
nucSur=75;
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/bsAll.mat');
bsAll.pos=round(bsAll.pos);
bsAll.type=promoterIDXvecF(bsAll.pos);
bsAll.dir=ismember(bsAll.dir,'+')*2-1;
c=0
clear mycCr haCr
% prom vs GB direct
nucSur=75;
bsSur=25;
zTh=2;
intTp=[0,60,90]
name={'not','bound'}
gtCols=lines(3);
typeStyle={'o','filled'}
nMax=[592,345]
c=0
subplot(1,1,1)
for t=1:numel(tfOrder)
    currBs=bsAll(contains(lower(bsAll.tf),tfOrder(t))&ismember(bsAll.type,[1,2]),:);
    smoothProf=movsum(cmbProfile,2*bsSur+1,1)-movsum(cmbProfile,2*5+1,1);
    smoothProf=smoothProf./(2*(bsSur-5));
    motScore2=smoothProf(currBs.pos,:);
    zSc=nan(size(smoothProf));
    zSc(promoterIDXvecF<3,:)=zscore(smoothProf(promoterIDXvecF<3,:));
    %figure('Color',[1 1 1],'Renderer','painters')
    motZsc=zSc(currBs.pos,:);
    h3.motScore=squeeze(mean(reshape(h3.data(currBs.pos+[-nucSur:nucSur],:),height(currBs),2*nucSur+1,width(h3.data)),2));
    %subplot(1,1,1)
    c=c+1;
    for g=1:numel(gtOrder)
        selSmp=find(ismember(cmbSmp.tf,tfOrder(t))&ismember(cmbSmp.gt,gtOrder(g))&ismember(cmbSmp.tp,intTp))
        clear haScore mycScore mycMean haMean %ncount
        for tp=1:numel(intTp)
            mycScore(:,tp)=mean(h3.motScore(:,ismember(h3.smeta.ab,'myc')&ismember(h3.smeta.tag2,'hht2')&ismember(h3.smeta.gt,gtOrder(g))&h3.smeta.tp==intTp(tp)&h3.smeta.bad==0),2);
            haScore(:,tp)=mean(h3.motScore(:,ismember(h3.smeta.ab,'ha')&ismember(h3.smeta.tag2,'hht2')&ismember(h3.smeta.gt,gtOrder(g))&h3.smeta.tp==intTp(tp)&h3.smeta.bad==0),2);
            for bound=[0,1]
                for gtype=[1,2]
                    selMots=currBs.type==gtype & (motZsc(:,selSmp(end))>=zTh)==bound;
                    mycMean(bound+2*(gtype-1)+1,tp)=median(mycScore(selMots,tp));
                    haMean(bound+2*(gtype-1)+1,tp)=median(haScore(selMots,tp));
                    nCount(bound+2*(gtype-1)+1,g)=sum(selMots);
                end
            end
        end
        %c=c+1;
%         subplot(2,numel(gtOrder),c)
%         hold off
%         subplot(2,numel(gtOrder),c+numel(gtOrder))
%         hold off
        for i=1:height(mycMean)/2
            %subplot(2,numel(gtOrder),c)
            subplot(2,numel(tfOrder),c)
            dotSize=rescale(nCount(4,g),20,200,'InputMin',0,'inputMax',nMax(t))
            scatter(mycMean(i,:),mycMean(i+2,:),dotSize,cell2mat(arrayfun(@(x)brighten(gtCols(g,:),x),[0.8,0.4,.0]','UniformOutput',false)),typeStyle{i})
            hold on
            plot(mycMean(i,:),mycMean(i+2,:),'Color',gtCols(g,:))
            subplot(2,numel(tfOrder),c+numel(tfOrder))
             scatter(haMean(i,:),haMean(i+2,:),dotSize,cell2mat(arrayfun(@(x)brighten(gtCols(g,:),x),[0.8,0.4,.0]','UniformOutput',false)),typeStyle{i})
            hold on
            plot(haMean(i,:),haMean(i+2,:),'Color',gtCols(g,:))           
        end
                subplot(2,numel(tfOrder),c)
        title(sprintf('%s %s',tfOrder{t},'myc'))
        xlim(quantile([xlim,ylim],[0 1]))
        ylim(quantile([xlim,ylim],[0 1])+[0 0.1])
        xlabel('promter')
        ylabel('GB')
        plot(xlim,xlim,'k--')
        subplot(2,numel(tfOrder),c+numel(tfOrder))
              title(sprintf('%s %s',tfOrder{t},'HA'))
        xlim(quantile([xlim,ylim],[0 1]))
        ylim(quantile([xlim,ylim],[0 1]))
        xlabel('promter')
        ylabel('GB')
        plot(xlim,xlim,'k--')
    end
end
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/Fig6-mycHA.fig')

%% Figure 5F
load('dataFS-single.mat', 'nucMeta')
[~,nucIdx]=min(abs(bsAll.pos-nucMeta.pos'),[],2);
bsAll.gene=nucMeta.gene(nucIdx);
tfOrder={'msn2','reb1'};

[typeVec,typeNames]=createTypeVec();
typeVec(promoterIDXvecF==1)=numel(typeNames)+1;
typeNames=[typeNames,'promAll'];
bsAll.type2=typeVec(bsAll.pos);
[~,intType]=ismember({'promAll','5','mid','3'},typeNames);
dotPos=[-0.602497883514248;-0.201539381931564;0.467383031397502;1.06369127860685;1.74244903405045]
expMat=repmat(dotPos,1,numel(intType));
nucMat=repmat(1:numel(intType),numel(dotPos),1);

gtOrder={'aaino80','aaspt6','aaspt16','aachd1','aaspt6spt16'};

bsSur=25;
nBins=5;
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/addData.mat','geneLvl','opnScore')

currPar=eval('geneLvl');
binEdge=quantile(currPar,[0.01,0.99]);
[geneId,~]=discretize(max(min(currPar,binEdge(2)),binEdge(1)),numel(dotPos));
geneId(isnan(currPar))=nan;
bsAll.gGroup(bsAll.gene>0)= geneId(bsAll.gene(bsAll.gene>0));

cLim={[],[]}
c=0
yDist=3.5;
subplot(1,1,1)
for t=1:numel(tfOrder)
    c=0;
    currBs=bsAll(contains(lower(bsAll.tf),tfOrder(t))&ismember(bsAll.type2,intType)&(bsAll.gene>0),:);
    motScore=squeeze(mean(reshape(cmbProfile(currBs.pos+[-bsSur:bsSur],:),height(currBs),2*bsSur+1,width(cmbProfile)),2));
    for g=1:numel(gtOrder)
        selSmp=find(ismember(cmbSmp.gt,gtOrder(g))&ismember(cmbSmp.tf,tfOrder(t))&cmbSmp.tp<=90);
        selMot=log2(motScore(:,selSmp)+1)-mean(log2(motScore(:,selSmp(cmbSmp.tp(selSmp)==0))+1),2);
        for b=1:max(geneId)
            for s=1:numel(intType)
                selRows=currBs.gGroup==b & currBs.type2==intType(s);
                if sum(selRows)>0
                    crColl(b,s)=corr(acol(selMot(selRows,:)),acol(repmat(cmbSmp.tp(selSmp)',sum(selRows),1)));
                    fitColl(b,s,:)=polyfit(acol(repmat(cmbSmp.tp(selSmp)',sum(selRows),1)),acol(selMot(selRows,:)),1);
                else
                    crColl(b,s)=nan;
                    fitColl(b,s,:)=nan;
                end
            end
        end
        c=c+1
        %subplot(numel(tfOrder),numel(gtOrder),)
        subplot(1,numel(tfOrder),t)
        scatter(nucMat(:),expMat(:)-c*yDist,rescale(acol(crColl),40,100,'InputMin',0,'InputMax',1),acol(fitColl(:,:,1).*90),'filled','MarkerEdgeColor',[1 1 1].*.65)      
        hold on
    end
    xlim([-1 numel(intType)]+.5)
    set(gca,'Xtick',1:numel(intType),'Xticklabel',typeNames(intType),'Ytick',fliplr(mean(expMat(:))-yDist*[1:numel(gtOrder)]),'YtickLabel',fliplr(gtOrder))
    title(tfOrder{t})
    ylabel(colorbar(),'log2 change')
    bL=mean(expMat(:))-yDist*[1:numel(gtOrder)]-yDist/2
    hold on
    plot(xlim',bL.*[1;1],'k-')
    ylim([numel(gtOrder),0].*-yDist+mean(expMat(:))-yDist/2)
end
saveas(gcf,sprintf('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/Fig6SumScatter.fig'))
sgtitle('GB motifs')
