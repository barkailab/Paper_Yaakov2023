%% Figure 4A
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/dataFScomb.mat')
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');

clearvars -except smeta data tcs nucMeta intChaps chapNames normProfile GP
segMat=nan(6701,3,height(smeta));
%segChap=nan(6701,3,width(normProfile));

for g=find(abs(diff(GP.gene_infoR64.felixTss(:,2:3),1,2))>700)'
    posi=GP.gene_infoR64.felixTss(g,:);
    posi=GP.chrIdx(posi(1,1))+posi(1,2:3);
    segMat(g,1,:)=mean(data(posi(1)+[0:350].*GP.gene_infoR64.dir(g),:),1);
    segMat(g,3,:)=mean(data(posi(2)+[-350:0].*GP.gene_infoR64.dir(g),:),1);
    if abs(diff(posi))>1050
        segMat(g,2,:)=mean(data(posi(1)+350.*GP.gene_infoR64.dir(g):GP.gene_infoR64.dir(g):posi(2)-350.*GP.gene_infoR64.dir(g),:));
    end
    
%     segChap(g,1,:)=mean(normProfile(posi(1)+[0:350].*GP.gene_infoR64.dir(g),:),1);
%     segChap(g,3,:)=mean(normProfile(posi(2)+[-350:0].*GP.gene_infoR64.dir(g),:),1);
%     if abs(diff(posi))>1050
%         segChap(g,2,:)=mean(normProfile(posi(1)+350.*GP.gene_infoR64.dir(g):GP.gene_infoR64.dir(g):posi(2)-350.*GP.gene_infoR64.dir(g),:));
%     end
end
close all
nBins=5;
clearvars crVal pVal
sSel={1,2,3};
segMat=segMat;
dotPos=[-0.602497883514248;-0.201539381931564;0.467383031397502;1.06369127860685;1.74244903405045]
nSegs= size(sSel,2);
ySeg=5;
allChaps=unique(smeta.gt);
%chapOrder=[14,15,2,11,12,1,16,13,3,4     5     6     7     8     9    10]
chapOrder=[14,11,16,15,12,13,2,1,3,5,8]

allChaps=allChaps(chapOrder)
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/addData.mat','geneLvl')
intParas={'geneLvl'}
for p=1
    currPar=eval(intParas{p});
    binEdge=quantile(currPar,[0.01,0.99]);
    [geneId,~]=discretize(max(min(currPar,binEdge(2)),binEdge(1)),nBins);
    geneId(isnan(currPar))=nan;
    expMat=repmat(dotPos,1,numel(sSel));
    nucMat=repmat([-1 0,1],nBins,1);
   % figure('Color',[1 1 1],'Renderer','painters','Position', [61 346 1692 560])
    subplot(1,1,1)
    subplot(1,5,1)
    hold off
    c=0
    for g=allChaps'
        c=c+1;
        sb=1;
        subplot(1,5,1)
%         [~,currCh]=ismember(extractAfter(g,2),lower(chapNames));
%         if currCh>0
%             binBind=cell2mat(cellfun(@(s)accumarray(geneId(geneId>0),chapEnr(geneId>0,s,currCh),[],@(x)mean(x,'omitnan')),sSel,'UniformOutput',false));
%             scatter(acol(nucMat),expMat(:)-(ySeg)*c,40,acol(binBind),'filled','MarkerEdgeColor',[1 1 1].*0.65)
%             hold on
%         end
        caxis([0 .75])
        xlim([-1 1].*1.5)
        ylim(quantile(expMat(:,1),[0 1])-[numel(allChaps),1]*ySeg)
        plot(xlim,[-1.93].*[1 1]-(ySeg)*c,'-','Color',[1 1 1].*0.65)
        set(gca,'Ytick',fliplr(mean(expMat(:,1))-[1:numel(allChaps)]*ySeg),'Yticklabel',flipud(extractAfter(allChaps,2)),'Xtick',-1:1,'Xticklabel',{'5p','mid','3p'})
        colormap(gca,brewermap(128,'YlOrRd'))
        clear crVal pVal crType
        ylabel('chapEnr')
        for t2={'hht2','h2b'}
            for ab={'myc','ha'}
                sb=sb+1;
                intSmp=find(ismember(smeta.tid,find(ismember(tcs.ab,ab)&ismember(tcs.tag2,t2)&ismember(tcs.gt,g))) & smeta.bad==0 &smeta.tp<=90);
                if numel(intSmp)>0
                    dSeg=log2(segMat(:,:,intSmp)+.1)-log2(mean(segMat(:,:,intSmp)+.1,3));
                    timeVec=repmat(permute(smeta.tp(intSmp),[3,2,1]),6701,width(dSeg));
                    for b=1:max(geneId)
                        sc=0;
                        for s=sSel;
                            sc=sc+1;
                            tempY=acol(dSeg(geneId==b,s{1},:));
                            tempX=acol(timeVec(geneId==b,s{1},:));
                            tempX=tempX(~isnan(tempY));
                            tempY=tempY(~isnan(tempY));
                            [crVal{c}(b,sc),pVal{c}(b,sc)]=corr(tempX,tempY);
                            pFit=polyfit(tempX,tempY,1);
                            fitVal{c}(b,sc)=pFit(1);
                        end
                    end
                    crType(c,:)=[t2,ab,g];
                    subplot(1,5,sb)
                    scatter(acol(nucMat),expMat(:)-(ySeg)*c,rescale(min(-log10(pVal{c}(:)),10).*abs(crVal{c}(:)).*2,20,100,'InputMax',10),fitVal{c}(:)*90,'filled','MarkerEdgeColor',[1 1 1].*0.4)
                    hold on
                    plot(xlim,[-1.93].*[1 1]-(ySeg)*c,'-','Color',[1 1 1].*0.65)
                    caxis([-1 1].*1.2)
                    xlim([-1 1].*1.5)
                    yticks([])
                    ylim(quantile(expMat(:,1),[0 1])-[numel(allChaps),1]*ySeg)
                    colormap(gca,flipud(brighten(brewermap(128,'RdBu'),0.5)));
                    ylabel([t2{1} '-' ab{1}])
                    xlabel('nucPosition')
                    set(gca,'Xtick',[-1:1],'Xticklabel',{'5p','mid','3p'})
                    %ylabel(colorbar(),'RapEffect-log2 in 90min')
                end
                if sb==2
                    set(gca,'Ytick',fliplr(mean(expMat(:,1))-[1:numel(allChaps)]*ySeg),'Yticklabel',flipud(extractAfter(allChaps,2)),'Xtick',-1:1,'Xticklabel',{'5p','mid','3p'})
                end
            end
        end
    end
    sgtitle(sprintf('%s - 5p,mid,3p',intParas{p}))
end

%% Figure  4B, S5C
clearvars -except smeta data tcs nucMeta intChaps chapNames normProfile  GP
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');

promLen=250;
fiveLen=350;
midLen=300;
threeLen=350;
termLen=100;
plotMat=nan(6701,promLen+1+fiveLen+midLen+threeLen+1+termLen,width(data));
%plotChap=nan(6701,promLen+1+fiveLen+midLen+threeLen+1+termLen,width(normProfile));

for g=find(abs(diff(GP.gene_infoR64.felixTss(:,2:3),1,2))>1500)'
    posi=GP.gene_infoR64.felixTss(g,:);
    posi=GP.chrIdx(posi(1,1))+posi(1,2:3);
    promBases=posi(1)+[-promLen:fiveLen].*GP.gene_infoR64.dir(g);
    termBases=posi(2)+[-threeLen:termLen].*GP.gene_infoR64.dir(g);
    
    plotMat(g,1:promLen+1+fiveLen,:)=data(promBases,:);
    plotMat(g,end-termLen-threeLen:end,:)=data(termBases,:);
    
%     plotChap(g,1:promLen+1+fiveLen,:)=normProfile(promBases,:);
%     plotChap(g,end-termLen-threeLen:end,:)=normProfile(termBases,:);
    
    if abs(diff(posi))>1050
        geneLen=abs(diff(posi))-fiveLen-threeLen;
        
        midBases=posi(1)+[fiveLen+1:fiveLen+geneLen].*GP.gene_infoR64.dir(g);
        
        tempMat=movmean(data(midBases,:),10);
        keepBases=round(linspace(1,height(tempMat),midLen));
        plotMat(g,promLen+fiveLen+2:promLen+fiveLen+midLen+1,:)=tempMat(keepBases,:);
        
%         tempMat=movmean(normProfile(midBases,:),10);
%         plotChap(g,promLen+fiveLen+2:promLen+fiveLen+midLen+1,:)=tempMat(keepBases,:);
        
        
    end
end

load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/addData.mat','geneLvl')
binEdge=quantile(geneLvl,[0.01,0.99]);
[geneId,~]=discretize(max(min(geneLvl,binEdge(2)),binEdge(1)),5);
geneId(isnan(geneLvl))=nan;
windowSize=51;
selChapColl={{'aaspt6','aaspt16','aaspt6spt16'};{'aaspt6','aa3spt6'};{ 'aaspt16','aa3spt16'};{'aaspt6spt16','aa5'}}
for ab={'myc','ha'}
    figure('Color',[1 1 1],'Renderer','painters','Position',[2127 288 560 361])
    subplot(1,1,1)
    for gg=1:height(selChapColl)
        selChaps=selChapColl{gg};
        for g=selChaps
            sel0=find(smeta.bad==0&ismember(smeta.gt,g)&ismember(smeta.ab,ab)&ismember(smeta.tag2,'hht2')&smeta.tp==0);
            sel90=find(smeta.bad==0&ismember(smeta.gt,g)&ismember(smeta.ab,ab)&ismember(smeta.tag2,'hht2')&smeta.tp==90);
            
            dMat=log2(movmean(mean(plotMat(geneId==5,:,sel90),[3,1],'omitnan'),windowSize)+0.25)-log2(movmean(mean(plotMat(geneId==5,:,sel0),[3,1],'omitnan'),51)+0.25);
            subplot(4,1,gg)
            plot([-promLen:fiveLen+midLen+threeLen+1+termLen],dMat,'Linewidth',2,'DisplayName',g{1})
            hold on
            xlim([-promLen,fiveLen+midLen+threeLen+1+termLen])
            ylabel('\Delta Log2')
            title('Chap Effect')
            %    [~,currCh]=ismember(extractAfter(g,2),lower(chapNames));
            %     subplot(2,1,2)
            %     plot([-promLen:fiveLen+midLen+threeLen+1+termLen],rescale(movmean(mean(plotChap(geneId==5,:,currCh),[3,1],'omitnan'),windowSize),0,1),'Linewidth',2,'DisplayName',g{1})
            %     hold on
            %     xlim([-promLen,fiveLen+midLen+threeLen+1+termLen])
            %     ylabel('relative  Bind')
            %     title('Chap loc')
        end
        legend()
    end
    pause(1)
    saveas(gcf,sprintf('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/Fig5Meta%s.fig',ab{1}))
end
%% Figure 4C
clearvars -except smeta data tcs nucMeta intChaps chapNames normProfile GP
segMat=nan(6701,3,height(smeta));
%segChap=nan(6701,3,width(normProfile));

for g=find(abs(diff(GP.gene_infoR64.felixTss(:,2:3),1,2))>700)'
    posi=GP.gene_infoR64.felixTss(g,:);
    posi=GP.chrIdx(posi(1,1))+posi(1,2:3);
    segMat(g,1,:)=mean(data(posi(1)+[0:350].*GP.gene_infoR64.dir(g),:),1);
    segMat(g,3,:)=mean(data(posi(2)+[-350:0].*GP.gene_infoR64.dir(g),:),1);
    if abs(diff(posi))>1050
        segMat(g,2,:)=mean(data(posi(1)+350.*GP.gene_infoR64.dir(g):GP.gene_infoR64.dir(g):posi(2)-350.*GP.gene_infoR64.dir(g),:));
    end
    
%     segChap(g,1,:)=mean(normProfile(posi(1)+[0:350].*GP.gene_infoR64.dir(g),:),1);
%     segChap(g,3,:)=mean(normProfile(posi(2)+[-350:0].*GP.gene_infoR64.dir(g),:),1);
%     if abs(diff(posi))>1050
%         segChap(g,2,:)=mean(normProfile(posi(1)+350.*GP.gene_infoR64.dir(g):GP.gene_infoR64.dir(g):posi(2)-350.*GP.gene_infoR64.dir(g),:));
%     end
end
intChaps={'aaspt6','aaspt6spt16','aaspt16'}
selSmp=find(ismember(smeta.ab,'myc')&ismember(smeta.gt,intChaps)&smeta.bad==0&ismember(smeta.tp,[0 20 40 60 90]));
pcaData=reshape(segMat(:,:,selSmp),[],numel(selSmp));
selSeg=all(~isnan(pcaData),2);
[c,s,l,e,k]=pca(pcaData(selSeg,:)','NumComponents',2)
cScales={'greens','blues','Reds','oranges','greys'}
hold off
for ch=1:numel(intChaps)
    plotSmp=ismember(smeta.gt(selSmp),intChaps(ch))
    cMap=brewermap(7,cScales{ch})
    meanS=cell2mat(arrayfun(@(tp)mean(s(plotSmp&smeta.tp(selSmp)==tp,:),1),[0 20 40 60 90]','UniformOutPut',false));
    stdS=cell2mat(arrayfun(@(tp)std(s(smeta.tp(selSmp)==tp&plotSmp,:),[],1)./sqrt(sum(smeta.tp(selSmp)==tp&plotSmp)),[0 20 40 60 90]','UniformOutPut',false));
    
    scatter(meanS(:,1),meanS(:,2),200,cMap(2:6,:),'filled','DisplayName',intChaps{ch})
    hold on
    errorbar(meanS(:,1),meanS(:,2),stdS(:,1),stdS(:,1),stdS(:,2),stdS(:,2),'Linestyle','none','Capsize',0,'Color',[1 1 1].*0.85)
end
axis tight
sgtitle('PCA analysis')
xlabel(sprintf('PC-1 %.2f',k(1)))
ylabel(sprintf('PC-2 %.2f',k(2)))
title('means+sterror')
saveas(gcf,sprintf('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/Fig5PCA.fig'))
