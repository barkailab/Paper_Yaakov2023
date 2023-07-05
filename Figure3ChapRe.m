load('chapReloc.mat')

%% Figure 3A (and S5A)
promLen=250;
fiveLen=350;
midLen=300;
threeLen=350;
termLen=100;
plotChap=nan(6701,promLen+1+fiveLen+midLen+threeLen+1+termLen,width(normProfile));

for g=find(abs(diff(GP.gene_infoR64.felixTss(:,2:3),1,2))>1500)'
    posi=GP.gene_infoR64.felixTss(g,:);
    posi=GP.chrIdx(posi(1,1))+posi(1,2:3);
    promBases=posi(1)+[-promLen:fiveLen].*GP.gene_infoR64.dir(g);
    termBases=posi(2)+[-threeLen:termLen].*GP.gene_infoR64.dir(g);
    
    plotChap(g,1:promLen+1+fiveLen,:)=normProfile(promBases,:);
    plotChap(g,end-termLen-threeLen:end,:)=normProfile(termBases,:);
    
    if abs(diff(posi))>1050
        geneLen=abs(diff(posi))-fiveLen-threeLen;
        
        midBases=posi(1)+[fiveLen+1:fiveLen+geneLen].*GP.gene_infoR64.dir(g);
        
        
        tempMat=movmean(normProfile(midBases,:),10);
         keepBases=round(linspace(1,height(tempMat),midLen));
        plotChap(g,promLen+fiveLen+2:promLen+fiveLen+midLen+1,:)=tempMat(keepBases,:);
        
        
    end
end

atf=unique(tcs.tf)
agt={'aaino80','aaspt6','aaspt16','aachd1'}';unique(tcs.gt)
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/addData.mat','geneLvl','opnScore')
[~,geneIdx]=sort(-geneLvl)

atf=unique(tcs.tf)
agt={'aaino80','aaspt6','aaspt16','aachd1'}';unique(tcs.gt)
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/addData.mat','geneLvl','opnScore')
[~,geneIdx]=sort(-geneLvl)

figure('Color',[1 1 1],'Renderer','painters','Position',[1 41 1920 962])
lWidth=1;
geneSel={geneIdx(1:300),geneIdx(301:end)};
geneName={'highGenes','rest'}


aaOrder=[2,4,3,1]
tfOrder=[1,2,4,5,3]
for g=1%:numel(geneSel)
    selGenes=geneSel{g};
    c=0;
    for aa=agt(aaOrder)'
        for tf=atf(tfOrder)'
            c=c+1;
            subplot(numel(aaOrder),numel(tfOrder),c)
            hold off
            for tp=[0,60]
                selSmp=find(ismember(sTable.tf,tf)&ismember(sTable.gt,aa)&sTable.tp==tp&~contains(sTable.rpt,{'$','*'}))'
                meanVal=mean(movmean(mean(plotChap(selGenes,:,selSmp),[1],'omitnan'),51),3);
                stdVal=std(movmean(mean(plotChap(selGenes,:,selSmp),[1],'omitnan'),51),[],3)./sqrt(numel(selSmp));
                if tp==0
                    %                 for s=selSmp
                    %                     plot([-promLen:fiveLen+midLen+threeLen+termLen+1],movmean(mean(plotChap(selGenes,:,s),[1,3],'omitnan'),51),'DisplayName',sprintf('%d',tp),'Linewidth',.25,'Color',[1 1 1].*0.25)
                    %                     hold on
                    %                 end
                    plot([-promLen:fiveLen+midLen+threeLen+termLen+1],meanVal,'DisplayName',sprintf('%d',tp),'Linewidth',lWidth,'Color',[1 1 1].*0.25)
                    hold on
                    fill([[-promLen:fiveLen+midLen+threeLen+termLen+1],fliplr([-promLen:fiveLen+midLen+threeLen+termLen+1])],[meanVal-stdVal,fliplr(meanVal+stdVal)],[1 1 1].*.25,'Linestyle','none','FaceAlpha',0.15)
                else
                    %                 for s=selSmp
                    %                     plot([-promLen:fiveLen+midLen+threeLen+termLen+1],movmean(mean(plotChap(selGenes,:,s),[1,3],'omitnan'),51),'DisplayName',sprintf('%d',tp),'Linewidth',.25,'Color',[1 1 1].*0.65)
                    %                     hold on
                    %                 end
                    plot([-promLen:fiveLen+midLen+threeLen+termLen+1],meanVal,'-','DisplayName',sprintf('%d',tp),'Linewidth',lWidth,'Color',[1 1 1].*0.65)
                    hold on
                    fill([[-promLen:fiveLen+midLen+threeLen+termLen+1],fliplr([-promLen:fiveLen+midLen+threeLen+termLen+1])],[meanVal-stdVal,fliplr(meanVal+stdVal)],[1 1 1].*0.65,'Linestyle','none','FaceAlpha',0.15)
                end
                hold on
                
            end
            %legend()
            title([tf{1},'-',aa{1},geneName{1}])
            axis tight
            xlabel('postiion relative to TSS (bp)')
        end
    end
    pause(1)
    saveas(gcf,sprintf('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/svg/Fig4Meta-%s.fig',geneName{g}))
end

%% Figure 3B

clearvars -except smeta data tcs nucMeta intChaps chapNames normProfile GP sTable

segChap=nan(6701,5,width(normProfile));

for g=find(abs(diff(GP.gene_infoR64.felixTss(:,2:3),1,2))>700)'
    posi=GP.gene_infoR64.felixTss(g,:);
    posi=GP.chrIdx(posi(1,1))+posi(1,2:3);
%     segMat(g,1,:)=mean(data(posi(1)+[0:350].*GP.gene_infoR64.dir(g),:),1);
%     segMat(g,3,:)=mean(data(posi(2)+[-350:0].*GP.gene_infoR64.dir(g),:),1);
%     if abs(diff(posi))>1050
%         segMat(g,2,:)=mean(data(posi(1)+350.*GP.gene_infoR64.dir(g):GP.gene_infoR64.dir(g):posi(2)-350.*GP.gene_infoR64.dir(g),:));
%     end
    segChap(g,1,:)=mean(normProfile(posi(1)+[-250:0].*GP.gene_infoR64.dir(g),:),1);
    segChap(g,2,:)=mean(normProfile(posi(1)+[0:350].*GP.gene_infoR64.dir(g),:),1);
    segChap(g,4,:)=mean(normProfile(posi(2)+[-350:0].*GP.gene_infoR64.dir(g),:),1);
     segChap(g,5,:)=mean(normProfile(posi(2)+[1:100].*GP.gene_infoR64.dir(g),:),1);
    if abs(diff(posi))>1050
        segChap(g,3,:)=mean(normProfile(posi(1)+350.*GP.gene_infoR64.dir(g):GP.gene_infoR64.dir(g):posi(2)-350.*GP.gene_infoR64.dir(g),:));
    end    
end


agt=unique(tcs.gt)
atf=unique(tcs.tf)

gtOrder=fliplr([2,4,3,1])
tfOrder=[1,2,5,4,3]

clearvars crVal pVal
sSel={1,2,3,4,5};
segMat=segChap;
dotPos=[-0.602497883514248;-0.201539381931564;0.467383031397502;1.06369127860685;1.74244903405045]
%dotPos=[-0.602497883514248;-0.201539381931564]
nDots=range(dotPos).*2.5;
nBins=numel(dotPos)  %how to bin the genes
%close all
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/addData.mat','geneLvl','opnScore')
intParas={'geneLvl'}
for p=1%:numel(intParas)
    %dSegKeep{1}=nan(nBins,width(segMat),width(data));
    currPar=eval(intParas{p});
    binEdge=quantile(currPar,[0.01,0.99]);
    [geneId,~]=discretize(max(min(currPar,binEdge(2)),binEdge(1)),5);
    geneId(isnan(currPar))=nan;
    %geneId=1+(currPar>11);
    %geneId(isnan(currPar))=nan;
    expMat=repmat(dotPos,1,numel(sSel));
    nucMat=repmat(1:numel(sSel),numel(dotPos),1);
    %figure('Color',[1 1 1],'Renderer','painters','Position',[2127 288 560 361])
    subplot(1,1,1)
    c=0
    for g=agt(gtOrder)'
        sb=0
        c=c+1;
        clear crVal pVal crType fVal
        for tf=atf(tfOrder)'
            sb=sb+1;
            intSmp=find(ismember(sTable.gt,g)&ismember(sTable.tf,tf)&sTable.tp<=60&~contains(sTable.rpt,{'*','$'}));
            if numel(intSmp)>0
                %dSeg=log2(segMat(:,:,intSmp)+1)-log2(mean(segMat(:,:,intSmp)+1,3));
                %dSeg=log2(segMat(:,:,intSmp)+.1)-log2(mean(segMat(:,:,intSmp(sTable.tp(intSmp)==0))+.1,3));
                timeVec=repmat(permute(sTable.tp(intSmp),[3,2,1]),6701,width(segMat));
                clear dSeg
                for b=1:max(geneId)
                    dSeg(b,:,:)=log2(mean(segMat(geneId==b,:,intSmp),1,'omitnan')+0.1)-mean(log2(mean(segMat(geneId==b,:,intSmp),1,'omitnan')+0.1),3)
                end
                timeVec=repmat(permute(sTable.tp(intSmp),[3,2,1]),max(geneId),width(dSeg));
                for b=1:max(geneId)
                    %dSegKeep{1}(b,:,intSmp)=mean(dSeg(geneId==b,:,:),1,'omitnan');
                    sc=0;
                    for s=sSel
                        sc=sc+1;
                        %[crVal{c}(b,sc),pVal{c}(b,sc)]=corr(acol(dSeg(geneId==b,s{1},:)),acol(timeVec(geneId==b,s{1},:)),'rows','complete');
                        [crVal{c}(b,sc),pVal{c}(b,sc)]=corr(acol(dSeg(b,s{1},:)),acol(timeVec(b,s{1},:)),'rows','complete');
                        %selGenes=all(~isnan(dSeg(:,s{1},:)),3)&geneId==b;
                        %pFit=polyfit(acol(timeVec(selGenes,s{1},:)),acol(dSeg(selGenes,s{1},:)),1);
                       % pFit=polyfit(acol(mean(timeVec(selGenes,s{1},:))),acol(median(dSeg(selGenes,s{1},:))),1);
                        pFit=polyfit(acol(timeVec(b,s{1},:)),acol(dSeg(b,s{1},:)),1);
                        fVal{c}(b,sc)=90*pFit(1);
                    end
                    
                end
                %crType(c,:)=[t2,ab,g];
                subplot(5,numel(tfOrder),sb+[1:4].*numel(tfOrder))
                %scatter(acol(nucMat)+nDots*c,expMat(:),rescale(min(-log10(pVal{c}(:)),10).*abs(crVal{c}(:)).*2,20,100,'InputMax',10),fVal{c}(:),'filled')
                if numel(intSmp)>2
                    sizeMat=rescale(min(-log10(pVal{c}(:)),10).*abs(crVal{c}(:)).*2,20,100,'InputMax',10)%50%;
                else
                    sizeMat=50;
                end
                scatter(acol(nucMat),expMat(:)+nDots*c,sizeMat,fVal{c}(:),'filled','MarkerEdgeColor',[1 1 1].*0.65)
                hold on
                caxis([-1 1].*0.8)
                xlim(quantile(nucMat,[0 1],'all')'+[-.5 .5])
                ylim(quantile(expMat(:),[0,1])+[-.1 .1]*nDots+[1 numel(gtOrder)]*nDots)
                yticks([])
                if sb==1
                    set(gca,'ytick',mean(dotPos)+[1:numel(gtOrder)]*nDots,'Yticklabel',agt(gtOrder))
                end
                %xlim([4.5 nDots*numel(intChaps)+3.5]
                xticks(1:size(segMat,2))
                xticklabels({'prom','5p','mid','3p','term'})
                
                %ylabel(colorbar(),'RapEffect')
                colormap(gca,flipud(brighten(brewermap(128,'RdBu'),0.5)));
                
            end
        end
    end
    
    sb=0
    c=c+1;
    clear crVal pVal crType fVal    
    for tf=atf(tfOrder)'
        sb=sb+1;
        selSmp=find(ismember(sTable.tf,tf)&sTable.tp==0&~contains(sTable.rpt,{'*','$'}));
        meanLvl(:,:,sb)=mean(segMat(:,:,selSmp),3)
    end
    relLvl=log2(meanLvl+1)-log2(mean(meanLvl,3)+1);
    
    sb=0
    for tf=atf(tfOrder)'
        valMat=[];
        sb=sb+1;
        for b=1:nBins
            valMat(b,:)=mean(relLvl(geneId==b,:,sb),1,'omitnan');
        end
        subplot(5,numel(tfOrder),sb)
        scatter(acol(nucMat),expMat(:),50,valMat(:),'filled','MarkerEdgeColor',[1 1 1].*0.65)
        ylim(quantile(expMat(:),[0,1])+[-.1 .1]*nDots)
        xlim([.5 5.5])
        title(sprintf('%s',tf{1}))
    end
    sgtitle(sprintf('%s - geneQuantiles - chap Zscore',intParas{p}))
end
saveas(gcf,sprintf('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/svg/Fig4SumMeansBins.fig'))


