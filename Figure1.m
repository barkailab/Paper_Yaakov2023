%%load data - due to size ocnstrains on the github 
clear all
load('dataFS-single.mat')
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/promoterIDXvecFimp.mat','promoterIDXvecF')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/addData.mat','geneLvl','opnScore')

%% Figure 1B
r=1;
intGT=unique(tcs.gt)%{'aaspt6','aacac1','aaino80','aanofrb'}
chapOrder=([5,11,1,12,16,4,3,13,14,6,15,7,9,8,2,10]) %optimalleaforder(chapLink,chapDist)
intGT=intGT(chapOrder)
intGT={'aaspt6','aaino80','aartt106','aavps75','aahir1'}
cols=min(numel(intGT),5);
midTss=GP.chrIdx(min(GP.gene_infoR64.felixTss(:,1),17))+mean(GP.gene_infoR64.felixTss(:,2:3),2);

intGene=GP.groups{23}{2}{45}(7)%[279,1867]%
regSize=2000;
binSize=20;
intRegions=arrayfun(@(x)GP.chrIdx(GP.gene_infoR64.felixTss(x,1))+GP.gene_infoR64.felixTss(x,2)+[-regSize:regSize-1],[intGene],'UNiformoutput',false);
region=intRegions{1};
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/addData.mat','geneLvl','opnScore')
orfLabel=bwlabel(ismember([1:height(normProfile)]',region)&promoterIDXvecF==2);
orfCent=regionprops('Table',orfLabel,'centroid')
[~,orfIdx]=min(abs(orfCent.Centroid(:,2)-midTss'),[],2);
expLabel=zeros(height(normProfile),1);
expLabel(orfLabel>0)=max(geneLvl(orfIdx(orfLabel(orfLabel>0))),5);

%close all
subplot(1,1,1)
for i=1:numel(intGT)%:numel(intRegions)
    c=mod(i-1,5)+1;
    [~,intTCs]=ismember(table({'myc','ha','myc','ha'}',{'hht2','hht2','h2b','h2b'}',repmat(intGT(i),4,1),'VariableNames',{'ab','tag2','gt'}),tcs(:,{'ab','tag2','gt'}));

    
    
    intChap=find(contains(chapNames,extractAfter(unique(tcs.gt(intTCs(intTCs>0))),'aa'),'ignoreCase',true))
    subplot(numel(intTCs)+2,cols,c)
    if numel(intChap)==1
        plot(region-GP.chrIdx(GP.gene_infoR64.position(intGene(r),1)),movmean(normProfile(region,intChap),21))
        title(chapNames{intChap})
    end
    axis tight
    colorbar()
    i
    for tc = intTCs'
        c=c+cols*4;
        if tc>0
            selSmp=ismember(smeta(:,{'ab','tag2','gt'}),tcs(tc,{'ab','tag2','gt'}))&smeta.bad==0 & smeta.tp<=90;
            imMat=movmean(cell2mat(arrayfun(@(x)mean(data(region,selSmp&smeta.tp==x),2),unique(smeta.tp(selSmp))','uniformoutput',false)),1,1);
            imMat=movmean(squeeze(mean(reshape(imMat,binSize,regSize*2/binSize,width(imMat)),1)),7)            
            subplot(4*(numel(intTCs)+2),cols,c+[0 1 2]*cols)
            imagesc(imMat','XData',quantile(region,[0 1])-GP.chrIdx(GP.gene_infoR64.position(intGene(r),1)))
            %caxis([-1.5 1.5])
            title(strjoin(tcs{tc,:}))
            colorbar()
            if strcmp(tcs.tag2(tc),'h2b')
                colormap(gca,brighten(brewermap(128,'purples'),0.5));
            else
                colormap(gca,brighten(brewermap(128,'blues'),0.5));
            end
            %
            xticks([])
        end
    end
    c=c+4*cols;
    subplot(4*(numel(intTCs)+2),cols,c+[0 1 2]*cols)
    imagesc(expLabel(region)','XData',quantile(region,[0 1])-GP.chrIdx(GP.gene_infoR64.position(intGene(r),1)))
    title('promoter vs.genebody')
    colormap(gca,flipud(gray()))
    caxis([4 14])
    xticks(orfCent.Centroid(:,2)-GP.chrIdx(GP.gene_infoR64.position(intGene(r),1)))
    xticklabels(GP.gene_infoR64.nameNew(orfIdx))
    %set(gca,'Position', [0.1300 0.1954 0.7372 0.0126])
    colorbar()
end
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/Fig1Raw.fig')

%% Figure 1D, S2E
clearvars -except smeta data tcs nucMeta intChaps chapNames normProfile  GP
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');

promLen=250;
fiveLen=350;
midLen=300;
threeLen=350;
termLen=100;
plotMat=nan(6701,promLen+1+fiveLen+midLen+threeLen+1+termLen,width(data));
plotChap=nan(6701,promLen+1+fiveLen+midLen+threeLen+1+termLen,width(normProfile));

for g=find(abs(diff(GP.gene_infoR64.felixTss(:,2:3),1,2))>700)' %was 1500
    posi=GP.gene_infoR64.felixTss(g,:);
    posi=GP.chrIdx(posi(1,1))+posi(1,2:3);
    promBases=posi(1)+[-promLen:fiveLen].*GP.gene_infoR64.dir(g);
    termBases=posi(2)+[-threeLen:termLen].*GP.gene_infoR64.dir(g);
    
    plotMat(g,1:promLen+1+fiveLen,:)=data(promBases,:);
    plotMat(g,end-termLen-threeLen:end,:)=data(termBases,:);
    
    plotChap(g,1:promLen+1+fiveLen,:)=normProfile(promBases,:);
    plotChap(g,end-termLen-threeLen:end,:)=normProfile(termBases,:);
    
    if abs(diff(posi))>1050
        geneLen=abs(diff(posi))-fiveLen-threeLen;
        
        midBases=posi(1)+[fiveLen+1:fiveLen+geneLen].*GP.gene_infoR64.dir(g);
        
        tempMat=movmean(data(midBases,:),10);
        keepBases=round(linspace(1,height(tempMat),midLen));
        plotMat(g,promLen+fiveLen+2:promLen+fiveLen+midLen+1,:)=tempMat(keepBases,:);    
        tempMat=movmean(normProfile(midBases,:),10);
        plotChap(g,promLen+fiveLen+2:promLen+fiveLen+midLen+1,:)=tempMat(keepBases,:);
    end
end
metaProfs=squeeze(mean(plotMat,1,'omitnan'));
chapProfs=squeeze(mean(plotChap,1,'omitnan'));
nGenes=sum(any(~isnan(plotMat),[2,3]));
clearvars plotMat plotChap
allChaps=unique(tcs.gt);
regBoarder=cumsum([1,promLen,fiveLen,midLen,threeLen+1,termLen])
chapOrder=[5,11,1,15,6,14,13,3,9,2,16,4,12,7,8,10] %optimalleaforder(chapLink,chapDist)
%chapOrder=([5,11,1,12,16,4,3,13,14,6,15,7,9,8,2,10]) %optimalleaforder(chapLink,chapDist)
allChaps=allChaps(chapOrder);
clear chapOrder

endTp=ones(size(allChaps))*90;
endTp(ismember(allChaps,'aanofrb'))=120;


%figure('Color',[1 1 1],'Renderer','painters','Position',[1 41 1920 962])
subplot(1,1,1)
hold off
c=0
for g=allChaps'
    c=c+1;
    subplot(16,1,c)
    hold off
    selSmp=ismember(smeta.gt,g)&ismember(smeta.ab,'myc')&ismember(smeta.tag2,'hht2')&smeta.bad==0&smeta.tp<=endTp(ismember(allChaps,g));
    imMat=cell2mat(arrayfun(@(x)mean(metaProfs(:,selSmp&smeta.tp==x),2),unique(smeta.tp(selSmp))','uniformoutput',false));
    imMat=log2(imMat+1);
    imMat=imMat-imMat(:,1);
    imagesc(imMat',[-1 1].*.55)
    %colormap(gca,brighten(gy2bl,0.5));
    colormap(gca,flipud(brighten(brewermap(128,'RdBu'),0.5)));
    ylabel(extractAfter(g{1},'aa'))
    xticks(regBoarder)
    xticklabels([])
    hold on
    intChap=find(contains(chapNames,extractAfter(g,'aa'),'ignoreCase',true))
    plot(movmean(rescale(-chapProfs(:,intChap),min(ylim),max(ylim)),21),'-','Color',[1 1 1].*0.65)
end
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/svg/Fig1Meta.fig')

types={'hht2','myc';'hht2','ha';'h2b','myc';'h2b','ha'}
for t=1:height(types)
    subplot(1,1,1)
    c=0
    for g=allChaps'
        c=c+1;
        subplot(16,1,c)
        hold off
        selSmp=ismember(smeta.gt,g)&ismember(smeta.ab,types(t,2))&ismember(smeta.tag2,types(t,1))&smeta.bad==0&smeta.tp<=endTp(ismember(allChaps,g));
        if sum(selSmp)>0
        imMat=cell2mat(arrayfun(@(x)mean(metaProfs(:,selSmp&smeta.tp==x),2),unique(smeta.tp(selSmp))','uniformoutput',false));
        imMat=log2(imMat+1);
        imMat=imMat-imMat(:,1);
        imagesc(imMat',[-1 1].*.55)
        end
        %colormap(gca,brighten(gy2bl,0.5));
        colormap(gca,flipud(brighten(brewermap(128,'RdBu'),0.5)));
        ylabel(extractAfter(g{1},'aa'))
        xticks(regBoarder)
        xticklabels([])
        hold on
        intChap=find(contains(chapNames,extractAfter(g,'aa'),'ignoreCase',true))
        plot(movmean(rescale(-chapProfs(:,intChap),min(ylim),max(ylim)),21),'-','Color',[1 1 1].*0.65)
        
    end
    sgtitle(sprintf('%s-%s (%d)',types{t,1},types{t,2},nGenes))
    saveas(gcf,sprintf('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/svg/Fig1SMeta%s%s.fig',types{t,1},types{t,2}))
end





%% Figure 1C, S2D
clearvars -except smeta data tcs nucMeta chapNames normProfile typeVec typeNames 
load('typeVec.mat')
newTypes=[1,2,3,4,5,6,6,7,8,9,10,11,12,13,14];
typeVec(typeVec>0)=newTypes(typeVec(typeVec>0));
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/promoterIDXvecFimp.mat','promoterIDXvecF')
typeNames=typeNames([1:6,8:end])
typeVec(promoterIDXvecF==1&ismember(typeVec,[0,7,8]))=8;
typeVec(typeVec==7)=0;
for t=1:max(typeVec)
    tempId=bwlabel(typeVec==t);
    typeColl{t}=nan(max(tempId),width(data));
    chapColl{t}=nan(max(tempId),width(normProfile));
    
    for i=1:max(tempId)
        if sum(tempId==i)>100
            typeColl{t}(i,:)=mean(data(tempId==i,:),1);
            chapColl{t}(i,:)=mean(normProfile(tempId==i,:),1);
        else
            typeColl{t}(i,:)=mean(data([],:));
            chapColl{t}(i,:)=mean(normProfile([],:),1);
        end
        
    end
end

for t=1:max(typeVec)
    typeScore(t,:)=mean(data(typeVec==t,:),'omitnan');
    typeChap(t,:)=mean(normProfile(typeVec==t,:),'omitnan');
end

clear crAll pAll pvAll tAll

allChaps=unique(tcs.gt)
endTp=ones(size(allChaps))*90;
endTp(ismember(allChaps,'aanofrb'))=120;

c=0
for t2={'hht2','h2b'}
    for ab={'myc','ha'}
        c=c+1;
        cc=0;
        for g=allChaps'
            cc=cc+1;
            intSmp=find(ismember(smeta.ab,ab)&ismember(smeta.tag2,t2)&ismember(smeta.gt,g)& smeta.bad==0 &smeta.tp<=endTp(ismember(allChaps,g)));
            %intSmp=find(ismember(smeta.ab,ab)&ismember(smeta.tag2,t2)&ismember(smeta.gt,g)& smeta.bad==0 &ismember(smeta.tp,[60,90]));
            conSmp=find(ismember(smeta.ab,ab)&ismember(smeta.tag2,t2)& smeta.bad==0&ismember(smeta.gt,g) &smeta.tp==0);
            conSmpAll=find(ismember(smeta.ab,ab)&ismember(smeta.tag2,t2)& smeta.bad==0 &smeta.tp==0);
            finalSample=find(ismember(smeta.ab,ab)&ismember(smeta.tag2,t2)&ismember(smeta.gt,g)& smeta.bad==0 &smeta.tp==90);
            crAll{c}(cc,:)=nan(1,numel(typeNames));
            pAll{c}(cc,:)=nan(1,numel(typeNames));
            tAll{c}(cc,:)=nan(1,numel(typeNames));
            pvAll{c}(cc,:)=nan(1,numel(typeNames));
            
            if numel(intSmp)>0
                for t=1:numel(typeColl)
                    if numel(typeColl{t})>0
                        %dTemp=log2(typeColl{t}(:,intSmp)+0.1)-mean(log2(typeColl{t}(:,intSmp)+0.1),2);
                        dTemp=log2(typeColl{t}(:,intSmp)+0.5)-mean(log2(typeColl{t}(:,intSmp)+0.5),2);
                        %dTemp=log2(typeColl{t}(:,intSmp)+0.1)-mean(log2(typeColl{t}(:,conSmp)+0.1),2);
                        %dTemp=dTemp./(std(dTemp,[],2)+0.1);
                        [crAll{c}(cc,t),pvAll{c}(cc,t)]=corr(dTemp(:),acol(repmat(smeta.tp(intSmp)',height(dTemp),1)),'rows','complete');
                        pFit=polyfit(smeta.tp(intSmp),mean(dTemp,1,'omitnan'),1);%robustfit(smeta.tp(intSmp),mean(dTemp,1,'omitnan'))
                        pAll{c}(cc,t)=pFit(1);
                        %pAll{c}(cc,t)=mean(mean(log2(typeColl{t}(:,intSmp)+0.1),2)-mean(log2(typeColl{t}(:,conSmp)+0.1),2),'omitnan')/90;
                        tAll{c}(cc,t)=mean(-log10(mattest(typeColl{t}(:,[finalSample,finalSample]),typeColl{t}(:,conSmpAll))),'omitnan');
                    end
                end
            end
        end
    end
end

chapBind=nan(numel(typeNames),numel(allChaps));
[~,idx]=ismember(extractAfter(allChaps,2),lower(chapNames));
chapBind(:,idx>0)=typeChap(:,idx(idx>0));
chapEnr=log2(chapBind)-mean(log2(chapBind),2,'omitnan')

nType=accumarray(typeVec(typeVec>0),1);
chapOrder=[5,11,1,15,6,14,13,3,9,2,16,4,12,7,8,10] %optimalleaforder(chapLink,chapDist)
typeOrder=[12,2,3,4,5,13,1,6,8,9:11];

xMat=repmat([1:numel(typeOrder)],numel(chapOrder),1)-0.125;
yMat=repmat([1:numel(chapOrder)]',1,numel(typeOrder));

sMat=repmat(rescale(log2(nType(typeOrder)),50,200)',numel(chapOrder),1)
load('cMaps.mat')
%figure('Color',[1 1 1],'Renderer','painters','Position',[2127 288 560 361])

cMap=brighten(brewermap(128,'YlOrRd'),0)
c=0
for t2={'hht2','h2b'}
    for ab={'myc','ha'}
        c=c+1;
        subplot(2,2,c)
        hold off
         colorIdx=acol(round(rescale(chapEnr(typeOrder,chapOrder)',1,128,'InputMin',0,'InputMax',1.5)))
         sizeVal=acol(round(rescale(chapEnr(typeOrder,chapOrder)',1,64,'InputMin',0,'InputMax',1.5)))
         colorBin=ones(height(colorIdx),3);%        
        colorBin(colorIdx>0,:)=cMap(colorIdx(colorIdx>0),:);
        %image(permute(reshape(colorBin,height(chapEnr),width(chapEnr),3),[2,1,3]))
        %imagesc(pAll{c}(chapOrder,newOrder)*90,[-1.1 1])
        tempMat=pAll{c}(chapOrder,typeOrder)*90;
        sizeMat=rescale(abs(crAll{c}(chapOrder,typeOrder)),10,150,'InputMax',1,'InputMin',0)
        %scatter(xMat(:),yMat(:),sMat(:),tempMat(:),'filled','MarkerEdgeColor',[1 1 1].*0.35)
        scatter(xMat(:),yMat(:),sizeMat(:),tempMat(:),'filled','MarkerEdgeColor',[1 1 1].*0.35)
        hold on
        scatter(xMat(:)+0.4,yMat(:),sizeVal,[1  1 1].*.45,'|','linewidth',2)
        scatter(xMat(isnan(tempMat)),yMat(isnan(tempMat)),sMat(isnan(tempMat))+10,[1 1 1],'filled')
        
        
        
        set(gca,'xtick',1:14,'Xticklabel',typeNames(typeOrder),'Ytick',1:numel(chapOrder),'Yticklabel',allChaps(chapOrder),'Ydir','reverse')
        title(sprintf('%s - %s',t2{1},ab{1}))
        ylabel(colorbar(),'mean(\Delta log2) in 90 mins')
        if  ismember(t2,{'h2b'})
            colormap(gca,(brighten(gy2pu,0.5)));
        else
            colormap(gca,(brighten(gy2bl,0.5)));
        end
        caxis([-.8 .8])
        plot([1.5:numel(typeOrder)].*[1;1],ylim','-','Color',[1 1 1].*.65)
        plot(xlim',[1.5:numel(chapOrder)].*[1;1],'-','Color',[1 1 1].*.65)
        xlim([.5 numel(typeOrder)+.5])
        ylim([.5 numel(chapOrder)+.5])
        axis square
    end
end
sgtitle('log change')
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/svg/Fig1SummaryCRScaleNONDR.fig')



%% S2A

%correlaion of chc data
[~,idx]=ismember(extractAfter(allChaps,'aa'),lower(chapNames));
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/promoterIDXvecFimp.mat','promoterIDXvecF')
crAll=corr(movmean(normProfile(promoterIDXvecF<3,:),150,1))
subplot(1,1,1)
imagesc(crAll(idx(idx>0),idx(idx>0)))
set(gca,'Ytick',1:sum(idx>0),'Yticklabel',chapNames(idx(idx>0)))
bL=[1:14]
hold on
plot([bL;bL]+.5,ylim','k-')
plot(ylim',[bL;bL]+.5,'k-')
sgtitle('chapCorr-150bp')
saveas(gcf,sprintf('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/svg/ChapcrMat.fig'))


%% Figure S2B
clear all
load('dataFS-single.mat')
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/promoterIDXvecFimp.mat','promoterIDXvecF')

allChaps=unique(tcs.gt);
chapOrder=[5,11,1,15,6,14,13,3,9,2,16,4,12,7,8,10] 
allChaps=allChaps(chapOrder)';
clear chapOrder;

[~,smeta.chapId]=ismember(smeta.gt,allChaps)
[uSmp,~,smeta.uid]=unique(smeta(:,{'ab','tag2','chapId','gt','tp'}))



uSmp
c=0;
for t2={'hht2','h2b'}
    for ab={'myc','ha'}
        c=c+1;
        
        allZero=find(uSmp.tp==0&ismember(uSmp.ab,ab)&ismember(uSmp.tag2,t2));
        selSmp=ismember(smeta.uid,allZero)&smeta.bad==0
        
        subplot(2,2,c)
        imagesc(corr(movmean(data(promoterIDXvecF<3,selSmp),250)))
        set(gca,'xtick',1:sum(selSmp),'Xticklabel',smeta.name(selSmp),'Ytick',1:sum(selSmp),'Yticklabel',smeta.name(selSmp))
        colorbar()
    end
end



figure('Color',[1 1 1],'Renderer','painters','Position',[2127 288 560 361])
allZero=find(uSmp.tp==0);
selSmp=ismember(smeta.uid,allZero)&smeta.bad==0
imagesc(corr(movmean(data(promoterIDXvecF<3,selSmp),250)))
set(gca,'xtick',1:sum(selSmp),'Xticklabel',smeta.name(selSmp),'Ytick',1:sum(selSmp),'Yticklabel',smeta.name(selSmp))
caxis([0 1])

allZero=find(uSmp.tp==0);
zeroData=cell2mat(arrayfun(@(u)mean(data(:,smeta.uid==u&smeta.bad==0),2),allZero','UniformOutput',false));
imagesc(corr(movmean(zeroData(promoterIDXvecF<3,:),250)),[0 1])
set(gca,'xtick',1:numel(allZero),'Xticklabel',strcat(uSmp.ab(allZero),'-',uSmp.tag2(allZero),'-',uSmp.gt(allZero)),'Ytick',1:numel(allZero),'Yticklabel',strcat(uSmp.ab(allZero),'-',uSmp.tag2(allZero),'-',uSmp.gt(allZero)))
colorbar()
title('correlation means smoothed(250bp)')
[~,bL]=unique(uSmp((allZero),{'ab','tag2'}))
bL=bL'-0.5;
hold on
plot([bL;bL],ylim','k-','Linewidth',1)
plot(ylim',[bL;bL],'k-','Linewidth',1)

bL=[1:numel(allZero)]+0.5;
hold on
plot([bL;bL],ylim','-','Color',[1 1 1].*0.75)
plot(ylim',[bL;bL],'-','Color',[1 1 1].*0.75)


saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/svg/FigS1crAllZero.fig')