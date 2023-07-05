clear all
load('dataFS-single.mat')
clearvars -except smeta data tcs nucMeta intChaps chapNames normProfile GP
GP=load('/home/labs/barkailab/felixj/Documents/MATLAB/scripts/gene/group_imp.mat');

%% Figure S3B
chapBind=squeeze(sum(metaProfilePromLenDivya(normProfile,'promEnd','felixTss','afterTss',150,'promLen',300),2));
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/addData.mat','geneLvl','opnScore')
chapEnr=(chapBind-mean(chapBind,'omitnan'))./std(chapBind,'omitnan');
promVal=nan(6701,4,width(data));
for g=nucMeta.gene(nucMeta.order==1)'
    plus1=find(nucMeta.order==1&nucMeta.gene==g);
    promVal(g,4,:)=mean(data(nucMeta.pos(plus1)+[-50:50],:));
    ndr=(nucMeta.pos(plus1-GP.gene_infoR64.dir(g))+50*GP.gene_infoR64.dir(g)):GP.gene_infoR64.dir(g):(nucMeta.pos(plus1)-50*GP.gene_infoR64.dir(g));
    promVal(g,3,:)=mean(data(ndr,:));
    minus1= find(nucMeta.order==-1&nucMeta.gene==g);
    if numel(minus1)>0
        promVal(g,2,:)=mean(data(nucMeta.pos(minus1)+[-50:50],:));
    end
    minus2= find(nucMeta.order==-2&nucMeta.gene==g);
    if numel(minus2)>0
        promVal(g,1,:)=mean(data(nucMeta.pos(minus2)+[-50:50],:));
    end
end

segMat=promVal;
load('holstege.mat')
geneVar=std(holstege.TsXMut,[],2,'omitnan');
opnScoreS=opnScore;
opnScoreS(opnScoreS<0)=opnScoreS(opnScoreS<0)/3;
intParas={'geneLvl','geneVar','opnScoreS','ndrWidth'}
geneSegs=0
nBins=5;
allChaps=unique(tcs.gt)
chapOrder=[5,11,1,15,6,14,13,3,9,2,16,4,12,7,8,10] %optimalleaforder(chapLink,chapDist)
allChaps=allChaps(chapOrder);
clearvars crVal pVal chapOrder
sSel={[1,2],3,4};
close all
nSegs= size(sSel,2);
ySeg=5;
endTp=ones(size(allChaps))*90;
endTp(ismember(allChaps,'aanofrb'))=120;

for p=3
    currPar=eval(intParas{p});
    binEdge=quantile(currPar,[0.01,0.99]);
    [geneId,~]=discretize(max(min(currPar,binEdge(2)),binEdge(1)),nBins);
    geneId(isnan(currPar))=nan;
    expMat=repmat(accumarray(geneId(geneId>0),currPar(geneId>0),[],@mean),1,numel(sSel));
    %nucMat=repmat([-2 -1 1:geneSegs geneSegs+2],10,1);
    %nucMat=repmat([-2 -1 0],nBins,1);
    nucMat=repmat([-1 0,1],nBins,1);
    figure('Color',[1 1 1],'Renderer','painters','Position',[2127 288 560 361])
    
    subplot(1,5,1)
    hold off
    c=0
    for g=allChaps'
        c=c+1;
        sb=1;
        subplot(1,5,1)
        [~,currCh]=ismember(extractAfter(g,2),lower(chapNames));
        if currCh>0
            binBind=accumarray(geneId(geneId>0),chapEnr(geneId>0,currCh),[],@(x)mean(x,'omitnan'));
            scatter(ones(nBins,1).*0,expMat(:,1)-(ySeg)*c,40,binBind,'filled','MarkerEdgeColor',[1 1 1].*0.65)
            hold on
        end
        caxis([0 .75])
        xlim([-1 1].*.5)
        ylim(quantile(expMat(:,1),[0 1])-[numel(allChaps),1]*ySeg)
        plot(xlim,[-1.93].*[1 1]-(ySeg)*c,'-','Color',[1 1 1].*0.65)
        set(gca,'Ytick',fliplr(mean(expMat(:,1))-[1:numel(allChaps)]*ySeg),'Yticklabel',flipud(extractAfter(allChaps,2)))
        colormap(gca,brewermap(128,'YlOrRd'))
        clear crVal pVal crType
        ylabel('chapEnr')
        for t2={'hht2','h2b'}
            for ab={'myc','ha'}
                sb=sb+1;
                intSmp=find(ismember(smeta.tid,find(ismember(tcs.ab,ab)&ismember(tcs.tag2,t2)&ismember(tcs.gt,g))) & smeta.bad==0 &smeta.tp<=endTp(ismember(allChaps,g)));
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
                    caxis([-.8 .8])
                    xlim([-1 1].*1.5)
                    yticks([])
                    ylim(quantile(expMat(:,1),[0 1])-[numel(allChaps),1]*ySeg)
                    colormap(gca,flipud(brighten(brewermap(128,'RdBu'),0.5)));
                    ylabel([t2{1} '-' ab{1}])
                    xlabel('nucPosition')
                    %ylabel(colorbar(),'RapEffect-log2 in 90min')
                end
                
            end
        end
    end
    sgtitle(sprintf('%s - (-2,,-1), NDR, +1',intParas{p}))
end
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/svg/Fig2Summary.fig')

%% Figure S3C
clearvars -except smeta data tcs nucMeta allChaps chapNames normProfile GP 
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Naama/extData/addData.mat','geneLvl','opnScore')
opnScoreS=opnScore;
opnScoreS(opnScoreS<0)=opnScoreS(opnScoreS<0)/3;

promLen=500;
fiveLen=150;
metaPro=metaProfilePromLenDivya(data,'promEnd','felixTss','afterTss',fiveLen,'promLen',promLen);
chapPro=metaProfilePromLenDivya(normProfile,'promEnd','felixTss','afterTss',fiveLen,'promLen',promLen);

nBins=5;
binEdge=quantile(opnScore,[0.01,0.99]);
[geneId,~]=discretize(max(min(opnScore,binEdge(2)),binEdge(1)),nBins);
geneId(isnan(opnScore))=nan;
allTps=unique(smeta.tp(smeta.tp<=90))
cMaps=load('cMaps2.mat')

uTPs=unique(smeta.tp(smeta.tp<=90))
cIdx=round(linspace(1,128,numel(uTPs)));
cMap{1}=cMaps.gy2bl(cIdx,:);
cMap{5}=cMaps.gy2rd(cIdx,:);
groupName{5}='OPN';
groupName{1}='DPN';


allChaps={'aartt106','aasth1','aaino80','aaspt6','aaspt16'}';

figure('Color',[1 1 1],'Renderer','painters','Position',[2127 288 560 361])
c=0
for g=allChaps'
    c=c+1;
    [~,currCh]=ismember(extractAfter(g,2),lower(chapNames));
    subplot(8,10,c)
    title(extractAfter(g{1},2))
    if currCh>0
        
        hold off
        plot([-promLen:fiveLen],movmean(mean(chapPro(geneId==1,:,currCh),'omitnan'),21))
        hold on
        plot([-promLen:fiveLen],movmean(mean(chapPro(geneId==nBins,:,currCh),'omitnan'),21))
        axis tight
        title(chapNames(currCh))
    end
    selSmp=find(smeta.bad==0&ismember(smeta.gt,g)&ismember(smeta.ab,'myc') &ismember(smeta.tag2,'hht2')&smeta.tp<=90)
    for b=[1,5]
        c=c+1;
        [uVal,u1,u2]=unique(smeta.tp(selSmp));
        [~,cIdx]=ismember(uVal,uTPs)
        subplot(8,10,c)
        for i=1:numel(u1)
            tempMat(i,:)=movmean(mean(metaPro(geneId==b,:,selSmp(u2==i)),[3,1]),51);
            plot([-promLen:fiveLen],tempMat(i,:),'-','linewidth',1,'Color',cMap{b}(cIdx(i),:))
            hold on
        end
        axis tight
        ylim([0 max(ylim)])
        c=c+1;
        %set(gca,'ColorOrder',cMap{b}(2:end,:),'Color',[1 1 1].*.9)
        title(sprintf('Myc Dynamics  %s',groupName{b}))
        subplot(8,10,c)
        dMat=log2(tempMat+0.25)-log2(tempMat(1,:)+0.25);
        for i=1:numel(u1)
            plot([-promLen:fiveLen],dMat(i,:),'-','linewidth',1)
            hold on
        end
        title(sprintf('Myc \\Delta log2  %s',groupName{b}))
        set(gca,'ColorOrder',cMap{b}(2:end,:))
        axis tight
        ylim([-1 1])
    end
end
saveas(gcf,'/home/labs/barkailab/felixj/Documents/MATLAB/projects/GIlad-AA/figureScripts/Fig2Inds.fig')