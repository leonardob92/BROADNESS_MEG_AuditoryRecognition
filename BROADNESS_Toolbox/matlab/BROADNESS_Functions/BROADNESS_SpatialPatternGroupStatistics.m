function [SPATIAL_GROUP_STATS,FIGURES] = BROADNESS_SpatialPatternGroupStatistics(BROADNESS,group_labels,varargin)
% BROADNESS_SPATIALPATTERNGROUPSTATISTICS Independent-group spatial maps.
%
% Participant maps must share ONE common BROADNESS decomposition.
% Input maps: sources x components x conditions x participants.
% Tests are between independent participants, separately per condition.
%
% Two groups: pairwise Welch t maps; three or more groups: additionally
% omnibus Welch one-way F maps. Pairwise effect size is Hedges' g.
% Spatial clusters use a heteroscedasticity-aware, group-centered wild
% bootstrap (Rademacher multipliers), NOT paired sign flips or ordinary
% label permutations. It is approximate, especially for small groups.
% Corrected cluster p-values control maximum cluster mass over sources
% and selected components WITHIN each condition and contrast; separate
% conditions and contrasts are separate families. Significant clusters
% do not establish significance of each constituent source.
%
% Options: 'principalcomps', 'conditionnames', 'groupnames', 'mni_coords',
% 'mapnormalization' ('none'/'rms'), 'permutations', 'clusteralpha',
% 'alpha', 'neighbourdistance', 'randomseed', 'figuremode', 'figurelayout',
% 'OutputPath', 'figureformats', 'figureprefix'.
% Requires Statistics and Machine Learning Toolbox (tcdf, tinv, fcdf,
% finv); also uses BROADNESS_FDRCorrection and figure helpers.
%
% For final inference, use >=5000 bootstrap replicates, examine model
% assumptions, and validate the procedure on representative/null data.

opts = struct('principalcomps',[1 2],'conditionnames',{{}}, ...
    'groupnames',{{}},'mni_coords',[],'mapnormalization','none', ...
    'permutations',1000,'clusteralpha',0.01,'alpha',0.05, ...
    'neighbourdistance',[],'randomseed',1,'figuremode','off', ...
    'figurelayout','both','outputpath',[],'figureformats',{{'png'}}, ...
    'figureprefix','');
if mod(numel(varargin),2) ~= 0
    error('Optional arguments must be name-value pairs.');
end
for k = 1:2:numel(varargin)
    key = lower(char(string(varargin{k})));
    if ~isfield(opts,key), error('Unknown option: %s',key); end
    opts.(key) = varargin{k+1};
end
if ~isstruct(BROADNESS) || ~isfield(BROADNESS,'ActivationPatterns_BrainNetworks_Individual')
    error('BROADNESS must contain ActivationPatterns_BrainNetworks_Individual.');
end
maps = BROADNESS.ActivationPatterns_BrainNetworks_Individual;
if ~isnumeric(maps) || ndims(maps) ~= 4
    error('Individual maps must be sources x components x conditions x participants.');
end
[nSources,nComponents,nConditions,nParticipants] = size(maps);
components = opts.principalcomps(:)';
if isempty(components) || any(components<1 | components>nComponents | fix(components)~=components) || numel(unique(components))~=numel(components)
    error('Invalid principalcomps.');
end
if ~isnumeric(group_labels) || numel(group_labels)~=nParticipants || any(~isfinite(group_labels(:))) || any(group_labels(:)<1) || any(fix(group_labels(:))~=group_labels(:))
    error('group_labels must contain one positive integer per participant.');
end
group_labels = group_labels(:);
nGroups = max(group_labels);
if ~isequal(unique(group_labels)',1:nGroups), error('Group labels must be consecutive: 1,2,...'); end
groupIndices = cell(nGroups,1);
for g=1:nGroups
    groupIndices{g} = find(group_labels==g);
    if numel(groupIndices{g})<3
        error('Each group needs at least three participants for this bootstrap example.');
    end
end
if nGroups<2, error('At least two groups are required.'); end
if ~ismember(lower(char(opts.mapnormalization)),{'none','rms'})
    error('mapnormalization must be none or rms.');
end
if opts.permutations<3 || fix(opts.permutations)~=opts.permutations
    error('permutations must be an integer >= 3.');
end
if opts.alpha<=0 || opts.alpha>=1 || opts.clusteralpha<=0 || opts.clusteralpha>=1
    error('alpha and clusteralpha must be between 0 and 1.');
end
if ~isscalar(opts.randomseed) || ~isfinite(opts.randomseed) || fix(opts.randomseed)~=opts.randomseed
    error('randomseed must be a finite integer.');
end
conditionNames = make_names(opts.conditionnames,nConditions,'Condition');
groupNames = make_names(opts.groupnames,nGroups,'Group');
if isempty(opts.mni_coords)
    f=which('MNI152_8mm_coord_dyi.mat');
    if isempty(f), error('Provide mni_coords or add the bundled MNI file to the path.'); end
    s=load(f,'MNI8'); coords=s.MNI8;
else
    coords=opts.mni_coords;
end
if ~isnumeric(coords) || ~isequal(size(coords),[nSources 3])
    error('mni_coords must be sources x 3.');
end
valid=all(isfinite(coords),2);
if nnz(valid)<2 || any(~isfinite(maps(valid,components,:,:)),'all')
    error('Insufficient valid MNI sources or nonfinite spatial maps.');
end
validCoords=coords(valid,:);
if isempty(opts.neighbourdistance)
    steps=[];
    for d=1:3
        z=diff(unique(validCoords(:,d)));
        steps=[steps;z(z>1e-6)]; %#ok<AGROW>
    end
    if isempty(steps), error('Cannot infer spatial grid spacing.'); end
    neighbourDistance=min(steps);
else
    neighbourDistance=opts.neighbourdistance;
end
if ~isscalar(neighbourDistance) || ~isfinite(neighbourDistance) || neighbourDistance<=0
    error('Invalid neighbourdistance.');
end
neighbours=make_neighbours(validCoords,neighbourDistance);
if ~any(cellfun(@(x)~isempty(x),neighbours))
    error('No spatial neighbours found.');
end
rawMaps=maps(valid,components,:,:);
analysisMaps=rawMaps;
if strcmpi(opts.mapnormalization,'rms')
    for c=1:numel(components)
        for cond=1:nConditions
            for p=1:nParticipants
                v=rawMaps(:,c,cond,p);
                scale=sqrt(mean(v.^2));
                if scale==0, error('Cannot RMS-normalize an all-zero map.'); end
                analysisMaps(:,c,cond,p)=v./scale;
            end
        end
    end
end
oldRng=rng;
cleanupRng=onCleanup(@()rng(oldRng)); %#ok<NASGU>
rng(opts.randomseed,'twister');
ns=nnz(valid); nc=numel(components);
pairs=nchoosek(1:nGroups,2); np=size(pairs,1);
% Pairwise arrays: sources x selected components x conditions x group pairs.
shape=[nSources nc nConditions np];
meanDiff=nan(shape); tMap=nan(shape); gMap=nan(shape); rawP=nan(shape);
clusterLabels=zeros(shape); clusterP=nan(shape); significant=false(shape);
clusters=cell(nc,nConditions,np);
nullMass=zeros(opts.permutations,nConditions,np);
% Omnibus arrays: sources x selected components x conditions.
omniEnabled=nGroups>2;
FMap=nan(nSources,nc,nConditions); omniP=nan(size(FMap));
omniLabels=zeros(size(FMap)); omniClusterP=nan(size(FMap));
omniSignificant=false(size(FMap)); omniClusters=cell(nc,nConditions);
omniNull=zeros(opts.permutations,nConditions);

fprintf('Computing spatial group comparisons (%d groups, %d conditions).\n',nGroups,nConditions);
for cond=1:nConditions
    for pair=1:np
        g1=pairs(pair,1); g2=pairs(pair,2);
        observed=cell(nc,1); null1=cell(nc,1); null2=cell(nc,1);
        for c=1:nc
            x=reshape(analysisMaps(:,c,cond,groupIndices{g1}),ns,[]);
            y=reshape(analysisMaps(:,c,cond,groupIndices{g2}),ns,[]);
            [t,df,delta,hedges]=welch_t(x,y);
            observed{c}=t;
            sourceIds=find(valid);
            meanDiff(sourceIds,c,cond,pair)=delta;
            tMap(sourceIds,c,cond,pair)=t;
            gMap(sourceIds,c,cond,pair)=hedges;
            rawP(sourceIds,c,cond,pair)=2*tcdf(-abs(t),df);
            % Group-centered, leverage-adjusted residuals; use the same
            % participant multipliers across sources and components.
            null1{c}=bsxfun(@minus,x,mean(x,2))*sqrt(size(x,2)/(size(x,2)-1));
            null2{c}=bsxfun(@minus,y,mean(y,2))*sqrt(size(y,2)/(size(y,2)-1));
        end
        maximum=zeros(opts.permutations,1);
        for b=1:opts.permutations
            e1=2*(rand(numel(groupIndices{g1}),1)>.5)-1;
            e2=2*(rand(numel(groupIndices{g2}),1)>.5)-1;
            for c=1:nc
                x=bsxfun(@times,null1{c},e1');
                y=bsxfun(@times,null2{c},e2');
                [tb,dfb]=welch_t(x,y);
                threshold=tinv(1-opts.clusteralpha/2,dfb);
                [~,mass]=label_clusters(tb,threshold,neighbours,true);
                if ~isempty(mass), maximum(b)=max(maximum(b),max(mass)); end
            end
        end
        nullMass(:,cond,pair)=maximum;
        for c=1:nc
            x=reshape(analysisMaps(:,c,cond,groupIndices{g1}),ns,[]);
            y=reshape(analysisMaps(:,c,cond,groupIndices{g2}),ns,[]);
            [~,df]=welch_t(x,y);
            threshold=tinv(1-opts.clusteralpha/2,df);
            [lab,masses,signs]=label_clusters(observed{c},threshold,neighbours,true);
            sourceIds=find(valid);
            clusterLabels(sourceIds,c,cond,pair)=lab;
            current=empty_clusters();
            for j=1:numel(masses)
                corrected=(1+sum(maximum>=masses(j)))/(opts.permutations+1);
                ids=sourceIds(lab==j);
                clusterP(ids,c,cond,pair)=corrected;
                significant(ids,c,cond,pair)=corrected<=opts.alpha;
                current(j)=make_cluster(j,signs(j),masses(j),corrected,ids,coords,opts.alpha); %#ok<AGROW>
            end
            clusters{c,cond,pair}=current;
        end
    end
    if omniEnabled
        observed=cell(nc,1); residuals=cell(nc,nGroups);
        for c=1:nc
            xGroups=cell(nGroups,1);
            for g=1:nGroups
                xGroups{g}=reshape(analysisMaps(:,c,cond,groupIndices{g}),ns,[]);
                residuals{c,g}=bsxfun(@minus,xGroups{g},mean(xGroups{g},2))* ...
                    sqrt(size(xGroups{g},2)/(size(xGroups{g},2)-1));
            end
            [f,df2]=welch_f(xGroups);
            observed{c}=f;
            sourceIds=find(valid);
            FMap(sourceIds,c,cond)=f;
            omniP(sourceIds,c,cond)=fcdf(f,nGroups-1,df2,'upper');
        end
        maximum=zeros(opts.permutations,1);
        for b=1:opts.permutations
            multipliers=cell(nGroups,1);
            for g=1:nGroups
                multipliers{g}=2*(rand(numel(groupIndices{g}),1)>.5)-1;
            end
            for c=1:nc
                boot=cell(nGroups,1);
                for g=1:nGroups
                    boot{g}=bsxfun(@times,residuals{c,g},multipliers{g}');
                end
                [fb,dfb]=welch_f(boot);
                threshold=finv(1-opts.clusteralpha,nGroups-1,dfb);
                [~,mass]=label_clusters(fb,threshold,neighbours,false);
                if ~isempty(mass), maximum(b)=max(maximum(b),max(mass)); end
            end
        end
        omniNull(:,cond)=maximum;
        for c=1:nc
            xGroups=cell(nGroups,1);
            for g=1:nGroups
                xGroups{g}=reshape(analysisMaps(:,c,cond,groupIndices{g}),ns,[]);
            end
            [~,df2]=welch_f(xGroups);
            threshold=finv(1-opts.clusteralpha,nGroups-1,df2);
            [lab,masses]=label_clusters(observed{c},threshold,neighbours,false);
            sourceIds=find(valid);
            omniLabels(sourceIds,c,cond)=lab;
            current=empty_clusters();
            for j=1:numel(masses)
                corrected=(1+sum(maximum>=masses(j)))/(opts.permutations+1);
                ids=sourceIds(lab==j);
                omniClusterP(ids,c,cond)=corrected;
                omniSignificant(ids,c,cond)=corrected<=opts.alpha;
                current(j)=make_cluster(j,1,masses(j),corrected,ids,coords,opts.alpha); %#ok<AGROW>
            end
            omniClusters{c,cond}=current;
        end
    end
    fprintf('Spatial condition %d / %d completed.\n',cond,nConditions);
end

% Six summary measures are computed from ORIGINAL, unnormalized maps.
summaryNames={'MapStrength','CentroidX','CentroidY','CentroidZ', ...
    'SpatialDispersion','Lateralization'};
summaryValues=nan(6,nc,nConditions,nParticipants);
left=validCoords(:,1)<0; right=validCoords(:,1)>0;
for c=1:nc
    for cond=1:nConditions
        for p=1:nParticipants
            v=rawMaps(:,c,cond,p);
            w=abs(v); total=sum(w);
            strength=sqrt(mean(v.^2));
            centroid=[NaN NaN NaN]; dispersion=NaN;
            if total>0
                centroid=sum(bsxfun(@times,validCoords,w),1)/total;
                dispersion=sqrt(sum(w.*sum(bsxfun(@minus,validCoords,centroid).^2,2))/total);
            end
            l=sum(w(left)); r=sum(w(right)); laterality=NaN;
            if l+r>0, laterality=(l-r)/(l+r); end
            summaryValues(:,c,cond,p)=[strength centroid dispersion laterality];
        end
    end
end
summaryShape=[6 nc nConditions np];
summaryT=nan(summaryShape); summaryG=nan(summaryShape);
summaryP=nan(summaryShape); summaryDiff=nan(summaryShape);
summaryF=nan(6,nc,nConditions); summaryOverallP=nan(size(summaryF));
for m=1:6
    for c=1:nc
        for cond=1:nConditions
            groups=cell(nGroups,1);
            for g=1:nGroups
                v=reshape(summaryValues(m,c,cond,groupIndices{g}),1,[]);
                groups{g}=v(isfinite(v));
            end
            if omniEnabled && all(cellfun(@numel,groups)>=2) && ...
                    all(cellfun(@(v) var(v)>0,groups))
                [f,df2]=welch_f(groups);
                summaryF(m,c,cond)=f;
                summaryOverallP(m,c,cond)=fcdf(f,nGroups-1,df2,'upper');
            end
            for pair=1:np
                a=groups{pairs(pair,1)}; b=groups{pairs(pair,2)};
                if numel(a)<2 || numel(b)<2 || var(a)==0 || var(b)==0
                    continue
                end
                [t,df,delta,hedges]=welch_t(a,b);
                summaryT(m,c,cond,pair)=t;
                summaryP(m,c,cond,pair)=2*tcdf(-abs(t),df);
                summaryDiff(m,c,cond,pair)=delta;
                summaryG(m,c,cond,pair)=hedges;
            end
        end
    end
end
[summarySig,summaryAdj,summaryCrit]=BROADNESS_FDRCorrection(summaryP,opts.alpha);
if omniEnabled
    [summaryOverallSig,summaryOverallAdj,summaryOverallCrit]= ...
        BROADNESS_FDRCorrection(summaryOverallP,opts.alpha);
else
    summaryOverallSig=false(size(summaryF));
    summaryOverallAdj=nan(size(summaryF)); summaryOverallCrit=NaN;
end

% Figures follow the original toolbox's spatial-map and summary style.
figureSettings=BROADNESS_FigureSettings(opts.figuremode,opts.figurelayout, ...
    opts.outputpath,opts.figureformats,opts.figureprefix,'SpatialPatternGroupStatistics');
handles=gobjects(0); files={};
if ~strcmpi(opts.figuremode,'off')
    if figureSettings.MakeIndividual
        for c=1:nc
            for cond=1:nConditions
                if omniEnabled
                    [h,f]=plot_map(coords,valid,FMap(:,c,cond), ...
                        omniSignificant(:,c,cond), ...
                        sprintf('Network %d: %s - overall groups',components(c),conditionNames{cond}), ...
                        false,figureSettings,sprintf('SpatialGroup_BN%d_C%d_Overall',components(c),cond));
                    handles(end+1,1)=h; files=[files;f]; %#ok<AGROW>
                end
                for pair=1:np
                    label=sprintf('Network %d: %s - %s vs %s',components(c), ...
                        conditionNames{cond},groupNames{pairs(pair,1)},groupNames{pairs(pair,2)});
                    [h,f]=plot_map(coords,valid,tMap(:,c,cond,pair), ...
                        significant(:,c,cond,pair),label,true,figureSettings, ...
                        sprintf('SpatialGroup_BN%d_C%d_G%d_vs_G%d',components(c),cond,pairs(pair,1),pairs(pair,2)));
                    handles(end+1,1)=h; files=[files;f]; %#ok<AGROW>
                end
            end
        end
    end
    if figureSettings.MakeSummary
        for c=1:nc
            for cond=1:nConditions
                [h,f]=plot_summary(summaryValues(:,c,cond,:),summaryNames, ...
                    groupIndices,groupNames,components(c),conditionNames{cond}, ...
                    figureSettings,sprintf('SpatialGroupSummary_BN%d_C%d',components(c),cond));
                handles(end+1,1)=h; files=[files;f]; %#ok<AGROW>
            end
        end
    end
end
FIGURES.Handles=handles; FIGURES.Files=files;

SPATIAL_GROUP_STATS=struct();
SPATIAL_GROUP_STATS.MNI_Coordinates=coords;
SPATIAL_GROUP_STATS.ValidSources=valid;
SPATIAL_GROUP_STATS.Voxelwise.Pairwise.GroupPairs=pairs;
SPATIAL_GROUP_STATS.Voxelwise.Pairwise.MeanDifference=meanDiff;
SPATIAL_GROUP_STATS.Voxelwise.Pairwise.TStatistics=tMap;
SPATIAL_GROUP_STATS.Voxelwise.Pairwise.HedgesG=gMap;
SPATIAL_GROUP_STATS.Voxelwise.Pairwise.PValues=rawP;
SPATIAL_GROUP_STATS.Voxelwise.Pairwise.ClusterLabels=clusterLabels;
SPATIAL_GROUP_STATS.Voxelwise.Pairwise.ClusterCorrectedP=clusterP;
SPATIAL_GROUP_STATS.Voxelwise.Pairwise.Significant=significant;
SPATIAL_GROUP_STATS.Voxelwise.Pairwise.Clusters=clusters;
SPATIAL_GROUP_STATS.Voxelwise.Pairwise.NullMaximumClusterMass=nullMass;
SPATIAL_GROUP_STATS.Voxelwise.Pairwise.CorrectionFamily= ...
    'Sources and selected components, separately for each condition and group pair';
SPATIAL_GROUP_STATS.Voxelwise.Omnibus.Enabled=omniEnabled;
SPATIAL_GROUP_STATS.Voxelwise.Omnibus.FStatistics=FMap;
SPATIAL_GROUP_STATS.Voxelwise.Omnibus.PValues=omniP;
SPATIAL_GROUP_STATS.Voxelwise.Omnibus.ClusterLabels=omniLabels;
SPATIAL_GROUP_STATS.Voxelwise.Omnibus.ClusterCorrectedP=omniClusterP;
SPATIAL_GROUP_STATS.Voxelwise.Omnibus.Significant=omniSignificant;
SPATIAL_GROUP_STATS.Voxelwise.Omnibus.Clusters=omniClusters;
SPATIAL_GROUP_STATS.Voxelwise.Omnibus.NullMaximumClusterMass=omniNull;
SPATIAL_GROUP_STATS.Voxelwise.Omnibus.CorrectionFamily= ...
    'Sources and selected components, separately for each condition';
SPATIAL_GROUP_STATS.Summary.Names=summaryNames;
SPATIAL_GROUP_STATS.Summary.Values=summaryValues;
for m=1:6
    SPATIAL_GROUP_STATS.Summary.(summaryNames{m})=squeeze(summaryValues(m,:,:,:));
end
SPATIAL_GROUP_STATS.SummaryTests.Pairwise.GroupPairs=pairs;
SPATIAL_GROUP_STATS.SummaryTests.Pairwise.MeanDifference=summaryDiff;
SPATIAL_GROUP_STATS.SummaryTests.Pairwise.TStatistics=summaryT;
SPATIAL_GROUP_STATS.SummaryTests.Pairwise.HedgesG=summaryG;
SPATIAL_GROUP_STATS.SummaryTests.Pairwise.PValues=summaryP;
SPATIAL_GROUP_STATS.SummaryTests.Pairwise.AdjustedPValues=summaryAdj;
SPATIAL_GROUP_STATS.SummaryTests.Pairwise.Significant=summarySig;
SPATIAL_GROUP_STATS.SummaryTests.Pairwise.FDRCriticalP=summaryCrit;
SPATIAL_GROUP_STATS.SummaryTests.Pairwise.FDRMethod='Benjamini-Hochberg';
SPATIAL_GROUP_STATS.SummaryTests.Pairwise.FDRFamily= ...
    'Measures, selected components, conditions and group pairs';
SPATIAL_GROUP_STATS.SummaryTests.OverallANOVA.Enabled=omniEnabled;
SPATIAL_GROUP_STATS.SummaryTests.OverallANOVA.FStatistics=summaryF;
SPATIAL_GROUP_STATS.SummaryTests.OverallANOVA.PValues=summaryOverallP;
SPATIAL_GROUP_STATS.SummaryTests.OverallANOVA.AdjustedPValues=summaryOverallAdj;
SPATIAL_GROUP_STATS.SummaryTests.OverallANOVA.Significant=summaryOverallSig;
SPATIAL_GROUP_STATS.SummaryTests.OverallANOVA.FDRCriticalP=summaryOverallCrit;
SPATIAL_GROUP_STATS.SummaryTests.OverallANOVA.FDRMethod='Benjamini-Hochberg';
SPATIAL_GROUP_STATS.SummaryTests.OverallANOVA.FDRFamily= ...
    'Measures, selected components and conditions';
SPATIAL_GROUP_STATS.Settings.SelectedComponents=components;
SPATIAL_GROUP_STATS.Settings.GroupLabels=group_labels;
SPATIAL_GROUP_STATS.Settings.GroupNames=groupNames;
SPATIAL_GROUP_STATS.Settings.GroupPairs=pairs;
SPATIAL_GROUP_STATS.Settings.ConditionNames=conditionNames;
SPATIAL_GROUP_STATS.Settings.MapNormalization=opts.mapnormalization;
SPATIAL_GROUP_STATS.Settings.Permutations=opts.permutations;
SPATIAL_GROUP_STATS.Settings.ClusterFormingAlpha=opts.clusteralpha;
SPATIAL_GROUP_STATS.Settings.Alpha=opts.alpha;
SPATIAL_GROUP_STATS.Settings.NeighbourDistance=neighbourDistance;
SPATIAL_GROUP_STATS.Settings.RandomSeed=opts.randomseed;
SPATIAL_GROUP_STATS.Settings.ResamplingMethod= ...
    'Group-centered Rademacher wild bootstrap, studentized Welch statistics';
SPATIAL_GROUP_STATS.Settings.InferenceData='Unthresholded individual spatial maps';
SPATIAL_GROUP_STATS.Figures=FIGURES;
fprintf('Spatial group statistics completed.\n');
end

function [t,df,delta,g]=welch_t(x,y)
% x and y: sources x participants (or 1 x participants).
n1=size(x,2); n2=size(y,2);
m1=mean(x,2); m2=mean(y,2);
v1=var(x,0,2); v2=var(y,0,2);
a=v1/n1; b=v2/n2;
delta=m1-m2;
t=delta./sqrt(a+b);
df=(a+b).^2./(a.^2/(n1-1)+b.^2/(n2-1));
pooled=sqrt(((n1-1)*v1+(n2-1)*v2)/(n1+n2-2));
% First-order Hedges small-sample correction (same as time-series script).
g=(1-3/(4*(n1+n2)-9))*delta./pooled;
t(~isfinite(t))=0; df(~isfinite(df))=1;
g(~isfinite(g))=NaN;
end

function [f,df2]=welch_f(groups)
ng=numel(groups); n=zeros(ng,1);
ns=size(groups{1},1);
m=zeros(ns,ng); v=zeros(ns,ng);
for g=1:ng
    n(g)=size(groups{g},2);
    m(:,g)=mean(groups{g},2);
    v(:,g)=var(groups{g},0,2);
end
% Nonpositive variances cannot support the standard Welch approximation.
invalid=any(v<=0 | ~isfinite(v),2);
v=max(v,realmin);
w=bsxfun(@rdivide,n',v);
ws=sum(w,2);
wm=sum(w.*m,2)./ws;
correction=sum(bsxfun(@rdivide,(1-bsxfun(@rdivide,w,ws)).^2,(n'-1)),2);
df1=ng-1;
df2=(ng^2-1)./(3*correction);
f=(sum(w.*bsxfun(@minus,m,wm).^2,2)/df1)./ ...
    (1+2*(ng-2)*correction/(ng^2-1));
f(invalid | ~isfinite(f))=0;
df2(invalid | ~isfinite(df2))=1;
end

function [labels,masses,signs]=label_clusters(stat,threshold,neighbours,twoSided)
stat=stat(:); threshold=threshold(:);
labels=zeros(size(stat)); masses=[]; signs=[];
if twoSided
    masks={stat>=threshold,stat<=-threshold}; directions=[1 -1];
else
    masks={stat>=threshold}; directions=1;
end
for k=1:numel(masks)
    active=masks{k}; visited=false(size(stat)); seeds=find(active);
    for s=seeds(:)'
        if visited(s), continue; end
        queue=s; visited(s)=true; head=1; members=[];
        while head<=numel(queue)
            current=queue(head); head=head+1;
            members(end+1,1)=current; %#ok<AGROW>
            adj=neighbours{current};
            new=adj(active(adj) & ~visited(adj));
            if ~isempty(new)
                visited(new)=true;
                queue=[queue;new(:)]; %#ok<AGROW>
            end
        end
        j=numel(masses)+1;
        labels(members)=j;
        masses(j,1)=sum(abs(stat(members))); %#ok<AGROW>
        signs(j,1)=directions(k); %#ok<AGROW>
    end
end
end

function neighbours=make_neighbours(coords,distance)
n=size(coords,1); neighbours=cell(n,1);
norms=sum(coords.^2,2)';
for start=1:250:n
    last=min(n,start+249);
    a=coords(start:last,:);
    dist2=bsxfun(@plus,sum(a.^2,2),norms)-2*a*coords';
    for j=1:size(a,1)
        idx=start+j-1;
        neighbours{idx}=find(dist2(j,:) <= (distance*1.001)^2 & dist2(j,:)>1e-8)';
    end
end
end

function names=make_names(value,n,prefix)
if isempty(value)
    names=arrayfun(@(x)sprintf('%s %d',prefix,x),1:n,'UniformOutput',false);
elseif isstring(value)
    names=cellstr(value(:));
elseif iscell(value)
    names=value(:);
else
    error('Names must be a cell array or string array.');
end
if numel(names)~=n, error('Provide exactly %d %s names.',n,lower(prefix)); end
names=cellfun(@char,names,'UniformOutput',false);
end

function s=empty_clusters()
s=struct('Cluster',{},'Sign',{},'Mass',{},'CorrectedP',{}, ...
    'Significant',{},'SourceIndices',{},'MNI_Coordinates',{});
end

function s=make_cluster(j,sgn,mass,p,ids,coords,alpha)
s=struct('Cluster',j,'Sign',sgn,'Mass',mass,'CorrectedP',p, ...
    'Significant',p<=alpha, 'SourceIndices',ids,'MNI_Coordinates',coords(ids,:));
end

function [h,files]=plot_map(coords,valid,stat,sig,label,twoSided,settings,fileName)
h=figure('Visible',settings.Visible,'Color','w','Position',[100 100 1150 500]);
layout=tiledlayout(h,1,2,'TileSpacing','compact','Padding','compact');
ax=nexttile(layout);
scatter3(ax,coords(valid,1),coords(valid,2),coords(valid,3),12,stat(valid),'filled');
colorbar(ax);
if twoSided
    maximum=max(abs(stat(valid)));
    if maximum>0, caxis(ax,[-maximum maximum]); end
    title(ax,'Welch t-statistic');
else
    title(ax,'Welch F-statistic');
end
axis(ax,'equal'); axis(ax,'off'); view(ax,[-90 10]);
ax=nexttile(layout);
positive=sig & stat>0;
negative=sig & stat<0;
if any(positive)
    scatter3(ax,coords(positive,1),coords(positive,2),coords(positive,3), ...
        22,[.76 .16 .14],'filled'); hold(ax,'on');
end
if twoSided && any(negative)
    scatter3(ax,coords(negative,1),coords(negative,2),coords(negative,3), ...
        22,[.14 .32 .72],'filled');
end
if ~any(sig)
    text(ax,.5,.5,'No cluster survived correction','Units','normalized', ...
        'HorizontalAlignment','center');
end
title(ax,'Cluster-corrected differences');
axis(ax,'equal'); axis(ax,'off'); view(ax,[-90 10]);
title(layout,label,'Interpreter','none','FontWeight','normal');
rotate3d(h,'on');
files=BROADNESS_FinalizeFigure(h,settings,fileName);
end

function [h,files]=plot_summary(values,names,groupIndices,groupNames,component,condition,settings,fileName)
values=reshape(values,6,[]);
h=figure('Visible',settings.Visible,'Color','w','Position',[100 80 1150 700]);
layout=tiledlayout(h,2,3,'TileSpacing','compact','Padding','compact');
labels={'RMS weight','MNI x (mm)','MNI y (mm)','MNI z (mm)', ...
    'Distance (mm)','(Left-right)/(left+right)'};
for m=1:6
    ax=nexttile(layout); hold(ax,'on');
    for g=1:numel(groupIndices)
        x=values(m,groupIndices{g}); x=x(isfinite(x));
        scatter(ax,repmat(g,numel(x),1),x,13,'filled', ...
            'MarkerFaceAlpha',.4,'HandleVisibility','off');
        if ~isempty(x)
            errorbar(ax,g,mean(x),std(x)/sqrt(numel(x)),'ko', ...
                'MarkerFaceColor','w','LineWidth',1.4);
        end
    end
    xticks(ax,1:numel(groupNames)); xticklabels(ax,groupNames);
    xtickangle(ax,25); xlim(ax,[.5 numel(groupNames)+.5]);
    title(ax,names{m}); ylabel(ax,labels{m}); grid(ax,'on'); box(ax,'off');
end
title(layout,sprintf('Spatial summaries - Network %d, %s',component,condition), ...
    'FontWeight','normal');
files=BROADNESS_FinalizeFigure(h,settings,fileName);
end
