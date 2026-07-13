function [u, info] = confidence_weighted_graph_tv(observed, confidence, seed, varargin)
%CONFIDENCE_WEIGHTED_GRAPH_TV Edge-preserving confidence-weighted smoothing.
%   U minimizes a robust graph objective by iteratively reweighted Jacobi
%   updates. OBSERVED supplies the data term, CONFIDENCE controls its weight,
%   and SEED defines fixed bilateral graph weights. NaNs denote missing grid
%   locations and remain NaN. No ground-truth or material labels are used.

p = inputParser;
addParameter(p,'Lambda',0.25,@(x)isnumeric(x)&&isscalar(x)&&x>=0);
addParameter(p,'Iterations',30,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'EdgeSigma',0.20,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'TvEpsilon',0.03,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'HighConfidenceThreshold',0.80,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'HighConfidenceWeight',50,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'LowConfidenceWeight',0.5,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'PhysicalRange',[0.5 10],@(x)isnumeric(x)&&numel(x)==2);
parse(p,varargin{:}); opt=p.Results;

assert(isequal(size(observed),size(confidence),size(seed)), ...
    'Observed, confidence, and seed maps must have equal dimensions.');
valid=isfinite(observed)&isfinite(confidence)&isfinite(seed);
u=observed; u(~valid)=NaN;
data_weight=opt.LowConfidenceWeight+confidence.^2;
data_weight(confidence>=opt.HighConfidenceThreshold)=opt.HighConfidenceWeight;
data_weight(~valid)=0;

% Fixed bilateral graph weights prevent propagation across strong seed edges.
right_valid=valid(:,1:end-1)&valid(:,2:end);
down_valid=valid(1:end-1,:)&valid(2:end,:);
wr=exp(-((seed(:,2:end)-seed(:,1:end-1))/opt.EdgeSigma).^2).*right_valid;
wd=exp(-((seed(2:end,:)-seed(1:end-1,:))/opt.EdgeSigma).^2).*down_valid;
wr(~right_valid)=0; wd(~down_valid)=0;

for iter=1:round(opt.Iterations)
    numerator=data_weight.*observed;
    denominator=data_weight;
    dr=u(:,2:end)-u(:,1:end-1);
    dd=u(2:end,:)-u(1:end-1,:);
    wr_iter=wr./sqrt(dr.^2+opt.TvEpsilon^2);
    wd_iter=wd./sqrt(dd.^2+opt.TvEpsilon^2);

    value=opt.Lambda*wr_iter.*u(:,2:end);
    numerator(:,1:end-1)=numerator(:,1:end-1)+value;
    denominator(:,1:end-1)=denominator(:,1:end-1)+opt.Lambda*wr_iter;
    value=opt.Lambda*wr_iter.*u(:,1:end-1);
    numerator(:,2:end)=numerator(:,2:end)+value;
    denominator(:,2:end)=denominator(:,2:end)+opt.Lambda*wr_iter;

    value=opt.Lambda*wd_iter.*u(2:end,:);
    numerator(1:end-1,:)=numerator(1:end-1,:)+value;
    denominator(1:end-1,:)=denominator(1:end-1,:)+opt.Lambda*wd_iter;
    value=opt.Lambda*wd_iter.*u(1:end-1,:);
    numerator(2:end,:)=numerator(2:end,:)+value;
    denominator(2:end,:)=denominator(2:end,:)+opt.Lambda*wd_iter;

    next=numerator./max(denominator,eps);
    next=min(max(next,opt.PhysicalRange(1)),opt.PhysicalRange(2));
    next(~valid)=NaN; u=next;
end
info=struct('iterations',round(opt.Iterations),'lambda',opt.Lambda, ...
    'mean_data_weight',mean(data_weight(valid),'omitnan'), ...
    'mean_horizontal_edge_weight',mean(wr(right_valid),'omitnan'), ...
    'mean_vertical_edge_weight',mean(wd(down_valid),'omitnan'));
end
