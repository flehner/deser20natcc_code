close all
clear all

% ------------------------------------------------------------------------------
% fig2.m
%
% Description: script to create Fig. 2 in Deser et al., 2020,
% https://doi.org/10.1038/s41558-020-0731-2
%
% Note: this is essentially just a plotting script.
%
% The full data is too large to share on github, but is available via:
% - Large Ensembles: http://www.cesm.ucar.edu/projects/community-projects/MMLEA/
% - CMIP5: https://esgf-node.llnl.gov/search/cmip5/  (optional)
%
% Required data:
% - the MMLEA archive
%   (data has been regridded to common 2.5x2.5 degree grid prior,
%   identified by "g025"; this step needs to be done by the user)
% - CMIP5: https://esgf-node.llnl.gov/search/cmip5/  (optional)
%
% Required functions:
% - am.m
% - seasmean.m
%
% Author: Flavio Lehner, May 2020, flehner@ucar.edu
%
% ------------------------------------------------------------------------------

pathin      = '/Users/flehner/Dropbox/work/';
vars        = {'tas','pr'};
comp        = 'Amon';
seasons     = {'DJF','JJA','annual'};
hist_or_kernel = 2; % histogram or kernel density function
% -- location to plot:
% -- Upper Colorado
lat         = [38.75 41.25];
lon         = [248.75 253.75];

% -- CLIVAR parameters
models      = {'canesm2_lens','csiro_mk36_lens','gfdl_cm3_lens','cesm_lens','gfdl_esm2m_lens','ec_earth_lens','mpi_lens'};
model_names = {'CanESM2','CSIRO-Mk3-6-0','GFDL-CM3','CESM1-CAM5','GFDL-ESM2M','EC-EARTH','MPI-ESM'};
ensmem      = [50,30,20,40,30,16,100];
start0      = [1950,1850,1920,1920,1950,1860,1850];
ende0       = [2100,2100,2100,2100,2100,2100,2099];

% -- CMIP5 parameters
models_cmip5 = {'ACCESS1-0','ACCESS1-3','bcc-csm1-1-m','bcc-csm1-1','BNU-ESM',...
                'CanESM2','CCSM4','CESM1-BGC','CESM1-CAM5','CMCC-CESM','CMCC-CM',...
                'CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','FGOALS-g2',...
                'FIO-ESM','GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','GISS-E2-H-CC',...
                'GISS-E2-H','GISS-E2-R-CC','GISS-E2-R','HadGEM2-AO','HadGEM2-CC',...
                'HadGEM2-ES','inmcm4','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR',...
                'MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR',...
                'MRI-CGCM3','MRI-ESM1','NorESM1-ME','NorESM1-M'};
start_cmip5  = 1870;
ende_cmip5   = 2100;

% -- general paramaters
start       = 1950;
ende        = 2099;
time        = start:ende;
refstart    = 1971; % 1950 2001 1971
refende     = 2000; % 1979 2030 2000
start_trend = 1951;
ende_trend  = 2099;
time_trend  = start_trend:ende_trend;
lati        = [-88.75:2.5:88.75];
loni        = [1.25:2.5:358.75];
% -- find grid cells to plot in models
for x = 1:length(lat)
  ii(x)          = find(abs(lati-lat(x))==min(abs(lati-lat(x))));
end
for x = 1:length(lon)
  jj(x)          = find(abs(loni-lon(x))==min(abs(loni-lon(x))));
end

wl          = 10; % running mean length in years

% -- load data
for v = 2:2%length(vars)
  vari = vars{v}
  if strcmp(vari,'tas')==1
    units       = 'K';
    f           = 1;
  elseif strcmp(vari,'pr')==1
    units       = '%';
    f           = 86400;
  end
  for s = 3:3%length(seasons)
    seas = seasons{s}
    % -- load CMIP5 data --
    start_cmip5 = 1870;
    ende_cmip5  = 2100;
    for m = 1:length(models_cmip5)
      m
      filein  = [pathin 'cmip5-ng/' vari '/' vari '_mon_' models_cmip5{m} '_rcp85_r1i1p1_g025.nc'];
      tmp0    = ncread([filein],[vari])*f;
      tmp     = squeeze(nanmean(nanmean(squeeze(tmp0(jj,ii,(start_trend-start_cmip5)*12+1:(ende_trend-start_cmip5+1)*12)),2),1));
      if strcmp(seas,'annual')==1
        tmp    = am(tmp);
      else
        tmp    = seasmean(tmp,seas);
      end
      cmip5(m,:) = tmp;
    end
    % -- load CLIVAR data
    for m = 1:length(models)
      m
      filein  = [pathin models{m} '/' comp '/' vari '/' vari '_' comp '_' model_names{m} '_' seas '_g025.nc'];
      tmp0    = ncread([filein],'data')*f;
      for e = 1:ensmem(m)
        tmp     = squeeze(nanmean(nanmean(squeeze(tmp0(jj,ii,e,start_trend-start0(m)+1:ende_trend-start0(m)+1)),2),1));
        eval([models{m} '(e,:) = tmp;']);
        model_raw_ts(m,e,:) = tmp;
      end
    end

  end % -- end of seasons loop
end % -- end of vars loop








%% --- PLOTTING ---------------------------------------------------------------
close all

% -- rainbow-ish colorbar
cols = [255 0 0;
        255 160 16;
        255 224 32;
        0 192 0;
        80 208 255;
        0 32 255;
        160 32 255]/255;
% -- list of line styles:
lstyle = {'-',':','--','-.','-','--','-.'};


if strcmp(vari,'tas')==1
  xlim = [-3 5];
elseif strcmp(vari,'pr')==1
  xlim = [-.8 .8];
end

ylim1 = [-.4 .6];
ylim2 = [-.04 .06];

figure1 = figure;
set(figure1, 'units', 'centimeters', 'pos', [10 10 18 20])

subplot(2,4,[1 3])
hold on
string1 = ['(a) ' num2str(wl) '-yr running mean'];
string2 = ['relative to ' num2str(refstart) '-' num2str(refende)];
title([string1 ' ' string2],'FontSize',10)
anom = NaN(length(models),max(ensmem),length(time_trend));
for m = 1:length(models)
  ref(m)      = nanmean(nanmean(squeeze(model_raw_ts(m,1:ensmem(m),refstart-start+1:refende-start))));
  anom(m,1:ensmem(m),:) = (rm(squeeze(model_raw_ts(m,1:ensmem(m),:)),wl,2)'-ref(m))';
  idx         = ~isnan(anom(m,1,:));
  jbfill(time_trend(idx),prctile(squeeze(anom(m,:,idx)),95),prctile(squeeze(anom(m,:,idx)),5),cols(m,:),'none',1,.5)
  hold on
end
for m = 1:length(models)
  h0(m) = plot(time_trend(idx),squeeze(nanmean(anom(m,:,idx),2)),'Color',cols(m,:),'LineWidth',1.5,'LineStyle',lstyle{m})
end
set(gca,'YLim',ylim1)
hline(0,'k')
box on
legend([h0],[model_names],'Location','NorthWest','FontSize',8)
legend boxoff
xlabel('Time (Year)')
ylabel('Precipitation anomaly (mm/day)')

subplot(2,4,[5 7])
hold on
string1 = ['(b) Standard deviation of ' num2str(wl) '-yr running means'];
string2 = ['relative to ' num2str(refstart) '-' num2str(refende)];
title([string1 char(10) string2],'FontSize',10)
for m = 1:length(models)
  ref = nanmean(std(rm(squeeze(model_raw_ts(m,1:ensmem(m),refstart-start+1:refende-start)),wl,2)));
  ref = nanmean(ref(~isnan(ref)));
  std_ts(m,:) = std(rm(squeeze(model_raw_ts(m,1:ensmem(m),:)),wl,2))-ref;
  std_ts_percent(m,:) = (std(rm(squeeze(model_raw_ts(m,1:ensmem(m),:)),wl,2))./ref)*100-100;
end
h0 = jbfill(time_trend(idx),max(std_ts(:,idx)),min(std_ts(:,idx)),[.7 .7 .7],'none')
hold on
h1 = plot(time_trend(idx),nanmean(std_ts(:,idx)),'k','LineWidth',2)
set(gca,'YLim',ylim2)
hline(0,'k')
box on
legend([h0 h1],'Range across LEs','Mean across LEs','Location','NorthWest')
legend boxoff
xlabel('Time (Year)')
ylabel('Standard deviation change (mm/day)')

subplot(2,4,4)
hold on
title('2070-2099')
pos = [1 6 2 3 4 5 7]; % hand-arranged position of models in side plot
for m = 1:length(models)
  % [f0,x0,u] = ksdensity(nanmean(squeeze(anom(m,:,end-29:end)),2));  %,binranges);
  % plot(f0*.25,x0, 'Color',cols(m,:),'LineWidth',1.5)
  tmp = anom(m,:,end-29:end);
  plot([pos(m) pos(m)],[min(tmp(:)) max(tmp(:))],'Color',cols(m,:))
  plot([pos(m) pos(m)],[prctile(tmp(:),5) prctile(tmp(:),95)],'Color',cols(m,:),'LineWidth',3)
  plot([pos(m) pos(m)],[nanmean(tmp(:)) nanmean(tmp(:))],'o','Color',cols(m,:),'LineWidth',2,'MarkerSize',8)
  % plot([m-.1 m+.1],[nanmedian(tmp(:)) nanmedian(tmp(:))],'-','Color',cols(m,:),'LineWidth',2)
end
for m = 1:length(models_cmip5)
  ref_cmip5(m)      = nanmean(nanmean(squeeze(cmip5(m,refstart-start+1:refende-start))));
  anom_cmip5(m,:) = (rm(squeeze(cmip5(m,:)),wl,2)'-ref_cmip5(m))';
end
plot([8 8],[min(anom_cmip5(:,end-4)) max(anom_cmip5(:,end-4))],'Color',[.5 .5 .5])
plot([8 8],[prctile(anom_cmip5(:,end-4),5) prctile(anom_cmip5(:,end-4),95)],'Color',[.5 .5 .5],'LineWidth',3)
plot([8 8],[mean(anom_cmip5(:,end-4)) mean(anom_cmip5(:,end-4))],'o','Color',[.5 .5 .5],'LineWidth',2,'MarkerSize',8)
% plot([7.9 8.1],[median(anom_cmip5(:,end-4)) median(anom_cmip5(:,end-4))],'-','Color',[.5 .5 .5],'LineWidth',2)
text(6.2,.53,'CMIP5','Color',[.5 .5 .5])
% set(gca,'YLim',ylim1,'XTick',[])
set(gca,'Ylim',ylim1,'XLim',[0 9],'ycolor',[.999 .999 .999],'xcolor',[.999 .999 .999],'Ytick',[],'Xtick',[])
hline(0,'k')
box on

subplot(2,4,8)
hold on
title('2070-2099')
tmp = std_ts(:,end-29:end);
plot([1 1],[min(tmp(:)) max(tmp(:))],'Color',[0 0 0])
plot([1 1],[prctile(tmp(:),5) prctile(tmp(:),95)],'Color',[0 0 0],'LineWidth',3)
plot([1 1],[nanmean(tmp(:)) nanmean(tmp(:))],'o','Color',[0 0 0],'LineWidth',2,'MarkerSize',8)
set(gca,'YLim',ylim2,'XLim',[0 9],'ycolor',[.999 .999 .999],'xcolor',[.999 .999 .999],'Ytick',[],'Xtick',[])
hline(0,'k')
box on

% -- save figure (adjust as needed)
% set(gcf,'PaperPositionMode','auto');
% fileo = ['/Users/flehner/Dropbox/publication/clivar19_perspective/clivar_lens_variability_change_plots_' vari '_' seas '_' num2str(wl) 'yr'];
% print('-r300','-loose', '-depsc', ['' fileo '.eps'])
% save2pdf(['' fileo '.pdf'])
% saveas(gcf,fileo,'jpg')
