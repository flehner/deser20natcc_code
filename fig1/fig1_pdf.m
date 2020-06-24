close all
clear all

% ------------------------------------------------------------------------------
% fig1_pdf.m
%
% Description: script to create the center panel of Fig. 1 in Deser et al., 2020,
% https://doi.org/10.1038/s41558-020-0731-2
%
% Note: this is essentially just a plotting script; spatially-averaged trends
% were calculated prior; the necessary time series are provided.
%
% The full data is too large to share on github, but is available via:
% - Large Ensembles: http://www.cesm.ucar.edu/projects/community-projects/MMLEA/
% - CMIP5: https://esgf-node.llnl.gov/search/cmip5/
%
% Required data:
% - fig1_pdf_data.tar
%   (spatially-averaged trend values for all Large Ensembles and CMIP5 models)
%
% Required functions:
% - am.m
% - jbfill.m
% - vline.m
%
% Author: Flavio Lehner, May 2020, flehner@ucar.edu
%
% ------------------------------------------------------------------------------

pathin      = '~/Dropbox/publication/clivar19_perspective/code_for_github/';
vars        = {'tas'};
comp        = 'Amon';
seasons     = {'annual'};
hist_or_kernel = 2; % histogram or kernel density function


% -- CLIVAR parameters
models      = {'cesm_lens','canesm2_lens','csiro_mk36_lens','gfdl_cm3_lens','gfdl_esm2m_lens','ec_earth_lens','mpi_lens'};
model_names = {'CESM1-CAM5','CanESM2','CSIRO-Mk3-6-0','GFDL-CM3','GFDL-ESM2M','EC-EARTH','MPI-ESM'};
ensmem      = [40,50,30,20,30,16,100];
start0      = [1920,1950,1850,1920,1950,1860,1850];
ende0       = [2100,2100,2100,2100,2100,2100,2099];
for m = 1:length(models)
  model_names_legend(m) = strcat(model_names(m),{' ('},num2str(ensmem(m)),{')'})
end

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
start_trend = 1951;
ende_trend  = 2010;
time_trend  = start_trend:ende_trend;


% -- load CLIVAR data
vari = vars{1};
if strcmp(vari,'tas')==1
  units       = 'K';
  f           = 1;
elseif strcmp(vari,'pr')==1
  units       = '%';
  f           = 86400;
end
seas = seasons{1};

% -- load domain-wide trends --
% -- LEs
for m = 1:length(models)
  m
  filein  = [pathin '/fig1_pdf_data/lens/' vari '_' comp '_' model_names{m} '_' seas '_trend_values_' num2str(start_trend) '-' num2str(ende_trend) '.nc'];
  tmp0    = ncread([filein],'trend_values');
  eval([models{m} ' = tmp0;']);
end
% -- CMIP5 (one member per model)
start_cmip5 = 1870;
ende_cmip5  = 2100;
for m = 1:length(models_cmip5)
  m
  filein  = [pathin 'fig1_pdf_data/cmip5/' vari '_mon_' models_cmip5{m} '_rcp85_r1i1p1_g025_ts.nc'];
  tmp     = ncread([filein],'NA')*f;
  if strcmp(seas,'annual')==1
    tmp    = am(tmp);
  else
    tmp    = seasmean(tmp,seas);
  end
  trend     = polyfit(time_trend,tmp(start_trend-start_cmip5+1:ende_trend-start_cmip5+1),1);
  cmip5(m)  = trend(1)*length(time_trend);
end
% -- CMIP5 (all members per model)
ec_count = ones(40,1);
cmip5_all = NaN(length(models_cmip5),11);
for m = 1:length(models_cmip5)
  m
  for e = 1:20
    filein = [pathin 'fig1_pdf_data/cmip5/' vari '_mon_' models_cmip5{m} '_rcp85_r' num2str(e) 'i1p1_g025_ts.nc'];
    if isfile(filein)
      ec_count(m) = ec_count(m)+1;
      tmp     = ncread([filein],'NA')*f;
      if strcmp(seas,'annual')==1
        tmp    = am(tmp);
      else
        tmp    = seasmean(tmp,seas);
      end
      trend     = polyfit(time_trend,tmp(start_trend-start_cmip5+1:ende_trend-start_cmip5+1),1);
      cmip5_all(m,e)  = trend(1)*length(time_trend);
    else
      cmip5_all(m,e) = NaN;
    end
  end
end

% -- obs --
'observations'
obs        = ncread([pathin 'fig1_pdf_data/observations/tas_' seas '_trend_values_' num2str(start_trend) '-' num2str(ende_trend) '.nc'],'trend_values');
obs_name  = {'BEST'};

% -- OLENS domain-wide trends --
% -- Observational Large Ensemble (McKinnon and Deser, 2018, https://doi.org/10.1175/JCLI-D-17-0901.1)
%    Note: not use in the final version of the figure, but potentially useful
'OLENS'
olens        = ncread([pathin 'fig1_pdf_data/olens_mckinnon/tas_' seas '_trend_values_' num2str(start_trend) '-' num2str(ende_trend) '.nc'],'trend_values');
olens_cesm_lens  = olens;
for m = 2:length(models)
  eval(['olens_' models{m} '   = olens_cesm_lens-mean(cesm_lens)+mean(' models{m} ');'])
end






%% --- PLOTTING ---------------------------------------------------------------
close all

cols = [0 192 0;
        255 160 16;
        255 224 32;
        255 0 0;
        80 208 255;
        0 32 255;
        160 32 255;...
        0 0 0]/255;

xlim = [-.25 2.75];

if start_trend == 1951
  xlim1 = [-.2 3];
  xlim2 = xlim1;
elseif start_trend == 1965
  xlim1 = [.3 3.5];
  xlim2 = xlim1;
end

figure3 = figure;
set(figure3, 'units', 'centimeters', 'pos', [10 10 20 9])

hold on
title(['North America annual land temperature trend ' num2str(start_trend) '-' num2str(ende_trend)])

bw = .18; % band width of kernel smoother; set to 0 if want to use default

% -- CMIP5
if bw == 0
  [f0,x0,u] = ksdensity(cmip5_all(:));
else
  [f0,x0,u] = ksdensity(cmip5_all(:),'bandwidth',bw);
end
h0 = jbfill(x0, f0*.25*length(time_trend), x0*0, [.8 .8 .8],'none');
hold on
% -- all models but CESM and MPI
for m = 2:length(models)-1
  % -- kernel
  tmp0 = eval([models{m}]);
  if bw == 0
    [f0,x0,u] = ksdensity(tmp0);
  else
    [f0,x0,u] = ksdensity(tmp0,'bandwidth',bw)%,'support',[min(tmp0)-.1 max(tmp0)+.1])
  end
  h3(m) = line(x0, f0*.25*length(time_trend), 'Color',[.2 .2 .2], 'LineWidth',1);
end
% -- CESM and MPI
for m = [1 7]
  % -- kernel
  tmp0 = eval([models{m}]);
  if bw == 0
    [f0,x0,u] = ksdensity(tmp0);
  else
    [f0,x0,u] = ksdensity(tmp0,'bandwidth',bw)%,'support',[min(tmp0)-.1 max(tmp0)+.1])
  end
  y = f0*.25*length(time_trend);
  h3(m) = line(x0, y, 'Color',cols(m,:), 'LineWidth',2);
end
set(gca,'Layer','top','YTickLabel',[],'YTick',[],'XLim',xlim1,'YLim',[0 30])
% -- Obs
h1 = vline(obs,'k')
set(h1,'LineWidth',3)
vline(0,'k--')
xlabel(['Trend (\circ C ' num2str(length(time_trend)) '^{-1} years)'])
ylabel('Relative density')
box on
legend([h3(1) h3(7) h3(2) h0(1) h1],[ model_names_legend(1) model_names_legend(7) 'other LEs' 'CMIP5' 'Observations'],'Location','Northeast','FontSize',8)
legend boxoff


% -- save figure (adjust as needed)
% set(gcf,'PaperPositionMode','auto');
% fileo = ['~/Dropbox/publication/clivar19_perspective/clivar_lens_pdf_plots_' vari '_' seas '_OLENS_1panel_with_CMIP5_' num2str(start_trend) '-' num2str(ende_trend)];
% print('-r300','-loose', '-depsc', ['' fileo '.eps'])
% save2pdf(['' fileo '.pdf'])
% saveas(gcf,fileo,'jpg')
