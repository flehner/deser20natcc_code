close all
clear all

% ------------------------------------------------------------------------------
% fig3.m
%
% Description: script to create Fig. 3 in Deser et al., 2020,
% https://doi.org/10.1038/s41558-020-0731-2
%
% Note: this is essentially just a plotting script.
%
% The full data is too large to share on github, but is available via:
% - Large Ensembles: http://www.cesm.ucar.edu/projects/community-projects/MMLEA/
%
% Required data:
% - the MMLEA archive
%   (data has been regridded to common 2.5x2.5 degree grid prior,
%   identified by "g025"; this step needs to be done by the user)
%
% Required functions:

%
% Author: Flavio Lehner, May 2020, flehner@ucar.edu
%
% ------------------------------------------------------------------------------


pathin      = '/project/cas/flehner/';
vars        = {'tas','pr'};
comp        = 'day';
seasons     = {'DJF','JJA','annual'};
hist_or_kernel = 2; % histogram or kernel density function

plot1   = 0; % PDF
plot2   = 0; % PDF
plot3   = 0; % PDF
plot4   = 0; % christoph schaer style plot
plot4b  = 1; % same as plot4 but with discrete boxplots
plot4c  = 0; % like plot4, just for visual picking of ensemble member

% -- location to plot:
% -- Dallas
city_name   = 'Dallas, TX';
file_name   = 'dallas'
lat         = 32.7;
lon         = -96.8+360;
% % -- Spokane, WA
% city_name   = 'Spokane, WA';
% file_name   = 'spokane'
% lat         = 47.6;
% lon         = -117.4+360;
% % -- W US
% city_name   = 'Western US';
% file_name   = 'western_us';
% lat         = [32.1 48.5];
% lon         = [-123.4+360 -115.2+360];

% -- month to plot:
% -- 1=Jan, 2=Feb, ...
month = 7; % 7 9
month = 'july'; % july september
ndays = 31; % 31 30
vari  = 'tas'; % tas pr

% -- area --
filein  = [pathin 'cmip5-ng/area_g025.nc'];
area    = ncread([filein],'AREA');

% -- CLIVAR parameters
% -- all
% models      = {'cesm_lens','canesm2_lens','csiro_mk36_lens','gfdl_cm3_lens','gfdl_esm2m_lens','ec_earth_lens','mpi_lens'};
% model_names = {'CESM1-CAM5','CanESM2','CSIRO-Mk3-6-0','GFDL-CM3','GFDL-ESM2M','EC-EARTH','MPI-ESM'};
% ensmem      = [40,50,30,20,30,16,100];
% start0      = [1920,1950,1850,1920,1950,1860,1850];
% ende0       = [2100,2100,2100,2100,2100,2100,2099];
% -- day
models      = {'canesm2_lens','csiro_mk36_lens','gfdl_cm3_lens','gfdl_esm2m_lens','cesm_lens','ec_earth_lens'};
model_names = {'CanESM2','CSIRO-Mk3-6-0','GFDL-CM3','GFDL-ESM2M','CESM1-CAM5','EC-EARTH'};
ensmem0     = [50,30,20,26,40,16];
start0      = [1950,1850,1920,1950,1920,1860];
ende0       = [2100,2100,2100,2100,2100,2100];
for m = 1:length(models)
  model_names_legend(m) = strcat(model_names(m),{' ('},num2str(ensmem0(m)),{')'})
end

% -- CMIP5 parameters
% models_cmip5 = {'ACCESS1-0','ACCESS1-3','bcc-csm1-1-m','bcc-csm1-1','BNU-ESM',...
%                 'CanESM2','CCSM4','CESM1-BGC','CESM1-CAM5','CMCC-CESM','CMCC-CM',...
%                 'CMCC-CMS','CNRM-CM5','CSIRO-Mk3-6-0','EC-EARTH','FGOALS-g2',...
%                 'FIO-ESM','GFDL-CM3','GFDL-ESM2G','GFDL-ESM2M','GISS-E2-H-CC',...
%                 'GISS-E2-H','GISS-E2-R-CC','GISS-E2-R','HadGEM2-AO','HadGEM2-CC',...
%                 'HadGEM2-ES','inmcm4','IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR',...
%                 'MIROC5','MIROC-ESM-CHEM','MIROC-ESM','MPI-ESM-LR','MPI-ESM-MR',...
%                 'MRI-CGCM3','MRI-ESM1','NorESM1-ME','NorESM1-M'};
% start_cmip5  = 1870;
% ende_cmip5   = 2100;

% -- general paramaters
start       = 1950;
ende        = 2099;
time        = start:ende;
refstart1   = 1980; % as in McKinnon et al. (2014, JGR)
refende1    = 1989; % as in McKinnon et al. (2014, JGR)
refstart2   = 2006; % as in McKinnon et al. (2014, JGR)
refende2    = 2015; % as in McKinnon et al. (2014, JGR)
start_trend = 1951;
ende_trend  = 2010;
time_trend  = start_trend:ende_trend;
lati        = [-88.75:2.5:88.75];
loni        = [1.25:2.5:358.75];
lati_obs    = [-88.75:2.5:88.75];
loni_obs    = [0:2.5:357.5];
msd         = [1 32 60 91 121 152 182 213 244 274 305 335]; % month start day
med         = [31 59 90 120 151 181 212 243 273 304 334 365]; % month end day
% -- find grid cells to plot in models
% jj          = find(abs(lati-lat)==min(abs(lati-lat)));
% ii          = find(abs(loni-lon)==min(abs(loni-lon)));
% jj          = jj(1);
% ii          = ii(1);

for x = 1:length(lat)
  tmp0    = find(abs(lati-lat(x))==min(abs(lati-lat(x))));
  tmp(x)  = tmp0(1);
end
ii = tmp(1):tmp(end);
for x = 1:length(lon)
  tmp0    = find(abs(loni-lon(x))==min(abs(loni-lon(x))));
  tmp(x)  = tmp0(1);
end
jj = tmp(1):tmp(end);
npoints = length(ii)*length(jj);
clear('tmp','tmp0')
% -- cut out area --
weights = area(jj,ii);
weights = weights/sum(sum(weights));

% -- load data
if strcmp(vari,'tas')==1
  units       = 'K';
  f           = 1;
elseif strcmp(vari,'pr')==1
  units       = '%';
  f           = 86400;
end

% -- load obs (ERA-interim) --
% -- load full data with month preselected with cdo selmon,7 (--> FASTER) --
start_obs = 1979;
ende_obs  = 2018;
pathin_obs = '/project/cas/flehner/observations/erai/day/';
filein    = ['t2m.day.19790101-20181231.remapbil_g025.' month '.nc']
tmp0      = ncread([pathin_obs filein],'t2m')*f;
for i = 1:length(tmp0)
  tmp(i)     = sum(sum(squeeze(tmp0(jj,ii,i)).*weights));
  % tmp99        = squeeze(tmp0(jj,ii,i));
  % tmp(i,:)     = tmp99(:);
end
tmp3      = squeeze(tmp((refstart1-start_obs)*ndays+1:((refende1-start_obs+1)*ndays)));
ref1      = nanmean(tmp3(:));
obs_gridpoint_ref1 = tmp3(:)-ref1;
tmp3      = squeeze(tmp((refstart2-start_obs)*ndays+1:((refende2-start_obs+1)*ndays)));
obs_gridpoint_ref2 = tmp3(:)-ref1;

% % -- load domain-wide trends
% for m = 1:length(models)
%   m
%   filein  = [pathin models{m} '/' comp '/' vari '/' vari '_' comp '_' model_names{m} '_' seas '_trend_values_' num2str(start_trend) '-' num2str(ende_trend) '.nc'];
%   tmp0    = ncread([filein],'trend_values');
%   eval([models{m} ' = tmp0;']);
% end

% -- load grid point data
for m = 1:length(models) %[1 2 4]
  m
  pathin2   = [pathin models{m} '/' comp '/' vari '/'];
  eval([models{m} '_gridpoint_ref1 = NaN(ensmem0(m),ndays*(refende1-refstart1+1));']);
  eval([models{m} '_gridpoint_ref2 = NaN(ensmem0(m),ndays*(refende2-refstart2+1));']);
  % eval([models{m} '_gridpoint_ref1 = NaN(ensmem0(m),ndays*(refende1-refstart1+1)*npoints);']);
  % eval([models{m} '_gridpoint_ref2 = NaN(ensmem0(m),ndays*(refende2-refstart2+1)*npoints);']);
  for e = 1:ensmem0(m)
    e
    datestr(now,'HH:MM:SS')

    % -- load full data and select grid cell AND month (--> TOO SLOW!!) --
    % filein    = [vari '_' comp '_' model_names{m} '_historical_rcp85_r' num2str(e) 'i1p1_' num2str(start0(m)) '0101-' num2str(ende0(m)) '1231_remapcon_g025.nc']
    % tmp0      = ncread([pathin2 filein],vari)*f;
    % tmp1      = squeeze(tmp0(ii,jj,((refstart1-start0(m))*365)+1:((refende1-start0(m)+1)*365)));
    % for i = 1:(refende1-refstart1+1)
    %   tmp2(i,:) = tmp1((i-1)*365+1:i*365);
    % end
    % tmp3      = tmp2(:,msd(month):med(month));
    % ref1      = nanmean(tmp3(:));
    % eval([models{m} '_gridpoint_ref1(e,:) = tmp3(:)-ref1;']);
    % tmp1      = squeeze(tmp0(ii,jj,((refstart2-start0(m))*365)+1:((refende2-start0(m)+1)*365)));
    % for i = 1:(refende1-refstart1+1)
    %   tmp2(i,:) = tmp1((i-1)*365+1:i*365);
    % end
    % tmp3      = tmp2(:,msd(month):med(month));
    % eval([models{m} '_gridpoint_ref2(e,:) = tmp3(:)-ref1;']);

    % -- load full data with month preselected with CDO (--> FASTER) --
    filein    = [vari '_' comp '_' model_names{m} '_historical_rcp85_r' num2str(e) 'i1p1_' num2str(start0(m)) '0101-' num2str(ende0(m)) '1231_remapcon_g025_' month '.nc']
    tmp0      = ncread([pathin2 filein],vari)*f;
    for i = 1:length(tmp0)
      tmp(i)     = sum(sum(squeeze(tmp0(jj,ii,i)).*weights));
      % tmp99        = squeeze(tmp0(jj,ii,i));
      % tmp(i,:)     = tmp99(:);
    end
    tmp3      = squeeze(tmp((refstart1-start0(m))*ndays+1:((refende1-start0(m)+1)*ndays)));
    % tmp3      = squeeze(tmp((refstart1-start0(m))*ndays+1:((refende1-start0(m)+1)*ndays),:));
    ref1      = nanmean(tmp3(:));
    eval([models{m} '_gridpoint_ref1(e,:) = tmp3(:)-ref1;']);
    tmp3      = squeeze(tmp((refstart2-start0(m))*ndays+1:((refende2-start0(m)+1)*ndays)));
    % tmp3      = squeeze(tmp((refstart2-start0(m))*ndays+1:((refende2-start0(m)+1)*ndays),:));
    eval([models{m} '_gridpoint_ref2(e,:) = tmp3(:)-ref1;']);
    tmp4      = squeeze(tmp((start-start0(m))*ndays+1:((ende-start0(m)+1)*ndays)));
    eval([models{m} '_gridpoint(e,:) = squeeze(tmp4)-ref1;']);
    eval([models{m} '_gridpoint_abs(e,:) = tmp4-273.15;']);
  end
end
clear('tmp1','tmp2','tmp3','tmp4')








%% --- PLOTTING ---------------------------------------------------------------
close all

cols_no = [1 2 3 4 5 6]; % to be consistent with Fig 1, 2, and 4

cols = [255 0 0;
        255 160 16;
        255 224 32;
        0 192 0;
        80 208 255;
        0 32 255;
        160 32 255]/255;
% cols = [0 0 0;
%         255 160 16;
%         255 224 32;
%         0 192 0;
%         80 208 255;
%         0 32 255;
%         160 32 255]/255;
cols_light1 = (1-cols)*.75+cols;

binranges = [-12:0.1:12];
xlim = [binranges(1) binranges(end)];
if length(lat) == 1
  ylim = [0 .3];
else
  ylim = [0 .5];
end
% bw        = 0.1;

tmp1 = obs_gridpoint_ref1(:);
tmp2 = obs_gridpoint_ref2(:);
[obs_f0,obs_x0,u] = ksdensity(obs_gridpoint_ref1(:),binranges);
[obs_f1,obs_x1,u] = ksdensity(obs_gridpoint_ref2(:),binranges);


if plot1 == 1
% ---------
close all
figure1 = figure;
set(figure1, 'units', 'centimeters', 'pos', [10 10 8 13])

mc = 1;
for m = [1 2 4]

  tmp1 = eval([models{m} '_gridpoint_ref1(:);']);
  tmp2 = eval([models{m} '_gridpoint_ref2(:);']);
  [f0,x0,u] = ksdensity(tmp1,binranges);
  [f1,x1,u] = ksdensity(tmp2,binranges);
  % [f0,x0,u] = ksdensity(tmp1,binranges,'Bandwidth',bw);
  % [f1,x1,u] = ksdensity(tmp2,binranges,'Bandwidth',bw);

  subplot(3,1,mc)
  hold on
  % title([{'July daily temperatures,'},{'grid cell centered on Boulder, CO'}])
  title([model_names_legend{m}],'Interpreter','none')
  [x_1,y0_1,y1_1] = pdf_shading_low(tmp1,5,x1,f1);
  [x_2,y0_2,y1_2] = pdf_shading_low(tmp1,10,x1,f1);
  [x_3,y0_3,y1_3] = pdf_shading_high(tmp1,90,x1,f1);
  [x_4,y0_4,y1_4] = pdf_shading_high(tmp1,95,x1,f1);
  h02 = patch([x_2 fliplr(x_2)],[y0_2 fliplr(y1_2)],cols(2,:),'Edgecolor','none');
  h01 = patch([x_1 fliplr(x_1)],[y0_1 fliplr(y1_1)],cols(1,:),'Edgecolor','none');
  h03 = patch([x_3 fliplr(x_3)],[y0_3 fliplr(y1_3)],cols(3,:),'Edgecolor','none');
  h04 = patch([x_4 fliplr(x_4)],[y0_4 fliplr(y1_4)],cols(4,:),'Edgecolor','none');
  h1 = line(x0, f0, 'Color',[.5 .5 .5], 'LineWidth',2);
  h2 = line(x1, f1, 'Color',[0 0 0], 'LineWidth',2);
  set(gca,'XLim',xlim,'YTickLabel',[])
  m1 = median(tmp1);
  m2 = median(tmp2);
  yl = ylim;
  plot([m1 m1],[0 yl(2)],'Color',[.5 .5 .5],'LineStyle','--','LineWidth',2)
  plot([m2 m2],[0 yl(2)],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
  box on
  if mc == 3
    xlabel('Daily temperature anomalies (degC)')
  end
  if mc == 1
    legend([h1 h2],[num2str(refstart1) '-' num2str(refende1)],[num2str(refstart2) '-' num2str(refende2)],'Location','NorthWest')
  end
  if mc == 3
    legend([h01 h02 h03 h04],['5th %-ile'],['10th %-ile'],['90th %-ile'],['95th %-ile'],'Location','NorthWest')
  end

  mc = mc + 1;
  clear('tmp1','tmp2','tmp3x','tmp3y','tmp4x','tmp4y')

end
% return
set(gcf,'PaperPositionMode','auto');
fileo = ['/home/flehner/publication/clivar19_perspective/clivar_lens_pdf_extremes_' vari '_mckinnon'];
% print('-r300','-loose', '-depsc', ['' fileo '.eps'])
% save2pdf(['' fileo '.pdf'])
% saveas(gcf,fileo,'jpg')
return
end




if plot2 == 1
% ---------
close all
figure2 = figure;
set(figure2, 'units', 'centimeters', 'pos', [10 10 8 13])

mc = 1;
for m = [1 2 4]

  tmp1 = eval([models{m} '_gridpoint_ref1(:);']);
  tmp2 = eval([models{m} '_gridpoint_ref2(:);']);
  [f0,x0,u] = ksdensity(tmp1,binranges);
  [f1,x1,u] = ksdensity(tmp2,binranges);
  % [f0,x0,u] = ksdensity(tmp1,binranges,'Bandwidth',bw);
  % [f1,x1,u] = ksdensity(tmp2,binranges,'Bandwidth',bw);

  subplot(3,1,mc)
  hold on
  % title([{'July daily temperatures,'},{'grid cell centered on Boulder, CO'}])
  % title([models{m}],'Interpreter','none')
  title([model_names_legend{m}],'Interpreter','none')
  for e = 1:ensmem0(m)
    [f0,x0,u] = eval(['ksdensity(' models{m} '_gridpoint_ref1(e,:),binranges);']);  %,binranges);
    tmp3y(e,:) = f0;
    tmp3x(e,:) = x0;
    [f0,x0,u] = eval(['ksdensity(' models{m} '_gridpoint_ref2(e,:),binranges);']);  %,binranges);
    tmp4y(e,:) = f0;
    tmp4x(e,:) = x0;
  end
  line(tmp3x', tmp3y', 'Color',cols(2,:), 'LineWidth',1);
  line(tmp4x', tmp4y', 'Color',cols(3,:), 'LineWidth',1);
  h2 = jbfill(binranges,prctile(tmp4y,95),prctile(tmp4y,5),cols(3,:),'none',1,.5)
  hold on
  h1 = jbfill(binranges,prctile(tmp3y,95),prctile(tmp3y,5),cols(2,:),'none',1,.5)
  hold on
  vline(0,'k--')
  set(gca,'XLim',xlim)
  set(gca,'XLim',xlim,'YTickLabel',[])
  box on
  if mc == 3
    xlabel('Daily temperature anomalies (degC)')
  end
  if mc == 1
    legend([h1 h2],[num2str(refstart1) '-' num2str(refende1)],[num2str(refstart2) '-' num2str(refende2)],'Location','NorthWest')
  end

  mc = mc + 1;
  clear('tmp1','tmp2','tmp3x','tmp3y','tmp4x','tmp4y')

end
% return
set(gcf,'PaperPositionMode','auto');
fileo = ['/home/flehner/publication/clivar19_perspective/clivar_lens_pdf_extremes_' vari '_fischer'];
print('-r300','-loose', '-depsc', ['' fileo '.eps'])
% save2pdf(['' fileo '.pdf'])
saveas(gcf,fileo,'jpg')
return
end




if plot3 == 1
% ------
close all

figure3 = figure;
set(figure3, 'units', 'centimeters', 'pos', [10 10 10 17])

mc = 1;
for m = [1 2 4]

  tmp1 = eval([models{m} '_gridpoint_ref1(:);']);
  tmp2 = eval([models{m} '_gridpoint_ref2(:);']);
  [f0,x0,u] = ksdensity(tmp1,binranges);
  [f1,x1,u] = ksdensity(tmp2,binranges);
  % [f0,x0,u] = ksdensity(tmp1,binranges,'Bandwidth',bw);
  % [f1,x1,u] = ksdensity(tmp2,binranges,'Bandwidth',bw);

  subplot(3,1,mc)
  hold on
  if mc == 1
    title([{'July daily temperatures,'},{city_name}])
  end
  ylabel([model_names_legend{m}],'Interpreter','none')
  for e = 1:ensmem0(m)
    [f0,x0,u] = eval(['ksdensity(' models{m} '_gridpoint_ref1(e,:),binranges);']);  %,binranges);
    tmp3y(e,:) = f0;
    tmp3x(e,:) = x0;
    [f0,x0,u] = eval(['ksdensity(' models{m} '_gridpoint_ref2(e,:),binranges);']);  %,binranges);
    tmp4y(e,:) = f0;
    tmp4x(e,:) = x0;
    % -- find ensemble members with highest/lowest change in 95th percentile:
    p1 = eval(['prctile(' models{m} '_gridpoint_ref1(e,:),95);']);
    p2 = eval(['prctile(' models{m} '_gridpoint_ref2(e,:),95);']);
    p(e) = p2-p1;
  end
  h2 = jbfill(binranges,prctile(tmp4y,95),prctile(tmp4y,5),cols(3,:),'none',1,.5)
  hold on
  h1 = jbfill(binranges,prctile(tmp3y,95),prctile(tmp3y,5),cols(2,:),'none',1,.5)
  hold on
  [vmax,imax] = max(p);
  [vmin,imin] = min(p);
  h3 = line(tmp3x(imax,:)', tmp3y(imax,:)','Color','b')% 'Color',cols(2,:));
  line(tmp4x(imax,:)', tmp4y(imax,:)','Color','r')% 'Color',cols(3,:));
  h4 = line(tmp3x(imin,:)', tmp3y(imin,:)','Color','b','LineStyle','--')% 'Color',cols(2,:));
  line(tmp4x(imin,:)', tmp4y(imin,:)','Color','r','LineStyle','--')%'Color',cols(3,:));
  % line(mean(tmp3x)', mean(tmp3y)', 'Color','b', 'LineWidth',1);
  % line(mean(tmp4x)', mean(tmp4y)', 'Color','r', 'LineWidth',1);
  h5 = line(obs_x0, obs_f0, 'Color','b', 'LineWidth',2);
  h5 = line(obs_x1, obs_f1, 'Color','r', 'LineWidth',2);
  set(gca,'XLim',xlim,'YLim',ylim,'YTickLabel',[],'YTick',[])
  if mc == 3
    xlabel('Daily temperature anomalies (degC)')
  end
  if mc == 1
    % lgd = legend([h1 h2 h3 h4 h5],[num2str(refstart1) '-' num2str(refende1) char(10) ' 5-95% range'],...
    % [num2str(refstart2) '-' num2str(refende2) char(10) ' 5-95% range'],...
    % ['Ens. mem. with' char(10) 'largest diff. in' char(10) '95th %-ile'],...
    % ['Ens. mem. with' char(10) 'smallest diff. in' char(10) '95th %-ile'],...
    % ['Observations (ERA-interim)'],...
    % 'Location','NorthWest')
    % lgd.FontSize = 6;
    % legend boxoff
    lgd = legend([h1 h2],[num2str(refstart1) '-' num2str(refende1) char(10) ' 5-95% range'],...
    [num2str(refstart2) '-' num2str(refende2) char(10) ' 5-95% range'],...
    'Location','NorthWest')
    lgd.FontSize = 6;
    legend boxoff
  end
  if mc == 2
    lgd = legend([h3 h4 h5],...
    ['Ens. mem. with' char(10) 'largest diff. in' char(10) '95th %-ile'],...
    ['Ens. mem. with' char(10) 'smallest diff. in' char(10) '95th %-ile'],...
    ['Observations' char(10) '(ERA-interim)'],...
    'Location','NorthWest')
    lgd.FontSize = 6;
    legend boxoff
  end
  vline(0,'k--')
  box on

  mc = mc + 1;
  clear('p','tmp1','tmp2','tmp3x','tmp3y','tmp4x','tmp4y')

end
% return
set(gcf,'PaperPositionMode','auto');
fileo = ['/home/flehner/publication/clivar19_perspective/clivar_lens_pdf_extremes_' file_name '_' vari '_fischer'];
print('-r300','-loose', '-depsc', ['' fileo '.eps'])
saveas(gcf,fileo,'jpg')
% save2pdf(['' fileo '.pdf'])
return
end





if plot4 == 1
% ------
close all

figure4 = figure;
set(figure4, 'units', 'centimeters', 'pos', [10 10 20 20])

xlim = [start ende];

time = NaN((ende-start+1)*ndays,1);
for i = 1:length(time)/ndays
  for d = 1:ndays
    time((i-1)*ndays+d) = start+i-1+(d/ndays);
  end
end

pp = 99.9; %99.9 99
em = 10;
wl = 10; % window length years
for m = 1:length(models) % [1 2 4]
  for e = 1:ensmem0(m)
    tmp0 = eval([models{m} '_gridpoint(e,1:50*ndays);']);
    th(e) = prctile(tmp0,pp);
  end
  th = nanmean(th);
  for e = 1:ensmem0(m)
    tmp1 = eval([models{m} '_gridpoint(e,:);']);
    events_tmp = tmp1>th;
    if m == 1 && e < em+1 % plot vertical lines just for first LE
      h(e) = subplot(em+3,1,e);
      plot([time(tmp1>th) time(tmp1>th)],[0 1],'Color',[.5 .5 .5])%cols_light1(1,:))
      pos = get(h(e),'Position');
      % pos(2) = 0.5 ;                         % horizontal pos
      pos(4) = 0.05 ;                        % height
      set( h(e), 'Position', pos ) ;
      set(gca,'Layer','top','XLim',xlim,'XTick',[],'YTick',[],'YTickLabel',[],'XTickLabel',[])
      ylabel(['CESM' char(10) '#' num2str(e)],'FontSize',8)
      % if e < em
      %   set(gca,'XLim',xlim,'XTick',[],'YTick',[],'YTickLabel',[],'XTickLabel',[])
      % else
      %   set(gca,'XLim',xlim,'YTick',[],'YTickLabel',[])
      %   xlabel('Time (Year)')
      % end
    end
    for i = 1:(length(time)/ndays)-wl
      freq_runsum1(m,e,i) = sum(events_tmp((i-1)*ndays+1:(i+wl-1)*ndays))/(wl*ndays);
    end
  end
  tmp = eval([models{m} '_gridpoint(:,:);']);
  events = tmp>th;
  freq(m,:) = sum(events,1)/ensmem0(m);
  tmpsum = sum(events,1);
  for i = 1:(length(events)/ndays)-wl
    freq_runsum2(m,i) = sum(tmpsum((i-1)*ndays+1:(i+wl-1)*ndays))/(ensmem0(m)*wl*ndays);
  end
end


he = subplot(em+3,1,[em+1 em+3])
hold on
title(['Probability to exceed ' num2str(pp) 'th percentile (of 1950-1999 ' month ' daily ' vari ')'])
pos = get(he,'Position');
pos(4) = 0.15 ;                        % height
set(he,'Position',pos);
% % plot(time,freq)
% % plot(time,rm(freq(1,:),wl),'Color',cols_light1(1,:))
% plot(start+wl/2:ende-wl/2,squeeze(freq_runsum1(1,:,:)),'Color',cols_light1(1,:))
for m = 1:length(models)
  jbfill(start+wl/2:ende-wl/2,prctile(squeeze(freq_runsum1(m,1:ensmem0(m),:)),95),prctile(squeeze(freq_runsum1(m,1:ensmem0(m),:)),5),cols_light1(cols_no(m),:),'none',1,.5)
  % plot(start+wl/2:ende-wl/2,prctile(squeeze(freq_runsum1(m,1:ensmem0(m),:)),95),'Color',cols_light1(cols_no(m),:),'LineWidth',.5)
  % plot(start+wl/2:ende-wl/2,prctile(squeeze(freq_runsum1(m,1:ensmem0(m),:)),5),'Color',cols_light1(cols_no(m),:),'LineWidth',.5)
end
hold on
for m = 1:length(models)
  h(m) = plot(start+wl/2:ende-wl/2,freq_runsum2(m,:),'Color',cols(cols_no(m),:),'LineWidth',1.5)
end
box on
set(gca,'Layer','top')
legend([h(2) h(3) h(4) h(1) h(5)],model_names([2 3 4 1 5]),'Location','NorthWest','Interpreter','none','FontSize',8)
legend boxoff
ylabel('Probability')
xlabel('Time (Year)')


% return
set(gcf,'PaperPositionMode','auto');
fileo = ['/home/flehner/publication/clivar19_perspective/clivar_lens_pdf_extremes_' file_name '_' vari '_schaer'];
print('-r300','-loose', '-depsc', ['' fileo '.eps'])
saveas(gcf,fileo,'jpg')
% save2pdf(['' fileo '.pdf'])
return
end







if plot4b == 1
% ------
close all

figure4 = figure;
set(figure4, 'units', 'centimeters', 'pos', [10 10 20 20])
% set(figure4, 'units', 'centimeters', 'pos', [10 10 20 12])

xlim = [start ende];

time = NaN((ende-start+1)*ndays,1);
for i = 1:length(time)/ndays
  for d = 1:ndays
    time((i-1)*ndays+d) = start+i-1+(d/ndays);
  end
end

pp = 99.9; %99.9 99
em = 2;%10;
wl = 10; % window length years
for m = 1:length(models) % [1 2 4]
  m
  for e = 1:ensmem0(m)
    tmp0 = eval([models{m} '_gridpoint_abs(e,1:50*ndays);']);
    th = prctile(tmp0,pp);
  end
  th = nanmean(th);
  for e = 1:ensmem0(m)
    tmp1 = eval([models{m} '_gridpoint_abs(e,:);']);
    events_tmp = tmp1>th;
    for i = 1:(length(time)/ndays)-wl
      freq_runsum1(m,e,i) = ( sum(events_tmp((i-1)*ndays+1:(i+wl-1)*ndays))/(wl*ndays) )*100;
    end
  end
  tmp = eval([models{m} '_gridpoint_abs(:,:);']);
  events = tmp>th;
  freq(m,:) = sum(events,1)/ensmem0(m);
  tmpsum = sum(events,1);
  for i = 1:(length(events)/ndays)-wl
    freq_runsum2(m,i) = ( sum(tmpsum((i-1)*ndays+1:(i+wl-1)*ndays))/(ensmem0(m)*wl*ndays) )*100;
  end
end

% -- "bar code" -- plot just for CESM
m = 4;
for e = 1:ensmem0(m)
  tmp0 = eval([models{m} '_gridpoint_abs(e,1:50*ndays);']);
  th(e) = prctile(tmp0,pp);
end
th = nanmean(th);
for e = 1:ensmem0(m)
  tmp1 = eval([models{m} '_gridpoint_abs(e,:);']);
  event_count(e) = sum(tmp1>th);
end
e_sel(1) = find(event_count==min(event_count));
e_sel(2) = find(event_count==max(event_count));
e_sel(3:5) = [4,8,17];
for e = 1:length(e_sel)
  tmp99 = eval([models{m} '_gridpoint_abs(e_sel(e),:);']);
  h(e) = subplot(8,1,e);
  hold on
  if e == 1
    title(['(a) Daily heat extreme occurrences in July for ' city_name],'FontSize',8)
  end
  hl = plot([time(tmp99>th) time(tmp99>th)],[20 40],'Color',[.5 .5 .5])%cols_light1(1,:))
  % hl = vline(time(tmp99>th),'k')
  set(hl,'Color',[.5 .5 .5])
  legend([hl(1)],['no of events = ' num2str(event_count(e_sel(e)))],'Location','NorthWest')
  legend boxoff
  set(gca,'Layer','top','XLim',xlim,'XTick',[],'YTick',[],'YTickLabel',[],'XTickLabel',[])
  ylabel(['CESM' char(10) '#' num2str(e_sel(e))],'FontSize',8)
  box on
end

% -- shading --
he = subplot(8,1,[6 8])
hold on
title(['(b) Probability to exceed historical (1950-1999) daily ' num2str(pp) 'th percentile'],'FontSize',8)
for m = 1:length(models)
  jbfill(start+wl/2+1:ende-wl/2+1,prctile(squeeze(freq_runsum1(m,1:ensmem0(m),:)),95),prctile(squeeze(freq_runsum1(m,1:ensmem0(m),:)),5),cols_light1(cols_no(m),:),'none',1,.5)
  % plot(start+wl/2:ende-wl/2,prctile(squeeze(freq_runsum1(m,1:ensmem0(m),:)),95),'Color',cols_light1(cols_no(m),:),'LineWidth',.5)
  % plot(start+wl/2:ende-wl/2,prctile(squeeze(freq_runsum1(m,1:ensmem0(m),:)),5),'Color',cols_light1(cols_no(m),:),'LineWidth',.5)
end
hold on
% -- lines --
for m = 1:length(models)
  h(m) = plot(start+wl/2+1:ende-wl/2+1,freq_runsum2(m,:),'Color',cols(cols_no(m),:),'LineWidth',1.5)
end
% -- box plots
offset = 0.6;
ds = [1970 1990 2010 2030 2050 2070 2090];
for d = ds
  for m = 1:length(models)
    xpos = d+wl/2 - 3*offset + m*offset;  % d+wl/2-1+offset*m-2.5*offset
    plot([xpos xpos],[prctile(squeeze(freq_runsum1(m,1:ensmem0(m),d-start)),5) prctile(squeeze(freq_runsum1(m,1:ensmem0(m),d-start)),95)],'Color',cols(cols_no(m),:),'LineWidth',1.5)
    plot([xpos xpos],[nanmean(squeeze(freq_runsum1(m,1:ensmem0(m),d-start))) nanmean(squeeze(freq_runsum1(m,1:ensmem0(m),d-start)))],'o','Color',cols(cols_no(m),:),'LineWidth',1.5,'MarkerSize',5)
  end
end
box on
set(gca,'Layer','top','XLim',xlim,'XTick',ds)
legend([h(1) h(2) h(3) h(4) h(5) h(6)],model_names,'Location','NorthWest','Interpreter','none','FontSize',8)
legend boxoff
ylabel('Probability (%)')
xlabel('Time (Year)')


return
set(gcf,'PaperPositionMode','auto');
fileo = ['/home/flehner/publication/clivar19_perspective/clivar_lens_pdf_extremes_' file_name '_' vari '_schaer'];
print('-r300','-loose', '-depsc', ['' fileo '.eps'])
% print('-r300','-loose', '-painters', ['' fileo '.eps'])
saveas(gcf,fileo,'jpg')
save2pdf(['' fileo '.pdf'])
return
end






if plot4c == 1
% ------
close all

figure4 = figure;
set(figure4, 'units', 'centimeters', 'pos', [10 10 40 20])

xlim = [start ende];

time = NaN((ende-start+1)*ndays,1);
for i = 1:length(time)/ndays
  for d = 1:ndays
    time((i-1)*ndays+d) = start+i-1+(d/ndays);
  end
end

pp = 99.9; %99.9 99
em = 40;
wl = 10; % window length years
for m = 4:4%length(models) % [1 2 4]
  for e = 1:ensmem0(m)
    tmp0 = eval([models{m} '_gridpoint_abs(e,1:50*ndays);']);
    th(e) = prctile(tmp0,pp);
  end
  th = nanmean(th);
  for e = 1:ensmem0(m)
    e
    tmp1 = eval([models{m} '_gridpoint_abs(e,:);']);
    event_count(e) = sum(tmp1>th);
    h(e) = subplot(em/2,2,e);
    plot([time(tmp1>th) time(tmp1>th)],[0 1],'Color',[.5 .5 .5])%cols_light1(1,:))
    % hl = vline(time(tmp1>th),'r');
    % set(hl,'Color',[.5 .5 .5])
    set(gca,'Layer','top','XLim',xlim,'XTick',[],'YTick',[],'YTickLabel',[],'XTickLabel',[])
    ylabel(['CESM' char(10) '#' num2str(e)],'FontSize',8)
    box on
  end
end
event_count
find(event_count==min(event_count))
find(event_count==max(event_count))


end
