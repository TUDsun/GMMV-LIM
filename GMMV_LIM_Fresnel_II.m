clearvars; restoredefaultpath; add_path

%% load data
% TM
str = 'FoamDielIntTM'; fre = 2 : 8; para.pt = PT.TM; NTX = 8; centre = [0, 0, 0]; radius = 0.06;
% str = 'FoamDielExtTM'; fre = 2 : 8; para.pt = PT.TM; NTX = 8; centre = [0 0, 0]; radius = 0.06;
% str = 'FoamTwinDielTM'; fre = 2 : 8; para.pt = PT.TM; NTX = 18; centre = [0 0, 0]; radius = 0.06;
% str = 'FoamMetExtTM'; fre = 2 : 8; para.pt = PT.TM; NTX = 18; centre = [0 0, 0]; radius = 0.06;
% TE
% str = 'FoamDielIntTE'; fre = 2 : 8; para.pt = PT.TE; NTX = 8; centre = [, 0, 0]; radius = 0.06;
% str = 'FoamDielExtTE'; fre = 2 : 8; para.pt = PT.TE; NTX = 8; centre = [0 0, 0]; radius = 0.06;
% str = 'FoamTwinDielTE'; fre = 2 : 8; para.pt = PT.TE; NTX = 18; centre = [0 0, 0]; radius = 0.06;
% str = 'FoamMetExtTE'; fre = 2 : 8; para.pt = PT.TE; NTX = 18; centre = [0 0, 0]; radius = 0.06;

xs              = 100; 
ys              = 100; 
para.m_unit     = 2.0e-3;
regSize         = [-3.0, 3.0; -3.0, 3.0]; % 12GHz
TxInterval      = 360 / NTX;
para.NTX        = NTX;
interval        = NTX;
rho             = 0.20066;
DTr             = rho / 0.6964;
% DTr             = zeros(1, 18);

% if para.pt == PT.TE; DTr(2) = 0; end

para.Tr      	= 1.67 + DTr / 2;
para.Rr       	= 1.67 + DTr / 2;
rawdat          = load([str '.txt']);
rawdat(:, 1)    = TxInterval * (rawdat(:, 1) - 1) + 0.5;

CVr             = 0.1;
% CV              = sort(randsample(NRX, round(CVr * NRX)));
CV              = 24 : 24 : 241; 
CV(end)         = [];
CV              = sort([CV, CV + 1, CV + 2, CV + 3]);

para.nK         = 0 * floor(2 * pi * radius ./ (0.3 ./ fre));

%% Configuration

[Emea, dat, Phi, ~, ~, ~, ~, ~, grid3d, pars, pars_LSM] ...
    = Pre_InvFresnel_Conf(rawdat, fre, regSize, centre, xs, ys, para);

%%

[ILSM, time]    = LSM(Emea, Phi, pars, para, pars_LSM);

ILSM            = ILSM / max(ILSM(:));

ILSMdB          = db(ILSM, 'power');

%%
[NRX, NTX]      = size(dat{1});
dat             = cell2mat(dat);
pars.interval   = interval;
pars.NR         = length(pars.Runiq);
pars.pt         = para.pt;
Nfre            = length(fre);


%% PEC Inversion
% opts            = spgSetParms('optTol', 1e-3, 'verbosity', 1); % Turn off the SPGL1 log output
IterJ           = 300;
opts            = spgSetParms('optTol', 1e-3, 'decTol', 1e-4, 'iterations', IterJ); % Turn off the SPGL1 log output

NP              = 1 : 3;
CVc             = setdiff(1 : (NRX / NP(para.pt)), CV);

X               = GMMV_LIM(dat, Phi, CV, CVc, opts, pars);

I             	= vec2scalar(X,        pars.Ninv, NTX, Nfre, para.pt);
IBP          	= vec2scalar(pars.XBP, pars.Ninv, NTX, Nfre, para.pt);

IdB         	= db(I, 'power');
IBPdB        	= db(IBP, 'power');

%% display
dBrange             = 25;
fontsize            = 8;
[Xh, Yv]            = ndgrid(grid3d{1}.l{1 : 2});
yy                  = Xh(:, 1) * grid3d{1}.unitvalue;
xx                  = Yv(1, :).' * grid3d{1}.unitvalue;
yinv                = yy(pars.ny) + para.m_unit / 2;
xinv                = xx(pars.nx) + para.m_unit / 2;
close all
% figure; imagesc(1e3 * yinv, 1e3 * xinv, I, [0 1]);
% xlabel('y / mm'); ylabel('x / mm'); axis equal tight;
% grid on; set(gca, 'layer', 'top'); colormap(flipud(hot)); colorbar; plottools('on')
% axis(1e3 * [invdom(3) invdom(4) invdom(1) invdom(2)])
% set(gca, 'fontsize', fontsize); 
% ax = gca; ax.XMinorGrid = 'on';ax.YMinorGrid = 'on';ax.ZMinorGrid = 'on';
% 
% 
% figure; imagesc(1e3 * yinv, 1e3 * xinv, ILSM, [0 1]);
% xlabel('y / mm'); ylabel('x / mm'); axis equal tight;
% grid on; set(gca, 'layer', 'top'); colormap(flipud(hot)); colorbar; plottools('on')
% axis(1e3 * [invdom(3) invdom(4) invdom(1) invdom(2)])
% set(gca, 'fontsize', fontsize); 
% ax = gca; ax.XMinorGrid = 'on';ax.YMinorGrid = 'on';ax.ZMinorGrid = 'on';


figure; 
subplot(1, 2, 1);
imagesc(1e3 * yinv, 1e3 * xinv, IdB, [-dBrange 0]);
colormap(subplot(1, 2, 1), flipud(hot)); colorbar; 
xlabel('$y$ / mm', 'interpreter', 'latex'); 
ylabel('$x$ / mm', 'interpreter', 'latex'); 
axis equal tight;
grid on; set(gca, 'layer', 'top'); 
% plottools('on')
% axis(1e3 * [invdom(3) invdom(4) invdom(1) invdom(2)])
set(gca, 'fontsize', fontsize); 
ax = gca; ax.XMinorGrid = 'on';ax.YMinorGrid = 'on';ax.ZMinorGrid = 'on';

subplot(1, 2, 2);
imagesc(1e3 * yinv, 1e3 * xinv, ILSMdB, [-dBrange 0]);
colormap(subplot(1, 2, 2), flipud(hot)); colorbar; 
xlabel('$y$ / mm', 'interpreter', 'latex'); 
ylabel('$x$ / mm', 'interpreter', 'latex'); 
axis equal tight;
grid on; set(gca, 'layer', 'top'); 
% plottools('on')
% axis(1e3 * [invdom(3) invdom(4) invdom(1) invdom(2)])
set(gca, 'fontsize', fontsize); 
ax = gca; ax.XMinorGrid = 'on';ax.YMinorGrid = 'on';ax.ZMinorGrid = 'on';

% figure(1); saveTightFigure(gcf, [str 'MMV'])
% figure(2); saveTightFigure(gcf, [str 'LSM'])
% figure(1); saveTightFigure(gcf, [str 'MMVdB'])
% figure(2); saveTightFigure(gcf, [str 'LSMdB'])
% figure; imagesc(1e3 * yinv, 1e3 * xinv, IBP, [0 1]);
% xlabel('y / mm'); ylabel('x / mm'); axis equal tight;
% grid on; set(gca, 'layer', 'top'); colormap(flipud(hot)); colorbar; plottools('on')
% axis(1e3 * [invdom(3) invdom(4) invdom(1) invdom(2)])
% % saveTightFigure(gcf, [str 'LSM'])
% 
% figure; imagesc(1e3 * yinv, 1e3 * xinv, IBPdB, [-dBrange 0]);
% xlabel('y / mm'); ylabel('x / mm'); axis equal tight;
% grid on; set(gca, 'layer', 'top'); colormap(flipud(hot)); colorbar; plottools('on')
% axis(1e3 * [invdom(3) invdom(4) invdom(1) invdom(2)])
% % saveTightFigure(gcf, [str 'LSM'])
time


