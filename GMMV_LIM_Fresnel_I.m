clearvars; restoredefaultpath; add_path

%% Load data, incident fields, measurement matrices, stiffness matrices

% str = 'rectTM_cent'; ptype = PT.TM; centre = [0, 0, 0]; xs = 40; ys = 40; m_unit = 1.0e-3; fre = 4 : 4 : 16; radius = 27.6e-3;

% str = 'rectTM_dece'; ptype = PT.TM; centre = [0, 40, 0]; xs = 40; ys = 40; m_unit = 1.0e-3; fre = 4 : 4 : 16; radius = 27.6e-3;

str = 'uTM_shaped'; para.pt = PT.TM; centre = [0, 0, 0]; xs = 70; ys = 70; para.m_unit = 1.5e-3; fre = 4 : 4 : 16; radius = 50e-3;

% str = 'rectTE_8f'; ptype = PT.TE; centre = [0, 0, 0]; xs = 40; ys = 40; m_unit = 1.0e-3; fre = 4 : 4 : 16; radius = 27.6e-3;

% str = 'dielTM_dec8f'; ptype = PT.TM;

% str = 'twodielTM_8f'; ptype = PT.TM; centre = [0, 0, 0]; xs = 70; ys = 100; m_unit = 1.5e-3; fre = 2 : 2 : 8;

regSize         = [-1.2, 1.2; -1.2, 1.2]; % 12GHz
para.NTX        = 36;
interval        = para.NTX;
% rho             = 0.20066;
% DTr             = rho/0.6964*ones(1,18);
DTr             = 0; % zeros(1, 18);
para.Tr      	= 0.72 + DTr / 2; % para.Tr      	= 0.72 * ones(1, 18) + DTr / 2;
para.Rr       	= 0.76 + DTr / 2;
rawdat          = load([str '.txt']);
TxInterval      = 360 / para.NTX;
rawdat(:, 1)    = TxInterval * (rawdat(:, 1) - 1) - 2.5;
rawdat(:, 2)    = 5 * (rawdat(:, 2) - 1);

NRx             = 49;
CVN             = round(NRx / 5);
CVI             = [CVN CVN + 1];
CV              = [];
for jj = 0 : 3
    CV  = [CV CVI + jj * CVN];
end
CV              = sort(CV);


para.nK         = floor(2 * pi * radius ./ (0.3 ./ fre));

%% Configuration

[Emea, dat, Phi, A, chiinv, vJ, Einc, eTotInv, grid3d, pars, pars_LSM] ...
    = Pre_InvFresnel_Conf(rawdat, fre, regSize, centre, xs, ys, para);

%% LSM

[ILSM, time]    = LSM(Emea, Phi, pars, para, pars_LSM);

% Emea            = cell2mat(Emea);
% [U, S, V]         = svd(Emea);
% Spse            = S';
% Spse            = Spse ./ (Spse.^2 + (0.01 * S(1))^2);
% ILSM            = sum(abs(Spse * U' * cell2mat(Phi)).^2, 1);

% ILSM            = ILSM(:, 1 : NP(ptype) : end) + ILSM(:, NP(ptype) : NP(ptype) : end);
% ILSM            = ILSM / Nfre;

ILSM            = ILSM / max(ILSM(:));

ILSMdB          = db(ILSM, 'power');


%%
[NRX, NTX]      = size(dat{1});
dat             = cell2mat(dat);
pars.interval   = interval;
pars.NR         = length(pars.Runiq);
pars.pt         = para.pt;
Nfre            = length(fre);


%% Inversion
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

%%
% close all
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
subplot(2, 2, 1.5); 
imagesc(1e3 * yinv, 1e3 * xinv, ILSMdB, [-dBrange 0]);
colormap(subplot(2, 2, 1.5), flipud(hot)); colorbar; 
xlabel('$y$ / mm', 'interpreter', 'latex'); 
ylabel('$x$ / mm', 'interpreter', 'latex'); 
axis equal tight; grid on; set(gca, 'layer', 'top'); 
% plottools('on')
% axis(1e3 * [invdom(3) invdom(4) invdom(1) invdom(2)])
set(gca, 'fontsize', fontsize); 
% ax = gca; ax.XMinorGrid = 'on';ax.YMinorGrid = 'on';ax.ZMinorGrid = 'on';

subplot(2, 2, 3);
imagesc(1e3 * yinv, 1e3 * xinv, I);
colormap(subplot(2, 2, 3), jet); colorbar;
xlabel('$y$ / mm', 'interpreter', 'latex'); 
ylabel('$x$ / mm', 'interpreter', 'latex'); 
axis equal tight; grid on; set(gca, 'layer', 'top');  
% plottools('on')
% axis(1e3 * [invdom(3) invdom(4) invdom(1) invdom(2)])
set(gca, 'fontsize', fontsize); 
% ax = gca; ax.XMinorGrid = 'on';ax.YMinorGrid = 'on';ax.ZMinorGrid = 'on';

subplot(2, 2, 4);
imagesc(1e3 * yinv, 1e3 * xinv, IdB, [-dBrange 0]);
colormap(subplot(2, 2, 4), flipud(hot)); colorbar;
xlabel('$y$ / mm', 'interpreter', 'latex'); 
ylabel('$x$ / mm', 'interpreter', 'latex'); 
axis equal tight; grid on; set(gca, 'layer', 'top');  
% plottools('on')
% axis(1e3 * [invdom(3) invdom(4) invdom(1) invdom(2)])
set(gca, 'fontsize', fontsize); 
% ax = gca; ax.XMinorGrid = 'on';ax.YMinorGrid = 'on';ax.ZMinorGrid = 'on';



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



