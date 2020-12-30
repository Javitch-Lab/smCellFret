%--------------------------------------------------------------------------
%
% theorGamma:
%   Calculation of the net emission spectrum and determination of the gamma 
%   and crosstalk correction factors. 
% 
% Description: 
%   The script calculates the net emission spectra of the single molecule 
%   FRET dyes LD555p-SNAP-CLIP and LD655-SNAP-CLIP taking into account the 
%   optical setup of the microscope used. The calculation is performed by
%   multiplying the net optical transmission with the standardized dye 
%   emission spectra. Additionally, the detection efficiencies etaD (donor
%   emission channel) and etaA (acceptor emission channel) are calculated  
%   by integration over the net emission spectra. Using the ratio of the 
%   detection efficiencies etaAD = etaA?etaD and the ratio of the 
%   fluorescence quantum yields phiAD=phiA?phiD the correction factor 
%   gamma = etaAD * phiAD can be calculated [1]. The integrated donor emission  
%   etaC that bleeds into the acceptor channel is a measure of the spectral  
%   donor crosstalk. An estimate of the crosstalk correction factor can be 
%   calculated from the detection efficiency ratio etaC?etaD.
%
%   [1] McCann, J. J., et al. (2010). "Optimizing methods to recover 
%   absolute FRET efficiency from immobilized single molecules." 
%   Biophysical Journal 99(3): 961-970.
%
% Syntax:  
%   theorGamma
% 
% Inputs:
%   code section 
%   1. 'load Data' - Update the path name where all microscope and dye 
%       spectra are located.
%   2. 'Emission Tier 1' - Update the file names of the microscope filter   
%       transmission spectra. 
%   3. 'Quantum Efficiency Curve Evolve 512' - Update the file name of the 
%       camera quantum efficiency curve
%   4. 'Transmittance of NA 1.7 Objective' - Update the file name of the 
%       objective transmission curve. 
%   5. 'Load Spectra' - Update the file names of the dye absorption and 
%       emission spectra.
%   6. 'Calculate Gamma' - Update the value of the fluorescence quantum 
%       yield for each dye. 
%
%   
% Outputs:
%   1. Figure1 - A graph with the filter curves, the camera quantum 
%      efficiency curve and the objective transmission curve. 
%   2. Figure2 - Dye absorption and emission spectra.
%   3. Figure3 - Net transmission spectrum of all optical elements
%   4. Figure4 - Total net emission spectrum, which is a combination of the 
%      net transmission spectrum of the optical elements and the LD  
%      Snap-Clip dye emission spectra ( cf. reference [1]).
%   5. Figure5 - Total net emission spectrum, which is a combination of the 
%      net transmission spectrum of the optical elements and the LD free 
%      dye emission spectra. 
%   6. gammaFree - Workspace variable with the gamma value calculated for 
%      the free LD dyes
%   7. gammaSnap - Workspace variable with the gamma value calculated for 
%      the LD Snap-Clip dyes
%   8. Figure6 - The figure shows the area of the crosstalk for the LD 
%      Snap-Clip dyes
%   9. crtRatioSnap - Workspace variable with the crosstalk value  
%      calculated for the LD Snap-Clip dyes.
%   10. Figure7 - The figure shows the area of the crosstalk for the LD 
%      free dyes
%   11. crtRatioFree - Workspace variable with the crosstalk value  
%      calculated for the LD Free dyes.
%
%
% Author: 
%   - P.G. Jul 2020
%
% Copyright:
%    2011-2020 Columbia Univiversity / RFMH All Rights Reserved.
%
% -------------------------------------------------------------------------

%% Load Data
clear;

% microscope transmission spectra
path = 'F:\Disk E\#BackupMicData\2020 Jan-Apr\200706_PG_Microscope Transmission Spectra';

%% Emission Tier 1
emTier1Raw=readmatrix([path filesep 'emissionTier1.txt']);
% only range from 400nm - 900nm, stepsize 1nm
idx=401:2:1401;
emTier1=emTier1Raw(idx,:);

% Dichroic Tier 1
dcTier1Raw=readmatrix([path filesep 'dichroicTier1.txt']);
idx=21:1:521;
dcTier1=dcTier1Raw(idx,:);

% Emission Green Tier 2
emGreenTier2Raw=readmatrix([path filesep 'emGreenTier2.txt']);
idx=101:1:601;
emGreenTier2=emGreenTier2Raw(idx,:);

% Emission Red Tier 2
emRedTier2Raw=readmatrix([path filesep 'emRedTier2.txt']);
idx=201:2:1201;
emRedTier2=emRedTier2Raw(idx,:);

% Dichroic Tier 2
dcTier2Raw=readmatrix([path filesep 'dichroicTier2.txt']);
idx=21:1:521;
dcTier2 = dcTier2Raw(idx,:);


%% Quantum Efficiency Curve Evolve 512
QE=readmatrix([path filesep 'evolve512-QE.txt']);

% interpolate QE curve over the filter range 
x=QE(:,1);
y=QE(:,2)/100; %normalize to 1
xQE=400:1:900; %range from 400 to 900
yQE=interp1(x,y,xQE,'spline');
% figure;
% scatter(x,y);hold on;
% plot(xq,yq);

%% Transmittance of NA 1.7 Objective
objTrm = readmatrix([path filesep 'Objective NA 1.7 Transmittance']);

% interpolate QE curve over the filter range 
x=objTrm(:,1);
y=objTrm(:,2)/100; %normalize to 1
xObjTrm=400:1:900; %range from 400 to 900
yObjTrm=interp1(x,y,xObjTrm,'spline');

% figure;
% scatter(x,y);hold on;
% plot(xObjTrm,yObjTrm);


%% Load Spectra 

%---- LD Snap-Clip Dyes
% Donor Abs
absSnapLD555p=readmatrix([path filesep 'LD555p labeled SNAP absorption spectrum']);
% Donor Em
emSnapLD555pRaw=readmatrix([path filesep 'LD555p labeled SNAP emission spectrum']);
idx=1:2:341; 
emSnapLD555p=emSnapLD555pRaw(idx,:);

% Acceptor Abs
absSnapLD655=readmatrix([path filesep 'LD655 labeled SNAP absorption spectrum']);
% Acceptor Em
emSnapLD655Raw=readmatrix([path filesep 'LD655 labeled SNAP emission spectrum']);
idx=1:2:341; 
emSnapLD655=emSnapLD655Raw(idx,:);

%---- Free LD Dyes (for comparison)
% Donor Abs
absFreeLD555p=readmatrix([path filesep 'Free LD555p dye absorption spectrum']);
% Donor Em
emFreeLD555pRaw=readmatrix([path filesep 'Free LD555p dye emission spectrum']);
idx=1:2:341; 
emFreeLD555p=emFreeLD555pRaw(idx,:);

% Acceptor Abs
absFreeLD655=readmatrix([path filesep 'Free LD655 dye absorption spectrum']);
% Acceptor Em
emFreeLD655Raw=readmatrix([path filesep 'Free LD655 dye emission spectrum']);
idx=1:2:341; 
emFreeLD655=emFreeLD655Raw(idx,:);


%% plot filter curves 
figure;
wavelength = 400:1:900;
filter = [emTier1(:,2) dcTier1(:,2) emGreenTier2(:,2) emRedTier2(:,2) dcTier2(:,2) yQE' yObjTrm'];
line(wavelength, filter); hold on

%% plot Spectra
figure;

% excitation Waverlength 532nm
x=[532 532];
y=[0 1.01];
line(x,y,'Color', [0 1 0.5]);

%--- LD SNAP-CLIP Dyes
%Spectra Donor
line(absSnapLD555p(:,1),absSnapLD555p(:,2),'LineStyle','--','Color','g');hold on
line(emSnapLD555p(:,1),emSnapLD555p(:,2),'Color','g');

%Spectra Acceptor
line(absSnapLD655(:,1),absSnapLD655(:,2),'LineStyle','--','Color','r');
line(emSnapLD655(:,1),emSnapLD655(:,2),'Color','r');

title('LD SNAP-CLIP Dyes');

%--- LD Dyes
figure;

% excitation Waverlength 532nm
line(x,y,'Color', [0 1 0.5]);

%Spectra Donor
line(absFreeLD555p(:,1),absFreeLD555p(:,2),'LineStyle','--','Color','g');hold on
line(emFreeLD555p(:,1),emFreeLD555p(:,2),'Color','g');

%Spectra Acceptor
line(absFreeLD655(:,1),absFreeLD655(:,2),'LineStyle','--','Color','r');
line(emFreeLD655(:,1),emFreeLD655(:,2),'Color','r');

title('Free LD Dyes');

%% reduced transmission 
dcTier2(1:225,2)=1;
emRedTier2(1:254,2)=1;
emGreenTier2(220:501,2)=1;

trsmFilter = [emTier1(:,2) dcTier1(:,2) emGreenTier2(:,2) dcTier2(:,2) emRedTier2(:,2) yQE' yObjTrm'];
trsmEffFlt = prod(trsmFilter,2);
netFilter = [wavelength' trsmEffFlt];
figure; line(netFilter(:,1), netFilter(:,2)); hold on

%% reduced final emission spectra for LD SNAP-CLIP Dyes

%---------- reduced Donor
idx=131:1:301; % only take the donor part of the net filter 
netFilterDon = netFilter(idx,:);
netFilterDon(128:171,:)=0;
trsmDonor  = [netFilterDon(:,2) emSnapLD555p(:,2)];
trsmEffDon = prod(trsmDonor,2);

% integrate Donor Emission
reducedDonor = [ emSnapLD555p(:,1) trsmEffDon];
%etaDon1=sum(trsmEffDon); the sum is the same because wavelength spacing is 1  
etaDonSnap=trapz(emSnapLD555p(:,1),trsmEffDon);

figure;line(reducedDonor(:,1), reducedDonor(:,2)); hold on
line(emSnapLD555p(:,1), emSnapLD555p(:,2));

%---------- reduced Acceptor
idx=231:1:401; % only take the acceptor part of the net filter 
netFilterAcc = netFilter(idx,:);
trsmAcc  = [netFilterAcc(:,2) emSnapLD655(:,2)];
trsmEffAcc = prod(trsmAcc,2);

% integrate Acceptor
reducedAcc = [ emSnapLD655(:,1) trsmEffAcc];
% etaAcc1=sum(trsmEffAcc)
% etaAcc2=trapz(trsmEffAcc)
etaAccSnap=trapz(emSnapLD655(:,1),trsmEffAcc);

line(reducedAcc(:,1), reducedAcc(:,2));
line(emSnapLD655(:,1), emSnapLD655(:,2));

title('LD SNAP-CLIP Dyes');

%% reduced final emission spectra for LD Free Dyes

%---------- reduced Donor
trsmDonorFree  = [netFilterDon(:,2) emFreeLD555p(:,2)];
trsmEffDonFree = prod(trsmDonorFree,2);

% integrate Donor Emission
reducedDonorFree = [ emFreeLD555p(:,1) trsmEffDonFree];
etaDonFree=trapz(emFreeLD555p(:,1),trsmEffDonFree);

figure;
line(reducedDonorFree(:,1), reducedDonorFree(:,2),'LineStyle','--'); hold on
line(reducedDonor(:,1), reducedDonor(:,2))
line(emSnapLD555p(:,1), emSnapLD555p(:,2));
line(emFreeLD555p(:,1), emFreeLD555p(:,2),'LineStyle','--');

%---------- reduced Acceptor
trsmAccFree  = [netFilterAcc(:,2) emFreeLD655(:,2)];
trsmEffAccFree = prod(trsmAccFree,2);

% integrate Acceptor
reducedAccFree = [ emFreeLD655(:,1) trsmEffAccFree];
etaAccFree=trapz(emFreeLD655(:,1),trsmEffAccFree);

line(reducedAccFree(:,1), reducedAccFree(:,2),'LineStyle','--');
line(reducedAcc(:,1), reducedAcc(:,2));
line(emSnapLD655(:,1), emSnapLD655(:,2));
line(emFreeLD655(:,1), emFreeLD655(:,2),'LineStyle','--');
title('Free LD Dyes');

%% Calculate Gamma
% fluorescence Quantum Yield SnapLD555p and SnapLD655p
phiDonSnap = 0.47;
phiAccSnap = 0.59;
phiRatioSnap = phiAccSnap/phiDonSnap;
% Gamma LD Snap
etaRatioSnap = etaAccSnap/etaDonSnap;
gammaSnap    = etaRatioSnap*phiRatioSnap;

% fluorescence Quantum Yield SnapLD555p and SnapLD655p
phiDonFree = 0.28;
phiAccFree = 0.39;
phiRatioFree = phiAccFree/phiDonFree;
% Gamma Free LD
etaRatioFree = etaAccFree/etaDonFree;
gammaFree    = etaRatioFree*phiRatioFree;

%% Calculate CRT

%--- LD SNAP-CLIP Dyes
figure;
idx=129:1:168;
xSnap = emSnapLD555p(idx,1);
xCrtSnap = xSnap(1):1:740;
ySnap = emSnapLD555p(idx,2);
yCrtSnap = interp1(xSnap,ySnap,xCrtSnap,'linear','extrap');

line(emSnapLD555p(:,1), emSnapLD555p(:,2));hold on
line(xCrtSnap, yCrtSnap,'Color','r');
line(reducedAcc(:,1), reducedAcc(:,2));

crtDonSnap=trapz(xCrtSnap,yCrtSnap);
crtRatioSnap=crtDonSnap/etaDonSnap;
area(xCrtSnap,yCrtSnap)

%% --- LD Free Dyes
figure;
idx=129:1:169;
xFree = emFreeLD555p(idx,1);
xCrtFree = xFree(1):1:736;
yFree = emFreeLD555p(idx,2);
yCrtFree = interp1(xFree,yFree,xCrtFree,'linear','extrap');

line(emFreeLD555p(:,1), emFreeLD555p(:,2));hold on
line(xCrtFree, yCrtFree,'Color','r');
line(reducedAccFree(:,1), reducedAccFree(:,2));

crtDonFree=trapz(xCrtFree,yCrtFree);
crtRatioFree=crtDonFree/etaDonSnap;
area(xCrtFree,yCrtFree)


%%

