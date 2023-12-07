%%%%%% MBOX 4 box ocean model BJW Mills 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Define parameters   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:100;
for j= 1:100;
try

ic = (i-20)/100; %  -0.1 2.3
jc = (j-20)/200; % -0.05 1.05

%%%%%%% set up global structures
global stepnumber
global pars
global workingstate

%%%%%% water reservoir sizes in m3 (m=margins, s=surface, h= hi-lat, d=deep)
pars.vol_p  = 2.6e15 ;  %%%% approx volume of all shelves and slope to depth 100m, area pecentage 5%
pars.vol_di = 5.4e15 ;  %%%% approx volume of all shelves and slope in depth 100-1000m, area pecentage 5%
pars.vol_s  = 2.75e16 ; %%%% approx volume of suface water to depth 100m, area pecentage 76.5%
pars.vol_h  = 1.22e16 ; %%%% approx volume of hi-lat to depth 250m, area pecentage 13.5%
pars.vol_d  = 1.35e18 ;
pars.vol_ocean = pars.vol_p + pars.vol_di + pars.vol_s + pars.vol_h + pars.vol_d;


%%%% mixing coefficient (Sv)
pars.mixcoeff_dip = 30.28;
pars.mixcoeff_ds  = 46.33;
pars.mixcoeff_dh  = 54.9;


%%%%%% initial inorganic carbon reservoirs in moles C
pars.CO2_a_0  = 5e16 ;
pars.DIC_p_0  = 5.2e15 ;
pars.DIC_di_0 = 1.08e16 ;
pars.DIC_s_0  = 5.37e16 ;
pars.DIC_h_0  = 2.71e16 ;
pars.DIC_d_0  = 3e18 ;
pars.ALK_p_0  = 5.2e15 ;
pars.ALK_di_0 = 1.08e16 ;
pars.ALK_s_0  = 5.37e16 ;
pars.ALK_h_0  = 2.71e16 ;
pars.ALK_d_0  = 3e18 ;

%%%%%% initial C isotope composition
pars.R13C12C_PDB   = 0.01124;
pars.d13c_atm_0    = -7 ;
pars.d13c_DIC_p_0  = 0.1 ;
pars.d13c_DIC_di_0 = 0.1 ;
pars.d13c_DIC_s_0  = 0.1 ;
pars.d13c_DIC_h_0  = 0.1 ;
pars.d13c_DIC_d_0  = 0.1 ;
pars.R13C12C_atm_0    = (pars.d13c_atm_0/1000+1)*pars.R13C12C_PDB ;
pars.R13C12C_DIC_p_0  = (pars.d13c_DIC_p_0/1000+1)*pars.R13C12C_PDB ;
pars.R13C12C_DIC_di_0 = (pars.d13c_DIC_di_0/1000+1)*pars.R13C12C_PDB ;
pars.R13C12C_DIC_s_0  = (pars.d13c_DIC_s_0/1000+1)*pars.R13C12C_PDB ;
pars.R13C12C_DIC_h_0  = (pars.d13c_DIC_h_0/1000+1)*pars.R13C12C_PDB ;
pars.R13C12C_DIC_d_0  = (pars.d13c_DIC_d_0/1000+1)*pars.R13C12C_PDB ;



%%%%%% initial POC reservoirs in moles C
pars.POC_p_0  = 630e12 ;
pars.POC_di_0 = 250e12 ;
pars.POC_s_0  = 2329e12 ;
pars.POC_h_0  = 1084e12 ;
pars.POC_d_0  = 56000e12 ;


%%%%%% initial dissolved phosphate in moles
pars.DP_p_0  = 1.82e12 ;
pars.DP_di_0 = 7.56e12 ;
pars.DP_s_0  = 0.55e12 ;
pars.DP_h_0  = 16.27e12 ;
pars.DP_d_0  = 2970e12 ;

%%%%%% initial d13C of POC
pars.d13c_POC_p_0  = -26;
pars.d13c_POC_di_0 = -26;
pars.d13c_POC_s_0  = -26 ;
pars.d13c_POC_h_0  = -26 ;
pars.d13c_POC_d_0  = -26 ;

pars.R13C12C_POC_p_0  = (pars.d13c_POC_p_0/1000+1)*pars.R13C12C_PDB ;
pars.R13C12C_POC_di_0  = (pars.d13c_POC_di_0/1000+1)*pars.R13C12C_PDB ;
pars.R13C12C_POC_s_0  = (pars.d13c_POC_s_0/1000+1)*pars.R13C12C_PDB ;
pars.R13C12C_POC_h_0  = (pars.d13c_POC_h_0/1000+1)*pars.R13C12C_PDB ;
pars.R13C12C_POC_d_0  = (pars.d13c_POC_d_0/1000+1)*pars.R13C12C_PDB ;

%%%%%% initial amount of O2 in moles C
pars.O2_a_0  = 3.7e19 ;
pars.O2_p_0  = 6.705e14 ;
pars.O2_di_0 = 8.964e14 ;
pars.O2_s_0  = 9.139e15 ;
pars.O2_h_0  = 4.02e15 ;
pars.O2_d_0  = 1.823e17 ;


%%%%%% initial amount of FeIII in moles
pars.FeIII_p_0  = 1.56e9 ;
pars.FeIII_di_0 = 3.24e9 ;
pars.FeIII_s_0  = 9.625e9 ;
pars.FeIII_h_0  = 3.66e9 ;
pars.FeIII_d_0  = 810e9 ;

%%%%%% initial amount of sulfate in moles

pars.SO4_p_0  = 7.28e16;
pars.SO4_di_0 = 1.512e17;
pars.SO4_s_0  = 7.7e17;
pars.SO4_h_0  = 3.416e17;
pars.SO4_d_0  = 3.78e19;

%%%%%% initial amount of FeII in moles
pars.FeII_p_0  = 0;
pars.FeII_di_0 = 0;
pars.FeII_s_0  = 0;
pars.FeII_h_0  = 0;
pars.FeII_d_0  = 0;

%%%%%% initial amount of H2S in moles
pars.H2S_p_0  = 0;
pars.H2S_di_0 = 0;
pars.H2S_s_0  = 0;
pars.H2S_h_0  = 0;
pars.H2S_d_0  = 0;

%%%%%% initial size of land S reservoirs in moles
pars.SO4_l_0 = 9.375e19;
pars.pyS_l_0 = 1.875e20;

%%%%%% present day rates, mol yr-1
pars.k_carbw = 12e12 ;
pars.k_sfw = 0 ; %%% Seafloor weathering
pars.k_mccb = 20e12 ;
pars.k_silw = pars.k_mccb - pars.k_carbw ;
basfrac = 0.3 ;
pars.k_granw = pars.k_silw * (1-basfrac) ;
pars.k_basw = pars.k_silw * basfrac ;

%%%%%% organic C cycle, mol yr-1
pars.k_ocdeg = 1.25e12;
pars.k_locb = 2.5e12;
pars.k_mocb = 7e12;
pars.k_oxidw = pars.k_mocb + pars.k_locb - pars.k_ocdeg;

%%%%%% present P, Fe, pyrite and sulfate weathering rate, mol yr-1
pars.k_phosw0 = 0.0967e12;
pars.k_FeIIIw0 = 2e9;
pars.k_pyritew = 1.85e12; %mol S yr-1;
pars.k_sulfatew = 1.25e12; %mol S yr-1;

%%%%%% present pyrite and sulfate burial rate
pars.k_sulfateb   = 1.25e12; %mol S yr-1;

%%%%%% Redfeild ratio
pars.Red_C_P = 106;
pars.Red_C_N = 106/16;
pars.Red_C_O = 106/138;
pars.Red_C_Fe = 106*2000;
pars.Red_Fe_P = pars.Red_C_P/pars.Red_C_Fe;

%%%%%% Monod constant; mol/m3
pars.KP = 0.1e-3;
pars.KFe = 0.1e-6;
pars.KmO2 = 13e-3;
pars.KmFeIII = 10;
pars.KmSO4 = 0.5;

%%%%%% reaction rate constant; m3 mol-1 yr-1
pars.kpy = 0.3708/1e3*24*365.25;
pars.kironO = 1.4e4;
pars.ksulfO = 1e2;
pars.kSironR = 8;
pars.kSide = 4e3;

%%%%%% Ksp
pars.Kspside = 10^(-8.4)*1e6;%%%mol2 m-6
pars.KspFeSaq = 10^(-5.08)*1e3;%%%mol m-3
pars.STFeSaq = 10^(-5.7)*1e3;%%%mol m-3



%%%%%% dust Fe fluxes,mol yr-1
pars.FeIIIa_s0 = 1.443e9; %%%
pars.FeIIIa_h = 0;  %%%% 

%%%% riverine solid FeIII fluxes, mol yr-1
pars.sFeIII_p0 = 190e9;
pars.sFeIII_di0 = 70e9;
pars.sFeIII_d0 =  50e9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Initialise   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% set maximum step size for solver
options = odeset('maxstep',1e5) ;

%%%% run beginning
fprintf('Beginning run: \n')

%%%% set stepnumber to 1
stepnumber = 1 ;

%%%%%%% model timeframe in years (0 = present day)
pars.whenstart = -10e7 ;
pars.whenend = 400e7 ;

%%%% model start state
pars.startstate(1) = pars.CO2_a_0 ;
pars.startstate(2) = pars.DIC_p_0 ;
pars.startstate(3) = pars.DIC_di_0 ;
pars.startstate(4) = pars.DIC_s_0 ;
pars.startstate(5) = pars.DIC_h_0 ;
pars.startstate(6) = pars.DIC_d_0 ;

pars.startstate(7) = pars.ALK_p_0 ;
pars.startstate(8) = pars.ALK_di_0 ;
pars.startstate(9) = pars.ALK_s_0 ;
pars.startstate(10) = pars.ALK_h_0 ;
pars.startstate(11) = pars.ALK_d_0 ;

pars.startstate(12) = pars.CO2_a_0 /(1+1/pars.R13C12C_atm_0);
pars.startstate(13) = pars.DIC_p_0 /(1+1/pars.R13C12C_DIC_p_0);
pars.startstate(14) = pars.DIC_di_0 /(1+1/pars.R13C12C_DIC_di_0);
pars.startstate(15) = pars.DIC_s_0 /(1+1/pars.R13C12C_DIC_s_0);
pars.startstate(16) = pars.DIC_h_0 /(1+1/pars.R13C12C_DIC_h_0);
pars.startstate(17) = pars.DIC_d_0 /(1+1/pars.R13C12C_DIC_d_0);

pars.startstate(18) = pars.POC_p_0
pars.startstate(19) = pars.POC_di_0
pars.startstate(20) = pars.POC_s_0
pars.startstate(21) = pars.POC_h_0
pars.startstate(22) = pars.POC_d_0

pars.startstate(23) = pars.DP_p_0
pars.startstate(24) = pars.DP_di_0
pars.startstate(25) = pars.DP_s_0
pars.startstate(26) = pars.DP_h_0
pars.startstate(27) = pars.DP_d_0

pars.startstate(28) = pars.POC_p_0/(1+1/pars.R13C12C_POC_p_0);
pars.startstate(29) = pars.POC_di_0/(1+1/pars.R13C12C_POC_di_0);
pars.startstate(30) = pars.POC_s_0/(1+1/pars.R13C12C_POC_s_0);
pars.startstate(31) = pars.POC_h_0/(1+1/pars.R13C12C_POC_h_0);
pars.startstate(32) = pars.POC_d_0/(1+1/pars.R13C12C_POC_d_0);

pars.startstate(33) = pars.O2_a_0
pars.startstate(34) = pars.O2_p_0
pars.startstate(35) = pars.O2_di_0
pars.startstate(36) = pars.O2_s_0
pars.startstate(37) = pars.O2_h_0
pars.startstate(38) = pars.O2_d_0

pars.startstate(39) = pars.FeIII_p_0
pars.startstate(40) = pars.FeIII_di_0
pars.startstate(41) = pars.FeIII_s_0
pars.startstate(42) = pars.FeIII_h_0
pars.startstate(43) = pars.FeIII_d_0

pars.startstate(44) = pars.SO4_p_0
pars.startstate(45) = pars.SO4_di_0
pars.startstate(46) = pars.SO4_s_0
pars.startstate(47) = pars.SO4_h_0
pars.startstate(48) = pars.SO4_d_0

pars.startstate(49) = pars.FeII_p_0
pars.startstate(50) = pars.FeII_di_0
pars.startstate(51) = pars.FeII_s_0
pars.startstate(52) = pars.FeII_h_0
pars.startstate(53) = pars.FeII_d_0

pars.startstate(54) = pars.H2S_p_0
pars.startstate(55) = pars.H2S_di_0
pars.startstate(56) = pars.H2S_s_0
pars.startstate(57) = pars.H2S_h_0
pars.startstate(58) = pars.H2S_d_0

pars.startstate(59) = pars.SO4_l_0
pars.startstate(60) = pars.pyS_l_0


%%% Vectors to store the results

statevector =1:1;

state.pO2_af = statevector
state.Atmospheric_CO2_ppmf = statevector

state.DIC_conc_pf = statevector
state.DIC_conc_dif = statevector
state.DIC_conc_sf = statevector
state.DIC_conc_hf = statevector
state.DIC_conc_df = statevector

state.ALK_conc_pf = statevector
state.ALK_conc_dif = statevector
state.ALK_conc_sf = statevector
state.ALK_conc_hf = statevector
state.ALK_conc_df = statevector

state.pH_pf = statevector
state.pH_dif = statevector
state.pH_sf = statevector
state.pH_hf = statevector
state.pH_df = statevector

state.T_sf = statevector
state.T_hf = statevector
state.T_df = statevector
state.T_contf = statevector
state.GASTf = statevector

state.ccdegf = statevector
state.baswf = statevector
state.granwf = statevector
state.silwf = statevector
state.carbwf = statevector

state.mccb_pf = statevector
state.mccb_dif = statevector
state.mccb_df = statevector

state.POC_pf = statevector
state.POC_dif = statevector
state.POC_sf = statevector
state.POC_hf = statevector
state.POC_df = statevector

state.DP_conc_pf = statevector
state.DP_conc_dif = statevector
state.DP_conc_sf = statevector
state.DP_conc_hf = statevector
state.DP_conc_df = statevector

state.O2_conc_pf = statevector
state.O2_conc_dif = statevector
state.O2_conc_sf = statevector
state.O2_conc_hf = statevector
state.O2_conc_df = statevector

state.FeIII_conc_pf = statevector
state.FeIII_conc_dif = statevector
state.FeIII_conc_sf = statevector
state.FeIII_conc_hf = statevector
state.FeIII_conc_df = statevector

state.SO4_conc_pf = statevector
state.SO4_conc_dif = statevector
state.SO4_conc_sf = statevector
state.SO4_conc_hf = statevector
state.SO4_conc_df = statevector

state.FeII_conc_pf = statevector
state.FeII_conc_dif = statevector
state.FeII_conc_sf = statevector
state.FeII_conc_hf = statevector
state.FeII_conc_df = statevector

state.H2S_conc_pf = statevector
state.H2S_conc_dif = statevector
state.H2S_conc_sf = statevector
state.H2S_conc_hf = statevector
state.H2S_conc_df = statevector

state.O2_conc_pf = statevector
state.O2_conc_dif = statevector
state.O2_conc_sf = statevector
state.O2_conc_hf = statevector
state.O2_conc_df = statevector

state.FeIIIwf = statevector

state.FeIIIscavenging_pf = statevector
state.FeIIIscavenging_dif = statevector
state.FeIIIscavenging_sf = statevector
state.FeIIIscavenging_hf = statevector
state.FeIIIscavenging_df = statevector

state.pyritewf = statevector
state.sulfatewf = statevector
%state.pyriteb_pf = statevector
%state.pyriteb_dif = statevector
state.sulfatebf = statevector

state.pyF_pf = statevector
state.pyF_dif = statevector
state.pyF_sf = statevector
state.pyF_hf = statevector
state.pyF_df = statevector

state.ironO_pf = statevector
state.ironO_dif = statevector
state.ironO_sf = statevector
state.ironO_hf = statevector
state.ironO_df = statevector

state.SironR_pf = statevector
state.SironR_dif = statevector
state.SironR_sf = statevector
state.SironR_hf = statevector
state.SironR_df = statevector

state.SideP_pf = statevector
state.SideP_dif = statevector
state.SideP_sf = statevector
state.SideP_hf = statevector
state.SideP_df = statevector

state.pripr_pf = statevector
state.pripr_sf = statevector
state.pripr_hf = statevector

state.remin_pf = statevector
state.remin_dif = statevector
state.remin_sf = statevector
state.remin_hf = statevector
state.remin_df = statevector

state.mocb_pf = statevector
state.mocb_dif = statevector
state.mocb_df = statevector
state.phoswf  = statevector

state.sulfatebentic_pf = statevector
state.sulfatebf = statevector
state.sulfR_pf = statevector
state.sulfO_pf  = statevector
state.SironR_pf = statevector
state.pyF_pf  = statevector
state.pyritewf  = statevector
state.sulfatewf  = statevector
state.SO4_lf  = statevector

state.mocb_p_FeIIIf  = statevector
state.mocb_di_FeIIIf  = statevector
state.mocb_d_FeIIIf  = statevector

state.ocbother_pf  = statevector
state.ocbother_dif  = statevector
state.ocbother_df  = statevector

state.water_sediment_pf  = statevector
state.water_sediment_dif  = statevector
state.water_sediment_df  = statevector

state.BEf_p  = statevector
state.BEf_di  = statevector
state.BEf_d  = statevector

state.AR_pf = statevector;
state.AR_dif = statevector;
state.AR_sf = statevector;
state.AR_hf = statevector;
state.AR_df = statevector;

state.ironR_pf = statevector;
state.ironR_dif = statevector;
state.ironR_sf = statevector;
state.ironR_hf = statevector;
state.ironR_df = statevector;

state.sulfR_pf = statevector;
state.sulfR_dif = statevector;
state.sulfR_sf = statevector;
state.sulfR_hf = statevector;
state.sulfR_df = statevector;

state.oxygenbentic_pf = statevector;
state.oxygenbentic_dif = statevector;
state.oxygenbentic_df = statevector;

state.sulfatebentic_pf = statevector;
state.sulfatebentic_dif = statevector;
state.sulfatebentic_df = statevector;

state.methanogenesis_pf = statevector;
state.methanogenesis_dif = statevector;
state.methanogenesis_df = statevector;

state.d13c_atmf = statevector;
%%%%%


pars.sFeIII_p = pars.sFeIII_p0*10^(ic-2);
pars.sFeIII_di = pars.sFeIII_di0*10^(ic-2);
pars.sFeIII_d =  pars.sFeIII_d0*10^(ic-2);
pars.k_FeIIIw = pars.k_FeIIIw0*10^(ic-2);
pars.FeIIIa_s = pars.FeIIIa_s0*10^(ic-2);

pars.k_phosw = pars.k_phosw0*10^(jc-1);

%%%%%%% run the system
[rawoutput.T,rawoutput.Y] = ode15s(@MBOX_equations,[pars.whenstart pars.whenend],pars.startstate,options);

%%%%% start time counter
tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Postprocessing   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% takes 'workingstate' from model and turns into 'state' %%%%%%%

%%%% size of output
pars.output_length = length(rawoutput.T) ;
%%%%%%%%%% model finished output to screen
fprintf('Integration finished \t') ; fprintf('Total steps: %d \t' , stepnumber ) ; fprintf('Output steps: %d \n' , pars.output_length )
toc

%%%%%%%%% print final model states using final state for each timepoint
%%%%%%%%% during integration
fprintf('assembling state vectors... \t')
tic

%%%% trecords is index of shared values between ode15s output T vector and
%%%% model recorded workingstate t vector
[sharedvals,trecords] = intersect(workingstate.time,rawoutput.T,'stable') ;

%%%%%% assemble output state vectors
field_names = fieldnames(workingstate) ;
for numfields = 1:length(field_names)
    eval([' state.' char( field_names(numfields) ) ' = workingstate.' char( field_names(numfields) ) '(trecords) ; '])
end

%%%%%% done message
fprintf('Done: ')
endtime = toc ;
fprintf('time (s): %d \n', endtime )



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Plotting script   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% load data

statefnumber = 41000
state.pO2_af = state.pO2_a(statefnumber)
state.Atmospheric_CO2_ppmf = state.Atmospheric_CO2_ppm(statefnumber)

state.DIC_conc_pf = state.DIC_conc_p(statefnumber)
state.DIC_conc_dif = state.DIC_conc_di(statefnumber)
state.DIC_conc_sf = state.DIC_conc_s(statefnumber)
state.DIC_conc_hf = state.DIC_conc_h(statefnumber)
state.DIC_conc_df = state.DIC_conc_d(statefnumber)

state.ALK_conc_pf = state.ALK_conc_p(statefnumber)
state.ALK_conc_dif = state.ALK_conc_di(statefnumber)
state.ALK_conc_sf = state.ALK_conc_s(statefnumber)
state.ALK_conc_hf = state.ALK_conc_h(statefnumber)
state.ALK_conc_df = state.ALK_conc_d(statefnumber)

state.pH_pf = state.pH_p(statefnumber)
state.pH_dif = state.pH_di(statefnumber)
state.pH_sf = state.pH_s(statefnumber)
state.pH_hf = state.pH_h(statefnumber)
state.pH_df= state.pH_d(statefnumber)

state.T_sf = state.T_s(statefnumber)
state.T_hf = state.T_h(statefnumber)
state.T_df = state.T_d(statefnumber)
state.T_contf = state.T_cont(statefnumber)
state.GASTf = state.GAST(statefnumber)

state.ccdegf = state.ccdeg(statefnumber)
state.baswf = state.basw(statefnumber)
state.granwf = state.granw(statefnumber)
state.silwf = state.silw(statefnumber)
state.carbwf = state.carbw(statefnumber)

state.mccb_pf = state.mccb_p(statefnumber)
state.mccb_dif = state.mccb_di(statefnumber)
state.mccb_df = state.mccb_d(statefnumber)

state.POC_pf = state.POC_p(statefnumber)
state.POC_dif = state.POC_di(statefnumber)
state.POC_sf = state.POC_s(statefnumber)
state.POC_hf = state.POC_h(statefnumber)
state.POC_df = state.POC_d(statefnumber)

state.DP_conc_pf = state.DP_conc_p(statefnumber)
state.DP_conc_dif = state.DP_conc_di(statefnumber)
state.DP_conc_sf = state.DP_conc_s(statefnumber)
state.DP_conc_hf = state.DP_conc_h(statefnumber)
state.DP_conc_df = state.DP_conc_d(statefnumber)

state.O2_conc_pf = state.O2_conc_p(statefnumber)
state.O2_conc_dif = state.O2_conc_di(statefnumber)
state.O2_conc_sf = state.O2_conc_s(statefnumber)
state.O2_conc_hf = state.O2_conc_h(statefnumber)
state.O2_conc_df = state.O2_conc_d(statefnumber)

state.FeIII_conc_pf = state.FeIII_conc_p(statefnumber)
state.FeIII_conc_dif = state.FeIII_conc_di(statefnumber)
state.FeIII_conc_sf = state.FeIII_conc_s(statefnumber)
state.FeIII_conc_hf = state.FeIII_conc_h(statefnumber)
state.FeIII_conc_df = state.FeIII_conc_d(statefnumber)

state.SO4_conc_pf = state.SO4_conc_p(statefnumber)
state.SO4_conc_dif = state.SO4_conc_di(statefnumber)
state.SO4_conc_sf = state.SO4_conc_s(statefnumber)
state.SO4_conc_hf = state.SO4_conc_h(statefnumber)
state.SO4_conc_df = state.SO4_conc_d(statefnumber)

state.FeII_conc_pf = state.FeII_conc_p(statefnumber)
state.FeII_conc_dif = state.FeII_conc_di(statefnumber)
state.FeII_conc_sf = state.FeII_conc_s(statefnumber)
state.FeII_conc_hf = state.FeII_conc_h(statefnumber)
state.FeII_conc_df = state.FeII_conc_d(statefnumber)

state.H2S_conc_pf = state.H2S_conc_p(statefnumber)
state.H2S_conc_dif = state.H2S_conc_di(statefnumber)
state.H2S_conc_sf = state.H2S_conc_s(statefnumber)
state.H2S_conc_hf = state.H2S_conc_h(statefnumber)
state.H2S_conc_df = state.H2S_conc_d(statefnumber)

state.FeIIIwf = state.FeIIIw(statefnumber)

state.FeIIIscavenging_pf = state.FeIIIscavenging_p(statefnumber)
state.FeIIIscavenging_dif = state.FeIIIscavenging_di(statefnumber)
state.FeIIIscavenging_sf = state.FeIIIscavenging_s(statefnumber)
state.FeIIIscavenging_hf = state.FeIIIscavenging_h(statefnumber)
state.FeIIIscavenging_df = state.FeIIIscavenging_d(statefnumber)

state.pyritewf = state.pyritew(statefnumber)
state.sulfatewf = state.sulfatew(statefnumber)
%state.pyriteb_pf = state.pyriteb_p(statefnumber)
%state.pyriteb_dif = state.pyriteb_di(statefnumber)
state.sulfatebf = state.sulfateb(statefnumber)

state.pyF_pf = state.pyF_p(statefnumber)
state.pyF_dif = state.pyF_di(statefnumber)
state.pyF_sf = state.pyF_s(statefnumber)
state.pyF_hf = state.pyF_h(statefnumber)
state.pyF_df = state.pyF_d(statefnumber)

state.ironO_pf = state.ironO_p(statefnumber)
state.ironO_dif = state.ironO_di(statefnumber)
state.ironO_sf = state.ironO_s(statefnumber)
state.ironO_hf = state.ironO_h(statefnumber)
state.ironO_df = state.ironO_d(statefnumber)

state.SironR_pf = state.SironR_p(statefnumber)
state.SironR_dif = state.SironR_di(statefnumber)
state.SironR_sf = state.SironR_s(statefnumber)
state.SironR_hf = state.SironR_h(statefnumber)
state.SironR_df = state.SironR_d(statefnumber)

state.SideP_pf = state.SideP_p(statefnumber)
state.SideP_dif = state.SideP_di(statefnumber)
state.SideP_sf = state.SideP_s(statefnumber)
state.SideP_hf = state.SideP_h(statefnumber)
state.SideP_df = state.SideP_d(statefnumber)

state.pripr_pf = state.pripr_p(statefnumber)
state.pripr_sf = state.pripr_s(statefnumber)
state.pripr_hf = state.pripr_h(statefnumber)

state.remin_pf = state.remin_p(statefnumber)
state.remin_dif = state.remin_di(statefnumber)
state.remin_sf = state.remin_s(statefnumber)
state.remin_hf = state.remin_h(statefnumber)
state.remin_df = state.remin_d(statefnumber)

state.mocb_pf = state.mocb_p(statefnumber)
state.mocb_dif = state.mocb_di(statefnumber)
state.mocb_df = state.mocb_d(statefnumber)
state.phoswf  = state.phosw(statefnumber)

state.sulfatebentic_pf = state.sulfatebentic_p(statefnumber)
state.sulfatebf = state.sulfateb(statefnumber)
state.sulfR_pf = state.sulfR_p(statefnumber)
state.sulfO_pf  = state.sulfO_p(statefnumber)
state.SironR_pf = state.SironR_p(statefnumber)
state.pyF_pf  = state.pyF_p(statefnumber)

state.pyritewf  = state.pyritew(statefnumber)
state.sulfatewf  = state.sulfatew(statefnumber)

state.SO4_lf  = state.SO4_l(statefnumber)

state.mocb_p_FeIIIf  = state.mocb_p_FeIII(statefnumber)
state.mocb_di_FeIIIf  = state.mocb_di_FeIII(statefnumber)
state.mocb_d_FeIIIf  = state.mocb_d_FeIII(statefnumber)

state.ocbother_pf  = state.ocbother_p(statefnumber)
state.ocbother_dif  = state.ocbother_di(statefnumber)
state.ocbother_df  = state.ocbother_d(statefnumber)

state.water_sediment_pf  = state.water_sediment_p(statefnumber)
state.water_sediment_dif  = state.water_sediment_di(statefnumber)
state.water_sediment_df  = state.water_sediment_d(statefnumber)

BEf_p = state.mocb_pf/state.water_sediment_pf;
BEf_di = state.mocb_dif/state.water_sediment_dif;
BEf_d = state.mocb_df/state.water_sediment_df;

state.AR_pf = state.AR_p(statefnumber);
state.AR_dif = state.AR_di(statefnumber);
state.AR_sf = state.AR_s(statefnumber);
state.AR_hf = state.AR_h(statefnumber);
state.AR_df = state.AR_d(statefnumber);

state.ironR_pf = state.ironR_p(statefnumber);
state.ironR_dif = state.ironR_di(statefnumber);
state.ironR_sf = state.ironR_s(statefnumber);
state.ironR_hf = state.ironR_h(statefnumber);
state.ironR_df = state.ironR_d(statefnumber);

state.sulfR_pf = state.sulfR_p(statefnumber);
state.sulfR_dif = state.sulfR_di(statefnumber);
state.sulfR_sf = state.sulfR_s(statefnumber);
state.sulfR_hf = state.sulfR_h(statefnumber);
state.sulfR_df = state.sulfR_d(statefnumber);

state.oxygenbentic_pf = state.oxygenbentic_p(statefnumber);
state.oxygenbentic_dif = state.oxygenbentic_di(statefnumber);
state.oxygenbentic_df = state.oxygenbentic_d(statefnumber);

state.sulfatebentic_pf = state.sulfatebentic_p(statefnumber);
state.sulfatebentic_dif = state.sulfatebentic_di(statefnumber);
state.sulfatebentic_df = state.sulfatebentic_d(statefnumber);

state.methanogenesis_pf = state.methanogenesis_p(statefnumber);
state.methanogenesis_dif = state.methanogenesis_di(statefnumber);
state.methanogenesis_df = state.methanogenesis_d(statefnumber);

state.d13c_atmf = state.d13c_atm(statefnumber);


%%% results processing
ARs_p = -state.oxygenbentic_pf*pars.Red_C_O;
ARs_di = -state.oxygenbentic_dif*pars.Red_C_O;
ARs_d  = -state.oxygenbentic_df*pars.Red_C_O;

sulfRs_p = -state.sulfatebentic_pf*2*pars.Red_C_O;
sulfRs_di = -state.sulfatebentic_dif*2*pars.Red_C_O;
sulfRs_d  = -state.sulfatebentic_df*2*pars.Red_C_O;

methans_p = state.methanogenesis_pf*2;
methans_di = state.methanogenesis_dif*2;
methans_d = state.methanogenesis_df*2;



%%%output
ColName = {'i';'j';'pO2_af';'O2_conc_df';'SO4_conc_df';'mocb_Fe';'ocbother';'pripr';'BEf_p';'BEf_di';'BEf_d';'AR_p';'AR_di';'AR_s';'AR_h';'AR_d';'ironR_p';'ironR_di';'ironR_s';'ironR_h';'ironR_d';'sulfR_p';'sulfR_di';'sulfR_s';'sulfR_h';'sulfR_d';'ARs_p';'ARs_di';'ARs_d';'sulfRs_p';'sulfRs_di';'sulfRs_d';'methans_p';'methans_di';'methans_d'};
mocb_Fe = state.mocb_p_FeIIIf+state.mocb_di_FeIIIf+state.mocb_d_FeIIIf;
ocbother = state.ocbother_pf+state.ocbother_dif+state.ocbother_df;
pripr = state.pripr_pf + state.pripr_sf + state.pripr_hf;

T = [i',j',state.pO2_af',state.O2_conc_df',state.SO4_conc_df',mocb_Fe',ocbother',pripr',BEf_p',BEf_di',BEf_d',state.AR_pf',state.AR_dif',state.AR_sf',state.AR_hf',state.AR_df',state.ironR_pf',state.ironR_dif',state.ironR_sf',state.ironR_hf',state.ironR_df',state.sulfR_pf',state.sulfR_dif',state.sulfR_sf',state.sulfR_hf',state.sulfR_df',ARs_p',ARs_di',ARs_d',sulfRs_p',sulfRs_di',sulfRs_d',methans_p',methans_di',methans_d'];


save('results.ascii','T','-ascii','-append','-tabs')




catch exception
    continue
end
end
end
