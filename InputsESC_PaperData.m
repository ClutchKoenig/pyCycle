SYS.SolverTolerance = 1e-3;
SYS.NumberOfStoredDataPoints = 20;

F = 96485; %As/mol
R_m = 8.3145; %J/(mol*K)
MR = 0;

M_h2 = 2e-3;
M_h2o = 18e-3;
M_o2 = 32e-3;
M_n2 = 28e-3;

dR_Hm_H2 = 2.479e5; %J/mol


C_cat_ESC = 62;
C_an_ESC = 62;

speedPEN_ESC = 1;
speedCAT_ESC = 1;
speedIC_ESC = 1;
speedAN_ESC = 1;

p0 = 1e5; %Pa
%% Cell
% delta_el_ESC  = 165e-6; %m
% delta_an_ESC  = 45e-6; %m
% delta_cat_ESC  = 45e-6; %m
delta_el_ESC  = 90e-6; %m
delta_an_ESC  = 35e-6; %m
delta_cat_ESC  = 35e-6; %m
delta_pen_ESC  = delta_an_ESC  + delta_el_ESC  + delta_cat_ESC ; %m

% delta_pen_ESC = 160e-6; %m
N_seg_ESC  = 10;
N_seg_ESC_passive  = 8;

%length_total_ESC  = 0.11; %m
%width_total_ESC  = 0.1; %m

length_total_ESC  = 0.162; %m
width_total_ESC  = 0.142; %m


% length_ESC = 0.1; %m
% width_ESC = 0.1; %m
length_ESC = 0.09; %m
width_ESC = 0.142; %m

length_passive_ESC = length_total_ESC - length_ESC;
width_passive_ESC = width_total_ESC - width_ESC;

thickness_IC_tot_ESC  = 850e-6; %m
A_cell_ESC = length_ESC *width_ESC ; %m2

Aic_loss_act_ESC = 2*(length_ESC/N_seg_ESC)*(thickness_IC_tot_ESC + width_total_ESC);
Aic_loss_pas_ESC = 2*(length_passive_ESC/N_seg_ESC_passive)*(thickness_IC_tot_ESC + width_total_ESC);
Aic_loss_face_ESC = 0;


No_cat_channels_ESC = 24;
width_cat_channel_ESC = 0.0042; %m
height_cat_channel_ESC =  1/1000; %m
d_h_cat_ESC = 2*width_cat_channel_ESC*height_cat_channel_ESC/(width_cat_channel_ESC + height_cat_channel_ESC); %m
A_cross_cat_ESC = width_cat_channel_ESC * height_cat_channel_ESC * No_cat_channels_ESC; %m2

No_an_channels_ESC = 24;
width_an_channel_ESC = 0.0042; %m
height_an_channel_ESC = 1/1000; %m
d_h_an_ESC =  2*width_an_channel_ESC*height_an_channel_ESC/(width_an_channel_ESC + height_an_channel_ESC); %m
A_cross_an_ESC = width_an_channel_ESC * height_an_channel_ESC * No_an_channels_ESC; %m2

%% Convective Cat Heat Transfer Active
Nu_cat_ESC = 4.9;

A_ht_cat_ic_ESC = No_cat_channels_ESC*(2*height_cat_channel_ESC*length_ESC + width_cat_channel_ESC*length_ESC); %m2
A_ht_cat_pen_ESC = No_cat_channels_ESC*width_cat_channel_ESC*length_ESC; %m2
%% Convective Cat Heat Transfer Passive
A_ht_cat_ic_ESC_pas = No_cat_channels_ESC*(2*height_cat_channel_ESC*length_passive_ESC + width_cat_channel_ESC*length_passive_ESC); %m2
A_ht_cat_pen_ESC_pas  = No_cat_channels_ESC*width_cat_channel_ESC*length_passive_ESC; %m2
%% Convective An Heat Transfer Active
Nu_an_ESC = 4.9;

A_ht_an_ic_ESC = No_an_channels_ESC*(2*height_an_channel_ESC*length_ESC + width_an_channel_ESC*length_ESC); %m2
A_ht_an_pen_ESC = No_an_channels_ESC*width_an_channel_ESC*length_ESC; %m2

%% Convective An Heat Transfer Passive
A_ht_an_ic_ESC_pas = No_an_channels_ESC*(2*height_an_channel_ESC*length_passive_ESC + width_an_channel_ESC*length_passive_ESC); %m2
A_ht_an_pen_ESC_pas = No_an_channels_ESC*width_an_channel_ESC*length_passive_ESC; %m2
%% Conduction Active Active
lambda_PEN_ESC = 2;

Ax_lambda_pen_ESC =  width_total_ESC *delta_pen_ESC ; %m2
Ax_lambda_ESC = width_total_ESC *thickness_IC_tot_ESC ; %m2

dx_ESC  = length_ESC/N_seg_ESC; 
%% Conduction Active Passive
Ax_lambda_pen_ESC_pas = Ax_lambda_pen_ESC; %m2

Ax_lambda_ESC_pas  = Ax_lambda_ESC; %m2


dx_ESC_pas  = length_passive_ESC/N_seg_ESC_passive; 
%% Radiation Active & Passive
A_rad_pen_cat_ESC = A_ht_cat_pen_ESC/N_seg_ESC;
A_rad_pen_an_ESC = A_ht_an_pen_ESC/N_seg_ESC;
A_rad_ic_cat_ESC = (A_ht_cat_pen_ESC + A_ht_cat_ic_ESC)/N_seg_ESC;
A_rad_ic_an_ESC = (A_ht_an_pen_ESC + A_ht_an_ic_ESC)/N_seg_ESC;

A_rad_pen_cat_ESC_pas = A_ht_cat_pen_ESC_pas/N_seg_ESC_passive;
A_rad_pen_an_ESC_pas = A_ht_an_pen_ESC_pas/N_seg_ESC_passive;
A_rad_ic_cat_ESC_pas = (A_ht_cat_pen_ESC_pas + A_ht_cat_ic_ESC_pas)/N_seg_ESC_passive;
A_rad_ic_an_ESC_pas = (A_ht_an_pen_ESC_pas + A_ht_an_ic_ESC_pas)/N_seg_ESC_passive;

epsilon_pen_ESC = 0.4;
epsilon_ic_ESC = 0.3;
% epsilon_pen_ESC = 0.0;
% epsilon_ic_ESC = 0.0;

sigma_rad = 0;%5.670374419e-8;
phi_pen_ic_ESC = 1;
phi_ic_pen_cat_ESC = A_rad_pen_cat_ESC/A_rad_ic_cat_ESC;
phi_ic_pen_an_ESC = A_rad_pen_an_ESC/A_rad_ic_an_ESC;

% sigma_rad = 0;
% phi_pen_ic_ESC = 0;
% phi_ic_pen_cat_ESC = A_rad_pen_cat_ESC/A_rad_ic_cat_ESC;
% phi_ic_pen_an_ESC = A_rad_pen_an_ESC/A_rad_ic_an_ESC;


epsilon_pen_ic_cat_ESC = 1/((1/phi_pen_ic_ESC) + ((1/epsilon_pen_ESC) - 1) + ((1/epsilon_ic_ESC) - 1)*(A_rad_pen_cat_ESC/A_rad_ic_cat_ESC));
epsilon_pen_ic_an_ESC = 1/((1/phi_pen_ic_ESC) + ((1/epsilon_pen_ESC) - 1) + ((1/epsilon_ic_ESC) - 1)*(A_rad_pen_an_ESC/A_rad_ic_an_ESC));
epsilon_ic_pen_cat_ESC = 1/((1/phi_ic_pen_cat_ESC) + ((1/epsilon_ic_ESC) - 1) + ((1/epsilon_pen_ESC) - 1)*(A_rad_ic_cat_ESC/A_rad_pen_cat_ESC));
epsilon_ic_pen_an_ESC = 1/((1/phi_ic_pen_an_ESC) + ((1/epsilon_ic_ESC) - 1) + ((1/epsilon_pen_ESC) - 1)*(A_rad_ic_an_ESC/A_rad_pen_an_ESC));

epsilon_pen_ic_cat_ESC_pas = 1/((1/phi_pen_ic_ESC) + ((1/epsilon_pen_ESC) - 1) + ((1/epsilon_ic_ESC) - 1)*(A_rad_pen_cat_ESC_pas/A_rad_ic_cat_ESC_pas));
epsilon_pen_ic_an_ESC_pas = 1/((1/phi_pen_ic_ESC) + ((1/epsilon_pen_ESC) - 1) + ((1/epsilon_ic_ESC) - 1)*(A_rad_pen_an_ESC_pas/A_rad_ic_an_ESC_pas));
epsilon_ic_pen_cat_ESC_pas = 1/((1/phi_ic_pen_cat_ESC) + ((1/epsilon_ic_ESC) - 1) + ((1/epsilon_pen_ESC) - 1)*(A_rad_ic_cat_ESC_pas/A_rad_pen_cat_ESC_pas));
epsilon_ic_pen_an_ESC_pas = 1/((1/phi_ic_pen_an_ESC) + ((1/epsilon_ic_ESC) - 1) + ((1/epsilon_pen_ESC) - 1)*(A_rad_ic_an_ESC_pas/A_rad_pen_an_ESC_pas));


%% Fluid property data
T_0 = 1000; %K
%Hydrogen 
%molar enthalpy %J/(mol*K)
cmp_H2_a = 0;
cmp_H2_b = 2.200E-03;
cmp_H2_c = 2.804E+01;

%viscosity %Pa*s

%thermal conductivity %W/(mK)
a_lambda_h2 = 3.799454545e-4;
b_lambda_h2 = 0.0804581818;

%Water
%molar enthalpy %J/(mol*K)
cmp_H2O_a = 0;
cmp_H2O_b = 1.228E-02;
cmp_H2O_c = 2.902E+01;

%thermal conductivity %W/(mK)
lambda_H2O_a = 0;
lambda_H2O_b = 1.042E-04;
lambda_H2O_c = -1.590E-02;
%a + b*T + c*T^2
%lambda_H2O_a = -0.0086593279;
%lambda_H2O_b = 0.0000747671;
%lambda_H2O_c = 0.0000000294;
%Oxygen
%molar enthalpy %J/(mol*K)
cmp_O2_a = 0;
cmp_O2_b = 7.058E-03;
cmp_O2_c = 2.778E+01;


%thermal conductivity %W/(mK)
lambda_O2_a = 0;
lambda_O2_b = 6.116E-05;
lambda_O2_c = 1.128E-02;

%Nitrogen
%molar enthalpy %J/(mol*K)
cmp_N2_a = 0;
cmp_N2_b = 6.006E-03;
cmp_N2_c = 2.661E+01;

%thermal conductivity %W/(mK)
lambda_N2_a = -0.0000000137;
lambda_N2_b = 0.0000740423;
lambda_N2_c = 0.0052068912;

%specific entropy %kJ/(kgK)
%hydrogen 250-1000K
%s = a*ln(T) + b
entr_H2_a = 14.5435245901;
entr_H2_b = - 29.4584462581;

%water 400 - 1200K
%s = a + b*T + c*T^2
entr_H2O_a = 5.7223384381;
entr_H2O_b = 0.0051575476;
entr_H2O_c = -0.0000014616;

%oxygen 250-1200K
%s = a*ln(T) + b
entr_O2_a = 1.0053005351;
entr_O2_b = 0.6621118402;

%nitrogen 250-1200K
%s = a*ln(T) + b
entr_N2_a = 1.0878937988;
entr_N2_b = 0.6264841382;

%heat capacity in kJ/(kgK)
%hydrogen 250 - 1000 K
%c_p = a + b*T + c*T^2 + d*T^3 + e*T^4 + f*T^5
c_p_h2_a = 8.75948223230951000000;
c_p_h2_b = 0.04510779844943950000;
c_p_h2_c = - 0.00013838909011465500;
c_p_h2_d = 0.00000020637895412073;
c_p_h2_e = - 0.00000000014909309395;
c_p_h2_f = 0.00000000000004223795;

%water 400 - 1200K
%c_p = a + b*T + c*T^2 + d*T^3
c_p_h2o_a = 2.2806788744;
c_p_h2o_b = -0.0016535600;
c_p_h2o_c = 0.0000026290;
c_p_h2o_d = -0.0000000010;

%oxygen 250 - 1200K
%c_p = a + b*T + c*T^2 + d*T^3 + e*T^4
c_p_o2_a = 0.962655389253072;
c_p_o2_b = -0.000595928076977;
c_p_o2_c = 0.000002028431244;
c_p_o2_d = -0.000000001896252;
c_p_o2_e = 0.000000000000591;

%oxygen 250 - 1200K
%c_p = a + b*T + c*T^2 + d*T^3
c_p_o2_a = 1.0913006006;
c_p_o2_b = -0.0003661460;
c_p_o2_c = 0.0000007501;
c_p_o2_d = -0.0000000003;

%Viscosity
A_H2O_visc = 4.076228959276E-08;
B_H2O_visc = -3.01414856711936E-06;

A_H2_visc = -3.77906E-12;
B_H2_visc = 2.1613E-08;
C_H2_visc = 2.85329E-06;

A_O2_visc = -8.75780474E-12;
B_O2_visc = 5.10759E-08;
C_O2_visc = 6.91769E-06;

A_N2_visc = -7.57643E-12;
B_N2_visc = 4.29767E-08;
C_N2_visc = 6.28023E-06;


%viscosity %Pa*s
A_eta_h2 = 0.18024e-5;
B_eta_h2 = 0.27174e-7;
C_eta_h2 = -0.13395e-10;
D_eta_h2 = 0.00585e-12;
E_eta_h2 = -0.00104e-15;
%viscosity %Pa*s
A_eta_o2 = -0.10257e-5;
B_eta_o2 = 0.92625e-7;
C_eta_o2 = -0.80657e-10;
D_eta_o2 = 0.05113e-12;
E_eta_o2 = -0.01295e-15;

%viscosity %Pa*s
A_eta_h2o = -3.01414856711936E-06;
B_eta_h2o = 4.076228959276E-08;


%viscosity %Pa*s
A_eta_n2 = -0.01020e-5;
B_eta_n2 = 0.74785e-7;
C_eta_n2 = -0.59037e-10;
D_eta_n2 = 0.03230e-12;
E_eta_n2 = -0.00673e-15;



%%init data
T_pen_init_ESC = [673.5544337	674.4399837	675.7570638	677.5881919	684.5881274	693.0196508	701.6518421	710.4865985	719.3129271	727.8716494	735.8360463	742.7762226	748.0355688	750.0059444	744.8966483	743.803425	743.2133053	742.9563506] +100 + 273.15;
T_ic_init_ESC = [676.9207426	677.3423327	678.0917783	679.1265057	682.5735055	689.9224196	697.9235636	706.1279341	714.2342288	721.9512508	728.948732	734.8427339	739.2059511	741.6351507	741.9157971	741.9737581	742.0050473	742.0188377] +100 + 273.15;
T_air_init_ESC = [653.6283496	656.8343228	659.7282302	662.4065972	672.9734048	682.182244	690.9518885	699.6074031	708.1746149	716.5216087	724.4110116	731.5145986	737.3988411	741.3893571	741.3893571	741.6462067	741.8045509	741.9067466] +100 + 273.15;
T_fuel_init_ESC = [674.8570903	676.1003037	677.0928832	678.4498942	683.3882676	691.1658998	699.4209717	707.879135	716.2766322	724.3358825	731.7291936	738.0559662	742.7988642	745.0681543	745.0681543	743.3177151	742.7848523	742.5316641] +100 + 273.15;

Ucell_init_ESC = 0.7;
Current_init_ESC =0.2*[9.728396147	10.21534569	10.83176665	11.5147179	12.20237754	12.82359239	13.28675987	13.47297832	13.22563192	12.29843357];
