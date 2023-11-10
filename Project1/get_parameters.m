function  [aa ai ar I0 n0 T f alpha] = get_parameters(birthdate)

rng(birthdate);

g = 9.81;          
kB = 1.38e-23;     
q = 1.60e-19;     
epsilon0 = 8.85e-12;   
me = 9.1e-31;     
mp = 1.67e-27;
mO = 16*mp; 
p = 1e5;

aa = random('unif',8.0e-22,9.0e-22);   
ai = random('unif',7.0e-23,11.0e-23);
ar = random('unif',0.5e-16,2.0e-16);    
I0 = random('unif',1e14,2e14);
T = random('unif',1000,1500);
n0 = random('unif',2e18,3e18);
alpha = random('unif',25,65);
 
H = kB*T/(mO*g);    
zmax = H*log(aa*n0*H);
nemax = sqrt(ai/ar*n0*exp(-zmax/H)*I0*exp(-H*aa*n0*exp(-zmax/H)));
fpemax = 1/(2*pi)*sqrt(nemax*q^2/(epsilon0*me));
factor = random('unif',0.6,0.9);
f = factor*fpemax;


disp(' ');
disp('====================================');
disp('======  Modelling parameters  ======');
disp('====================================');
disp('');
str = ['aa = ' num2str(aa,'%.2e') '    (m^2)'];
disp(str);
str = ['ai = ' num2str(ai,'%.2e') '    (m^2)'];
disp(str);
str = ['ar = ' num2str(ar,'%.2e') '    (m^3 s^-1)'];
disp(str);
str = ['I0 = ' num2str(I0,'%.2e') '    (photons m^-2 s^-1)'];
disp(str);
str = ['n0 = ' num2str(n0,'%.2e') '    (m^-3)'];
disp(str);
str = ['T = ' num2str(T,'%.2e') '     (K)'];
disp(str);
str = ['f = ' num2str(f,'%.2e') '     (Hz)'];
disp(str);
str = ['alpha = ' num2str(alpha,'%2.0f') '       (degrees)'];
disp(str);
disp(' ');


return;