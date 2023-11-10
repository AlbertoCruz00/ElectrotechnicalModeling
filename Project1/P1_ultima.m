close all; clear; clc;

k      = 1.38*10^-23;  % Boltzmann
m_e    = 9.1e-31;      % Electron mass
m_p    = 1.67*10^-27;  % Proton mass
m_02   = 16*m_p;       % O2 mass
q      = 1.6e-19;      % Unit charge
ep_0   = 8.85e-12;     % Vacuum permitivity
g      = 9.8;          % Gravity

% Step 1

[aa, ai, ar, I0, n0, T, f, alpha] = get_parameters(001126);

alpha = deg2rad(alpha);

% Step 2

T_1 = T+100;

z = 1e3*(0:1:1000); % Altitude vector (m)

H = k*T/(g*m_02);  % Atm scale height (m)
H_1 = k*T_1/(g*m_02);  % Atm scale height (m)

ne = sqrt( ai/ar .* n0 .* exp(- z/H) .* I0 .* exp(-H .* aa .* n0 .* exp(-z/H) ) );
ne_1 = sqrt( ai/ar .* n0 .* exp(- z/H_1) .* I0 .* exp(-H_1 .* aa .* n0 .* exp(-z/H_1) ) );


hfig = figure;
fname = 'Ex1_1';
semilogx(ne,z*1e-3,'LineWidth',1.5)
xlabel('Electron number density ($ number \; of \; electrons/m^3$)'); xlim([1e-8 1e13]);
ylabel('Height (Km)'); ylim([0 1e3]);
grid on
title('Height vs Electron density')
Figures

% hfig = figure;
% fname = 'Ex1_2';
% semilogx(ne,z*1e-3,ne_1,z*1e-3,'LineWidth',1.5)
% xlabel('Electron number density (number of electrons/m^3)'); xlim([1e-8 1e13]);
% ylabel('Height (km)'); ylim([0 1e3]);
% legend(['T = ',num2str(T),' K'],['T = ',num2str(T_1),' K'])
% grid on
% title('Height vs Electron density')
% Figures

%% PART 2

f_pe = 1/(2*pi) * sqrt( (ne * q^2) /( ep_0 * m_e)  );
f_pe_1 = 1/(2*pi) * sqrt( (ne_1 * q^2) /( ep_0 * m_e)  );


hfig = figure;
fname = 'Ex2_1';
plot(f_pe*1e-6,z*1e-3,'LineWidth',1.5)
xlabel('Plasma frequency (MHz)')
ylabel('Height (Km)')
grid on
title('Height vs Plasma frequency')
Figures

% hfig = figure;
% fname = 'Ex2_2';
% plot(f_pe*1e-6,z*1e-3,f_pe_1*1e-6,z*1e-3,'LineWidth',1.5)
% legend(['T = ',num2str(T),' K'],['T = ',num2str(T_1),' K'])
% xlabel('Plasma frequency (MHz)')
% ylabel('Height (km)')
% grid on
% title('Height vs Plasma frequency')
% Figures

%% Step 3 again

% Compute z reflection

p=1;
ind = 0;
while p<length(f_pe)
    if f <= f_pe(p)/cos(alpha)
        ind = p;
        break
    end
    p = p+1;
end

f_pe(ind);
z_ref = z(ind);

x_1 = z(1:ind)*sin(alpha);
z_1 = z(1:ind);

x_n = [x_1 x_1+x_1(end)];
z_n = [z_1,flip(z_1)];

hfig = figure;
fname = 'Ex3_1';
hold on
grid on
title('Simplified model for the radio wave propagation')
plot(x_n*1e-3,z_n*1e-3,'-','LineWidth',1.5)
yline(z_ref*1e-3,'-','Zref')
ylabel('Height (km)')
xlabel('Horizontal distance (km)')
Figures



%% Step 4

dz = 0.1;

zi = (0:dz:1000)*1e3; % Altitude vector (m)

H = k*T/(g*m_02);  % Atm scale height (m)
ne = sqrt( ai/ar .* n0 .* exp(- zi/H) .* I0 .* exp(-H .* aa .* n0 .* exp(-zi/H) ) );
f_pe = 1/(2*pi) * sqrt( (ne * q^2) /( ep_0 * m_e)  );

xi = zeros(length(zi),1);

ni = sqrt(1 - f_pe.^2./f^2 );

alpha_i = zeros(length(zi),1);

flag_1 = 0;

p = 1;

while flag_1 == 0
    if p == 1
        alpha_i(p) = alpha;
        xi(p) = 0;
    else
        alpha_i(p) = asin(sin(alpha_i(p-1))*ni(p-1)/ni(p));  
        if alpha_i(p) >= pi/2
            ind = p;
            flag_1 = 1;
            z_ref = zi(ind);
        end
        xi(p) = xi(p-1)+(zi(p)-zi(p-1))*tan(alpha_i(p));
    end
    % xi(p) = zi(p)*tan(alpha_i(p));
    p = p+1;

end

% xi = xi';
% 
% z_f = [zi(1:ind) flip(zi(1:ind))];
% x_f = [xi(1:ind) flip(xi(1:ind))+abs(xi(ind))];
% 
% figure
% plot(x_f,z_f)

z_final = [zi(1:ind), zi(ind:-1:1)];

xi_rev = zeros(ind, 1);
xi_end = xi(ind);

for i=1:1:ind
    if i == 1
        xi_rev(i) = xi_end;
    else
        xi_rev(i) = xi_rev(i-1)+(xi(ind-i+2)-xi(ind-i+1));
    end
end
% xi(ind)+ xi(ind:-1:1
x_final = [xi(1:ind)',xi_rev'];


hfig = figure;
fname = 'Ex4_1';
hold on, grid on
plot(x_final*1e-3,z_final*1e-3,'-','LineWidth',1.5), 
yline(zi(ind)*1e-3,'-','Zref')
xlabel('Horizontal distance (Km)')
ylabel('Height (km)')
title('More realistic model of the ray-tracing path')
Figures

hfig = figure;
fname = 'Ex4_2';
hold on, grid on
plot(x_n*1e-3,z_n*1e-3,'-','LineWidth',1.5,'DisplayName','Simple model')
plot(x_final*1e-3,z_final*1e-3,'-','LineWidth',1.5,'DisplayName','Realistic model'), 
xlabel('Horizontal distance (Km)')
ylabel('Height (km)')
legend;legend('boxoff')
title('Comparison of models')
Figures



%% Step 5 - evaluating the effects of changing the atmospheric temperature

dz = 0.001;

zi = (0:dz:1000)*1e3; % Altitude vector (m)

T_1 = T + 100; % New temperature (K)
H = k*T_1/(g*m_02);  % Atm scale height (m)
ne = sqrt( ai/ar .* n0 .* exp(- zi/H) .* I0 .* exp(-H .* aa .* n0 .* exp(-zi/H) ) );
f_pe = 1/(2*pi) * sqrt( (ne * q^2) /( ep_0 * m_e)  );

xi = zeros(length(zi),1);

ni = sqrt(1 - f_pe.^2./f^2 );

alpha_i = zeros(length(zi),1);

flag_1 = 0;

p = 1;

while flag_1 == 0
    if p == 1
        alpha_i(p) = alpha;
        xi(p) = 0;
    else
        alpha_i(p) = asin(sin(alpha_i(p-1))*ni(p-1)/ni(p));  
        if alpha_i(p) >= pi/2
            ind = p;
            flag_1 = 1;
            z_ref = zi(ind);
        end
        xi(p) = xi(p-1)+(zi(p)-zi(p-1))*tan(alpha_i(p));
    end
    p = p+1;
end

z_5 = [zi(1:ind), zi(ind:-1:1)];

xi_rev = zeros(ind, 1);
xi_end = xi(ind);

for i=1:1:ind
    if i == 1
        xi_rev(i) = xi_end;
    else
        xi_rev(i) = xi_rev(i-1)+(xi(ind-i+2)-xi(ind-i+1));
    end
end

x_5 = [xi(1:ind)',xi_rev'];


hfig = figure;
fname = 'Ex5_1';
hold on, grid on
%plot(x_n*1e-3,z_n*1e-3,'-','LineWidth',1.5)
plot(x_final*1e-3,z_final*1e-3,'-','LineWidth',1.5,'DisplayName',['T = ',num2str(T),' K'])
plot(x_5*1e-3,z_5*1e-3,'-','LineWidth',1.5,'DisplayName',['T = ',num2str(T_1),' K']);
legend;legend('boxoff')
xlabel('Horizontal distance (Km)')
ylabel('Height (km)')
title('Comparison of models')
Figures


%% Step 6

I_6=I0*10;

dz = 0.001;

zi = (0:dz:1000)*1e3; % Altitude vector (m)

H = k*T/(g*m_02);  % Atm scale height (m)
ne = sqrt( ai/ar .* n0 .* exp(- zi/H) .* I_6 .* exp(-H .* aa .* n0 .* exp(-zi/H) ) );
f_pe = 1/(2*pi) * sqrt( (ne * q^2) /( ep_0 * m_e)  );

xi = zeros(length(zi),1);

ni = sqrt(1 - f_pe.^2./f^2 );

alpha_i = zeros(length(zi),1);

flag_1 = 0;

p = 1;

while flag_1 == 0
    if p == 1
        alpha_i(p) = alpha;
        xi(p) = 0;
    else
        alpha_i(p) = asin(sin(alpha_i(p-1))*ni(p-1)/ni(p));  
        if alpha_i(p) >= pi/2
            ind = p;
            flag_1 = 1;
            z_ref = zi(ind);
        end
        xi(p) = xi(p-1)+(zi(p)-zi(p-1))*tan(alpha_i(p));
    end
    p = p+1;
end

z_6 = [zi(1:ind), zi(ind:-1:1)];

xi_rev = zeros(ind, 1);
xi_end = xi(ind);

for i=1:1:ind
    if i == 1
        xi_rev(i) = xi_end;
    else
        xi_rev(i) = xi_rev(i-1)+(xi(ind-i+2)-xi(ind-i+1));
    end
end

x_6 = [xi(1:ind)',xi_rev'];


hfig = figure;
fname = 'Ex6_1';
hold on, grid on
%plot(x_n*1e-3,z_n*1e-3,'-','LineWidth',1.5)
plot(x_final*1e-3,z_final*1e-3,'-','LineWidth',1.5,'DisplayName','Original case')
%plot(x_5*1e-3,z_5*1e-3,'-','LineWidth',1.5,'DisplayName','Larger temperature');
plot(x_6*1e-3,z_6*1e-3,'-','LineWidth',1.5,'DisplayName','Larger intensity');
legend;legend('boxoff','')
xlabel('Horizontal distance (Km)')
ylabel('Height (km)')
title('Comparison of models')
Figures


