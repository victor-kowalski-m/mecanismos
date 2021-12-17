clc
clear
close all

%% Polinômio interpolador

T = 1;
Rpin1 = [330 -105 30]';
Rpfin2 = [240 180 30]';
[Rpfin1, Rpin2] = deal((Rpfin2 + Rpin1)/2 + [0 0 40]');

Rp = @(t, Rpin, Rpfin) ...
    Rpin + (6*(t/T)^5 - 15*(t/T)^4 + 10*(t/T)^3)*(Rpfin - Rpin);
Vp = @(t, Rpin, Rpfin) ...
    (30*(t/T)^4 - 60*(t/T)^3 + 30*(t/T)^2)*(Rpfin - Rpin)/T;
Ap = @(t, Rpin, Rpfin) ...
    (120*(t/T)^3 - 180*(t/T)^2 + 60*(t/T))*(Rpfin - Rpin)/T^2;

tempo  = linspace(0,1,21);
RP = zeros(3, 2*length(tempo) - 1);
VP = zeros(3, 2*length(tempo) - 1);
AP = zeros(3, 2*length(tempo) - 1);

for k = 1:length(tempo)-1
       RP(:, k) = Rp(tempo(k), Rpin1, Rpfin1);
       VP(:, k) = Vp(tempo(k), Rpin1, Rpfin1);
       AP(:, k) = Ap(tempo(k), Rpin1, Rpfin1);
end
for k = 1:length(tempo)
       RP(:, length(tempo)-1+k) = Rp(tempo(k), Rpin2, Rpfin2);
       VP(:, length(tempo)-1+k) = Vp(tempo(k), Rpin2, Rpfin2);
       AP(:, length(tempo)-1+k) = Ap(tempo(k), Rpin2, Rpfin2);
end

%% Cinemática inversa

syms RP0
q = [];
for i=1:3
   var_str = strcat('q',num2str(i),'(t)');
   syms(var_str); 
   q = [q str2sym(var_str)];
end

l1 = 250; 
l2 = 350; 
l3 = 175; 

T10 = @(q1)...
    [cos(q1) -sin(q1) 0 0  ;...
     sin(q1) cos(q1)  0 0  ;...
     0         0      1 0 ;...
     0         0      0 1  ];

T21 = @(q2)...
    [cos(q2) 0 sin(q2)  0 ;...
     0       1  0       0  ;...
    -sin(q2) 0  cos(q2) l1  ;...
     0       0  0       1  ];

T32 = @(q3)...
    [cos(q3) 0  sin(q3) l2 ;...
     0       1  0       0  ;...
    -sin(q3) 0  cos(q3) 0  ;...
     0       0  0       1  ];

F = @(q, RP0)...
    T10(q(1))*...
    T21(q(2))*...
    T32(q(3))*...
    [l3; 0; 0; 1] - RP0;

strFp = char(diff(F(q, RP0), t));
strFpp = char(diff(F(q, RP0), t, 2));
olds = {'q_(t)', 'diff(q(_), t)', 'diff(q(_), t, t)'};
news = {'q(_)', 'qp(_)', 'qpp(_)'};
for j=1:length(olds)
    for i=1:length(q)
        old = strrep(olds(j), '_', num2str(i));
        new = strrep(news(j), '_', num2str(i));
        strFp = strrep(strFp, old, new);
        strFpp = strrep(strFpp, old, new);
    end
end        
Fp = eval(['@(q, qp, VP0)' char(strFp) '- VP0']);
Fpp = eval(['@(q, qp, qpp, AP0)' char(strFpp) '- AP0']);

pontos = 2*length(tempo) - 1;
theta = zeros(3, pontos);
thetap = zeros(3, pontos);
thetapp = zeros(3, pontos);

for k = 1:pontos   
    theta(:, k) = fsolve(@(q) F(q, [RP(:, k); 1]), [0 0 pi/2]);
    thetap(:, k) = fsolve(@(qp) Fp(theta(:, k), qp, [VP(:, k); 0]), [0 0 0]);
    thetapp(:, k) = fsolve(@(qpp) Fpp(theta(:, k), thetap(:, k), qpp, [AP(:, k); 0]), [0 0 0]);
end

%% Plot 3D

figure("Name", "Trajetória da garra", 'NumberTitle','off');
[x, y, z] = meshgrid(-50:1:500, -200:-50, 0);
surf(x,y,z, 'LineStyle', 'none', 'FaceColor', [50 50 50]/155); hold on
[x, y, z] = meshgrid(-50:1:500, 50:1:200, 0);
surf(x,y,z, 'LineStyle', 'none', 'FaceColor', [50 50 50]/155);
grid on
xlabel("x (mm)")
ylabel("y (mm)")
zlabel("z (mm)")
title("Trajetória da garra")
plot_circulo([0 0 20], 20, 20, [150 150 150]);

i = 1;
while i < 5

    if i > 1, delete(h2), end
    
    for k=1:pontos
        
        ang = theta(:, k);
        R0 = RP(:, k);
        
        O = [0 0 0]';
        A = T10(ang(1))*[0 0 l1 1]';
        A = A(1:3);
        B = T10(ang(1))*T21(ang(2))*[l2 0 0 1]';
        B = B(1:3);
        C = T10(ang(1))*T21(ang(2))*T32(ang(3))*[l3 0 0 1]';
        C = C(1:3);

        if k > 1 || i > 1, delete(h1), end
        
        h1(1) = plot3([O(1), A(1)],...
             [O(2), A(2)],...
             [O(3), A(3)], 'k-', 'Color', [150 150 150]/255, 'LineWidth', 3);
        h1(2) = plot3([A(1), B(1)],...
             [A(2), B(2)],...
             [ A(3), B(3)], 'k-o', 'Color', [150 150 150]/255, 'LineWidth', 3);
        h1(3) = plot3([B(1), C(1)],...
             [B(2), C(2)],...
             [B(3), C(3)], 'k', 'Color', [150 150 150]/255, 'LineWidth', 3);
        h1(4) = plot3(R0(1), R0(2), R0(3)-15, 'o', 'LineWidth', 10, 'Color', [156, 119, 65]/255);
                
        h1(5) = plot3([R0(1), R0(1)],...
        [R0(2)-10, R0(2)+10],...
        [R0(3), R0(3)], 'k', 'Color', [150 150 150]/255, 'LineWidth', 3);
    
        h1(6) = plot3([R0(1), R0(1)],...
        [R0(2)-10, R0(2)-10],...
        [R0(3), R0(3)-20], 'k', 'Color', [150 150 150]/255, 'LineWidth', 3);
    
        h1(7) = plot3([R0(1), R0(1)],...
        [R0(2)+10, R0(2)+10],...
        [R0(3), R0(3)-20], 'k', 'Color', [150 150 150]/255, 'LineWidth', 3);
        
        axis equal
        axis([-50 500 -300 300 0 300])
        
        h2(k) = plot3(R0(1), R0(2), R0(3)-15, 'r.');
        pause(0.1)

    end
    
i = i + 1;

end

%% Função auxiliar para plot da base do mecanismo

function plot_circulo(c, R, altura, cor)

    for r = linspace(0, R, 100)
        x = linspace(c(1)-r, c(1)+r, 100); 
        y = sqrt(r^2-(x-c(1)).^2) + c(2);
        plot3(x, y, ones(length(x), 1)*c(3), "k-", 'Color', cor/255);
        hold on
        plot3(x, 2*c(2)-y, ones(length(x), 1)*c(3), "k-", 'Color', cor/255);
    end
    for i = linspace(0, altura, 100)
        plot3(x, y, ones(length(x), 1)*c(3)-i, "k-", 'Color', cor/255);
        plot3(x, 2*c(2)-y, ones(length(x), 1)*c(3)-i, "k-", 'Color', cor/255);
    end
end
