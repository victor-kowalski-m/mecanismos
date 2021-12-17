clc
clear
close

% Elipse da trajetória
elipse = [30, 10];
y = @(x) sqrt(elipse(2)^2*(1-x^2/elipse(1)^2)); % Equação de elipse

% Pontos e direções pré-definidos
P1=[0; 0; 0]; % coordenadas do ponto P1
P2=[15; 5; 0]; % coordenadas do ponto P2
P3=[-7.5; y(-7.5); 0]; % coordenadas do ponto P3
P = [P1 P2 P3];
alpha21=20*pi/180; % valor predeterminado de alpha21
alpha31=20*pi/180; % valor predeterminado de alpha31

% Lado direito %%%%%%%%%%%%%%%%%%%

W=12; % valor atribuido a W
phi=175*pi/180; % valor atribuido a phi

% Cálculo do sistema de equações, sendo x = [Z, beta1, beta21, beta31]
E = @(x) [...
W*cos(x(2)+x(3)) + x(1)*cos(phi+alpha21) + P1(1)-P2(1) - x(1)*cos(phi) - W*cos(x(2));...
W*sin(x(2)+x(3)) + x(1)*sin(phi+alpha21) + P1(2)-P2(2) - x(1)*sin(phi) - W*sin(x(2));...
W*cos(x(2)+x(4)) + x(1)*cos(phi+alpha31) + P1(1)-P3(1) - x(1)*cos(phi) - W*cos(x(2));...
W*sin(x(2)+x(4)) + x(1)*sin(phi+alpha31) + P1(2)-P3(2) - x(1)*sin(phi) - W*sin(x(2))];
x0=[20; 240*pi/180; 100*pi/180; 250*pi/180]; %valores do chute inicial para uso do metodo de Newton-Raphson
YE=fsolve(E,x0); %resultado;
YEgraus = [YE(1); 180/pi*YE(2:4);];
Z = YE(1);
beta1 = YE(2);
beta21 = YE(3);
beta31 = YE(4);

% Lado esquerdo %%%%%%%%%%%%%%%%%%%
U=40; % valor atribuido a W
sigma=5*pi/180; % valor atribuido a phi

% Cálculo do sistema de equações, sendo q = [S, gamma1, gamma21, gamma31]
D = @(q) [...
U*cos(q(2)+q(3)) + q(1)*cos(sigma+alpha21) + P1(1)-P2(1) - q(1)*cos(sigma) - U*cos(q(2));...
U*sin(q(2)+q(3)) + q(1)*sin(sigma+alpha21) + P1(2)-P2(2) - q(1)*sin(sigma) - U*sin(q(2));...
U*cos(q(2)+q(4)) + q(1)*cos(sigma+alpha31) + P1(1)-P3(1) - q(1)*cos(sigma) - U*cos(q(2));...
U*sin(q(2)+q(4)) + q(1)*sin(sigma+alpha31) + P1(2)-P3(2) - q(1)*sin(sigma) - U*sin(q(2))];
angs=[20; 230*pi/180; 60*pi/180;20*pi/180]; %valores do chute inicial para uso do metodo de Newton-Raphson
YD=fsolve(D,angs); %resultado;
YDgraus = [YD(1); 180/pi*YD(2:4);];
S = YD(1);
gamma1 = YD(2);
gamma21 = YD(3);
gamma31 = YD(4);

% Cálculo de posições dos pontos A e B
A1=[P1(1)-Z*cos(phi); P1(2)-Z*sin(phi);0];
A0=[A1(1)-W*cos(beta1); A1(2)-W*sin(beta1);0];
A2=[A0(1)+W*cos(beta1+beta21); A0(2)+W*sin(beta1+beta21);0];
A3=[A0(1)+W*cos(beta1+beta31); A0(2)+W*sin(beta1+beta31);0];
A = [A1 A2 A3];
B1=[P1(1)-S*cos(sigma); P1(2)-S*sin(sigma);0];
B0=[B1(1)-U*cos(gamma1); B1(2)-U*sin(gamma1);0];
B2=[B0(1)+U*cos(gamma1+gamma21);B0(2)+U*sin(gamma1+gamma21);0];
B3=[B0(1)+U*cos(gamma1+gamma31);B0(2)+U*sin(gamma1+gamma31);0];
B = [B1 B2 B3];

A1B1 = norm(B1-A1); % comprimento da barra A1B1
A0B0 = norm(B0-A0); % comprimento da barra A0B0
A1A0 = norm(A1-A0); % comprimento da barra A1B0
B1B0 = norm(B1 - B0); % comprimento da barra B1B0

% Vetor que de A0 a B0
vecA0B0 = B0 - A0; % vetor
ang_hor_vecA0B0 = atan2(abs(vecA0B0(2)), abs(vecA0B0(1)));

% Posição do ponto P em relação à barra A1B1
anguloP_BA=subspace((P1-A1),(B1-A1));
dist_P_BA = Z*sin(anguloP_BA);
dist_projP_BA_A = Z*cos(anguloP_BA);

% Matrizes de transformação
T21 = @(theta)...
    [cos(theta) -sin(theta) 0 A1A0  ;...
     sin(theta) cos(theta)  0 0  ;...
     0          0           1 0 ;...
     0          0           0 1 ];
 
 T10 = @(theta)...
    [cos(theta) -sin(theta) 0 0  ;...
     sin(theta) cos(theta)  0 0  ;...
     0          0           1 0 ;...
     0          0           0 1 ];
 
 T40 = @(theta)...
     [cos(theta) -sin(theta) 0 vecA0B0(1)  ;...
      sin(theta) cos(theta)  0 vecA0B0(2)  ;...
      0          0           1 0 ;...
      0          0           0 1 ];

 % Vetor de angulos de entrada
num_angs = 100;
angs = linspace(0, 2*pi, num_angs);

% Insere betas dos três pontos P se não estiverem no vetor de angulos
betas = mod([beta1 beta1+beta21 beta1+beta31], 2*pi);
angs_com_betas = angs;
j = 0; % contador de betas inseridos
for i=1:num_angs-1
    for beta=betas
        if angs(i) < beta && angs(i+1) > beta
            angs_com_betas = [angs_com_betas(1:i+j) beta angs_com_betas(i+1+j:end)];
            j = j + 1;
        end
    end
end
angs = angs_com_betas;
num_angs = num_angs+j;

% Matrizes para guardar posição dos pontos a cada angulo de entrada
As = zeros(3, num_angs);
Bs = zeros(3, num_angs);
Cs = zeros(3, num_angs);
Ds = zeros(3, num_angs);

% Calcula a posição de pontos do mecanismo para cada angulo de entrada 
for k=1:num_angs  
    F = @(q)...
    [cos(angs(k)) -sin(angs(k)) 0 0  ;...
     sin(angs(k)) cos(angs(k))  0 0  ;...
     0           0            1 0 ;...
     0           0            0 1  ]*...
    [cos(q(1)) -sin(q(1)) 0 A1A0  ;...
     sin(q(1)) cos(q(1))  0 0  ;...
     0         0          1 0 ;...
     0         0          0 1  ]*...
    [A1B1; 0; 0; 1] - ...
    [cos(q(2)) -sin(q(2)) 0 vecA0B0(1)  ;...
     sin(q(2)) cos(q(2))  0 vecA0B0(2)  ;...
     0         0          1 0 ;...
     0         0          0 1  ]*...
    [B1B0; 0; 0; 1];
    if (mod(angs(k), 2*pi) > pi/2 && mod(angs(k), 2*pi) <= pi)
        theta0 = [pi/6 3*pi/2];
    elseif mod(angs(k), 2*pi) > pi && mod(angs(k), 2*pi) <= 7*pi/6
        theta0 = [2*pi 3*pi/2];
    elseif mod(angs(k), 2*pi) > 3*pi/2+pi/6
        theta0 = [pi 3*pi/2];
    else
        theta0 = [3*pi/2 3*pi/2];
    end
    theta0grad = theta0*180/pi;
    thetas = fsolve(F, theta0);
    thetasgrad = thetas*180/pi;
    
    Ak=[A0; 1]+T10(angs(k))*[A1A0; 0; 0; 1];
    Bk=[A0; 1]+T40(thetas(2))*[B1B0; 0; 0; 1];
    Ck=[A0; 1]+T40(thetas(2))*[-10; 0; 0; 1];
    Dk=[A0; 1]+T40(thetas(2))*[-20; 0; 0; 1];
    
    As(:, k) = Ak(1:3);
    Bs(:, k) = Bk(1:3);
    Cs(:, k) = Ck(1:3);
    Ds(:, k) = Dk(1:3);
      
end

% Função de um círculo
r = 5;
x_roda = linspace(A0(1)-r, A0(1)+r, 100); 
y_roda = sqrt(r^2-(x_roda-A0(1)).^2) + A0(2);

% Plots
count = 1;
max = 3;
marcacoes_P = 0;
while count<=max
    
    for k = 1:num_angs
        
        % Posição calculada dos pontos para o angulo
        Ak = As(:, k);
        Bk = Bs(:, k);
        Ck = Cs(:, k);
        Dk = Ds(:, k);
        
        % Versores paralelos e ortogonais à sapata
        versor_BA = (Bk-Ak)/norm(Bk-Ak);
        ortog_sapata = [-versor_BA(2); versor_BA(1); 0];        
                
        % Estrutura do mecanismo
        h(1) = plot([A0(1),Ak(1),Bk(1),B0(1), Ck(1)], ...
          [A0(2),Ak(2),Bk(2),B0(2), Ck(2)], "k-", 'LineWidth', 5);
        hold on
        h(end+1)= plot([A0(1),Ak(1),Bk(1),B0(1)], ...
          [A0(2),Ak(2),Bk(2),B0(2)], "ko", 'LineWidth', 3);
        h(end+1)=plot([Ck(1), Dk(1)], ...
          [Ck(2), Dk(2)], "k-", 'LineWidth', 10);
        h(end+1)=plot(A0(1), A0(2), "k-o", 'LineWidth', 5);
        h(end+1)=plot(x_roda, y_roda, "k-");
        h(end+1)=plot(x_roda, 2*A0(2)-y_roda, "k-");

        % Sapata
        Pk = Ak+dist_projP_BA_A*versor_BA-ortog_sapata*dist_P_BA;
        ini_sapata = Pk+7*versor_BA;
        fim_sapata = Pk-7*versor_BA;
        h(end+1)=plot([ini_sapata(1), fim_sapata(1)],[ini_sapata(2), fim_sapata(2)],"-k");
        h(end+1)=plot([ini_sapata(1), ini_sapata(1)-dist_P_BA*ortog_sapata(1)], ...
            [ini_sapata(2), ini_sapata(2)-dist_P_BA*ortog_sapata(2)],"-k");
        h(end+1)=plot([fim_sapata(1), fim_sapata(1)-dist_P_BA*ortog_sapata(1)], ...
            [fim_sapata(2), fim_sapata(2)-dist_P_BA*ortog_sapata(2)],"-k");
        h(end+1)=plot([ini_sapata(1), ini_sapata(1)+dist_P_BA*ortog_sapata(1)], ...
            [ini_sapata(2), ini_sapata(2)+dist_P_BA*ortog_sapata(2)],"-k");
        h(end+1)=plot([fim_sapata(1), fim_sapata(1)+dist_P_BA*ortog_sapata(1)], ...
            [fim_sapata(2), fim_sapata(2)+dist_P_BA*ortog_sapata(2)],"-k");
        
        % Trajetória
        h2(k) = plot(Pk(1), Pk(2), "r.");          
        if ((norm(Pk-P1) < 10e-8) || (norm(Pk-P2) < 10e-8) || (norm(Pk-P3) < 10e-8)) && marcacoes_P < 3
            marcacoes_P = marcacoes_P + 1;
            plot(Pk(1), Pk(2), "kx", 'MarkerSize', 10);
            txt = sprintf("(%.2f, %.2f)", Pk(1), Pk(2)) ;
            htext(marcacoes_P) = text(Pk(1)-5, Pk(2)-3, txt, "FontSize", 10);
            h2(k) = plot(Pk(1), Pk(2), "r.");    
        end
        
        % Configurações do plot
        grid on
        axis("equal")
        axis([-40,50,-10,60])
        xlabel("Coordenada X (m)")
        ylabel("Coordenada Y (m)")
        pause(0.001)
        
        % Se for ultima iteração, não deleta gráficos
        if count==max && k==num_angs
            break
        end
        delete(h);
        
    end
    
    % Se for ultima iteração, não deleta gráficos e traz textos para a frente
    if count==max
        delete(htext)
        for i=1:length(P)
            Pk = P(:, i);
            txt = sprintf("(%.2f, %.2f)", Pk(1), Pk(2)) ;
            text(Pk(1)-5, Pk(2)-3, txt, "FontSize", 10);
        end
        break
    end
    delete(h);
    delete(h2);
    count = count+1;    
    
end


