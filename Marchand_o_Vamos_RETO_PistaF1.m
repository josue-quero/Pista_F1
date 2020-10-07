
function Marchand_o_Vamos_RETO_Simulacion()
    f = figure;
    grid on;
    close all;
    syms x real;
    desicion = input("1) Colocacion segura gradas 2) Recorrido de la pista");
    if desicion == 1
        v = 69.44;
    elseif desicion == 2
        v = 26;
    end
    miu = 0.8;
    radio = 100;
    vMax = calcularV(radio, miu);
    xi = 100;
    xf = 2900;
    x = [100, 1100, 1800, 2800];
    y = [2900, 3500, 600, 1200];
    coeficientes = polyfit(x, y, 3);%Ajuste a un polinomio de grado n
    z = @(xx) spline(x, y, xx);
    hold on;
    fplot(z,[x(1),x(end)], '--b')

    
    coeficientesFinal = vpa(coeficientes);
    ecuacion = polinomial(coeficientesFinal);
    primeraDerivada = diff(ecuacion);
    segundaDerivada = diff(primeraDerivada);
    
    [maximo_x, maximo_y, minimo_x, minimo_y] = obtenerMaximoMinimo(primeraDerivada, segundaDerivada, ecuacion);
    fprintf("Maximo(x): %.6f Maximo(y): %.6f \n", maximo_x, maximo_y);
    fprintf("Minimo(x): %.6f Minimo(y): %.6f \n", minimo_x, minimo_y);
    
    paso1 = 1;
    [punto_derrape_maximo_x, punto_derrape_maximo_y] = calcular_punto_derrape(primeraDerivada, segundaDerivada, paso1, maximo_x, ecuacion);
    [punto_derrape_minimo_x, punto_derrape_minimo_y] = calcular_punto_derrape(primeraDerivada, segundaDerivada, paso1, minimo_x, ecuacion);
    fprintf("Punto derrape en zona maxima(x): %.6f Punto derrape en zona minima(x): %.6f \n", punto_derrape_maximo_x, punto_derrape_minimo_x);
    
    [recta_tangente_1, pendiente_1] = obtener_recta_tangente(primeraDerivada, punto_derrape_maximo_x, punto_derrape_maximo_y);
    [recta_tangente_2, pendiente_2] = obtener_recta_tangente(primeraDerivada, punto_derrape_minimo_x, punto_derrape_minimo_y);
    distancia = 20;
    paralela_derrape_maxima = obtener_paralela(recta_tangente_1, pendiente_1, distancia);
    paralela_derrape_minima = obtener_paralela(recta_tangente_2, pendiente_2, distancia);
    distancia = 10;
    paralela2_derrape_maxima = obtener_paralela(paralela_derrape_maxima, pendiente_1, distancia);
    paralela2_derrape_minima = obtener_paralela(paralela_derrape_minima, pendiente_2, distancia);
    distancia = 80;
    [perpendicular_derrape_maxima, pendiente_per_max] = obtener_perpendicular(paralela_derrape_maxima, maximo_y + 20, pendiente_1);
    perpendicular2_derrape_maxima = obtener_paralela(perpendicular_derrape_maxima, pendiente_per_max, distancia);
    [perpendicular_derrape_minima, pendiente_per_min] = obtener_perpendicular(paralela_derrape_minima, minimo_y - 20, pendiente_2);
    perpendicular2_derrape_minima = obtener_paralela(perpendicular_derrape_minima, pendiente_per_min, distancia);
   
    coor_max = obtener_coordenadas_gradas(paralela_derrape_maxima, paralela2_derrape_maxima, perpendicular_derrape_maxima, perpendicular2_derrape_maxima);
    coor_min = obtener_coordenadas_gradas(paralela_derrape_minima, paralela2_derrape_minima, perpendicular_derrape_minima, perpendicular2_derrape_minima);
    
    longitud = obtenerLongitud(primeraDerivada, xi, xf);
    disp("Longitud de la curva")
    disp(longitud)
    barrera = str2sym("-4.15*x+8055");
    graficar_gradas(paralela_derrape_maxima, paralela2_derrape_maxima, perpendicular_derrape_maxima, perpendicular2_derrape_maxima, coor_max);
    graficar_gradas(paralela_derrape_minima, paralela2_derrape_minima, perpendicular_derrape_minima, perpendicular2_derrape_minima, coor_min);
    graficar_rectas_tangentes(recta_tangente_1, recta_tangente_2, barrera);
    if v > vMax
        graficar_carro_derrape(ecuacion, recta_tangente_1, punto_derrape_maximo_x);
    else
        graficar_carro(ecuacion, maximo_y, minimo_y);
    end
    
    
end

function vMax = calcularV(radio, miu)
    vMax = sqrt(9.81*radio*miu);
end

function ecuacion = polinomial(coeficientes)
    syms x
    ecuacion = coeficientes(1)*x^3+coeficientes(2)*x^2+coeficientes(3)*x+coeficientes(4);
end

function evaluada = evaluacion(ecuacion, valor_x)
    syms x real;
    evaluada = subs(ecuacion, {x}, valor_x);
end

function [maximo_x, maximo_y, minimo_x, minimo_y] = obtenerMaximoMinimo(primeraDerivada, segundaDerivada, ecuacion)
    syms x
    ecuacion_derivada = primeraDerivada == 0;
    maximoMinimo = solve(ecuacion_derivada, x);
    [maximo_x, minimo_x] = evaluarMaximoMinimo(maximoMinimo, segundaDerivada);
    maximo_y = evaluacion(ecuacion, maximo_x);
    minimo_y = evaluacion(ecuacion, minimo_x);
end

function [maximo, minimo] = evaluarMaximoMinimo(maximoMinimo, segundaDerivada)
    syms x
    for i = 1:1:2
        x = maximoMinimo(i);
        evaluacion = subs(segundaDerivada);
        if evaluacion < 0
            maximo = maximoMinimo(i);
        elseif evaluacion > 0
            minimo = maximoMinimo(i);
        end
    end
end

function [radio, y] = obtenerRadio(primeraDerivada, segundaDerivada, punto, ecuacion)
    primeraEvaluada = evaluacion(primeraDerivada, punto);
    segundaEvaluada = evaluacion(segundaDerivada, punto);
    radio = ((1 + (primeraEvaluada)^2)^1.5) / (abs(segundaEvaluada));
    y = evaluacion(ecuacion, punto);
end

function longitud = obtenerLongitud(funcion, a, b)
    syms x real;
    ecuacion = sqrt(1 + funcion^2);
    n = 20;
    fx = evaluacion(ecuacion, a);
    h = (b - a) / n;
    
    for k = 0:1:((n / 2) - 1)
        fx = fx + 4 * evaluacion(ecuacion, a + ((2 * k + 1)*h));
    end
    
    for k = 1:1:((n / 2) - 1)
        fx = fx + 2 * evaluacion(ecuacion, a + (2 * k * h));
    end 
    fx = fx + evaluacion(ecuacion, b);
    longitud = fx * (h / 3);
end

function [punto_derrape_x, punto_derrape_y] = calcular_punto_derrape(primeraDerivada, segundaDerivada, paso1, maximo, ecuacion)
    for i = maximo - 60:paso1:maximo + 60
        [radio, y] = obtenerRadio(primeraDerivada, segundaDerivada, i, ecuacion);
        if radio < 100
            punto_derrape_x = i;
            punto_derrape_y = y;
            break
        end
    end
end

function graficar_carro(ecuacion, maximo, minimo)
    syms x real;
    longitudinstante=0;
    vi=83.33;
    y=4;
    peso=823;
    coefFricc=1.7;
    g=9.81;
    t=0;
    cont = 1;
    ax = 0;
    ay = 0;
    fplot(recta_tangente_1);
    for i = 100:10:2800
            y(i) = subs(ecuacion, {x}, i);
            h(cont) = plot(i, y(i), 'b*');
            pbaspect([1 1 1]);
            xlim([100 2800]);
            ylim([-600 5000]);
            pause(0.0000000001);
            cont = 1 + cont;
            delete(findall(gcf,'type','text'));
            if cont > 2
                delete(h(cont - 2));
            end
            vx(cont-1) = encontrarVx(ecuacion, vi,i,y);
            vy(cont-1) = encontrarVy(ecuacion, vi,i,y);
            t=t+calcularTiempo(vx(cont-1),vi);
            energia = vi^2*823;
        if i > 110
            [ax, ay] = calcularAceleracion(vx(cont-1),vx(cont-2),vy(cont-1),vy(cont-2),t);
        end
            text(2000,4000,['Energia cinetica = ' num2str(double(energia))]);
            text(2000,3800,['Tiempo transcurrido = ' num2str(double(t))]);
            text(2000,3600,['Velocidad en X = ' num2str(double(vx(cont-1)))]);
            text(2000,3400,['Velocidad en Y = ' num2str(double(vy(cont-1)))]);
            text(2000,3200,['Aceleracion en X = ' num2str(double(ax))]);
            text(2000,3000,['Aceleracion en Y = ' num2str(double(ay))]);
            text(2000,2800,['Posicion en x = ' num2str(double(i))]);
            text(2000,2600,['Posicion en Y = ' num2str(double(y(i)))]);
        end
end
    

function [recta_tangente, pendiente] = obtener_recta_tangente(ecuacion, punto_x, punto_y)
    syms x real;
    pendiente = subs(ecuacion, {x}, punto_x);
    b = punto_y - pendiente*punto_x;
    recta_tangente = pendiente*x + b;
end

function paralela = obtener_paralela(recta_tangente, pendiente, distancia)
    syms x real;
    b = subs(recta_tangente, {x}, 0);
    ecuacion = ((abs(-b + x))/(sqrt(pendiente^2 + 1)))==distancia;
    b_paralela = solve(ecuacion, x);
    for i = 1:1:2
        if (b_paralela(i) > b) && (pendiente > 0) 
            paralela = pendiente*x + b_paralela(i);
        
        elseif (b_paralela(i) < b) && (pendiente < 0)
            paralela = pendiente*x + b_paralela(i);
        end
    end
    
end

function [perpendicular, pendiente_perpendicular] = obtener_perpendicular(recta_paralela, recta_maxima, pendiente)
    syms x real;
    pendiente_perpendicular = -(1/pendiente);
    coordenada_x = solve(recta_paralela == recta_maxima);
    
    b_perpendicular = solve(pendiente_perpendicular*coordenada_x + x == recta_maxima);
    perpendicular = pendiente_perpendicular*x + b_perpendicular; 
end

function coordenadas = obtener_coordenadas_gradas(paralela1, paralela2, perpendicular1, perpendicular2)
       syms x real;
       coordenadas(1) = solve(paralela1 == perpendicular2);
       coordenadas(2) = solve(paralela1 == perpendicular1);
       coordenadas(3) = solve(paralela2 == perpendicular2);
       coordenadas(4) = solve(paralela2 == perpendicular1);
end

function graficar_gradas(funcion1, funcion2, funcion3, funcion4, coor)
    syms x real;
    fplot(funcion1, [double(coor(3)) double(coor(2))]);
    fplot(funcion2, [double(coor(1)) double(coor(4))]);
    fplot(funcion3, [double(coor(4)) double(coor(2))]);
    fplot(funcion4, [double(coor(3)) double(coor(1))]);
end

function graficar_carro_derrape(ecuacion, recta_tangente_1, punto_derrape_maximo_x)
    syms x real;
    
    longitudinstante=0;
    vi=69.44;
    y=4;
    peso=823;
    coefFricc=1.7;
    g=9.81;
    t=0;
    ax = 0;
    ay = 0;
    cont = 1;
    
    for i = 100:10:punto_derrape_maximo_x
        y(i) = subs(ecuacion, {x}, i);
        h(cont) = plot(i, y(i), 'b*');
        xlim([100 2800]);
        ylim([-600 5000]);
        pause(0.0000000001);
        cont = 1 + cont;
        delete(findall(gcf,'type','text'));
        if cont > 2
            delete(h(cont - 2));
        end
        vx(cont-1) = encontrarVx(ecuacion, vi,i,y);
        vy(cont-1) = encontrarVy(ecuacion, vi,i,y);
        energia = vi^2*823;
        t=t+calcularTiempo(vx(cont-1),vi);
        if i > 110
            [ax, ay] = calcularAceleracion(vx(cont-2),vx(cont-1),vy(cont-2),vy(cont-1),t);
        end
        text(2000,4000,['Energia cinetica = ' num2str(double(energia))]);
        text(2000,3800,['Tiempo transcurrido = ' num2str(double(t))]);
        text(2000,3600,['Velocidad en X = ' num2str(double(vx(cont-1)))]);
        text(2000,3400,['Velocidad en Y = ' num2str(double(vy(cont-1)))]);
        text(2000,3200,['Aceleracion en X = ' num2str(double(ax))]);
        text(2000,3000,['Aceleracion en Y = ' num2str(double(ay))]);
        text(2000,2800,['Posicion en x = ' num2str(double(i))]);
        text(2000,2600,['Posicion en Y = ' num2str(double(y(i)))]);
    end
    delete(h(cont - 1));
    cont1 = 1;
    fplot(recta_tangente_1);
    for u = punto_derrape_maximo_x:10:835.4
        y1(cont1) = subs(recta_tangente_1, {x}, u);
        h1(cont1) = plot(u, y1(cont1), 'b*');
        xlim([500 1088]);
        ylim([4450 4600]);
        pause(0.0000000001);
        cont1 = 1 + cont1;
        delete(findall(gcf,'type','text'));
        if cont1 > 2
            delete(h1(cont1 - 2));
        end
        distancia = distanciaPunto(587.9452, 4515.6861, u, y1(cont1-1));
        energia2(cont1) = energia_perdida(distancia);
        [v] = velocidad_energia(energia2(cont1));
        vx1(cont1-1) = encontrarVx(recta_tangente_1, v,u,y1);
        vy1(cont1-1) = encontrarVy(recta_tangente_1, v,u,y1);
        t=t+calcularTiempo(vx1(cont1-1), vy1);
        text(950,4560,['Energia perdida = ' num2str(double(energia2(cont1)))]);
        text(950,4550,['Tiempo transcurrido = ' num2str(double(t))]);
        text(950,4540,['Velocidad en X = ' num2str(double(vx1(cont1-1)))]);
        text(950,4530,['Velocidad en Y = ' num2str(double(vy1(cont1-1)))]);
        text(950,4500,['Posicion en x = ' num2str(double(u))]);
        text(950,4490,['Posicion en Y = ' num2str(double(y1(cont1-1)))]);
        text(700, 4580, ('SE HA SALIDO DE LA PISTA'));
    end
    delete(h1(cont1 - 1));
    cont2 = 1;
    distancia = 720;
    for l = 835.4:-10:distancia
        y2(cont2) = subs(recta_tangente_1, {x}, l);
        h2(cont2) = plot(l, y2(cont2), 'b*');
        xlim([500 1088]);
        ylim([4450 4600]);
        pause(0.0000000001);
        cont2 = 1 + cont2;
        delete(findall(gcf,'type','text'));
        if cont2 > 2
            delete(h2(cont2 - 2));
        end
        t=t+calcularTiempo(vx(cont2-1), vy);
        distancia = distanciaPunto(835.4, 4588, l, y2(cont2-1));
        energia3(cont2) = energia_perdida(distancia);
        [v1] = velocidad_energia(energia3(cont2));
        vx2(cont2-1) = encontrarVx(recta_tangente_1, v1,l,y2);
        vy2(cont2-1) = encontrarVy(recta_tangente_1, v1,l,y2);
        text(700, 4580, ('HA COLISIONADO'));
        text(950,4560,['Energia perdida = ' num2str(double(energia3(cont2)))]);
        text(950,4550,['Tiempo transcurrido = ' num2str(double(t))]);
        text(950,4540,['Velocidad en X = ' num2str(double(vx2(cont2-1)))]);
        text(950,4530,['Velocidad en Y = ' num2str(double(vy2(cont2-1)))]);
        text(950,4500,['Posicion en x = ' num2str(double(l))]);
        text(950,4490,['Posicion en Y = ' num2str(double(y2(cont2-1)))]);
    end
end

function [vx] = encontrarVx(ecuacion, vi, xx, y)
    syms x;
    nueva_ecuacion=diff(ecuacion);
    m=subs(nueva_ecuacion,xx);
    b=y-m*xx;
    ecua=m*x+b;
    xxx=[0,10];
    yy=subs(ecua, xxx);
    toa=(yy(2)-yy(1))/(xxx(2)-xxx(1));
    angulo=vpa(atan(toa));
    vx=vpa(vi*cos(angulo));
end

function [vy] = encontrarVy(ecuacion, vi, xx, y)
    syms x;
    nueva_ecuacion=diff(ecuacion);
    m=subs(nueva_ecuacion,xx);
    b=y-m*xx;
    ecua=m*x+b;
    xxx=[0,10];
    yy=subs(ecua, xxx);
    toa=(yy(2)-yy(1))/(xxx(2)-xxx(1));
    angulo=vpa(atan(toa));
    vy=vpa(vi*sin(angulo));
end

function [t] = calcularTiempo(vx, vy)
    t=10/vx;
end

function [energia]=calcularEnergia(ecuacion, x1, x2, peso, g, coefFricc,longitudinstante)
    syms x;
    primeraDerivada=diff(ecuacion);
    expresion = sqrt(1+((primeraDerivada)^2));
    longitud =vpa(int(expresion,x1,x2));
    longitudinstante = longitudinstante + longitud;
    fuerza=peso*g*coefFricc;
    energia=fuerza*longitudinstante;
end

function [ax, ay] = calcularAceleracion(v1x,v2x,v1y,v2y,t)
    ax = (v2x-v1x)/t;
    ay = (v2y-v1y)/t;
end

function [distancia] = distancia_despues()
    v = velocidad_antes();
    coefFricc=0.9;
    m = 823;
    g=9.81 ;
    e = 0.15;
    v2 = v * e;
    ek = 1/2 * m * v2^2;
    deltat = 1;
    distancia=0;
    while ek>0
        deltad = v2 * deltat;
        ek = ek-m * coefFricc * g * deltad;
        distancia = distancia + v2 * deltat;
        v2 = vpa(sqrt( 2 * ek/m));
    end
end

function [distancia] = distanciaPunto(x1,y1,x2,y2)
    distancia = vpa(sqrt((x2-x1)^2+(y2-y1)^2));
end

function [v2] = velocidad_antes()
    vi=69.44;
    v2=vi;
    coefFricc=0.9 ;
    m=823;
    g=9.81;
    ek=1/2*m*v2^2;
    deltat=.1;
    distancia=distanciaPunto(587.9452,4515.6861,835.4,4588);
    while distancia>0
        deltad=v2*deltat;
        ek=ek-m*coefFricc*g*deltad;
        distancia=distancia-deltad;
        v2=vpa(sqrt(2*ek/m));
    end
end

function [energia]=energia_perdida(distancia)
    coefFricc=0.9;
    m=823;
    g=9.81;
    energia=coefFricc*m*g*distancia;
end

function [v] = velocidad_energia(energia)
    m=823;
    v=vpa(sqrt(2*energia/m));
    fprintf("La velocidad es %.6f m/s",v)
end

function graficar_rectas_tangentes(recta_tangente_1, recta_tangente_2, barrera)
    fplot(recta_tangente_1);
    fplot(recta_tangente_2);
    fplot(barrera);
end