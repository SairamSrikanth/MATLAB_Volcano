global rho A P0 h0 L0 Lp0 Patm mu rc rs Us dt G hslim;

Patm=100000;
rho=2600; 
mu=1000;
G=10^10;
rc=3;
Nf=Nf_fun(rho, mu, rc);
Fr=Fr_fun(Nf);
lbd=rc*lbd_p(Fr, Nf);
rs=rc-lbd;
A=1-(rs/rc)^2;
Us=Fr*sqrt(2*9.81*rc);
h0=200-5; %1458*461*M/3.14/rs^3/2/2400/9.81; %4475;
P0=rho*9.81*h0+Patm;
L0=5;%2*rs;
gamma=1.0;
Lp0=rho*9.81*Us*L0/gamma/P0;
dt=0.1;
hslim=1; %sqrt(A*P0*L0/rho/9.81)-Patm/rho/9.81;

sythetic_data = volcano_data_generator()

function sythetic_data = volcano_data_generator()
    global h0 L0 Lp0 rc rs dt Us hslim;

    k = zeros([4, 1]);
    j = zeros([4, 1]);
       
    t  = [0];
    L  = [L0];
    Lp = [Lp0];
    z1 = [50];
    z2 = [50+h0];
    z3 = [50+h0+2*rs];
    z4 = [50+h0+2*rs];
    signal = [0];

    for i = 1:6000
        k(1) = f2(t(i), L(i), Lp(i));
        j(1) = f1(Lp(i));

        k(2) = f2(t(i)+dt/2, L(i)+dt/2*j(1), Lp(i)+dt/2*k(1));
        j(2) = f1(Lp(i)+dt/2*k(1));

        k(3) = f2(t(i)+dt/2, L(i)+dt/2*j(2), Lp(i)+dt/2*k(2));
        j(3) = f1(Lp(i)+dt/2*k(2));

        k(4) = f2(t(i)+dt, L(i)+dt*j(3), Lp(i)+dt*k(3));
        j(4) = f1(Lp(i)+dt*k(3));

        t(i+1) = t(i) + dt;
        L(i+1) = L(i) + 1/6*dt*(j(1)+2*j(2)+2*j(3)+j(4));
        Lp(i+1) = Lp(i) + 1/6*dt*(k(1)+2*k(2)+2*k(3)+k(4));
        z1(i+1) = z1(1)-(L(i+1)-L(1))/rc^2*rs^2;
        z3(i+1) = z3(i)-Us*dt;
        z2(i+1) = z3(i+1)-L(i+1); 
        z4(i+1) = z4(1);        
        
        signal(i+1) = -tilt(z1(1), z3(1), z1(i+1), z2(i+1), z3(i+1), z4(i+1), Lp(i+1));
        
        if z2(i+1)-z1(i+1) < hslim
            break;
        end
        
%         plot(t, -signal);
%         pause(0.05);
        hold on; 
        plot(t, -z1, "g");
        plot(t, -z2, "r");
        plot(t, -z3, "b");
        plot(t, -z1-hslim, "c")
        hold off;
        pause(0.05);
        
    end
    
    sythetic_data = [t; L; Lp; z1; z2; z3; z4; signal];
end

function rhs = f1(Lp)
    rhs = Lp;
end

function rhs = f2(t, L, Lp)
    global rho A P0 L0 Patm mu rc;
    rhs = 2/(rho*(2-A))*(P0*L0/L/head(t, L) - rho*9.81 - Patm/head(t, L) - 8*mu/rc^2*Lp);
end

function output = Nf_fun(rho, mu, rc)
    output=rho/mu*sqrt(9.81*(2*rc)^3);
end

function output = Fr_fun(Nf)
    if Nf < 2
        output = 0.01*Nf;
    elseif Nf > 300
        output = 0.345;
    else 
        output = 0.345*(1-exp(-Nf/34.5));
    end
end

function output = lbd_p(Fr, Nf)
    output = (6*Fr/Nf)^(1/3);
end

function output = head(t, L)
    global h0 Us L0 rs rc;
    output = h0-Us*t-(L-L0)*(1-(rs/rc)^2);
end

function output = tilt(z1p, z3p, z1, z2, z3, z4, Lp)
    output = tilt_sigma(z1p, z3p, z1, z2, z3, z4)+tilt_tau(z1, z2, z3, Lp);
end

function output = tilt_sigma(z1p, z3p, z1, z2, z3, z4)
    output = tilt_sigma_s1(z1, z1p);
    output = output + tilt_sigma_s2(z1p, z2, z1);
    output = output + tilt_sigma_s3(z2, z3, z1p, z1);
    output = output + tilt_sigma_s4(z3, z3p, z1p, z1, z2);
    output = output + tilt_sigma_s5(z3p, z4);
end

function output = tilt_sigma_s1(z1, z1p)
    global rc rho G;
    
    output = rc^2*rho*9.81/G*3*500*int_exp2(z1, z1p);
    output = output - rc^2*rho*9.81/G*15*500/2*int_exp4(z1, z1p);
    output = output + rc^2*rho*9.81/G*3*0.25*int_exp2(z1, z1p);

    output = output - rc^2*rho*9.81*z1/G*3*500*int_exp1(z1, z1p);
    output = output + rc^2*rho*9.81*z1/G*15*500/2*int_exp3(z1, z1p);
    output = output - rc^2*rho*9.81*z1/G*3*0.25*int_exp1(z1, z1p);
end 

function output = tilt_sigma_s2(z1p, z2, z1)
    global rc rho G;
    
    output = rc^2*rho*9.81*(z1p-z1)/G*3*500*int_exp1(z1p, z2);
    output = output - rc^2*rho*9.81*(z1p-z1)/G*15*500/2*int_exp3(z1p, z2);
    output = output + rc^2*rho*9.81*(z1p-z1)/G*3*0.25*int_exp1(z1p, z2);
end 

function output = tilt_sigma_s3(z2, z3, z1p, z1)
    global rc rho A G;

    output = -rc^2*rho*9.81*(1-A)/G*3*500*int_exp2(z2,z3);
    output = output + rc^2*rho*9.81*(1-A)/G*15*500/2*int_exp4(z2,z3);
    output = output - rc^2*rho*9.81*(1-A)/G*3*0.25*int_exp2(z2,z3);

    output = output + rc^2*rho*9.81*(z1p-z1+(1-A)*z2)/G*3*500*int_exp1(z2,z3);
    output = output - rc^2*rho*9.81*(z1p-z1+(1-A)*z2)/G*15*500/2*int_exp3(z2,z3);
    output = output + rc^2*rho*9.81*(z1p-z1+(1-A)*z2)/G*3*0.25*int_exp1(z2,z3);
end 

function output = tilt_sigma_s4(z3, z3p, z1p, z1, z2)
    global rc rho A G;
    output = rc^2*rho*9.81*(z1p-z1-(1-A)*(z3-z2))/G*3*500*int_exp1(z3,z3p);
    output = output - rc^2*rho*9.81*(z1p-z1-(1-A)*(z3-z2))/G*15*500/2*int_exp3(z3,z3p);
    output = output + rc^2*rho*9.81*(z1p-z1-(1-A)*(z3-z2))/G*3*0.25*int_exp1(z3,z3p);
end

function output = tilt_sigma_s5(z3p, z4)
    global rc rho A G;

    output = rc^2*rho*9.81*(1-A)/G*3*500*int_exp2(z3p, z4);
    output = output - rc^2*rho*9.81*(1-A)/G*15*500/2*int_exp4(z3p, z4);
    output = output + rc^2*rho*9.81*(1-A)/G*3*0.25*int_exp2(z3p, z4);

    output = output - rc^2*rho*9.81*(1-A)*z3p/G*3*500*int_exp1(z3p, z4);
    output = output + rc^2*rho*9.81*(1-A)*z3p/G*15*500/2*int_exp3(z3p, z4);
    output = output - rc^2*rho*9.81*(1-A)*z3p/G*3*0.25*int_exp1(z3p, z4);
end

function output = tilt_tau(z1, z2, z3, Lp)
    global rho rc rs G mu;
    
    output = 3*rho*9.81*(rs^2/2/rc-rc/2)*500/4/3.14/G*(1-2*0.25)/(1-0.25)*2/5*int_exp3(z2, z3);
    output = output+3*2*mu*Lp/rc^2/4/3.14/G*(1-2*0.25)/(1-0.25)*int_exp5(z1, z2);
    output = output+3*2*mu*Lp/rc^2/4/3.14/G*(1-2*0.25)/(1-0.25)*2*500^2/5*int_exp3(z1, z2);
end

function output = int_exp1(x1, x2)
    output = -1/3/(x2^2+500^2)^(3/2)+1/3/(x1^2+500^2)^(3/2);
end 

function output = int_exp2(x1, x2)
    output = x2^3/(3*500^2*(x2^2+500^2)^(3/2))-x1^3/(3*500^2*(x1^2+500^2)^(3/2));
end 

function output = int_exp3(x1, x2)
    output = (-x2^2-100000)/3/(x2^2+500^2)^(5/2)-(-x1^2-100000)/3/(x1^2+500^2)^(5/2);
end 

function output = int_exp4(x1, x2)
    output = x2^5/5/500^2/(x2^2+500^2)^(5/2)-x1^5/5/500^2/(x1^2+500^2)^(5/2);
end 

function output = int_exp5(x1, x2)
    output = (-(3*x2^2+500000)/3/(x2^2+500^2)^(3/2))-(-(3*x1^2+500000)/3/(x1^2+500^2)^(3/2));
end
