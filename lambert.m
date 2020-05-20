%% Lambert
function [v1, v2] = lambert(r1, r2, dt)
global mue
R1 = norm(r1);
R2 = norm(r2);
deltat = dt;
rcross = cross(r1,r2);
if rcross(3) >= 0
    dtheta = acosd(dot(r1,r2)/(R1*R2));
else
    dtheta = 360 - acosd(dot(r1,r2)/(R1*R2));
end
A = sind(dtheta)*sqrt(R1*R2/(1-cosd(dtheta)));
z = -100;
while F(z,deltat) < 0
    z = z + 0.1;
end
error = 10^-8;
ratio = 1;
i = 0;
while abs(ratio)>error
    i = i + 1;
    ratio = F(z,deltat)/dFdz(z);
    z = z - ratio;
    if i == 10000
        break
    end
end
f = 1 - y(z)/R1;
g = A*sqrt(y(z)/mue);
gdot = 1 - y(z)/R2;
v1 = (1/g)*(r2-f*r1);
v2 = (1/g)*(gdot*r2-r1);
    function dum = y(z)
        dum = R1 + R2 + A*(z*S(z) - 1)/sqrt(C(z));
    end
    function dum = F(z,t)
        dum = (y(z)/C(z))^1.5*S(z) + A*sqrt(y(z)) - sqrt(mue)*t;
    end
    function dum = dFdz(z)
        if z == 0
            dum = sqrt(2)/40*y(0)^1.5 + A/8*(sqrt(y(0)) + A*sqrt(1/2/y(0)));
        else
            dum = (y(z)/C(z))^1.5*(1/2/z*(C(z) - 3*S(z)/2/C(z)) ...
                + 3*S(z)^2/4/C(z)) + A/8*(3*S(z)/C(z)*sqrt(y(z)) ...
                + A*sqrt(C(z)/y(z)));
        end
    end
    function dum = C(z)
        dum = stumpC(z);
    end
    function c = stumpC(z)
        % wwwwwwwwwwwwwwwwwwwwww
        %{
This function evaluates the Stumpff function C(z) according
to Equation 3.53.
z - input argument
c - value of C(z)
User M-functions required: none
        %}
        % ----------------------------------------------
        if z > 0
            c = (1 - cos(sqrt(z)))/z;
        elseif z < 0
            c = (cosh(sqrt(-z)) - 1)/(-z);
        else
            c = 1/2;
        end
    end
    function dum = S(z)
        dum = stumpS(z);
    end
    function s = stumpS(z)
        % wwwwwwwwwwwwwwwwwwwwww
        %{
This function evaluates the Stumpff function S(z) according
to Equation 3.52.
z - input argument
s - value of S(z)
User M-functions required: none
        %}
        % ----------------------------------------------
        if z > 0
            s = (sqrt(z) - sin(sqrt(z)))/(sqrt(z))^3;
        elseif z < 0
            s = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z))^3;
        else
            s = 1/6;
        end
    end
end