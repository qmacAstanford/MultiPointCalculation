function [A,B,C] = offset(a,b,c,MIN)
% Insures that a,b,c are separated in value by at least MIN
% If they are too close it moves them appart
% Maintains the order of a,b,c and the average of them
    out=[a,b,c];
    [Amin,Imin]=min(out);
    [Amax,Imax]=max(out);
    if (Amax-Amin)<2*MIN
        avj=(a+b+c)/3.0;
        out=[avj,avj,avj];
        out(Imin)=avj-MIN;
        out(Imax)=avj+MIN;
    elseif abs(a-b)<MIN
        avj=(a+b)/2.0;
        if a > b
            out=[avj+0.5*MIN,avj-0.5*MIN,out(3)];
        else
            out=[avj-0.5*MIN,avj+0.5*MIN,out(3)];
        end
    elseif abs(b-c)<MIN
        avj=(b+c)/2.0;
        if b>c
            out=[out(1),0.5*MIN+avj,-0.5*MIN+avj];
        else
            out=[out(1),-0.5*MIN+avj,0.5*MIN+avj];
        end
    elseif abs(a-c)<MIN
        avj=(a+c)/MIN;
        if a<c
            out=[-0.5*MIN+avj,out(2),0.5*MIN+avj];
        else
            out=[0.5*MIN+avj,out(2),-0.5*MIN+avj];
        end
    end
    A=out(1);
    B=out(2);
    C=out(3);
end