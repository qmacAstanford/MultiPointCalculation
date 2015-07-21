function S4=s4gaussian2(N,NM,LAM,FA,Q1,Q2,Q3,Q4,d)

S4=zeros(2,2,2,2);
MIN=1e-3;

% Reset Qs to column vectors if entered as rows
if isrow(Q1)==1
    Q1=transpose(Q1);
    Q2=transpose(Q2);
    Q3=transpose(Q3);
    Q4=transpose(Q4);
end

% Begin calculation of s4
if sum(power(Q1+Q2+Q3+Q4,2)) <= MIN

    % Evaluate the quantities for s4 calculation
    Q1MAG=sqrt(sum(power(Q1,2)));
    Q2MAG=sqrt(sum(power(Q2,2)));
    Q3MAG=sqrt(sum(power(Q3,2)));
    Q4MAG=sqrt(sum(power(Q4,2)));
    
    Q12MAG=sqrt(sum(power(Q1+Q2,2)));
    Q13MAG=sqrt(sum(power(Q1+Q3,2)));
    Q14MAG=sqrt(sum(power(Q1+Q4,2)));
    Q23MAG=sqrt(sum(power(Q2+Q3,2)));
    Q24MAG=sqrt(sum(power(Q2+Q4,2)));
    Q34MAG=sqrt(sum(power(Q3+Q4,2)));
    
    S2_1=s2gaussian(N,NM,FA,LAM,Q1MAG,d);
    S2_2=s2gaussian(N,NM,FA,LAM,Q2MAG,d);
    S2_3=s2gaussian(N,NM,FA,LAM,Q3MAG,d);
    S2_4=s2gaussian(N,NM,FA,LAM,Q4MAG,d);
    
    S2_12=s2gaussian(N,NM,FA,LAM,Q12MAG,d);
    S2_13=s2gaussian(N,NM,FA,LAM,Q13MAG,d);
    S2_14=s2gaussian(N,NM,FA,LAM,Q14MAG,d);
    S2_23=s2gaussian(N,NM,FA,LAM,Q23MAG,d);
    S2_24=s2gaussian(N,NM,FA,LAM,Q24MAG,d);
    S2_34=s2gaussian(N,NM,FA,LAM,Q34MAG,d);
    
    for I1=1:2
            for I2=1:2
                    for I3=1:2
                            for I4=1:2
                                
                                % PERMUTATIONS
                                % 1 12 4 112 443
                                % 1 12 3 112 334
                                % 1 13 4 113 442
                                % 1 13 2 113 224
                                % 1 14 3 114 332
                                % 1 14 2 114 223
                                % 2 23 4 223 441
                                % 2 23 1 223 114
                                % 2 24 3 224 331
                                % 2 24 1 224 113
                                % 3 34 2 334 221
                                % 3 34 1 334 112
    
                                S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+S2_1(I1,I2)*S2_12(I2,I3)*S2_4(I3,I4);
                                S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+S2_1(I1,I2)*S2_12(I2,I3)*S2_3(I3,I4);
                                S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+S2_1(I1,I2)*S2_13(I2,I3)*S2_4(I3,I4);
                                S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+S2_1(I1,I2)*S2_13(I2,I3)*S2_2(I3,I4);
                                S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+S2_1(I1,I2)*S2_14(I2,I3)*S2_3(I3,I4);
                                S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+S2_1(I1,I2)*S2_14(I2,I3)*S2_2(I3,I4);
                                S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+S2_2(I1,I2)*S2_23(I2,I3)*S2_4(I3,I4);
                                S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+S2_2(I1,I2)*S2_23(I2,I3)*S2_1(I3,I4);
                                S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+S2_2(I1,I2)*S2_24(I2,I3)*S2_3(I3,I4);
                                S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+S2_2(I1,I2)*S2_24(I2,I3)*S2_1(I3,I4);
                                S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+S2_3(I1,I2)*S2_34(I2,I3)*S2_2(I3,I4);
                                S4(I1,I2,I3,I4)=S4(I1,I2,I3,I4)+S2_3(I1,I2)*S2_34(I2,I3)*S2_1(I3,I4);
    
                            end
                    end
            end
    end
    S4 = S4/12;
    
end
end
