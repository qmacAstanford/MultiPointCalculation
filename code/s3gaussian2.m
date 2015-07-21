function S3=s3gaussian2(N,NM,LAM,FA,Q1,Q2,Q3,d)

S3=zeros(2,2,2);
MIN=1e-5;

% Reset Qs to column vectors if entered as rows
if isrow(Q1)==1
    Q1=transpose(Q1);
    Q2=transpose(Q2);
    Q3=transpose(Q3);
end

if sum(power(Q1+Q2+Q3,2)) <= MIN
    
    % Evaluate the quantities for s3 calculation
    Q1MAG=sqrt(sum(power(Q1,2)));
    Q2MAG=sqrt(sum(power(Q2,2)));
    
    Q12MAG=sqrt(sum(power(Q1+Q2,2)));
    Q13MAG=sqrt(sum(power(Q1+Q3,2)));
    Q23MAG=sqrt(sum(power(Q2+Q3,2)));

    S2_1=s2gaussian(N,NM,FA,LAM,Q1MAG,d);
    S2_2=s2gaussian(N,NM,FA,LAM,Q2MAG,d);
    
    S2_12=s2gaussian(N,NM,FA,LAM,Q12MAG,d);
    S2_13=s2gaussian(N,NM,FA,LAM,Q13MAG,d);
    S2_23=s2gaussian(N,NM,FA,LAM,Q23MAG,d);
    
    for I1=1:2
            for I2=1:2
                    for I3=1:2
                            S3(I1,I2,I3)=S3(I1,I2,I3)+S2_1(I1,I2)*S2_12(I2,I3);
                            S3(I1,I2,I3)=S3(I1,I2,I3)+S2_1(I1,I2)*S2_13(I2,I3);
                            S3(I1,I2,I3)=S3(I1,I2,I3)+S2_2(I1,I2)*S2_23(I2,I3);
                    end
            end
    end
    S3 = S3/3;

end
end