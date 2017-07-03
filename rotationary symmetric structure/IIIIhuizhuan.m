E=2.1e11;
A=1e-4;
density=7.3e3; 
node_number=30;
elment_number=54;
asd = sqrt(3);
 
nc=[0.0,      2.0;
            asd,      1.0;
            asd,      3.0;
            2*asd,      2.0;
            2*asd,      4.0;
            3*asd,      3.0;
            asd,      -1.0;
            2*asd,      0.0;
           2*asd,      -2.0;
           3*asd,      -1.0;
           3*asd,      -3.0;
           0.0,     -2.0;
           asd,      -3.0;
           0.0,      -4.0;
           asd,      -5.0;
           0.0,      -6.0;
           -asd,      -1.0;
           -asd,      -3.0;
           -2*asd,      -2.0;
           -2*asd,      -4.0;
           -3*asd,      -3.0;
          -asd,   1.0;
          -2*asd,    0.0;
          -2*asd,    2.0;
          -3*asd,    1.0;
          -3*asd,    3.0;
          -asd,    3.0;
          0.0,    4.0;
          -asd,    5.0;
          0.0,     6.0];

 
en=[1,          2;
                  1,          3;
                  2,          3;
                  2,          4;
                  3,          4;
                  3,         5;
                  4,          5;
                  4,          6;
                  5,         6;
                2,         7;
                2,         8;
          7,   8;
          7,   9;
          8,   9;
          8,   10;
          9,   10;
          9,   11;
          10,   11;
          7,   12;
          7,   13;
          12,   13;
          12,   14;
          13,   14;
          13,   15;
          14,   15;
          14,  16;
          15,  16;
          12,  17;
          12,   18;
          17,   18;
          17,   19;
          18,   19;
          18,   20;
          19,   20;
          19,  21;
          20,  21;
          17,  22;
          17,  23;
          22,  23;
          22,  24;
          23,  24;
          23,  25;
          24,  25;
          24,  26;
          25,  26;
          1,  22;
          22,  27;
          1,  27;
          1,  28;
          27,  28;
          27,  29;
          28,  29;
          28,  30;
          29,  30];

ed(1:node_number,1:2)=1;

constraint=[6,1;6,2;11,1;11,2;16,1;16,2;21,1;21,2;26,1;26,2;30,1;30,2];
 
for loopi=1:length(constraint);
    ed(constraint(loopi,1),constraint(loopi,2))=0;
end
dof=0;
for loopi=1:node_number
     for loopj=1:2
        if ed(loopi,loopj)~=0
            dof=dof+1;
            ed(loopi,loopj)=dof;
        end
    end
end
 


k = kk1;
m = mm1;
[v,d] = eig(k,m); 
 
frequency=sqrt(diag(d))/(2*pi);

[frequency,indexf]=sort(frequency);
d=d(:,indexf);