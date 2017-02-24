%module testing Q and ISR_EF and ISR_WA for varying gaussian noise levels,
%G


%generating gaussian noise values

num = 300

Gau = zeros(num,1);
Gau(1,1) = 0.001;
for i=2:num
    Gau(i,1) = Gau(i-1,1) + 0.001;
end


ISREF = zeros(0,1);
ISRWA = zeros(0,1);
Q = zeros(0,1);
COR = zeros(0,1);
PC = zeros(0,1);

%running pipeline
for i=1:num
    load('DSTREAM_good')
    load('PRE_good')
    load('Amix_good')
    ISREFtmp = zeros(0,1);
    ISRWAtmp = zeros(0,1);
    Qtmp = zeros(0,1);
    CORtmp = zeros(0,1);
    PCtmp = zeros(0,1);
    
    for j=1:10
        i
        j
        %generating simulations 
        [X,PRE,A,PAR] = bss_simu('A',A,'G',Gau(i,1));
        %running pipeline
        bss_pipeline
        
        Qtmp = cat(1,Qtmp,max(max(DSTREAM.INFO.FIND.LBQstats)));
        CORtmp = cat(1,CORtmp,max(max(DSTREAM.INFO.FIND.corrstats)));
        ISREFtmp = cat(1,ISREFtmp,mean(mean(DSTREAM.INFO.MCOMBI.ISRef1)));
        ISRWAtmp = cat(1,ISRWAtmp,mean(mean(DSTREAM.INFO.MCOMBI.ISRwa1)));
        PCtmp = cat(1,PCtmp,max(DSTREAM.DATA.PCA.latent));
    end
  
    %extracting parameters
    Q = cat(1,Q,mean(Qtmp));
    COR = cat(1,COR,mean(CORtmp));
    ISREF = cat(1,ISREF,mean(ISREFtmp));
    ISRWA = cat(1,ISRWA,mean(ISRWAtmp));
    PC = cat(1,PC,mean(PCtmp));
end


    
    

