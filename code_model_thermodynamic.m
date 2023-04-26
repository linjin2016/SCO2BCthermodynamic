%This code is the thermodynamic model used in the article "Thermodynamic–economic evaluation and multi-objective optimization of supercritical carbon dioxide recompression Brayton cycle considering printed circuit recuperator types and designs
%This code can be available online when the article i
clear
clc
%known conditions
intMC=0.89;%Isentropic efficiency of main compressor
intRC=0.89;%Isentropic efficiency of recompressor
intTB=0.93;%Isentropic efficiency of turbine
P1=7615;%Inlet pressure of main compressor 
P2=25000;%Outlet pressure of main compressor
T1=32+273;%Inlet temmerature of main compressor,K
S1=refpropm('S','T',T1,'P',P1,'co2');%Inlet entropy of main compressor J/kg K
Han1=refpropm('H','T',T1,'P',P1,'co2');%Inlet enthalpy of main compressor J/kg
Han2isen=refpropm('H','P',P2,'S',S1,'co2');%isentropic exit enthalpy of main compression , isentropic compression condition K
Wcomp=(Han2isen-Han1)/intMC;
Han2=Han1+Wcomp;
T2=refpropm('T','P',P2,'H',Han2,'co2');

T6=700+273;%Inlet temmerature of turbine，K

FRC=0.3682; % Recompression fraction: variable
ksaiLTR=0.95; % Enthalpy efficiency of low temperature recuperator: variable
ksaiHTR=0.95; % Enthalpy efficiency of high temperature recuperator: variable
GLTR=330; % Mass flow flux of low temperature recuperator
GHTR=330; % Mass flow flux of high temperature recuperator

%initial value definition 
P5=P2-500;
P9=P1+77;
P8=P9+30;
P3=P2;P4=P3;P10=P3;
P6=P5-245;

P7=P1+200;


T8=470; %  T8 initial value definition 

%low temperature recuperator zigzag channel
NumLTR=401;
twLTR=0.0005;
conwLTR=16.3; 
DiamLTR=2.75;
Jh=12;Jc=12;
DUh=36;DUc=36; %bending angle

ALTR=(pi*(DiamLTR/2)^2)/2;% cross-sectional area，mm2
DLTR=pi*DiamLTR/2/(1+pi/2);% hydraulic diameter,mm
AzhouchangLTR=pi*DiamLTR/2+DiamLTR;% perimeter,mm

mhLTR=GLTR*ALTR/1000000;% mass flow rate of channel
mcLTR=mhLTR*(1-FRC);% mass flow rate of channel

dPHLTR=zeros(1,NumLTR);
dPCLTR=zeros(1,NumLTR);
HanhLTR=zeros(1,NumLTR);
HancLTR=zeros(1,NumLTR);
LengthLTR=dPCLTR;
PHLTR=zeros(1,NumLTR+1);
PCLTR=zeros(1,NumLTR+1);
THLTR=PHLTR;TCLTR=PHLTR;
LTRaveh=dPCLTR;

THLTR(NumLTR+1)=T8;
TCLTR(1)=T2;
PHLTR(NumLTR+1)=P8;
PHLTR(1)=P9;
PCLTR(1)=P2;
PCLTR(NumLTR+1)=P3;

%high temperature recuperator airfoil channel
NumHTR=401;
twHTR=0.0005;
conwHTR=16.3; 

Jzhyi=12;Jzcyi=12; % Periodic longitudinal pitch mm
Wchyi=4.2;Wccyi=4.2; % Channel widch mm
Hhyi=0.95; Hcyi=0.95; % Channel height mm

Sa=4.94; %The area of the top surface of the airfoil mm2
Vh=(Jzhyi*Wchyi-Sa)*Hhyi; %The volume of the repeating unit
Pa=12.57; %The side length of the airfoil mm
S=Pa*Hhyi+2*(Jzhyi*Wchyi-Sa); %The wetting area of the repeating unit

AHTR=Wchyi*Hhyi;% cross-sectional area
DHTR=4*Vh/S;
AzhouchangHTR=S/Jzhyi;

mhHTR=GHTR*AHTR/1000000;
mcHTR=mhHTR;

dPHHTR=zeros(1,NumHTR);
dPCHTR=zeros(1,NumHTR);
HanhHTR=zeros(1,NumHTR);
HancHTR=zeros(1,NumHTR);
LengthHTR=dPCHTR;
PHHTR=zeros(1,NumHTR+1);
PCHTR=zeros(1,NumHTR+1);
THHTR=PHHTR;TCHTR=PHHTR;
HTRaveh=dPCHTR;

    
%Initial value definition of convergence
mLTR=1;
nLTR=1;
pLTR=1;
mT8star=1;
nP5star=1;
nP7star=1;

    while abs((nP5star-P5)/P5)>0.00005|abs((nP7star-P7)/P7)>0.00005
     
    
     
     S6=refpropm('S','T',T6,'P',P6,'co2');
     Han6=refpropm('H','T',T6,'P',P6,'co2');
     Han7isen=refpropm('H','P',P7,'S',S6,'co2');
     Wturb=intTB*(Han6-Han7isen);
     Han7=Han6-Wturb;
     T7=refpropm('T','P',P7,'H',Han7,'co2');
     
    while abs((mT8star-T8)/T8)>0.00001
         
     THLTR(NumLTR+1)=T8;
     TCLTR(1)=T2;
     PHLTR(NumLTR+1)=P8;
     PHLTR(1)=P9;
     PCLTR(1)=P2;
     PCLTR(NumLTR+1)=P3;
        
     HanThinPHLTR=refpropm('H','T',THLTR(NumLTR+1),'P',PHLTR(NumLTR+1),'co2');
     HanTcinPHLTR=refpropm('H','T',TCLTR(1),'P',PHLTR(1),'co2');
     HanThinPCLTR=refpropm('H','T',THLTR(NumLTR+1),'P',PCLTR(NumLTR+1),'co2');
     HanTcinPCLTR=refpropm('H','T',TCLTR(1),'P',PCLTR(1),'co2');
     delHmaxLTR=min(mhLTR*(1-FRC)*(HanThinPCLTR-HanTcinPCLTR),mhLTR*(HanThinPHLTR-HanTcinPHLTR));
     QLTR=ksaiLTR*delHmaxLTR;     
     Han8=refpropm('H','T',T8,'P',P8,'co2');
     Han9=Han8-QLTR/mhLTR;
     T9=refpropm('T','P',P9,'H',Han9,'co2');
     THLTR(1)=T9;
   
    for i=1:NumLTR 
    %Thermal physical properties of hot side 
    CphLTR=refpropm('C','T',THLTR(i),'P',PHLTR(i),'co2'); 
    MiuhLTR=refpropm('V','T',THLTR(i),'P',PHLTR(i),'co2');
    DenhLTR=refpropm('D','T',THLTR(i),'P',PHLTR(i),'co2'); 
    ConhLTR=refpropm('L','T',THLTR(i),'P',PHLTR(i),'co2'); 
    RehLTR=1000*mhLTR*DLTR/(MiuhLTR*ALTR); 
    PrhLTR=refpropm('^','T',THLTR(i),'P',PHLTR(i),'co2'); 
    HhLTR=refpropm('H','T',THLTR(i),'P',PHLTR(i),'co2'); 
    HanhLTR(i)=HhLTR;
    %Thermal physical properties of cold side
    CpcLTR=refpropm('C','T',TCLTR(i),'P',PCLTR(i),'co2'); 
    MiucLTR=refpropm('V','T',TCLTR(i),'P',PCLTR(i),'co2');
    DencLTR=refpropm('D','T',TCLTR(i),'P',PCLTR(i),'co2'); 
    ConcLTR=refpropm('L','T',TCLTR(i),'P',PCLTR(i),'co2'); 
    RecLTR=1000*mcLTR*DLTR/(MiucLTR*ALTR); 
    PrcLTR=refpropm('^','T',TCLTR(i),'P',PCLTR(i),'co2');
    HcLTR=refpropm('H','T',TCLTR(i),'P',PCLTR(i),'co2');  
    HancLTR(i)=HcLTR;
    
    %Nu and f of hot side
    fhLTR=0.1924*(RehLTR^(-0.091)); 
    NuhLTR=0.184*(RehLTR^0.629)*(PrhLTR^0.317);
    %Nu and f of cold side
    fcLTR=0.1924*(RecLTR^(-0.091)); 
    NucLTR=0.184*(RecLTR^0.629)*(PrcLTR^0.317);
    % heat transfer coefficient
    hhLTR=1000*NuhLTR*ConhLTR/DLTR; 
    hcLTR=1000*NucLTR*ConcLTR/DLTR; 
    hallLTR=1/(1/hhLTR+1/hcLTR+twLTR/conwLTR); 
    LTRaveh(i)=hallLTR; 
    
    delTQLTR=QLTR/NumLTR;%J
    delTLLTR=delTQLTR/hallLTR/(THLTR(i)-TCLTR(i))/(AzhouchangLTR/1000);%m
    LengthLTR(i)=delTLLTR;% m
    
    %pressure drop of hot side 
    deltaPhLTR=fhLTR*2*(LengthLTR(i))*mhLTR^2/(DenhLTR*(DLTR/1000)*(ALTR/1000000)^2)/1000; %根据f公式计算节点间压�?�?��除以1000是转化为kpa
    dPHLTR(i)=deltaPhLTR; 
    PHLTR(i+1)=PHLTR(i)+deltaPhLTR; 
    % pressure drop of cold side
    deltaPcLTR=fcLTR*2*LengthLTR(i)*i*mcLTR^2/(DencLTR*(DLTR/1000)*(ALTR/1000000)^2)/1000;%根据f公式计算节点间压�?,�?��除以1000是转化为kpa
    dPCLTR(i)=deltaPcLTR; 
    PCLTR(i+1)=PCLTR(1)-deltaPcLTR; 
    
    % Calculate node parameters according to enthalpy value
    AFLTR=(AzhouchangLTR/1000)*LengthLTR(i);%m2
    HanhLTR(i+1)=HanhLTR(i)+hallLTR*AFLTR*(THLTR(i)-TCLTR(i))/mhLTR;
    HancLTR(i+1)=HancLTR(i)+hallLTR*AFLTR*(THLTR(i)-TCLTR(i))/mcLTR;
    THLTR(i+1)=THLTR(i)+(HanhLTR(i+1)-HanhLTR(i))/CphLTR;
    TCLTR(i+1)=TCLTR(i)+(HancLTR(i+1)-HancLTR(i))/CpcLTR;  
    end
    
    
    P8=PHLTR(NumLTR+1);
 
    T3=TCLTR(NumLTR+1);
    P3=PCLTR(NumLTR+1);
    
    % T10
   P10=P3;
   S9=refpropm('S','T',T9,'P',P9,'co2'); 
   Han9=refpropm('H','T',T9,'P',P9,'co2'); 
   Han10isen=refpropm('H','P',P10,'S',S9,'co2'); 
   Wrcomp=(Han10isen-Han9)/intRC;
   Han10=Han9+Wrcomp;
   T10=refpropm('T','P',P10,'H',Han10,'co2');
   
   % T4
   P4=P3;
   Han3=refpropm('H','T',T3,'P',P3,'co2');
   Han4=(1-FRC)*Han3+FRC*Han10;
   T4=refpropm('T','P',P4,'H',Han4,'co2');
   
   THHTR(NumHTR+1)=T7;
   TCHTR(1)=T4;
   PHHTR(1)=P8;
   PCHTR(1)=P4;
   PCHTR(NumHTR+1)=P5;
   PHHTR(NumHTR+1)=P7;
   
   HanThinPHHTR=refpropm('H','T',THHTR(NumHTR+1),'P',PHHTR(NumHTR+1),'co2');
   HanTcinPHHTR=refpropm('H','T',TCHTR(1),'P',PHHTR(1),'co2');
   HanThinPCHTR=refpropm('H','T',THHTR(NumHTR+1),'P',PCHTR(NumHTR+1),'co2');
   HanTcinPCHTR=refpropm('H','T',TCHTR(1),'P',PCHTR(1),'co2');
   
   delHmaxHTR=min(mcHTR*(HanThinPCHTR-HanTcinPCHTR),mhHTR*(HanThinPHHTR-HanTcinPHHTR));
   QHTR=ksaiHTR*delHmaxHTR;     
   Han8=Han7-QHTR/mhHTR;
   T8star=refpropm('T','P',PHHTR(1),'H',Han8,'co2');
   
   mT8star=T8star;
   if abs((mT8star-T8)/T8)>0.00001
    aT8=T8star-T8;
    dqT8=abs(0.5*aT8);
    if aT8>0
        dT8=T8+dqT8;
    elseif aT8<0
        dT8=T8-dqT8;
    end
    T8=dT8;
   end
   
    end
    
    PHHTR(1)=P8;
    THHTR(1)=T8;
    PCHTR(1)=P4;
    TCHTR(1)=T4;
    
    % High temperature recuperator
    for i=1:NumHTR     
    %Thermal physical properties of hot side
    CphHTR=refpropm('C','T',THHTR(i),'P',PHHTR(i),'co2');
    MiuhHTR=refpropm('V','T',THHTR(i),'P',PHHTR(i),'co2'); 
    DenhHTR=refpropm('D','T',THHTR(i),'P',PHHTR(i),'co2'); 
    ConhHTR=refpropm('L','T',THHTR(i),'P',PHHTR(i),'co2'); 
    RehHTR=1000*mhHTR*DHTR/(MiuhHTR*AHTR); 
    PrhHTR=refpropm('^','T',THHTR(i),'P',PHHTR(i),'co2'); 
    HhHTR=refpropm('H','T',THHTR(i),'P',PHHTR(i),'co2'); 
    HanhHTR(i)=HhHTR;
    %Thermal physical properties of cold side
    CpcHTR=refpropm('C','T',TCHTR(i),'P',PCHTR(i),'co2'); 
    MiucHTR=refpropm('V','T',TCHTR(i),'P',PCHTR(i),'co2'); 
    DencHTR=refpropm('D','T',TCHTR(i),'P',PCHTR(i),'co2'); 
    ConcHTR=refpropm('L','T',TCHTR(i),'P',PCHTR(i),'co2'); 
    RecHTR=1000*mcHTR*DHTR/(MiucHTR*AHTR); 
    PrcHTR=refpropm('^','T',TCHTR(i),'P',PCHTR(i),'co2'); 
    HcHTR=refpropm('H','T',TCHTR(i),'P',PCHTR(i),'co2');  
    HancHTR(i)=HcHTR;
    % Nu and f of hot side 
    fhHTR=0.357*(1.2/(Wchyi/2))^(-0.255)*(RehHTR)^((-0.173)*(1.2/(Wchyi/2))^(-0.274));
    NuhHTR=0.0318*(RehHTR)^(0.78)*PrhHTR^(0.4);
    %Nu and f of cold side
    fcHTR=0.357*(1.2/(Wchyi/2))^(-0.255)*(RecHTR)^((-0.173)*(1.2/(Wchyi/2))^(-0.274)); %冷侧f
    NucHTR=0.0318*(RecHTR)^(0.78)*PrcHTR^(0.4);
    
    % Heat Transfer Coefficient
    hhHTR=1000*NuhHTR*ConhHTR/DHTR; 
    hcHTR=1000*NucHTR*ConcHTR/DHTR;
    hallHTR=1/(1/hhHTR+1/hcHTR+twHTR/conwHTR); 
    HTRaveh(i)=hallHTR; 
    
    delTQHTR=QHTR/NumHTR;%J
    delTLHTR=delTQHTR/hallHTR/(THHTR(i)-TCHTR(i))/(AzhouchangHTR/1000);%m
    LengthHTR(i)=delTLHTR;% m
    
    % pressure drop of hot side
    deltaPhHTR=fhHTR/2*LengthHTR(i)*mhHTR^2/(DenhHTR*(DHTR/1000)*(AHTR/1000000)^2)/1000; %根据f公式计算节点间压�?    dPHHTR(i)=deltaPhHTR; 
    PHHTR(i+1)=PHHTR(i)+deltaPhHTR; 
    % pressure drop of cold side
    deltaPcHTR=fcHTR/2*LengthHTR(i)*i*mcHTR^2/(DencHTR*DHTR/1000*(AHTR/1000000)^2)/1000;%根据f公式计算节点间压�?    dPCHTR(i)=deltaPcHTR; 
    PCHTR(i+1)=PCHTR(1)-deltaPcHTR; 
    
    
    AFHTR=(AzhouchangHTR/1000)*LengthHTR(i);
    HanhHTR(i+1)=HanhHTR(i)+hallHTR*AFHTR*(THHTR(i)-TCHTR(i))/mhHTR;
    HancHTR(i+1)=HancHTR(i)+hallHTR*AFHTR*(THHTR(i)-TCHTR(i))/mcHTR;
    THHTR(i+1)=THHTR(i)+(HanhHTR(i+1)-HanhHTR(i))/CphHTR;
    TCHTR(i+1)=TCHTR(i)+(HancHTR(i+1)-HancHTR(i))/CpcHTR;  
    end
    
    
    P7star=PHHTR(NumHTR+1);
    P5star=PCHTR(NumHTR+1);
    T5=TCHTR(NumHTR+1);
    T7=THHTR(NumHTR+1);
    
    nP5star=P5star;
    nP7star=P7star;
    
    if abs((nP5star-P5)/P5)>0.00005
    P5=nP5star;
    P6=P5-250;
    end
    
    if abs((nP7star-P7)/P7)>0.00005
    P7=nP7star;
    end
    
    end
    
   T0=293.15;% environment temperature 
   P0=101.33;%  environment pressure
   Han0=refpropm('H','T',T0,'P',P0,'co2');   
   Han1=refpropm('H','T',T1,'P',P1,'co2');
   Han2=refpropm('H','T',T2,'P',P2,'co2');
   Han3=refpropm('H','T',T3,'P',P3,'co2');
   Han4=refpropm('H','T',T4,'P',P4,'co2');
    
   Han7=refpropm('H','T',T7,'P',P7,'co2');
   Han8=refpropm('H','T',T8,'P',P8,'co2');
   Han10=refpropm('H','T',T10,'P',P10,'co2');
    
   Han6=refpropm('H','T',T6,'P',P6,'co2');
   Han5=refpropm('H','T',T5,'P',P5,'co2');
    
   mCO2=74.968;
   QPHX=mCO2*(Han6-Han5);
    
   
   Han9=refpropm('H','T',T9,'P',P9,'co2');
   QPC=mCO2*(1-FRC)*(Han9-Han1);
   
   intTH=1-(QPC/QPHX);
   
   % cost calculation
   % compressor and turbine
   WTturb=mCO2*Wturb; 
   WTcomp=mCO2*(1-FRC)*Wcomp+mCO2*FRC*Wrcomp; 
   Costcomp=6898*(WTcomp/1000)^(0.7865); 
   Costturb=7790*(WTturb/1000)^(0.6842); 
   % heat and cold source 
   UAPHE=QPHX/(T6-T5);
   UAPC=QPC/(T9-T1);
   CostPHE=3500*(UAPHE/1000);
   CostPC=2300*(UAPC/1000);
   % recuperator
   % low temperature recuperator
   LengthL=sum(LengthLTR(:));
   Wch=3.3;
   Hch=4;
   ncnLTR=mCO2/mhLTR;
   VolumeLTR=Wch*Hch*LengthL*ncnLTR/1000000;
   CostLTR=7940*120*VolumeLTR;
   % high temperature recuperator
   LengthH=sum(LengthHTR(:));
   Wfyi=4.2;
   tpyi=3.26;
   ncnHTR=mCO2/mhHTR;
   VolumeHTR=Wfyi*tpyi*LengthH*ncnHTR/1000000;
   CostHTR=7940*120*VolumeHTR;
   CostRecup=CostLTR+CostHTR;
   % total cost 
   CostTotol=CostRecup+CostPC+CostPHE+Costcomp+Costturb;
  
   % Exergy efficiency
   Shang0=refpropm('S','T',T0,'P',P0,'co2');
   Shang1=refpropm('S','T',T1,'P',P1,'co2');
   Shang2=refpropm('S','T',T2,'P',P2,'co2');
   Shang3=refpropm('S','T',T3,'P',P3,'co2');
   Shang4=refpropm('S','T',T4,'P',P4,'co2');
    
   Shang7=refpropm('S','T',T7,'P',P7,'co2');
   Shang8=refpropm('S','T',T8,'P',P8,'co2');
   Shang9=refpropm('S','T',T9,'P',P9,'co2');
   Shang10=refpropm('S','T',T10,'P',P10,'co2');
    
   Shang6=refpropm('S','T',T6,'P',P6,'co2');
   Shang5=refpropm('S','T',T5,'P',P5,'co2');
     
   Wnet=WTturb-WTcomp; 
   
   Ex1=((Han1-Han0)-T0*(Shang1-Shang0));%1 point
   Ex2=((Han2-Han0)-T0*(Shang2-Shang0));%2 point
   Ex3=((Han3-Han0)-T0*(Shang3-Shang0));%3 point
   Ex4=((Han4-Han0)-T0*(Shang4-Shang0));%4 point
   Ex5=((Han5-Han0)-T0*(Shang5-Shang0));%5 point
   Ex6=((Han6-Han0)-T0*(Shang6-Shang0));%6 point
   Ex7=((Han7-Han0)-T0*(Shang7-Shang0));%7 point
   Ex8=((Han8-Han0)-T0*(Shang8-Shang0));%8 point
   Ex9=((Han9-Han0)-T0*(Shang9-Shang0));%9 point
   Ex10=((Han10-Han0)-T0*(Shang10-Shang0));%10 point
   
   Exsunturbine=mCO2*(Ex6-Ex7)-WTturb; 
   ExsunHTR=mCO2*(Ex4+Ex7-Ex5-Ex8); 
   ExsunLTR=mCO2*(1-FRC)*Ex2+mCO2*Ex8-mCO2*(1-FRC)*Ex3-mCO2*Ex9; 
   ExsunmainCom=mCO2*(1-FRC)*Wcomp-mCO2*(1-FRC)*(Ex2-Ex1); 
   ExsunreCom=mCO2*FRC*Wrcomp-mCO2*FRC*(Ex10-Ex9); 
   ExsunPHE=(1-1/0.8)*(Ex5-Ex6);%PHE
   ExsunPC=(1-0.5)*(Ex9-Ex1);%PC
   Exsunbujian=Exsunturbine+ExsunHTR+ExsunLTR+ExsunmainCom+ExsunreCom+ExsunPC+ExsunPHE;% total
   AME=[intTH*100,CostTotol/1000000,Exsunbujian/1000000];