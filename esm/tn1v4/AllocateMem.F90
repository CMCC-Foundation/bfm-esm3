!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! SUBROUTINE: AllocateMem
!
! DESCRIPTION
!   Allocation of memory for Global State variables and other Variables
!
! COPYING
!
!   Copyright (C) 2023 BFM System Team (bfm_st@cmcc.it)
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation.
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU General Public License for more details.
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!
! INCLUDE
#include"cppdefs.h"

  subroutine AllocateMem

! USES
  use global_mem
  use mem
  use mem_Param

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Implicit typing is never allowed
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  IMPLICIT NONE

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Local Variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  integer:: status
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  !Start the allocation of pelagic state global
  ! matrix and pointers
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  LEVEL1 "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-"
  LEVEL1 "# Allocating State Variables and Rates array ..."

#ifndef NOT_STANDALONE
    
   allocate(D3STATE(1:NO_BOXES,1:NO_D3_BOX_STATES),stat=status)
   if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D3STATE")
   D3STATE = ZERO
     allocate(D3SOURCE(1:NO_BOXES,1:NO_D3_BOX_STATES),stat=status)
     if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D3SOURCE")
     D3SOURCE = ZERO
   allocate(D3STATETYPE(1:NO_D3_BOX_STATES ),stat=status)
   if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem","D3STATETYPE")
   D3STATETYPE = ZERO

    
   allocate(D3DIAGNOS(1:NO_BOXES,1:NO_D3_BOX_DIAGNOSS),stat=status)
   if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D3DIAGNOS")
   D3DIAGNOS = ZERO

    
   allocate(D2DIAGNOS(1:NO_BOXES_XY,1:NO_D2_BOX_DIAGNOSS),stat=status)
   if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D2DIAGNOS")
   D2DIAGNOS = ZERO

#ifdef INCLUDE_SEAICE


#endif
    
   allocate(D2STATE_BEN(1:NO_BOXES_XY,1:NO_D2_BOX_STATES_BEN),stat=status)
   if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D2STATE_BEN")
   D2STATE_BEN = ZERO
     allocate(D2SOURCE_BEN(1:NO_BOXES_XY,1:NO_D2_BOX_STATES_BEN),stat=status)
     if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D2SOURCE_BEN")
     D2SOURCE_BEN = ZERO
   allocate(D2STATETYPE_BEN(1:NO_D2_BOX_STATES_BEN ),stat=status)
   if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem","D2STATETYPE_BEN")
   D2STATETYPE_BEN = ZERO

    
   allocate(D2DIAGNOS_BEN(1:NO_BOXES_XY,1:NO_D2_BOX_DIAGNOSS_BEN),stat=status)
   if (status /= 0) call error_msg_prn(ALLOC,"AllocateMem", "D2DIAGNOS_BEN")
   D2DIAGNOS_BEN = ZERO

#endif



#ifdef INCLUDE_SEAICE


#endif



  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocation of Pelagic variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

O2o => D3STATE(:,ppO2o); O2o=ZERO
N1p => D3STATE(:,ppN1p); N1p=ZERO
N3n => D3STATE(:,ppN3n); N3n=ZERO
N4n => D3STATE(:,ppN4n); N4n=ZERO
N5s => D3STATE(:,ppN5s); N5s=ZERO
N7f => D3STATE(:,ppN7f); N7f=ZERO
B1c => D3STATE(:,ppB1c); B1c=ZERO
P1c => D3STATE(:,ppP1c); P1c=ZERO
P1n => D3STATE(:,ppP1n); P1n=ZERO
P1p => D3STATE(:,ppP1p); P1p=ZERO
P1l => D3STATE(:,ppP1l); P1l=ZERO
P1s => D3STATE(:,ppP1s); P1s=ZERO
P1f => D3STATE(:,ppP1f); P1f=ZERO
P2c => D3STATE(:,ppP2c); P2c=ZERO
P2n => D3STATE(:,ppP2n); P2n=ZERO
P2p => D3STATE(:,ppP2p); P2p=ZERO
P2l => D3STATE(:,ppP2l); P2l=ZERO
P2f => D3STATE(:,ppP2f); P2f=ZERO
Z4c => D3STATE(:,ppZ4c); Z4c=ZERO
Z5c => D3STATE(:,ppZ5c); Z5c=ZERO
O3c => D3STATE(:,ppO3c); O3c=ZERO
O3h => D3STATE(:,ppO3h); O3h=ZERO
O5c => D3STATE(:,ppO5c); O5c=ZERO
R1c => D3STATE(:,ppR1c); R1c=ZERO
R1n => D3STATE(:,ppR1n); R1n=ZERO
R1p => D3STATE(:,ppR1p); R1p=ZERO
R1f => D3STATE(:,ppR1f); R1f=ZERO
R6c => D3STATE(:,ppR6c); R6c=ZERO
R6n => D3STATE(:,ppR6n); R6n=ZERO
R6p => D3STATE(:,ppR6p); R6p=ZERO
R6s => D3STATE(:,ppR6s); R6s=ZERO
R6f => D3STATE(:,ppR6f); R6f=ZERO
R3c => D3STATE(:,ppR3c); R3c=ZERO
R2c => D3STATE(:,ppR2c); R2c=ZERO
N6r => D3STATE(:,ppN6r); N6r=ZERO
O4n => D3STATE(:,ppO4n); O4n=ZERO


#ifdef INCLUDE_SEAICE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocation of SeaIce variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Allocation of Benthic variables
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

Q1c => D2STATE_BEN(:,ppQ1c); Q1c=ZERO
Q1n => D2STATE_BEN(:,ppQ1n); Q1n=ZERO
Q1p => D2STATE_BEN(:,ppQ1p); Q1p=ZERO
Q11c => D2STATE_BEN(:,ppQ11c); Q11c=ZERO
Q11n => D2STATE_BEN(:,ppQ11n); Q11n=ZERO
Q11p => D2STATE_BEN(:,ppQ11p); Q11p=ZERO
Q6c => D2STATE_BEN(:,ppQ6c); Q6c=ZERO
Q6n => D2STATE_BEN(:,ppQ6n); Q6n=ZERO
Q6p => D2STATE_BEN(:,ppQ6p); Q6p=ZERO
Q6s => D2STATE_BEN(:,ppQ6s); Q6s=ZERO
Q16c => D2STATE_BEN(:,ppQ16c); Q16c=ZERO
Q16n => D2STATE_BEN(:,ppQ16n); Q16n=ZERO
Q16p => D2STATE_BEN(:,ppQ16p); Q16p=ZERO
Q16s => D2STATE_BEN(:,ppQ16s); Q16s=ZERO


  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Start the allocation of other pelagic variables which can be outputted
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  LEVEL1 "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-"
  LEVEL1 "# Allocating Other Global Variables .."

ppqpcPPY(iiP1)=42
ppqpcPPY(iiP2)=43
ppqncPPY(iiP1)=44
ppqncPPY(iiP2)=45
ppqscPPY(iiP1)=46
ppqscPPY(iiP2)=47
ppqlcPPY(iiP1)=48
ppqlcPPY(iiP2)=49
ppqfcPPY(iiP1)=50
ppqfcPPY(iiP2)=51
ppqccPPY(iiP1)=52
ppqccPPY(iiP2)=53
ppqpcMEZ(iiZ4)=54
ppqncMEZ(iiZ4)=55
ppqpcMIZ(iiZ5)=56
ppqncMIZ(iiZ5)=57
ppqpcOMT(iiR1)=58
ppqpcOMT(iiR6)=59
ppqpcOMT(iiR3)=60
ppqpcOMT(iiR2)=61
ppqncOMT(iiR1)=62
ppqncOMT(iiR6)=63
ppqncOMT(iiR3)=64
ppqncOMT(iiR2)=65
ppqscOMT(iiR1)=66
ppqscOMT(iiR6)=67
ppqscOMT(iiR3)=68
ppqscOMT(iiR2)=69
ppqfcOMT(iiR1)=70
ppqfcOMT(iiR6)=71
ppqfcOMT(iiR3)=72
ppqfcOMT(iiR2)=73
ppqpcPBA(iiB1)=74
ppqncPBA(iiB1)=75
ppsediPPY(iiP1)=76
ppsediPPY(iiP2)=77
ppsediMIZ(iiZ5)=78
ppsediMEZ(iiZ4)=79
ppsunPPY(iiP1)=80
ppsunPPY(iiP2)=81
ppeiPPY(iiP1)=82
ppeiPPY(iiP2)=83
ppELiPPY(iiP1)=84
ppELiPPY(iiP2)=85

Depth => D3DIAGNOS(:,ppDepth); Depth=ZERO
Volume => D3DIAGNOS(:,ppVolume); Volume=ZERO
Area => D3DIAGNOS(:,ppArea); Area=ZERO
ETW => D3DIAGNOS(:,ppETW); ETW=ZERO
ESW => D3DIAGNOS(:,ppESW); ESW=ZERO
ERHO => D3DIAGNOS(:,ppERHO); ERHO=ZERO
EIR => D3DIAGNOS(:,ppEIR); EIR=ZERO
ESS => D3DIAGNOS(:,ppESS); ESS=ZERO
EPR => D3DIAGNOS(:,ppEPR); EPR=ZERO
DIC => D3DIAGNOS(:,ppDIC); DIC=ZERO
CO2 => D3DIAGNOS(:,ppCO2); CO2=ZERO
pCO2 => D3DIAGNOS(:,pppCO2); pCO2=ZERO
HCO3 => D3DIAGNOS(:,ppHCO3); HCO3=ZERO
CO3 => D3DIAGNOS(:,ppCO3); CO3=ZERO
ALK => D3DIAGNOS(:,ppALK); ALK=ZERO
pH => D3DIAGNOS(:,pppH); pH=ZERO
OCalc => D3DIAGNOS(:,ppOCalc); OCalc=ZERO
OArag => D3DIAGNOS(:,ppOArag); OArag=ZERO
ffCO2 => D3DIAGNOS(:,ppffCO2); ffCO2=ZERO
totpelc => D3DIAGNOS(:,pptotpelc); totpelc=ZERO
totpeln => D3DIAGNOS(:,pptotpeln); totpeln=ZERO
totpelp => D3DIAGNOS(:,pptotpelp); totpelp=ZERO
totpels => D3DIAGNOS(:,pptotpels); totpels=ZERO
totpelf => D3DIAGNOS(:,pptotpelf); totpelf=ZERO
cxoO2 => D3DIAGNOS(:,ppcxoO2); cxoO2=ZERO
eO2mO2 => D3DIAGNOS(:,ppeO2mO2); eO2mO2=ZERO
Chla => D3DIAGNOS(:,ppChla); Chla=ZERO
flPTN6r => D3DIAGNOS(:,ppflPTN6r); flPTN6r=ZERO
flN3O4n => D3DIAGNOS(:,ppflN3O4n); flN3O4n=ZERO
flN4N3n => D3DIAGNOS(:,ppflN4N3n); flN4N3n=ZERO
fscavN7f => D3DIAGNOS(:,ppfscavN7f); fscavN7f=ZERO
IRONDEPO => D3DIAGNOS(:,ppIRONDEPO); IRONDEPO=ZERO
sediR2 => D3DIAGNOS(:,ppsediR2); sediR2=ZERO
sediR3 => D3DIAGNOS(:,ppsediR3); sediR3=ZERO
sediR6 => D3DIAGNOS(:,ppsediR6); sediR6=ZERO
sediO5 => D3DIAGNOS(:,ppsediO5); sediO5=ZERO
xEPS => D3DIAGNOS(:,ppxEPS); xEPS=ZERO
ABIO_eps => D3DIAGNOS(:,ppABIO_eps); ABIO_eps=ZERO
BGE => D3DIAGNOS(:,ppBGE); BGE=ZERO
BIOALK => D3DIAGNOS(:,ppBIOALK); BIOALK=ZERO
DMS => D3DIAGNOS(:,ppDMS); DMS=ZERO

qpcPPY => D3DIAGNOS(:,ppqpcPPY(iiP1): ppqpcPPY(iiP2))
qpcPPY=ZERO
qncPPY => D3DIAGNOS(:,ppqncPPY(iiP1): ppqncPPY(iiP2))
qncPPY=ZERO
qscPPY => D3DIAGNOS(:,ppqscPPY(iiP1): ppqscPPY(iiP2))
qscPPY=ZERO
qlcPPY => D3DIAGNOS(:,ppqlcPPY(iiP1): ppqlcPPY(iiP2))
qlcPPY=ZERO
qfcPPY => D3DIAGNOS(:,ppqfcPPY(iiP1): ppqfcPPY(iiP2))
qfcPPY=ZERO
qccPPY => D3DIAGNOS(:,ppqccPPY(iiP1): ppqccPPY(iiP2))
qccPPY=ZERO
qpcMEZ => D3DIAGNOS(:,ppqpcMEZ(iiZ4): ppqpcMEZ(iiZ4))
qpcMEZ=ZERO
qncMEZ => D3DIAGNOS(:,ppqncMEZ(iiZ4): ppqncMEZ(iiZ4))
qncMEZ=ZERO
qpcMIZ => D3DIAGNOS(:,ppqpcMIZ(iiZ5): ppqpcMIZ(iiZ5))
qpcMIZ=ZERO
qncMIZ => D3DIAGNOS(:,ppqncMIZ(iiZ5): ppqncMIZ(iiZ5))
qncMIZ=ZERO
qpcOMT => D3DIAGNOS(:,ppqpcOMT(iiR1): ppqpcOMT(iiR2))
qpcOMT=ZERO
qncOMT => D3DIAGNOS(:,ppqncOMT(iiR1): ppqncOMT(iiR2))
qncOMT=ZERO
qscOMT => D3DIAGNOS(:,ppqscOMT(iiR1): ppqscOMT(iiR2))
qscOMT=ZERO
qfcOMT => D3DIAGNOS(:,ppqfcOMT(iiR1): ppqfcOMT(iiR2))
qfcOMT=ZERO
qpcPBA => D3DIAGNOS(:,ppqpcPBA(iiB1): ppqpcPBA(iiB1))
qpcPBA=ZERO
qncPBA => D3DIAGNOS(:,ppqncPBA(iiB1): ppqncPBA(iiB1))
qncPBA=ZERO
sediPPY => D3DIAGNOS(:,ppsediPPY(iiP1): ppsediPPY(iiP2))
sediPPY=ZERO
sediMIZ => D3DIAGNOS(:,ppsediMIZ(iiZ5): ppsediMIZ(iiZ5))
sediMIZ=ZERO
sediMEZ => D3DIAGNOS(:,ppsediMEZ(iiZ4): ppsediMEZ(iiZ4))
sediMEZ=ZERO
sunPPY => D3DIAGNOS(:,ppsunPPY(iiP1): ppsunPPY(iiP2))
sunPPY=ZERO
eiPPY => D3DIAGNOS(:,ppeiPPY(iiP1): ppeiPPY(iiP2))
eiPPY=ZERO
ELiPPY => D3DIAGNOS(:,ppELiPPY(iiP1): ppELiPPY(iiP2))
ELiPPY=ZERO


Area2d => D2DIAGNOS(:,ppArea2d); Area2d=ZERO
ThereIsLight => D2DIAGNOS(:,ppThereIsLight); ThereIsLight=ZERO
SUNQ => D2DIAGNOS(:,ppSUNQ); SUNQ=ZERO
EWIND => D2DIAGNOS(:,ppEWIND); EWIND=ZERO
EICE => D2DIAGNOS(:,ppEICE); EICE=ZERO
ETAUB => D2DIAGNOS(:,ppETAUB); ETAUB=ZERO
EPCO2air => D2DIAGNOS(:,ppEPCO2air); EPCO2air=ZERO
CO2airflux => D2DIAGNOS(:,ppCO2airflux); CO2airflux=ZERO
dpco2 => D2DIAGNOS(:,ppdpco2); dpco2=ZERO
dpo2 => D2DIAGNOS(:,ppdpo2); dpo2=ZERO
totsysc => D2DIAGNOS(:,pptotsysc); totsysc=ZERO
totsysn => D2DIAGNOS(:,pptotsysn); totsysn=ZERO
totsysp => D2DIAGNOS(:,pptotsysp); totsysp=ZERO
totsyss => D2DIAGNOS(:,pptotsyss); totsyss=ZERO
totsysf => D2DIAGNOS(:,pptotsysf); totsysf=ZERO


   PELSURFACE => D2DIAGNOS(:,15+1:15+36); PELSURFACE=ZERO
   jsurO2o => D2DIAGNOS(:,15+ppO2o); jsurO2o=ZERO
   jsurN1p => D2DIAGNOS(:,15+ppN1p); jsurN1p=ZERO
   jsurN3n => D2DIAGNOS(:,15+ppN3n); jsurN3n=ZERO
   jsurN4n => D2DIAGNOS(:,15+ppN4n); jsurN4n=ZERO
   jsurN5s => D2DIAGNOS(:,15+ppN5s); jsurN5s=ZERO
   jsurN7f => D2DIAGNOS(:,15+ppN7f); jsurN7f=ZERO
   jsurB1c => D2DIAGNOS(:,15+ppB1c); jsurB1c=ZERO
   jsurP1c => D2DIAGNOS(:,15+ppP1c); jsurP1c=ZERO
   jsurP1n => D2DIAGNOS(:,15+ppP1n); jsurP1n=ZERO
   jsurP1p => D2DIAGNOS(:,15+ppP1p); jsurP1p=ZERO
   jsurP1l => D2DIAGNOS(:,15+ppP1l); jsurP1l=ZERO
   jsurP1s => D2DIAGNOS(:,15+ppP1s); jsurP1s=ZERO
   jsurP1f => D2DIAGNOS(:,15+ppP1f); jsurP1f=ZERO
   jsurP2c => D2DIAGNOS(:,15+ppP2c); jsurP2c=ZERO
   jsurP2n => D2DIAGNOS(:,15+ppP2n); jsurP2n=ZERO
   jsurP2p => D2DIAGNOS(:,15+ppP2p); jsurP2p=ZERO
   jsurP2l => D2DIAGNOS(:,15+ppP2l); jsurP2l=ZERO
   jsurP2f => D2DIAGNOS(:,15+ppP2f); jsurP2f=ZERO
   jsurZ4c => D2DIAGNOS(:,15+ppZ4c); jsurZ4c=ZERO
   jsurZ5c => D2DIAGNOS(:,15+ppZ5c); jsurZ5c=ZERO
   jsurO3c => D2DIAGNOS(:,15+ppO3c); jsurO3c=ZERO
   jsurO3h => D2DIAGNOS(:,15+ppO3h); jsurO3h=ZERO
   jsurO5c => D2DIAGNOS(:,15+ppO5c); jsurO5c=ZERO
   jsurR1c => D2DIAGNOS(:,15+ppR1c); jsurR1c=ZERO
   jsurR1n => D2DIAGNOS(:,15+ppR1n); jsurR1n=ZERO
   jsurR1p => D2DIAGNOS(:,15+ppR1p); jsurR1p=ZERO
   jsurR1f => D2DIAGNOS(:,15+ppR1f); jsurR1f=ZERO
   jsurR6c => D2DIAGNOS(:,15+ppR6c); jsurR6c=ZERO
   jsurR6n => D2DIAGNOS(:,15+ppR6n); jsurR6n=ZERO
   jsurR6p => D2DIAGNOS(:,15+ppR6p); jsurR6p=ZERO
   jsurR6s => D2DIAGNOS(:,15+ppR6s); jsurR6s=ZERO
   jsurR6f => D2DIAGNOS(:,15+ppR6f); jsurR6f=ZERO
   jsurR3c => D2DIAGNOS(:,15+ppR3c); jsurR3c=ZERO
   jsurR2c => D2DIAGNOS(:,15+ppR2c); jsurR2c=ZERO
   jsurN6r => D2DIAGNOS(:,15+ppN6r); jsurN6r=ZERO
   jsurO4n => D2DIAGNOS(:,15+ppO4n); jsurO4n=ZERO

   PELBOTTOM => D2DIAGNOS(:,51+1:51+36); PELBOTTOM=ZERO
   jbotO2o => D2DIAGNOS(:,51+ppO2o); jbotO2o=ZERO
   jbotN1p => D2DIAGNOS(:,51+ppN1p); jbotN1p=ZERO
   jbotN3n => D2DIAGNOS(:,51+ppN3n); jbotN3n=ZERO
   jbotN4n => D2DIAGNOS(:,51+ppN4n); jbotN4n=ZERO
   jbotN5s => D2DIAGNOS(:,51+ppN5s); jbotN5s=ZERO
   jbotN7f => D2DIAGNOS(:,51+ppN7f); jbotN7f=ZERO
   jbotB1c => D2DIAGNOS(:,51+ppB1c); jbotB1c=ZERO
   jbotP1c => D2DIAGNOS(:,51+ppP1c); jbotP1c=ZERO
   jbotP1n => D2DIAGNOS(:,51+ppP1n); jbotP1n=ZERO
   jbotP1p => D2DIAGNOS(:,51+ppP1p); jbotP1p=ZERO
   jbotP1l => D2DIAGNOS(:,51+ppP1l); jbotP1l=ZERO
   jbotP1s => D2DIAGNOS(:,51+ppP1s); jbotP1s=ZERO
   jbotP1f => D2DIAGNOS(:,51+ppP1f); jbotP1f=ZERO
   jbotP2c => D2DIAGNOS(:,51+ppP2c); jbotP2c=ZERO
   jbotP2n => D2DIAGNOS(:,51+ppP2n); jbotP2n=ZERO
   jbotP2p => D2DIAGNOS(:,51+ppP2p); jbotP2p=ZERO
   jbotP2l => D2DIAGNOS(:,51+ppP2l); jbotP2l=ZERO
   jbotP2f => D2DIAGNOS(:,51+ppP2f); jbotP2f=ZERO
   jbotZ4c => D2DIAGNOS(:,51+ppZ4c); jbotZ4c=ZERO
   jbotZ5c => D2DIAGNOS(:,51+ppZ5c); jbotZ5c=ZERO
   jbotO3c => D2DIAGNOS(:,51+ppO3c); jbotO3c=ZERO
   jbotO3h => D2DIAGNOS(:,51+ppO3h); jbotO3h=ZERO
   jbotO5c => D2DIAGNOS(:,51+ppO5c); jbotO5c=ZERO
   jbotR1c => D2DIAGNOS(:,51+ppR1c); jbotR1c=ZERO
   jbotR1n => D2DIAGNOS(:,51+ppR1n); jbotR1n=ZERO
   jbotR1p => D2DIAGNOS(:,51+ppR1p); jbotR1p=ZERO
   jbotR1f => D2DIAGNOS(:,51+ppR1f); jbotR1f=ZERO
   jbotR6c => D2DIAGNOS(:,51+ppR6c); jbotR6c=ZERO
   jbotR6n => D2DIAGNOS(:,51+ppR6n); jbotR6n=ZERO
   jbotR6p => D2DIAGNOS(:,51+ppR6p); jbotR6p=ZERO
   jbotR6s => D2DIAGNOS(:,51+ppR6s); jbotR6s=ZERO
   jbotR6f => D2DIAGNOS(:,51+ppR6f); jbotR6f=ZERO
   jbotR3c => D2DIAGNOS(:,51+ppR3c); jbotR3c=ZERO
   jbotR2c => D2DIAGNOS(:,51+ppR2c); jbotR2c=ZERO
   jbotN6r => D2DIAGNOS(:,51+ppN6r); jbotN6r=ZERO
   jbotO4n => D2DIAGNOS(:,51+ppO4n); jbotO4n=ZERO








#ifdef INCLUDE_SEAICE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Start the allocation of other seaIce vairables which can be outputted 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=







#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Start the allocation of other benthic vairables which can be outputted 
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
ppqpcBOM(iiQ1)=5
ppqpcBOM(iiQ11)=6
ppqpcBOM(iiQ6)=7
ppqpcBOM(iiQ16)=8
ppqncBOM(iiQ1)=9
ppqncBOM(iiQ11)=10
ppqncBOM(iiQ6)=11
ppqncBOM(iiQ16)=12
ppqscBOM(iiQ1)=13
ppqscBOM(iiQ11)=14
ppqscBOM(iiQ6)=15
ppqscBOM(iiQ16)=16

totbenc => D2DIAGNOS_BEN(:,pptotbenc); totbenc=ZERO
totbenn => D2DIAGNOS_BEN(:,pptotbenn); totbenn=ZERO
totbenp => D2DIAGNOS_BEN(:,pptotbenp); totbenp=ZERO
totbens => D2DIAGNOS_BEN(:,pptotbens); totbens=ZERO

qpcBOM => D2DIAGNOS_BEN(:,ppqpcBOM(iiQ1): ppqpcBOM(iiQ16))
qpcBOM=ZERO
qncBOM => D2DIAGNOS_BEN(:,ppqncBOM(iiQ1): ppqncBOM(iiQ16))
qncBOM=ZERO
qscBOM => D2DIAGNOS_BEN(:,ppqscBOM(iiQ1): ppqscBOM(iiQ16))
qscBOM=ZERO






  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Start the allocation of Pelagic vars for calculation of combined fluxes for output
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  LEVEL1 "#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-"
  LEVEL1 "# Allocating Fluxes .."

  allocate( D3FLUX_MATRIX(1:NO_D3_BOX_STATES, 1:NO_D3_BOX_STATES),stat=status )
  allocate( D3FLUX_FUNC(1:NO_BOXES, 1:NO_D3_BOX_FLUX),stat=status )
  D3FLUX_FUNC = 0

  allocate( D3FLUX_MATRIX(ppZ5c, ppO3c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppP1c, ppR6c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppZ4c, ppR3c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppR1c, ppB1c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppP1c, ppR3c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppN4n, ppP2n)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppZ4c, ppR6c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppP2c, ppR6c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppP2c, ppR3c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppR2c, ppB1c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppP2c, ppR1c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppZ5c, ppZ5c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppZ5c, ppZ5c)%dir( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppN3n, ppP2n)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppP1c, ppR1c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppZ4c, ppR1c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppR3c, ppB1c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppB1c, ppR1c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppO3c, ppO5c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppZ5c, ppR2c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppO3c, ppP2c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppP1s, ppR6s)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppN1p, ppP1p)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppN5s, ppP1s)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppB1c, ppR6c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppO2o, ppO2o)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppO2o, ppO2o)%dir( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppZ5c, ppZ4c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppR6s, ppN5s)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppN7f, ppP2f)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppB1c, ppR3c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppZ5c, ppR3c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppN1p, ppP2p)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppZ4c, ppO3c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppB1c, ppZ4c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppO3c, ppP1c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppP1c, ppO3c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppZ5c, ppR6c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppN7f, ppP1f)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppP2c, ppO3c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppB1c, ppR2c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppZ4c, ppZ5c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppP1c, ppZ5c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppZ5c, ppR1c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppP2c, ppZ5c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppP2p, ppR1p)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppP2c, ppR2c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppN3n, ppP1n)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppP1c, ppR2c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppP1p, ppR1p)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppB1c, ppZ5c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppR6c, ppB1c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppZ4c, ppR2c)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppP2c, ppZ4c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppP1p, ppR6p)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppN4n, ppP1n)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppP2p, ppR6p)%p( 1:1 ) )
  allocate( D3FLUX_MATRIX(ppP1c, ppZ4c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppZ4c, ppZ4c)%p( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppZ4c, ppZ4c)%dir( 1:2 ) )
  allocate( D3FLUX_MATRIX(ppB1c, ppO3c)%p( 1:1 ) )


  ! exPPYc = (P.c -> R1c+R2c+R3c+R6c)
  D3FLUX_MATRIX(ppP1c, ppR1c)%p(1)=+12
  D3FLUX_MATRIX(ppP1c, ppR2c)%p(1)=+12
  D3FLUX_MATRIX(ppP1c, ppR3c)%p(1)=+12
  D3FLUX_MATRIX(ppP1c, ppR6c)%p(1)=+12
  D3FLUX_MATRIX(ppP2c, ppR1c)%p(1)=+12
  D3FLUX_MATRIX(ppP2c, ppR2c)%p(1)=+12
  D3FLUX_MATRIX(ppP2c, ppR3c)%p(1)=+12
  D3FLUX_MATRIX(ppP2c, ppR6c)%p(1)=+12
  
  ! exPPYp = (P.p -> R1p+R6p)
  D3FLUX_MATRIX(ppP2p, ppR1p)%p(1)=+24
  D3FLUX_MATRIX(ppP2p, ppR6p)%p(1)=+24
  D3FLUX_MATRIX(ppP1p, ppR1p)%p(1)=+24
  D3FLUX_MATRIX(ppP1p, ppR6p)%p(1)=+24
  
  ! resPPYc = P.c -> O3c
  D3FLUX_MATRIX(ppP2c, ppO3c)%p(1)=+4
  D3FLUX_MATRIX(ppP1c, ppO3c)%p(1)=+4
  
  ! ruPZOOc = Z.c <- P.c
  D3FLUX_MATRIX(ppP2c, ppZ4c)%p(2)=+16
  D3FLUX_MATRIX(ppP1c, ppZ4c)%p(2)=+16
  D3FLUX_MATRIX(ppP2c, ppZ5c)%p(2)=+16
  D3FLUX_MATRIX(ppP1c, ppZ5c)%p(2)=+16
  
  ! ruPPYp = P.p <- N1p
  D3FLUX_MATRIX(ppN1p, ppP1p)%p(1)=+26
  D3FLUX_MATRIX(ppN1p, ppP2p)%p(1)=+26
  
  ! ruPPYc = P.c  <- O3c
  D3FLUX_MATRIX(ppO3c, ppP1c)%p(1)=+1
  D3FLUX_MATRIX(ppO3c, ppP2c)%p(1)=+1
  
  ! ruZZOOc = Z.c <- Z.c
  D3FLUX_MATRIX(ppZ4c, ppZ4c)%p(2)=+18
  D3FLUX_MATRIX(ppZ5c, ppZ4c)%p(2)=+18
  D3FLUX_MATRIX(ppZ4c, ppZ5c)%p(2)=+18
  D3FLUX_MATRIX(ppZ5c, ppZ5c)%p(2)=+18
  D3FLUX_MATRIX(ppZ4c, ppZ4c)%dir(2)=0
  D3FLUX_MATRIX(ppZ5c, ppZ5c)%dir(2)=0
  
  ! resPPYc2 = P2c -> O3c
  D3FLUX_MATRIX(ppP2c, ppO3c)%p(2)=+6
  
  ! resPBAc = B.c -> O3c
  D3FLUX_MATRIX(ppB1c, ppO3c)%p(1)=+7
  
  ! ruPPYn = (P.n <- N3n+N4n)
  D3FLUX_MATRIX(ppN3n, ppP2n)%p(1)=+25
  D3FLUX_MATRIX(ppN4n, ppP2n)%p(1)=+25
  D3FLUX_MATRIX(ppN3n, ppP1n)%p(1)=+25
  D3FLUX_MATRIX(ppN4n, ppP1n)%p(1)=+25
  
  ! prdO5c = O3c -> O5c
  D3FLUX_MATRIX(ppO3c, ppO5c)%p(1)=+11
  
  ! exPBAc = (B.c -> R1c+R2c+R3c+R6c)
  D3FLUX_MATRIX(ppB1c, ppR1c)%p(1)=+22
  D3FLUX_MATRIX(ppB1c, ppR2c)%p(1)=+22
  D3FLUX_MATRIX(ppB1c, ppR3c)%p(1)=+22
  D3FLUX_MATRIX(ppB1c, ppR6c)%p(1)=+22
  
  ! ruR1PBAc = R1c -> B1c
  D3FLUX_MATRIX(ppR1c, ppB1c)%p(2)=+20
  
  ! ruPPYf = P.f <- N7f
  D3FLUX_MATRIX(ppN7f, ppP2f)%p(1)=+27
  D3FLUX_MATRIX(ppN7f, ppP1f)%p(1)=+27
  
  ! ruPBAc = (B.c <- R1c+R2c+R3c+R6c)
  D3FLUX_MATRIX(ppR1c, ppB1c)%p(1)=+19
  D3FLUX_MATRIX(ppR2c, ppB1c)%p(1)=+19
  D3FLUX_MATRIX(ppR3c, ppB1c)%p(1)=+19
  D3FLUX_MATRIX(ppR6c, ppB1c)%p(1)=+19
  
  ! ruPPYc1 = P1c  <- O3c
  D3FLUX_MATRIX(ppO3c, ppP1c)%p(2)=+2
  
  ! ruPPYs = P.s <- N5s
  D3FLUX_MATRIX(ppN5s, ppP1s)%p(1)=+28
  
  ! exPPYc1 = (P1c -> R1c+R2c+R3c+R6c)
  D3FLUX_MATRIX(ppP1c, ppR1c)%p(2)=+13
  D3FLUX_MATRIX(ppP1c, ppR2c)%p(2)=+13
  D3FLUX_MATRIX(ppP1c, ppR3c)%p(2)=+13
  D3FLUX_MATRIX(ppP1c, ppR6c)%p(2)=+13
  
  ! exPPYs = (P.s -> R6s)
  D3FLUX_MATRIX(ppP1s, ppR6s)%p(1)=+29
  
  ! resPELo = (O2o->*)
  D3FLUX_MATRIX(ppO2o, ppO2o)%p(1)=+31
  D3FLUX_MATRIX(ppO2o, ppO2o)%dir(1)=1
  
  ! remR6s = (R6s -> N5s)
  D3FLUX_MATRIX(ppR6s, ppN5s)%p(1)=+30
  
  ! resZOOc4 = Z4c -> O3c
  D3FLUX_MATRIX(ppZ4c, ppO3c)%p(2)=+9
  
  ! ruPPYc2 = P2c  <- O3c
  D3FLUX_MATRIX(ppO3c, ppP2c)%p(2)=+3
  
  ! exPPYc2 = (P2c -> R1c+R2c+R3c+R6c)
  D3FLUX_MATRIX(ppP2c, ppR1c)%p(2)=+14
  D3FLUX_MATRIX(ppP2c, ppR2c)%p(2)=+14
  D3FLUX_MATRIX(ppP2c, ppR3c)%p(2)=+14
  D3FLUX_MATRIX(ppP2c, ppR6c)%p(2)=+14
  
  ! resZOOc = Z.c -> O3c
  D3FLUX_MATRIX(ppZ5c, ppO3c)%p(1)=+8
  D3FLUX_MATRIX(ppZ4c, ppO3c)%p(1)=+8
  
  ! exZOOc = (Z.c -> R1c+R2c+R3c+R6c)
  D3FLUX_MATRIX(ppZ5c, ppR1c)%p(1)=+23
  D3FLUX_MATRIX(ppZ5c, ppR2c)%p(1)=+23
  D3FLUX_MATRIX(ppZ5c, ppR3c)%p(1)=+23
  D3FLUX_MATRIX(ppZ5c, ppR6c)%p(1)=+23
  D3FLUX_MATRIX(ppZ4c, ppR1c)%p(1)=+23
  D3FLUX_MATRIX(ppZ4c, ppR2c)%p(1)=+23
  D3FLUX_MATRIX(ppZ4c, ppR3c)%p(1)=+23
  D3FLUX_MATRIX(ppZ4c, ppR6c)%p(1)=+23
  
  ! resPPYc1 = P1c -> O3c
  D3FLUX_MATRIX(ppP1c, ppO3c)%p(2)=+5
  
  ! ruZOOc = (Z.c<-P.c+B.c+Z.c)
  D3FLUX_MATRIX(ppP2c, ppZ5c)%p(1)=+15
  D3FLUX_MATRIX(ppP1c, ppZ5c)%p(1)=+15
  D3FLUX_MATRIX(ppB1c, ppZ5c)%p(1)=+15
  D3FLUX_MATRIX(ppZ5c, ppZ5c)%p(1)=+15
  D3FLUX_MATRIX(ppZ4c, ppZ5c)%p(1)=+15
  D3FLUX_MATRIX(ppP2c, ppZ4c)%p(1)=+15
  D3FLUX_MATRIX(ppP1c, ppZ4c)%p(1)=+15
  D3FLUX_MATRIX(ppB1c, ppZ4c)%p(1)=+15
  D3FLUX_MATRIX(ppZ5c, ppZ4c)%p(1)=+15
  D3FLUX_MATRIX(ppZ4c, ppZ4c)%p(1)=+15
  D3FLUX_MATRIX(ppZ5c, ppZ5c)%dir(1)=0
  D3FLUX_MATRIX(ppZ4c, ppZ4c)%dir(1)=0
  
  ! resZOOc5 = Z5c -> O3c
  D3FLUX_MATRIX(ppZ5c, ppO3c)%p(2)=+10
  
  ! ruR6PBAc = R6c -> B1c
  D3FLUX_MATRIX(ppR6c, ppB1c)%p(2)=+21
  
  ! ruBZOOc = Z.c <- B.c
  D3FLUX_MATRIX(ppB1c, ppZ4c)%p(2)=+17
  D3FLUX_MATRIX(ppB1c, ppZ5c)%p(2)=+17
  


#ifdef INCLUDE_SEAICE
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Start the allocation of Seaice vars for calculation of combined fluxes for output
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  allocate( D2FLUX_MATRIX_ICE(1:NO_D2_BOX_STATES_ICE, 1:NO_D2_BOX_STATES_ICE),stat=status )


#endif

  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! Start the allocation of Benthic vars for calculation of combined fluxes for output
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  allocate( D2FLUX_MATRIX_BEN(1:NO_D2_BOX_STATES_BEN, 1:NO_D2_BOX_STATES_BEN),stat=status )




  end subroutine AllocateMem

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
! MODEL  BFM - Biogeochemical Flux Model
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
