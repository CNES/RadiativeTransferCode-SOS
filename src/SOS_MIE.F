C******************************************************************************
C* Copyright 2017, Centre National d'Etudes Spatiales (CNES)
C* 
C* This file is part of SOS.
C* 
C* SOS is free software: you can redistribute it and/or modify
C* it under the terms of the GNU General Public License as published by
C* the Free Software Foundation, either version 3 of the License, or
C* (at your option) any later version.
C* 
C* SOS is distributed in the hope that it will be useful,
C* but WITHOUT ANY WARRANTY; without even the implied warranty of
C* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C* GNU General Public License for more details.
C* 
C* You should have received a copy of the GNU General Public License
C* along with SOS. If not, see <http://www.gnu.org/licenses/>.
C******************************************************************************


C******************************************************************************
C* FICHIER: SOS_MIE.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Ce module gere les procedures permettant les calculs de la theorie
C*       de MIE.
C*
C* AUTEUR: 
C*   Code initial: Laboratoire d'Optique Atmospherique (LOA, Villeneuve d'Ascq).
C*   Reprise du codage et developpement : CS SI, B.Lafrance, C. Casteras.
C* DATE: 30/04/02
C*
C* MOD:VERSION:1.0 : Mise a niveau du code des OS du LOA
C*             1.1 : Ajustement des commentaires
C*
C* MOD:VERSION:2.0: 10/01/2010 
C*   - Gestion des angles : 
C*         * Suppression de la lecture d'un fichier d'angles de Gauss
C*         * Passage des arguments d'entree : 
C*                 MIE_NBMU : nombre d'angles utiles
C*                 RMU      : tableau des cosinus des angles
C*                 CHR      : tableau des poids des angles 
C*
C*   - Suppression de la constante SOS_LENDIR  
C*
C*   - Modification de la taille des tableaux des fonctions de phase :
C*     ==> Passage de SOS_MIE_NBMU a SOS_MIE_NBMU_MAX.
C*     ==> Boucle sur les angles jusqu'a MIE_NBMU egal au nombre effectif d'angles utilises.
C*
C*   - SOS_FPHASE_MIE : 
C*         * Ajout de l'argument d'entree : MIE_NBMU (nb angles utiles)
C*         * Modification de la ligne d'enregistrement des tableaux  
C*           de fonctions de phase par valeur de parametre de taille Alpha
C*           ==> limitation aux bornes utiles -MIE_NBMU:MIE_NBMU
C* 
C* MOD:VERSION:2.1: 24/06/2010 
C*     Les etiquettes de FORMAT en ,X, passent en ,1X, pour eviter des erreurs de compilation.
C*
C* MOD:VERSION:2.2: 08/12/2015 
C*    - Ajustements mineurs pour le respect strict des 72 colonnes 
C*      (incluant au décodage des constante) : requis pour compilation gfortran                   
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    SOS_PI : Valeur de PI
C    SOS_MIE_NBMU_MAX : Nombre maximal d'angles de Gauss positifs
C    SOS_MIE_DIM : Dimension des tableaux
C                              
C----------------------------------------------------------------------------
#include "SOS.h"

C==============================================================================
C PROCEDURE: SOS_MIE
C ==========      
C      Cette procedure applique la theorie de MIE pour des particules spheriques.
C
C      La nature chimique des particules est caracterisee par leur indice de 
C      refraction relativement a l'environnement : m.
C      m = mr + i * mi, avec "mr" la partie reelle et "mi" la partie imaginaire.
C      La partie imaginaire est entree negative. Elle est responsable des effets
C      d'absorption.
C
C      Les calculs de Mie sont realises en fonction du parametre de taille : 
C      alpha = 2*pi * r / wa, ou r est le rayon de la particule et wa est la 
C      longueur d'onde du rayonnement. 
C
C      Le pas sur alpha est:
C         si        alpha  <= 0.1  ===> pas = 0.0001
C         si 0.1  < alpha  <= 1.0  ===> pas = 0.001
C         si 1.0  < alpha  <= 10   ===> pas = 0.01
C         si 10   < alpha  <= 30   ===> pas = 0.05
C         si 30   < alpha  <= 100  ===> pas = 0.1
C         si 100  < alpha          ===> pas = 1.0
C
C
C      Ce programme calcule pour l'indice m :      
C      
C         - Le coefficient d'efficacite en extinction: Qext(alpha) 
C           de alpha0 a alphaF
C         - Le coefficient d'efficacite en diffusion:  Qsca(alpha) 
C           de alpha0 a alphaF
C         - Le facteur d'anisotropie: G
C         - Les fonctions de phases : Imie(alpha,mu), 
C                                     Qmie(alpha,mu),
C				      Umie(alpha,mu)
C 	    de alpha0 a alphaF et pour mu, cosinus de l'angle de diffusion, des
C           angles de Gauss entre 0 et pi.
C
C Description des parametres
C ----------------------------
C       MIE_NBMU (I4) : (E) Nombre d'angles (positifs) effectivement utiles
C	RMU(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX) (double) : (E) Cosinus des angles.
C	CHR(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX) (double) : (E) Poids des angles.
C	RN (double) (E) : partie reelle de l'indice de refraction.
C	IN (double) (E) : partie imaginaire de l'indice de refraction.
C	ALPHAO (double) (E) : valeur de alpha a partir de laquelle 
C			      on commence le calcul.
C	ALPHAF (double) (E) : valeur de alpha pour laquelle on arrete le calcul.			  
C	FICMIE (CHARACTER*SOS_LENFIC2) (E) : chemin et nom du fichier MIE.
C	FICLOG (CHARACTER*SOS_LENFIC2) (E) : nom du fichier TRACE de MIE .
C	IER (I4) (S) : code d'erreur = 0 si pas d'erreur, sinon =1	
C
C
C Resultat fourni par la procedure
C --------------------------------
C   Fichier de MIE
C   Fichier de trace (optionnel)
C
C
C Variable d'environnement
C ------------------------
C   Aucune
C 
C Description des fichiers utilises
C ---------------------------------
C   -->Contenu du fichier de MIE (cree)
C       Fichier binaire non formatte, cree par SOS_MIE
C       1 enregistrement d'entete avec
C                Indice de refraction (partie reelle)
C                Indice de refraction (partie imaginaire)
C                Parametre de taille limite
C                Nombre d'angles de Gauss positif
C
C      1 enregistrement pour chaque valeur du parametre de taille avec:
C        - Parametre de taille (reel),
C        - Coefficient d'efficacite en extinction (reel),
C        - Coefficient d'efficacite en diffusion (reel),
C        - Facteur d'anisotropie (double precision),
C        - Tableau de la fonction de phase Imie(mu) (reel),
C        - Tableau de la fonction de phase polarisee Qmie(mu) (reel),
C        - Tableau de la fonction de phase polarisee Umie(mu) (reel).
C
C
C Common utilise:
C --------------
C  Aucun
C
C
C Cas d'erreur :
C ------------
C     -Decodage des parametres incorrect
C     -Parametres du programme incorrects ALPHA0 >ALPHAF
C     -Ouverture d'un fichier
C      Erreur si le fichier de MIE existe deja
C     -Lecture/ecriture dans un fichier
C
C     Affichage d'un message a l'ecran, arret du programme et 
C     retour du status 1 au shell
C        
C   
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel aux routines:
C      - SOS_XERO
C      - SOS_FPHASE_MIE
C
C==============================================================================

      SUBROUTINE SOS_MIE(MIE_NBMU,RMU,CHR,RN,IN,ALPHAO,ALPHAF,
     &                   FICMIE,FICLOG,IER)


      IMPLICIT NONE


C* Definition des variables                                                   
C----------------------------------------------------------------------------        
      DOUBLE PRECISION RN	! partie reelle de l'indice de refraction.
      DOUBLE PRECISION IN	! partie imaginaire de l'indice de refraction.
      DOUBLE PRECISION ALPHAO   ! valeur de alpha a partir de laquelle 
      				! on commence le calcul.				
      DOUBLE PRECISION ALPHAF   ! valeur de alpha pour laquelle on arrete le calcul.
                                          
      DOUBLE PRECISION RMU(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)   !Cosinus des angles.
      DOUBLE PRECISION CHR(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)   !Poids des angles.

      DOUBLE PRECISION CNA(-1:SOS_MIE_DIM)   !Fonction reelle Cn(alpha). 
      DOUBLE PRECISION SNA(-1:SOS_MIE_DIM)   !Fonction reelle Sn(alpha).
           
      DOUBLE PRECISION RGNA(-1:SOS_MIE_DIM)  !Partie reelle de la fonction Gn(alpha).	
      DOUBLE PRECISION IGNA(-1:SOS_MIE_DIM)  !Partie imaginaire de Gn(alpha).
      
      DOUBLE PRECISION RDNA(0:SOS_MIE_DIM)   !Fonction reelle Dn(alpha).
      DOUBLE PRECISION RDNB(0:SOS_MIE_DIM)   !Partie reelle de la fonction Dn(alpha*m).
      DOUBLE PRECISION IDNB(0:SOS_MIE_DIM)   !Partie imaginaire de  Dn(alpha*m).
      
      DOUBLE PRECISION RA(0:SOS_MIE_DIM)     !Partie reelle de la fonction An(alpha). 
      DOUBLE PRECISION IA(0:SOS_MIE_DIM)     !Partie imaginaire de An(alpha). 
      DOUBLE PRECISION RB(0:SOS_MIE_DIM)     !Partie reelle de la fonction Bn(alpha). 
      DOUBLE PRECISION IB(0:SOS_MIE_DIM)     !Partie imaginaire de Bn(alpha).
      
      

      
      DOUBLE PRECISION QEXT		!Coefficient d'efficacite en extinction.
      DOUBLE PRECISION QSCA		!Coefficient d'efficacite en diffusion.
      DOUBLE PRECISION G		!Facteur d'anisotropie.
      
      DOUBLE PRECISION ALPHA		!Valeur courante du parametre de taille.
      DOUBLE PRECISION PAS		!Pas entre deux valeurs de ALPHA.
      
      DOUBLE PRECISION RBETA		!Partie reelle de m*alpha.
      DOUBLE PRECISION IBETA            !Partie imaginaire de m*alpha.
      

      DOUBLE PRECISION X
      DOUBLE PRECISION Y
      DOUBLE PRECISION Z
      DOUBLE PRECISION W
      
      DOUBLE PRECISION X1
      DOUBLE PRECISION X2
      DOUBLE PRECISION X3
      DOUBLE PRECISION X4
      DOUBLE PRECISION X5
      DOUBLE PRECISION X6
      DOUBLE PRECISION X7
           
      DOUBLE PRECISION Y1
      DOUBLE PRECISION Y2
      DOUBLE PRECISION Y3
      DOUBLE PRECISION Y4
      DOUBLE PRECISION Y5
      DOUBLE PRECISION Y6
      DOUBLE PRECISION Y7
      DOUBLE PRECISION Y8
      
      DOUBLE PRECISION Z3
      DOUBLE PRECISION Z4
      DOUBLE PRECISION Z5
      DOUBLE PRECISION Z6
      DOUBLE PRECISION Z7
      DOUBLE PRECISION Z8
      DOUBLE PRECISION Z9
      
      DOUBLE PRECISION T
      DOUBLE PRECISION Q
      
      DOUBLE PRECISION W6		!Valeur de 2 / ALPHA**2
      DOUBLE PRECISION XX
      DOUBLE PRECISION YY
      DOUBLE PRECISION ZZ
      DOUBLE PRECISION TT
      
      DOUBLE PRECISION A2

      INTEGER*4 MIE_NBMU   		! Nombre d'angles effectifs des simulations      
      INTEGER*4 UN			!Valeur 1 ou -1
      INTEGER*4 TEST
      
      INTEGER*4 N1			! Ordre maximal de calcul des suites Dn et Sn
      INTEGER*4 N2			! Ordre maximal de calcul des suites Gn et Cn
      INTEGER*4 N2P1			
      INTEGER*4 M
       
      INTEGER*4 J				
      INTEGER*4 I
      INTEGER*4 N
      INTEGER*4 DIM                   ! dimension
      INTEGER*4 IFIN                  ! nb caracteres significatifs d'une chaine
      
      LOGICAL TRACE                   ! = vrai,si ecriture dans le fichier trace
       
      CHARACTER*SOS_LENFIC2 FICMIE    ! chemin d'acces au fichier MIE, nom du fichier MIE
      CHARACTER*SOS_LENFIC2 FICLOG    ! nom du fichier TRACE de MIE
      
      INTEGER*4 IER       	      ! code d'erreur =0 si pas d'erreur, =-1 sinon      
      
                       
C* Initialisation du code d'erreur
C------------------------------------------------------
      IER=0     

C* Test de compatibilite entre alphaF et l'ordre maximal 
C* des developpements pour le calcul des fonctions de 
C* Ricatti-Bessel (et de leurs derivees)
C------------------------------------------------------
      IF (DINT(ALPHAF+ALPHAF+20).GT.SOS_MIE_DIM) GOTO 997



C* Lecture des parametres/arguments d'entree et decodage 
C------------------------------------------------------      
      IF ((FICLOG(1:1).EQ.'0').AND.(FICLOG(2:2).EQ.' ')) THEN
	TRACE=.FALSE.
      ELSE
        TRACE=.TRUE.
      END IF     
      
C* Ouverture des fichiers 
C------------------------------------------------------
      OPEN (1,FILE=FICMIE,FORM='UNFORMATTED',
     *      STATUS='NEW',ERR=911)
            
            
C     Fichier trace
      IF (TRACE) THEN
         OPEN (99,FILE=FICLOG,ERR=913)
	 WRITE(99,*,err=923) ' '
         WRITE(99,*,err=923) 'Calcul du fichier MIE'     
      END IF
 
               
C* Initialisation des tableaux a 0   
C------------------------------------------------------  
      DIM=SOS_MIE_DIM 
      CALL SOS_XERO(CNA,-1,DIM)
      CALL SOS_XERO(SNA,-1,DIM)
      CALL SOS_XERO(RGNA,-1,DIM)  
      CALL SOS_XERO(IGNA,-1,DIM)
      CALL SOS_XERO(RDNA,0,DIM) 
      CALL SOS_XERO(RDNB,0,DIM)           
      CALL SOS_XERO(IDNB,0,DIM)    
      
                                                       
C* Fichier trace 
C------------------------------------------------------  
      IF (TRACE) THEN
        WRITE(99,*,err=923)'PARAMETRES DE SIMULATION '
        WRITE(99,*,err=923)
     *      ' Indice de refraction (partie reelle): ',RN
        WRITE(99,*,err=923)
     *      ' Indice de refraction (partie imaginaire): ',IN
        WRITE(99,*,err=923)' Parametres de tailles limites: ',
     *           ALPHAO,' et ',ALPHAF
        WRITE(99,*,err=923)' '
        WRITE(99,*,err=923)' Liste des angles utilises : '
	WRITE(99,*,err=923) "   Indice, mu, angle (degres), poids"
	DO J=-MIE_NBMU,MIE_NBMU
	    WRITE(99,222,ERR=923) J, RMU(J),
     &	          DACOS(RMU(J))*180./SOS_PI, CHR(J)
	ENDDO
	WRITE(99,*,err=923)' '
	WRITE(99,*,err=923)'RESULTATS'
	WRITE(99,*,err=923)' Valeur du parametre de taille alpha'
	WRITE(99,*,err=923)
     *      ' Valeur du coefficient d''efficacite en extinction Qext'
	WRITE(99,*,err=923)
     *      ' Valeur du coefficient d''efficacite en diffusion Qsca'
	WRITE(99,*,err=923)' Valeur du facteur d''anisotropie G'
	WRITE(99,*,err=923)' '
        WRITE(99,10,err=923)'  alpha','Qext','Qsca','G'
      END IF

C* Ecriture de l'entete du fichier de sortie  
C------------------------------------------------------                         
      WRITE(1,err=921)RN,IN,ALPHAF,MIE_NBMU
      
  
  
C* Boucle sur les valeurs du parametre de taille ALPHA de
C* ALPHA0 a ALPHAF par pas variable. (etiquette 7777)
C------------------------------------------------------          
      ALPHA=ALPHAO
           
 7777 CONTINUE   
 
C* Estimation du pas pour la valeur de ALPHA suivante  
C------------------------------------------------------
         PAS=0.0001
         IF(ALPHA.gt.0.1)  PAS =0.001
         IF(ALPHA.gt.1.00) PAS =0.01
         IF(ALPHA.gt.10.)  PAS =0.05
         IF(ALPHA.gt.30.)  PAS =0.10
         IF(ALPHA.gt.100.) PAS =1.00
 
 
      
C* Calcul des fonctions derivees de Ricatti-Bessel : 
C* Gn(alpha), Dn(alpha), Dn(alpha*m) - parties reelle et imaginaire
C*
C* Calcul des fonctions reelles Sn(alpha) et Cn(alpha)
C------------------------------------------------------    

C*  Ordres de developpement limites : N1 pour Dn et Sn
C*                                    N2 pour Gn et Cn
                 
         N1   = DINT(ALPHA+ALPHA+20)
	 IF(N1.GT.SOS_MIE_DIM) GOTO 997  
         N2   = DINT(ALPHA+ALPHA+5)  
         N2P1 = N2+1  
	 
	 

C*  Recurrence sur Gn(alpha) : partie reelle RGNA et imaginaire IGNA, 
C*  et sur Cn(alpha) de n=1 a N2.
C------------------------------------------------------  
      
C*   Initialisation des suites Cn (CNA reel)
C*   et Gn (partie reelle RGNA et imaginaire IGNA)    

         CNA(-1) = -DSIN(ALPHA)   
	 CNA(0)  =  DCOS(ALPHA)
	     
         RGNA(-1) =  0.D+00
         RGNA(0)  =  0.D+00         
         IGNA(-1) =  0.D+00
         IGNA(0)  = -1.D+00


C*  Boucle de recurrence
    
         DO 25 I=1,N2
            X=RGNA(I-1)
            Z=I/ALPHA
            Y=IGNA(I-1)
            W=((Z-X)*(Z-X)+(Y*Y)) 
            RGNA(I)=(Z-X)/W-Z 
            IGNA(I)=Y/W

            CNA(I)=(2*I-1.D+00)*CNA(I-1)/ALPHA-CNA(I-2)	
	    
C*    Controle : 
C*      Si Cn(alpha) ne diverge pas => on poursuit la recurrence.   
C*      Sinon, on redefinit les ordres limites N2 et N1, puis on quitte
C*      la boucle de recurrence.	    
            IF(CNA(I).LT.1.D+304) GO TO 25	    
            N2=I 
            N2P1=I+1
            N1=I+15
            GOTO 100
	    
   25    CONTINUE
  100    CONTINUE



C*  Recurrence sur Dn(alpha)   : fonction reelle RDNA
C*                 Dn(alpha*m) : partie reelle RDNB et imaginaire IDNB
C*                 Sn(alpha)   : fonction reelle SNA
C*  Boucle de N1 a 0
C------------------------------------------------------         

         RBETA=RN*ALPHA       
         IBETA=IN*ALPHA    
         X1=RBETA*RBETA+IBETA*IBETA 
         X2=RBETA/X1
         X3=IBETA/X1
	 
C*   Initialisation a N1            
         RDNA(N1)=0.D+00
         RDNB(N1)=0.D+00
         IDNB(N1)=0.D+00
      
         SNA(N1)=0.D+00 
         SNA(N1-1)=1.D+00

C*   Boucle de recurrence de N1-1 a 0    
         DO 30 I=N1-1,0,-1
	    
            X=RDNB(I+1) 
            Y=IDNB(I+1)
            Z=X+(I+1.D+00)*X2
            W=Y-(I+1.D+00)*X3 
            X4=Z*Z+W*W 
            RDNB(I)=(I+1.D+00)*X2-Z/X4  
            IDNB(I)=-(I+1.D+00)*X3+W/X4
	    
            Z=(I+1.D+00)/ALPHA
            X=RDNA(I+1)  
            RDNA(I)=Z-1.D+00/(X+Z)
	    
            SNA(I-1)=(2.D+00*I+1.D+00)*SNA(I)/ALPHA-SNA(I+1)

C*   Si les valeurs de Sn deviennent trop grandes, on normalise
C*   l'ensemble des termes calcules par Sn-1. 	    
            IF (SNA(I-1).GT.1D+304) THEN    	     	
                TEST=I-1 
                X=SNA(TEST)
                DO 35 J=TEST,N2
                   SNA(J)=SNA(J)/X 
   35           CONTINUE
            ENDIF
	    		  
   30    CONTINUE 
 
 
C*   Normalisation de Sn(alpha) par sin(alpha)/S0(alpha)
C*   pour verifier la condition n=0 : S0(alpha) = sin(alpha)                 
         Q=SNA(0)/DSIN(ALPHA)    
                                                           
         DO 45 I=0,N2       
             SNA(I)=SNA(I)/Q 
   45    CONTINUE
   
         
          


C* Calcul des fonctions An(alpha) et Bn(alpha) :
C* parties reelle (RA, RB) et imaginaire (IA, IB)
C------------------------------------------------------    
         UN=1 
	                                        
         DO 50 I=1,N2     
            X1=SNA(I)
            X2=CNA(I)
            X3=RDNB(I)
            X4=IDNB(I)                
            X5=RDNA(I)
            X6=RGNA(I)         
            X7=IGNA(I)           
            Y1=X3-RN*X5          
            Y2=X4-IN*X5         
            Y3=X3-RN*X6+IN*X7         
            Y4=X4-RN*X7-IN*X6       
            Y5=RN*X3-IN*X4-X5        
            Y6=IN*X3+RN*X4        
            Y7=RN*X3-IN*X4-X6       
            Y8=IN*X3+RN*X4-X7       
            Z4=Y2*Y3-Y1*Y4       
            Z3=Y1*Y3+Y2*Y4        
            Z5=X1*X1+X2*X2       
            Z6=Y3*Y3+Y4*Y4
            Z7=Y5*Y7+Y6*Y8       
            Z8=Y6*Y7-Y5*Y8         
            Z9=Y7*Y7+Y8*Y8          
            Q=(I+I+1.D+00)/I/(I+1.D+00)*UN

C*   Test sur X2 (fonction Cn(alpha)) : tend a diverger pour n important.
C*   Z5 tend vers l'infini en X2**2. La valeur de X2 / Z5 tend alors vers 0.
C*   Ce test evite un probleme numerique de rapport Inf / Inf.          	    
	    IF (X2.GT.1D+300) THEN
		Y1=0.D+00
		Y2=0.D+00
		Y3=0.D+00
		Y4=0.D+00		
	    ELSE
                Y1=X1*(X1*Z3+X2*Z4)/Z5/Z6          
                Y2=X1*(X1*Z4-X2*Z3)/Z5/Z6         
                Y3=X1*(X1*Z7+X2*Z8)/Z5/Z9          
                Y4=X1*(X1*Z8-X2*Z7)/Z5/Z9  
	    ENDIF
            
            RA(I)=Y2*Q               
            IB(I)=Y3*Q           
            Q=-Q              
            RB(I)=Y4*Q         
            IA(I)=Y1*Q        
            UN=-UN    
   50    CONTINUE
    	 
                                                                 
         RA(0)=0.D+00
         IA(0)=0.D+00
         RB(0)=0.D+00
         IB(0)=0.D+00
         RA(N2P1)=0.D+00 
         IA(N2P1)=0.D+00
         RB(N2P1)=0.D+00
         IB(N2P1)=0.D+00


C* Calcul des coefficients d'efficacite en diffusion (Qsca) 
C* et en extinction (Qext).
C* Calcul du facteur d'anisotropie (G).
C------------------------------------------------------
         QEXT=0.D+00
         QSCA=0.D+00
         G=0.D+00
         J=-1
         X=RA(1)      
         Y=IA(1)       
         Z=RB(1)        
         T=IB(1)
                                                                         
         DO 55 N=1,N2         
            M=N+1
            XX=RA(M) 
            YY=IA(M)
            ZZ=RB(M) 
            TT=IB(M) 
            A2=(N+1.D+00)
            QEXT=QEXT+N*A2*J*(Y-T) 
            QSCA=QSCA+N*N*A2*A2/(N+A2)*(X*X+Y*Y+Z*Z+T*T)    
            J=-J  
            G=G-A2*N/(A2+N)*(N*(A2+1.)*(A2+1.)/(2.*N+3.)*
     1        (Y*YY+X*XX+T*TT+Z*ZZ)+Y*T+X*Z)
            X=XX
            Y=YY
            Z=ZZ
            T=TT
   55    CONTINUE 
                                                             
         W6=2.D+00/ALPHA/ALPHA  
         QEXT=W6*QEXT         
         QSCA=W6*QSCA             
         G=4.D+00*G/QSCA/ALPHA/ALPHA

C* Appel de la routine permettant le calcul des fonctions de phase de Mie
C* et enregistrant les resultats pour ALPHA.       
C------------------------------------------------------
         CALL SOS_FPHASE_MIE(MIE_NBMU,RMU,RA,IA,RB,IB,ALPHA,
     *                      QEXT,QSCA,G,N2,TRACE,IER) 
         IF (IER.NE.0) GOTO 996   ! cas d'erreur
      
C* Passage a la valeur suivante du parametre de taille ALPHA (arret si la 
C* limite ALPHAF est depassee).                 
C------------------------------------------------------
         ALPHA=ALPHA+PAS 
      IF (ALPHA.LE.ALPHAF) GOTO 7777 
                                                   
 
C* Fermeture fichier
C-------------------
      CLOSE(1)
      CLOSE(4)

      
C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur et retour du status 1 au shell 
C----------------------------------------------
  911 WRITE(6,*) 'Erreur a l''ouverture du fichier MIE'
      IER=-1  
      GOTO 9999
  913 WRITE(6,*) 'Erreur a l''ouverture du fichier trace MIE'
      IER=-1  
      GOTO 9999        
  921 WRITE(6,*) 'Erreur d''ecriture sur le fichier MIE'
      IER=-1  
      GOTO 9999
  922 WRITE(6,*) 'Erreur de lecture sur le fichier angles de Gauss'
      IER=-1  
      GOTO 9999
  923 WRITE(6,*) 'Erreur d''ecriture sur le fichier trace MIE'
      IER=-1  
      GOTO 9999    
  996 WRITE(6,*) 'Erreur dans la routine SOS_PHASE_MIE'
      IER=-1
      GOTO 9999
  997 WRITE(6,*) 'Valeur AlphaMax trop grande devant SOS_MIE_DIM'
      IER=-1
      GOTO 9999      
  
C* Format
C---------   
   10 FORMAT(A7,17X,A4,16X,A4,15X,A1)          
  222 FORMAT(I4,1X,D21.14,2X,F7.3,2X,D21.14,2X)   
                                                   
 9999 CONTINUE
 
      RETURN
      END    !FIN DE LA PROCEDURE SOS_MIE   
      

 
                                                                      

C==============================================================================
C PROCEDURE: SOS_XERO
C ==========      
C  Cette procedure affecte zero aux elements d'un tableau IT(N:M)
C
C
C Description des parametres
C --------------------------
C    N,M (I4) (E)  Dimensions du tableau IT           
C    IT(N:M) (double) (E/S)  Tableau a initialiser
C 
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     Aucun
C
C==============================================================================
      SUBROUTINE SOS_XERO(IT,N,M) 
   
      IMPLICIT NONE
      
C* Definition des variables 
C---------------------------
      INTEGER*4 N,M			! Dimensions du tableau   
      DOUBLE PRECISION IT(N:M)		! Tableau 
               
      INTEGER*4 I			!Indice du tableau
      
C*   Initialisation du tableau 
C-----------------------------           
      DO 1 I=N,M  
    1 IT(I)=0.D+00
    
      RETURN 
      END  	!FIN DE LA PROCEDURE SOS_XERO       
      
 
                         
                                                        
C==============================================================================
C PROCEDURE: SOS_FPHASE_MIE
C ==========      
C      Cette procedure calcule les fonctions de phase Imie(mu), Qmie(mu),
C      Umie(mu) pour une valeur fixee du parametre de taille ALPHA et 
C      pour les mu de Gauss du tableau RMU.
C
C      Il y a enregistrement (non formatte) des resultats:
C        - Parametre de taille (reel),
C        - Coefficient d'efficacite en extinction (reel),
C        - Coefficient d'efficacite en diffusion (reel),
C        - Facteur d'anisotropie (double precision),
C        - Tableau de la fonction de phase Imie(mu) (reel),
C        - Tableau de la fonction de phase polarisee Qmie(mu) (reel),
C        - Tableau de la fonction de phase polarisee Umie(mu) (reel).
C
C    --> Les tableaux ont une taille utile limitee a MIE_NBMU : 
C        ils sont ecris avec cette limitation.
C
C     
C Description des parametres
C --------------------------
C    MIE_NBMU (I4)   (E)   Nombre d'angles (positifs) effectivement utiles
C    RMU(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX) (double) (E) Cosinus des angles de Gauss  
C    RA(0:SOS_MIE_DIM)  (double) (E)   Partie reelle de la fonction An(alpha).
C    IA(0:SOS_MIE_DIM)  (double) (E)   Partie imaginaire de An(alpha).  
C    RB(0:SOS_MIE_DIM)  (double) (E)   Partie reelle de la fonction Bn(alpha). 
C    IB(0:SOS_MIE_DIM)  (double) (E)   Partie imaginaire de Bn(alpha).
C    ALPHA (double)  (E)   Valeur du parametre de taille. 
C    KMA1  (double)  (E)   Coefficient d'efficacite en extinction.
C    KMA2  (double)  (E)   Coefficient d'efficacite en diffusion.
C    G     (double)  (E)   Facteur d'anisotropie.
C    N2    (I4)      (E)   Ordre maximal de sommation sur N.
C    TRACE (logical) (E)   vrai,si ecriture dans le fichier trace
C    IER (I4) (S) : code d'erreur = 0 si pas d'erreur, sinon =1	
C
C
C Resultat stocke dans le fichier MIE:
C ------------------------------------
C    AY1   (R4)	        Valeur du parametre de taille ALPHA.
C    AY2   (R4)	        Coefficient d'efficacite en extinction.
C    AY3   (R4)	        Coefficient d'efficacite en diffusion.
C    AY4   (double)     Facteur d'anisotropie.
C    IMIE(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)    (R4)	      
C                 fonction de phase : Imie(mu) pour ALPHA.
C    QMIE(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)    (R4)	     
C                  fonction de phase polarisee: Qmie(mu).
C    UMIE(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)    (R4)     
C                  fonction de phase polarisee: Umie(mu).       
C      
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     Affichage d'un message a l'ecran, arret du programme et 
C     retour du status 1 au shell
C
C==============================================================================
      SUBROUTINE SOS_FPHASE_MIE(MIE_NBMU,RMU,RA,IA,RB,IB,ALPHA,
     *                          KMA1,KMA2,G,N2,TRACE,IER)
 

      IMPLICIT NONE
      
C* Definition des variables
C----------------------------------------------------------------------------
      DOUBLE PRECISION RMU(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX) 
                                        !Cosinus des angles de Gauss 	
      
      DOUBLE PRECISION ALPHA       !Valeur du parametre de taille. 
      DOUBLE PRECISION KMA1	   !Coefficient d'efficacite en extinction.
      DOUBLE PRECISION KMA2        !Coefficient d'efficacite en diffusion.
      DOUBLE PRECISION G 	   !Facteur d'anisotropie.


      DOUBLE PRECISION RA(0:SOS_MIE_DIM)     !Partie reelle de la fonction An(alpha). 
      DOUBLE PRECISION IA(0:SOS_MIE_DIM)     !Partie imaginaire de An(alpha).
      DOUBLE PRECISION RB(0:SOS_MIE_DIM)     !Partie reelle de la fonction Bn(alpha).
      DOUBLE PRECISION IB(0:SOS_MIE_DIM)     !Partie imaginaire de Bn(alpha).
      DOUBLE PRECISION AR	!Valeur du tableau RA pour n donne.
      DOUBLE PRECISION AI       !Valeur du tableau IA pour n donne.
      DOUBLE PRECISION BR       !Valeur du tableau RB pour n donne.
      DOUBLE PRECISION BI       !Valeur du tableau IB pour n donne.   	      
            
      DOUBLE PRECISION COEF
      DOUBLE PRECISION X        !Valeur mu = -cos(tetad) 
      				     !avec tetad l'angle de diffusion.
				           
      DOUBLE PRECISION PIM      !Fonction PI(mu) d'ordre n-1.
      DOUBLE PRECISION PIV      !Fonction PI(mu) d'ordre n.
      DOUBLE PRECISION PIP     	!Fonction PI(mu) d'ordre n+1.	
      DOUBLE PRECISION TAU      !Fonction TAU(mu) d'ordre n.
      
      DOUBLE PRECISION RES1     !Partie reelle de la fonction S1(alpha,m,mu).
      DOUBLE PRECISION IMS1     !Partie imaginaire de S1(alpha,m,mu).
      DOUBLE PRECISION RES2     !Partie reelle de la fonction S2(alpha,m,mu).
      DOUBLE PRECISION IMS2     !Partie imaginaire de S2(alpha,m,mu).      
      
      DOUBLE PRECISION Y1
      DOUBLE PRECISION Y2
      DOUBLE PRECISION Y3
      DOUBLE PRECISION Y4
      DOUBLE PRECISION AY4	!Facteur d'anisotropie.
      
      REAL*4 AY1		!Valeur du parametre de taille ALPHA.
      REAL*4 AY2		!Coefficient d'efficacite en extinction.
      REAL*4 AY3		!Coefficient d'efficacite en diffusion.
      REAL*4 IMIE(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)     
                             !fonction de phase : Imie(mu) pour ALPHA.
      REAL*4 QMIE(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)     
                             !fonction de phase polarisee: Qmie(mu).
      REAL*4 UMIE(-SOS_MIE_NBMU_MAX:SOS_MIE_NBMU_MAX)
                             !fonction de phase polarisee: Umie(mu).       

      LOGICAL TRACE          ! = vrai,si ecriture dans le fichier trace

      INTEGER*4 MIE_NBMU     ! Nombre d'angles effectifs des simulations         
      INTEGER*4 N2      !Ordre maximal de sommation sur N 
      INTEGER*4 J
      INTEGER*4 N
      
      INTEGER*4 IER        ! code d'erreur =0 si pas d'erreur, =-1 sinon           


      AY1=ALPHA           
      AY2=KMA1      
      AY3=KMA2    
      AY4=G
      COEF=2.D+00/KMA2/ALPHA**2
      
      DO 1 J=-MIE_NBMU,MIE_NBMU              
         X=-RMU(J)            
         PIM=0.D+00 
         PIV=1.D+00  
         TAU=X          
         RES1=0.D+00
         RES2=0.D+00
         IMS1=0.D+00
         IMS2=0.D+00
	 
         DO 2 N=1,N2           
            AI=IA(N)
            BI=IB(N)
            AR=RA(N)
            BR=RB(N)
            RES1=RES1-AI*PIV-BI*TAU  
            RES2=RES2+AI*TAU+BI*PIV  
            IMS1=IMS1+AR*PIV+BR*TAU 
            IMS2=IMS2-AR*TAU-BR*PIV  
            PIP=((2.D+00*N+1.D+00)*X*PIV-(N+1.D+00)*PIM)/N 
            PIM=PIV
            PIV=PIP  
            TAU=(N+1.D+00)*X*PIV-(N+2.D+00)*PIM 
    2    CONTINUE 
                                                                     
         Y1=RES1*RES1+IMS1*IMS1  
         Y2=RES2*RES2+IMS2*IMS2 
         Y3=2.D+00*RES2*RES1   
         Y4=2.D+00*IMS2*IMS1
         IMIE(J)=COEF*(Y1+Y2)      
         QMIE(J)=COEF*(Y2-Y1)
         UMIE(J)=COEF*(Y3+Y4)
    1 CONTINUE    
    
C* Enregistrement des resultats pour ALPHA
C*    Enregistrement des tableaux de fonctions de phase 
C*    pour des bornes d'enregistrement limitees au dimensionnement 
C*    utile -SOS_MIE_NBMU:SOS_MIE_NBMU
C------------------------------------------------------
      WRITE(1,err=921) AY1,AY2,AY3,AY4,
     &                (IMIE(J),J=-MIE_NBMU,MIE_NBMU),
     &                (QMIE(J),J=-MIE_NBMU,MIE_NBMU),
     &                (UMIE(J),J=-MIE_NBMU,MIE_NBMU)
      
C* Fichier trace :alpha,Qext,Qsca,G
C------------------------------------------------------  
      IF (TRACE) THEN
	   WRITE(99,20,err=923) AY1,AY2,AY3,AY4
      END IF


C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur et retour du status -1 au progamme appelant
C----------------------------------------------------------
  921 WRITE(6,*) 'Erreur d''ecriture sur le fichier MIE'
      IER=-1 
      GOTO 9999
  923 WRITE(6,*) 'Erreur d''ecriture sur le fichier trace MIE'
      IER=-1 
      GOTO 9999

C* Format 
C----------        
20    FORMAT(1X,E15.7,5X,E15.7,5X,E15.7,5X,E15.7)  

 9999 CONTINUE
 
      RETURN                                                                    
      END         !FIN DE LA PROCEDURE SOS_FPHASE_MIE       

