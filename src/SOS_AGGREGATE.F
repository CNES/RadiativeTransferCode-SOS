C******************************************************************************
C* Copyright 2023, Centre National d'Etudes Spatiales (CNES)
C* 
C* This file is part of the SOS-ABS radiative transfer code.
C* 
C* SOS-ABS is free software: you can redistribute it and/or modify
C* it under the terms of the GNU General Public License as published by
C* the Free Software Foundation, either version 3 of the License, or
C* (at your option) any later version.
C* 
C* SOS-ABS is distributed in the hope that it will be useful,
C* but WITHOUT ANY WARRANTY; without even the implied warranty of
C* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C* GNU General Public License for more details.
C* 
C* You should have received a copy of the GNU General Public License
C* along with SOS-ABS. If not, see <http://www.gnu.org/licenses/>.
C******************************************************************************


C******************************************************************************
C* FICHIER: SOS_AGGREGATE.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Ponderation et aggregation des sorties de SOS.F pour chaque exponentielle de la méthode CKD
C*      - Aggregation des paramètres I Q U pour chaque ordre de la decomposition en serie de Fourier 
C*	 - Aggregation des transmissions atmosphériques
C*	 - Aggregation des épaisseurs optiques
C*
C* AUTEUR: 
C*   CS GROUP France, X.Lenot, B.Lafrance
C* DATE: 06/06/2013
C*
C* MOD:VERSION:6.0: Creation du fichier pour adapter les OS à la méthode CKD (absorption atmosphérique)
C*
C* MOD:VERSION:6.1: 01/09/2020
C*     Passage du nom du fichier temporaire FICOS_AGG_TMP pour permettre
C*     plusieurs lancement simultanés depuis le même répertoire.
C*
C* MOD:VERSION:6.2: 31/08/2022
C*    - Renommage des constantes de SOS_* en CTE_*
C*    - Suppression des arguments TRACE et IDLOG (non utilisés)
C*
C******************************************************************************

C----------------------------------------------------------------------------
C Definition des constantes  
C---------------------------------------------------------------------------- 
C Constantes utilisees :
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C    CTE_OS_NBMU_MAX : Nombre maximal d'angles positifs pour la resolution de l'ETR
C----------------------------------------------------------------------------
#include "SOS.h"
#define INCTE_PI DACOS(-1.D+00)

C==============================================================================
C PROGRAMME: SOS_AGGREGATE
C ==========      
C      Programme qui pondère les paramètres I, Q, U, les transmissions et épaisseurs optiques
C      calculées pour une exponentielle de la méthode CKD, et qui les aggrege avec les memes paramètres
C      calcules precedemment.
C
C Donnees en entrée du programme 
C ------------------------------
C    -->Les paramètres/arguments du programme
C
C             NBMU (I4) : (E) Nombre d'angles effectifs pour les calculs de luminances      			               
C
C             AIK (I4) : (E) Coefficient de ponderationde l'exponentielle CKD traitee      			               
C
C             FICOS_TMP (CHARACTER) : (E) Nom complet du fichier resultat des OS obtenus precedemment
C                                     (repertoire + nom fichier + extension)
C
C             TTOT_TRONC_TMP (DOUBLE) : (E) Epaisseur optique totale au sol 
C				            (obtenue precedemment) (tronquee)    
C
C             TTOT_VRAI_TMP (DOUBLE) : (E) Epaisseur optique totale au sol 
C				           (obtenue precedemment) (non tronquee)    
C
C             TAUOUT_TMP (DOUBLE) : (E) Epaisseur optique au niveau de sortie pour le champ montant
C			               (obtenue precedemment) (tronquee)    
C
C             TDIFMUS_TMP (DOUBLE) : (E) Transmission diffuse pour la direction d'incidence solaire
C      				      de l'atmosphere equivalente (obtenue precedemment)  
C
C             TDIFMUG_TMP (DOUBLE) : (E) Tableau 1D de 2*CTE_OS_NBMU_MAX + 1 elements
C			     	             contenant la transmission diffuse de l'atmosphere equivalente
C				             pour les angles de Gauss et les angles fournis par l'utilisateur
C   				         (obtenue precedemment)
C
C             EMOINS_TMP (DOUBLE) : (E) Flux diffus eclairant le sol, normalisé par le flux solaire incident au TOA.
C
C             EPLUS_TMP (DOUBLE) : (E)  Flux diffus emergeant au TOA, normalisé par le flux solaire incident au TOA.
C        
C             FICOS (CHARACTER) : (E) Nom complet du fichier resultat des OS (resultats aggreges)
C                                     (repertoire + nom fichier + extension)
C
C             TTOT_TRONC (DOUBLE) : (S) Epaisseur optique totale au sol 
C				           (resultats aggreges) (tronquee)    
C
C             TTOT_VRAI (DOUBLE) : (S) Epaisseur optique totale au sol 
C				          (resultats aggreges) (non tronquee)    
C
C             TAUOUT (DOUBLE) : (S) Epaisseur optique au niveau de sortie pour le champ montant
C			               (resultats aggreges) (tronquee)    
C
C             TDIFMUS (DOUBLE) : (S) Transmission diffuse pour la direction d'incidence solaire
C      				  de l'atmosphere equivalente (resultats aggreges)  
C
C             TDIFMUG (DOUBLE) : (S) Tableau 1D de 2*CTE_OS_NBMU_MAX + 1 elements
C			     	         contenant la transmission diffuse de l'atmosphere equivalente
C				         pour les angles de Gauss et les angles fournis par l'utilisateur
C   				         (resultats aggreges)
C
C             EMOINS (DOUBLE) : (S) Flux diffus eclairant le sol, normalisé par le flux solaire incident au TOA
C   				        (resultats aggreges)
C
C             EPLUS (DOUBLE) : (S) Flux diffus emergeant au TOA, normalisé par le flux solaire incident au TOA
C   				       (resultats aggreges)
C    
C      	IER (I4) : (S) code d'erreur =0 si pas d'erreur, =-1 sinon
C
C    -->Fichiers
C   Fichier resultat des OS (decomposition en series de Fourier)
C
C
C Resultats fournis par le programme
C ----------------------------------
C   Fichier resultat des OS (decomposition en series de Fourier)
C   Transmissions diffuses
C   Epaisseurs optiques (sol et au niveau de sortie)
C
C   Si un ou des fichiers de meme noms existent deja, ils sont ecrases
C
C
C 
C Description des fichiers utilises
C ---------------------------------
C
C   -->Contenu du fichier resultat des OS
C       Fichier binaire non formatte, cree par SOS_OS
C       Il contient les paramètres de Stokes Q , U et I pour un ordre S
C       du développement en series de Fourier, en fonction de l'angle mu
C       (valeurs en double précision).
C
C   -->Contenu du fichier temporaire FICOS_TMP
C       Identique au fichier de resultat des OS.
C       Fichier fourni en entrée avec les paramètres calcules precedemment
C
C
C Common utilise:
C --------------
C  Aucun
C
C
C Cas d'erreur :
C ------------        
C     -Decodage des paramètres incorrect
C     -Erreur dans un sous-programme
C     -Erreur a l'ouverture du fichier Trace
C     -Paramètres non valables
C
C     Affichage d'un message à l'écran, arret du programme et 
C     retour du status 1 au shell
C   
C
C  Sous programmes utilises:
C --------------------------
C  Aucun
C
C============================================================================

      SUBROUTINE SOS_AGGREGATE(NBMU, AIK, FICOS_TMP, TTOT_TRONC_TMP,
     &	 	    	      TTOT_VRAI_TMP, TAUOUT_TMP, TDIFMUS_TMP,
     &		    	      TDIFMUG_TMP, EMOINS_TMP, EPLUS_TMP,
     &                     FICOS_AGG_TMP, FICOS,
     &		    	      TTOT_TRONC, TTOT_VRAI, TAUOUT,            
     &	 	    	      TDIFMUS, TDIFMUG, EMOINS, EPLUS,
     &                     IER)
     

C* Definition des variables   
C*------------------------------------------------------
				
     
      DOUBLE PRECISION TDIFMUS  !Transmission diffuse pour la direction d'incidence solaire
      				    !de l'atmosphere equivalente
				
      DOUBLE PRECISION TDIFMUG(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				!Transmission diffuse pour les angles de Gauss de l'atmosphere equivalente

      DOUBLE PRECISION EMOINS   ! Flux diffus eclairant le sol, normalisé par
      				    ! le flux solaire incident au TOA.
                                
      DOUBLE PRECISION EPLUS    ! Flux diffus emergeant au TOA, normalisé par
      				    ! le flux solaire incident au TOA.
                                
      DOUBLE PRECISION TDIFMUS_TMP  !Transmission diffuse pour la direction d'incidence solaire
      				        !de l'atmosphere equivalente (calcul pour un septuor ik1*ik2*ik3*ik4*ik5*ik6*ik7)
				
      DOUBLE PRECISION TDIFMUG_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)  
      				!Transmission diffuse pour les angles de Gauss de l'atmosphere equivalente
				!  (calcul pour un octuor ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8 de décomposition CKD)

      DOUBLE PRECISION EMOINS_TMP   ! Flux diffus eclairant le sol, normalisé par
      				        ! le flux solaire incident au TOA.
                                    ! (calcul pour un octuor ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8)

      DOUBLE PRECISION EPLUS_TMP    ! Flux diffus emergeant au TOA, normalisé par
      				        ! le flux solaire incident au TOA.
                                    ! (calcul pour un octuor ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8)

      DOUBLE PRECISION TTOT_TRONC ! Epaisseur optique totale au sol (tronquee).
      DOUBLE PRECISION TTOT_VRAI  ! Epaisseur optique totale au sol (non tronquee).

      DOUBLE PRECISION TTOT_TRONC_TMP ! Epaisseur optique totale au sol (tronquee).
                                      ! (calcul pour un octuor ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8)

      DOUBLE PRECISION TTOT_VRAI_TMP  ! Epaisseur optique totale au sol (non tronquee).	
					   ! (calcul pour un octuor ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8)

      DOUBLE PRECISION TAUOUT	! Epaisseur optique au niveau de
      					! sortie des resultats pour le champ montant (tronquee).
      DOUBLE PRECISION TAUOUT_TMP 	! Epaisseur optique au niveau de
      				  	! sortie des resultats pour le champ montant (tronquee).
				 	! (calcul pour un octuor ik1*ik2*ik3*ik4*ik5*ik6*ik7*ik8)


      ! Paramètres de Stokes resultant de l'aggregation
      DOUBLE PRECISION I3_RES(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Paramètre de Stokes I pour un ordre S
      				! du développement en series de Fourier,
				! en fonction de l'angle mu.
      DOUBLE PRECISION Q3_RES(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Paramètre de Stokes Q.
      DOUBLE PRECISION U3_RES(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Paramètre de Stokes U.

      ! Paramètres de Stokes correspondant a l'iteration en cours
      DOUBLE PRECISION I3_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Paramètre de Stokes I pour un ordre S
      				! du développement en series de Fourier,
				! en fonction de l'angle mu.
      DOUBLE PRECISION Q3_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Paramètre de Stokes Q.
      DOUBLE PRECISION U3_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Paramètre de Stokes U.

      ! Paramètres de Stokes issus des iterations precedentes
      DOUBLE PRECISION I3_RES_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Paramètre de Stokes I pour un ordre S
      				! du développement en series de Fourier,
				! en fonction de l'angle mu.
      DOUBLE PRECISION Q3_RES_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Paramètre de Stokes Q.
      DOUBLE PRECISION U3_RES_TMP(-CTE_OS_NBMU_MAX:CTE_OS_NBMU_MAX)	
      				! Paramètre de Stokes U.

      DOUBLE PRECISION AIK  !Coefficient de ponderation

      DOUBLE PRECISION TRANS  !Transmission atmosphérique

      INTEGER*4 NBMU    ! nombre d'angles utiles		

      CHARACTER*CTE_LENFIC2 FICOS_TMP      ! Nom du fichier temporaire resultat des OS (calcul pour une itération)

      CHARACTER*CTE_LENFIC2 FICOS_AGG_TMP   ! Nom du fichier temporaire resultat des OS (calcul aggrégé)

      CHARACTER*CTE_LENFIC2 FICOS       ! Nom du fichier résultat final des OS (incluant toutes les aggrégations)

      INTEGER*4 IER        !code d'erreur =0 si pas d'erreur, =-1 sinon

      LOGICAL EX	! Verifie l'existence d'un fichier

      INTEGER*4 I,J   	! Indice de tableau

      INTEGER*4 IOS16,IOS17   	! Statuts de fichiers

      INTEGER*4 LFICOS ! longueur de la chaine de caractere FICOS
      INTEGER*4 LFICOSTMP ! longueur de la chaine de caractere FICOS_TMP
      INTEGER*4 LFICOSAGGTMP ! longueur de la chaine de caractere FICOS_AGG_TMP

      INTEGER*4 LCOM     	        ! longueur de la chaine de caractere COMMAND
      
      CHARACTER*CTE_LENCOM  COMMAND     ! Chaine de commande pour appel system.

 	  
      	        		 

C* Initiatlisation de valeurs
C ---------------------------

C*    Init des valeurs des statuts de fichiers
      IOS16=0
      IOS17=0


C*    Longueur des chaines de caracteres
      LFICOS=INDEX(FICOS,' ')-1 
      LFICOSTMP=INDEX(FICOS_TMP,' ')-1 
      LFICOSAGGTMP=INDEX(FICOS_AGG_TMP,' ')-1 

C* Controle de l'existence du fichier de Matrice de reflexion
C  ----------------------------------------------------------	
      INQUIRE(FILE=FICOS(1:LFICOS),EXIST=EX)

	  
C--------------------------------------------------------------
C
C       ******  AGGREGATION DES RESULTATS   ******
C
C--------------------------------------------------------------
      
     	 	
 
C************************************************************************ 
C* Paramètres de Stokes
C* Les valeurs IQU des decompositions en serie de Fourier sont lues
C* dans le fichier temporaire, pondérées par le coeff AIK de pondération
C* et additionnées aux valeurs IQU stockées dans le fichier résultat
C************************************************************************ 
 
C* Ouverture des fichiers a lire
C------------------------------------------------     
C     Ouverture du fichier temporaire des OS
      OPEN(UNIT=16,FILE=FICOS_TMP,FORM='UNFORMATTED',
     &     STATUS='OLD',ERR=1101) 

      IF (EX.EQV..FALSE.) THEN  

C* Cas d'une premiere iteration: le fichier resultat n'existe pas encore
C* On affecte la valeur IOSTAT17 du fichier de sortie a une valeur 
C* différente de 0 => le fichier ne sera pas lu
	
	  IOS17=-1

      ELSE

C* Ouverture du fichier resultat des OS
          OPEN(UNIT=17,FILE=FICOS,FORM='UNFORMATTED',
     &         STATUS='OLD',ERR=1201) 
      
      ENDIF 
      

C* Ouverture du fichier de sortie
C------------------------------------------------     
C     Ouverture du fichier temporaire de sortie
      OPEN(UNIT=18,FILE=FICOS_AGG_TMP,FORM='UNFORMATTED',
     &     STATUS='NEW',ERR=1301) 
     
C* Aggregation des paramètres IQU pour chaque ordre de la serie de Fourier 
C-------------------------------------------------------------------------
    
C*    Boucle sur l'ordre IS
      IS=0

C*    Initialisation des valeurs lues dans les fichiers
      DO J=-NBMU,NBMU 
      
	     Q3_TMP(J)=0.
	     U3_TMP(J)=0.
	     I3_TMP(J)=0.
	     Q3_RES_TMP(J)=0.
	     U3_RES_TMP(J)=0.
	     I3_RES_TMP(J)=0.
	  
      ENDDO
      
C*    Lecture du fichier TMP qui a ete enregistre avec des 
C*    bornes de tableaux limitées au dimensionnement 
C*    utile -NBMU:NBMU
   10 IF (IOS16.EQ.0) THEN
   
   	    READ(16,IOSTAT=IOS16,err=1102) (Q3_TMP(J),J=-NBMU,NBMU),
     & 		                           (U3_TMP(J),J=-NBMU,NBMU),
     &		                           (I3_TMP(J),J=-NBMU,NBMU)
 
      ENDIF
      
C*    Lecture du fichier Resultat qui a ete enregistre avec des 
C*    bornes de tableaux limitées au dimensionnement 
C*    utile -NBMU:NBMU
      IF (IOS17.EQ.0) THEN
      
      	    READ(17,IOSTAT=IOS17,err=1202) (Q3_RES_TMP(J),J=-NBMU,NBMU),
     & 		                     (U3_RES_TMP(J),J=-NBMU,NBMU),
     &		                            (I3_RES_TMP(J),J=-NBMU,NBMU)

      ELSE
             ! Cas où les deux fichiers sont termines
	     IF (IOS16.NE.0) GOTO 9
	  
      ENDIF
      
C*    Calcul de l'aggregation
      	
      DO J=-NBMU,NBMU 
      
            ! On ajoute la valeur pondérée du nouveau calcul 
	     ! au resultat à l'itération précédente
	     Q3_RES(J)=Q3_RES_TMP(J)+AIK*Q3_TMP(J)
	     U3_RES(J)=U3_RES_TMP(J)+AIK*U3_TMP(J)
	     I3_RES(J)=I3_RES_TMP(J)+AIK*I3_TMP(J)
	  
	     !Re-initialisation des tableaux au cas ou le fichier est fini
	     Q3_TMP(J)=0.
	     U3_TMP(J)=0.
	     I3_TMP(J)=0.
	     Q3_RES_TMP(J)=0.
	     U3_RES_TMP(J)=0.
	     I3_RES_TMP(J)=0.
	  
      ENDDO

C*    Ecriture des resultats

      WRITE(18,err=1302) (Q3_RES(J),J=-NBMU,NBMU),
     & 		    (U3_RES(J),J=-NBMU,NBMU),
     &		           (I3_RES(J),J=-NBMU,NBMU)
        
      IS=IS+1
      GOTO 10
      
 9    CLOSE(16) 
      CLOSE(17)
      CLOSE(18)
      
C --- Destruction du fichier resultat des OS et remplacement par le nouveau

C*    Destruction de l'ancien fichier

      OPEN(UNIT=17,FILE=FICOS,FORM='UNFORMATTED') 
      CLOSE(17,STATUS='DELETE')    

C*    Changement de nom du fichier resultat

      COMMAND='mv '//FICOS_AGG_TMP(1:LFICOSAGGTMP)//
     &	      ' '//FICOS(1:LFICOS)
      LCOM = INDEX(COMMAND(10:CTE_LENCOM),' ') - 1
      IF((LCOM.LE.0).OR.(LCOM.GT.CTE_LENCOM)) GOTO 1001
      CALL SYSTEM(COMMAND)
      

C*    Test d'existence du fichier
      INQUIRE(FILE=FICOS(1:LFICOS),EXIST=EX)
      IF (EX.EQV..FALSE.) GOTO 1002
      
C************************************************************************ 
C* Transmission atmosphérique
C* Simple aggregation de contenu avec ponderation de la transmission
C************************************************************************ 
      TDIFMUS=TDIFMUS+AIK*TDIFMUS_TMP

      EMOINS=EMOINS+AIK*EMOINS_TMP
      EPLUS=EPLUS+AIK*EPLUS_TMP
      
      DO J=-NBMU,NBMU 
      	TDIFMUG(J)=TDIFMUG(J)+AIK*TDIFMUG_TMP(J)
      ENDDO
      
C***************************************************************************************** 
C* Epaisseurs optiques
C* Il faut passer en transmissions (avec mus=1), faire l'aggregation, puis revenir en épaisseur optique
C***************************************************************************************** 
  
C*    Epaisseur optique totale tronquee
      IF (TTOT_TRONC.NE.0) THEN
      	   TRANS = AIK*DEXP(-TTOT_TRONC_TMP)+DEXP(-TTOT_TRONC)
      ELSE
          TRANS = AIK*DEXP(-TTOT_TRONC_TMP)
      ENDIF	
      TTOT_TRONC = -DLOG(TRANS)
      
C*    Epaisseur optique totale vraie
      IF (TTOT_VRAI.NE.0) THEN
      	   TRANS = AIK*DEXP(-TTOT_VRAI_TMP)+DEXP(-TTOT_VRAI)
      ELSE
          TRANS = AIK*DEXP(-TTOT_VRAI_TMP)
      ENDIF	
      TTOT_VRAI = -DLOG(TRANS)

C*    Epaisseur optique au niveau de sortie pour le champ montant
      IF (TAUOUT.NE.0) THEN
      	   TRANS = AIK*DEXP(-TAUOUT_TMP)+DEXP(-TAUOUT)
      ELSE
          TRANS = AIK*DEXP(-TAUOUT_TMP)
      ENDIF	
      TAUOUT = -DLOG(TRANS)

C* Fin nominale 
C-------------------
      GOTO 9999
      
C* Cas d'erreur
C-------------------

 1001 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1001 : '
      WRITE(6,*) '  => File copy command maybe too long:' 
      WRITE(6,*) '  => Check if CTE_LENCOM has to be adjusted.' 
      IER=-1      
      GOTO 9999                  	

 1002 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1002 : '
      WRITE(6,*) '  => Ouput OS file has note been created: ',FICOS 
      IER=-1      
      GOTO 9999                  	

 1101 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1101 : '
      WRITE(6,*) '  => Error opening file: ',FICOS_TMP 
      IER=-1      
      GOTO 9999                  	

 1102 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1102 : '
      WRITE(6,*) '  => Error reading file: ',FICOS_TMP 
      IER=-1      
      GOTO 9999                  	

 1201 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1201 : '
      WRITE(6,*) '  => Error opening file: ',FICOS 
      IER=-1      
      GOTO 9999                  	

 1202 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1202 : '
      WRITE(6,*) '  => Error reading file: ',FICOS 
      IER=-1      
      GOTO 9999                  	

 1301 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1301 : '
      WRITE(6,*) '  => Error opening file: ',FICOS_AGG_TMP 
      IER=-1      
      GOTO 9999                  	


 1302 WRITE(6,*) '  SOS_AGGREGATE : ERROR_1302 : '
      WRITE(6,*) '  => Error writing file: ',FICOS_AGG_TMP 
      IER=-1      
      GOTO 9999                  	

C* Sortie de la fonction

 9999 RETURN
 
      END	!FIN DU PROGRAMME SOS_AGGREGATE
