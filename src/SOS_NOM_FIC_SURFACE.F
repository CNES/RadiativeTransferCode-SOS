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
C* FICHIER: SOS_NOM_FIC_SURFACE.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Ce module gere la nomenclature des fichiers de BRDF / BPDF pour les OS.
C*
C* AUTEUR: 
C*    CS SI, B.Lafrance.
C*
C* DATE: 30/04/02
C*
C* MOD:VERSION:1.0: Creation
C* MOD:VERSION:1.1: Correction du champ pour l'angle limite solaire des calculs
C*                  du modele de Roujean (VS devient SL : angle Solaire Limite).
C* MOD:VERSION:1.2: Mise a jour des commentaires.
C*
C* MOD:VERSION:2.0: 10/01/2010 
C*     - Creation de la procedure SOS_COMPLETE_NOM_FIC_SURF
C*     - Ajout et gestion des arguments d'entree NBMU_GAUSS et FICANGLES_USER 
C*       et appel de la procedure SOS_COMPLETE_NOM_FIC_SURF par les procedures
C*          SOS_NOM_FICGLITTER, SOS_NOM_FICROUJ, SOS_NOM_FICRH, SOS_NOM_FICBREON, 
C*          SOS_NOM_FICNADAL
C*   
C* MOD:VERSION:2.1: 08/12/2015
C*    - Ajustements mineurs pour le respect strict des 72 colonnes 
C*      (incluant au décodage des constante) : requis pour compilation gfortran
C******************************************************************************


C==============================================================================
C PROCEDURE: SOS_NOM_FICGLITTER
C ==========
C      Cette procedure genere le nom du fichier de Glitter en fonction de :
C	- l'indice de l'eau : IND
C      	- la vitesse du vent : WIND
C       - le nombre d'angles de Gauss utilises NBMU_GAUSS
C       - le nom du fichier d'angles utilisateur FICANGLES_USER
C          (en supprimant son arborescence) 
C
C      On impose d'avoir :
C         IND au format F5.3.
C         WIND au format F4.1, 
C
C      Le nom du fichier de Glitter doit etre sous la forme:
C
C      GLITTER-IND-WIND-MUnb-NomFicUser avec 
C      - IND au format F5.3 (ex : 1.500) 
C      - WIND au format F4.1 (ex : 12.5) 
C      - MUnb le nombre d'angles de Gauss (ex : MU24)
C      - NomFicUser le nom du fichier utilisateur (sans arborescence) (ex : MyAngles)
C
C      On convertit WIND en une chaine de caracteres CWIND(1:4) 
C      Pour cela on passe du Format F4.1 * 10 => A4
C           Le dernier caractere CRN(4) est la decimale de WIND,
C           etc.
C
C Description des parametres entree / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C       FICANGLES_USER (CHARACTER*SOS_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C	WIND (double)	: Vitesse du vent m/s
C	IND (double)	: Indice de l'eau
C	FICS (CHARACTER*SOS_LENFIC1) : (S) nom du fichier de GLITTER (sans arborescence)	
C
C 
C Constantes utilisees:
C --------------------
C    SOS_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C    SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     Aucun
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel a la routine:
C     - SOS_COMPLETE_NOM_FIC_SURF
C
C==============================================================================
      SUBROUTINE SOS_NOM_FICGLITTER(NBMU_GAUSS,FICANGLES_USER,
     &                              IND,WIND,FICS) 

#include "SOS.h"

      IMPLICIT NONE

C* Declaration des variables
C  -------------------------         
      DOUBLE PRECISION WIND 	  ! (E) Vitesse du vent m/s
      DOUBLE PRECISION IND	  ! (E) Indice de l'eau.
				
      CHARACTER*SOS_LENFIC2 FICANGLES_USER  ! (E) chemin et nom du fichier des angles 
                                            ! "utilisateurs"
      
      CHARACTER*SOS_LENFIC1 FICS  ! (S) nom du fichier GLITTER 	
          
      CHARACTER*3 CWIND		  ! Chaines de caracteres
      CHARACTER*4 CIND
      
      INTEGER*2 I   
      INTEGER*4 NBMU_GAUSS  ! (E) Nombre d'angles de Gauss utilises 
      			      	 
	    
C --------- Conversion en chaines de caracteres des parametres d'entree	   
      WRITE(CIND,'(I4)') int(IND*1000)		
      WRITE(CWIND,'(I3)') int(WIND*10)  		

C --------- Les blancs sont remplaces par des 0	      
      DO I=1,4
         IF( CIND(I:I).EQ.' ') CIND(I:I) = '0'
      ENDDO   

      DO I=1,3
	 IF( CWIND(I:I).EQ.' ') CWIND(I:I) = '0'
      ENDDO
      
     
      FICS = 'GLITTER-'
           
C -------- Construction de la chaine de caracteres pour l'indice IND	    
      FICS(9:9) = CIND(1:1)     ! Chiffre des unités
      FICS(10:10) = '.'	        ! Virgule
      FICS(11:13) = CIND(2:4)   ! Trois premieres decimales   	    
      FICS(14:14) = '-'

C -------- Construction de la chaine de caracteres pour l'indice WIND	    
      FICS(15:16) = CWIND(1:2)     ! Chiffres des unites
      FICS(17:17) = '.'	           ! Virgule
      FICS(18:18) = CWIND(3:3)     ! premiere decimale   
      FICS(19:19) = '-'

C -------- Completude de la fin du fichier avec MUnb-NomFicUser	              
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)
      	  
     	      
      END !SOS_NOM_FICGLITTER
 



C==============================================================================
C PROCEDURE: SOS_NOM_FICROUJ
C ==========
C      Cette procedure genere le nom du fichier du modele de ROUJEAN 
C      en fonction des parametres K0,K1 et K2
C
C      On impose d'avoir : 
C	  kx au format F7.3.
C
C      Le nom du fichier de Roujean doit etre sous la forme:
C
C      ROUJ-K0-K1-K2-VL{SOS_TETAV_LIM_ROUJEAN}-
C          -SL{SOS_TETAS_LIM_ROUJEAN}-MUnb-NomFicUser avec  
C      - kx au format F5.3 (ex : 0.500) 
C      - SOS_TETAV_LIM_ROUJEAN et SOS_TETAS_LIM_ROUJEAN respectivement
C        l'angle zenithal de visee et solaire limites pour les calculs de BRDF
C        (en degres).
C      - MUnb le nombre d'angles de Gauss (ex : MU24)
C      - NomFicUser le nom du fichier utilisateur (sans arborescence) (ex : MyAngles)
C
C      On convertit Kx en une chaine de caracteres CKx(1:6) 
C      Pour cela on passe du Format F7.3 * 1000 => A6
C        Le premier caractere CKx(1) est le chiffre des centaines de Kx,
C        Le second caractere CKx(2) est le chiffre des dizaines de Kx,
C        Le troisieme caractere CKx(3) est le chiffre des unites de Kx,
C           etc.
C
C
C Description des parametres entree / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C       FICANGLES_USER (CHARACTER*SOS_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C	K0  (double)	: Parametre du modele de Roujean
C	K1  (double)	: Parametre du modele de Roujean
C	K2  (double)	: Parametre du modele de Roujean
C	FICS (CHARACTER*SOS_LENFIC1) : (S) nom du fichier de Roujean (sans arborescence)	
C
C 
C Constantes utilisees:
C --------------------
C     SOS_TETAV_LIM_ROUJEAN : Angle zenithal de visee maximal pour les calculs
C                             de BRDF.
C     SOS_TETAS_LIM_ROUJEAN : Angle zenithal solaire maximal pour les calculs
C                             de BRDF.
C     SOS_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C     SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     Aucun
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel a la routine:
C     - SOS_COMPLETE_NOM_FIC_SURF
C
C==============================================================================
      SUBROUTINE SOS_NOM_FICROUJ(NBMU_GAUSS,FICANGLES_USER,K0,K1,K2,
     &                           FICS) 

#include "SOS.h"

      IMPLICIT NONE

C* Declaration des variables
C  -------------------------              
      DOUBLE PRECISION K0,K1,K2    ! (E) Parametres du modele de Roujean

      CHARACTER*SOS_LENFIC2 FICANGLES_USER  ! (E) chemin et nom du fichier  
                                            ! des angles "utilisateurs"					
      CHARACTER*SOS_LENFIC1 FICS   ! (S) Nom du fichier de ROUJEAN 	
          
      CHARACTER*6 CK0,CK1,CK2	   ! Chaines de caracteres
      CHARACTER*4 CVLIM,CSLIM	   ! Chaines de caracteres
      
      INTEGER*2 I
      INTEGER*4 NBMU_GAUSS  ! (E) Nombre d'angles de Gauss utilises 

	    
C --------- Conversion en chaines de caracteres des parametres d'entree	
      WRITE(CK0,'(I6)') int(K0*1000)		
      WRITE(CK1,'(I6)') int(K1*1000)
      WRITE(CK2,'(I6)') int(K2*1000)
      WRITE(CVLIM,'(I2)') int(SOS_TETAV_LIM_ROUJEAN)
      WRITE(CSLIM,'(I2)') int(SOS_TETAS_LIM_ROUJEAN)
      
C --------- Les blancs sont remplaces par des 0	      
      DO I=1,6
         IF( CK0(I:I).EQ.' ') CK0(I:I) = '0'
         IF( CK1(I:I).EQ.' ') CK1(I:I) = '0'   
         IF( CK2(I:I).EQ.' ') CK2(I:I) = '0'
      ENDDO   
      DO I=1,2
         IF( CVLIM(I:I).EQ.' ') CVLIM(I:I) = '0'
         IF( CSLIM(I:I).EQ.' ') CSLIM(I:I) = '0'
      ENDDO 
      
      FICS = 'ROUJ-'

C -------- Construction de la chaine de caracteres pour K0	    
      FICS(6:8) = CK0(1:3)       ! Chiffres des unités
      FICS(9:9) = '.'	         ! Virgule
      FICS(10:12) = CK0(4:6)     ! Trois premieres decimales   	    
      FICS(13:13) = '-'

C -------- Construction de la chaine de caracteres pour K1	    
      FICS(14:16) = CK1(1:3)     ! Chiffres des unités
      FICS(17:17) = '.'	         ! Virgule
      FICS(18:20) = CK1(4:6)     ! Trois premieres decimales   	    
      FICS(21:21) = '-'

C -------- Construction de la chaine de caracteres pour K2	    
      FICS(22:24) = CK2(1:3)     ! Chiffres des unités
      FICS(25:25) = '.'	         ! Virgule
      FICS(26:28) = CK2(4:6)     ! Trois premieres decimales   	    
      FICS(29:29) = '-'

C -------- Construction de la chaine de caracteres pour la limite
C -------- sur l'angle de visee.
      FICS(30:31) = 'VL'	    
      FICS(32:33) = CVLIM(1:2)     ! Chiffre des unités
      FICS(34:34) = '-'	        
 
C -------- Construction de la chaine de caracteres pour la limite
C -------- sur l'angle solaire.
      FICS(35:36) = 'SL'	    
      FICS(37:38) = CSLIM(1:2)     ! Chiffre des unités
      FICS(39:39) = '-'	        
                        
C -------- Completude de la fin du fichier avec MUnb-NomFicUser	                  
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)
      
      END !SOS_NOM_FICROUJ
      


      
C==============================================================================
C PROCEDURE: SOS_NOM_FICRH
C ==========
C      Cette procedure genere le nom du fichier de Rondeaux - Herman
C      en fonction de :
C	- l'indice de vegetation : IND
C       - Du modele de BRDF associe pour la composante d'intensite :
C           modele de ROUJEAN en fonction des parametres K0,K1 et K2.
C
C      On impose d'avoir :
C         IND au format F5.3.
C	  kx au format F7.3.
C
C      Le nom du fichier de RH doit etre sous la forme:
C
C      RH-IND-ROUJ-K0-K1-K2-VL{SOS_TETAV_LIM_ROUJEAN}-
C        -SL{SOS_TETAS_LIM_ROUJEAN}-MUnb-NomFicUser 
C      avec
C      - IND au format F5.3 (ex : 1.500) 
C      - kx au format F7.3 (ex : 0.500) 
C      - SOS_TETAV_LIM_ROUJEAN et SOS_TETAS_LIM_ROUJEAN respectivement
C        l'angle zenithal de visee et solaire limites pour les calculs de BRDF
C        (en degres).
C      - MUnb le nombre d'angles de Gauss (ex : MU24)
C      - NomFicUser le nom du fichier utilisateur (sans arborescence) (ex : MyAngles)
C
C      On convertit IND en une chaine de caracteres IND(1:4) 
C      Pour cela on passe du Format 5.3 * 1000 => A4
C           Le premier caractere CIND(1) est le chiffre unite de IND,
C           etc.
C      On convertit Kx en une chaine de caracteres CKx(1:6) 
C      Pour cela on passe du Format F7.3 * 1000 => A6
C        Le premier caractere CKx(1) est le chiffre des centaines de Kx,
C        Le second caractere CKx(2) est le chiffre des dizaines de Kx,
C        Le troisieme caractere CKx(3) est le chiffre des unites de Kx,
C           etc.
C
C
C Description des parametres entree / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C       FICANGLES_USER (CHARACTER*SOS_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C	IND (double)	: Indice de vegetation
C	K0  (double)	: Parametre du modele de Roujean
C	K1  (double)	: Parametre du modele de Roujean
C	K2  (double)	: Parametre du modele de Roujean
C	FICS (CHARACTER*SOS_LENFIC1) : (S) nom du fichier de RH (sans arborescence)	
C
C 
C Constantes utilisees:
C --------------------
C     SOS_TETAV_LIM_ROUJEAN : Angle zenithal de visee maximal pour les calculs
C                             de BRDF.
C     SOS_TETAS_LIM_ROUJEAN : Angle zenithal solaire maximal pour les calculs
C                             de BRDF.
C     SOS_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C     SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     Aucun
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel a la routine:
C     - SOS_COMPLETE_NOM_FIC_SURF
C
C==============================================================================
      SUBROUTINE SOS_NOM_FICRH(NBMU_GAUSS,FICANGLES_USER,IND,K0,K1,K2,
     &                         FICS) 

#include "SOS.h"

      IMPLICIT NONE

C* Declaration des variables
C  -------------------------         
      DOUBLE PRECISION IND		! (E) Indice de vegetation.
      DOUBLE PRECISION K0,K1,K2 	! (E) Parametres du modele de Roujean.	
      CHARACTER*SOS_LENFIC2 FICANGLES_USER  ! (E) chemin et nom du fichier des angles 
                                            ! "utilisateurs"			
      CHARACTER*SOS_LENFIC1 FICS	! (S) nom du fichier RH. 	

      CHARACTER*6 CK0,CK1,CK2	! Chaines de caracteres
      CHARACTER*4 CVLIM,CSLIM	! Chaines de caracteres          
      CHARACTER*4 CIND 		! Chaine de caractere de l'indice de vegetation.
      
      INTEGER*2 I
      INTEGER*4 NBMU_GAUSS  ! (E) Nombre d'angles de Gauss utilises 
      			      	 

	    
C --------- Conversion en chaines de caracteres des parametres d'entree	   		
      WRITE(CIND,'(I4)') int(IND*1000)		
      WRITE(CK0,'(I6)') int(K0*1000)		
      WRITE(CK1,'(I6)') int(K1*1000)
      WRITE(CK2,'(I6)') int(K2*1000)
      WRITE(CVLIM,'(I2)') int(SOS_TETAV_LIM_ROUJEAN)
      WRITE(CSLIM,'(I2)') int(SOS_TETAS_LIM_ROUJEAN)
      
C --------- Les blancs sont remplaces par des 0	      
      DO I=1,4
	 IF( CIND(I:I).EQ.' ') CIND(I:I) = '0'
      ENDDO   
      DO I=1,6
         IF( CK0(I:I).EQ.' ') CK0(I:I) = '0'
         IF( CK1(I:I).EQ.' ') CK1(I:I) = '0'   
         IF( CK2(I:I).EQ.' ') CK2(I:I) = '0'
      ENDDO   
      DO I=1,2
         IF( CVLIM(I:I).EQ.' ') CVLIM(I:I) = '0'
         IF( CSLIM(I:I).EQ.' ') CSLIM(I:I) = '0'
      ENDDO 
  
      
      FICS = 'RH-'
      
C -------- Construction de la chaine de caracteres pour l'indice IND	    
      FICS(4:4) = CIND(1:1)     ! Chiffre des unités
      FICS(5:5) = '.'	        ! Virgule
      FICS(6:8) = CIND(2:4)     ! Trois premieres decimales   	    

C -------- Ajout du nom du fichier de BRDF associe
      FICS(9:14) = '-ROUJ-'

C -------- Construction de la chaine de caracteres pour K0	    
      FICS(15:17) = CK0(1:3)     ! Chiffres des unités
      FICS(18:18) = '.'	         ! Virgule
      FICS(19:21) = CK0(4:6)     ! Trois premieres decimales   	    
      FICS(22:22) = '-'

C -------- Construction de la chaine de caracteres pour K1	    
      FICS(23:25) = CK1(1:3)     ! Chiffres des unités
      FICS(26:26) = '.'	         ! Virgule
      FICS(27:29) = CK1(4:6)     ! Trois premieres decimales   	    
      FICS(30:30) = '-'

C -------- Construction de la chaine de caracteres pour K2	    
      FICS(31:33) = CK2(1:3)     ! Chiffres des unités
      FICS(34:34) = '.'	         ! Virgule
      FICS(35:37) = CK2(4:6)     ! Trois premieres decimales   	    
      FICS(38:38) = '-'

C     -------- Construction de la chaine de caracteres pour la limite
C     -------- sur l'angle de visee.
      FICS(39:40) = 'VL'	    
      FICS(41:42) = CVLIM(1:2)     ! Chiffre des unités
      FICS(43:43) = '-'	        
 
C     -------- Construction de la chaine de caracteres pour la limite
C     -------- sur l'angle solaire.
      FICS(44:45) = 'SL'	    
      FICS(46:47) = CSLIM(1:2)     ! Chiffre des unités
      FICS(48:48) = '-'	 
                   
C -------- Completude de la fin du fichier avec MUnb-NomFicUser	                  
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)

      END !SOS_NOM_FICRH      
      


      
      
C==============================================================================
C PROCEDURE: SOS_NOM_FICBREON
C ==========
C      Cette procedure genere le nom du fichier de BREON en fonction de :
C	- l'indice de sol : IND
C       - Du modele de BRDF associe pour la composante d'intensite :
C           modele de ROUJEAN en fonction des parametres K0,K1 et K2.
C
C      On impose d'avoir :
C         IND au format F5.3.
C	  kx au format F7.3.
C
C      Le nom du fichier de Breon doit etre sous la forme:
C
C      BREON-IND-ROUJ-K0-K1-K2-VL{SOS_TETAV_LIM_ROUJEAN}-
C        -SL{SOS_TETAS_LIM_ROUJEAN}-MUnb-NomFicUser
C      avec
C      - IND au format F5.3 (ex : 1.500) 
C      - kx au format F7.3 (ex : 0.500) 
C      - SOS_TETAV_LIM_ROUJEAN et SOS_TETAS_LIM_ROUJEAN respectivement
C        l'angle zenithal de visee et solaire limites pour les calculs de BRDF
C        (en degres).
C      - MUnb le nombre d'angles de Gauss (ex : MU24)
C      - NomFicUser le nom du fichier utilisateur (sans arborescence) (ex : MyAngles)
C
C      On convertit IND en une chaine de caracteres IND(1:4) 
C      Pour cela on passe du Format 5.3 * 1000 => A4
C           Le premier caractere CIND(1) est le chiffre unite de IND,
C           etc.
C      On convertit Kx en une chaine de caracteres CKx(1:6) 
C      Pour cela on passe du Format F7.3 * 1000 => A6
C        Le premier caractere CKx(1) est le chiffre des centaines de Kx,
C        Le second caractere CKx(2) est le chiffre des dizaines de Kx,
C        Le troisieme caractere CKx(3) est le chiffre des unites de Kx,
C           etc.
C
C Description des parametres entree / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C       FICANGLES_USER (CHARACTER*SOS_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C	IND (double)	: Indice de sol
C	K0  (double)	: Parametre du modele de Roujean
C	K1  (double)	: Parametre du modele de Roujean
C	K2  (double)	: Parametre du modele de Roujean
C	FICS (CHARACTER*SOS_LENFIC1) : (S) nom du fichier de Breon (sans arborescence)	
C
C 
C Constantes utilisees:
C --------------------
C     SOS_TETAV_LIM_ROUJEAN : Angle zenithal de visee maximal pour les calculs
C                             de BRDF.
C     SOS_TETAS_LIM_ROUJEAN : Angle zenithal solaire maximal pour les calculs
C                             de BRDF.
C     SOS_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C     SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     Aucun
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel a la routine:
C     - SOS_COMPLETE_NOM_FIC_SURF
C
C==============================================================================
      SUBROUTINE SOS_NOM_FICBREON(NBMU_GAUSS,FICANGLES_USER,IND,
     &                            K0,K1,K2,FICS) 

#include "SOS.h"

      IMPLICIT NONE

C* Declaration des variables
C  -------------------------         
      DOUBLE PRECISION IND		! (E) Indice de sol.
      DOUBLE PRECISION K0,K1,K2 	! (E) Parametres du modele de Roujean.	
     
      CHARACTER*SOS_LENFIC2 FICANGLES_USER  ! (E) chemin et nom du fichier des angles "utilisateurs"				
      CHARACTER*SOS_LENFIC1 FICS	! (S) nom du fichier Breon. 	

      CHARACTER*6 CK0,CK1,CK2	! Chaines de caracteres
      CHARACTER*4 CVLIM,CSLIM	! Chaines de caracteres          
      CHARACTER*4 CIND 		! Chaine de caractere de l'indice de vegetation.
      
      INTEGER*2 I
      INTEGER*4 NBMU_GAUSS  ! (E) Nombre d'angles de Gauss utilises 
	    
C --------- Conversion en chaines de caracteres des parametres d'entree	   		
      WRITE(CIND,'(I4)') int(IND*1000)		
      WRITE(CK0,'(I6)') int(K0*1000)		
      WRITE(CK1,'(I6)') int(K1*1000)
      WRITE(CK2,'(I6)') int(K2*1000)
      WRITE(CVLIM,'(I2)') int(SOS_TETAV_LIM_ROUJEAN)
      WRITE(CSLIM,'(I2)') int(SOS_TETAS_LIM_ROUJEAN)
      
C --------- Les blancs sont remplaces par des 0	   
      DO I=1,4
	 IF( CIND(I:I).EQ.' ') CIND(I:I) = '0'
      ENDDO      
      DO I=1,6
         IF( CK0(I:I).EQ.' ') CK0(I:I) = '0'
         IF( CK1(I:I).EQ.' ') CK1(I:I) = '0'   
         IF( CK2(I:I).EQ.' ') CK2(I:I) = '0'
      ENDDO   
      DO I=1,2
         IF( CVLIM(I:I).EQ.' ') CVLIM(I:I) = '0'
         IF( CSLIM(I:I).EQ.' ') CSLIM(I:I) = '0'
      ENDDO 
  
      
      FICS = 'BREON-'
      
C -------- Construction de la chaine de caracteres pour l'indice IND	    
      FICS(7:7) = CIND(1:1)      ! Chiffre des unités
      FICS(8:8) = '.'	         ! Virgule
      FICS(9:11) = CIND(2:4)     ! Trois premieres decimales   	    

C -------- Ajout du nom du fichier de BRDF associe
      FICS(12:17) = '-ROUJ-'

C -------- Construction de la chaine de caracteres pour K0	    
      FICS(18:20) = CK0(1:3)     ! Chiffres des unités
      FICS(21:21) = '.'	         ! Virgule
      FICS(22:24) = CK0(4:6)     ! Trois premieres decimales   	    
      FICS(25:25) = '-'

C -------- Construction de la chaine de caracteres pour K1	    
      FICS(26:28) = CK1(1:3)     ! Chiffres des unités
      FICS(29:29) = '.'	         ! Virgule
      FICS(30:32) = CK1(4:6)     ! Trois premieres decimales   	    
      FICS(33:33) = '-'

C -------- Construction de la chaine de caracteres pour K2	    
      FICS(34:36) = CK2(1:3)     ! Chiffres des unités
      FICS(37:37) = '.'	         ! Virgule
      FICS(38:40) = CK2(4:6)     ! Trois premieres decimales   	    
      FICS(41:41) = '-'

C     -------- Construction de la chaine de caracteres pour la limite
C     -------- sur l'angle de visee.
      FICS(42:43) = 'VL'	    
      FICS(44:45) = CVLIM(1:2)     ! Chiffre des unités
      FICS(46:46) = '-'	        
 
C     -------- Construction de la chaine de caracteres pour la limite
C     -------- sur l'angle solaire.
      FICS(47:48) = 'SL'	    
      FICS(49:50) = CSLIM(1:2)     ! Chiffre des unités
      FICS(51:51) = '-'	 
                   
C -------- Completude de la fin du fichier avec MUnb-NomFicUser	                  
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)


      END !SOS_NOM_FICBREON     
      


 
 
      
C==============================================================================
C PROCEDURE: SOS_NOM_FICNADAL
C ==========
C      Cette procedure genere le nom du fichier du modele de NADAL
C      en fonction de :
C	- la valeur de saturation alpha : ALPHA
C	- la valeur de beta : BETA
C       - l'indice de refraction surface / atmosphere : IND
C
C      On impose d'avoir : 
C         IND au format F5.3.
C         ALPHA au format F6.4,
C	  BETA au format F4.1.
C
C      Le nom du fichier de Nadal doit etre sous la forme:
C
C      NADAL-IND-K0-K1-K2-MUnb-NomFicUser avec  
C      - ALPHA au format F6.4 (ex : 0.0229) 
C      - BETA au format F3.1 (ex: 38)
C      - MUnb le nombre d'angles de Gauss (ex : MU24)
C      - NomFicUser le nom du fichier utilisateur (sans arborescence) (ex : MyAngles)
C
C      On convertit ALPHA en une chaine de caracteres CALPHA(1:6) 
C      Pour cela on passe du Format F6.4 * 10000 => A6
C      Le premier caractere CALPHA(1) est le chiffre unite de ALPHA,
C           etc.
C
C      On convertit BETA en une chaine de caracteres CBETA(1:4) 
C      Pour cela on passe du Format F4.1 * 10 => A4
C      Le premier caractere CBETA(1) est le chiffre des dizaines de BETA,
C           etc.
C
C Description des parametres entree / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C       FICANGLES_USER (CHARACTER*SOS_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C       IND (double)	: Indice de refraction surface / atmosphere
C	ALPHA (double)	: Parametre alpha du modelde de Nadal
C	BETA  (double)	: Parametre beta du modelde de Nadal
C	K0  (double)	: Parametre du modele de Roujean
C	K1  (double)	: Parametre du modele de Roujean
C	K2  (double)	: Parametre du modele de Roujean
C	FICS (CHARACTER*SOS_LENFIC1) : (S) nom du fichier de Nadal (sans arborescence)	
C
C 
C Constantes utilisees:
C --------------------
C     SOS_TETAV_LIM_ROUJEAN : Angle zenithal de visee maximal pour les calculs
C                             de BRDF.
C     SOS_TETAS_LIM_ROUJEAN : Angle zenithal solaire maximal pour les calculs
C                             de BRDF.
C     SOS_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C     SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     Aucun
C
C  Sous programmes utilises:
C --------------------------
C  Ce programme fait appel a la routine:
C     - SOS_COMPLETE_NOM_FIC_SURF
C
C==============================================================================
      SUBROUTINE SOS_NOM_FICNADAL(NBMU_GAUSS,FICANGLES_USER,IND,
     &                            ALPHA,BETA,K0,K1,K2,FICS) 

#include "SOS.h"

      IMPLICIT NONE

C* Declaration des variables
C  -------------------------  
      DOUBLE PRECISION IND	 ! (E) Indice de sol       
      DOUBLE PRECISION ALPHA	 ! (E) Parametre alpha du modelde de Nadal
      DOUBLE PRECISION BETA	 ! (E) Parametre beta du modelde de Nadal
      DOUBLE PRECISION K0,K1,K2  ! (E) Parametres du modele de Roujean.				
      
      CHARACTER*SOS_LENFIC2 FICANGLES_USER  ! (E) chemin et nom du fichier des angles "utilisateurs"					
      CHARACTER*SOS_LENFIC1 FICS ! (S) nom du fichier NADAL	
          
      CHARACTER*5 CALPHA	! Chaines de caracteres
      CHARACTER*3 CBETA 
      CHARACTER*4 CIND 
      CHARACTER*6 CK0,CK1,CK2	
      CHARACTER*4 CVLIM,CSLIM	
             
      INTEGER*2 I
      INTEGER*4 NBMU_GAUSS  ! (E) Nombre d'angles de Gauss utilises 
      


C --------- Conversion en chaines de caracteres des parametres d'entree
      WRITE(CIND,'(I4)') int(IND*1000)				
      WRITE(CALPHA,'(I5)') int(ALPHA*10000)
      WRITE(CBETA,'(I3)') int(BETA*10)	
      WRITE(CK0,'(I6)') int(K0*1000)		
      WRITE(CK1,'(I6)') int(K1*1000)
      WRITE(CK2,'(I6)') int(K2*1000)
      WRITE(CVLIM,'(I2)') int(SOS_TETAV_LIM_ROUJEAN)
      WRITE(CSLIM,'(I2)') int(SOS_TETAS_LIM_ROUJEAN)
            
    
C --------- Les blancs sont remplaces par des 0	 
      DO I=1,4
         IF( CIND(I:I).EQ.' ') CIND(I:I) = '0'
      ENDDO   
      	     
      DO I=1,5
         IF( CALPHA(I:I).EQ.' ') CALPHA(I:I) = '0'
      ENDDO         
      DO I=1,3
	 IF( CBETA(I:I).EQ.' ') CBETA(I:I) = '0'
      ENDDO   
      
      DO I=1,6
         IF( CK0(I:I).EQ.' ') CK0(I:I) = '0'
         IF( CK1(I:I).EQ.' ') CK1(I:I) = '0'   
         IF( CK2(I:I).EQ.' ') CK2(I:I) = '0'
      ENDDO   
      DO I=1,2
         IF( CVLIM(I:I).EQ.' ') CVLIM(I:I) = '0'
         IF( CSLIM(I:I).EQ.' ') CSLIM(I:I) = '0'
      ENDDO 



      
      FICS = 'NADAL-'
  
C -------- Construction de la chaine de caracteres pour l'indice IND	    
      FICS(7:7) = CIND(1:1)       ! Chiffre des unités
      FICS(8:8) = '.'	          ! Virgule
      FICS(9:11) = CIND(2:4)      ! Trois premieres decimales   	    
      FICS(12:12) = '-'
            
C -------- Construction de la chaine de caracteres pour alpha	    
      FICS(13:13) = CALPHA(1:1)     ! Chiffre des unites
      FICS(14:14) = '.'	            ! Virgule
      FICS(15:18) = CALPHA(2:5)     ! Quatre premieres decimales   	    
      FICS(19:19) = '-'

C -------- Construction de la chaine de caracteres pour beta	    
      FICS(20:21) = CBETA(1:2)      ! Chiffre des unites 
      FICS(22:22) = '.'
      FICS(23:23) = CBETA(3:3)      ! Chiffre des unites 

C -------- Ajout du nom du fichier de BRDF associe
      FICS(24:29) = '-ROUJ-'

C -------- Construction de la chaine de caracteres pour K0	    
      FICS(30:32) = CK0(1:3)     ! Chiffres des unités
      FICS(33:33) = '.'	         ! Virgule
      FICS(34:36) = CK0(4:6)     ! Trois premieres decimales   	    
      FICS(37:37) = '-'
		
C -------- Construction de la chaine de caracteres pour K1	    
      FICS(38:40) = CK1(1:3)     ! Chiffres des unités
      FICS(41:41) = '.'	         ! Virgule
      FICS(42:44) = CK1(4:6)     ! Trois premieres decimales   	    
      FICS(45:45) = '-'
		
C -------- Construction de la chaine de caracteres pour K2	    
      FICS(46:48) = CK2(1:3)     ! Chiffres des unités
      FICS(49:49) = '.'	         ! Virgule
      FICS(50:52) = CK2(4:6)     ! Trois premieres decimales   	    
      FICS(53:53) = '-'
		
C     -------- Construction de la chaine de caracteres pour la limite
C     -------- sur l'angle de visee.
      FICS(54:55) = 'VL'	    
      FICS(56:57) = CVLIM(1:2)     ! Chiffre des unités
      FICS(58:58) = '-'	        
 
C     -------- Construction de la chaine de caracteres pour la limite
C     -------- sur l'angle solaire.
      FICS(59:60) = 'SL'	    
      FICS(61:62) = CSLIM(1:2)     ! Chiffre des unités
      FICS(63:63) = '-'	 
      
    			                 
C -------- Completude de la fin du fichier avec MUnb-NomFicUser	          
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)
   
      END !SOS_NOM_FICNADAL      


 
      
C==============================================================================
C PROCEDURE: SOS_COMPLETE_NOM_FIC_SURF
C ==========

C
C Description des parametres entree / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C       FICANGLES_USER (CHARACTER*SOS_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C	FICS (CHARACTER*SOS_LENFIC1) : (E/S) nom du fichier de BRDF/BPDF (sans arborescence)
C 
C Constantes utilisees:
C --------------------
C    SOS_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C    SOS_LENFIC2 : Longueur des noms de fichiers avec arborescence.
C
C Common utilise:
C --------------
C     Aucun
C
C Cas d'erreur :
C ------------
C     Aucun
C==============================================================================
      SUBROUTINE SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,
     &                                     FICS)

#include "SOS.h"

      IMPLICIT NONE

C* Declaration des variables
C  -------------------------  
      CHARACTER*SOS_LENFIC2 FICANGLES_USER  ! (E) chemin et nom du fichier des angles "utilisateurs"	
      CHARACTER*SOS_LENFIC1 FICS            ! (E/S) nom du fichier GLITTER 	

      CHARACTER*4 CNBMU     ! Chaine de caracteres associee a NBMU_GAUSS
      
      INTEGER*2 I
      
      INTEGER*4 NBMU_GAUSS  ! (E) Nombre d'angles de Gauss utilises 
      INTEGER*4 IDEB, IDEB_CH, IFIN, IFIN_CH  ! position de debut et fin d'une chaine de caracteres  
      
	    
C --------- Conversion en chaines de caracteres de NBMU_GAUSS   
      WRITE(CNBMU,'(I4)') NBMU_GAUSS      
      

C -------- Construction de la chaine de caracteres pour MUgauss
C --------Structuration pour avoir NBMU_GAUSS a gauche (sans "blancs")	            
      IFIN=INDEX(FICS,' ')
      
      FICS(IFIN:IFIN+1) = 'MU'

      IF (NBMU_GAUSS.LT.10) THEN
          FICS(IFIN+2:IFIN+2) = CNBMU(4:4)	!MUn
      ELSE 
          IF (NBMU_GAUSS.LT.100) THEN
              FICS(IFIN+2:IFIN+3) = CNBMU(3:4)  !MUnn
	  ELSE 
	     IF (NBMU_GAUSS.LT.1000) THEN
                 FICS(IFIN+2:IFIN+4) = CNBMU(2:4) !MUnnn
	     ELSE !NBMU_GAUSS < 10000
	     	 FICS(IFIN+2:IFIN+5) = CNBMU(1:4) !MUnnnn
	     ENDIF !NBMU_GAUSS < 1000
	  ENDIF !NBMU_GAUSS < 100
      ENDIF !NBMU_GAUSS < 10

       
C -------- Construction de la chaine de caracteres pour NomFicUser
C -------- Structuration pour avoir NBMU_GAUSS a gauche (sans "blancs")	            

      IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN

C       Determination de la position du nom de fichier sans son arborescence          
	  IFIN_CH=INDEX(FICANGLES_USER,' ')
          IFIN_CH=IFIN_CH-1
	  IF (IFIN_CH.LE.0) IFIN_CH=SOS_LENDIR
	  
	  IDEB_CH=0
	  	  
10	  IDEB=INDEX(FICANGLES_USER(IDEB_CH:IFIN_CH),'/')
	  IDEB_CH=IDEB_CH+IDEB
          IF (IDEB.NE.0) GOTO 10	  

		
C       Ajout du nom de fichier d'angles "utilisateur" au nom du fichier     
          IFIN=INDEX(FICS,' ')  
	  FICS(IFIN:IFIN) = '-'
	  
	  IFIN=IFIN+1
	  FICS(IFIN:IFIN+IFIN_CH-IDEB_CH+1) 
     &       = FICANGLES_USER(IDEB_CH:IFIN_CH)
	  
      ENDIF 
   
      END !SOS_COMPLETE_NOM_FIC_SURF   
