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
C* FICHIER: SOS_NOM_FIC_SURFACE.F
C* PROJET: Ordres successifs de diffusion
C* ROLE: Ce module gère la nomenclature des fichiers de BRDF / BPDF pour les OS.
C*
C* AUTEUR: 
C*    CS GROUP France, B.Lafrance.
C*
C* DATE: 30/04/2002
C* MOD:VERSION:1.0: Creation
C*
C* MOD:VERSION:1.1: Correction du champ pour l'angle limite solaire des calculs
C*                  du modèle de Roujean (VS devient SL : angle Solaire Limite).
C* MOD:VERSION:1.2: Mise à jour des commentaires.
C*
C* MOD:VERSION:2.0: 10/01/2010 
C*     - Creation de la procédure SOS_COMPLETE_NOM_FIC_SURF
C*     - Ajout et gestion des arguments d'entrée NBMU_GAUSS et FICANGLES_USER 
C*       et appel de la procédure SOS_COMPLETE_NOM_FIC_SURF par les procedures
C*          SOS_NOM_FICGLITTER, SOS_NOM_FICROUJ, SOS_NOM_FICRH, SOS_NOM_FICBREON, 
C*          SOS_NOM_FICNADAL
C*
C* MOD:VERSION:2.1: 24/02/2016 
C*    - Ajustements mineurs pour le respect strict des 72 colonnes 
C*
C* MOD:VERSION:2.2: 17/07/2020 
C*    - Modification de l'identifiant "ROUJ" en "ROUJEAN" 
C*      dans le nom des fichiers produits
C*    - Ajout de la routine SOS_NOM_FICMAIGNAN 
C*
C* MOD:VERSION:2.3: 07/07/2022
C*    - Renommage des constantes de SOS_* en CTE_*
C*
C*
C******************************************************************************
#include "SOS.h"

C==============================================================================
C PROCEDURE: SOS_NOM_FICGLITTER
C ==========
C      Cette procédure génère le nom du fichier de Glitter en fonction de :
C	- l'indice de réfraction de l'eau / air : IND
C      - la vitesse du vent : WIND
C      - le nombre d'angles de Gauss utilises NBMU_GAUSS
C      - le nom du fichier d'angles utilisateur FICANGLES_USER
C          (en supprimant son arborescence) 
C
C      On impose d'avoir :
C         IND au format F5.3.
C         WIND au format F4.1, 
C
C      Le nom du fichier de Glitter doit être sous la forme:
C
C      GLITTER-IND-WIND-MUnb-NomFicUser avec 
C      - IND au format F5.3 (ex : 1.500) 
C      - WIND au format F4.1 (ex : 12.5) 
C      - MUnb le nombre d'angles de Gauss (ex : MU24)
C      - NomFicUser le nom du fichier utilisateur (sans arborescence) (ex : MyAngles)
C
C      On convertit WIND en une chaine de caractères CWIND(1:4) 
C      Pour cela on passe du Format F4.1 * 10 => A4
C      Le dernier caractère CRN(4) est la decimale de WIND,
C      etc.
C
C Description des paramètres entrée / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C      FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C	WIND (double)	: Vitesse du vent m/s
C	IND (double)	: Indice de l'eau
C	FICS (CHARACTER*CTE_LENFIC1) : (S) nom du fichier de GLITTER (sans arborescence)	
C
C 
C Constantes utilisees:
C --------------------
C    CTE_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
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


      IMPLICIT NONE

C* Declaration des variables
C  -------------------------         
      DOUBLE PRECISION WIND 	  ! (E) Vitesse du vent m/s
      DOUBLE PRECISION IND	  ! (E) Indice de l'eau.
				
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! (E) chemin et nom du fichier des angles 
                                            ! "utilisateurs"
      
      CHARACTER*CTE_LENFIC1 FICS  ! (S) nom du fichier GLITTER 	
          
      CHARACTER*3 CWIND	! Chaines de caractères
      CHARACTER*4 CIND
      
      INTEGER*2 I   
      INTEGER*4 NBMU_GAUSS  ! (E) Nombre d'angles de Gauss utilises 
      			      	 
	    
C --------- Conversion en chaines de caractères des paramètres d'entrée	   
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
           
C -------- Construction de la chaine de caractères pour l'indice IND	    
      FICS(9:9) = CIND(1:1)     ! Chiffre des unités
      FICS(10:10) = '.'	    ! Virgule
      FICS(11:13) = CIND(2:4)   ! Trois premières décimales   	    
      FICS(14:14) = '-'

C -------- Construction de la chaine de caractères pour l'indice WIND	    
      FICS(15:16) = CWIND(1:2)     ! Chiffres des unites
      FICS(17:17) = '.'	       ! Virgule
      FICS(18:18) = CWIND(3:3)     ! premiere decimale   
      FICS(19:19) = '-'

C -------- Completude de la fin du fichier avec MUnb-NomFicUser	                  
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)
      	  
     	      
      END !SOS_NOM_FICGLITTER
 



C==============================================================================
C PROCEDURE: SOS_NOM_FICROUJ
C ==========
C      Cette procédure génère le nom du fichier du modèle de ROUJEAN 
C      en fonction des paramètres K0,K1 et K2
C
C      On impose d'avoir : 
C	  kx au format F7.3.
C
C      Le nom du fichier de Roujean doit être sous la forme:
C
C      ROUJEAN-K0-K1-K2-VL{CTE_TETAV_LIM_ROUJEAN}-
C          -SL{CTE_TETAS_LIM_ROUJEAN}-MUnb-NomFicUser avec  
C      - kx au format F5.3 (ex : 000.500) 
C      - CTE_TETAV_LIM_ROUJEAN et CTE_TETAS_LIM_ROUJEAN respectivement
C        l'angle zénithal de visée et solaire limites pour les calculs de BRDF
C        (en degrés).
C      - MUnb le nombre d'angles de Gauss (ex : MU24)
C      - NomFicUser le nom du fichier utilisateur (sans arborescence) (ex : MyAngles)
C
C      On convertit Kx en une chaine de caractères CKx(1:6) 
C      Pour cela on passe du Format F7.3 * 1000 => A6
C        Le premier caractère CKx(1) est le chiffre des centaines de Kx,
C        Le second caractère CKx(2) est le chiffre des dizaines de Kx,
C        Le troisième caractère CKx(3) est le chiffre des unites de Kx,
C           etc.
C
C
C Description des paramètres entrée / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C       FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C	K0  (double)	: Paramètre du modèle de Roujean
C	K1  (double)	: Paramètre du modèle de Roujean
C	K2  (double)	: Paramètre du modèle de Roujean
C	FICS (CHARACTER*CTE_LENFIC1) : (S) nom du fichier de Roujean (sans arborescence)	
C
C 
C Constantes utilisees:
C --------------------
C     CTE_TETAV_LIM_ROUJEAN : Angle zénithal de visée maximal pour les calculs
C                             de BRDF.
C     CTE_TETAS_LIM_ROUJEAN : Angle zénithal solaire maximal pour les calculs
C                             de BRDF.
C     CTE_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C     CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
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


      IMPLICIT NONE

C* Declaration des variables
C  -------------------------              
      DOUBLE PRECISION K0,K1,K2    ! (E) Paramètres du modèle de Roujean

      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! (E) chemin et nom du fichier des angles 
                                            ! "utilisateurs"					
      CHARACTER*CTE_LENFIC1 FICS   ! (S) Nom du fichier de ROUJEAN 	
          
      CHARACTER*6 CK0,CK1,CK2	   ! Chaines de caractères
      CHARACTER*4 CVLIM,CSLIM	   ! Chaines de caractères
      
      INTEGER*2 I
      INTEGER*4 NBMU_GAUSS  ! (E) Nombre d'angles de Gauss utilises 

	    
C --------- Conversion en chaines de caractères des paramètres d'entrée	
      WRITE(CK0,'(I6)') int(K0*1000)		
      WRITE(CK1,'(I6)') int(K1*1000)
      WRITE(CK2,'(I6)') int(K2*1000)
      WRITE(CVLIM,'(I2)') int(CTE_TETAV_LIM_ROUJEAN)
      WRITE(CSLIM,'(I2)') int(CTE_TETAS_LIM_ROUJEAN)
      
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
      
      FICS = 'ROUJEAN-'

C -------- Construction de la chaine de caractères pour K0	    
      FICS(9:11) = CK0(1:3)      ! Chiffres des unités
      FICS(12:12) = '.'	     ! Virgule
      FICS(13:15) = CK0(4:6)     ! Trois premières décimales   	    
      FICS(16:16) = '-'

C -------- Construction de la chaine de caractères pour K1	    
      FICS(17:19) = CK1(1:3)     ! Chiffres des unités
      FICS(20:20) = '.'	     ! Virgule
      FICS(21:23) = CK1(4:6)     ! Trois premières décimales   	    
      FICS(24:24) = '-'

C -------- Construction de la chaine de caractères pour K2	    
      FICS(25:27) = CK2(1:3)     ! Chiffres des unités
      FICS(28:28) = '.'	     ! Virgule
      FICS(29:31) = CK2(4:6)     ! Trois premières décimales   	    
      FICS(32:32) = '-'

C -------- Construction de la chaine de caractères pour la limite
C -------- sur l'angle de visée.
      FICS(33:34) = 'VL'	    
      FICS(35:36) = CVLIM(1:2)     ! Chiffre des unités
      FICS(37:37) = '-'	        
 
C -------- Construction de la chaine de caractères pour la limite
C -------- sur l'angle solaire.
      FICS(38:39) = 'SL'	    
      FICS(40:41) = CSLIM(1:2)     ! Chiffre des unités
      FICS(42:42) = '-'	        
                        
C -------- Completude de la fin du fichier avec MUnb-NomFicUser	                  
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)
      
      END !SOS_NOM_FICROUJ
      


      
C==============================================================================
C PROCEDURE: SOS_NOM_FICRH
C ==========
C      Cette procédure génère le nom du fichier de Rondeaux - Herman
C      en fonction de :
C	- l'indice de réfraction vegetation / air : IND
C       - Du modèle de BRDF associé pour la composante d'intensité :
C           modèle de ROUJEAN en fonction des paramètres K0,K1 et K2.
C
C      On impose d'avoir :
C         IND au format F5.3.
C	  kx au format F7.3.
C
C      Le nom du fichier de RH doit être sous la forme:
C
C      RH-IND-ROUJEAN-K0-K1-K2-VL{CTE_TETAV_LIM_ROUJEAN}-
C        -SL{CTE_TETAS_LIM_ROUJEAN}-MUnb-NomFicUser 
C      avec
C      - IND au format F5.3 (ex : 1.500) 
C      - kx au format F7.3 (ex : 000.500) 
C      - CTE_TETAV_LIM_ROUJEAN et CTE_TETAS_LIM_ROUJEAN respectivement
C        l'angle zénithal de visée et solaire limites pour les calculs de BRDF
C        (en degrés).
C      - MUnb le nombre d'angles de Gauss (ex : MU24)
C      - NomFicUser le nom du fichier utilisateur (sans arborescence) (ex : MyAngles)
C
C      On convertit IND en une chaine de caractères IND(1:4) 
C      Pour cela on passe du Format 5.3 * 1000 => A4
C           Le premier caractère CIND(1) est le chiffre unite de IND,
C           etc.
C      On convertit Kx en une chaine de caractères CKx(1:6) 
C      Pour cela on passe du Format F7.3 * 1000 => A6
C        Le premier caractère CKx(1) est le chiffre des centaines de Kx,
C        Le second caractère CKx(2) est le chiffre des dizaines de Kx,
C        Le troisième caractère CKx(3) est le chiffre des unites de Kx,
C           etc.
C
C
C Description des paramètres entrée / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C       FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C	IND (double)	: Indice de vegetation
C	K0  (double)	: Paramètre du modèle de Roujean
C	K1  (double)	: Paramètre du modèle de Roujean
C	K2  (double)	: Paramètre du modèle de Roujean
C	FICS (CHARACTER*CTE_LENFIC1) : (S) nom du fichier de RH (sans arborescence)	
C
C 
C Constantes utilisees:
C --------------------
C     CTE_TETAV_LIM_ROUJEAN : Angle zénithal de visée maximal pour les calculs
C                             de BRDF.
C     CTE_TETAS_LIM_ROUJEAN : Angle zénithal solaire maximal pour les calculs
C                             de BRDF.
C     CTE_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C     CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
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


      IMPLICIT NONE

C* Declaration des variables
C  -------------------------         
      DOUBLE PRECISION IND		! (E) Indice de vegetation.
      DOUBLE PRECISION K0,K1,K2 	! (E) Paramètres du modèle de Roujean.	
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! (E) chemin et nom du fichier des angles 
                                            ! "utilisateurs"			
      CHARACTER*CTE_LENFIC1 FICS	! (S) nom du fichier RH. 	

      CHARACTER*6 CK0,CK1,CK2	! Chaines de caractères
      CHARACTER*4 CVLIM,CSLIM	! Chaines de caractères          
      CHARACTER*4 CIND 		! Chaine de caractère de l'indice de vegetation.
      
      INTEGER*2 I
      INTEGER*4 NBMU_GAUSS  ! (E) Nombre d'angles de Gauss utilises 
      			      	 

	    
C --------- Conversion en chaines de caractères des paramètres d'entrée	   		
      WRITE(CIND,'(I4)') int(IND*1000)		
      WRITE(CK0,'(I6)') int(K0*1000)		
      WRITE(CK1,'(I6)') int(K1*1000)
      WRITE(CK2,'(I6)') int(K2*1000)
      WRITE(CVLIM,'(I2)') int(CTE_TETAV_LIM_ROUJEAN)
      WRITE(CSLIM,'(I2)') int(CTE_TETAS_LIM_ROUJEAN)
      
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
      
C -------- Construction de la chaine de caractères pour l'indice IND	    
      FICS(4:4) = CIND(1:1)     ! Chiffre des unités
      FICS(5:5) = '.'	    ! Virgule
      FICS(6:8) = CIND(2:4)     ! Trois premières décimales   	    

C -------- Ajout du nom du fichier de BRDF associe
      FICS(9:17) = '-ROUJEAN-'

C -------- Construction de la chaine de caractères pour K0	    
      FICS(18:20) = CK0(1:3)     ! Chiffres des unités
      FICS(21:21) = '.'	     ! Virgule
      FICS(22:24) = CK0(4:6)     ! Trois premières décimales   	    
      FICS(25:25) = '-'

C -------- Construction de la chaine de caractères pour K1	    
      FICS(26:28) = CK1(1:3)     ! Chiffres des unités
      FICS(29:29) = '.'	     ! Virgule
      FICS(30:32) = CK1(4:6)     ! Trois premières décimales   	    
      FICS(33:33) = '-'

C -------- Construction de la chaine de caractères pour K2	    
      FICS(34:36) = CK2(1:3)     ! Chiffres des unités
      FICS(37:37) = '.'	     ! Virgule
      FICS(38:40) = CK2(4:6)     ! Trois premières décimales   	    
      FICS(41:41) = '-'

C     -------- Construction de la chaine de caractères pour la limite
C     -------- sur l'angle de visée.
      FICS(42:43) = 'VL'	    
      FICS(44:45) = CVLIM(1:2)     ! Chiffre des unités
      FICS(46:46) = '-'	        
 
C     -------- Construction de la chaine de caractères pour la limite
C     -------- sur l'angle solaire.
      FICS(47:48) = 'SL'	    
      FICS(49:50) = CSLIM(1:2)     ! Chiffre des unités
      FICS(51:51) = '-'	 
                   
C -------- Completude de la fin du fichier avec MUnb-NomFicUser	                  
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)

      END !SOS_NOM_FICRH      
      


      
      
C==============================================================================
C PROCEDURE: SOS_NOM_FICBREON
C ==========
C      Cette procédure génère le nom du fichier de BREON en fonction de :
C	- l'indice de réfraction sol / air : IND
C       - Du modèle de BRDF associé pour la composante d'intensité :
C           modèle de ROUJEAN en fonction des paramètres K0,K1 et K2.
C
C      On impose d'avoir :
C         IND au format F5.3.
C	  kx au format F7.3.
C
C      Le nom du fichier de Breon doit être sous la forme:
C
C      BREON-IND-ROUJEAN-K0-K1-K2-VL{CTE_TETAV_LIM_ROUJEAN}-
C        -SL{CTE_TETAS_LIM_ROUJEAN}-MUnb-NomFicUser
C      avec
C      - IND au format F5.3 (ex : 1.500) 
C      - kx au format F7.3 (ex : 0.500) 
C      - CTE_TETAV_LIM_ROUJEAN et CTE_TETAS_LIM_ROUJEAN respectivement
C        l'angle zénithal de visée et solaire limites pour les calculs de BRDF
C        (en degrés).
C      - MUnb le nombre d'angles de Gauss (ex : MU24)
C      - NomFicUser le nom du fichier utilisateur (sans arborescence) (ex : MyAngles)
C
C      On convertit IND en une chaine de caractères IND(1:4) 
C      Pour cela on passe du Format 5.3 * 1000 => A4
C           Le premier caractère CIND(1) est le chiffre unite de IND,
C           etc.
C      On convertit Kx en une chaine de caractères CKx(1:6) 
C      Pour cela on passe du Format F7.3 * 1000 => A6
C        Le premier caractère CKx(1) est le chiffre des centaines de Kx,
C        Le second caractère CKx(2) est le chiffre des dizaines de Kx,
C        Le troisième caractère CKx(3) est le chiffre des unites de Kx,
C           etc.
C
C Description des paramètres entrée / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C       FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C	IND (double)	: Indice de sol
C	K0  (double)	: Paramètre du modèle de Roujean
C	K1  (double)	: Paramètre du modèle de Roujean
C	K2  (double)	: Paramètre du modèle de Roujean
C	FICS (CHARACTER*CTE_LENFIC1) : (S) nom du fichier de Breon (sans arborescence)	
C
C 
C Constantes utilisees:
C --------------------
C     CTE_TETAV_LIM_ROUJEAN : Angle zénithal de visée maximal pour les calculs
C                             de BRDF.
C     CTE_TETAS_LIM_ROUJEAN : Angle zénithal solaire maximal pour les calculs
C                             de BRDF.
C     CTE_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C     CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
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


      IMPLICIT NONE

C* Declaration des variables
C  -------------------------         
      DOUBLE PRECISION IND		! (E) Indice de sol.
      DOUBLE PRECISION K0,K1,K2 	! (E) Paramètres du modèle de Roujean.	
     
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! (E) chemin et nom du fichier des angles "utilisateurs"				
      CHARACTER*CTE_LENFIC1 FICS	! (S) nom du fichier Breon. 	

      CHARACTER*6 CK0,CK1,CK2	! Chaines de caractères
      CHARACTER*4 CVLIM,CSLIM	! Chaines de caractères          
      CHARACTER*4 CIND 		! Chaine de caractère de l'indice de vegetation.
      
      INTEGER*2 I
      INTEGER*4 NBMU_GAUSS  ! (E) Nombre d'angles de Gauss utilises 
	    
C --------- Conversion en chaines de caractères des paramètres d'entrée	   		
      WRITE(CIND,'(I4)') int(IND*1000)		
      WRITE(CK0,'(I6)') int(K0*1000)		
      WRITE(CK1,'(I6)') int(K1*1000)
      WRITE(CK2,'(I6)') int(K2*1000)
      WRITE(CVLIM,'(I2)') int(CTE_TETAV_LIM_ROUJEAN)
      WRITE(CSLIM,'(I2)') int(CTE_TETAS_LIM_ROUJEAN)
      
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
      
C -------- Construction de la chaine de caractères pour l'indice IND	    
      FICS(7:7) = CIND(1:1)      ! Chiffre des unités
      FICS(8:8) = '.'	     ! Virgule
      FICS(9:11) = CIND(2:4)     ! Trois premières décimales   	    

C -------- Ajout du nom du fichier de BRDF associe
      FICS(12:20) = '-ROUJEAN-'

C -------- Construction de la chaine de caractères pour K0	    
      FICS(21:23) = CK0(1:3)     ! Chiffres des unités
      FICS(24:24) = '.'	     ! Virgule
      FICS(25:27) = CK0(4:6)     ! Trois premières décimales   	    
      FICS(28:28) = '-'

C -------- Construction de la chaine de caractères pour K1	    
      FICS(29:31) = CK1(1:3)     ! Chiffres des unités
      FICS(32:32) = '.'	     ! Virgule
      FICS(33:35) = CK1(4:6)     ! Trois premières décimales   	    
      FICS(36:36) = '-'

C -------- Construction de la chaine de caractères pour K2	    
      FICS(37:39) = CK2(1:3)     ! Chiffres des unités
      FICS(40:40) = '.'	     ! Virgule
      FICS(41:43) = CK2(4:6)     ! Trois premières décimales   	    
      FICS(44:44) = '-'

C     -------- Construction de la chaine de caractères pour la limite
C     -------- sur l'angle de visée.
      FICS(45:46) = 'VL'	    
      FICS(47:48) = CVLIM(1:2)     ! Chiffre des unités
      FICS(49:49) = '-'	        
 
C     -------- Construction de la chaine de caractères pour la limite
C     -------- sur l'angle solaire.
      FICS(50:51) = 'SL'	    
      FICS(52:53) = CSLIM(1:2)     ! Chiffre des unités
      FICS(54:54) = '-'	 
                   
C -------- Completude de la fin du fichier avec MUnb-NomFicUser	                  
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)


      END !SOS_NOM_FICBREON     
      


 
 
      
C==============================================================================
C PROCEDURE: SOS_NOM_FICNADAL
C ==========
C      Cette procédure génère le nom du fichier du modèle de NADAL
C      en fonction de :
C	- la valeur de saturation alpha : ALPHA
C	- la valeur de beta : BETA
C       - l'indice de réfraction surface / atmosphere : IND
C
C      On impose d'avoir : 
C         IND au format F5.3.
C         ALPHA au format F6.4,
C	  BETA au format F4.1.
C
C      Le nom du fichier de Nadal doit être sous la forme:
C
C      NADAL-IND-K0-K1-K2-MUnb-NomFicUser avec  
C      - ALPHA au format F6.4 (ex : 0.0229) 
C      - BETA au format F3.1 (ex: 38)
C      - MUnb le nombre d'angles de Gauss (ex : MU24)
C      - NomFicUser le nom du fichier utilisateur (sans arborescence) (ex : MyAngles)
C
C      On convertit ALPHA en une chaine de caractères CALPHA(1:6) 
C      Pour cela on passe du Format F6.4 * 10000 => A6
C      Le premier caractère CALPHA(1) est le chiffre unite de ALPHA,
C           etc.
C
C      On convertit BETA en une chaine de caractères CBETA(1:4) 
C      Pour cela on passe du Format F4.1 * 10 => A4
C      Le premier caractère CBETA(1) est le chiffre des dizaines de BETA,
C           etc.
C
C Description des paramètres entrée / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C       FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C       IND (double)	: Indice de réfraction surface / atmosphere
C	ALPHA (double)	: Paramètre alpha du modelde de Nadal
C	BETA  (double)	: Paramètre beta du modelde de Nadal
C	K0  (double)	: Paramètre du modèle de Roujean
C	K1  (double)	: Paramètre du modèle de Roujean
C	K2  (double)	: Paramètre du modèle de Roujean
C	FICS (CHARACTER*CTE_LENFIC1) : (S) nom du fichier de Nadal (sans arborescence)	
C
C 
C Constantes utilisees:
C --------------------
C     CTE_TETAV_LIM_ROUJEAN : Angle zénithal de visée maximal pour les calculs
C                             de BRDF.
C     CTE_TETAS_LIM_ROUJEAN : Angle zénithal solaire maximal pour les calculs
C                             de BRDF.
C     CTE_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C     CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
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


      IMPLICIT NONE

C* Declaration des variables
C  -------------------------  
      DOUBLE PRECISION IND	 ! (E) Indice de sol       
      DOUBLE PRECISION ALPHA ! (E) Paramètre alpha du modelde de Nadal
      DOUBLE PRECISION BETA	 ! (E) Paramètre beta du modelde de Nadal
      DOUBLE PRECISION K0,K1,K2  ! (E) Paramètres du modèle de Roujean.				
      
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! (E) chemin et nom du fichier des angles "utilisateurs"					
      CHARACTER*CTE_LENFIC1 FICS ! (S) nom du fichier NADAL	
          
      CHARACTER*5 CALPHA	! Chaines de caractères
      CHARACTER*3 CBETA 
      CHARACTER*4 CIND 
      CHARACTER*6 CK0,CK1,CK2	
      CHARACTER*4 CVLIM,CSLIM	
             
      INTEGER*2 I
      INTEGER*4 NBMU_GAUSS  ! (E) Nombre d'angles de Gauss utilises 
      


C --------- Conversion en chaines de caractères des paramètres d'entrée
      WRITE(CIND,'(I4)') int(IND*1000)				
      WRITE(CALPHA,'(I5)') int(ALPHA*10000)
      WRITE(CBETA,'(I3)') int(BETA*10)	
      WRITE(CK0,'(I6)') int(K0*1000)		
      WRITE(CK1,'(I6)') int(K1*1000)
      WRITE(CK2,'(I6)') int(K2*1000)
      WRITE(CVLIM,'(I2)') int(CTE_TETAV_LIM_ROUJEAN)
      WRITE(CSLIM,'(I2)') int(CTE_TETAS_LIM_ROUJEAN)
            
    
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
  
C -------- Construction de la chaine de caractères pour l'indice IND	    
      FICS(7:7) = CIND(1:1)       ! Chiffre des unités
      FICS(8:8) = '.'	      ! Virgule
      FICS(9:11) = CIND(2:4)      ! Trois premières décimales   	    
      FICS(12:12) = '-'
            
C -------- Construction de la chaine de caractères pour alpha	    
      FICS(13:13) = CALPHA(1:1)     ! Chiffre des unites
      FICS(14:14) = '.'	        ! Virgule
      FICS(15:18) = CALPHA(2:5)     ! Quatre premières décimales   	    
      FICS(19:19) = '-'

C -------- Construction de la chaine de caractères pour beta	    
      FICS(20:21) = CBETA(1:2)      ! Chiffre des unites 
      FICS(22:22) = '.'
      FICS(23:23) = CBETA(3:3)      ! Chiffre des unites 

C -------- Ajout du nom du fichier de BRDF associe
      FICS(24:32) = '-ROUJEAN-'

C -------- Construction de la chaine de caractères pour K0	    
      FICS(33:35) = CK0(1:3)     ! Chiffres des unités
      FICS(36:36) = '.'	     ! Virgule
      FICS(37:39) = CK0(4:6)     ! Trois premières décimales   	    
      FICS(40:40) = '-'
		
C -------- Construction de la chaine de caractères pour K1	    
      FICS(41:43) = CK1(1:3)     ! Chiffres des unités
      FICS(44:44) = '.'	     ! Virgule
      FICS(45:47) = CK1(4:6)     ! Trois premières décimales   	    
      FICS(48:48) = '-'
		
C -------- Construction de la chaine de caractères pour K2	    
      FICS(49:51) = CK2(1:3)     ! Chiffres des unités
      FICS(52:52) = '.'	     ! Virgule
      FICS(53:55) = CK2(4:6)     ! Trois premières décimales   	    
      FICS(56:56) = '-'
		
C     -------- Construction de la chaine de caractères pour la limite
C     -------- sur l'angle de visée.
      FICS(57:58) = 'VL'	    
      FICS(59:60) = CVLIM(1:2)     ! Chiffre des unités
      FICS(61:61) = '-'	        
 
C     -------- Construction de la chaine de caractères pour la limite
C     -------- sur l'angle solaire.
      FICS(62:63) = 'SL'	    
      FICS(64:65) = CSLIM(1:2)     ! Chiffre des unités
      FICS(66:66) = '-'	 
      
    			                 
C -------- Completude de la fin du fichier avec MUnb-NomFicUser	                  
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)
   
      END !SOS_NOM_FICNADAL      


 

      
C==============================================================================
C PROCEDURE: SOS_NOM_FICMAIGNAN
C ==========
C      Cette procédure génère le nom du fichier de MAIGNAN en fonction de :
C	- L'indice de réfraction sol / air : IND
C      - Le coefficient C du modèle de Maignan
C      - Du modèle de BRDF associé pour la composante d'intensité :
C         modèle de ROUJEAN en fonction des paramètres K0,K1 et K2.
C
C      On impose d'avoir :
C        IND au format F5.3.
C        C au format F7.3
C	  kx au format F7.3.
C
C      Le nom du fichier de Maignan doit être sous la forme:
C
C      MAIGNAN-IND-C-ROUJEAN-K0-K1-K2-VL{CTE_TETAV_LIM_ROUJEAN}-
C        -SL{CTE_TETAS_LIM_ROUJEAN}-MUnb-NomFicUser
C      avec
C      - IND au format F5.3 (ex : 1.500) 
C      - C au format F7.3 (ex : 007.290) 
C      - kx au format F7.3 (ex : 000.500) 
C      - CTE_TETAV_LIM_ROUJEAN et CTE_TETAS_LIM_ROUJEAN respectivement
C        l'angle zénithal de visée et solaire limites pour les calculs de BRDF
C        (en degrés).
C      - MUnb le nombre d'angles de Gauss (ex : MU24)
C      - NomFicUser le nom du fichier utilisateur (sans arborescence) (ex : MyAngles)
C
C      On convertit IND en une chaine de caractères IND(1:4) 
C      Pour cela on passe du Format 5.3 * 1000 => A4
C           Le premier caractère CIND(1) est le chiffre unite de IND,
C           etc.
C      On convertit Kx en une chaine de caractères CKx(1:6) 
C      Pour cela on passe du Format F7.3 * 1000 => A6
C        Le premier caractère CKx(1) est le chiffre des centaines de Kx,
C        Le second caractère CKx(2) est le chiffre des dizaines de Kx,
C        Le troisième caractère CKx(3) est le chiffre des unites de Kx,
C           etc.
C      De même pour le coefficient C.
C
C Description des paramètres entrée / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C      FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C	IND (double)	  : (E) Indice de réfraction surface / air
C      COEF_C (double) : (E) Paramètre du modèle de Maignan
C	K0  (double)	: (E) Paramètre du modèle de Roujean
C	K1  (double)	: (E) Paramètre du modèle de Roujean
C	K2  (double)	: (E) Paramètre du modèle de Roujean
C	FICS (CHARACTER*CTE_LENFIC1) : (S) nom du fichier de Maignan (sans arborescence)	
C
C 
C Constantes utilisees:
C --------------------
C     CTE_TETAV_LIM_ROUJEAN : Angle zénithal de visée maximal pour les calculs
C                             de BRDF.
C     CTE_TETAS_LIM_ROUJEAN : Angle zénithal solaire maximal pour les calculs
C                             de BRDF.
C     CTE_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C     CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
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
      SUBROUTINE SOS_NOM_FICMAIGNAN(NBMU_GAUSS,FICANGLES_USER,IND,
     &                              COEF_C,K0,K1,K2,FICS) 


      IMPLICIT NONE

C* Declaration des variables
C  -------------------------         
      DOUBLE PRECISION IND		! (E) Indice de sol.
      DOUBLE PRECISION COEF_C      ! (E) Paramètre C*exp(-NDVI) du modèle de Maignan.
      DOUBLE PRECISION K0,K1,K2 	! (E) Paramètres du modèle de Roujean.	
     
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! (E) chemin et nom du fichier des angles "utilisateurs"				
      CHARACTER*CTE_LENFIC1 FICS	! (S) nom du fichier Breon. 	

      CHARACTER*6 CK0,CK1,CK2	! Chaines de caractères
      CHARACTER*4 CVLIM,CSLIM	! Chaines de caractères          
      CHARACTER*4 CIND 		! Chaine de caractère de l'indice de vegetation.
      CHARACTER*6 CC               ! Chaines de caractères pour coefficient C
      
      INTEGER*2 I
      INTEGER*4 NBMU_GAUSS  ! (E) Nombre d'angles de Gauss utilises 
	    
C --------- Conversion en chaines de caractères des paramètres d'entrée	   		
      WRITE(CIND,'(I4)') int(IND*1000)		
      WRITE(CC,'(I6)')  int(COEF_C*1000)		
      WRITE(CK0,'(I6)') int(K0*1000)		
      WRITE(CK1,'(I6)') int(K1*1000)
      WRITE(CK2,'(I6)') int(K2*1000)
      WRITE(CVLIM,'(I2)') int(CTE_TETAV_LIM_ROUJEAN)
      WRITE(CSLIM,'(I2)') int(CTE_TETAS_LIM_ROUJEAN)
      
C --------- Les blancs sont remplaces par des 0	   
      DO I=1,4
	 IF( CIND(I:I).EQ.' ') CIND(I:I) = '0'
      ENDDO      
      DO I=1,6
         IF( CC(I:I).EQ.' ') CC(I:I) = '0'
         IF( CK0(I:I).EQ.' ') CK0(I:I) = '0'
         IF( CK1(I:I).EQ.' ') CK1(I:I) = '0'   
         IF( CK2(I:I).EQ.' ') CK2(I:I) = '0'
      ENDDO   
      DO I=1,2
         IF( CVLIM(I:I).EQ.' ') CVLIM(I:I) = '0'
         IF( CSLIM(I:I).EQ.' ') CSLIM(I:I) = '0'
      ENDDO 
  
      
      FICS = 'MAIGNAN-'
      
C -------- Construction de la chaine de caractères pour l'indice IND	    
      FICS(9:9) = CIND(1:1)       ! Chiffre des unités
      FICS(10:10) = '.'	      ! Virgule
      FICS(11:13) = CIND(2:4)     ! Trois premières décimales   	    

C -------- Construction de la chaine de caractères pour le coefficient C 
      FICS(14:15) = '-C'
      FICS(16:18) = CC(1:3)      ! Chiffres des unités
      FICS(19:19) = '.'	     ! Virgule
      FICS(20:22) = CC(4:6)      ! Trois premières décimales   


C -------- Ajout du nom du fichier de BRDF associe
      FICS(23:31) = '-ROUJEAN-'

C -------- Construction de la chaine de caractères pour K0	    
      FICS(32:34) = CK0(1:3)     ! Chiffres des unités
      FICS(35:35) = '.'	     ! Virgule
      FICS(36:38) = CK0(4:6)     ! Trois premières décimales   	    
      FICS(39:39) = '-'

C -------- Construction de la chaine de caractères pour K1	    
      FICS(40:42) = CK1(1:3)     ! Chiffres des unités
      FICS(43:43) = '.'	     ! Virgule
      FICS(44:46) = CK1(4:6)     ! Trois premières décimales   	    
      FICS(47:47) = '-'

C -------- Construction de la chaine de caractères pour K2	    
      FICS(48:50) = CK2(1:3)     ! Chiffres des unités
      FICS(51:51) = '.'	     ! Virgule
      FICS(52:54) = CK2(4:6)     ! Trois premières décimales   	    
      FICS(55:55) = '-'

C     -------- Construction de la chaine de caractères pour la limite
C     -------- sur l'angle de visée.
      FICS(56:57) = 'VL'	    
      FICS(58:59) = CVLIM(1:2)     ! Chiffre des unités
      FICS(60:60) = '-'	        
 
C     -------- Construction de la chaine de caractères pour la limite
C     -------- sur l'angle solaire.
      FICS(61:62) = 'SL'	    
      FICS(63:64) = CSLIM(1:2)     ! Chiffre des unités
      FICS(65:65) = '-'	 
                   
C -------- Completude de la fin du fichier avec MUnb-NomFicUser	                  
      CALL SOS_COMPLETE_NOM_FIC_SURF(NBMU_GAUSS,FICANGLES_USER,FICS)


      END !SOS_NOM_FICMAIGNAN     
      



      
C==============================================================================
C PROCEDURE: SOS_COMPLETE_NOM_FIC_SURF
C ==========

C
C Description des paramètres entrée / sortie
C ------------------------------------------
C	NBMU_GAUSS (I4)	: (E) Nombre d'angles de Gauss utilises
C       FICANGLES_USER (CHARACTER*CTE_LENFIC2)	: (E) Nom du fichier des angles 
C                                                 "utilisateurs" (chemin complet)
C	FICS (CHARACTER*CTE_LENFIC1) : (E/S) nom du fichier de BRDF/BPDF (sans arborescence)
C 
C Constantes utilisees:
C --------------------
C    CTE_LENFIC1 : Longueur des noms de fichiers sans arborescence.
C    CTE_LENFIC2 : Longueur des noms de fichiers avec arborescence.
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

      IMPLICIT NONE

C* Declaration des variables
C  -------------------------  
      CHARACTER*CTE_LENFIC2 FICANGLES_USER  ! (E) chemin et nom du fichier des angles "utilisateurs"	
      CHARACTER*CTE_LENFIC1 FICS            ! (E/S) nom du fichier GLITTER 	

      CHARACTER*4 CNBMU     ! Chaine de caractères associee a NBMU_GAUSS
      
      INTEGER*2 I
      
      INTEGER*4 NBMU_GAUSS  ! (E) Nombre d'angles de Gauss utilises 
      INTEGER*4 IDEB, IDEB_CH, IFIN, IFIN_CH  ! position de debut et fin d'une chaine de caractères  
      
      			      	 

	    
C --------- Conversion en chaines de caractères de NBMU_GAUSS   
      WRITE(CNBMU,'(I4)') NBMU_GAUSS      
      

C -------- Construction de la chaine de caractères pour MUgauss
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
  
       
C -------- Construction de la chaine de caractères pour NomFicUser
C -------- Structuration pour avoir NBMU_GAUSS a gauche (sans "blancs")	            


      IF (FICANGLES_USER.NE.'NO_USER_ANGLES') THEN

C       Determination de la position du nom de fichier sans son arborescence          
	  IFIN_CH=INDEX(FICANGLES_USER,' ')
          IFIN_CH=IFIN_CH-1
	  IF (IFIN_CH.LE.0) IFIN_CH=CTE_LENDIR
	  
	  IDEB_CH=0
	  	  
10	  IDEB=INDEX(FICANGLES_USER(IDEB_CH:IFIN_CH),'/')
	  IDEB_CH=IDEB_CH+IDEB
          IF (IDEB.NE.0) GOTO 10	  

C       Ajout du nom de fichier d'angles "utilisateur" au nom du fichier     
          IFIN=INDEX(FICS,' ')  
	  FICS(IFIN:IFIN) = '-'
	  
	  IFIN=IFIN+1
	  FICS(IFIN:IFIN+IFIN_CH-IDEB_CH+1) 
     &	  = FICANGLES_USER(IDEB_CH:IFIN_CH)
	  
      ENDIF 
   
      END !SOS_COMPLETE_NOM_FIC_SURF   
