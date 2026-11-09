# %% [markdown]
# # Campagnes thermohydrauliques turbulentes — contrats de validation
#
# Cette fiche définit quatre campagnes physiques. Les statuts initiaux sont
# PLANNED/BLOCKED : les données de référence et la théorie ne constituent pas
# une exécution CFDX ni une qualification.
#
# ## Matrice synthétique
#
# | ID | Cas | Référence | Modèle de départ | Quantités principales |
# |---|---|---|---|---|
# | TUR-CHANNEL | Canal plan DNS | Moser–Kim–Mansour (1999) | RANS SST puis RSM; DNS/LES si disponible | U+, contraintes, Cf, Re_b |
# | SEP-BFS-DRIVER | Marche descendante | Driver–Seegmiller (1985) | SST, realizable k-epsilon | x_r/h, Cp, Cf, U/V, R_ij |
# | TH-MIXED-CONVECTION-PIPE | Conduite verticale chauffée | Expérience upflow (1998) | SST avec flottabilité | U, T, Nu, Cf, Ri |
# | TH-TJUNCTION | Mélange en T | OECD/NEA–Vattenfall | SST/URANS puis LES | profils U/T, RMS T, spectres |
#
# ## 1. TUR-CHANNEL — Moser, Kim & Mansour
#
# ### Définition physique
#
# Canal plan incompressible pleinement développé, demi-hauteur \\(h\\),
# périodicité longitudinale et spanwise, parois fixes sans glissement.
# La DNS fournit trois nombres de Reynolds frictionnels :
# \\(Re_\\tau=u_\\tau h/\\nu=180,395,590\\).
# Le domaine de référence est périodique selon x et z ; les parois sont en
# y=±h. La base de données UT Austin contient profils moyens, contraintes de
# Reynolds, skewness/flatness, bilans et corrélations.
#
# ### Données et provenance à geler
#
# - Source primaire : https://turbulence.oden.utexas.edu/MKM_1999.html
# - Index miroir : https://torroja.dmt.upm.es/turbdata/agard/chapter5/PCH10/
# - Télécharger séparément chan180, chan395 et chan590.
# - Archiver URL, date de récupération, licence/conditions, nom exact du fichier
#   et SHA-256. Ne jamais reconstituer les profils depuis une figure.
#
# ### Grandeurs de comparaison
#
# \\(U^+=U/u_\\tau\\), \\(y^+=yu_\\tau/\\nu\\),
# \\(u_{rms}^+,v_{rms}^+,w_{rms}^+\\),
# \\(-\\overline{u'v'}^+=-\\overline{u'v'}/u_\\tau^2\\),
# \\(k^+=k/u_\\tau^2\\), \\(C_f\\), \\(Re_b\\) et gradient de pression.
# Vérifier aussi la contrainte totale, somme contrainte visqueuse et contrainte
# turbulente, qui doit suivre le profil linéaire imposé par l'équilibre.
#
# ### Limite des fermetures
#
# Le SST à viscosité turbulente peut être comparé à U+ et à la contrainte de
# cisaillement, mais ne constitue pas une prédiction libre des six composantes
# de R_ij. Les contraintes normales DNS sont à traiter comme diagnostic
# d'anisotropie ; une revendication sur ces composantes nécessite un RSM.
# Un RANS stationnaire ne reproduit pas les spectres DNS.
#
# ### Plan de calcul
#
# - Lot 1 : Re_tau=180, vérifier la convention de pression/force motrice,
#   le débit et l'équilibre de contrainte.
# - Lot 2 : répéter à 395 puis 590, sans recalibrer les seuils après résultat.
# - Au moins trois résolutions en y, et résolution x/z cohérente avec le modèle.
# - Pour RANS, rapporter y+ effectif et traitement de paroi ; pour LES/DNS,
#   démontrer la résolution des échelles et la convergence statistique.
#
# ## 2. SEP-BFS-DRIVER — Driver & Seegmiller
#
# ### Référence et géométrie
#
# Driver, D. M. & Seegmiller, H. L. (1985), “Features of a Reattaching
# Turbulent Shear Layer in Divergent Channel Flow”, AIAA Journal 23, 163–171.
# Utiliser le cas à paroi opposée droite (straight-wall), pas le cas divergent.
# La marche vaut h=12.7 mm ; le rapport d'expansion est 1.125 ; la hauteur
# amont est 8h et la largeur rapportée est 12h. Le Reynolds de référence est
# Re_h=U_ref h/nu≈37,500, avec U_ref≈44.2 m/s.
#
# La couche limite amont est mesurée : à 4h en amont de la marche,
# δ≈19 mm, δ*≈2.7 mm, θ≈1.9 mm, shape factor H≈1.42 et Cf≈0.0029.
# Utiliser le profil expérimental d'entrée et les stations LDV du jeu de
# données. Ne pas substituer un profil uniforme.
#
# ### Quantités et extraction
#
# - Position de recollement x_r/h définie par le changement de signe de τ_w.
# - Cp(x), Cf(x), profils U/U_ref et V/U_ref aux stations expérimentales.
# - Contraintes u'u', v'v', -u'v' aux stations correspondantes.
# - Bilan massique et sens du cisaillement pariétal dans la zone séparée.
#
# La longueur de recollement ne doit pas être estimée visuellement sur une
# image de vitesse. La définition expérimentale et la méthode d'interpolation
# du zéro de τ_w doivent être archivées.
#
# ### Modèles et distinction des cas
#
# SST est le premier candidat RANS pour le gradient de pression adverse et
# la séparation ; realizable k-epsilon est une comparaison utile. Tester
# sensibilité au y+, à la résolution du cisaillement détaché et au profil
# d'entrée. Le cas SEP-BFS déjà enregistré pour Armaly/VMFL064 (Re_D=200)
# reste un autre benchmark : ne pas réutiliser ses paramètres ou sa référence.
#
# IMPORTANT — Le fichier TH-07.json existant est la DNS Wagner-Wetzel 2022
# d'un CANAL PLAN VERTICAL (parois differentiellement chauffees). Ce n'est
# PAS l'experience 1998 de conduite CIRCULAIRE ascendante ci-dessous.
# L'experience de conduite possede maintenant TH-PIPE-UPFLOW-1998.json,
# une fiche de preparation sans statut de simulation ni valeurs inventees.
# Ne jamais la verifier a partir des sorties ou des references TH-07.
#
# ## 3. TH-MIXED-CONVECTION-PIPE — conduite verticale chauffée
#
# ### Référence expérimentale
#
# “Upflow turbulent mixed convection heat transfer in vertical pipes”,
# International Journal of Heat and Mass Transfer 41 (1998), 4037–4054,
# DOI 10.1016/S0017-9310(98)00177-X. L'expérience porte sur l'eau en
# écoulement ascendant dans une conduite verticale chauffée, du régime forcé
# vers le régime de convection mixte, et documente la réduction du transfert
# thermique associée à l'effet de flottabilité et à la laminarisation proche
# paroi dans certains régimes.
#
# Les valeurs précises de géométrie, propriétés, flux, Re et Gr doivent être
# transcrites depuis les tables de l'article avant de geler les fichiers
# d'entrée. Ne pas inventer de valeurs intermédiaires ni utiliser les plages
# d'une autre publication comme si elles appartenaient à cette expérience.
#
# ### Théorie et balayage
#
# Utiliser Re=Ub D/nu, Pr=nu/alpha, Gr=g beta ΔT D^3/nu^2 et
# Ri=Gr/Re^2, avec ΔT et propriétés de référence explicitement définies.
# Conserver la convention de signe de l'écoulement ascendant et de la gravité.
# Construire trois familles : forced-convection baseline, buoyancy-aided et
# buoyancy-opposed lorsque la référence le permet. Pour chaque point, documenter
# le paramètre de flottabilité effectivement couvert par l'expérience.
#
# ### Quantités
#
# Profils U(r), T(r), k et contraintes disponibles ; τ_w, Cf, Nu local et
# moyen, température de mélange, perte de charge et bilan énergétique.
# Examiner l'évolution du profil de vitesse et du frottement avec Ri.
# La comparaison Boussinesq n'est admissible que si l'hypothèse de faibles
# variations de densité est vérifiée ; sinon employer les propriétés variables
# de l'expérience.
#
# ## 4. TH-TJUNCTION — OECD/NEA–Vattenfall
#
# ### Configuration de référence
#
# Utiliser le rapport officiel OECD/NEA, NEA/CSNI/R(2011)5 et l'article
# Smith, Mahaffy & Angele (2013), Nuclear Engineering and Design 264, 80–88,
# DOI 10.1016/j.nucengdes.2013.02.030. Il s'agit d'une expérience à l'échelle
# 1:2, réalisée à Vattenfall R&D, conçue pour l'étude du mélange turbulent et
# de la fatigue thermique à grand nombre de cycles.
#
# Le dossier public décrit de l'eau froide à 19 °C dans la conduite principale
# et de l'eau chaude à 36 °C dans la branche, soit ΔT=17 K pour la campagne
# thermique. Les débits de référence rapportés sont 9 L/s et 6 L/s,
# correspondant approximativement à Re=80,000 et 100,000. La conduite froide
# possède une longue section amont (environ 80D) ; la branche chaude dispose
# d'environ 20 diamètres hydrauliques. Les dimensions exactes, repères de
# diamètres, coordonnées des thermocouples et incertitudes doivent être copiés
# depuis le rapport officiel et non déduits d'un schéma.
#
# ### Deux campagnes obligatoirement séparées
#
# A. Isotherme : deux entrées à même température, comparaison LDV des profils
# de vitesse dans les branches et aux stations aval. Cela isole la prédiction
# hydrodynamique et le partage de débit.
#
# B. Non-isotherme : entrées à 19 °C et 36 °C ; comparaison des températures
# moyennes et des signaux thermocouples aux positions publiées. La température
# affecte densité/viscosité selon la plage et le modèle retenu.
#
# ### Quantités et exigences temporelles
#
# Comparer profils U/V, température moyenne, RMS de température, amplitude,
# fréquence dominante et densité spectrale lorsque les signaux bruts et leur
# fréquence d'échantillonnage sont disponibles. Les stations aval rapportées
# comprennent x/D=2.6 et 6.6 pour les mesures LDV. Le placement exact de chaque
# thermocouple doit venir du rapport.
#
# Un RANS stationnaire peut être évalué sur les moyennes et les bilans, mais
# ne peut pas valider les fluctuations thermiques, RMS, spectres ni les
# extrema temporels. Ces objectifs nécessitent un calcul instationnaire
# (URANS ou LES) et une durée d'échantillonnage suffisante.
#
# ## 5. Critères transversaux et statut
#
# Pour chaque modèle/maillage : résidus, stabilité des QoI, bilan massique,
# bilan de quantité de mouvement si exploitable, bilan énergétique, temps
# physique et statistique, y+, schémas, sous-relaxations, propriétés, version
# des données et SHA du cas doivent être conservés.
#
# Utiliser au moins trois maillages et, pour les cas instationnaires, trois
# pas de temps lorsque cela est réalisable. Calculer les erreurs avec les
# incertitudes de référence ; définir les seuils avant l'exécution.
#
# Statut initial des quatre campagnes : PLANNED/BLOCKED. La disponibilité
# d'une référence publique n'implique pas que le chemin de production CFDX
# sache encore exécuter le modèle physique demandé.
#
# ## Références et données
#
# - Moser–Kim–Mansour DNS : https://turbulence.oden.utexas.edu/MKM_1999.html
# - Driver–Seegmiller (1985) : AIAA Journal 23, 163–171 ; fiche de données
#   AGARD-AR-345, cas CMP32.
# - OECD/NEA Vattenfall : https://www.oecd-nea.org/upload/docs/application/pdf/2021-03/csni-r2011-5.pdf
# - T-junction benchmarking paper : DOI 10.1016/j.nucengdes.2013.02.030
# - Upflow mixed convection (1998) : DOI 10.1016/S0017-9310(98)00177-X
# - High-fidelity mixed-convection pipe dataset (2026), DOI 10.57745/6CGPI4,
#   Recherche Data Gouv; use as an additional numerical dataset, not as
#   experimental evidence.
