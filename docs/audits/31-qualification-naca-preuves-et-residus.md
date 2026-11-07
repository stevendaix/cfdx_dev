# FICHE MÉTHODE 31 — Intégrité des preuves NACA et qualification sur trois grilles

**En-tête d'évidence**
- Scope : lecture des historiques physiques, certificats de forces et résumés de raffinement N64/N128/N256.
- Statut : vérification logicielle ciblée ; validation aérodynamique et GCI qualifié non obtenus.
- Base : archive locale CFDX du 2026-11-07 ; aucune équivalence supposée avec les PR parallèles.
- Tests : contre-exemples négatifs, comparaison de traînée pression/viscosité et non-régression Python.

## 1. Références sources
- NASA NPARC, *Examining Spatial (Grid) Convergence*, https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html ; approche Roache/indices de convergence de grille.
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer 2016, DOI 10.1007/978-3-319-16874-6.
- NASA NPARC, *Verification Assessment*, https://www.grc.nasa.gov/www/wind/valid/tutorial/verassess.html.

## 2. Formulation continue
Le débit de masse normal à une face est $\rho\mathbf{U}\cdot\mathbf{S}$. À équilibre, le flux net d'un volume fermé doit être nul. Les résidus de continuité et de momentum sont des **normes**, donc ne peuvent jamais être négatifs. Pour une même référence de normalisation aérodynamique, la traînée totale est la somme de ses composantes pression et visqueuse : $C_d=C_{d,p}+C_{d,v}$.

## 3. Discrétisation Volume Fini
Les historiques doivent représenter le résidu physique réévalué à l'état retourné, pas uniquement le résidu algébrique préconditionné. Une condition nécessaire de convergence est $0\leq R/R_0\leq\varepsilon$ pour chaque résidu physique requis, avec une définition non dégénérée de l'échelle $R_0$. Les certificats doivent satisfaire la fermeture des efforts dans la précision de leur représentation numérique.

Sur trois grilles géométriquement comparables, un ordre observé et un GCI peuvent être estimés **seulement** après convergence itérative de tous les niveaux. Le raffinement de l'épaisseur de première couche et des dimensions radiales doit être contrôlé indépendamment du nombre de faces sur le profil. Un calcul des différences de traînée seul ne suffit pas.

## 4. Traitement des frontières
- Les intégrales pression/viscosité sont calculées sur les faces physiques du profil, avec des normales et une orientation cohérentes.
- Les frontières farfield et de symétrie ne doivent pas contaminer la force de paroi.
- Un fichier tronqué ou non fini, un nombre d'itérations nul/fractionnaire, une norme négative, une fermeture des forces fausse ou un certificat incompatible avec l'historique provoquent un refus explicite.

## 5. Pièges connus
- `-1e-12 < 1e-8` ne démontre pas une convergence : une norme négative traduit une preuve corrompue, même si elle Â« passe Â» une comparaison naïve.
- Le triplet `(N64,N128,N256)` ne garantit ni des familles géométriquement semblables ni un ordre observé asymptotique.
- `naca_three_grid_gci.py` est **DIAGNOSTIC_ONLY** : il ne lit pas les historiques physiques ni le rapport de qualité des trois maillages. La voie de certification forte est `naca_force_mesh_qualification.py`, qui contrôle aussi les certificats d'efforts et la similarité des grilles.
- Un marqueur `PASS` isolé ou des champs JSON cohérents entre eux ne sont pas une validation expérimentale.

## 6. Ordre de précision attendu et validation exigée
- Ces contrôles ne changent pas le schéma FV, son ordre ou les seuils physiques. La tolérance sur la somme des coefficients de traînée couvre uniquement les erreurs d'arrondi de leurs représentations en sortie ; elle n'est pas un critère de convergence.
- Tests red-green : normes négatives dans certificat/historique, continuité négative, nombre d'itérations invalide et incohérence de décomposition de traînée.
- Les campagnes physiques complètes, la convergence itérative stricte N64/N128/N256, les profils $C_p/C_f$, le GCI par quantité, les conditions de sortie et l'incertitude expérimentale doivent être prouvés indépendamment avant toute promotion.

## 7. Statut d'implémentation
- Lot local 2026-11-07 : gardes de qualification renforcés dans cinq scripts NACA, régressions négatives ajoutées.
- Tests Python NACA élargis : 109 pass, 3 skipped ; les campagnes longues N128/N256 et la référence expérimentale ne sont pas exécutées ici.
- Statut : **partiel / non qualifié physiquement**. Aucun coefficient de traînée ou de cisaillement n'est promu sur ce seul lot.
