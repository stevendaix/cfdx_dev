# FICHE MÃTHODE 31 â IntÃ©gritÃ© des preuves NACA et qualification sur trois grilles

**En-tÃªte d'Ã©vidence**
- Scope : lecture des historiques physiques, certificats de forces et rÃ©sumÃ©s de raffinement N64/N128/N256.
- Statut : vÃ©rification logicielle ciblÃ©e ; validation aÃ©rodynamique et GCI qualifiÃ© non obtenus.
- Base : archive locale CFDX du 2026-11-07 ; aucune Ã©quivalence supposÃ©e avec les PR parallÃ¨les.
- Tests : contre-exemples nÃ©gatifs, comparaison de traÃ®nÃ©e pression/viscositÃ© et non-rÃ©gression Python.

## 1. RÃ©fÃ©rences sources
- NASA NPARC, *Examining Spatial (Grid) Convergence*, https://www.grc.nasa.gov/www/wind/valid/tutorial/spatconv.html ; approche Roache/indices de convergence de grille.
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer 2016, DOI 10.1007/978-3-319-16874-6.
- NASA NPARC, *Verification Assessment*, https://www.grc.nasa.gov/www/wind/valid/tutorial/verassess.html.

## 2. Formulation continue
Le dÃ©bit de masse normal Ã  une face est $\rho\mathbf{U}\cdot\mathbf{S}$. Ã Ã©quilibre, le flux net d'un volume fermÃ© doit Ãªtre nul. Les rÃ©sidus de continuitÃ© et de momentum sont des **normes**, donc ne peuvent jamais Ãªtre nÃ©gatifs. Pour une mÃªme rÃ©fÃ©rence de normalisation aÃ©rodynamique, la traÃ®nÃ©e totale est la somme de ses composantes pression et visqueuse : $C_d=C_{d,p}+C_{d,v}$.

## 3. DiscrÃ©tisation Volume Fini
Les historiques doivent reprÃ©senter le rÃ©sidu physique rÃ©Ã©valuÃ© Ã  l'Ã©tat retournÃ©, pas uniquement le rÃ©sidu algÃ©brique prÃ©conditionnÃ©. Une condition nÃ©cessaire de convergence est $0\leq R/R_0\leq\varepsilon$ pour chaque rÃ©sidu physique requis, avec une dÃ©finition non dÃ©gÃ©nÃ©rÃ©e de l'Ã©chelle $R_0$. Les certificats doivent satisfaire la fermeture des efforts dans la prÃ©cision de leur reprÃ©sentation numÃ©rique.

Sur trois grilles gÃ©omÃ©triquement comparables, un ordre observÃ© et un GCI peuvent Ãªtre estimÃ©s **seulement** aprÃ¨s convergence itÃ©rative de tous les niveaux. Le raffinement de l'Ã©paisseur de premiÃ¨re couche et des dimensions radiales doit Ãªtre contrÃ´lÃ© indÃ©pendamment du nombre de faces sur le profil. Un calcul des diffÃ©rences de traÃ®nÃ©e seul ne suffit pas.

## 4. Traitement des frontiÃ¨res
- Les intÃ©grales pression/viscositÃ© sont calculÃ©es sur les faces physiques du profil, avec des normales et une orientation cohÃ©rentes.
- Les frontiÃ¨res farfield et de symÃ©trie ne doivent pas contaminer la force de paroi.
- Un fichier tronquÃ© ou non fini, un nombre d'itÃ©rations nul/fractionnaire, une norme nÃ©gative, une fermeture des forces fausse ou un certificat incompatible avec l'historique provoquent un refus explicite.

## 5. PiÃ¨ges connus
- `-1e-12 < 1e-8` ne dÃ©montre pas une convergence : une norme nÃ©gative traduit une preuve corrompue, mÃªme si elle Â« passe Â» une comparaison naÃ¯ve.
- Le triplet `(N64,N128,N256)` ne garantit ni des familles gÃ©omÃ©triquement semblables ni un ordre observÃ© asymptotique.
- `naca_three_grid_gci.py` est **DIAGNOSTIC_ONLY** : il ne lit pas les historiques physiques ni le rapport de qualitÃ© des trois maillages. La voie de certification forte est `naca_force_mesh_qualification.py`, qui contrÃ´le aussi les certificats d'efforts et la similaritÃ© des grilles.
- Un marqueur `PASS` isolÃ© ou des champs JSON cohÃ©rents entre eux ne sont pas une validation expÃ©rimentale.

## 6. Ordre de prÃ©cision attendu et validation exigÃ©e
- Ces contrÃ´les ne changent pas le schÃ©ma FV, son ordre ou les seuils physiques. La tolÃ©rance sur la somme des coefficients de traÃ®nÃ©e couvre uniquement les erreurs d'arrondi de leurs reprÃ©sentations en sortie ; elle n'est pas un critÃ¨re de convergence.
- Tests red-green : normes nÃ©gatives dans certificat/historique, continuitÃ© nÃ©gative, nombre d'itÃ©rations invalide et incohÃ©rence de dÃ©composition de traÃ®nÃ©e.
- Les campagnes physiques complÃ¨tes, la convergence itÃ©rative stricte N64/N128/N256, les profils $C_p/C_f$, le GCI par quantitÃ©, les conditions de sortie et l'incertitude expÃ©rimentale doivent Ãªtre prouvÃ©s indÃ©pendamment avant toute promotion.

## 7. Statut d'implÃ©mentation
- Lot local 2026-11-07 : gardes de qualification renforcÃ©s dans cinq scripts NACA, rÃ©gressions nÃ©gatives ajoutÃ©es.
- Tests Python NACA Ã©largis : 109 pass, 3 skipped ; les campagnes longues N128/N256 et la rÃ©fÃ©rence expÃ©rimentale ne sont pas exÃ©cutÃ©es ici.
- Statut : **partiel / non qualifiÃ© physiquement**. Aucun coefficient de traÃ®nÃ©e ou de cisaillement n'est promu sur ce seul lot.
