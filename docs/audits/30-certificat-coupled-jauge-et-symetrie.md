# FICHE MÃTHODE 30 â Certification FV COUPLED, symÃ©trie et ancrage de pression

**En-tÃªte d'Ã©vidence**
- Scope : contrat de vÃ©rification des Ã©quations gelÃ©es de Picard en COUPLED, lignes de jauge du bloc CSR, existence effective des conditions aux limites de pression FIXED_VALUE en SIMPLE, COUPLED et PCD.
- Statut : QUALIFICATION NUMÃRIQUE CIBLÃE ; la qualification physique NACA N64/N128/N256 et la validation-total restent ouvertes.
- Base : archive utilisateur `cfdx_dev-master (1)(1).zip` reÃ§ue le 2026-11-07, et non HEAD distant supposÃ© Ã©quivalent.
- Sources : `src/cfdx/physics/{steady_incompressible_solver,coupled_frozen_fv_certificate,coupled_assembly_certificate}.h`.
- Preuve : Ã©checs reproduits puis corrigÃ©s dans les deux premiers tests (`test_n9_s9_rhie_chow_consistency`, `test_coupled_assembly_certificate`) ; rÃ©gression de prÃ©validation ajoutÃ©e dans `test_steady_incompressible_solver`.

## 1. RÃ©fÃ©rences sources
- C. M. Rhie et W. L. Chow, *AIAA Journal*, vol. 21, 1983 : DOI 10.2514/3.8284.
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016 : DOI 10.1007/978-3-319-16874-6 ; ISBN 978-3-319-16873-9.
- PETSc KSP : https://petsc.org/main/manual/ksp/ (normes des rÃ©sidus prÃ©conditionnÃ© et non prÃ©conditionnÃ©, modes nuls).
- NASA : https://www.grc.nasa.gov/www/wind/valid/tutorial/verassess.html (vÃ©rification du code, itÃ©ration, conservation, Ã©tude de maillage).

## 2. Formulation continue
Pour un Ã©coulement incompressible, la conservation de masse impose
\[\nabla\!\cdot\!\mathbf U=0.\]
La pression n'a de sens qu'Ã  une constante prÃ¨s en l'absence de donnÃ©e de pression absolue. Une condition de symÃ©trie impose le flux normal nul : \(\mathbf U\cdot\mathbf n=0\) Ã  la face, mÃªme lorsque la vitesse arbitraire du vecteur de test ne satisfait pas encore cette contrainte.

## 3. DiscrÃ©tisation Volume Fini
Pour chaque cellule \(P\), le rÃ©sidu de continuitÃ© non ancrÃ© est la somme signÃ©e des flux de masse de ses faces :
\[R_{m,P}=\sum_{f\in\partial P}\phi_f, \qquad \phi_f=\rho\,\mathbf U_f\cdot\mathbf S_f+\phi_{\mathrm{RC},f}.\]
Le bloc monolithique gelÃ© par Picard se prÃ©sente comme
\[\begin{bmatrix}M&G\\D&C\end{bmatrix}\begin{bmatrix}U\\p\end{bmatrix}=\begin{bmatrix}b_U\\b_p\end{bmatrix}.\]
Le certificat `coupled_frozen_fv_certificate.h` compare \(Ax-b\) Ã  une rÃ©Ã©valuation des Ã©quations de Picard sur les faces. Il ne constitue **pas** une comparaison du Jacobien de Newton au rÃ©sidu non linÃ©aire complet et ne prouve pas la justesse d'une Ã©quation de momentum erronÃ©e Ã  la source.

Une ligne de jauge imposÃ©e Ã  l'indice \(g\) doit Ãªtre \(A_{gg}=1\), \(A_{gj}=0\) si \(j\ne g\). Le certificat de structure exige dÃ©sormais aussi que la diagonale de Schur ajustÃ©e Ã  cet indice soit Ã©gale Ã  1 dans la tolÃ©rance dÃ©clarÃ©e. Ce contrÃ´le porte sur les mÃ©tadonnÃ©es utilisÃ©es par les prÃ©conditionneurs ; il ne modifie pas les coefficients physiques.

## 4. Traitement des frontiÃ¨res
- `EMPTY` et `SYMMETRY` : flux normal imposÃ© nul dans le certificat FV, en accord avec le bloc de continuitÃ© effectivement assemblÃ©. Une vitesse de test arbitraire ne doit pas crÃ©er un faux dÃ©faut de bilan sur ces faces.
- `FIXED_VALUE` pour \(p\) : cette condition ne supprime le mode nul que si son nom correspond Ã  un patch rÃ©el et non vide et que le patch n'est ni `EMPTY` ni `SYMMETRY`.
- Patch mal orthographiÃ©, sans faces, vide ou impermÃ©able : rejet explicite de cette configuration incompatible avant l'initialisation des champs utilisateur ; ne jamais dÃ©sactiver silencieusement le gauge sous prÃ©texte de prÃ©sence d'une entrÃ©e dans la table des CL.
- En l'absence d'une pression effectivement fixÃ©e : fixer une rÃ©fÃ©rence scalaire dans la matrice pression/Schur et dans la formulation COUPLED selon le contrat documentÃ© par le solveur.

## 5. PiÃ¨ges connus
- Une certification qui traite la face de symÃ©trie comme une sortie Ã  flux libre peut produire un **faux nÃ©gatif** de l'identitÃ© \(Ax-b=R_{FV}\).
- Une matrice dont la ligne de jauge est correcte peut prÃ©senter des mÃ©tadonnÃ©es de diagonale de Schur incohÃ©rentes et tromper les diagnostics de conditionnement.
- Une entrÃ©e `FIXED_VALUE` pour un patch absent peut crÃ©er un systÃ¨me singulier si elle supprime Ã  tort la jauge.
- Ne pas assimiler le rÃ©sidu linÃ©aire prÃ©conditionnÃ©, le vrai \(b-Ax\), la divergence de masse conservatrice et le dÃ©faut non linÃ©aire rÃ©Ã©valuÃ©. Toutes les mÃ©triques doivent Ãªtre vÃ©rifiÃ©es indÃ©pendamment.
- Ne pas annoncer une qualification physique sur la base d'un smoke test ou d'un certificat d'intÃ©gritÃ©, mÃªme si ces tests passent en Debug.

## 6. Ordre de prÃ©cision attendu et validation exigÃ©e
- Ces correctifs sont des corrections de **contrat et de contrÃ´le de validitÃ©** : ils ne changent ni l'ordre formel des opÃ©rateurs FV ni les seuils de convergence du solveur.
- Test red-green de la face `SYMMETRY` avec vitesse arbitrairement non tangentielle : le certificat doit Ãªtre cohÃ©rent avec la ligne de continuitÃ© de la matrice gelÃ©e.
- Test red-green de la diagonale de Schur de jauge corrompue Ã  zÃ©ro : le certificat doit refuser la preuve.
- Test de CL fantÃ´me et impermÃ©able : exception explicite, champs utilisateur U et p intacts si les contrÃ´les Ã©chouent.
- Tests de rÃ©gression de non-orthogonalitÃ©, RhieâChow, pression pure Neumann, AMG/PCD, SST, NACA court, thermique/rayonnement ; les suites longue N64/N128/N256, conservation, GCI, performances, MPI et sanitizers restent obligatoires avant une promotion globale.
- Aucune modification des critÃ¨res d'acceptation, budgets, tolÃ©rances ou valeurs cibles physiques n'est autorisÃ©e pour faire passer une campagne.

## 7. Statut d'implÃ©mentation
- 2026-11-07, base ZIP locale : correctifs intÃ©grÃ©s dans les trois entÃªtes et tests nÃ©gatifs rÃ©gressifs. Voir le rapport de qualification livrÃ© avec le patch exact.
- La suite CTest Debug ciblÃ©e et les jeux rapides NACA/thermiques sont verts sur cette copie. Il n'existe pas Ã  cette date, dans cette campagne locale, de preuve de convergence non linÃ©aire complÃ¨te N256 ni d'incertitude aÃ©rodynamique GCI.
- Pour fusion, vÃ©rifier la compatibilitÃ© des prÃ©images avec le HEAD exact et refaire les gates Release, DebugSanitizers et les campagnes physiques. Le prÃ©sent statut n'est **pas** une promotion en Â« parfaitement qualifiÃ© Â».
