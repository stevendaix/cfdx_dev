# FICHE MÉTHODE 30 — Certification FV COUPLED, symétrie et ancrage de pression

**En-tête d'évidence**
- Scope : contrat de vérification des équations gelées de Picard en COUPLED, lignes de jauge du bloc CSR, existence effective des conditions aux limites de pression FIXED_VALUE en SIMPLE, COUPLED et PCD.
- Statut : QUALIFICATION NUMÉRIQUE CIBLÉE ; la qualification physique NACA N64/N128/N256 et la validation-total restent ouvertes.
- Base : archive utilisateur `cfdx_dev-master (1)(1).zip` reçue le 2026-11-07, et non HEAD distant supposé équivalent.
- Sources : `src/cfdx/physics/{steady_incompressible_solver,coupled_frozen_fv_certificate,coupled_assembly_certificate}.h`.
- Preuve : échecs reproduits puis corrigés dans les deux premiers tests (`test_n9_s9_rhie_chow_consistency`, `test_coupled_assembly_certificate`) ; régression de prévalidation ajoutée dans `test_steady_incompressible_solver`.

## 1. Références sources
- C. M. Rhie et W. L. Chow, *AIAA Journal*, vol. 21, 1983 : DOI 10.2514/3.8284.
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016 : DOI 10.1007/978-3-319-16874-6 ; ISBN 978-3-319-16873-9.
- PETSc KSP : https://petsc.org/main/manual/ksp/ (normes des résidus préconditionné et non préconditionné, modes nuls).
- NASA : https://www.grc.nasa.gov/www/wind/valid/tutorial/verassess.html (vérification du code, itération, conservation, étude de maillage).

## 2. Formulation continue
Pour un écoulement incompressible, la conservation de masse impose
\[\nabla\!\cdot\!\mathbf U=0.\]
La pression n'a de sens qu'à une constante près en l'absence de donnée de pression absolue. Une condition de symétrie impose le flux normal nul : \(\mathbf U\cdot\mathbf n=0\) à la face, même lorsque la vitesse arbitraire du vecteur de test ne satisfait pas encore cette contrainte.

## 3. Discrétisation Volume Fini
Pour chaque cellule \(P\), le résidu de continuité non ancré est la somme signée des flux de masse de ses faces :
\[R_{m,P}=\sum_{f\in\partial P}\phi_f, \qquad \phi_f=\rho\,\mathbf U_f\cdot\mathbf S_f+\phi_{\mathrm{RC},f}.\]
Le bloc monolithique gelé par Picard se présente comme
\[\begin{bmatrix}M&G\\D&C\end{bmatrix}\begin{bmatrix}U\\p\end{bmatrix}=\begin{bmatrix}b_U\\b_p\end{bmatrix}.\]
Le certificat `coupled_frozen_fv_certificate.h` compare \(Ax-b\) à une réévaluation des équations de Picard sur les faces. Il ne constitue **pas** une comparaison du Jacobien de Newton au résidu non linéaire complet et ne prouve pas la justesse d'une équation de momentum erronée à la source.

Une ligne de jauge imposée à l'indice \(g\) doit être \(A_{gg}=1\), \(A_{gj}=0\) si \(j\ne g\). Le certificat de structure exige désormais aussi que la diagonale de Schur ajustée à cet indice soit égale à 1 dans la tolérance déclarée. Ce contrôle porte sur les métadonnées utilisées par les préconditionneurs ; il ne modifie pas les coefficients physiques.

## 4. Traitement des frontières
- `EMPTY` et `SYMMETRY` : flux normal imposé nul dans le certificat FV, en accord avec le bloc de continuité effectivement assemblé. Une vitesse de test arbitraire ne doit pas créer un faux défaut de bilan sur ces faces.
- `FIXED_VALUE` pour \(p\) : cette condition ne supprime le mode nul que si son nom correspond à un patch réel et non vide et que le patch n'est ni `EMPTY` ni `SYMMETRY`.
- Patch mal orthographié, sans faces, vide ou imperméable : rejet explicite de cette configuration incompatible avant l'initialisation des champs utilisateur ; ne jamais désactiver silencieusement le gauge sous prétexte de présence d'une entrée dans la table des CL.
- En l'absence d'une pression effectivement fixée : fixer une référence scalaire dans la matrice pression/Schur et dans la formulation COUPLED selon le contrat documenté par le solveur.

## 5. Pièges connus
- Une certification qui traite la face de symétrie comme une sortie à flux libre peut produire un **faux négatif** de l'identité \(Ax-b=R_{FV}\).
- Une matrice dont la ligne de jauge est correcte peut présenter des métadonnées de diagonale de Schur incohérentes et tromper les diagnostics de conditionnement.
- Une entrée `FIXED_VALUE` pour un patch absent peut créer un système singulier si elle supprime à tort la jauge.
- Ne pas assimiler le résidu linéaire préconditionné, le vrai \(b-Ax\), la divergence de masse conservatrice et le défaut non linéaire réévalué. Toutes les métriques doivent être vérifiées indépendamment.
- Ne pas annoncer une qualification physique sur la base d'un smoke test ou d'un certificat d'intégrité, même si ces tests passent en Debug.

## 6. Ordre de précision attendu et validation exigée
- Ces correctifs sont des corrections de **contrat et de contrôle de validité** : ils ne changent ni l'ordre formel des opérateurs FV ni les seuils de convergence du solveur.
- Test red-green de la face `SYMMETRY` avec vitesse arbitrairement non tangentielle : le certificat doit être cohérent avec la ligne de continuité de la matrice gelée.
- Test red-green de la diagonale de Schur de jauge corrompue à zéro : le certificat doit refuser la preuve.
- Test de CL fantôme et imperméable : exception explicite, champs utilisateur U et p intacts si les contrôles échouent.
- Tests de régression de non-orthogonalité, Rhie–Chow, pression pure Neumann, AMG/PCD, SST, NACA court, thermique/rayonnement ; les suites longue N64/N128/N256, conservation, GCI, performances, MPI et sanitizers restent obligatoires avant une promotion globale.
- Aucune modification des critères d'acceptation, budgets, tolérances ou valeurs cibles physiques n'est autorisée pour faire passer une campagne.

## 7. Statut d'implémentation
- 2026-11-07, base ZIP locale : correctifs intégrés dans les trois entêtes et tests négatifs régressifs. Voir le rapport de qualification livré avec le patch exact.
- La suite CTest Debug ciblée et les jeux rapides NACA/thermiques sont verts sur cette copie. Il n'existe pas à cette date, dans cette campagne locale, de preuve de convergence non linéaire complète N256 ni d'incertitude aérodynamique GCI.
- Pour fusion, vérifier la compatibilité des préimages avec le HEAD exact et refaire les gates Release, DebugSanitizers et les campagnes physiques. Le présent statut n'est **pas** une promotion en Â« parfaitement qualifié Â».
