# FICHE MÉTHODE — Restart du FGMRES couplé et plateaux de résidu (chemin COUPLED+PCD)

> Fiche 16. Créée à l’occasion de l’audit du blocage Ghia 32×32 COUPLED+PCD identifié par la branche de diagnostic #894.

**En-tête d’évidence (gouvernance §16.5)**
- Scope : restart Krylov du système couplé pression-vélocité (FGMRES) du solveur incompressible stationnaire pour les requêtes explicites (COUPLED + Schur PCD) ; test de qualification N8-PCD (Ghia, Re=100, 32×32). Non couvert : qualité intrinsèque de l’approximation PCD au-delà du restart, chemins SIMPLE/SIMPLEC, tolérances et gates de qualification (inchangés).
- Statut : DRAFT (2026-10-14)
- Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json (INC-GHIA) ; NUMERICAL_METHOD_CAPABILITY_MATRIX.json (pcd : partial — non avancé par cette fiche)
- Issue / PR liées : #894 (branche de diagnostic — ne jamais merger), #896 (fiche 15), #897 (canal)
- Dernier commit vérifié : a22db9cc (2026-10-14)
- Limitations connues : verdict quantitatif 32×32 en attente du run validation-total au HEAD mergé (test derrière CFDX_ENABLE_LONG_VALIDATION).

## 1. Références sources (obligatoire : DOI ou ISBN)
- Y. Saad, M. H. Schultz, « GMRES: A Generalized Minimal Residual Algorithm for Solving Nonsymmetric Linear Systems », SIAM J. Sci. Stat. Comput. 7(3), 1986, pp. 856–869 — DOI: 10.1137/0907058 (méthode GMRES, restart et troncature de l’espace de Krylov).
- Y. Saad, « A Flexible Inner-Outer Preconditioned GMRES Algorithm », SIAM J. Sci. Comput. 14(2), 1993, pp. 461–469 — DOI: 10.1137/0914028 (FGMRES, préconditionnement variable imbriqué).
- D. Silvester, H. Elman, D. Kay, A. Wathen, « Efficient preconditioning of the linearized Navier–Stokes equations for incompressible flow », J. Comput. Appl. Math. 128(1–2), 2001, pp. 261–274 — DOI: 10.1016/S0377-0427(00)00515-X (approximation PCD du complément de Schur).

## 2. Formulation continue
Système stationnaire couplé (momentum + continuité assemblés en un système de point de selle) :

$$ \begin{pmatrix} A_{uu} & G \\ D & C \end{pmatrix} \begin{pmatrix} \mathbf{u} \\ p \end{pmatrix} = \begin{pmatrix} b_u \\ b_p \end{pmatrix} $$

avec bloc de pression Rhie–Chow implicite (anti-checkerboarding). Le cas Ghia est tout-parois : toutes les frontières de pression sont Neumann, la pression n’est définie qu’à une constante près — la jauge est imposée par élimination de la ligne ET de la colonne de la cellule de référence avec diagonale unité (assemblage couplé réduit).

## 3. Discrétisation Volume Fini
La discrétisation Volume Fini (trois composantes de vitesse + une pression par cellule) produit le système couplé de point de selle à 4·nx·ny inconnues (32×32 → 4096 ; 16×16 → 1024). FGMRES restarté appliqué à ce système construit m vecteurs de base par cycle :
- espace complet (m ≥ dimension) : convergence garantie en au plus n itérations en arithmétique exacte ;
- m < n : le taux de réduction atteignable par cycle est borné par le polynôme de degré m sur le spectre du système préconditionné ; un plateau au-dessus de la tolérance devient possible, même avec un bon préconditionneur.

Sélection du restart dans CFDX :
- Chemin automatique (krylov=Auto, preconditioner=Auto) : espace complet pour n ≤ 4096 — min(n, max(1024, budget)) — et cap mémoire 512 au-delà. Politique instaurée après deux plateaux mesurés et documentés dans steady_incompressible_solver.h : Poiseuille 768 inconnues (cap 512 → plateau ~1e-7) ; cavité 16×32, 4096 inconnues (cap 512 → résidu figé à 3.66e-12 contre une gate 1e-12, trois retries identiques).
- Chemin explicite (FGMRES et PCD demandés, cas du test N8-PCD) : restart = min(gmres_restart, dimension). La valeur par défaut de LinearSolverRequest::gmres_restart est 40 : une requête explicite qui ne fixe pas le restart hérite silencieusement d’un FGMRES(40), tronquant l’espace de Krylov très en dessous de la dimension du système — défaut observé sur la qualification 32×32 (#894 : budget de 5000 itérations épuisé en 125 cycles, résidu vrai relatif figé à 4.803775e-09 contre la gate 1e-10, trois retries identiques).
- Correctif : le test applique au chemin explicite la même politique que le chemin automatique (espace complet ≤ 4096 inconnues, borné par la dimension). Le budget d’itérations couplées (coupled_max_iterations=5000) et les gates ne sont pas modifiés.

## 4. Traitement des frontières
Le cas Ghia est tout-parois (lid entraîné) : aucune frontière à pression imposée, toutes les frontières de pression sont Neumann (flux nul). La jauge de pression est traitée à deux niveaux, non modifiés par ce correctif :
- Système couplé : élimination de la ligne ET de la colonne de la cellule de référence avec diagonale unité (assemblage couplé réduit) — la composante de jauge est imposée par construction.
- Préconditionneur PCD : les opérateurs de pression internes de l’action imbriquée portent une rangée unité à la cellule de référence (impose_pcd_reference_row) ; épingler la référence à l’entrée/sortie de l’action fige le résidu de jauge à sa valeur initiale (défaut documenté et corrigé dans pcd_schur.h).
- Conséquence pour le diagnostic : un plateau du résidu vrai relatif ne peut pas être attribué au traitement de la jauge — il relève de la troncature de l’espace de Krylov (cf. section 3) ou de la qualité de l’approximation PCD, pas des frontières.

## 5. Pièges connus (Gotchas)
- Signature du plateau de troncature : résidu vrai relatif figé juste au-dessus de la gate après épuisement du budget, identique sur les retries de relaxation (les alpha de relaxation ne changent pas le système linéaire couplé : les trois retries reproduisent le même résidu au bit près). Ne pas conclure à un défaut du préconditionneur sans avoir d’abord éliminé le restart.
- Le chemin explicite contourne la politique du chemin automatique : toute requête couplée explicite doit fixer gmres_restart consciemment (ou hériter de la politique miroir appliquée par ce correctif).
- L’action PCD imbriquée doit retourner une composante de jauge : les opérateurs de pression internes portent une rangée unité à la cellule de référence (impose_pcd_reference_row) ; épingler la référence à l’entrée/sortie de l’action fige le résidu de jauge à sa valeur initiale (défaut documenté et corrigé dans pcd_schur.h).
- Diagnostics : DiagnosticsControls::coupled_matrix_summary imprime le préconditionneur effectif et le restart (ligne COUPLED_PRECONDITIONER) — activé par ce correctif sur le chemin COUPLED du test pour toute preuve future.

## 6. Ordre de précision attendu et validation exigée
- Sans objet pour l’ordre spatial : ce correctif ne touche pas la discrétisation. La validation exigée est le contrat N8-PCD : convergence du solveur linéaire couplé (gate 1e-10 inchangée), profil Ghia Re=100 32×32 aux tolérances pré-déclarées (max 0.15 / RMS 0.075), continuité et résidu momentum (1e-7), le tout inchangé.
- Le test reste derrière CFDX_ENABLE_LONG_VALIDATION : la CI standard ne couvre que la compilation et les régressions SIMPLE (quick, order gate), non affectées par le restart couplé. Verdict quantitatif attendu du run validation-total au HEAD mergé.
- Risques résiduels : mémoire du restart complet (deux bases de 4096² double ≈ 270 Mo, déjà pratiqué par le chemin automatique à cette taille) et budget temps du test (TIMEOUT 720 s ; le run défaillant a pris 569 s pour 5000 itérations restart-40 — la convergence attendue en moins d’itérations doit rester dans le budget, sinon constat à consigner au registre avant tout ajustement).

## 7. Statut d’implémentation CFDX
- Fichiers source concernés : tests/validation/test_ghia_cavity.cpp (restart du chemin explicite + diagnostic couplé) ; src/cfdx/physics/steady_incompressible_solver.h (politique du chemin automatique — NON modifiée) ; src/cfdx/core/linalg/pcd_schur.h (contrat de jauge — non modifié).
- Tests exécutables : test_ghia_cavity_n8_pcd_preflight (16×16, TIMEOUT 120 s) et test_ghia_cavity_n8_pcd (32×32, TIMEOUT 720 s, dépend du preflight), labels n8;qualification;validation;long, derrière CFDX_ENABLE_LONG_VALIDATION.
- Écarts connus code ↔ fiche : le correctif est une correction de configuration Krylov fondée sur des preuves documentées (plateaux mesurés) ; si le run validation-total confirme la convergence 32×32, lever DRAFT. Si le plateau persiste avec le restart complet, le défaut est alors dans la qualité de l’approximation PCD (solveurs imbriqués, tolérances internes) — audit suivant, hors périmètre de cette fiche.
- Gaps et prochaines actions : qualification de l’ordre Fp vs Kp (note PROVISIONAL de pcd_schur.h) ; statut pcd du NUMERICAL_METHOD_CAPABILITY_MATRIX.json inchangé (partial) jusqu’au verdict validation-total.

## 8. Journal des modifications de la fiche
- 2026-10-14 — PR fix/ghia-pcd-restart-coupled (diagnostic #894) — création.
- 2026-10-14 — restructuration pour la gate d’audit : section 3 renommée « Discrétisation Volume Fini » (la sélection du restart y est intégrée comme algèbre du système discret) ; section 4 consacrée au « Traitement des frontières » (jauge tout-Neumann et rangée de référence PCD). Fond inchangé.
