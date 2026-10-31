# FICHE 24 — Microscope de diagnostic NACA0012 / solveurs scalaires

**Statut : INSTRUMENTATION AJOUTÉE — CAUSE RACINE EN ATTENTE DES ARTEFACTS CI**

## 1. Symptômes observés sur la tête précédente

Campagne `CFDX total validation`, job `Long validation campaigns`, commit
`ca0180472e7b57e2f616f406d96d580124754693` (31 octobre 2026).

- `test_naca0012_qualification` rapporte une valeur de qualification provisoire
  `NACA0012_RE1000_A0 level=64` : Cd=0.0901496, Cl=-9.15e-11,
  Cd pression=0.0344606, Cd visqueux=0.055689, continuité=8.84e-11,
  456 itérations. Ce message ne constitue pas à lui seul un PASS final.
- Une correction de pression échoue à l'itération non linéaire 20 :
  status=1, 457 itérations, résidu=3.021189409e2, relatif=1.181456e-9.
- À l'itération 22, le solve Ux échoue :
  status=3 (DIVERGED), itérations=1, résidu rapporté=1.329e25.
  Les retries suivants rapportent des relatifs 5.482 puis 2.252 avant l'abandon.
- Le test de cisaillement rapporte provisoirement Cd_visqueux GG=0.055689,
  WLS=0.0580625, 1D=0.106308 (spread relatif=0.690069), puis échoue
  sur le même problème Ux. Ce spread doit être réévalué sur des solutions
  convergées et comparables ; il ne permet pas encore d'isoler l'erreur de gradient.

Ces traces établissent une défaillance linéaire/non linéaire réelle, mais ne
permettent pas de conclure si l'origine est la matrice assemblée, le choix
Krylov, le préconditionneur, le RHS ou l'état initial.

## 2. Modifications de diagnostic

### CG

- Toutes les sorties numériques `DIVERGED` après le calcul du résidu initial
  recalculent désormais le résidu vrai FP64 du système courant.
- Le statut public reste `DIVERGED` ; aucun nouvel enum, critère, tolérance ou
  budget n'est introduit.
- Un `SolverDiagnosticSample` conserve l'itération et les scalaires
  `r^T M^{-1}r`, `p^T A p`, `alpha`, ainsi que le résidu vrai et le résidu
  récursif lorsqu'ils sont disponibles.
- Le chemin `NOT_APPLICABLE` garde son contrat et utilise le même calcul de
  résidu vrai.

### Solve scalaire / prédicteur momentum

Sur tout échec final d'un solve scalaire, le journal émet :

- `SCALAR_LINEAR_ATTEMPT` pour le solve demandé et chaque fallback réellement
  exécuté (GMRES, Gauss-Seidel, CG), avec statut, itérations, résidus et dernier
  échantillon Krylov disponible ;
- `SCALAR_LINEAR_FAILURE` avec dimensions, nnz, norme infinie de matrice,
  défaut de symétrie absolu, bornes diagonales, norme RHS, nombre de coefficients
  et RHS non finis, solveur/préconditionneur demandés et politique de fallback ;
- `SCALAR_LINEAR_KRYLOV_SAMPLE` avec les scalaires de l'échantillon de breakdown.

Cette instrumentation est déclenchée uniquement à l'échec et ne change pas
l'ordre des solveurs, les vecteurs candidats, les critères de convergence ni
les portes physiques.

## 3. Test de non-régression ajouté

`cg_nonfinite_initial_preconditioned_dot_reports_true_residual` injecte une
action de préconditionneur finie mais suffisamment grande pour faire déborder
le produit scalaire initial. Le test exige :

- statut `DIVERGED`, zéro itération ;
- solution inchangée ;
- résidu vrai FP64 et résidu relatif correctement renseignés ;
- échantillon diagnostique conservant le produit préconditionné négatif ou non fini.

## 4. Prochaine lecture des artefacts

Sur le nouveau HEAD, examiner dans cet ordre :

1. `SCALAR_LINEAR_ATTEMPT` : déterminer quel solveur a été demandé et lesquels
   ont réellement été essayés ;
2. `SCALAR_LINEAR_FAILURE` : vérifier l'échelle, les diagonales, les valeurs
   non finies et la symétrie de la matrice Ux ;
3. `SCALAR_LINEAR_KRYLOV_SAMPLE` : localiser le premier scalaire qui diverge
   (`rMz`, `pAp`, `alpha`) ;
4. comparer le résidu vrai au résidu récursif et le RHS à la matrice ;
5. seulement ensuite décider d'une correction de formulation, de préconditionnement
   ou de sélection de Krylov.

La cause racine reste **NON DÉTERMINÉE** tant que ces nouvelles traces n'ont pas
été obtenues. Ne pas assouplir les tolérances pour faire disparaître l'échec.
