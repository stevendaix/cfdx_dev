# diag-tools - Outils d'analyse offline (branche de diagnostic, NE JAMAIS MERGER)

Ce dossier rassemble les outils d'analyse developpes pour diagnostiquer l'echec
inattendu de test_vmfl036_axisymmetric dans "Long validation campaigns" au head
3475826 de la PR #895. Branche ci/long-campaign-diagnostics uniquement.
Ces outils seront supprimes lors du nettoyage final une fois le diagnostic termine.

## Contenu

- port3.ts : port TypeScript (jumeau numerique) du solveur
  src/cfdx/physics/axisymmetric_vmfl036_solver.h au commit 3475826.
  Transcription operation par operation (maillage polaire revolution, SIMPLE
  sous-relaxe, GS momentum 8 sweeps, SOR omega=1.8 x400 sweeps, correction
  conservative de flux). Zero caractere backslash dans le fichier pour garantir
  un aller-retour transport/stockage exact.
- README.md : ce fichier.

## Verification du port (15/10/2026)

1. Identites de fermeture du maillage EXACTES :
   - sum_f s*A*n_x ~ 0 : pire residu 1.1e-13 (faces internes et de bord)
   - sum_f s*A*n_r = 2*pi*A_mer : pire residu 2.2e-16
   Le maillage (radii = 0.5 + 49.5*(i/nr)^1.65) et la topologie sphere/outer/axis
   sont donc conformes au C++.

2. Deux bugs de transcription corriges dans port3 vs la premiere ebauche :
   a. gradients() : la fermeture radiale (sgr -= q*2*pi*A_mer AVANT la division
      par le volume) etait ecrasee par l'affectation suivante. Corrige pour
      coller exactement au C++ (lignes "Axisymmetric radial closure").
   b. Stencil sphere : la boucle de construction de sphereNb sautait les faces
      INTERNES au lieu des faces de bord (condition inversee par rapport au C++
      "if(f.neighbour==(-1)) continue;"). Corrige : sphereNb/popule pour les
      nt cellules ring-0.

3. Diagnostic it=1 (cas Re=50 medium 44x54, mu=0.02) : mom=1.543e+2,
   concentre sur l'anneau externe 43 (cellule 2333, j=11). C'est le residu
   convectif REEL du champ GS apres seulement 8 sweeps (flux equilibres,
   somme des F = 0 par construction, D ~ 43 devant F ~ 400). La valeur
   "mom=0.0305" du precedent run d'evidence (38743667702) N'EST PAS une
   reference valide : code d'avant le correctif de gradients (bug b30125e),
   et probablement cas 36x54.

## Cibles de validation (issues du port precedent, port2, perdu dans un reset d'environnement)

- Re=50 medium 44x54 upwind : Cd = 1.5992 (+1.31 pourcents vs 1.5785),
  cd_p ~ 0.985, cd_v ~ 0.615
- Re=100 coarse 24x36 : Cd ~ 0.81 ; medium 44x54 : 1.1869 ; fine 48x72 : 1.1918
- mom ~ 5.1e-5 a it=120 (medium) / 1.0e-5 a it=300 (fine)
- La convergence AUX GATES (mom 2e-6, pcorr 2e-7, cont 2e-6) n'a JAMAIS ete
  confirmeE ni au port ni en CI locale : c'est l'objet du diagnostic.

## Hypothese principale sur l'echec CI

NOT_CONVERGED d'un niveau (budget) : la contraction SOR mesuree dans les
commentaires du source C++ (48x72 : relres 6.5e-4 a 400 sweeps, 8.2e-6 a 800)
rend le gate pcorr 2e-7 difficilement atteignable a 400 sweeps, et la metrique
de continuite herite du residu de ligne du solveur de pression (identite des
flux conservatifs). Verdict attendu du run d'evidence 38779413390 (workflow
Long-campaign diagnostics, branche ci/long-campaign-diagnostics, logs committes
sous diag-evidence/).

## Discipline

- Aucune tolerance/gate modifie sans accord explicite de Steven.
- Les extensions de budget (pressure_sweeps, max_outer_iterations) sont
  legitimes si justifiees mathematiquement (contraction SOR mesuree).
- Campagnes lancees par chunks avec checkpoints (budget temps sandbox ~35 s).
