# Outils d analyse VMFL036 axisymetrique (branche JAMAIS a merger)

Branche de sauvegarde des outils de diagnostic de la PR #895 (VMFL036 axisymetrique).
A SUPPRIMER au nettoyage une fois la solution merged. Ne jamais merger dans main.

## Contenu
- port3.ts : port TypeScript du solveur (jumeau numerique du commit 3475826).
  Deux bugs de transcription corriges le 15/10 vs la premiere version du port :
  (1) fermeture radiale du gradient ecrasee par l affectation suivante
      (sgr -= q*2*pi*A_mer AVANT la division, comme en C++) ;
  (2) stencil sphere inverse (le voisin ring-1 est sur une face INTERNE ;
      C++ saute les faces de bord).
- port3dbg.ts : port + instrumentation de debug (identites de fermeture du
  maillage, cellule de residu max, dump des faces).
- solver_3475826.h, axisymmetric_fvm.h : source C++ de reference (commit 3475826),
  copies byte-identiques (blob SHA verifies).
- test_axisymmetric_vmfl036_kernel.cpp : tests kernel C++ (reference de fidelite
  du port). ATTENTION : contient la corruption latente connue (backslash-n
  litteraux dans un commentaire qui avalent l assert visc.radial) - a reparer
  dans la PR d amelioration.
- test_vmfl036_axisymmetric.cpp : test de campagne C++ (gates + bandes).
- test_vmfl036_reference_case.cpp : test du cas de reference.

## Etat du diagnostic au 15/10 ~18:35Z
- Etabli (port) : 2 bugs de transcription corriges ; maillage du port verifie
  (identites sum_f s*A*n = (0, 2*pi*A_mer) exactes).
- Hypothese principale (non confirmee C++) : sous-budget du solveur de pression
  (SOR omega=1.8, 400 sweeps ; commentaires C++ : relres 48x72 6.5e-4@400, 8.2e-6@800
  vs gate pcorr 2e-7) -> NOT_CONVERGED probable au niveau pcorr/continuite.
- En attente : evidence du run 38779413390 (branche ci/long-campaign-diagnostics,
  logs complets de test_vmfl036_axisymmetric au head exact 3475826).
