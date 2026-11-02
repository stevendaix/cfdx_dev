# FICHE MÉTHODE 26 — Globalisation non linéaire SIMPLE (garde de croissance physique) et reconstruction pression/traction aux faces de paroi

**En-tête d'évidence**
- Scope : (a) garde opt-in de rejet des itérations non linéaires SIMPLE dont les résidus physiques (QDM relatif RHS, continuité normalisée) croissent sévèrement, avant publication au moniteur de convergence, avec restauration transactionnelle et rejeu par la reduction de sous-relaxation existante ; (b) reconstruction affine de la pression et de la traction visqueuse au centre des faces de paroi (post-traitement d'audit, hors solveur).
- Statut : DRAFT (2026-10-28, PR #1146 — port byte-exact du patch audité P0/P1 du 2026-10-28)
- Issues liées : audit cinq chantiers (garde SIMPLE, reconstruction aux faces, NACA --single, GCI fail-closed, diagnostic 3D) ; PRs #1145 (fixtures, base du port), #1146 (cette PR).
- Dernier commit vérifié : 015c5097ff02c6315998f5653c026cbca8e2e4f6 (branche PR #1146) ; base master ecf1d022.
- Limitations : aucune qualification physique n'est revendiquée (r48/r72 non convergés, Cd non validé, conduction linéaire 3D FAIL documenté). La traction pariétale du composant de face est d'ordre 1 en champ quadratique.

## 1. Références sources
- F. Moukalled, L. Mangani, M. Darwish (2016), ISBN 978-3-319-16873-9, ch. 10–12 (SIMPLE, sous-relaxation, critères de convergence des schémas segregés).
- J. H. Ferziger, M. Perić (2002), DOI 10.1007/978-3-642-56026-2, § 5.4 et § 8.4 (procédés itératifs imbriqués, rejet d'itération et globalisation).
- P. L. Roe (1986), « Characteristic-based schemes for the Euler equations », Annual Review of Fluid Mechanics 18, DOI 10.1146/annurev.fl.18.010186.002233 (motivation : ne jamais publier un état aux grandeurs physiques non physiquement croissantes sans diagnostic).

## 2. Formulation continue
Le point fixe de SIMPLE n'est pas modifié : la correction de pression et la sous-relaxation déterminent le chemin d'itération, pas l'équilibre continu. Le garde est un critère d'acceptation de candidat d'itération : un candidat dont les résidus physiques (équation de quantité de mouvement relative, continuité normalisée) croissent sévèrement entre deux itérations acceptées successives est rejeté avant toute publication au moniteur de convergence, et l'itération est rejouée avec une relaxation d'équation réduite par le contrôleur de retry existant. La reconstruction aux faces de paroi évalue la traction continue tau = mu (grad u + grad u^T) n et la pression au centre géométrique de face, exactes pour des champs affines.

## 3. Discrétisation Volume Fini
Le garde s'insère entre les gates de résidus internes et la mise à jour du moniteur : comparaison previous vs candidat (history.back()), rejet par NonlinearRetryableFailure si croissance sévère au-delà de la tolérance relative du contrat de convergence, restauration par la transaction NonlinearStateRollback (U, p, mass_flux, moniteur, historique, échantillons, SST) puis rejeu avec alpha_u/alpha_p réduits. SIMPLE uniquement ; aucune sous-relaxation indépendante U/p post-solve n'est réintroduite pour COUPLED. Le composant de face reconstruit p_f = p_P + grad p · (x_f − x_P) et la traction par correction rank-one de la normale ; il est exact sur champs affines (MMS) et d'ordre 1 en champ quadratique (correction de Hessien hors scope).

## 4. Traitement des frontières
La reconstruction est locale à chaque face du patch de paroi (owner uniquement, géométrie de face exacte) : pas de stencil transverse, pas de valeur limite plombée requise. Le rejet d'itération préserve les états déjà acceptés : un candidat rejeté ne peut ni écraser l'historique de convergence, ni publier un mass_flux appartenant à un état roulé back (autoritative_mass_flux inchangé jusqu'à acceptation).

## 5. Pièges connus
- Un garde sur le résidu *d'équation* (résidu linéaire de la matrice) ne détecte pas les croissances physiques : le critère porte sur le résidu physique relatif RHS et la continuité normalisée.
- Le bruit d'arrondi près de la convergence ne doit pas être rejeté : le seuil est relatif à la tolérance du contrat de convergence (tolérance relative × facteur du contrôleur de retry), testé explicitement (blowup réel rejeté, bruit d'arrondi accepté).
- Un candidat non fini (Inf/NaN) doit être rejeté avant le moniteur : le garde est fail-closed sur les non-finis.
- Le moniteur et l'historique sont tentatifs jusqu'au commit de la transaction : rejouer une itération avec le même index ne doit jamais voir d'échantillons périmés.
- La traction pariétale d'ordre 1 en quadratique : NE PAS promouvoir le composant de face comme post-processeur d'ordre 2 sans Hessien/WLS cohérent en paroi ; le test quadratique documente l'ordre mesuré, il ne le gate pas.

## 6. Ordre de précision attendu et validation exigée
- Garde SIMPLE : propriétés testées en unitaire (tests/unit/test_timestep_control.cpp) — rejet d'un blowup réel, non-rejet du bruit d'arrondi, rejet des non-finis ; comportement par défaut strictement inchangé (opt-in via NonlinearRetryControls).
- Reconstruction aux faces : MMS affine sur le maillage NACA réel (tests/unit/test_wall_face_reconstruction_mms.cpp, --face-mms de l'audit de cisaillement) — erreur < 1e-11 sur pression et traction ; géométrie dégénérée rejetée. Ordre mesuré quadratique : pression 2, traction 1 (tests/unit/test_wall_face_reconstruction_quadratic.cpp, ratios exacts 4/2 documentés, la promotion d'ordre 2 en traction échoue volontairement).
- Toute qualification physique (GCI NACA, Cd, campagnes r*) reste hors scope de cette fiche : l'estimateur GCI (scripts/qualify_naca_grid.py) est fail-closed (3 grilles convergées, ratios quasi uniformes, Cd = Cdp + Cdv) et aucune campagne n'est revendiquée ici.

## 7. Statut d'implémentation CFDX
- Sources : src/cfdx/physics/timestep_control.h (options + severe_physical_residual_growth), src/cfdx/physics/steady_incompressible_solver.h (garde avant convergence_monitor.update), src/cfdx/physics/wall_face_reconstruction.h (reconstruct_wall_face_linear, post-traitement).
- Tests : tests/unit/test_timestep_control.cpp, tests/unit/test_wall_face_reconstruction_mms.cpp, tests/unit/test_wall_face_reconstruction_quadratic.cpp, tests/validation/test_naca0012_wall_shear_audit.cpp (--face-mms), tests/python/test_qualify_naca_grid.py.
- Écarts code ↔ fiche : aucun — port byte-exact du patch audité P0/P1 du 2026-10-28 (chaque blob SHA-1 vérifié contre l'index du patch).
- Gaps assumés : qualification r48/r72 et n128/n256 non produites (r48 : MAX_ITERATIONS documenté) ; conduction linéaire 3D T=300+10z FAIL documenté (résidu matriciel ~290, diagnostic opt-in --linear-conduction fail-closed, exit 2) ; pas de promotion du composant d'ordre 2 en traction.

## 8. Journal
- 2026-10-28 — création (PR #1146) : synthèse le l'audit cinq chantiers P0/P1 ; statut DRAFT en attente de review, aucune qualification revendiquée.
