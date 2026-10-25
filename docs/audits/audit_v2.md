## Reprise opérationnelle — 25 octobre 2026

Cette section est le point d'entrée courant. Les statuts ci-dessous sont basés sur les HEAD et workflows observés le 25/10/2026 ; une campagne `skipped` n'est jamais une validation physique réussie.

| Sujet | HEAD / PR | Constat vérifié | Action obligatoire |
|---|---|---|---|
| F03 Rhie–Chow / COUPLED | #1078 | N9 physique réussi ; campagne N9 qualification et validation totale encore en cours lors du contrôle. | Attendre les deux campagnes et auditer le flux de continuité avant toute promotion. |
| F06 diffusion non orthogonale | #1103, HEAD `080e87d60749b95c0424d151afecfd585a74d7a0` | Nouveau HEAD poussé ; tous les jobs de cette révision étaient encore en file d'attente. Un ancien essai COUPLED skew Couette plafonnait à `continuity_linf=0.001547` après 33 itérations. | Examiner les logs du HEAD exact. Si l'échec se reproduit, corriger le flux/gradient ou le couplage ; ne pas relâcher les gates. |
| F11 contrainte transposée à viscosité variable | #1102, HEAD `5afd56a664fd1a81cb0f5fd705a8fb18e273b94f` | La branche cible est `fix/audit-1069-coupled-rhie-chow-contract`, divergente de master (25 commits en avance, 13 en retard). CI sur cette base a échoué sur `test_dual_time_navier_stokes` et `test_n9_negative_diagnostics`; ces résultats ne permettent pas d'isoler F11. | Porter F11 proprement sur master, puis relancer CI exact-HEAD avant d'évaluer le code physique. |
| F25 SST API | #1113 | CI classique et audit verts sur le HEAD observé ; maturité numérique encore en attente. N9 physique, qualification et total sont `skipped`. | Attendre maturité et ne pas présenter les qualifications ignorées comme réussies. |
| F26 source P1 | #1111 | Audit et garde-fous verts ; CI principale/maturité encore en cours au contrôle. N9/total `skipped`. | Attendre les workflows et garder l'oracle indépendant. |
| F27 courbure VOF | #1114 | Audit et N10 verts ; CI, GUI, documentation et maturité encore en file/en cours. Qualification géométrique complète non couverte. | Attendre CI puis garder le statut PARTIEL tant que l'intégration de production n'est pas qualifiée. |
| F30 facteurs de vue S2S | #1115 | Audit, GUI, N10 et garde-fous verts ; CI principale/maturité encore en cours. Qualification physique `skipped`. | Attendre les jobs restants ; la diagonale nulle ne qualifie pas une enceinte fermée. |
| Contrôle du registre | #1108 | L'ancien run a échoué car l'audit gate traitait `audit_v2.md` et `turbulence.md` comme des fiches méthode ; la correction correspondante a été fusionnée via #1112, après ce run. | Relancer les workflows sur une révision actuelle et vérifier les contrôles de contenu sur le HEAD courant. |

Règle de clôture : chaque ligne ne passe à CORRIGÉE que si le code courant a été audité, les tests ciblés passent et les checks obligatoires passent sur le HEAD exact. Une campagne physique ignorée, en cours ou échouée reste ouverte.

### Mise à jour — 25 octobre 2026 ~09:05Z (suite) : F30 mergée, GUI palette mergée, F31/F33 mergées, F32 dirty

Cette suite poursuit la section ~08:20Z ci-dessous (les merges du matin y restent exacts).

- **Push docs master `71858a20` vérifié** : le check-run du commit (« Build numerical test suite ») est `success` (54 s) ; l'unique annotation est le warning de dépréciation Node.js 20 des actions officielles, sans défaut.
- **F30 #1115 MERGÉE** : squash `b7ced399`, `expectedHeadSha` épinglé sur `92af08d2` ; CI du head exact 37/37 check-runs terminés, 0 échec / 0 pending (Pre-commit, Validate, Audit gate vérifiés verts ; validation-total `skipped` sans label — la qualification physique de l'enceinte fermée reste ouverte). **#1086 fermée comme supersédée** (delta intégralement couvert par #1115, commentaire posté).
- **#1120 MERGÉE** : squash `976d8b11`, head `20df917d` ; 36/36 check-runs terminés, 0 échec / 0 pending — les 6 jobs absents du canal check-runs dégradé (Pre-commit, Validate, Release, DebugSanitizers, GUI/Qt, Build numerical) ont été confirmés verts via l'API check-runs paginée. L'incident GUI noté à ~08:20Z est résolu par ce merge.
- **F31 #1080 MERGÉE** : squash `f90de151`, head `755d7ac0` ; 37/37 check-runs terminés, 0 échec / 0 pending (Pre-commit, Validate, Audit gate, Release, DebugSanitizers, GUI, E2E inclus).
- **F33 #1089 MERGÉE** : squash `78e2724b`, head `13dbbd71` ; 25/25 check-runs terminés, 0 échec / 0 pending.
- **F32 #1088 NON mergée** : `mergeable_state=dirty` (conflit réel apparu après les merges du jour — ce n'est pas un échec CI). Aucune réconciliation automatique effectuée ; rebase/réconciliation à instruire séparément, aucun critère numérique modifié.
- **#1118 (CHT géométrie) toujours draft** : aucune action ; #1109 reste ouverte en attendant le passage ready de #1118 (sa réconciliation manuelle reste annulée).

Règles inchangées : aucune tolérance, oracle ou gate modifié ; les campagnes `validation-total` restent `skipped` sans label et ne comptent jamais comme validations physiques réussies ; INCONCLUSIVE ≠ PASS.

### Mise à jour — 25 octobre 2026 ~09:55Z (suite) : F32 #1088 réconcilée byte-exact puis mergée — vague transport F31/F32/F33 complète

- **F32 #1088 réconciliée sur master** : le `mergeable_state=dirty` noté à ~09:05Z provenait de l'évolution de master dans les mêmes fichiers — #1080 (F31, squash `f90de151af`) et #1089 (F33, squash `78e2724b`) dans `finite_volume_transport.h`/`test_finite_volume_transport.cpp`, et `038e3aef` (distances asymétriques) dans `thermal_zones.h`/`test_thermal_zones_vv.cpp`, le nouvel oracle asymétrique occupant exactement le point d'insertion du test F32. Résolution en 3 étapes : `0f588204` alignement des 4 fichiers sur master, `548c11f89` update-branch sans conflit, `0a455abe` ré-application du delta F32 complet. 5 fichiers re-fetchés byte-identiques (FNV-1a 64 : `9c31fe78fa46c44a` solver_status.h, `2c697c534266c004` finite_volume_transport.h, `a4d73a02f18b8910` thermal_zones.h, `245dc70cbb7b56f4` test_finite_volume_transport.cpp, `8ac99ace3e50c945` test_thermal_zones_vv.cpp). Aucun oracle/tolérance/gate modifié ; le test `thermal_conduction_rejects_underrelaxed_unconverged_state` est réinséré à sa position d'origine (après `composite_wall`, avant le nouvel oracle asymétrique). Commentaire de consolidation : 6395799481.
- **F32 #1088 MERGÉE** : squash `1c39147c`, `expectedHeadSha` épinglé sur `0a455abe` ; CI du head exact terminée : 26 ✅ / 0 échec / 0 pending / 11 skipped (validation-total sans label = normal), jobs démarrés 09:18Z+ (frais du head réconcilié, pas de cache) ; mergeable_state `clean` au moment du merge (09:53:19Z).
- **Conséquence registre** : F31/F32/F33 sont portées par master au commit `1c39147c` (statut CORRIGÉE au sens « code audité + tests ciblés + checks obligatoires verts sur le HEAD exact » UNIQUEMENT). La qualification physique (campagnes `validation-total`, enceinte fermée S2S, ordres de convergence) reste NON acquise ; INCONCLUSIVE ≠ PASS.

### Mise à jour — 25 octobre 2026 ~10:10Z (suite) : CHT #1118 mergée (voie parallèle), #1121 GUI mergée, doublons SST/harnais fermés, résidus adjudiqués

- **CHT #1118 MERGÉE (voie parallèle)** : squash `dcd19b2d` (09:57Z) — interface CHT géométrie/normales dans `cht_solver.h`, `test_cht_validation.cpp`, `test_level_c_coupled_verification.cpp`. **#1109 fermée en supersession** (09:58Z) conformément au verdict byte-level déjà documenté (commentaire 6395268529). #1119 fermée, remplacée par le draft docs #1122 (tracker `turbulence.md`, voie parallèle — non touché).
- **#1121 MERGÉE (GUI Lot 3 G1a)** : squash `27c79784`, `expectedHeadSha` épinglé sur `6e00b23c` ; CI du head exact : 26 ✅ / 0 échec / 0 pending / 1 skipped (deploy, attendu), jobs 09:45–09:49Z frais ; `mergeable_state` clean au moment du merge. Fichiers GUI Python uniquement, orthogonaux aux commits master intervenus depuis sa base.
- **Fermées pour supersession byte-level** (commentaires 6396263512, 6396263992, 6396264434) : **#1094** — le limiteur omega est couvert par #1116 (`9fa9bfb9`), master étant plus riche (garde zéro-nut, oracle actif `551.2258064516129`, cas zéro-nut testé) ; **#1093** — F34 couvert (garde isfinite de `expect_close` + self-test `expect_comparator_rejects`, 5 cas) ; **#1098** — garde isfinite complète `cfdx_expect_near_a_…` du harnais déjà sur master.
- **#1099 FERMÉE, résidu docs sauvé** : code supersédé par #1116 ; fiche 09 §2 couvre déjà SST-2003 (avec branche ν_t=0). Le paragraphe « SST-2003 production contract » absent de `docs/validation/TURBULENCE_CLOSURE_EQUATIONS.md` a été intégré byte-exact au commit `470ddb81` (+607 octets, P̃k=min(ρνtS², 10β*ρkω), P_ω=γP̃k/νt avec garde ν_t=0, réf NASA TMR tmbwg.github.io). Commentaire 6396287415.
- **#1095 RE-SCOPÉE (reste ouverte)** : contrat R2 géométrique intégralement couvert sur master (conformité vertex, opposition des normales, aires, imperméabilité WALL/INTERFACE, conductivités, gardes isfinite ; le delta `test_thermal_challenge_vv.cpp` `e2.conductivity=1` est déjà appliqué). **Résidu unique : time-consistency R2** — la garde Δt identique et fini entre régions dans `solve_two_region_cht` et son test négatif sont absents de master (`cht_solver.h` ne référence aucun `.dt`). La PR reste le véhicule de ce résidu ; aucun critère numérique (tolérance/oracle/gate) concerné. Commentaire 6396287982.
- **Tâches planifiées 08:50Z et 09:35Z : completed** — merges #1115/#1120/#1080/#1089 et consolidation #1120 faits ; plus de risque de collision.

Règles inchangées : aucune tolérance, oracle ou gate modifié ; INCONCLUSIVE ≠ PASS ; campagnes `validation-total` restent `skipped` sans label et ne comptent jamais comme validations physiques réussies ; jamais de push sur les branches de la voie parallèle (#1103, #1102, #1092, #1095, #1096, #1101, #1078, #1122).

### Mise à jour — 25 octobre 2026 ~08:20Z : vagues de merges, doublons fermés, reste à faire

Cette mise à jour supplante les lignes F25/F26/F27 et « Contrôle du registre » du tableau ci-dessus (rédigées avant les merges de la matinée).

**Merges effectués dans master (squash, SHA de head épinglés, CI 100 % terminée et verte sur le head exact au moment du merge ; campagnes validation-total restant `skipped` sans label)** :

- `7fec3d97` #1112 — classification des rapports consolidés (gate d'audit) ; `72ddaab9` #1106 — GUI undo/redo borné ; `a30c95b7` — rejet des références de qualification non finies (turbulence) ;
- `d720662e` #1110 — CHT conductivité interface ; `9fa9bfb9` #1116 — SST production limiter appliqué à la source omega ;
- `e9d27214` #1104 — tracker turbulence.md ; `a761292f` #1113 — F25 : suppression de la surcharge SST ambiguë (26 checks verts sur le head fixé, N1–N16 inclus) ;
- `a672e021` #1111 — F26 : normalisation source P1 (23 checks verts, 0 échec, 0 en attente sur le head exact) ;
- `dd863b43` #1114 — F27 : contrat de courbure VOF axe-par-axe + fiche 25 (21 checks verts, 0 échec, 0 en attente) ;
- `be74325b` #1108 — contrôle du registre : l'action « relancer les workflows sur une révision actuelle » est satisfaite (merge postérieur à la correction #1112 de la gate).

**Doublons fermés sans merge (supersession vérifiée par lecture byte-level du code master, commentaires sur chaque PR)** : #1105, #1107 (matin) ; #1079 → couverte par #1111 ; #1081 → couverte par #1113 ; #1083 → couverte par #1114 (fiches en plus) ; #1097 → couverte par #1110. #1086 (S2S, branche empilée) sera fermée après le merge de #1115.

**En cours** :

- F30 #1115 : conflit avec #1111/#1113 résolu par réconciliation byte-exact sur master courant (head `92af08d2`, alignement → update-branch → re-application du delta, FNV-1a vérifié) ; CI relancée ~08:00Z sur ce head exact — merge uniquement CI verte.
- CHT géométrie/normales : #1118 (draft, base master courant) couvre intégralement le delta restant de #1109 (vérifié byte-level : 106 lignes incluses dans les 121 de #1109, les 15 excédentaires étant déjà dans master via #1110 ; test level-C identique). La réconciliation manuelle de #1109 est ANNULÉE ; #1109 sera fermée au merge de #1118.
- #1120 (robustesse du test GUI palette aux actions internes Qt sans nom) : CI en vol ; un échec du job GUI observé à un instant — à instruire sur les checks stabilisés.
- #1119 (docs audit) : draft.

**Reste à faire (ordre proposé)** :

1. F31/F32/F33 transport : #1080 / #1088 / #1089 (non-draft, non supersédées) — relever l'état CI exact-HEAD, merger si verts, sinon diagnostiquer.
2. F30 : merger #1115 dès CI verte, puis fermer #1086.
3. CHT : suivre #1118 jusqu'au merge, puis fermer #1109.
4. F06 #1103 : attribuer l'échec « N9 Physical Validation » (exit 8, run 40172900788) — pré-existant vs causé par F06 — sur le HEAD exact en évolution (draft).
5. F11 #1102 : portage propre sur master (branche divergente, échecs non isolables sur la base actuelle), puis CI exact-HEAD.
6. F03 #1078 : attendre les campagnes N9 qualification + validation totale, puis auditer le flux de continuité avant promotion.

Aucune tolérance, oracle ou gate n'a été modifié ; aucune qualification physique ignorée n'est comptée comme réussie (règle de clôture ci-dessus inchangée).

## Suivi d'exécution courant — 24 octobre 2026 (mise à jour depuis R2)

Cette section est le journal vivant de reprise. Elle complète les statuts historiques R2 ci-dessous sans réécrire les preuves anciennes. Un PR ouvert signifie « correction candidate à auditer », jamais « problème corrigé ». Un item n'est clos qu'après lecture du diff actuel, vérification de l'absence de doublon sur `master`, tests ciblés et CI sur le HEAD exact.

### État des reprises sur GitHub

| ID de l'audit | État constaté à la reprise | Action / gate suivante |
|---|---|---|
| F03 | PR #1078 ouverte ; contrat Rhie–Chow/COUPLED à qualifier | Auditer le flux de face, le RHS de continuité, la jauge et les cas pression affine ; exiger CI exact-HEAD et comparaison COUPLED/ségrégé indépendante. |
| F06 | PR draft #1103 ouverte | Auditer intégration réelle dans tous les appelants de transport, reconstruction et diffusion ; ajouter MMS polyédrique et vérifier Release/DebugSanitizers avant promotion. |
| F08 | **Code candidat déjà présent sur master** : `geometric_interpolation_weight` calcule une projection longitudinale et `test_finite_volume_transport.cpp` contient un contre-exemple distance euclidienne vs projection | Pas de nouveau correctif avant vérification ciblée : confirmer le test est enregistré/exécuté en CI et que le chemin `compute_gradient_gauss` consomme bien ce helper ; si oui, passer à F12/F13. |
| F11 | PR draft #1102 ouverte | Vérifier les trois prédicteurs, le résidu final indépendant, le signe et les flux opposés ; compléter MMS non orthogonal/polyédrique avant qualification. |
| F12/F13 | Aucun PR dédié visible dans l'inventaire ouvert consulté | Auditer les contrats dimensionnels des helpers momentum/énergie et leurs appelants de production avant de proposer un patch. |
| F18 | Aucun PR dédié visible dans l'inventaire ouvert consulté | Établir les identités thermodynamiques attendues pour `beta != 0`, puis tester dérivées et invariants dimensionnels. |
| F22–F24 | Aucun PR dédié visible dans l'inventaire ouvert consulté | Séparer les variantes physiques nommées (WALE, compressibilité, Kato–Launder) ; pas de correction générique sans formule et variables nécessaires. |
| F25 | PR #1081 ouverte | Auditer le diff et les appels de l'API supprimée ; vérifier les builds/CI avant clôture. |
| F27 | PR #1083 ouverte | Vérifier les dérivées par axe et les oracles plan/cercle/sphère ; garder la qualification VOF/Laplace distincte. |
| F30 | PR #1086 ouverte | Vérifier producteur→consommateur S2S, réciprocité, fermeture et fraction ouverte ; la cohérence algébrique seule ne qualifie pas la géométrie. |
| F32 | PR #1088 ouverte | Recalculer le résidu du champ effectivement retourné après relaxation et auditer le gate thermique R3. |
| F34 | PR #1100 ouverte, avec variantes historiques #1093/#1098 encore ouvertes | Éviter les doublons ; vérifier le comparateur courant, puis compléter les oracles indépendants manquants sans compter plusieurs sorties d'un même défaut comme plusieurs bugs. |
| F02 | PR #1096 ouverte ; limiteur documenté comme partiel | Distinguer la correction du limiteur de diffusivité de la qualification du modèle k-omega 2006. |
| F28 | PR #1110 créée depuis le master courant ; anciennes #1095/#1097/#1107 restent des branches candidates périmées/concurrentes | Auditer la CI de #1110 sur son HEAD exact ; fusionner uniquement si verte et mergeable, puis fermer/annoter les variantes obsolètes sans perdre les tests uniques. |
| F26 | PR #1079 ouverte | Le correctif local R2 n'est pas une preuve de la PR ; vérifier formule, oracle et CI exact-HEAD. |
| F01 / SST production | PR #1105 ouverte ; intégration de boucle/restart distincte (#1045/#1068) | Séparer source omega, intégration de production, persistance/restart et qualification physique ; aucune promotion globale SST sans les quatre gates. |

### Ordre d'exécution maintenu

1. **F08** : audit initial du code courant effectué. Le helper utilise déjà la projection longitudinale et un test compare explicitement projection et distance euclidienne. Vérifier l'enregistrement/CI du test et les appelants ; si la preuve est complète, marquer F08 comme corrigé dans le périmètre démontré et passer à F12/F13.
2. **F03 / F06 / F11** : auditer les PR existantes avant toute réimplémentation ; pas de cherry-pick aveugle depuis les branches historiques.
3. **F12 / F13 / F18 / F22–F24** : spécifier d'abord les contrats physiques et dimensionnels, puis ajouter des oracles indépendants.
4. **F25 / F27 / F30 / F32 / F34 / F02 / F26 / F28** : examiner les PR existantes et leur CI exact-HEAD, éliminer les doublons et ne fermer que sur preuve.
5. Après chaque étape, actualiser la ligne correspondante du registre de la section 5 et cette table ; garder les limitations de validation explicites.


### Journal de revue des PR — campagne de clôture (mise à jour)

| PR | Verdict de revue | Décision actuelle |
|---|---|---|
| #1100 — comparateur de qualification F34 | Correctif ciblé : rejet explicite NaN/Inf, tolérance non finie/négative ; CFDX CI, maturité numérique, Fast Guards et audit maturité verts sur HEAD `e0b4ca94a7292189e2bdf9eca6153cb34402a374`. Les campagnes N9/total étaient ignorées, pas réussies. | **Fusionnée** par squash, commit `a30c95b7248cb86dbcd3142a41830798e294d201`. La validation physique reste hors portée de ce test de comparateur. |
| #1103 — diffusion non orthogonale F06 | **Échec reproductible** sur HEAD `6fbf8be1e490f4e0c361462d91b311277d3b2ad8`. `test_nonorthogonal_laplacian_campaign` échoue dans `transport_nonorthogonal_linear_campaign` : le résidu corrigé ne satisfait pas l'oracle affine sur le maillage shear + warp non uniforme. | **Ne pas fusionner.** Corriger la reconstruction du gradient/flux sur maillage déformé ; maintenir l'oracle strict, ajouter preuve locale de gradient affine et relancer Release + DebugSanitizers + N3 sur le HEAD corrigé. |
| #1097 / #1110 — contrat de conductivité CHT F28 | #1097 a des checks verts sur un ancien HEAD mais n'est plus mergeable. J'ai porté le correctif audité sur le master courant dans #1110 : validation finie/positive/cohérente avant mutation, test négatif sans mutation, fixture thermique cohérente et mise à jour de la fiche. CI #1110 en attente sur HEAD `963fb7d745de44a9064b77f625c4c1291f26481b`. | **#1110 en validation** ; ne pas fusionner avant CI exacte. Éviter de fusionner #1097/#1107 en parallèle. |
| #1107 — contrat CHT empilé | Head `25b6d58c536e920ea676d2f1d49fb6253d1f1653`, PR draft et `mergeable=false`; aucun run CI PR déclenché trouvé sur ce SHA lors de la revue. Diff empile les changements CHT avec des modifications SST et de gate documentaire. | **Ne pas fusionner en bloc.** Rebaser/extraire uniquement les tests CHT supplémentaires sur la correction retenue ; séparer le SST et la classification de rapport si ces changements ne sont pas strictement nécessaires. |
| #1098 — comparateur F34 / régression F08 | Variante ouverte qui chevauche #1100 sur les comparateurs non finis, mais ajoute également des changements de test transport/F08. CI verte sur son ancien HEAD, pas une preuve de compatibilité avec le master après #1100. | Garder uniquement les tests F08 encore utiles ; supprimer la duplication du comparateur déjà fusionné, puis rebaser et relancer la CI. Ne pas fermer avant d'avoir isolé le delta non redondant. |
| #1099 / #1094 — source omega SST | Variantes historiques concurrentes ; #1099 présente des échecs CFDX CI, maturité numérique et Fast Guards sur son HEAD `7ad610f9e3d189c6db8912c0098d12318982523d`. #1094 est une variante plus ancienne avec CI verte mais campagnes physiques ignorées. | Ne pas fusionner ces branches. Auditer #1105 comme candidat consolidé, attendre tous les runs en cours, et qualifier l'intégration de production séparément de la formule de source. |
| #1078 — Rhie–Chow/COUPLED F03 | Les suites principales sont vertes sur l'ancien HEAD, mais Fast Guards est annulé, la PR reste draft et GitHub signale `mergeable=false`. | Ne pas fusionner : rebaser sur master courant, compléter Fast Guards et vérifier l'égalité des flux de face avec le bloc de continuité et l'oracle pression affine. |
| #1079, #1081, #1083, #1086, #1088 — F26/F25/F27/F30/F32 | Les CI principales sont vertes sur leurs anciens HEAD, mais GitHub signale `mergeable=false` pour les cinq PR. Les campagnes N9 et validation totale sont ignorées. | Aucun merge aveugle. Rebaser chaque changement indépendamment, confirmer qu'il n'est pas déjà sur master, exécuter les tests ciblés et la CI sur le nouveau HEAD ; garder la qualification physique explicitement distincte. |
| #1102 — contrainte visqueuse transposée F11 | CFDX CI, Fast Guards, audit gate, documentation et maturité numérique verts sur le HEAD examiné ; la PR est draft et sa propre fiche indique que l'oracle est seulement un maillage orthogonal uniforme. | Ne pas fusionner comme correction « terminée ». Ajouter MMS sur maillage déformé/polyédrique et contrôler le résidu final réassemblé ; promouvoir uniquement après ces gates. |
| #1096 / #1101 — dossiers consolidés de corrections | CI verte sur les HEAD consultés, mais les scopes se chevauchent partiellement (tests CHT/énergie et éléments d'audit) ; #1096 indique explicitement que des tests ne sont pas enregistrés tant que les PR dépendantes n'ont pas été intégrées. | Comparer les arbres et bases exacts avant fusion. Ne pas fusionner les deux dossiers ensemble ni traiter le nombre de tests historiques comme une campagne sur master courant. |

**Règle ajoutée à la clôture :** un statut vert sur un ancien SHA n'autorise pas une fusion après déplacement de `master`. La base, le diff et les checks doivent être vérifiés de nouveau sur le HEAD candidat final. « skipped » pour N9, validation totale ou campagne physique ne signifie jamais « passed ».

---


### Reprise après création de #1110

- #1100 a été fusionnée dans `master` (commit `a30c95b7248cb86dbcd3142a41830798e294d201`).
- #1110 porte le correctif F28 sur la base courante et attend CI ; ne pas conclure avant résultats sur `963fb7d745de44a9064b77f625c4c1291f26481b`.
- #1103 reste bloquée par un contre-exemple N3 sur le maillage transformé shear+warp ; l'erreur ne doit pas être masquée en modifiant le seuil.
- Les PR #1078, #1079, #1081, #1083, #1086 et #1088 ont des CI historiques majoritairement vertes, mais leur base est ancienne et GitHub signale `mergeable=false`. La reprise doit réconcilier chaque delta sur le master courant avant toute fusion.


### Nouveaux correctifs portés sur le master courant

| PR | Objet | État / gate |
|---|---|---|
| #1110 | F28 — cohérence conductivité CHT, porté sur master courant | CI du HEAD `963fb7d745de44a9064b77f625c4c1291f26481b` en cours ; l'audit gate initial a échoué car il traite `audit_v2.md` et `turbulence.md` comme des fiches méthode. Ce défaut de classification est séparé du code CHT. |
| #1111 | F26 — normalisation du terme source P1, porté sur master courant | CI en cours sur HEAD `e68b2056ec910d2fbc757dbc05137e11766136c7`. Oracle indépendant : source sans irradiation égale à (4\sigma T^4), et source nulle lorsque (J=I_b). Aucune fusion avant CI verte. |



| #1115 | F30 — estimateur de facteurs de forme S2S | Port courant : diagonale strictement nulle, réciprocité surfacique conservée, tests sur surfaces d'aires égales/inégales et validation par le consommateur S2S. CI en attente sur le HEAD final ; les rangées non fermées sont interprétées comme échange avec l'espace ouvert non modélisé, pas comme une qualification d'enceinte fermée. |
| #1114 | F27 — contrat de courbure VOF | Port courant de la correction : l'API accepte explicitement les dérivées (\partial_x n_x,\partial_y n_y,\partial_z n_z), supprime l'ancien stencil ambigu et teste les cas plan/cercle/sphère. CI en attente sur HEAD `4c42092a8bdb9e3871ace2ce8a1aa06b78e6f8c3`. Limite déclarée : ne qualifie pas une chaîne de reconstruction VOF multidimensionnelle complète. |
| #1113 | F25 — suppression de la surcharge SST ambiguë | Correctif porté sur master courant : surcharge historique supprimée au niveau compilation, appels tests/benchmarks migrés vers l'API avec strain-rate explicite, oracle numérique indépendant. CI exacte en attente sur HEAD `6f23e7dbd42d2b91d3d19b47d95ba43e0910625f`. |
| #1112 | Gate d'audit — rapports vs fiches méthode | PR créée sur master courant ; liste d'exclusion explicite limitée aux rapports `audit_v2.md` et `turbulence.md`, sans relâcher les exigences des fiches méthode. CI en attente sur HEAD `9dceab6781871e70a5fcd69c126e85e26c47648b`. |

### Garde-fous de reprise

- L'état R2 est fondé sur une archive ZIP figée et un patch cumulatif ; il ne doit pas être appliqué tel quel à `master` ou empilé sur une PR sans comparaison exacte des sources.
- L'inventaire GitHub courant montre des PR candidates, mais la liste des PR ouvertes ne suffit pas à établir leur diff, leur base/head ou leur CI. Aucune correction n'est déclarée fusionnée par cette section.
- Aucun seuil de convergence, budget, oracle ou test ne doit être affaibli pour obtenir un vert. Les validations longues restent séparées et explicitement étiquetées.

---

CFDX — Audit consolidé, correctifs appliqués et tests — R2
24 octobre 2026. Fichier de travail unique, révision R2.
0. Lire, appliquer et prolonger ce document
Les corrections sont maintenant dans ce fichier, pas seulement décrites ni renvoyées vers un ZIP externe. L'annexe D contient le patch cumulatif complet de 43 fichiers, y compris les sources des cinq tests d'audit. L'annexe C contient le script permettant d'extraire ce patch du Markdown, de vérifier le ZIP d'origine, de créer une copie et de lancer les tests. Il faut le ZIP source exact, mais pas l'ancien ZIP de correctifs.
Le registre de la section 5 et les fiches de la section 6 constituent l'état courant des problèmes. La section 10 explique le travail effectivement ajouté dans cette révision : F28 corrigé dans le périmètre testé, H05 nouvellement identifié et corrigé, 113 nouveaux contrôles, trois fixtures de validation rectifiées, et réexécution des 167 contrôles antérieurs.
Convention d'historique : les résultats V1 des sections 3–4 et des annexes A–B sont conservés et étiquetés comme antérieurs. Les 214/220 résultats CTest de V1 ne sont pas une campagne relancée sur R2. Les statuts de PR de V1 ne sont pas une vérification de leurs nouveaux commits. Aucune PR n'a été créée, modifiée ou fusionnée dans cette révision.
Pour une prochaine révision, conserver les identifiants F/G/A/H, mettre à jour la fiche et la ligne du registre, ajouter une preuve avant/après, puis régénérer ensemble le patch et son manifeste. Ne pas promouvoir un état sur le seul titre d'une PR.
1. Verdict et contenu courant
Le registre comprend 53 constats distincts : 34 corrigés et testés dans leur périmètre explicite, 2 partiellement traités et 17 non corrigés ou restant à qualifier. Ce comptage concerne les constats, pas 53 équations physiques distinctes ; plusieurs constats partagent un opérateur.
Le patch cumulatif modifie ou ajoute 43 fichiers par rapport à l'archive de référence. La copie R2 contient 1181 fichiers ordinaires : les 1174 fichiers du ZIP et sept nouveaux fichiers de tests/fixtures. L'archive source reste intacte. Le script refuse un autre ZIP ou un dossier de destination existant.
280/280 contrôles internes des cinq suites d'audit passent sur R2 ; les 10/10 tests CTest sélectionnés passent. Les cinq suites font partie des dix exécutables : ces comptes ne s'additionnent pas. Les deux nouveaux exécutables, soit 113 contrôles, passent aussi avec ASan/UBSan. La bibliothèque statique de base reste compilée en Release non instrumentée : la vérification sanitizer porte sur les nouvelles unités de test et les fonctions inline qu'elles compilent, pas sur l'intégralité de CFDX.
Les quatre blocages principaux restent F03, F06, F08 et F11 : COUPLED/Rhie–Chow, diffusion non orthogonale du transport, gradient Gauss décalé et tenseur visqueux variable. Les corrections R2 ne les ferment pas. MPI/GPU, maillages CHT non conformes et qualification physique globale restent hors du périmètre validé.
2. Provenance, versions et sens des statuts
•	Archive source : cfdx_dev-master.zip.
•	Commit porte par le commentaire du ZIP : 3f7a764799ba4db5fbb5e1495b85e8fc5c31c6bb.
•	SHA-256 du ZIP : 11f92e9f842d16b184432051208c40262a0961085d86ac933946b54613b955bb.
•	Compilateur : GCC 14.2.0 ; CMake 3.31.6 ; C++17 ; CPU Release ; MPI/GPU/bindings Python desactives ; HDF5 serial detecte.
•	CFDX_BUILD_RUNTIME_TESTS=OFF ne signifie pas disparition du test CPU test_performance_runtime : celui-ci est conserve et sa dependance MPI injustifiee est corrigee.
•	Les lignes citees dans le registre sont celles de l'archive originale. Les positions apres correction sont visibles dans le patch cumulatif embarque (annexe D).
Corrige et teste = cause du contre-exemple supprimee et tests precis executes sur la copie. Ce n'est ni une preuve pour tous les maillages, ni une validation d'une PR GitHub, ni une validation physique d'un modele complet. Partiel = portion explicite corrigee, autre exigence ouverte. Non corrige = aucun correctif de cette livraison ne doit etre suppose traiter le point.
Les comptes des campagnes precedentes (25 CTest, 578 pytest, 19 CTest, 38 comparaisons, 14 CTest) ne sont pas additionnes a cette campagne : plusieurs tests se recouvrent et leurs versions different. Le nombre historique de 1418 correspondait a des entrees ZIP incluant des repertoires ; le nombre controle ici est 1174 fichiers ordinaires.
3. PR : etat documentaire V1 conserve, non actualise par R2
3.1 Liste publique consultee pendant V1, le 24 octobre 2026
Lors de V1, une consultation avait retourne 19 PR ouvertes, dont les suivantes directement liees a nos audits. Les tentatives d'ouvrir les pages individuelles, les fichiers modifies, les patches et l'API ont echoue. Aucun SHA de tete recent n'a donc ete telecharge, compile ni teste. Un titre de PR n'est pas une preuve de correction.
PR affichee ouverte
Anomalies visees par le titre
Verdict de cette revue
#1092
G01–G09 — Reconstruction et Krylov
Intention visible ; correction et CI exact-HEAD non verifiees.
#1078
F03 — Compensation Rhie–Chow et bloc de continuite
Intention visible ; correction et CI exact-HEAD non verifiees.
#1079
F26 — Normalisation source P1
Intention visible ; correction et CI exact-HEAD non verifiees.
#1080
F31 — Separation des desequilibres global/local
Intention visible ; correction et CI exact-HEAD non verifiees.
#1081
F25 — Suppression de la surcharge SST ambigue
Intention visible ; correction et CI exact-HEAD non verifiees.
#1083
F27 — Courbure VOF par axe
Intention visible ; correction et CI exact-HEAD non verifiees.
#1084
F29 — Distances asymetriques de resistance de contact
Intention visible ; correction et CI exact-HEAD non verifiees.
#1086
F30 — Contrat estimateur de facteurs de vue S2S
Intention visible ; correction et CI exact-HEAD non verifiees.
#1088
F32 — Residu apres sous-relaxation
Intention visible ; correction et CI exact-HEAD non verifiees.
#1089
F33 — Regression de bornage convection
Intention visible ; correction et CI exact-HEAD non verifiees.
#1045
Integration SST — Boucle de production transactionnelle
Intention visible ; correction et CI exact-HEAD non verifiees.
#1068
Integration SST — Raccordement production et restart, WIP
Intention visible ; correction et CI exact-HEAD non verifiees.
#1055
Etage energie — Echec solveur distinct de qualification
Intention visible ; correction et CI exact-HEAD non verifiees.
Source publique : https://github.com/stevendaix/cfdx_dev/pulls ; capture logique conservee dans remote/latest_pr_inventory.json. Les acces en echec ne sont pas interpretes comme une absence de correctif. Le patch local est independant de ces branches et risque de les recouvrir ; ne pas l'empiler aveuglement sur #1092 ou sur les PR Fxx. Porter d'abord les tests sur leurs SHA exacts, puis comparer les differences restantes.
3.2 Corrections et limites deja documentees dans l'archive
Element
Preuve dans l'archive
Ce que cela permet de conclure
#511
QUALIFICATION_REGISTRY.md, phase 2 : fermee le 12 octobre, remplacee par #863.
La presentation anterieure de #511 comme PR encore en cours etait obsolete. Pas de preuve sur le HEAD courant de #863.
#571
Meme registre : fermee sans fusion le 12 octobre.
Sa reduction d'horizon temporel n'est pas un correctif accepte de cette archive.
#855 / #856
Campagnes N9/N8 et SHA consignes dans le registre.
Preuves documentaires de campagnes historiques ; elles ne testent pas automatiquement nos contre-exemples nouveaux.
#887 / #888
Etages T1/T2 READY, corrections de bilans locaux et de contrat energetique consignees.
READY n'est pas PASS. Le solveur energie sait recevoir T^n ; le caller CHT oubliait encore de le figer.
Gradient quadratique principal
Le noyau normalise donne deja une reponse correcte au changement d'unite.
QUICK utilisait un autre fit ; une correction dans le premier ne ferme pas G07/G08.
BDF2 du noyau
Correct sur pas variable dans le controle independant.
Le pilote double temps generique etait encore incorrect, F19.
P1 principal
Normalisation correcte dans la branche principale.
Le helper auxiliaire restait faux, F26.
3.3 Conditions de cloture a exiger pour les PR importantes
#1092 : executer les 87 controles numeriques sur le SHA de tete exact, ajouter notre cas CG a RHS 1e-200, vrai residu FP64 en precision mixte, zero RHS avec x0 non nul, toutes les politiques LS/WLS, QUICK petit stencil et changement d'unite. Il faut aussi examiner H01–H04, qui ne sont pas couverts par le seul titre G01–G09.
#1078 : verifier ensemble la semantique U/HbyA, le flux de face et le bloc de continuite. Un changement de post-traitement seul ne ferme pas F03. Exiger pression affine avec U constant, invariance a la jauge, bilan local et comparaison COUPLED/segregue sur le meme operateur. Les deux calculs compares ne doivent pas appeler le meme helper fautif.
#1088 : distinguer statut du predicteur lineaire et statut du champ accepte sous-relaxe ; recalculer le vrai residu de ce dernier sans casser la logique de la boucle non lineaire. Ne pas transformer arbitrairement toute sous-relaxation en echec lineaire fatal.
#1089 : combiner tests de bornitude physique et identite assembleur/flux pour les 7 schemas, modes bornes/non bornes et signes de flux. Les 56 tests de cette livraison apportent l'identite discrete, pas a eux seuls une preuve de monotonie globale.
#1083, #1086, #1084 : respecter respectivement les axes de divergence, reciprocite/fermeture/visibilite S2S, et les deux distances de conduction. Refuser les cas non pris en charge plutot que fabriquer une valeur plausible.
4. Verification V1 conservee et limites historiques
Cette section decrit V1. Pour les executions effectuees sur R2, voir la section 10.
4.1 Comparaisons avant/apres
Cas
Code original
Copie corrigee
Reference independante
Suite numerique nouvelle
58 echecs sur 87 controles
87/87
Analytique, lois constitutives et invariances
Reconstruction flux CENTRAL
UBSan : lecture de pointeur nul, ligne 629
56/56 combinaisons schemas/bornage sans anomalie ASan/UBSan
Flux local vs matrice, sources incluses
CG Laplace, Gamma=1e24
Converge en 0 iteration, vrai residu relatif 1
10 iterations, vrai residu 2.54e-16
phi=x
CG Laplace, Gamma=1e12
Vrai residu 1.46e-5 malgre tolerance 1e-10
3.88e-16
Invariance a Gamma constant
CG, A=b=1e-200
x=0 annonce converge
x=1
Solution scalaire exacte
Precision mixte, (1+1e-8)x=1
Residu annonce nul, erreur reelle 1e-8
Residu FP64 <=1e-12
Equation d'origine FP64
WLS, phi=x, exclusion frontieres
gx=0.5 au bord
gx=1
Gradient analytique
QUICK, maillage 2^3
Repli lineaire nul
Champ affine retrouve
phi=x
QUICK, echelle 1e-8
Courbure perdue
Erreur normalisee ~2.7e-14
phi=(x/L)^2
Conductance non uniforme
3.2
4.571428571
Resistances en serie
CHT un pas dt=0.1
T1=315.277778 ; T2=384.722222
308.333333 ; 391.666667
Systeme thermique 2x2
CHT asymetrique, 6 cas
Echecs et mauvais temps physique
24/24 controles
Deux capacites, sources et conductance distinctes
Source omega SST limitee
5555.55556
161.290323
SST-2003 avec Pk limite
Wilcox, forte deformation
nu_t=1
0.00342857143
Limiteur Wilcox, hors correction Pope
Enthalpie a Cp borne
Incoherente avec Cp utilise
Integrale du Cp borne
Primitive/quadrature independante par morceaux
RK vectoriel
Seule composante 0 avancee
3 composantes et ordres 2/3 retrouves
u'=u et derivee constante
Newton, un pas permis
Racine atteinte, statut echec
Racine et statut converge
R(u)=u-1
Diagonale matrix-free gamma=0
Diagonale 1, action 0
Diagonale 0, action 0
e_i^T A e_i
Les anciens programmes de sondage impriment aussi les anomalies laissees ouvertes (par exemple RC_PHYSICAL_VELOCITY_INPUT, WALE auxiliaire et VOF). Leur execution n'est pas comptee comme un succes global. Le commentaire historique de CHT_INTERFACE_ITERATIONS decrit le bug initial ; la preuve du temps correct reside dans les temperatures et le systeme 2x2, pas dans le nombre d'iterations d'interface.
Les grandes echelles numeriques sont des tests d'invariance, pas des parametres recommandes pour un tutoriel. Un schema peut encore etre mal conditionne meme quand sa certification de residu est correcte.
4.2 Campagne CTest finale
220 tests selectionnes ; 214 reussites ; 6 echecs ; 0 ignores. Journal complet : logs/final_ctest.log, XML : final_ctest.xml.
Echecs non dissimules :
•	test_n11_conservation_boundedness
•	test_cross_import_equivalence
•	test_openfoam_results_cli_production_e2e
•	test_production_solver_e2e
•	test_meshio_bridge
•	test_mesh_geometry_matrix
Le statut de ces tests est analyse dans l'annexe de campagne ; aucun echec n'est automatiquement classe environnemental sans sa trace.
Les cinq campagnes non executees sont test_n7_adaptive_relaxation_re400, test_ghia_cavity_order, test_polyhedral_poisson_convergence, test_n9_physical_matrix et test_n9_poiseuille_body_force_production. Leur configuration et leurs tolerances restent inchangees. Les timeouts eventuels des autres tests sont comptabilises en echecs, pas en succes ni en invalidation physique automatique.
4.3 Evolution des tests, sans masquer les erreurs
La premiere passe de 88 tests existants a donne 86 reussites et 2 echecs : interpolation et conservation transport. La campagne plus large a ensuite expose un appel WLS qui comptait sur le defaut de politique implicite ; le caller de test demande maintenant explicitement ZERO_GRADIENT_GHOST et verifie separement le rejet de rang sous EXCLUDE_BOUNDARY. Le vecteur des erreurs de raffinement, auparavant recree a chaque maillage, est aussi conserve pour rendre la verification de decroissance effective. L'interpolation imposait une demi-somme sur une pseudo-geometrie ouverte et coplanaire ; elle utilise maintenant deux cubes fermes, avec les memes egalites attendues et une inversion d'orientation en plus. Le test de flux sortant utilisait la valeur de Dirichlet au lieu de l'amont owner, contrairement a la matrice ; son attendu est devenu 6=3*2, tolerance 1e-15 inchangee. Le test SST/Wilcox a ete corrige d'apres la definition physique, avec les anciennes traces conservees. Aucun seuil de precision existant n'a ete elargi.
La campagne CHT plus exigeante a expose H04, qui n'etait pas connu au debut : certaines solutions 1x1 correctement arrondies etaient refusees parce qu'une tolerance relative etait comparee a un residu absolu. Le correctif a ete applique puis les six cas relances a tolerance 1e-12 identique.
Les nouvelles reductions de maxima d'erreur rejettent explicitement NaN/Inf : un std::max ne doit pas masquer une comparaison non finie.
MPI multirang, CUDA, champs distribues reels, performances de grande taille, toutes les variantes physiques et les derniers HEAD des PR ne sont pas qualifies par cette livraison. Les stencils deficitaires n'ont pas ete transformes en stencils de plein rang.

5. Registre unique des problemes - etat R2
ID stable
Priorite
Probleme
Etat du correctif local
PR publique ciblee, non testee
F01
P1
Production de omega non limitee dans SST-2003
CORRIGE ET TESTE (perimetre indique)
—
F02
P1
La fermeture nommee komega2006 omet le limiteur de Wilcox
PARTIEL
—
F03
P1
La branche COUPLED utilise une vitesse physique comme H/A dans le flux de pression
NON CORRIGE / A QUALIFIER
#1078
F04
P1
Le couplage CHT transitoire avance le temps a chaque iteration d'interface
CORRIGE ET TESTE (perimetre indique)
—
F05
P1
Moyenne harmonique non ponderee sur cellules non uniformes
CORRIGE ET TESTE (perimetre indique)
—
F06
P1
La diffusion du solveur n'utilise pas le laplacien non orthogonal corrige
NON CORRIGE / A QUALIFIER
—
F07
P1
Interpolation LINEAR centrale 1/2 non exacte sur grille etiree
CORRIGE ET TESTE (perimetre indique)
—
F08
P1
Le poids du gradient de Gauss confond distance et projection
NON CORRIGE / A QUALIFIER
—
F09
P1
Le flux convectif sortant reconstruit ne correspond pas a la matrice
CORRIGE ET TESTE (perimetre indique)
—
F10
P0
Lecture d'un gradient vide pendant la reconstruction de flux
CORRIGE ET TESTE (perimetre indique)
—
F11
P1
Diffusion de vitesse scalaire incomplete lorsque la viscosite varie
NON CORRIGE / A QUALIFIER
—
F12
P1
Helper momentum: convection integree melangee a des termes volumiques
NON CORRIGE / A QUALIFIER
—
F13
P1
Helper thermique: convection en W et conduction en W/kg
NON CORRIGE / A QUALIFIER
—
F14
P1
Integrale d'enthalpie tabulee incorrecte hors table
CORRIGE ET TESTE (perimetre indique)
—
F15
P1
L'enthalpie n'integre pas le Cp effectivement borne
CORRIGE ET TESTE (perimetre indique)
—
F16
P1
Cp par defaut et conductivite de l'air errones dans une surcharge transport
CORRIGE ET TESTE (perimetre indique)
—
F17
P2
Deux IdealGasEOS renvoient des derivees incompatibles
CORRIGE ET TESTE (perimetre indique)
—
F18
P2
Identites incompatibles dans IncompressibleEOS pour beta non nul
NON CORRIGE / A QUALIFIER
—
F19
P1
BDF2 a pas variable incorrect dans le pilote generique double temps
CORRIGE ET TESTE (perimetre indique)
—
F20
P1
Une iteration implicite peut accepter un NaN comme convergence
CORRIGE ET TESTE (perimetre indique)
—
F21
P1
Les integrateurs low-storage ne mettent a jour que la composante 0
CORRIGE ET TESTE (perimetre indique)
—
F22
P2
Une fonction nommee WALE calcule la loi de Smagorinsky
NON CORRIGE / A QUALIFIER
—
F23
P2
Correction de compressibilite turbulente non adimensionnelle
NON CORRIGE / A QUALIFIER
—
F24
P2
Correction Kato-Launder ajoutee avec des dimensions incompatibles
NON CORRIGE / A QUALIFIER
—
F25
P2
Surcharge historique SST sans invariant de deformation
NON CORRIGE / A QUALIFIER
#1081
F26
P2
Facteur pi manquant dans une source P1 auxiliaire
CORRIGE ET TESTE (perimetre indique)
#1079
F27
P2
Courbure VOF: derivees de composantes sur le meme axe
NON CORRIGE / A QUALIFIER
#1083
F28
P2
Interface CHT: distances et conductivites de deux contrats differents
CORRIGE ET TESTE R2 (interface conforme, conduction)
—
F29
P1
Resistance thermique de zone suppose des demi-distances egales
CORRIGE ET TESTE (perimetre indique)
#1084
F30
P2
Le producteur de facteurs de vue viole le contrat du solveur S2S
NON CORRIGE / A QUALIFIER
#1086
F31
P2
max_imbalance est un bilan global, pas le maximum local
CORRIGE ET TESTE (perimetre indique)
#1080
F32
P2
Statut lineaire rapporte avant la sous-relaxation du champ retourne
NON CORRIGE / A QUALIFIER
#1088
F33
P2
Bornage haute resolution different entre assemblage et reconstruction
CORRIGE ET TESTE (perimetre indique)
#1089
F34
P1
Plusieurs tests recodent la formule fautive au lieu d'un oracle independant
PARTIEL
—
G01
P1
Politiques de frontière incohérentes entre LS et WLS
CORRIGE ET TESTE (perimetre indique)
#1092
G02
P1
CG interne : convergence dépendante de l’échelle et résidu mal normalisé
CORRIGE ET TESTE (perimetre indique)
#1092
G03
P1
Précision mixte : le vrai résidu est calculé avec l’opérateur FP32
CORRIGE ET TESTE (perimetre indique)
#1092
G04
P1
BiCGStab : seuils de rupture absolus et diagnostic nul sur échec
CORRIGE ET TESTE (perimetre indique)
#1092
G05
P2
GMRES : second membre nul refusé avec un itéré initial non nul
CORRIGE ET TESTE (perimetre indique)
#1092
G06
P1
GMRES : produits scalaires Arnoldi incompatibles avec la politique de réduction
CORRIGE ET TESTE (perimetre indique)
#1092
G07
P1
QUICK : repli linéaire sur un stencil vide
CORRIGE ET TESTE (perimetre indique)
#1092
G08
P1
QUICK : base quadratique non normalisée et rang dépendant des unités
CORRIGE ET TESTE (perimetre indique)
#1092
G09
P2
Oracle ExactSchur : NaN traité comme un second membre nul
CORRIGE ET TESTE (perimetre indique)
#1092
A01
P2
Bilan energetique global susceptible de compenser des erreurs locales
NON CORRIGE / A QUALIFIER
—
A02
P1
Compilation CPU sans MPI interrompue par une inclusion MPI non conditionnelle
CORRIGE ET TESTE (perimetre indique)
—
A03
P2
Dependance meshio non declaree pour les chemins de conversion qui l'exigent
NON CORRIGE / A QUALIFIER
—
A04
P1
Une diffusivite NaN peut etre transformee en zero dans le flux
CORRIGE ET TESTE (perimetre indique)
—
A05
P2
Qualification inachevee de familles numeriques et physiques
NON CORRIGE / A QUALIFIER
—
H01
P1
Normes instables et fausses convergences par overflow/underflow
CORRIGE ET TESTE (perimetre indique)
—
H02
P2
Newton retourne un echec apres son dernier pas pourtant convergent
CORRIGE ET TESTE (perimetre indique)
—
H03
P2
Diagonale matrix-free artificielle et coefficient NaN accepte
CORRIGE ET TESTE (perimetre indique)
—
H04
P1
Tolerance relative traitee comme absolue dans le solveur scalaire 1x1
CORRIGE ET TESTE (perimetre indique)
—
| H05 | P2 | Diagnostic energetique de frontiere calcule avec la mauvaise distance | CORRIGE ET TESTE R2 (perimetre indique) | — |
Traçabilite des anciens identifiants
Le premier rapport utilisait F01–F13 avec un sens different du rapport equations. Pour eviter les collisions, les F01–F34 du rapport equations sont conserves. La correspondance du premier rapport est :
Ancien audit technique
Registre unique
F01
F10
F02
F07
F03
F08
F04
F06
F05
F05
F06
A01, avec restriction de portee
F07
A02
F08
A03
F09
F31
F10
F32
F11
A04
F12
F34
F13
A05
6. Fiches detaillees : equations, integration, solution et preuve
F01 — Production de omega non limitee dans SST-2003
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/sst_solver.h:47-67; src/cfdx/physics/sst_production_coupling.h; src/cfdx/physics/production_turbulence_contract.h.
Portee : Composant SST; raccordement C++ SST disponible, pas de RANS autorise par la CLI production actuelle
Equation / identite attendue :
[ \widetilde P_k=\min(\rho\nu_t S^2,10\beta^*\rho k\omega),\quad P_\omega=\gamma\widetilde P_k/\nu_t. ]
Constat sur le code initial : Le code limite Pk dans k mais emploie gammarhoS^2 dans omega. Avec rho=k=omega=1, S=100 et F1=F2=1 : nu_t=0.0031, Pk_lim=0.9, source omega=5555.55556 au lieu de 161.290323 (facteur 34.4444). La TMR signale explicitement cette forme comme la coquille corrigee de la publication SST-2003.
Consequence : Modification du rapport production/destruction et donc de nu_t, du frottement, des separations et des echanges thermiques dans les chemins utilisant ce SST. Aucun facteur 34.4 sur un tutoriel complet n'est affirme.
Modification effectuee : Une production limitee unique alimente k et omega ; source omega = gammaPk_lim/nu_t. Le test SST existant qui recopiait gammarho*S^2 a ete corrige selon la definition SST-2003, et non en augmentant une tolerance.
Tests et preuves : test_audit_numerics: SST/high-strain et SST/low-strain ; equation_probes.tsv ; test_turbulence_qualification_equations.
Limites de cette correction : Correctif de fermeture locale. Pas de qualification d'une campagne RANS complete ni du raccordement production/restart des PR #1045/#1068.
F02 — La fermeture nommee komega2006 omet le limiteur de Wilcox
Etat : PARTIEL ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/turbulence_solver.h:421-441; src/cfdx/physics/turbulence_transport.h.
Portee : Composant k-omega
Equation / identite attendue :
[ \nu_t=k/\widehat\omega,\quad\widehat\omega=\max(\omega,C_{lim}\sqrt{2\overline S:\overline S/\beta^*}). ]
Constat sur le code initial : komega2006_cell_sources renvoie k/omega. komega_clim=0.875 est declare mais n'intervient pas dans cette viscosite. Pour k=omega=1 et S=100 : 1 au lieu de 0.00342857143, facteur 291.6667. La correction tridimensionnelle de beta (Pope) n'est pas representee non plus par cette signature scalaire.
Consequence : Ce composant ne correspond pas au Wilcox2006 standard. Le cas 2D annule la correction de Pope mais pas le limiteur de viscosite.
Modification effectuee : Le limiteur nu_t=k/max(omega,C_lim*S/sqrt(beta_star)) est implemente. Les coefficients diffusifs restent proportionnels a k/omega, conformement a la definition de Wilcox ; ils ne sont pas remplaces par la viscosite limitee.
Tests et preuves : test_audit_numerics: Wilcox/limited-viscosity et Wilcox/unlimited-diffusivity ; test_turbulence_qualification_equations.
Limites de cette correction : PARTIEL : la correction 3D de beta/Pope et la qualification de la variante complete ne sont pas realisees.
F03 — La branche COUPLED utilise une vitesse physique comme H/A dans le flux de pression
Etat : NON CORRIGE / A QUALIFIER ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/steady_incompressible_solver.h:871-974; src/cfdx/physics/steady_incompressible_solver.h:1590-1733; src/cfdx/physics/steady_incompressible_solver.h:2976-3042.
Portee : Solveur incompressible COUPLED; chemin principal
Equation / identite attendue :
[ \phi_f=\rho{\overline{\mathbf U+d_A\nabla p}_f\cdot\mathbf S_f-(d_A\nabla p)_f\cdot\mathbf S_f},\quad d_A=V/A_P. ]
Constat sur le code initial : make_rhie_chow_mass_flux documente et implemente un premier argument H/A. La branche separee construit HbyA. La branche COUPLED fournit U deja issu du systeme quantite de mouvement-pression. Le bloc continuite contient aussi div(interp(U)) moins le flux du gradient complet, sans le terme compensateur interp(d_A grad(p)). Avec U=1, p=x, rho=d_A=1, le flux interne renvoye est 0; avec le H/A correct egal a 2, il est 1.
Consequence : Les inconnues, l'equation de continuite stabilisee et le flux transporte n'ont pas le meme contrat. L'audit qui rappelle deux fois la meme fonction ne detecte pas cette erreur de formulation. Le contre-exemple est local; aucune campagne complete COUPLED n'est annoncee ici.
Solution a implementer, non appliquee dans le patch : Separer deux API explicites: flux_from_hbya et flux_from_physical_velocity. Pour U physique, inclure la compensation du gradient cellulaire et assembler sa contribution de maniere coherente dans le bloc D/C; ne pas corriger seulement le post-traitement.
Critere de cloture : Patch test U constant/p affine avec flux exact; checkerboard; identite bloc/flux verifiee par un oracle distinct; solution manufacturee u-p sur maillages non uniformes.
PR affichee correspondante : #1078. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
F04 — Le couplage CHT transitoire avance le temps a chaque iteration d'interface
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/cht_solver.h:98-177; src/cfdx/physics/energy_solver.h:350-390.
Portee : solve_two_region_cht, API deux regions
Equation / identite attendue :
[ \rho c_p V(T^{n+1}-T^n)/\Delta t+\sum_f q_f^{n+1}=Q V;\quad T^n\ \text{fixe pendant toutes les iterations internes}. ]
Constat sur le code initial : Les deux appels solve_energy ne recoivent pas de temperature au temps n figee. Chaque appel prend le champ courant comme nouvel ancien temps. Deux cellules adiabatiques, Tn=(300,400) K, rho=cp=k=V=A=1, dt=0.1 : Euler implicite exact (308.333333,391.666667) K; code (315.277778,384.722222) K apres deux iterations, avec converged=true.
Consequence : Erreur sur le temps physique meme avec bilan energetique total exact: la somme des temperatures reste 700. Ce defaut explique pourquoi un controle de conservation seul ne suffit pas.
Modification effectuee : Les champs physiques T1^n et T2^n sont copies avant la boucle d'interface et passes explicitement a solve_energy a chaque iteration. Les temperatures iteratives ne deviennent plus un nouvel ancien temps.
Tests et preuves : test_audit_cht_time: 6 combinaisons dt={0.01,0.1,1} et relaxation={0.5,1}, 24 controles ; integration_probes.tsv.
Limites de cette correction : Couplage transitoire a proprietes constantes, interfaces orthogonales concordantes. F28 reste ouvert.
F05 — Moyenne harmonique non ponderee sur cellules non uniformes
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/finite_volume_transport.h:307-327.
Portee : Assembleur transport partage par momentum, energie, turbulence et certains rayonnements
Equation / identite attendue :
[ D_f=\frac{A_f}{d_{Pf}/\Gamma_P+d_{Nf}/\Gamma_N}. ]
Constat sur le code initial : Le code applique Gamma_f=2 GammaP GammaN/(GammaP+GammaN), puis D=Gamma_f*A/|dPN|. Pour dPf=.125, dNf=.375, GammaP=1, GammaN=4 : D=3.2 au lieu de 4.571428571 (-30%).
Consequence : Flux faux aux interfaces a coefficients discontinus, meme sur maillage orthogonal et meme si la matrice conserve exactement les flux internes.
Modification effectuee : La conductance interne vaut A/(dP/GammaP+dN/GammaN), avec distances longitudinales projetees. L'assemblage et le flux reconstruit emploient la meme resistance en serie.
Tests et preuves : test_audit_numerics: diffusion/material-series ; les 56 tests de flux emploient un maillage etire et Gamma variable.
Limites de cette correction : Valide pour le chemin orthogonal a deux points. Ce correctif n'ajoute pas le terme non orthogonal manquant F06.
F06 — La diffusion du solveur n'utilise pas le laplacien non orthogonal corrige
Etat : NON CORRIGE / A QUALIFIER ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/finite_volume_transport.h:307-327; src/cfdx/core/numerics/laplacian.h:146-208; src/cfdx/core/solvers/scalar_diffusion.h:78-139.
Portee : Transport principal et solveur Poisson deux-points
Equation / identite attendue :
[ \Gamma\mathbf S_f\cdot\nabla\phi=\Gamma[\alpha_f(\phi_N-\phi_P)+\mathbf T_f\cdot(\nabla\phi)f],\quad \mathbf T_f=\mathbf S_f-\alpha_f\mathbf d{PN}. ]
Constat sur le code initial : Le transport emploie Gamma*A/|dPN| sans correction tangentielle. Un autre operateur, compute_laplacian corrige, existe mais ne prouve pas la justesse de cet assembleur. Le Poisson assembleur utilise aussi des distances euclidiennes, dont celle centre-face a la frontiere.
Consequence : La qualification du laplacien autonome ne se transfere pas automatiquement aux equations effectivement resolues. L'hypothese orthogonale doit etre explicite.
Solution a implementer, non appliquee dans le patch : Mutualiser le calcul de flux diffusif, sa partie implicite, sa correction explicite et son audit. Refuser les maillages hors domaine qualifie tant que la correction n'est pas integree.
Critere de cloture : Champ affine et solution manufacturee sur grilles inclinees puis familles polyedriques; verifier l'equation assemblee, pas seulement un operateur d'evaluation.
F07 — Interpolation LINEAR centrale 1/2 non exacte sur grille etiree
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/core/numerics/interpolation.h:204-213.
Portee : Interpolation commune; CENTRAL/BLENDED et plusieurs helpers
Equation / identite attendue :
[ \phi_f=(1-w_f)\phi_P+w_f\phi_N,\quad w_f=d_{Pf}/(d_{Pf}+d_{Nf})\ \text{sur grille orthogonale}. ]
Constat sur le code initial : Deux cellules [0,.25] et [.25,1], phi=x : la face x=.25 recoit .375 au lieu de .25. Le nom LINEAR ne garantit donc pas l'exactitude affine sur maillage non uniforme.
Consequence : Erreur de convection, de reconstruction et de plusieurs couplages selon le chemin selectionne.
Modification effectuee : LINEAR et BLENDED utilisent la projection du centre de face sur l'axe joignant les centres : w=((xf-xP).(xN-xP))/|xN-xP|^2. Le test a geometrie ouverte coplanaire a ete remplace par deux cubes fermes de meme largeur, avec inversion owner/neighbour.
Tests et preuves : test_audit_numerics: linear/nonuniform ; test_interpolation ; integration_probes.tsv.
Limites de cette correction : Exactitude affine longitudinale sur grille orthogonale non uniforme. Ne constitue pas une correction de skewness transversale.
F08 — Le poids du gradient de Gauss confond distance et projection
Etat : NON CORRIGE / A QUALIFIER ; priorite P1.
Localisation dans le snapshot source : src/cfdx/core/numerics/gradient.h:44-59; src/cfdx/core/numerics/gradient.h:97-104.
Portee : Gradient Gauss et corrections qui le reutilisent
Equation / identite attendue :
[ \nabla\phi_P=V_P^{-1}\sum_f\phi_f\mathbf S_{Pf};\quad w_{proj}=((\mathbf x_f-\mathbf x_P)\cdot\mathbf d)/|\mathbf d|^2. ]
Constat sur le code initial : Le poids utilise |xf-xP|/|xN-xP|. Pour P=(0,0), N=(1,0), f=(.5,.5), il vaut .7071068 plutot que .5 pour la projection longitudinale. La reconstruction requiert aussi une correction du decalage transversal.
Consequence : Perte d'exactitude affine sur faces decalees; contamination des gradients de pression, des invariants turbulents et des corrections diffuses.
Solution a implementer, non appliquee dans le patch : Definir un schema de reconstruction geometrique coherent, avec valeurs frontiere et correction de skew; tester independamment Gauss et least-squares.
Critere de cloture : MMS affine avec conditions de Dirichlet exactes sur toutes les faces; distinguer erreur de poids, de BC et de fermeture geometrique.
F09 — Le flux convectif sortant reconstruit ne correspond pas a la matrice
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/finite_volume_transport.h:447-457; src/cfdx/physics/finite_volume_transport.h:686-688; src/cfdx/physics/energy_solver.h:263-344.
Portee : Audit de transport et bilans energie
Equation / identite attendue :
[ F_f>0:\ \phi^{conv}_f=F_f\psi_P;\qquad F_f<0:\ \phi^{conv}_f=F_f\psi_b\ \text{pour upwind}. ]
Constat sur le code initial : L'assembleur utilise psiP en sortie; reconstruct_scalar_transport_flux utilise psib meme en sortie Dirichlet. Une cellule, flux entrant/sortant -1/+1, psiP=2, psib_sortie=5 : residu matriciel 0, flux sortant reconstruit 5 au lieu de 2, bilan reconstruit 3 au lieu de 0.
Consequence : Un bilan dit independant peut refuser une solution exacte de la matrice, ou diagnostiquer la mauvaise equation. Ne pas desactiver le bilan pour masquer l'ecart.
Modification effectuee : A une frontiere sortante, la convection reconstruite utilise la valeur owner ; a une frontiere entrante elle utilise la valeur imposee. Le flux diffusif conserve sa condition aux limites. FIXED_GRADIENT conserve le contrat distinct de son assembleur : convection sur la valeur extrapolee pour les deux signes de flux. L'oracle ancien a 3 au lieu de 6 a ete corrige : F=3 et phi_owner=2 donnent F*phi_owner=6.
Tests et preuves : test_audit_flux_contract, test_transport_conservation et integration_probes.tsv.
Limites de cette correction : Ne change pas le choix upwind des frontieres de l'assembleur ; aligne son diagnostic sur ce choix.
F10 — Lecture d'un gradient vide pendant la reconstruction de flux
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P0.
Localisation dans le snapshot source : src/cfdx/physics/finite_volume_transport.h:575-650.
Portee : CENTRAL, BLENDED, QUICK, QUICK_BOUNDED
Equation / identite attendue :
[ \text{Tout champ utilise dans }\nabla\psi\text{ doit etre construit avant acces.} ]
Constat sur le code initial : grad n'est cree que pour SECOND_ORDER_UPWIND/TVD; les autres branches lisent grad aux lignes 629-631. UBSan signale load of null pointer of type double. Le nouveau reproducer s'arrete avec UBSAN_OPTIONS=halt_on_error=1; le comportement sans arret depend de l'optimisation.
Consequence : Comportement indefini dans un chemin de bilan; un plantage n'est pas garanti a chaque compilation. Ce n'est pas un defaut d'equation, mais il empeche sa verification fiable.
Modification effectuee : Le gradient de reconstruction n'est lu que pour SECOND_ORDER_UPWIND/TVD, apres calcul. CENTRAL/BLENDED/QUICK utilisent leur valeur de face deja reconstruite.
Tests et preuves : Reproduction baseline avec ASan/UBSan: lecture de pointeur nul ligne 629 ; candidate: 56/56 contrats de flux sans diagnostic sanitizer.
Limites de cette correction : Couverture des 7 schemas avec bornage active/desactive sur le cas fourni ; pas une preuve d'absence de toute erreur memoire dans CFDX.
F11 — Diffusion de vitesse scalaire incomplete lorsque la viscosite varie
Etat : NON CORRIGE / A QUALIFIER ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/steady_incompressible_solver.h:1008-1051; src/cfdx/physics/steady_incompressible_solver.h:2530-2556; src/cfdx/physics/steady_incompressible_solver.h:2845-2892.
Portee : Quantite de mouvement avec mu_eff spatial, notamment raccordement SST C++
Equation / identite attendue :
[ \nabla\cdot\tau=\nabla\cdot[\mu_{eff}(\nabla\mathbf U+\nabla\mathbf U^T)]\quad(\nabla\cdot\mathbf U=0). ]
Constat sur le code initial : Chaque composante recoit seulement div(mu_eff grad(Ui)). Le terme transpose n'est pas retrouve dans le source body-grad(p). Pour U=(y,x,0), mu=1+x : div(tau)=(0,2,0), diffusion scalaire=(0,1,0). A viscosite constante et div(U)=0, le terme transpose disparait: ne pas generaliser ce defaut au laminaire constant.
Consequence : RANS avec viscosite turbulente variable peut resoudre une contrainte differente de celle du modele, meme avec une bonne formule de nu_t.
Solution a implementer, non appliquee dans le patch : Ajouter la divergence explicite/coherente du terme transpose et, si necessaire, du terme deviatorique. Documenter la pression modifiee absorbant 2rhok/3 selon la variante.
Critere de cloture : MMS vectorielle divergence nulle a mu(x) variable, comparaison tenseur/flux par face; verification de l'effort visqueux. Reference [R4].
F12 — Helper momentum: convection integree melangee a des termes volumiques
Etat : NON CORRIGE / A QUALIFIER ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/incompressible.h:19-53; src/cfdx/core/numerics/convection.h:89-102; src/cfdx/core/numerics/laplacian.h:205-208.
Portee : evaluate_momentum_terms, API auxiliaire; ne pas confondre avec l'assembleur principal
Equation / identite attendue :
[ \partial_t\mathbf U=-V^{-1}\sum_f F_f\mathbf U_f+\nu\nabla^2\mathbf U-\nabla p_{cin}. ]
Constat sur le code initial : compute_convection renvoie une somme de flux sans division par V, alors que compute_laplacian et le gradient sont normalises. Pour un flux volumique, convection integree en m^4/s^2, diffusion en m/s^2. Les champs retournes ne sont pas des termes homogennes d'une meme equation locale.
Consequence : Toute integration directe de ces trois champs produit une dependance artificielle au volume. Le resultat sur cellules unitaires peut masquer le defaut.
Solution a implementer, non appliquee dans le patch : Choisir explicitement des termes tous integres, ou tous par unite de volume; preciser si p est cinematique ou dynamique.
Critere de cloture : Meme champ manufacture sur domaine redimensionne et sur cellules non unitaires; controle des dimensions.
F13 — Helper thermique: convection en W et conduction en W/kg
Etat : NON CORRIGE / A QUALIFIER ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/thermal.h:19-40; src/cfdx/core/numerics/convection.h:89-102.
Portee : evaluate_energy_terms, API auxiliaire
Equation / identite attendue :
[ \partial_t T=-(\rho V)^{-1}\sum_f\dot m_f T_f+k/(\rho c_p)\nabla^2T+Q/(\rho c_p). ]
Constat sur le code initial : La convection renvoyee est cpsum(mdotT), donc une puissance W. La conduction est (k/rho)*laplacien(T), donc W/kg. Ni un residu integre d'energie ni un taux de temperature coherent n'est obtenu en les additionnant.
Consequence : Une integration utilisant ce helper peut diverger de solve_energy, dont l'assemblage a cp constant est coherent.
Solution a implementer, non appliquee dans le patch : Decider si l'API renvoie W, W/m3, W/kg ou K/s; appliquer les facteurs rho, cp et V correspondants aux deux termes.
Critere de cloture : Diffusion thermique et convection identiques via le helper et l'assembleur energie, avec rho, cp et V non unitaires.
F14 — Integrale d'enthalpie tabulee incorrecte hors table
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/thermophysical_models.h:114-191.
Portee : ThermophysicalProperties::enthalpy
Equation / identite attendue :
[ h(T,T_r)=\int_{T_r}^{T}c_p(\theta)d\theta;\quad h(T,T)=0. ]
Constat sur le code initial : Cp=1 tabule sur [300,400], extrapolation CLAMP : h(200,100)=200 au lieu de 100; h(600,500)=200 au lieu de 100; h(200,200)=100 au lieu de 0. Les queues d'integration ne sont pas intersectees avec l'intervalle demande.
Consequence : Energie sensible incorrecte dans les chemins utilisant cette primitive, independamment du solveur lineaire.
Modification effectuee : L'integrale est decoupee sur l'intervalle effectivement demande, les noeuds de table et les changements de bornage. CLAMP et extrapolation LINEAR sont integres sur leurs bons sous-intervalles ; antisymetrie et intervalle nul sont respectes.
Tests et preuves : test_audit_numerics: enthalpy/clamp, linear, zero, antisymmetry, additivity ; equation_probes.tsv.
Limites de cette correction : Tables lineaires traitees exactement ; autres lois bornees via quadrature adaptative controlee, non une campagne thermophysique exhaustive.
F15 — L'enthalpie n'integre pas le Cp effectivement borne
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/thermophysical_models.h:63-65; src/cfdx/physics/thermophysical_models.h:114-135.
Portee : ThermophysicalProperties avec enforce_bounds
Equation / identite attendue :
[ \frac{\partial h}{\partial T}=c_p^{utilise}(T). ]
Constat sur le code initial : Cp de reference 1000 borne a 500 : cp(350)=500 mais h(400,300)=100000 au lieu de 50000 J/kg. L'integrale analytique ne tient pas compte des bornes actives.
Consequence : Etat calorique et equation de temperature incompatibles; un Newton thermodynamique pourrait employer une derivee differente de la fonction integree.
Modification effectuee : L'integrande est le Cp effectivement retourne, apres bornage. Constant/lineaire/tabule : integration exacte par morceaux ; polynome non borne : primitive analytique ; autres cas : Simpson adaptatif avec rejet si la precision ne peut etre atteinte.
Tests et preuves : test_audit_numerics: clipped-linear, clipped-polynomial, clipped-constant et NaN-interval-count.
Limites de cette correction : La coherence dh/dT=Cp est restauree pour ces chemins. Ne qualifie pas une equation d'energie enthalpique a Cp variable dans tous les solveurs.
F16 — Cp par defaut et conductivite de l'air errones dans une surcharge transport
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/transport_models.h:73-88; src/cfdx/physics/transport_models.h:125-135.
Portee : compute_transport(T,rho), pas la surcharge avec EOS explicite
Equation / identite attendue :
[ c_p=\gamma R/(\gamma-1),\qquad k=\mu c_p/Pr. ]
Constat sur le code initial : Le defaut vaut AIR_R*1005/1.4. Pour R du code et gamma=1.4 : Cp=206059.091737 au lieu de 1004.666218 J/kg/K. A 300 K, k=5.357293326 au lieu de .02612013661 W/m/K. La surcharge EOS explicite suit une autre voie.
Consequence : Diffusivite/energie et nombre de Prandtl effectif incoherents pour les utilisateurs de cette surcharge. Ce n'est pas une mesure dans tous les cas de production.
Modification effectuee : La surcharge air utilise Cp=gamma*R/(gamma-1) avec gamma=1.4, et non une expression multiplier erronee.
Tests et preuves : test_audit_numerics: transport/air-Cp ; equation_probes.tsv.
Limites de cette correction : Surcharge identifiee uniquement ; pas de changement des proprietes des autres gaz.
F17 — Deux IdealGasEOS renvoient des derivees incompatibles
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P2.
Localisation dans le snapshot source : src/cfdx/thermodynamics/equation_of_state.h:190-192; src/cfdx/physics/equation_of_state.h:249-250.
Portee : Ancienne/reorganisee API thermodynamics vs API physics
Equation / identite attendue :
[ (\partial p/\partial T)_\rho=\rho R=p/T. ]
Constat sur le code initial : La version thermodynamics evalue density(0,T)*R et renvoie 0. La version physics renvoie 337.75 a p=101325 Pa,T=300 K, valeur correcte. Aucun appel actif de la premiere n'a ete etabli dans le pilote principal.
Consequence : Mauvaise Jacobienne thermodynamique si ce doublon public est utilise; pas de preuve que le solveur incompressible actuel l'emploie.
Modification effectuee : IdealGasEOS du module thermodynamics calcule dp/dT|rho = density(p,T)*R, avec la pression fournie, et non density(0,T).
Tests et preuves : test_audit_numerics: EOS/dp-dT ; equation_probes.tsv.
Limites de cette correction : Ne corrige pas le modele IncompressibleEOS de F18.
F18 — Identites incompatibles dans IncompressibleEOS pour beta non nul
Etat : NON CORRIGE / A QUALIFIER ; priorite P2.
Localisation dans le snapshot source : src/cfdx/physics/equation_of_state.h:110-150; src/cfdx/thermodynamics/equation_of_state.h.
Portee : API EOS incompressible a beta optionnel
Equation / identite attendue :
[ a^2=(\partial p/\partial\rho)_s,\qquad (\partial\rho/\partial T)_p=d\rho(T)/dT. ]
Constat sur le code initial : Pour rho=1000,beta=.001, a^2=1e6 mais dp_drho_s=1. La densite retournee est constante tandis que drho_dT_p=-1 et sa difference finie vaut 0. Le meme parametre beta melange des interpretations qui ne sont pas un contrat EOS coherent.
Consequence : Une precondition ou une Jacobienne basee sur ces derivees ne correspond pas a la loi d'etat. En incompressible strict, le son fini doit etre un artifice numerique explicitement distinct.
Solution a implementer, non appliquee dans le patch : Separer dilatation thermique beta_T, compressibilite kappa et compressibilite artificielle; imposer les identites de derivees.
Critere de cloture : Differences finies de rho(p,T), identite acoustique et verification dimensionnelle de chaque parametre.
F19 — BDF2 a pas variable incorrect dans le pilote generique double temps
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/dual_time_field_integrator.h:177-185; src/cfdx/physics/dual_time_stepping.h:100-117; src/cfdx/core/numerics/temporal.h; src/cfdx/physics/dual_time_navier_stokes.h:113-125.
Portee : DualTimeFieldIntegrator; le noyau BDF2 et le chemin NS possedent la bonne variante
Equation / identite attendue :
[ \dot u^{n+1}=\frac{1}{h}[\frac{1+2r}{1+r}u^{n+1}-(1+r)u^n+\frac{r^2}{1+r}u^{n-1}],\quad r=h/h_{prec}. ]
Constat sur le code initial : Le pilote conserve dt_previous mais utilise le residu BDF2 a coefficients constants 1.5,-2,.5. Pour u'=1, u0=0, pas .1 puis .2 : .2666666667 au lieu de .3. Le meme cas via core::TimeIntegrationContext donne .3.
Consequence : L'adaptation du pas modifie la solution meme pour une solution lineaire en temps. Ne pas accuser toutes les implementations BDF2: le chemin NS lit le pas precedent.
Modification effectuee : Le dt physique precedent est transmis du stockage d'historique au pilote double temps. Les coefficients BDF2 sont a0=(1+2r)/(1+r), a1=-(1+r), a2=r^2/(1+r), r=dt/dt_precedent ; demarrage BE conserve.
Tests et preuves : test_audit_numerics: dual-time/variable-step sur 0.1,0.2,0.05,0.08 ; equation_probes.tsv.
Limites de cette correction : Pilote generique uniquement ; aucune modification du noyau BDF2 qui etait deja correct.
F20 — Une iteration implicite peut accepter un NaN comme convergence
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/core/numerics/temporal.h:146-200.
Portee : core::advance_time, iterations de Picard
Equation / identite attendue :
[ |u^{k+1}-u^k|_\infty<tol\quad\text{et tous les champs/residus sont finis}. ]
Constat sur le code initial : Avec u'=-2 sqrt(u), u0=1, dt=1, le premier iterant devient negatif puis NaN. La reduction std::max avec NaN conserve une valeur nulle et le resultat NaN est renvoye sans exception. L'equation implicite possede pourtant la racine positive (sqrt(2)-1)^2.
Consequence : Fausse convergence. La convergence d'un point fixe n'est par ailleurs pas garantie par la stabilite de la discretisation Euler implicite.
Modification effectuee : Les valeurs nouvelles et anciennes de l'iteration implicite doivent etre finies ; NaN/Inf declenchent une exception plutot qu'un changement maximal artificiellement nul.
Tests et preuves : test_audit_numerics: implicit/reject-nonfinite ; equation_probes.tsv.
Limites de cette correction : Le cas non contractant est rejete, pas resolu miraculeusement. Aucun NaN retourne n'est compte comme solution numerique correcte.
F21 — Les integrateurs low-storage ne mettent a jour que la composante 0
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/low_storage_time_integration.h:9-73.
Portee : RK2, SSP-RK3 et RK3 Williamson pour Field vectoriel
Equation / identite attendue :
[ \mathbf u^{n+1}=\operatorname{RK}(\mathbf u^n,\mathbf f)\quad\text{sur chaque composante}. ]
Constat sur le code initial : Avec u=(1,2,3), RHS=(1,1,1), dt=.1, RK2 renvoie (1.1,2,3) au lieu de (1.1,2.1,3.1). Les temporaires sont scalaires et les boucles utilisent component_data(0). La signature n'interdit pas les champs vectoriels.
Consequence : Vitesse vectorielle partiellement integree si l'API est employee directement. Les coefficients du RK scalaire ne sont pas en cause.
Modification effectuee : Les buffers low-storage ont la dimension du champ et les 3 schemas mettent a jour toutes les composantes. La forme et la finitude de chaque RHS sont verifiees.
Tests et preuves : test_audit_numerics: RK/all-components et ordre sur u'=u, 10/20/40 pas, RK2/WilliamsonRK3/SSPRK3.
Limites de cette correction : Les ordres temporels testes ne constituent pas une preuve de stabilite pour toute discretisation spatiale raide.
F22 — Une fonction nommee WALE calcule la loi de Smagorinsky
Etat : NON CORRIGE / A QUALIFIER ; priorite P2.
Localisation dans le snapshot source : src/cfdx/physics/turbulence_models.h:158-164; src/cfdx/physics/turbulence_models.h:411-416.
Portee : Helper wale_nu_t; une autre implementation tensorielle correcte existe
Equation / identite attendue :
[ \nu_t=(C_w\Delta)^2\frac{(S^d:S^d)^{3/2}}{(S:S)^{5/2}+(S^d:S^d)^{5/4}}. ]
Constat sur le code initial : wale_nu_t(Cw,delta,S) renvoie (Cwdelta)^2S. En cisaillement simple pur, grad(U)^2=0 donc Sd=0 et WALE doit donner 0; le helper donne .025 pour Cw=.5,delta=.1,S=10. wale_eddy_viscosity avec invariants donne bien 0.
Consequence : Deux contrats physiques differents sous le meme nom. L'implementation tensorielle ne doit pas etre corrigee par remplacement par le helper fautif.
Solution a implementer, non appliquee dans le patch : Supprimer/renommer le helper a trois scalaires; WALE requiert les invariants du tenseur grad(U). Unifier les tests et la selection du modele.
Critere de cloture : Cisaillement simple, rotation solide, invariance par rotation de repere, limite au mur. Reference [R3].
F23 — Correction de compressibilite turbulente non adimensionnelle
Etat : NON CORRIGE / A QUALIFIER ; priorite P2.
Localisation dans le snapshot source : src/cfdx/physics/turbulence_models.h:81-88.
Portee : Helpers de correction turbulence
Equation / identite attendue :
[ M_t=\sqrt{2k}/a\quad\text{pour une correction basee sur le Mach turbulent}. ]
Constat sur le code initial : Le code forme mach=strain_rate/sound_speed. S est en s^-1 et a en m/s : leur rapport est en m^-1, pas un Mach. Une longueur manquante ne peut pas etre implicitement egale a 1 pour des maillages generaux.
Consequence : Correction dependante des unites et de l'echelle geometrique. Cela concerne le helper, pas une preuve d'un solveur compressible complet actif.
Solution a implementer, non appliquee dans le patch : Choisir une correction publiee et fournir ses variables physiques (k, a, eventuellement dilatation).
Critere de cloture : Changer d'unites m->mm sans modifier la physique; le facteur adimensionnel doit rester identique.
F24 — Correction Kato-Launder ajoutee avec des dimensions incompatibles
Etat : NON CORRIGE / A QUALIFIER ; priorite P2.
Localisation dans le snapshot source : src/cfdx/physics/turbulence_models.h:101-108; src/cfdx/physics/turbulence_models.h:128-147.
Portee : turbulence_models::production corrigee
Equation / identite attendue :
[ P_{KL}=\mu_t S\Omega\quad\text{(volumique), ou }\nu_t S\Omega\text{ (par masse)}. ]
Constat sur le code initial : Le helper renvoie CSOmega puis ce terme est ajoute a une production portant nu_t/rho selon le contrat. Le facteur de viscosite manque et la correction est additive, alors que la forme Kato-Launder remplace le modele de production. QCR, dans ce meme ensemble, est un facteur scalaire et non une relation de contrainte tensorielle complete.
Consequence : Les noms de corrections publiees recouvrent des approximations algebraiques non equivalentes.
Solution a implementer, non appliquee dans le patch : Implementer la formule et les tenseurs de la variante choisie; sinon nommer ces fonctions heuristiques et ne pas les qualifier sous le nom du modele.
Critere de cloture : Controle des unites, cas S=0, Omega=0, invariance du tenseur; specification exacte de la variante. Reference [R1], section Kato-Launder.
F25 — Surcharge historique SST sans invariant de deformation
Etat : NON CORRIGE / A QUALIFIER ; priorite P2.
Localisation dans le snapshot source : src/cfdx/physics/turbulence.h:31-40.
Portee : turbulent_kinematic_viscosity_komega_sst a quatre arguments
Equation / identite attendue :
[ \nu_t=a_1 k/\max(a_1\omega,S F_2). ]
Constat sur le code initial : La surcharge historique remplace S par omega, tandis que la surcharge a cinq arguments utilise S. Pour F2=1, elle calcule a1*k/omega independamment du cisaillement.
Consequence : API conservant un nom SST mais une formule non generale. Les appelants historiques peuvent obtenir des viscosites differentes sans changement apparent du modele.
Solution a implementer, non appliquee dans le patch : Deprecier cette surcharge; exiger S ou un objet d'invariants. Ne pas assimiler compatibilite source a compatibilite physique.
Critere de cloture : Tester egalite des surcharges uniquement dans le cas particulier S=omega, et rejet/documentation pour l'ancienne API.
PR affichee correspondante : #1081. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
F26 — Facteur pi manquant dans une source P1 auxiliaire
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P2.
Localisation dans le snapshot source : src/cfdx/physics/m1_m4_models.h:276-281; src/cfdx/physics/radiation.h:135-145; src/cfdx/physics/radiation_solver.h:58-119.
Portee : m1m4::p1_source; le P1 principal utilise la bonne convention
Equation / identite attendue :
[ \nabla\cdot\mathbf q_r=4\pi\kappa(I_b-J)=\kappa(4\sigma T^4-G),\quad I_b=\sigma T^4/\pi,\ G=4\pi J. ]
Constat sur le code initial : La fonction m1m4 utilise 4kappa(Ib-J), tandis que p1_radiative_source utilise 4pikappa*(Ib-J). T=1000 K,kappa=1,J=0 : 72197.4494 au lieu de 226814.9768 W/m3, facteur pi.
Consequence : Mauvaise normalisation angulaire si le helper est utilise. Les equations et signes du solveur P1 principal ne presentent pas ce defaut dans le chemin examine.
Modification effectuee : Le helper P1 utilise 4pikappa*(Ib-J), compatible avec Ib=sigma*T^4/pi. Le solveur P1 principal, deja coherent, n'a pas ete modifie.
Tests et preuves : test_audit_numerics: P1/full-solid-angle ; equation_probes.tsv.
Limites de cette correction : Correction auxiliaire distincte de la PR #1079 : sa tete n'a pas pu etre telechargee ni compilee.
PR affichee correspondante : #1079. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
F27 — Courbure VOF: derivees de composantes sur le meme axe
Etat : NON CORRIGE / A QUALIFIER ; priorite P2.
Localisation dans le snapshot source : src/cfdx/physics/vof.h:67-77.
Portee : interface_curvature_from_normals; prototype VOF
Equation / identite attendue :
[ \kappa=-\nabla\cdot\mathbf n=-(\partial_x n_x+\partial_y n_y+\partial_z n_z). ]
Constat sur le code initial : La routine differencie nx et ny avec les memes voisins i-1,i+1 et le meme pas. Pour n=(cos(x),sin(x)), independant de y, elle donne environ -1 en x=0 alors que -div(n)=0.
Consequence : Une force sigmakappagrad(alpha) construite a partir de ce resultat peut produire une courbure et des courants parasites artificiels.
Solution a implementer, non appliquee dans le patch : Exiger un stencil multidimensionnel/maillage et calculer chaque derivee dans sa propre direction; ou renommer la fonction comme approximation strictement 1D definie.
Critere de cloture : Plan (kappa=0), cercle et sphere avec signe fixe, changement de repere, test de Laplace pression sigma*kappa.
PR affichee correspondante : #1083. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
F28 — Interface CHT : distances et conductivités de deux contrats différents
État R2 : CORRIGÉ ET TESTÉ dans le périmètre ci-dessous ; priorité P2.
Localisation initiale : src/cfdx/physics/cht_solver.h:125-145 et src/cfdx/physics/finite_volume_transport.h:197-219. Le patch cumulatif donne les positions après modification.
Contrat discret attendu : pour une interface conductrice, conforme et imperméable,
[ d_{n,i}=|(\mathbf x_f-\mathbf x_i)\cdot\mathbf n_i|,\qquad D_i=\frac{k_i A_i}{d_{n,i}},\qquad T_I=\frac{D_1T_1+D_2T_2}{D_1+D_2}. ]
La conservation porte sur les puissances opposées : [ Q_{1\to I}=D_1(T_1-T_I),\qquad Q_{I\to2}=D_2(T_I-T_2). ]
Défaut avant R2 : l'interface utilisait les distances euclidiennes et les conductivités de ChtInterfaceControls, alors que les régions assemblaient avec les distances normales et les matériaux de EnergySolverControls. Le diagnostic d'interface réutilisait le premier contrat et pouvait donc converger sans conserver l'énergie des matrices.
Correction appliquée : les conductances sont calculées une fois avec les matériaux des équations d'énergie et leur fonction boundary_normal_distance. Les anciens paramètres conductivity1/2 sont conservés comme assertions de compatibilité, pas comme une deuxième source de matériaux. Un désaccord est rejeté avant modification des champs. La pondération est mise à l'échelle pour éviter de sommer directement deux grandes conductances. Les diagnostics comparent les puissances effectivement utilisées, en W.
Le code rejette aussi : normales non opposées, interfaces vides ou dupliquées, formes non conformes (comparaison des sommets), géométrie non finie, distances normales dégénérées, flux de masse non nul à l'interface, champs non scalaires/non finis et pas physiques différents entre les régions. Le couplage ne prétend donc pas gérer un contact perméable, non conforme ou sous-cyclé.
Preuves : test_audit_cht_contract comprend 13 configurations admissibles (65 contrôles) et 18 configurations invalides (36 contrôles : rejet et absence de mutation). Avant R2 : 71 échecs/101. Après : 101/101. Les cas varient les longueurs, cisaillements géométriques, pas de temps, relaxations et conductivités ; le contraste maximal testé vaut 1000. L'erreur thermique maximale passe de 5,7686548084 K à 5,8094×10^-11 K. Le bilan énergétique normalisé maximal passe de 0,8856282323 à 1,7508×10^-11.
Limites : l'oracle est le système discret indépendant de deux capacités, pas la solution continue générale d'un domaine incliné. Les essais de solveur utilisent une cellule de chaque côté et une interface commune. Les contrôles de conformité ne remplacent pas une campagne multi-face/multi-région, ni la correction F06 de non-orthogonalité dans le transport.
F29 — Resistance thermique de zone suppose des demi-distances egales
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/thermal_zones.h:236-240.
Portee : thermal_zones, contact sur maillage multi-zone
Equation / identite attendue :
[ G=A/(d_{Pf}/k_P+R_c+d_{Nf}/k_N). ]
Constat sur le code initial : La resistance utilise distance/(2ko)+Rc+distance/(2kn). Cela suppose une face a mi-chemin. La routine cht_interface_conductance de thermal.h accepte, elle, deux distances et porte la bonne structure de resistances.
Consequence : Le meme probleme de geometrie que F05 se retrouve dans le contact thermique entre materiaux.
Modification effectuee : Resistance en serie A/(dP/kP+Rc+dN/kN), distances normales separees au lieu de la demi-distance moyenne.
Tests et preuves : test_audit_numerics: thermal/contact-series ; test_thermal_zones_vv.
Limites de cette correction : Geometrie orthogonale testee ; ne generalise pas le traitement des interfaces CHT F28.
PR affichee correspondante : #1084. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
F30 — Le producteur de facteurs de vue viole le contrat du solveur S2S
Etat : NON CORRIGE / A QUALIFIER ; priorite P2.
Localisation dans le snapshot source : src/cfdx/physics/radiation_advanced.h:235-289; src/cfdx/physics/radiation_s2s.h:48-85.
Portee : estimate_view_factor_matrix -> validate_s2s_view_factors
Equation / identite attendue :
[ A_iF_{ij}=A_jF_{ji},\quad\sum_jF_{ij}\le1;\quad F_{ii}=0\ \text{pour une face plane}. ]
Constat sur le code initial : L'estimateur force la fermeture en ajoutant le deficit sur Fii; le solveur S2S rejette tout Fii non nul. Deux surfaces planes paralleles unite distantes de 2 : Fii=.9204225285 et sortie de l'estimateur rejetee par le consommateur.
Consequence : Raccordement impossible et auto-vue artificielle. La conservation algebraique n'est pas une preuve de justesse geometrique des facteurs de vue.
Solution a implementer, non appliquee dans le patch : Garder une fraction ouverte vers l'ambiance, ou utiliser une vraie enceinte geometrique et une integration convergente. Aligner les hypotheses convexite/auto-vue du producteur et consommateur.
Critere de cloture : Plans paralleles, enceinte fermee, surfaces d'aires differentes, occultation, convergence statistique pour ray tracing.
PR affichee correspondante : #1086. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
F31 — max_imbalance est un bilan global, pas le maximum local
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P2.
Localisation dans le snapshot source : src/cfdx/physics/finite_volume_transport.h:535-538.
Portee : Diagnostic ScalarEquation
Equation / identite attendue :
[ r_{global}=|\sum_c r_c|,\quad r_{local}=\max_c|r_c|. ]
Constat sur le code initial : Le champ max_imbalance est abs(sum(divergence)). Pour deux cellules +10 et -10, le resultat est 0 et non 10. Les controles locaux recents du solveur energie existent: ce constat ne signifie pas que tous les solveurs ne verifient qu'un total.
Consequence : Nom trompeur et compensation si ce diagnostic est utilise comme preuve locale.
Modification effectuee : max_imbalance est desormais le maximum absolu des sommes locales de flux. Deux cellules +1/-1 renvoient 1, pas 0.
Tests et preuves : test_audit_numerics: transport/local-imbalance.
Limites de cette correction : La PR #1080 annonce en plus la separation des metriques globale et locale. Notre patch corrige la metrique existante, sans inventer une qualification du diff de cette PR.
PR affichee correspondante : #1080. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
F32 — Statut lineaire rapporte avant la sous-relaxation du champ retourne
Etat : NON CORRIGE / A QUALIFIER ; priorite P2.
Localisation dans le snapshot source : src/cfdx/physics/finite_volume_transport.h:885-898.
Portee : solve_scalar_equation, API generique
Equation / identite attendue :
[ x_{ret}=x_{old}+\alpha(x_{solve}-x_{old});\quad r_{ret}=b-Ax_{ret}. ]
Constat sur le code initial : Le SolverResult est celui du candidat resolu; ensuite x est sous-relaxe. Le residu renvoye ne mesure donc pas necessairement x retourne. Les pilotes qui recalculent le residu sur l'etat accepte evitent cette confusion.
Consequence : Une API appelee directement peut annoncer convergence sur un champ encore eloigne du systeme. Ce n'est pas une preuve que tous les pilotes font cette erreur.
Solution a implementer, non appliquee dans le patch : Renvoyer deux statuts/residus nommes: solution lineaire candidate et etat accepte; recalculer le second.
Critere de cloture : Systeme scalaire A=1,b=1,xold=0,alpha=.1; residu de xret=.1 doit valoir .9.
PR affichee correspondante : #1088. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
F33 — Bornage haute resolution different entre assemblage et reconstruction
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P2.
Localisation dans le snapshot source : src/cfdx/physics/finite_volume_transport.h:339-407; src/cfdx/physics/finite_volume_transport.h:629-650.
Portee : SECOND_ORDER_UPWIND/TVD avec bounded_convection=false
Equation / identite attendue :
[ F_f^{assemble}(\psi)=F_f^{audit}(\psi)\quad\text{pour un meme schema et les memes controles}. ]
Constat sur le code initial : L'assemblage borne la valeur haute resolution entre owner/neighbour, alors que la reconstruction conditionne ce bornage a bounded_convection. Le comportement differe quand cette option est fausse et qu'une extrapolation depasse ces bornes.
Consequence : Bilan reconstruit et equation assemblee potentiellement differents, independamment du defaut Dirichlet F09.
Modification effectuee : Les extrapolations SOU/TVD de l'assembleur ne sont bornees que lorsque bounded_convection est actif, comme la reconstruction. QUICK_BOUNDED garde son contrat propre.
Tests et preuves : test_audit_flux_contract: 7 schemas x 2 modes x 2 familles de frontieres x 2 signes de flux ; accord flux/matrice jusqu'a 1.42e-14.
Limites de cette correction : Ce test verifie l'identite discrete. Il ne prouve pas a lui seul la bornitude physique d'une solution transitoire.
PR affichee correspondante : #1089. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
F34 — Plusieurs tests recodent la formule fautive au lieu d'un oracle independant
Etat : PARTIEL ; priorite P1.
Localisation dans le snapshot source : tests/unit/test_turbulence_qualification_equations.cpp:53-85; tests/unit/test_thermophysical_models.cpp; tests/validation/test_cht_validation.cpp; tests/validation/test_n9_s9_rhie_chow_consistency.cpp.
Portee : Qualification mathematique
Equation / identite attendue :
[ \text{Test de verification}\ne\text{copie de la fonction attendue}. ]
Constat sur le code initial : Les attentes komega2006 reprennent k/omega; celles de SST reprennent gammarhoS^2. La campagne existante ciblee passe 19/19 alors que les sondes independantes montrent 28 ecarts sur 38 sorties (certaines sorties concernent le meme defaut). Le test CHT execute ne couvre pas la contre-preuve temporelle nouvelle.
Consequence : Un tableau de tests verts ne qualifie ni le modele publie ni son raccordement. Ne pas convertir les 28 sorties en 28 bugs distincts.
Modification effectuee : Ajout de 3 executables de regression a verifications actives en Release. Correction explicite de quatre fichiers de tests : SST/Wilcox, campagne WLS, geometrie d'interpolation et flux upwind de sortie. Les anciens journaux en echec sont conserves, y compris une regression Neumann de notre premiere version : 28/56 cas en echec, puis 56/56 apres restauration du contrat Neumann. Aucun echec de cette version intermediaire n'est dissimule.
Tests et preuves : 167 verifications nouvelles dans 3 tests ; anciens reproducteurs et campagne CTest.
Limites de cette correction : PARTIEL : les 18 entrees historiques non qualifiees et les campagnes completes ne sont pas automatiquement closes.
G01 — Politiques de frontière incohérentes entre LS et WLS
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/core/numerics/gradient.h:139-195; src/cfdx/core/numerics/gradient.h:802-901.
Portee : API de reconstruction. WLS ignore EXCLUDE_BOUNDARY ; LS ignore la demande ZERO_GRADIENT_GHOST. Ne prouve pas que toutes les reconstructions de production passent par WLS.
Pour les moindres carrés, les échantillons doivent correspondre à la politique demandée :
[ \min_{\mathbf g_P}\sum_{N\in\mathcal S_P}w_{PN} \left(\mathbf g_P\cdot(\mathbf x_N-\mathbf x_P)-(\phi_N-\phi_P)\right)^2. ]
Dans WLS, la branche frontière teste REJECT, mais pas EXCLUDE_BOUNDARY. Elle ajoute des cellules fantômes de gradient normal nul en l'absence de conditions explicites. Dans LS, une frontière non rejetée est au contraire toujours ignorée, même lorsque ZERO_GRADIENT_GHOST est demandé.
Contre-exemple : cube cartésien, champ exact (\phi=x), gradient exact ((1,0,0)). WLS + EXCLUDE_BOUNDARY donne (g_x=0,5) dans les cellules des faces x, au lieu de 1. Pour les maillages (3^3,6^3,12^3), l'erreur maximale reste 0,5. Les erreurs L2 globales sont respectivement 0,4082482905 ; 0,2886751346 ; 0,2041241452. Cette famille présente donc un ordre global observé de 0,5 et aucune convergence en norme maximale.
Contrôles positifs : LS avec exclusion donne l'exactitude affine ; WLS avec la politique fantôme et les bonnes valeurs de Dirichlet donne aussi l'exactitude affine. La réflexion géométrique Neumann/Dirichlet n'est pas mise en cause ici. Le défaut est le non-respect de l'option choisie, pas le fait qu'une condition de Neumann nulle soit incompatible avec (\phi=x).
Test de non-régression : champ affine avec chaque politique, BC explicites, cellules frontalières incluses dans les normes ; contrôler aussi que REJECT rejette effectivement.
Modification effectuee : WLS ignore effectivement les faces frontieres sous EXCLUDE_BOUNDARY ; LS construit les fantomes quand ZERO_GRADIENT_GHOST est demande. REJECT reste un rejet explicite.
Tests et preuves : gradient_repro.cpp sur cubes 3^3,6^3,12^3 ; test_audit_numerics: WLS/exclude/all-cells, LS/ghost/corner.
Limites de cette correction : La reconstruction d'un stencil 3D de rang deficient reste une limitation distincte ; pas de preuve d'invariance generale en rotation.
PR affichee correspondante : #1092. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
G02 — CG interne : convergence dépendante de l’échelle et résidu mal normalisé
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/core/linalg/cg_solver.h:154-178; src/cfdx/core/linalg/cg_solver.h:245-269; src/cfdx/core/solvers/scalar_diffusion.h:336-338; src/cfdx/physics/finite_volume_transport.h:780-785.
Portee : Chemin CG sans objet préconditionneur fourni, utilisant Jacobi interne. Appel direct par le solveur scalaire de diffusion et dernier repli du transport. Le chemin avec objet préconditionneur externe a un contrôle FP64 distinct.
Le chemin interne compare essentiellement (\sqrt{r^TM^{-1}r}) à (\mathrm{tol},|b|_2). Il mélange une norme naturelle préconditionnée avec une référence non préconditionnée.
Si (A,b,M) sont multipliés par (s>0), la solution ne change pas, mais le membre gauche croît comme (\sqrt{s}) et le seuil comme (s). Le critère devient donc artificiellement plus facile quand les coefficients augmentent.
Cas matriciel : (A=s\left[\begin{smallmatrix}4&1\1&3\end{smallmatrix}\right]), (b=s(6,7)^T), (x_*=(1,2)^T), tolérance (10^{-10}). À (s=10^{24}), le solveur retourne CONVERGED en zéro itération avec (x=0) et un vrai résidu relatif de 1.
Cas intégré : appel réel à solve_laplace_dirichlet, cube (4^3), BC (\phi=x), diffusivité constante (\Gamma), source nulle, tolérance (10^{-10}). La solution exacte (\phi=x) ne dépend pas de (\Gamma).
Gamma
Itérations
Résidu relatif annoncé
Vrai résidu relatif FP64
Erreur maximale du champ
1
10
4,1253e-18
2,6735e-16
2,2204e-16
1e8
9
5,4884e-11
7,7635e-7
3,8437e-7
1e12
8
1,1016e-11
1,4638e-5
1,1343e-5
1e24
0
1
1
0,875
Tous ces retours ont le statut CONVERGED. Ce sont des tests d'invariance à l'échelle, pas l'affirmation qu'un tutoriel utilise une diffusivité de 1e24.
Le calcul du vrai résidu effectué dans le reproducteur utilise directement la matrice et le second membre renvoyés par le solveur. Le défaut n'est donc pas limité à une reconstruction analytique externe de l'équation.
Précision de portée : linear_solver_dispatch.h crée un objet Identity/Jacobi, même pour certaines sélections sans préconditionnement. Ce chemin passe dans la branche externe et ne doit pas être assimilé au cas interne défectueux.
Modification effectuee : Toutes les acceptations CG reposent sur le vrai residu FP64 non preconditionne. La normalisation n'emploie plus un plancher absolu 1e-15 pour ||b||. Les valeurs non finies sont rejetees ; precision mixte peut redemarrer en FP64.
Tests et preuves : Laplace 4^3 a Gamma=1,1e8,1e12,1e24 ; CG scalaire a echelles 1e-100/1e-200 ; test_audit_numerics.
Limites de cette correction : Les changements de norme n'eliminent pas tout probleme de conditionnement ou d'overflow dans les produits scalaires a des echelles extremes.
PR affichee correspondante : #1092. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
G03 — Précision mixte : le vrai résidu est calculé avec l’opérateur FP32
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/core/linalg/gmres_solver.h:82-85; src/cfdx/core/linalg/gmres_solver.h:383-400; src/cfdx/core/linalg/gmres_solver.h:438-453; src/cfdx/core/linalg/bicgstab_solver.h:68-75; src/cfdx/core/linalg/bicgstab_solver.h:97-105; src/cfdx/core/linalg/mixed_precision.h:36-40.
Portee : Options de précision mixte GMRES, FGMRES et BiCGStab ; option désactivée par défaut. Ne constitue pas une preuve que la CLI de production active FP32.
La politique de précision promet un contrôle du vrai résidu FP64. Or la fonction nommée true_residual réutilise l'application de matrice FP32. Calculer ensuite la norme en double précision ne restaure pas l'information perdue dans le produit matrice-vecteur.
Contre-exemple : ((1+10^{-8})x=1), (x_0=0), tolérance (10^{-12}). Les trois solveurs retournent CONVERGED avec (x=1) et un résidu annoncé nul. Le résidu recalculé avec la matrice d'origine est (9,99999994\times10^{-9}), environ 10 000 fois la tolérance.
Contrôles positifs : en FP64, les trois solveurs renvoient (x=0,99999999000000017) et un résidu nul à la précision mesurée. Le mauvais résultat est donc bien lié au raccordement de la précision mixte.
Pour une API matrix-free, prévoir une application haute précision distincte ou annoncer explicitement qu'une certification par rapport à l'opérateur original n'est pas disponible. Ne pas prétendre qu'une tolérance arbitraire sera atteinte par un seul calcul FP32.
Modification effectuee : L'operateur de travail FP32 et l'operateur de certification FP64 sont distincts. GMRES/FGMRES reviennent en FP64 avec le budget restant si la precision mixte ne satisfait pas le vrai residu ; CG remplace son residu et redemarre en FP64 dans les cas controles.
Tests et preuves : test_audit_numerics: 4 solveurs sur (1+1e-8)x=1, tolerance 1e-12 ; precision_schur_repro.cpp.
Limites de cette correction : Preuve pour API matrices creuses ; pas de qualification GPU ni de tous les operateurs matriciels externes.
PR affichee correspondante : #1092. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
G04 — BiCGStab : seuils de rupture absolus et diagnostic nul sur échec
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/core/linalg/bicgstab_solver.h:108-124; src/cfdx/core/linalg/bicgstab_solver.h:165-168.
Portee : API BiCGStab ; la panne est reproduite sans préconditionneur. Les choix de repli du solveur appelant peuvent en masquer l’effet final.
Des produits scalaires sont comparés à 1e-30 sans mise à l'échelle. Les quantités testées n'ont pas toutes la même dimension, ce qui empêche un seuil absolu unique d'avoir une signification générale.
Contre-exemple : (sx=s), donc (x_*=1), conditionnement 1. Pour (s=1), convergence en une itération. Pour (s=10^{-16}), le premier produit (\widehat r^Tr) vaut (10^{-32}) et est déclaré en rupture avant la première mise à jour. Retour DIVERGED, (x=0), résidu relatif réel 1.
Le SolverResult indique en outre un résidu nul, resté à sa valeur par défaut. Il ne faut pas confondre ce zéro avec une solution calculée.
Test de non-régression : multiplication commune de (A,b) par plusieurs puissances de dix, matrice bien conditionnée, statut et vrai résidu recalculés. Tester séparément les vraies ruptures algébriques.
Modification effectuee : Les seuils rho et r.v de BiCGStab sont relatifs aux normes des vecteurs ; le statut d'echec rapporte le residu effectivement mesure. tt et omega sont testes pour nullite/finitude.
Tests et preuves : solver_repro.cpp et test_audit_numerics: sx=s, s=1 et 1e-16.
Limites de cette correction : Le controle d'echelle n'est pas une garantie de convergence de BiCGStab sur toute matrice non normale.
PR affichee correspondante : #1092. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
G05 — GMRES : second membre nul refusé avec un itéré initial non nul
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P2.
Localisation dans le snapshot source : src/cfdx/core/linalg/gmres_solver.h:88-94.
Portee : Noyau GMRES partagé avec FGMRES. Exécution reproduite sur GMRES ; portée FGMRES déduite du noyau commun.
Le code retourne immédiatement DIVERGED quand (|b|=0) et (|b-Ax_0|\ne0).
Contre-exemple : (A=[1]), (b=[0]), (x_0=[1]). La solution est exactement 0 et une itération suffit en arithmétique exacte. Le solveur retourne DIVERGED, zéro itération, (x=1).
Ce cas est pertinent pour un solveur réutilisant un itéré antérieur lorsque le second membre s'annule. Cette pertinence ne démontre pas qu'un tutoriel particulier rencontre déjà ce chemin.
En cas de système singulier, une solution non nulle peut appartenir au noyau : ne pas remplacer aveuglément toute solution par zéro en présence d'une jauge ou d'un contrat sur le noyau. Traiter la convergence du système effectivement demandé.
Modification effectuee : Un second membre nul ne signifie pas que le residu initial est nul. GMRES utilise le residu initial comme reference quand ||b||=0, puis resout A*x=0.
Tests et preuves : test_audit_numerics: GMRES/zero-rhs ; solver_repro.cpp.
Limites de cette correction : Cas x0 non nul couvert ; noyau commun GMRES/FGMRES, cas direct teste en GMRES.
PR affichee correspondante : #1092. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
G06 — GMRES : produits scalaires Arnoldi incompatibles avec la politique de réduction
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/core/linalg/gmres_solver.h:145-172; src/cfdx/core/linalg/krylov_reductions.h:60-89.
Portee : API de plage locale/ghosts, reproduite en CPU monoprocessus. Incohérence MPI déduite du code ; aucune campagne multirang exécutée.
Les normes passent par krylov_norm2(..., controls.reduction), qui respecte local_begin/local_end et la réduction MPI. Les coefficients de Hessenberg sont au contraire calculés par des boucles locales sur toutes les entrées du vecteur. La réorthogonalisation reproduit la même incohérence.
Contre-exemple CPU : matrice identité (2\times2), (b=(1,1)), deux entrées dont la seconde représente une copie fantôme. Avec la plage complète, une itération donne un résidu relatif de (2,22\times10^{-16}). Avec local_end=1, le solveur atteint sa limite de dix itérations avec un résidu relatif owned de (1,71143995\times10^{-4}), au lieu de (10^{-10}).
Il s'agit d'un test du contrat de réduction de l'API, pas d'une preuve d'exécution d'un calcul CFD distribué. Une qualification MPI doit en plus vérifier plusieurs partitions, le rafraîchissement des halos et les réductions collectives.
Modification effectuee : Les projections d'Arnoldi et la reorthogonalisation utilisent le meme intervalle owned et la meme reduction que les normes.
Tests et preuves : solver_repro.cpp: matrice identite avec intervalle local partiel ; test_gmres_solver et test_fgmres_solver.
Limites de cette correction : TEST CPU DU CONTRAT ONLY : pas de lancement MPI multirang, pas de preuve de synchronisation des halos.
PR affichee correspondante : #1092. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
G07 — QUICK : repli linéaire sur un stencil vide
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/core/numerics/gradient.h:397-519; src/cfdx/core/numerics/gradient.h:503-506; src/cfdx/core/numerics/interpolation.h:178-185; src/cfdx/core/numerics/interpolation.h:214-230.
Portee : Fit quadratique utilisé par QUICK/QUICK_BOUNDED dans la bibliothèque. Test exécuté sur QUICK. La CLI de production de cet instantané refuse QUICK.
Lorsque le nombre d'échantillons est inférieur à neuf, compute_quadratic_fit_least_squares appelle least_squares_gradient avec deux vecteurs vides. Le gradient de repli est donc nul, même quand le stencil permet une reconstruction affine correcte.
Contre-exemple : cube (2^3), (\phi=x), flux positif. Le gradient reconstruit est 0 au lieu de 1 ; l'erreur maximale aux faces internes est 0,25. Une face à (x=0,5) reçoit la valeur amont 0,25 au lieu de 0,5.
Contrôle positif : cube (3^3), même champ et même reconstruction, erreur de face (2,78\times10^{-17}).
La nécessité d'un repli lorsque le stencil quadratique est insuffisant n'est pas un bug. Le bug est de perdre les données disponibles et donc l'exactitude du premier ordre au lieu d'effectuer ce repli.
Modification effectuee : Le repli lineaire du fit QUICK reutilise les positions et valeurs deja collectees, au lieu de deux listes vides.
Tests et preuves : test_audit_numerics et quadratic_repro.cpp sur cube 2^3, champ phi=x.
Limites de cette correction : Les stencils qui ne permettent meme pas une reconstruction affine 3D ne deviennent pas pleinement determines.
PR affichee correspondante : #1092. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
G08 — QUICK : base quadratique non normalisée et rang dépendant des unités
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/core/numerics/gradient.h:440-502; src/cfdx/core/numerics/gradient.h:531-802.
Portee : Implémentation du fit QUICK distincte du gradient quadratique normalisé. Ne pas attribuer le défaut au second noyau, qui réussit le contrôle.
La base du fit QUICK juxtapose (dx,dy,dz) et (dx^2,dy^2,dz^2,dxdy,\ldots) sans normalisation. Le seuil de pivot est déduit du plus grand coefficient diagonal, puis des pivots quadratiques peuvent être ignorés. have_quad=true est pourtant attribué sans exiger un rang quadratique complet.
Contre-exemple : même cube (4^3), coordonnées multipliées par (L), champ (\phi=(x/L)^2). On compare les coefficients adimensionnés : (L g_x) à (2x/L), et (L^2 a_{xx}) à 1 (coefficient polynomial, pas la Hessienne qui comporte un facteur 2).
L
Erreur max de L g, fit QUICK
Erreur max de L2 a_xx
Erreur max de L g, gradient quadratique normalisé
1
2,4425e-15
5,7732e-15
6,2172e-15
1e-4
5,9952e-15
1,6320e-14
6,8834e-15
1e-6
0,3375
1
6,2172e-15
1e-8
0,3375
1
5,3291e-15
Le coefficient de courbure attendu non nul est perdu aux petites échelles. Le noyau quadratique normalisé reste exact sur la même géométrie : ce contrôle distingue le défaut de base et de rang d'un problème général de calcul des centres.
Les équations normales méritent aussi une revue de conditionnement. Une factorisation QR avec diagnostic de rang, ou une SVD lorsque le rang est déficient, constitue une piste de correction ; ce rapport ne prétend pas qu'une telle correction a été implémentée.
Modification effectuee : Les inconnues du fit sont adimensionnees par une longueur locale h : termes lineaires/h et quadratiques/h^2. Tous les pivots requis doivent etre disponibles, sinon repli lineaire.
Tests et preuves : quadratic_repro.cpp : echelles 1,1e-4,1e-6,1e-8 ; test_audit_numerics: QUICK/scale/curvature.
Limites de cette correction : Il reste une resolution sur equations normales et un test de pivots, pas une SVD/QR revelatrice de rang sur tout stencil pathologique.
PR affichee correspondante : #1092. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
G09 — Oracle ExactSchur : NaN traité comme un second membre nul
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P2.
Localisation dans le snapshot source : src/cfdx/core/linalg/exact_schur.h:83-99.
Portee : Outil de référence Schur utilisé par des tests. La recherche de raccordement ne démontre pas son usage dans le chemin principal de production.
Le code utilise if (!(r0 > 0.0)) return true; pour détecter un second membre nul. Cette expression est également vraie pour NaN.
Contre-exemple : blocs scalaires (A_{uu}=1, G=1, D=-1, C=0). Le complément (S=C-DA_{uu}^{-1}G) vaut 1. Avec un second membre égal à 1, le solveur retourne correctement 1. Avec NaN, apply() retourne true et une pression nulle.
Ce défaut menace l'intégrité d'un oracle numérique. Il ne suffit pas à conclure que tous les tests Schur actuels acceptent les NaN : ils peuvent avoir des contrôles supplémentaires. L'API de cet oracle accepte néanmoins l'entrée invalide dans le cas exécuté.
Modification effectuee : ExactSchur rejette RHS non fini et norme initiale non finie ; seule une norme exactement nulle produit la solution nulle immediate.
Tests et preuves : test_audit_numerics: Schur/reject-NaN ; precision_schur_repro.cpp.
Limites de cette correction : Ne prouve pas la correction de tous les oracles Schur ni du couplage F03.
PR affichee correspondante : #1092. Titre observe, diff et CI exact-HEAD non verifies ; ne pas confondre avec ce correctif local.
A01 — Bilan energetique global susceptible de compenser des erreurs locales
Etat : NON CORRIGE / A QUALIFIER ; priorite P2.
Localisation dans le snapshot source : src/cfdx/physics/energy_solver.h:138-160; src/cfdx/physics/radiation_solver.h:576.
Constat sur le code initial : La somme signee des residus peut etre nulle pour deux cellules ayant des residus opposes. Le premier audit etait trop general : le solveur energie principal possede aussi des controles locaux et independants ; le helper global ne suffit pas a le condamner.
Solution a implementer, non appliquee dans le patch : Conserver le bilan global mais exposer egalement max_i |r_i|/scale_i et L1 ; inspecter chaque consommateur, notamment les couplages radiatifs, avant de conclure a une fausse acceptation.
Critere de cloture : Deux cellules +q/-q ; le global doit etre nul et le local strictement positif. Ne pas supprimer le bilan global.
A02 — Compilation CPU sans MPI interrompue par une inclusion MPI non conditionnelle
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/core/linalg/communication_avoiding.h; tests/unit/test_performance_runtime.cpp:7.
Constat sur le code initial : Reproduit dans la compilation complete, puis corrige sans exclure le test.
Modification effectuee : Les inclusions et helpers MPI dans communication_avoiding.h sont conditionnes par CFDX_HAS_MPI. Le test_performance_runtime reste compile et execute en mode serial.
Tests et preuves : Compilation CMake complete CPU, test_performance_runtime.
Limites de cette correction : Aucune execution MPI/GPU dans cette campagne.
A03 — Dependance meshio non declaree pour les chemins de conversion qui l'exigent
Etat : NON CORRIGE / A QUALIFIER ; priorite P2.
Localisation dans le snapshot source : pyproject.toml; python/cfdx/io/pipeline.py:144.
Constat sur le code initial : Les campagnes precedentes signalaient des echecs de conversion faute de meshio. Ce point de packaging est conserve, sans le compter comme defaut d'equation.
Solution a implementer, non appliquee dans le patch : Declarer un extra I/O documente ou une dependance obligatoire selon le contrat de la CLI ; message explicite quand l'extra n'est pas installe. Tester un wheel en environnement vierge.
Critere de cloture : Installation du wheel avec/sans extra ; conversions OpenFOAM/HDF5/VTU et messages d'erreur attendus.
A04 — Une diffusivite NaN peut etre transformee en zero dans le flux
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/finite_volume_transport.h:318-325.
Constat sur le code initial : Le comportement venait des comparaisons avec zero avant moyenne harmonique. Une entree materiau invalide pouvait etre masquee.
Modification effectuee : Les diffusivites cellulaires non finies ou negatives sont rejetees avant le calcul harmonique au lieu de devenir silencieusement nulles.
Tests et preuves : test_audit_numerics: diffusion/reject-NaN-material.
Limites de cette correction : Validation locale des valeurs ; la coherence thermodynamique globale reste du ressort des modeles consommateurs.
A05 — Qualification inachevee de familles numeriques et physiques
Etat : NON CORRIGE / A QUALIFIER ; priorite P2.
Localisation dans le snapshot source : docs/validation/NUMERICAL_METHOD_MATURITY_AUDIT.json; QUALIFICATION_REGISTRY.md.
Constat sur le code initial : Il s'agit d'une dette de qualification, pas d'une 52e equation fausse. READY ou presence d'un noyau ne signifie pas PASS sur le solveur integre.
Solution a implementer, non appliquee dans le patch : Relier chaque promotion a un SHA exact, un oracle independant et un chemin de production ; garder les limitations MPI/GPU/multiphysiques visibles.
Critere de cloture : Campagnes de verification et validation et artefacts reproductibles sur le commit fusionne.
H01 — Normes instables et fausses convergences par overflow/underflow
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/core/linalg/vector.h:63-74; src/cfdx/core/linalg/krylov_reductions.h; src/cfdx/core/linalg/jacobi_solver.h; src/cfdx/core/linalg/gauss_seidel_solver.h.
Equation / identite attendue :
[ |x|_2=s\sqrt{\sum_i(x_i/s)^2},\quad s=\max_i|x_i|. ]
Constat sur le code initial : La norme de (1e200) devenait infinie et celle de (1e-200) devenait zero ; norm_inf(NaN) pouvait devenir zero. Jacobi/Gauss-Seidel pouvaient accepter x=0 pour x=1e200 parce que inf<=tol*inf.
Modification effectuee : Les normes vectorielles utilisent hypot ou une mise a l'echelle avant somme des carres ; NaN est propage/rejete. Jacobi/Gauss-Seidel refusent les normes initiales non finies.
Tests et preuves : additional_repro.cpp ; test_audit_numerics: norm/, stationary/huge/.
Limites de cette correction : La robustesse de toutes les reductions/dot products de CFDX n'est pas revendiquee ; MPI non execute.
H02 — Newton retourne un echec apres son dernier pas pourtant convergent
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P2.
Localisation dans le snapshot source : src/cfdx/core/numerics/newton.h:77-146.
Equation / identite attendue :
[ J(u^k),\delta u=-R(u^k),\quad u^{k+1}=u^k+\alpha\delta u. ]
Constat sur le code initial : R(u)=u-1, u0=0, un pas autorise : u=1 mais statut MAX_ITER_REACHED. Le critere etait teste au debut de l'iteration suivante, qui n'existait pas.
Modification effectuee : Newton teste le residu de l'etat accepte apres chaque pas, y compris le dernier pas autorise, avant de retourner MAX_ITER_REACHED.
Tests et preuves : test_audit_numerics: newton/last-iteration/status et newton/solution ; additional_repro.cpp.
Limites de cette correction : Test lineaire a budget un ; pas de qualification de tous les Jacobiennes non lineaires.
H03 — Diagonale matrix-free artificielle et coefficient NaN accepte
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P2.
Localisation dans le snapshot source : src/cfdx/core/linalg/matrix_free_fv_operator.h:12-18.
Equation / identite attendue :
[ D_{ii}=\mathbf e_i^T A\mathbf e_i. ]
Constat sur le code initial : Pour gamma=0, apply(e_i)=0 mais diagonal(i)=1. Le constructeur gamma<0 ne rejetait pas NaN.
Modification effectuee : La diagonale rendue par l'operateur est sa vraie diagonale, y compris zero. Un coefficient de diffusion NaN est rejete. Une protection de preconditionneur ne peut pas changer silencieusement le contrat de l'operateur.
Tests et preuves : test_audit_numerics: matrixfree/zero-diagonal, matrixfree/reject-NaN ; additional_repro.cpp.
Limites de cette correction : Un operateur nul reste singulier ; cette correction ne le rend pas inversible.
H04 — Tolerance relative traitee comme absolue dans le solveur scalaire 1x1
Etat : CORRIGE ET TESTE (perimetre indique) ; priorite P1.
Localisation dans le snapshot source : src/cfdx/physics/finite_volume_transport.h:711-743.
Equation / identite attendue :
[ \eta=|a,x-b|/|b|\le\mathrm{tol}\quad(b\ne0). ]
Constat sur le code initial : La campagne CHT asymetrique a revele des refus avec residu absolu 1.46e-11 mais erreur relative proche de 1e-16. Le champ residual_relative contenait le residu absolu. Le defaut s'ajoutait a F04 et empechait quatre cas de converger a tolerance 1e-12.
Modification effectuee : Le solveur direct 1x1 rapporte |Ax-b|/|b| et applique la tolerance relative annoncee. Il rejette les valeurs non finies et ne remplace pas l'etat accepte par un predicteur refuse.
Tests et preuves : test_audit_cht_time : quatre cas initialement rejetes malgre une erreur retrograde ~6e-17 passent sans changer la tolerance 1e-12.
Limites de cette correction : Le contrat de sous-relaxation multivariable F32 demeure ouvert ; le chemin exact 1x1 conserve sa semantique anterieure.
H05 — Diagnostic d'énergie : distance euclidienne incompatible avec l'assemblage
État : CORRIGÉ ET TESTÉ R2 ; priorité P2.
Localisation : src/cfdx/physics/energy_solver.h, fonction energy_balance_relative, instruction const double d=dvec.mag() avant correction.
Défaut : la matrice et la reconstruction principale utilisent la distance normale à la frontière, tandis que ce diagnostic global utilisait la distance euclidienne pour le flux de Dirichlet. Un champ satisfaisant exactement l'équation discrète pouvait donc recevoir un bilan diagnostic fortement non nul.
Oracle indépendant : une cellule de volume et aire unitaires, (\rho c_p=1), (\Delta t=0,1), (k=2), (d_n=0,5), (T^n=300) et (T_b=400) donne [ 10(T-300)+4(T-400)=0,\qquad T=\frac{4600}{14}. ] Les autres faces sont adiabatiques. Le cisaillement géométrique conserve le volume, l'aire de l'interface et sa distance normale.
Correction : const double d=boundary_normal_distance(geometry,f,o);. Aucun seuil de convergence ni coefficient de la matrice n'a été modifié pour fermer ce constat.
Preuves : trois cisaillements (0, 1, 3), avec condition de Dirichlet déclarée puis valeur de face explicite. Avant : matrice satisfaite, mais diagnostic 0,2928932188 pour le cisaillement 1 et 0,6837722340 pour le cisaillement 3 ; 4 échecs/12. Après : 12/12, diagnostic maximal 7,9581×10^-16. ASan/UBSan ne signale pas d'anomalie sur ce chemin instrumenté.
Portée : correction d'un helper public de bilan, pas correction de toutes les formulations d'énergie. solve_energy n'utilisait pas ce helper comme unique critère de convergence. Les limites de bilan global A01 et de diffusion non orthogonale F06 restent ouvertes.
7. Intégrer R2 depuis le Markdown
Utiliser le script complet de l'annexe C. Il vérifie le SHA-256 du ZIP source et du patch embarqué, refuse d'écraser la destination, extrait une copie, applique git apply --check puis le patch, et vérifie les empreintes des 43 fichiers modifiés/ajoutés. Il ne fait ni commit ni push.
python rejouer_depuis_md.py \
  --md CFDX_audit_consolide_correctifs_testes_2026-10-24.md \
  --archive cfdx_dev-master.zip \
  --output cfdx-revision-r2 \
  --test --sanitizers --jobs 2
Sans --test ni --sanitizers, le script extrait et vérifie seulement la copie corrigée. Les tests ont été exécutés ici sous Linux, GCC 14.2.0 et CMake 3.31.6 ; les sanitizers du script nécessitent GCC/Clang dans un environnement compatible (par exemple Linux/WSL). Ce n'est pas une qualification Windows/MSVC.
Ne pas appliquer ce patch sur une branche de PR différente. Il vise uniquement le ZIP original identifié en section 2. Pour une PR actuelle, porter d'abord les reproducteurs sur son commit exact, comparer ses modifications et relancer sa campagne.
Fichiers modifiés ou ajoutés dans le patch cumulatif
Fichier
Nature
CMakeLists.txt
Modification
src/cfdx/core/linalg/bicgstab_solver.h
Modification
src/cfdx/core/linalg/cg_solver.h
Modification
src/cfdx/core/linalg/communication_avoiding.h
Modification
src/cfdx/core/linalg/exact_schur.h
Modification
src/cfdx/core/linalg/gauss_seidel_solver.h
Modification
src/cfdx/core/linalg/gmres_solver.h
Modification
src/cfdx/core/linalg/jacobi_solver.h
Modification
src/cfdx/core/linalg/krylov_reductions.h
Modification
src/cfdx/core/linalg/matrix_free_fv_operator.h
Modification
src/cfdx/core/linalg/vector.h
Modification
src/cfdx/core/numerics/gradient.h
Modification
src/cfdx/core/numerics/interpolation.h
Modification
src/cfdx/core/numerics/newton.h
Modification
src/cfdx/core/numerics/temporal.h
Modification
src/cfdx/physics/cht_solver.h
Modification
src/cfdx/physics/dual_time_driver.h
Modification
src/cfdx/physics/dual_time_field_integrator.h
Modification
src/cfdx/physics/dual_time_stepping.h
Modification
src/cfdx/physics/energy_solver.h
Modification
src/cfdx/physics/finite_volume_transport.h
Modification
src/cfdx/physics/low_storage_time_integration.h
Modification
src/cfdx/physics/m1_m4_models.h
Modification
src/cfdx/physics/sst_solver.h
Modification
src/cfdx/physics/thermal_zones.h
Modification
src/cfdx/physics/thermophysical_models.h
Modification
src/cfdx/physics/transport_models.h
Modification
src/cfdx/physics/turbulence_solver.h
Modification
src/cfdx/thermodynamics/equation_of_state.h
Modification
tests/audit/affine_mesh_fixture.h
Ajout
tests/audit/chain_fixture.h
Ajout
tests/audit/test_audit_cht_contract.cpp
Ajout
tests/audit/test_audit_cht_time.cpp
Ajout
tests/audit/test_audit_energy_boundary_metric.cpp
Ajout
tests/audit/test_audit_flux_contract.cpp
Ajout
tests/audit/test_audit_numerics.cpp
Ajout
tests/unit/test_interpolation.cpp
Modification
tests/unit/test_transport_conservation.cpp
Modification
tests/unit/test_turbulence_qualification_equations.cpp
Modification
tests/validation/test_cht_validation.cpp
Modification
tests/validation/test_level_c_coupled_verification.cpp
Modification
tests/validation/test_polyhedral_gradient_campaign.cpp
Modification
tests/validation/test_thermal_challenge_vv.cpp
Modification
8. Priorites restantes
Premier lot restant : F03, F06, F08 et F11. Ce sont des changements d'equation ou de couplage, pas de simples gardes numeriques. Il faut unifier la definition du flux de pression, integrer la correction non orthogonale dans le transport consomme par les solveurs, traiter la skewness et le rang de reconstruction, et assembler le tenseur visqueux complet lorsque mu varie.
Deuxieme lot : F12/F13/F18/F22–F25/F27/F30/F32 et A01/A03/A05. Les modeles ou aides incomplets doivent etre corriges, explicitement renommes ou refuses. Une API qui ne dispose pas de la variable physique necessaire ne peut pas reconstruire cette physique par une constante implicite. Les PR correspondantes doivent etre evaluees a leur SHA de tete, avec les contre-exemples qui motivent leur correction.
Ne pas annoncer une qualification globale de CFDX tant que ces lots restent ouverts. Le resultat de cette livraison est une base corrigee et testee sur un perimetre explicite, pas une certification des 121 formulations recensees anterieurement.
9. Sources scientifiques et documentaires
•	Source primaire du code et reference versionnee : archive fournie par l'utilisateur, SHA et empreinte en section 2 ; aucune source secondaire substituee au code.
•	SST : TMR, definition SST-2003 et erratum sur la production omega : https://tmbwg.github.io/turbmodels/sst.html .
•	Wilcox : TMR, limiteur de viscosite et coefficients diffusifs distincts : https://tmbwg.github.io/turbmodels/wilcox.html .
•	Normes de convergence Krylov : PETSc KSPConvergedDefault : https://petsc.org/release/manualpages/KSP/KSPConvergedDefault/ .
•	Moindres carres et rang : LAPACK Users' Guide : https://netlib.org/lapack/lug/node27.html .
•	Statuts historiques internes : QUALIFICATION_REGISTRY.md, lignes 21–44, et registres JSON designes par ce document.
•	Liste publique consulte pendant V1, non revalidee par R2 : https://github.com/stevendaix/cfdx_dev/pulls ; pas de confirmation des derniers diffs.
Les formules de conduction en serie, BDF2 variable, primitive d'enthalpie et normes mises a l'echelle sont egalement verifiees directement par les contre-exemples analytiques livres. Les journaux signalant les limites et echecs sont conserves ; aucune erreur ouverte n'est remplacee par une simple promesse de correction.

10. Travail ajouté et exécuté dans R2
10.1 Comparaison des versions
La base du comparatif R2 est V1 : le ZIP source plus le patch de la livraison précédente, et non une branche GitHub courante. L'ancien fix F04 de l'état temporel est donc déjà présent dans le test « avant » ; les nouveaux écarts ne lui sont pas imputables. Les deux nouveaux tests sont compilés avec les mêmes sources de test contre V1 puis R2.
Suite
Avant R2, sur V1
Après R2
Référence
test_audit_cht_contract
30 réussites, 71 échecs / 101
101/101
Réseau de résistances + deux capacités + rejets de contrats
test_audit_energy_boundary_metric
8 réussites, 4 échecs / 12
12/12
Équation scalaire analytique assemblée indépendamment
test_audit_numerics réexécuté
Voir historique V1
87/87
Oracles des corrections antérieures
test_audit_flux_contract réexécuté
Voir historique V1
56/56
Identité flux/matrice
test_audit_cht_time réexécuté
Voir historique V1
24/24
Un pas physique, six configurations
Total interne des cinq suites R2
Non additionné aux comptes CTest
280/280
Périmètre défini ci-dessus
10.2 Un cas F28 concret, sans changer de tolérance
Deux cellules ont les longueurs normales (L_1=0,4), (L_2=1,6), des cisaillements géométriques 2 et -0,5, les conductivités (k_1=1), (k_2=4), les capacités volumiques 2 et 3 et les sources volumiques 10 et -5. Le pas vaut 0,1 et la relaxation d'interface 1. Les températures initiales sont 300 et 400 K.
Avec (G=A/(L_1/(2k_1)+L_2/(2k_2))), (M_i=\rho_i c_{p,i}V_i/\Delta t) et (b_i=M_iT_i^n+Q_iV_i), l'oracle est : [ \begin{pmatrix}M_1+G&-G\-G&M_2+G\end{pmatrix} \begin{pmatrix}T_1\T_2\end{pmatrix} =\begin{pmatrix}b_1\b_2\end{pmatrix}. ]
Valeur (K)
Avant R2
Après R2
Oracle discret
(T_1)
329,016746411481
323,248091603047
323,248091603053
(T_2)
397,456937799043
396,041984732823
396,041984732824
Ce cas illustre un défaut de raccordement entre conductances. Il n'est pas une validation continue générale du schéma sur maillage non orthogonal ; F06 demeure ouvert.
10.3 Trois fixtures de validation corrigées, pas trois seuils assouplis
Après durcissement du contrat, les trois tests ci-dessous ont d'abord échoué. Le journal de ce premier passage est conservé dans l'annexe F. Les causes ont été corrigées dans les données de test, sans modification de leurs seuils.
Test
Incohérence démontrée
Changement dans R2
test_cht_validation
Deux cubes superposés, toutes les faces utilisées comme interfaces ; normales sortantes de même sens.
Deux cubes adjacents, une face commune, autres faces adiabatiques ; aire de la vraie interface dans le bilan.
test_level_c_coupled_verification
Les deux régions occupaient le même volume et partageaient une face de même orientation.
Translation de la région froide et appariement de sa face inférieure avec la face supérieure de la région chaude ; noms des frontières adaptés.
test_thermal_challenge_vv
energy2.conductivity=2 par copie du matériau chaud, alors que le couplage et l'oracle demandaient 1.
energy2.conductivity=1, aligné sur les températures analytiques déjà attendues.
Ces trois cas passent après correction. Le constat transversal F34 reste partiel : il n'est pas déclaré résolu pour tous les tests du dépôt.
10.4 Régression et sanitizers : périmètre exact
Les dix cibles compilées et exécutées sont les cinq suites test_audit_* du manifeste, test_level_c_coupled_verification, test_cht_validation, test_thermal_vv, test_thermal_zones_vv et test_thermal_challenge_vv. 10/10 passent. Aucun test ignoré dans cette sélection. Ce n'est pas une compilation de tous les exécutables ou une réexécution des 220 tests V1.
Les deux nouvelles suites passent avec -fsanitize=address,undefined, UBSAN_OPTIONS=halt_on_error=1:print_stacktrace=1 et ASAN_OPTIONS=detect_leaks=0. La détection des fuites est donc désactivée. La bibliothèque libcfdx_core.a est le build CPU Release ; il ne s'agit pas d'une instrumentation complète de cette bibliothèque.
Une première compilation R2 a signalé l'absence d'un opérateur Vec3 / double dans l'API du dépôt ; les normales sont maintenant normalisées composante par composante. Une construction regroupée a ensuite été interrompue par la limite de durée de l'appel d'outil, puis relancée jusqu'à son terme. Seuls les résultats finaux ci-dessus sont comptés comme validations.
10.5 Source primaire et limites de publication
La pondération par conductances d'interface est également cohérente avec la structure kappa * deltaCoeffs et les poids de la condition couplée documentée dans OpenFOAM Foundation v10 : https://cpp.openfoam.org/v10/turbulentTemperatureCoupledBaffleMixedFvPatchScalarField_8C_source.html . Cette référence a été consultée pour la structure de l'équation ; aucun code de cette source n'a été copié dans le patch. L'oracle numérique est calculé indépendamment dans le test fourni.
Les corrections ont été effectuées sur une copie locale. Les PR actuelles et leur CI n'ont pas été relues dans R2. Les fichiers et tests sont fournis pour examen et intégration, sans allégation de fusion ni de qualification de master.
Annexe A. Historique V1 - analyse des regressions de la campagne
La premiere campagne large (avant la derniere restauration Neumann et le correctif du caller WLS de validation) a execute 220 tests : 212 ont reussi, 8 ont echoue. Les journaux ne sont pas ecrases : logs/before_neumann_final_ctest.log et before_neumann_ctest.xml.
Les deux echecs numeriques apparents de cette passe etaient test_polyhedral_gradient_campaign et son agregateur test_n10_qualification. Le caller demandait implicitement EXCLUDE_BOUNDARY, alors que son fonctionnement anterieur dependait de fantomes ajoutes malgre cette demande. Apres respect de l'enumeration par G01, les cellules tetraedriques de bord revelaient un rang deficient. La campagne qui mesure l'erreur interieure demande maintenant explicitement ZERO_GRADIENT_GHOST et ajoute un test de rejet sous exclusion. Les controles de rang et les seuils d'ordre ne sont pas reduits.
Une deuxieme faiblesse du test a ete corrigee : smooth_errors etait recree a chaque taille de maillage, donc le controle de decroissance ne comparait jamais deux maillages. L'historique est maintenant conserve ; les mesures de decroissance sont effectives. Ce point appartient a F34, pas a un nouveau modele physique.
Cinq echecs de cette passe atteignaient explicitement un ModuleNotFoundError: No module named 'meshio' : test_cross_import_equivalence, test_openfoam_results_cli_production_e2e, test_production_solver_e2e, test_meshio_bridge, test_mesh_geometry_matrix. C'est un pre-requis manquant dans l'environnement et une dette de dependance A03 ; ces tests ne sont pas comptes comme reussis ou qualifies.
Enfin, test_n11_conservation_boundedness refuse de produire une preuve exact-HEAD a partir d'une extraction ZIP sans repertoire .git. Le SHA du ZIP et les empreintes de la copie modifiee sont connus, mais ils ne sont pas substitues frauduleusement a un commit construit/fusionne. Ce refus est garde.
La campagne finale est relancee apres reconstruction de tous les executables sur les sources finales. Ses chiffres officiels sont ceux des sections 1 et 4 et de final_ctest.xml, pas ceux de cette premiere passe.
Annexe B. Historique V1 - resultat final de cette campagne anterieure
La compilation complete de la configuration CPU Release a termine avec le code de sortie 0. Les 220 tests de la selection finale ont ensuite ete executes sur les sources finales : 214 reussites, 6 echecs, 0 ignore. Les trois nouveaux executables sont inclus dans ces 220 tests ; leurs 167 controles internes ne sont pas additionnes artificiellement au nombre CTest.
Les deux tests de validation polyedrique qui echouaient dans la premiere passe (test_polyhedral_gradient_campaign et test_n10_qualification) reussissent maintenant avec la politique frontiere explicite et les controles de rang et de convergence conserves/renforces.
Les cinq echecs ci-dessous ont chacun une trace explicite de meshio manquant :
•	test_cross_import_equivalence
•	test_openfoam_results_cli_production_e2e
•	test_production_solver_e2e
•	test_meshio_bridge
•	test_mesh_geometry_matrix
Le sixieme echec, test_n11_conservation_boundedness, retourne INCOMPLETE et git_sha=null car la source provient d'un ZIP sans .git. Aucun faux commit, HEAD ou certificat de provenance n'a ete cree pour contourner ce refus. Ces six tests ne sont pas valides, meme si leur cause actuelle n'est pas une erreur numerique demontree. La dependance et la provenance doivent etre retablies avant de qualifier les chemins concernes.
Les 56 contrats de flux passent avec ASan/UBSan. L'erreur maximale flux/matrice mesuree est 1.4210855e-14. Les six cas CHT donnent une erreur de temperature maximale de 3.4674486e-11 K et un defaut maximal d'energie de 1.5194246e-10 dans les unites du probleme.
Les maximas des tests sont a rejet explicite de NaN/Inf. La suite numerique recompilee avec ces memes gardes echoue toujours sur 58/87 controles du code initial et passe 87/87 sur la copie corrigee.

Annexe C. Script autonome de reproduction depuis ce Markdown
Copier ce bloc dans rejouer_depuis_md.py, puis exécuter la commande de la section 7. Le script lit les annexes D et E directement dans ce document. Il refuse les données ambiguës ou non conformes et ne tente aucune action GitHub.
#!/usr/bin/env python3
"""Extract the embedded CFDX patch, apply it to an exact ZIP, optionally test it.
Python 3.9+, git, CMake and a C++17 compiler. Never modifies the input archive.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import stat
import subprocess
import zipfile


def embedded(text: str, name: str) -> str:
    lines = text.splitlines()
    begin = '<!-- CFDX_' + name + '_BEGIN -->'
    end = '<!-- CFDX_' + name + '_END -->'
    if lines.count(begin) != 1 or lines.count(end) != 1:
        raise ValueError('Missing or ambiguous embedded block: ' + name)
    first, last = lines.index(begin), lines.index(end)
    block = lines[first + 1:last]
    if len(block) < 3:
        raise ValueError('Empty embedded block: ' + name)
    fence_len = len(block[0]) - len(block[0].lstrip('`'))
    if fence_len < 3 or block[-1] != '`' * fence_len:
        raise ValueError('Invalid code fence: ' + name)
    return '\n'.join(block[1:-1]) + '\n'


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(command: list[str], cwd: Path | None = None, env: dict | None = None) -> None:
    print('+', ' '.join(map(str, command)), flush=True)
    subprocess.run(command, cwd=cwd, env=env, check=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--md', required=True, type=Path)
    p.add_argument('--archive', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path, help='New, non-existing directory')
    p.add_argument('--test', action='store_true', help='Build and run the 10 R2 CTest targets')
    p.add_argument('--sanitizers', action='store_true', help='Also run the two new suites under ASan/UBSan')
    p.add_argument('--jobs', type=int, default=2)
    args = p.parse_args()
    if args.jobs < 1:
        p.error('--jobs must be positive')
    report = args.md.read_text(encoding='utf-8')
    manifest = json.loads(embedded(report, 'MANIFEST'))
    patch = embedded(report, 'PATCH').encode('utf-8')
    if hashlib.sha256(patch).hexdigest() != manifest['patch_sha256']:
        raise ValueError('Embedded patch checksum mismatch')
    if digest(args.archive) != manifest['archive_sha256']:
        raise ValueError('Wrong source ZIP. Do not force onto a current PR/master.')
    out = args.output.resolve()
    if out.exists():
        p.error('--output already exists; refusing to overwrite')
    needed = ['git'] + (['cmake', 'ctest'] if args.test or args.sanitizers else [])
    for program in needed:
        if not shutil.which(program):
            p.error('Required program not found: ' + program)
    with zipfile.ZipFile(args.archive) as z:
        for member in z.infolist():
            name = PurePosixPath(member.filename)
            if (name.is_absolute() or '..' in name.parts or not name.parts
                    or name.parts[0] != 'cfdx_dev-master'
                    or '\\' in member.filename or ':' in member.filename
                    or stat.S_ISLNK(member.external_attr >> 16)):
                raise ValueError('Unsafe or unexpected ZIP member: ' + member.filename)
        out.mkdir(parents=True)
        z.extractall(out)
    src = out / 'cfdx_dev-master'
    for item in manifest['files']:
        f = src / item['path']
        if item['new']:
            if f.exists():
                raise ValueError('Unexpected existing new file: ' + str(f))
        elif digest(f) != item['sha256_before']:
            raise ValueError('Source checksum mismatch: ' + str(f))
    patch_path = out / 'correctifs.patch'
    patch_path.write_bytes(patch)
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n', encoding='utf-8')
    run(['git', 'apply', '--check', str(patch_path)], src)
    run(['git', 'apply', str(patch_path)], src)
    for item in manifest['files']:
        if digest(src / item['path']) != item['sha256_after']:
            raise ValueError('Patched file checksum mismatch: ' + item['path'])
    if args.test or args.sanitizers:
        build = out / 'build-audit'
        run(['cmake', '-S', str(src), '-B', str(build), '-DCMAKE_BUILD_TYPE=Release',
             '-DCFDX_ENABLE_MPI=OFF', '-DCFDX_ENABLE_GPU=OFF', '-DCFDX_BUILD_TESTS=ON',
             '-DCFDX_BUILD_RUNTIME_TESTS=OFF', '-DCFDX_BUILD_PYTHON=OFF',
             '-DCFDX_ENABLE_LONG_VALIDATION=OFF'])
        targets = manifest['ctest_targets']
        run(['cmake', '--build', str(build), '-j', str(args.jobs), '--target'] + targets)
        run(['ctest', '--test-dir', str(build), '--output-on-failure',
             '--output-junit', str(out / 'r2-results.xml'),
             '-R', '^(' + '|'.join(targets) + ')$'])
        if args.sanitizers:
            compiler = shutil.which('g++') or shutil.which('clang++')
            if not compiler:
                raise RuntimeError('GCC or Clang required for sanitizer replay')
            for target in ['test_audit_cht_contract', 'test_audit_energy_boundary_metric']:
                exe = build / (target + '-sanitized')
                run([compiler, '-std=c++17', '-O1', '-g1', '-DCFDX_HAS_MPI=0',
                     '-fsanitize=address,undefined', '-fno-omit-frame-pointer',
                     '-I', str(src / 'src'), '-I', str(src / 'tests/audit'),
                     str(src / ('tests/audit/' + target + '.cpp')),
                     str(build / 'libcfdx_core.a'), '-o', str(exe)])
                env = dict(os.environ, ASAN_OPTIONS='detect_leaks=0',
                           UBSAN_OPTIONS='halt_on_error=1:print_stacktrace=1')
                run([str(exe)], env=env)
    print('Done: exact source copy patched; no GitHub push or PR validation performed.')


if __name__ == '__main__':
    main()

Annexe D. Patch cumulatif integral - code et tests
Base exclusive : ZIP SHA-256 11f92e9f842d16b184432051208c40262a0961085d86ac933946b54613b955bb.
Patch SHA-256 : 8df8479192656cf47155ace39666ea245b4dcfb6d2bccb03fecb69697a12be59.
Ce bloc contient les 37 fichiers de la livraison V1 et les changements R2, soit 43 fichiers distincts au total. Il remplace le patch V1 ; ne pas appliquer les deux successivement.
--- a/CMakeLists.txt
+++ b/CMakeLists.txt
@@ -418,6 +418,19 @@
         target_include_directories(${test_name} PRIVATE ${CMAKE_SOURCE_DIR}/tests)
         add_test(NAME ${test_name} COMMAND ${test_name})
     endmacro()
+
+    # Audit counterexamples: active checks in Release as well as Debug.
+    add_cfdx_test(test_audit_numerics tests/audit/test_audit_numerics.cpp)
+    target_include_directories(test_audit_numerics PRIVATE ${CMAKE_SOURCE_DIR}/tests/validation)
+    add_cfdx_test(test_audit_flux_contract tests/audit/test_audit_flux_contract.cpp)
+    add_cfdx_test(test_audit_cht_time tests/audit/test_audit_cht_time.cpp)
+    set_tests_properties(test_audit_numerics test_audit_flux_contract test_audit_cht_time PROPERTIES
+        LABELS "audit;regression" TIMEOUT 60)
+
+    add_cfdx_test(test_audit_cht_contract tests/audit/test_audit_cht_contract.cpp)
+    add_cfdx_test(test_audit_energy_boundary_metric tests/audit/test_audit_energy_boundary_metric.cpp)
+    set_tests_properties(test_audit_cht_contract test_audit_energy_boundary_metric PROPERTIES
+        LABELS "audit;regression;thermal" TIMEOUT 60)
 
     # N8 physical/production qualification reuses the existing validation
     # executables and records a machine-readable campaign report. It is
--- a/src/cfdx/core/linalg/bicgstab_solver.h
+++ b/src/cfdx/core/linalg/bicgstab_solver.h
@@ -91,23 +91,26 @@
         result.residual_relative=bnorm>0?res/bnorm:0.0;return result;
     }
 
+    result.residual=res; result.residual_relative=bnorm>0 ? res/bnorm : res;
     double rho=1.0,alpha=1.0,omega=1.0;
     for(std::size_t i=0;i<n;++i){w.v(i)=0.0;w.p(i)=0.0;}
 
     auto true_residual = [&]() {
-        matvec(x, Ax);
+        mixed_precision_matvec(A, x, Ax, SolverPrecision::FP64);
         if (!finite_vector(Ax)) return std::numeric_limits<double>::infinity();
         for (std::size_t i = 0; i < n; ++i) {
             const double ri = b(i) - Ax(i);
             if (!std::isfinite(ri)) return std::numeric_limits<double>::infinity();
             w.r(i) = ri;
         }
-        return krylov_norm2(w.r, redp, reduction);
+        const double norm = krylov_norm2(w.r, SolverPrecision::FP64, reduction);
+        result.residual=norm; result.residual_relative=bnorm>0 ? norm/bnorm : norm;
+        return norm;
     };
 
     for(std::size_t iter=1;iter<=max_iter;++iter){
         const double rho_new=krylov_dot(w.r_hat, w.r, redp, reduction);
-        if(!std::isfinite(rho_new)||std::abs(rho_new)<1e-30){result.status=SolverStatus::DIVERGED;result.iterations=iter-1;return result;}
+        if(!std::isfinite(rho_new)||std::abs(rho_new) <= 64*std::numeric_limits<double>::epsilon()*krylov_norm2(w.r_hat,SolverPrecision::FP64,reduction)*krylov_norm2(w.r,SolverPrecision::FP64,reduction)){result.status=SolverStatus::DIVERGED;result.iterations=iter-1;return result;}
         const double beta=(rho!=0.0)?(rho_new/rho)*(alpha/omega):0.0;
         if(!std::isfinite(beta)){result.status=SolverStatus::DIVERGED;result.iterations=iter-1;return result;}
         for(std::size_t i=0;i<n;++i)w.p(i)=w.r(i)+beta*(w.p(i)-omega*w.v(i));
@@ -121,7 +124,7 @@
         if(!finite_vector(w.v)){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
 
         const double rv=krylov_dot(w.r_hat, w.v, redp, reduction);
-        if(!std::isfinite(rv)||std::abs(rv)<1e-30){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
+        if(!std::isfinite(rv)||std::abs(rv) <= 64*std::numeric_limits<double>::epsilon()*krylov_norm2(w.r_hat,SolverPrecision::FP64,reduction)*krylov_norm2(w.v,SolverPrecision::FP64,reduction)){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
         alpha=rho_new/rv;
         if(!std::isfinite(alpha)){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
         for(std::size_t i=0;i<n;++i)w.s(i)=w.r(i)-alpha*w.v(i);
@@ -164,8 +167,8 @@
 
         const double ts=krylov_dot(w.t, w.s, redp, reduction);
         const double tt=krylov_dot(w.t, w.t, redp, reduction);
-        if(!std::isfinite(ts)||!std::isfinite(tt)||tt<=1e-30){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
-        omega=ts/tt;if(!std::isfinite(omega)||std::abs(omega)<1e-30){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
+        if(!std::isfinite(ts)||!std::isfinite(tt)||!(tt>0.0)){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
+        omega=ts/tt;if(!std::isfinite(omega)||omega==0.0){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
         for(std::size_t i=0;i<n;++i){x(i)+=alpha*w.z(i)+omega*w.z_s(i);w.r(i)=w.s(i)-omega*w.t(i);}
         if(!finite_vector(x)||!finite_vector(w.r)){result.status=SolverStatus::DIVERGED;result.iterations=iter;return result;}
 
--- a/src/cfdx/core/linalg/cg_solver.h
+++ b/src/cfdx/core/linalg/cg_solver.h
@@ -92,8 +92,8 @@
         }
     if (null_space) null_space->remove(x);
     const double* Av = A.values_data();
-    const bool mp32 = precision.enabled && precision.operator_precision == SolverPrecision::FP32;
-    const SolverPrecision redp = precision.enabled ? precision.reduction_precision : SolverPrecision::FP64;
+    bool mp32 = precision.enabled && precision.operator_precision == SolverPrecision::FP32;
+    SolverPrecision redp = precision.enabled ? precision.reduction_precision : SolverPrecision::FP64;
     const auto* Ac = A.columns_data();
     const auto* Ar = A.row_offsets_data();
 
@@ -156,19 +156,20 @@
     double rsold = krylov_dot(r_vector, z_vector, redp, reduction);
 
     const double b_norm = krylov_norm2(b, SolverPrecision::FP64, reduction);
-    const double tol_abs = tolerance * std::max(b_norm, 1e-15);
-
-    if (!std::isfinite(rsold)) { result.status=SolverStatus::DIVERGED; return result; }
-
     const double initial_true_residual =
         krylov_norm2(r_vector, SolverPrecision::FP64, reduction);
+    const double rhs_scale=b_norm>0.0 ? b_norm : initial_true_residual;
+    const double tol_abs=tolerance*rhs_scale;
+    if (!std::isfinite(b_norm) || !std::isfinite(initial_true_residual) ||
+        !std::isfinite(tolerance) || !(tolerance>0.0) || !std::isfinite(rsold)) {
+        result.status=SolverStatus::DIVERGED; return result;
+    }
     if (preconditioner) {
         result.diagnostics.push_back(
             {0, initial_true_residual, initial_true_residual,
              rsold, 0.0, 0.0, 0.0});
     }
-    if ((preconditioner && initial_true_residual <= tol_abs) ||
-        (!preconditioner && rsold < tol_abs * tol_abs)) {
+    if (std::isfinite(initial_true_residual) && initial_true_residual <= tol_abs) {
         result.status = SolverStatus::CONVERGED;
         result.iterations = 0;
         result.residual = preconditioner
@@ -249,9 +250,8 @@
         // A preconditioned norm is not a true residual norm. Re-evaluate
         // b-A*x in FP64 for the external-preconditioner path before deciding
         // convergence, just as the reported final residual is evaluated.
-        const double res = preconditioner
-            ? mixed_precision_true_residual(A, b, x, rv)
-            : std::sqrt(std::abs(rsnew));
+        mixed_precision_true_residual(A, b, x, rv);
+        const double res = krylov_norm2(rv, SolverPrecision::FP64, reduction);
         if (preconditioner &&
             (iter == 1 || iter % 100 == 0 || iter == max_iter)) {
             const double recursive_residual =
@@ -262,7 +262,7 @@
                 {iter, res, recursive_residual, rsnew, pAp, alpha,
                  rsold != 0.0 ? rsnew / rsold : 0.0});
         }
-        if (res < tol_abs) {
+        if (res <= tol_abs) {
             result.status = SolverStatus::CONVERGED;
             result.iterations = iter;
             result.residual = res;
@@ -270,6 +270,20 @@
             return result;
         }
 
+        if (mp32 && (rsnew==0.0 ||
+                     krylov_norm2(r_vector,SolverPrecision::FP64,reduction)<=tol_abs)) {
+            // FP32 exhausted its attainable accuracy, but the original
+            // equation has not converged. Restart the remaining budget in FP64.
+            mp32=false;
+            redp=SolverPrecision::FP64;
+            for(std::size_t i=0;i<n;++i) r[i]=rv(i);
+            if(null_space) null_space->remove(r);
+            if(!apply_preconditioner()) {result.status=SolverStatus::NOT_APPLICABLE;return result;}
+            for(std::size_t i=0;i<n;++i) {r_vector(i)=r[i];z_vector(i)=z[i];}
+            rsold=krylov_dot(r_vector,z_vector,SolverPrecision::FP64,reduction);
+            p=z;
+            continue;
+        }
         // The recursive CG residual can drift away from b-Ax in finite
         // precision, especially once an effective AMG preconditioner has
         // reduced the true residual by many orders of magnitude. Continuing
--- a/src/cfdx/core/linalg/communication_avoiding.h
+++ b/src/cfdx/core/linalg/communication_avoiding.h
@@ -1,6 +1,8 @@
 #pragma once
 #include "cfdx/core/linalg/vector.h"
+#if defined(CFDX_HAS_MPI) && CFDX_HAS_MPI
 #include "cfdx/core/parallel/mpi_utils.h"
+#endif
 #include <algorithm>
 #include <cmath>
 #include <cstddef>
@@ -25,6 +27,7 @@
     return r;
 }
 
+#if defined(CFDX_HAS_MPI) && CFDX_HAS_MPI
 // Batch the two additive Krylov reductions into one MPI collective.
 // max_abs is reduced separately because it requires MPI_MAX rather than SUM.
 // This API is deliberately independent of a solver so distributed Krylov
@@ -72,4 +75,6 @@
     return mpi_fused_reduction(fused_reduction(a, b), comm);
 }
 
+#endif // CFDX_HAS_MPI
+
 } // namespace cfdx::core
--- a/src/cfdx/core/linalg/exact_schur.h
+++ b/src/cfdx/core/linalg/exact_schur.h
@@ -83,7 +83,7 @@
     bool apply(const Vector& rhs_p, Vector& pressure) const override {
         if (!blocks_ || !auu_solve_) return false;
         const std::size_t n = blocks_->pressure_size();
-        if (rhs_p.size() != n) return false;
+        if (rhs_p.size() != n || !rhs_p.is_valid()) return false;
         if (pressure_null_space_ && !pressure_null_space_->is_compatible(rhs_p)) return false;
 
         if (pressure.size() != n) pressure = Vector(n, 0.0);
@@ -96,7 +96,8 @@
         for (std::size_t i = 0; i < n; ++i) p(i) = r(i);
         double r_dot_r = r.dot(r);
         const double r0 = std::sqrt(r_dot_r);
-        if (!(r0 > 0.0)) return true;  // rhs is zero, pressure=0 is the solution
+        if (!std::isfinite(r0)) return false;
+        if (r0 == 0.0) return rhs_p.norm2() == 0.0;
 
         const double abs_tol = controls_.tolerance * r0;
 
--- a/src/cfdx/core/linalg/gauss_seidel_solver.h
+++ b/src/cfdx/core/linalg/gauss_seidel_solver.h
@@ -15,8 +15,8 @@
  for(std::size_t i=0;i<n;++i){bool f=false;for(std::size_t k=Ar[i];k<Ar[i+1];++k)if(Ac[k]==static_cast<std::uint32_t>(i)){d[i]=Av[k];f=true;break;}if(!f||d[i]==0||!std::isfinite(d[i])){r.status=SolverStatus::NOT_APPLICABLE;return r;}}
  for(std::size_t i=0;i<n;++i)if(!std::isfinite(b(i))||!std::isfinite(x(i))){r.status=SolverStatus::DIVERGED;return r;}
  const double bn=b.norm2(),tol=tolerance*std::max(bn,1e-15);
- auto residual=[&](){double s=0;for(std::size_t i=0;i<n;++i){double q=0;for(std::size_t k=Ar[i];k<Ar[i+1];++k)q+=Av[k]*x(Ac[k]);double v=q-b(i);if(!std::isfinite(v))return std::numeric_limits<double>::infinity();s+=v*v;}return std::sqrt(s);};
- double res=residual(); if(res<=tol){r.status=SolverStatus::CONVERGED;r.residual=res;r.residual_relative=bn>0?res/bn:0;return r;}
+ auto residual=[&](){double s=0;for(std::size_t i=0;i<n;++i){double q=0;for(std::size_t k=Ar[i];k<Ar[i+1];++k)q+=Av[k]*x(Ac[k]);double v=q-b(i);if(!std::isfinite(v))return std::numeric_limits<double>::infinity();s=std::hypot(s,v);}return s;};
+ double res=residual(); if(!std::isfinite(bn)||!std::isfinite(res)){r.status=SolverStatus::DIVERGED;return r;} if(res<=tol){r.status=SolverStatus::CONVERGED;r.residual=res;r.residual_relative=bn>0?res/bn:0;return r;}
  for(std::size_t it=1;it<=max_iter;++it){
   for(std::size_t i=0;i<n;++i){double rhs=b(i);for(std::size_t k=Ar[i];k<Ar[i+1];++k)if(Ac[k]!=static_cast<std::uint32_t>(i))rhs-=Av[k]*x(Ac[k]);x(i)=rhs/d[i];if(!std::isfinite(x(i))){r.status=SolverStatus::DIVERGED;r.iterations=it;return r;}}
   res=residual();if(!std::isfinite(res)){r.status=SolverStatus::DIVERGED;r.iterations=it;return r;}if(res<=tol){r.status=SolverStatus::CONVERGED;r.iterations=it;r.residual=res;r.residual_relative=bn>0?res/bn:0;return r;}
--- a/src/cfdx/core/linalg/gmres_solver.h
+++ b/src/cfdx/core/linalg/gmres_solver.h
@@ -20,6 +20,8 @@
 struct LinearOperator {
     std::size_t size = 0;
     std::function<void(const Vector&, Vector&)> apply;
+    // Optional original-precision action used only for certification.
+    std::function<void(const Vector&, Vector&)> apply_true;
 };
 
 inline SolverResult solve_gmres(
@@ -56,8 +58,8 @@
     // Scale the stopping criterion with the actual RHS norm. Using max(||b||,1)
     // makes physically small systems appear converged at x=0 and is incorrect
     // for nondimensionalized or otherwise small-scale operators.
-    const double rhs_scale = std::max(b_norm, 1e-300);
-    const double tol = tolerance * b_norm;
+    double rhs_scale = std::max(b_norm, 1e-300);
+    double tol = tolerance * b_norm;
 
     // Diagnostic: track the gauge-row value of the initial solution vector.
     // The gauge row is conventionally the last pressure unknown (index n-1),
@@ -80,18 +82,17 @@
         op.apply(w.vin, out);
     };
     auto true_residual = [&]() {
-        apply_operator(x.data(), w.ax);
+        if (op.apply_true) op.apply_true(x, w.ax);
+        else apply_operator(x.data(), w.ax);
         for (std::size_t i = 0; i < n; ++i) w.r(i) = b(i) - w.ax(i);
         return krylov_norm2(w.r, SolverPrecision::FP64, controls.reduction);
     };
 
     double beta = true_residual();
+    if (!std::isfinite(beta)) { result.status = SolverStatus::DIVERGED; return result; }
     if (b_norm == 0.0) {
-        result.status = beta == 0.0 ? SolverStatus::CONVERGED : SolverStatus::DIVERGED;
-        result.iterations = 0;
-        result.residual = beta;
-        result.residual_relative = beta;
-        return result;
+        rhs_scale = std::max(beta, 1e-300);
+        tol = tolerance * rhs_scale;
     }
     if (beta <= tol) {
         result.status = SolverStatus::CONVERGED;
@@ -152,14 +153,16 @@
             for (int i = 0; i <= j; ++i) {
                 const double* vi = w.v(static_cast<std::size_t>(i));
                 double h = 0.0;
-                for (std::size_t k = 0; k < n; ++k) h += vi[k] * w.w(k);
+                for (std::size_t k = std::min(controls.reduction.local_begin,n); k < controls.reduction.end_for(n); ++k) h += vi[k] * w.w(k);
+                h = krylov_sum(h, controls.reduction);
                 w.H(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) = h;
                 for (std::size_t k = 0; k < n; ++k) w.w(k) -= h * vi[k];
             }
             for (int i = 0; i <= j; ++i) {
                 const double* vi = w.v(static_cast<std::size_t>(i));
                 double correction = 0.0;
-                for (std::size_t k = 0; k < n; ++k) correction += vi[k] * w.w(k);
+                for (std::size_t k = std::min(controls.reduction.local_begin,n); k < controls.reduction.end_for(n); ++k) correction += vi[k] * w.w(k);
+                correction = krylov_sum(correction, controls.reduction);
                 w.H(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) += correction;
                 for (std::size_t k = 0; k < n; ++k) w.w(k) -= correction * vi[k];
             }
@@ -384,6 +387,10 @@
     bool setup_preconditioner = true) {
     LinearOperator op;
     op.size = A.n_rows();
+    op.apply_true = [&A](const Vector& in, Vector& out) {
+        if(out.size()!=A.n_rows()) out.resize(A.n_rows());
+        mixed_precision_matvec(A,in,out,SolverPrecision::FP64);
+    };
     op.apply = [&A, precision](const Vector& in, Vector& out) {
         if (out.size() != A.n_rows()) out.resize(A.n_rows());
         if (precision.enabled && precision.operator_precision == SolverPrecision::FP32) {
@@ -405,8 +412,20 @@
         result.status = SolverStatus::NOT_APPLICABLE;
         return result;
     }
-    return solve_fgmres(
-        op, b, x, restart, max_iter, tolerance, preconditioner, controls);
+    auto result = solve_fgmres(op, b, x, restart, max_iter, tolerance, preconditioner, controls);
+    if (precision.enabled && precision.operator_precision == SolverPrecision::FP32 &&
+        result.status != SolverStatus::CONVERGED && result.status != SolverStatus::NOT_APPLICABLE &&
+        result.iterations < max_iter && x.is_valid()) {
+        // Certification uses the original operator; finish in FP64 rather
+        // than accepting an FP32 residual or silently increasing the budget.
+        op.apply = op.apply_true;
+        op.apply_true = {};
+        auto refined = solve_fgmres(op, b, x, restart, max_iter-result.iterations,
+                                      tolerance, preconditioner, controls);
+        refined.iterations += result.iterations;
+        return refined;
+    }
+    return result;
 }
 
 inline SolverResult solve_fgmres(
@@ -439,6 +458,10 @@
     bool setup_preconditioner = true) {
     LinearOperator op;
     op.size = A.n_rows();
+    op.apply_true = [&A](const Vector& in, Vector& out) {
+        if(out.size()!=A.n_rows()) out.resize(A.n_rows());
+        mixed_precision_matvec(A,in,out,SolverPrecision::FP64);
+    };
     op.apply = [&A, precision](const Vector& in, Vector& out) {
         if (out.size() != A.n_rows()) out.resize(A.n_rows());
         if (precision.enabled && precision.operator_precision == SolverPrecision::FP32) { mixed_precision_matvec(A, in, out, SolverPrecision::FP32); return; }
@@ -457,7 +480,20 @@
         result.status = SolverStatus::NOT_APPLICABLE;
         return result;
     }
-    return solve_gmres(op, b, x, restart, max_iter, tolerance, preconditioner, controls);
+    auto result = solve_gmres(op, b, x, restart, max_iter, tolerance, preconditioner, controls);
+    if (precision.enabled && precision.operator_precision == SolverPrecision::FP32 &&
+        result.status != SolverStatus::CONVERGED && result.status != SolverStatus::NOT_APPLICABLE &&
+        result.iterations < max_iter && x.is_valid()) {
+        // Certification uses the original operator; finish in FP64 rather
+        // than accepting an FP32 residual or silently increasing the budget.
+        op.apply = op.apply_true;
+        op.apply_true = {};
+        auto refined = solve_gmres(op, b, x, restart, max_iter-result.iterations,
+                                      tolerance, preconditioner, controls);
+        refined.iterations += result.iterations;
+        return refined;
+    }
+    return result;
 }
 
 inline SolverResult solve_gmres(
--- a/src/cfdx/core/linalg/jacobi_solver.h
+++ b/src/cfdx/core/linalg/jacobi_solver.h
@@ -15,8 +15,8 @@
  for(std::size_t i=0;i<n;++i){bool f=false;for(std::size_t k=Ar[i];k<Ar[i+1];++k)if(Ac[k]==static_cast<std::uint32_t>(i)){d[i]=Av[k];f=true;break;}if(!f||d[i]==0||!std::isfinite(d[i])){r.status=SolverStatus::NOT_APPLICABLE;return r;}}
  for(std::size_t i=0;i<n;++i)if(!std::isfinite(b(i))||!std::isfinite(x(i))){r.status=SolverStatus::DIVERGED;return r;}
  const double bn=b.norm2(),tol=tolerance*std::max(bn,1e-15); Vector next(n,0.0);
- auto residual=[&](){double s=0;for(std::size_t i=0;i<n;++i){double q=0;for(std::size_t k=Ar[i];k<Ar[i+1];++k)q+=Av[k]*x(Ac[k]);double v=q-b(i);if(!std::isfinite(v))return std::numeric_limits<double>::infinity();s+=v*v;}return std::sqrt(s);};
- double res=residual(); if(res<=tol){r.status=SolverStatus::CONVERGED;r.residual=res;r.residual_relative=bn>0?res/bn:0;return r;}
+ auto residual=[&](){double s=0;for(std::size_t i=0;i<n;++i){double q=0;for(std::size_t k=Ar[i];k<Ar[i+1];++k)q+=Av[k]*x(Ac[k]);double v=q-b(i);if(!std::isfinite(v))return std::numeric_limits<double>::infinity();s=std::hypot(s,v);}return s;};
+ double res=residual(); if(!std::isfinite(bn)||!std::isfinite(res)){r.status=SolverStatus::DIVERGED;return r;} if(res<=tol){r.status=SolverStatus::CONVERGED;r.residual=res;r.residual_relative=bn>0?res/bn:0;return r;}
  for(std::size_t it=1;it<=max_iter;++it){for(std::size_t i=0;i<n;++i){double rhs=b(i);for(std::size_t k=Ar[i];k<Ar[i+1];++k)if(Ac[k]!=static_cast<std::uint32_t>(i))rhs-=Av[k]*x(Ac[k]);next(i)=rhs/d[i];if(!std::isfinite(next(i))){r.status=SolverStatus::DIVERGED;r.iterations=it;return r;}}for(std::size_t i=0;i<n;++i)x(i)=next(i);res=residual();if(!std::isfinite(res)){r.status=SolverStatus::DIVERGED;r.iterations=it;return r;}if(res<=tol){r.status=SolverStatus::CONVERGED;r.iterations=it;r.residual=res;r.residual_relative=bn>0?res/bn:0;return r;}}
  r.status=SolverStatus::MAX_ITER_REACHED;r.iterations=max_iter;r.residual=res;r.residual_relative=bn>0?res/bn:0;return r;
 }}
--- a/src/cfdx/core/linalg/krylov_reductions.h
+++ b/src/cfdx/core/linalg/krylov_reductions.h
@@ -77,16 +77,32 @@
                            const KrylovReductionPolicy& policy = {}) {
     const std::size_t begin = std::min(policy.local_begin, a.size());
     const std::size_t end = policy.end_for(a.size());
+    double scale = 0.0;
+    double invalid = 0.0;
+    for (std::size_t i = begin; i < end; ++i) {
+        if (!std::isfinite(a(i))) invalid = 1.0;
+        scale = std::max(scale, std::abs(a(i)));
+    }
+    if (krylov_sum(invalid, policy) > 0.0)
+        return std::numeric_limits<double>::quiet_NaN();
+#if CFDX_HAS_MPI
+    if (policy.mpi_enabled) {
+        double global_scale = 0.0;
+        MPI_Allreduce(&scale, &global_scale, 1, MPI_DOUBLE, MPI_MAX, policy.comm);
+        scale = global_scale;
+    }
+#endif
+    if (scale == 0.0) return 0.0;
     double local = 0.0;
     for (std::size_t i = begin; i < end; ++i) {
+        const double q = a(i) / scale;
         local += reduction_precision == SolverPrecision::FP32
-            ? static_cast<double>(static_cast<float>(a(i)) * static_cast<float>(a(i)))
-            : a(i) * a(i);
+            ? static_cast<double>(static_cast<float>(q) * static_cast<float>(q))
+            : q * q;
     }
     const double sum = krylov_sum(local, policy);
-    if (!std::isfinite(sum))
-        return std::numeric_limits<double>::quiet_NaN();
-    return std::sqrt(std::max(0.0, sum));
+    if (!std::isfinite(sum)) return std::numeric_limits<double>::quiet_NaN();
+    return scale * std::sqrt(std::max(0.0, sum));
 }
 
 } // namespace cfdx::core
--- a/src/cfdx/core/linalg/matrix_free_fv_operator.h
+++ b/src/cfdx/core/linalg/matrix_free_fv_operator.h
@@ -9,12 +9,12 @@
 namespace cfdx::core {
 class FvDiffusionOperator final:public LinearOperatorBase{
 public:
- FvDiffusionOperator(const Mesh& mesh,const std::vector<Vec3>& centres,const std::vector<Vec3>& Sf,double gamma):mesh_(mesh),cc_(centres),Sf_(Sf),gamma_(gamma){if(cc_.size()!=mesh_.n_cells()||Sf_.size()!=mesh_.n_faces()||gamma_<0.0)throw std::invalid_argument("invalid FvDiffusionOperator geometry");}
+ FvDiffusionOperator(const Mesh& mesh,const std::vector<Vec3>& centres,const std::vector<Vec3>& Sf,double gamma):mesh_(mesh),cc_(centres),Sf_(Sf),gamma_(gamma){if(cc_.size()!=mesh_.n_cells()||Sf_.size()!=mesh_.n_faces()||gamma_<0.0||!std::isfinite(gamma_))throw std::invalid_argument("invalid FvDiffusionOperator geometry");}
  FvDiffusionOperator(const Mesh& mesh,const GeometryCache& geometry,double gamma):FvDiffusionOperator(mesh,geometry.cell_centres,geometry.face_Sf,gamma){if(!is_valid(geometry,mesh))throw std::invalid_argument("invalid GeometryCache");}
  std::size_t rows()const noexcept override{return mesh_.n_cells();}
  std::size_t cols()const noexcept override{return mesh_.n_cells();}
  bool has_diagonal()const noexcept override{return true;}
- void diagonal(Vector& d)const override{d.resize(mesh_.n_cells());d.fill(0.0);const auto& own=mesh_.ownership();for(std::size_t f=0;f<mesh_.n_faces();++f){const std::size_t o=own.owner(f);const auto nr=own.neighbour(f);if(nr<0)continue;const std::size_t nb=static_cast<std::size_t>(nr);const double dist=(cc_[nb]-cc_[o]).mag();if(!(dist>0.0))throw std::runtime_error("degenerate diffusion face");const double a=gamma_*Sf_[f].mag()/dist;d(o)+=a;d(nb)+=a;}for(std::size_t i=0;i<d.size();++i)if(std::abs(d(i))<1e-30)d(i)=1.0;}
+ void diagonal(Vector& d)const override{d.resize(mesh_.n_cells());d.fill(0.0);const auto& own=mesh_.ownership();for(std::size_t f=0;f<mesh_.n_faces();++f){const std::size_t o=own.owner(f);const auto nr=own.neighbour(f);if(nr<0)continue;const std::size_t nb=static_cast<std::size_t>(nr);const double dist=(cc_[nb]-cc_[o]).mag();if(!(dist>0.0))throw std::runtime_error("degenerate diffusion face");const double a=gamma_*Sf_[f].mag()/dist;d(o)+=a;d(nb)+=a;}}
  void apply(const Vector& x,Vector& y)const override{const std::size_t n=mesh_.n_cells();if(x.size()!=n)throw std::invalid_argument("FvDiffusionOperator size mismatch");if(y.size()!=n)y.resize(n);y.fill(0.0);const auto& own=mesh_.ownership();for(std::size_t f=0;f<mesh_.n_faces();++f){const std::size_t o=own.owner(f);const auto nr=own.neighbour(f);const double area=Sf_[f].mag();if(!(area>0.0))continue;if(nr>=0){const std::size_t nb=static_cast<std::size_t>(nr);const double d=(cc_[nb]-cc_[o]).mag();if(!(d>0.0))throw std::runtime_error("degenerate diffusion face");const double a=gamma_*area/d;const double q=a*(x(o)-x(nb));y(o)+=q;y(nb)-=q;}}}
 private:const Mesh& mesh_;const std::vector<Vec3>& cc_;const std::vector<Vec3>& Sf_;double gamma_;
 };
--- a/src/cfdx/core/linalg/vector.h
+++ b/src/cfdx/core/linalg/vector.h
@@ -62,12 +62,15 @@
     }
     Value norm2() const {
         Value s = 0.0;
-        for (const auto v : data_) s += v * v;
-        return std::sqrt(s);
+        for (const auto v : data_) s = std::hypot(s, v);
+        return s;
     }
     Value norm_inf() const {
         Value s = 0.0;
-        for (const auto v : data_) s = std::max(s, std::abs(v));
+        for (const auto v : data_) {
+            if (std::isnan(v)) return v;
+            s = std::max(s, std::abs(v));
+        }
         return s;
     }
 
--- a/src/cfdx/core/numerics/gradient.h
+++ b/src/cfdx/core/numerics/gradient.h
@@ -177,6 +177,13 @@
                 if (boundary_policy == BoundaryGradientPolicy::REJECT_BOUNDARY_STENCIL)
                     throw std::runtime_error("compute_gradient_least_squares: boundary face in rejected stencil at cell "
                                              + std::to_string(c));
+                if (boundary_policy == BoundaryGradientPolicy::ZERO_GRADIENT_GHOST) {
+                    const Vec3 normal=geometry.face_normals[f];
+                    const double dn=(geometry.face_centres[f]-geometry.cell_centres[c]).dot(normal);
+                    if (!(dn>0.0) || !std::isfinite(dn)) throw std::runtime_error("invalid ghost distance");
+                    neighbours.push_back(geometry.cell_centres[c]+normal*(2.0*dn));
+                    neighbour_values.push_back(values[c]);
+                }
                 continue;
             }
             if (nb >= n_cells || nb == c)
@@ -435,12 +442,18 @@
         double a[9][9] = {};
         double b[9] = {};
         int npts = 0;
+        std::vector<Vec3> fit_centres;
+        std::vector<double> fit_values;
+        double h = 0.0;
         const auto add = [&](std::size_t nb) {
             if (nb == c || seen[nb]) return;
             seen[nb] = 1;
             const Vec3 d = geometry.cell_centres[nb] - P;
             const double r2 = d.mag2();
             if (!(r2 > 0.0) || !std::isfinite(r2)) return;
+            fit_centres.push_back(geometry.cell_centres[nb]);
+            fit_values.push_back(values[nb]);
+            h=std::max(h,std::sqrt(r2));
             const double w = 1.0 / r2;
             const double basis[9] = {d.x, d.y, d.z,
                                      d.x*d.x, d.y*d.y, d.z*d.z,
@@ -464,6 +477,11 @@
         double x[9] = {};
         bool have_quad = false;
         if (npts >= 9) {
+            const double scales[9]={h,h,h,h*h,h*h,h*h,h*h,h*h,h*h};
+            for(int i=0;i<9;++i) {
+                b[i]/=scales[i];
+                for(int j=0;j<9;++j) a[i][j]/=scales[i]*scales[j];
+            }
             // Solve the 9x9 normal equations with partial pivoting.
             double m[9][10] = {};
             for (int i = 0; i < 9; ++i) {
@@ -498,12 +516,13 @@
                 for (int j = i + 1; j < 9; ++j) rhs -= m[i][j] * x[j];
                 if (std::abs(m[i][i]) > rank_tol) x[i] = rhs / m[i][i];
             }
-            g = {x[0], x[1], x[2]};
-            have_quad = true;
-        } else {
-            g = least_squares_gradient(P, phi_c,
-                std::vector<Vec3>{}, std::vector<double>{});
-        }
+            have_quad = std::all_of(std::begin(pivoted),std::end(pivoted),[](bool p){return p;});
+            if(have_quad) {
+                for(int i=0;i<9;++i) x[i]/=scales[i];
+                g = {x[0], x[1], x[2]};
+            }
+        }
+        if(!have_quad) g=least_squares_gradient(P,phi_c,fit_centres,fit_values);
         if (!std::isfinite(g.x) || !std::isfinite(g.y) || !std::isfinite(g.z))
             throw std::runtime_error("compute_quadratic_fit_least_squares: non-finite gradient");
         fit(c, 0) = g.x;
@@ -839,6 +858,7 @@
                 if (boundary_policy == BoundaryGradientPolicy::REJECT_BOUNDARY_STENCIL)
                     throw std::runtime_error("compute_gradient_weighted_least_squares: boundary face in rejected stencil at cell "
                                              + std::to_string(c));
+                if (boundary_policy == BoundaryGradientPolicy::EXCLUDE_BOUNDARY) continue;
 
                 // ZERO_GRADIENT_GHOST is a real policy even when no explicit
                 // boundary-condition vector is supplied. In that case the
--- a/src/cfdx/core/numerics/interpolation.h
+++ b/src/cfdx/core/numerics/interpolation.h
@@ -200,17 +200,25 @@
             const double vo = cell_field.component_data(d)[owner];
             const double vn = cell_field.component_data(d)[neighbour];
             double v = vo;
+            const auto linear_value = [&]() {
+                const Vec3 delta=geometry.cell_centres[neighbour]-geometry.cell_centres[owner];
+                const double d2=delta.dot(delta);
+                if (!(d2>0.0) || !std::isfinite(d2))
+                    throw std::runtime_error("linear interpolation: degenerate centre distance");
+                const double w=(geometry.face_centres[f]-geometry.cell_centres[owner]).dot(delta)/d2;
+                return (1.0-w)*vo+w*vn;
+            };
 
             if (!internal) {
                 v = vo;
             } else if (scheme == InterpScheme::LINEAR) {
-                v = 0.5 * (vo + vn);
+                v = linear_value();
             } else if (scheme == InterpScheme::UPWIND) {
                 v = (flux[f] >= 0.0) ? vo : vn;
             } else if (scheme == InterpScheme::BLENDED) {
                 // Central/upwind blending: v = blend*linear + (1-blend)*upwind.
                 const double v_up = (flux[f] >= 0.0) ? vo : vn;
-                v = blend_linear * 0.5 * (vo + vn) + (1.0 - blend_linear) * v_up;
+                v = blend_linear * linear_value() + (1.0 - blend_linear) * v_up;
             } else if (scheme == InterpScheme::QUICK || scheme == InterpScheme::QUICK_BOUNDED) {
                 const bool owner_upwind = flux[f] >= 0.0;
                 const std::size_t up = owner_upwind ? owner : neighbour;
--- a/src/cfdx/core/numerics/newton.h
+++ b/src/cfdx/core/numerics/newton.h
@@ -138,6 +138,14 @@
             result.residual_relative = norm / scale;
             return result;
         }
+        const double accepted_norm=r.norm2();
+        if (std::isfinite(accepted_norm) && accepted_norm<=controls.tolerance*scale) {
+            result.status=SolverStatus::CONVERGED;
+            result.iterations=iteration+1;
+            result.residual=accepted_norm;
+            result.residual_relative=accepted_norm/scale;
+            return result;
+        }
     }
 
     result.status = SolverStatus::MAX_ITER_REACHED;
--- a/src/cfdx/core/numerics/temporal.h
+++ b/src/cfdx/core/numerics/temporal.h
@@ -182,6 +182,8 @@
                         candidate = (a0 * phi_data[cell] - a1 * prev[cell] +
                              b * rhs_data[cell]) / denom;
                     }
+                    if (!std::isfinite(candidate) || !std::isfinite(new_data[cell]))
+                        throw std::runtime_error("advance_time: non-finite implicit iterate");
                     max_delta = std::max(max_delta, std::abs(candidate - new_data[cell]));
                     new_data[cell] = candidate;
                 }
@@ -249,7 +251,9 @@
                         double* new_data = phi_new.component_data(d);
                         for (std::size_t cell = 0; cell < n_cells; ++cell) {
                             const double candidate = phi_data[cell] + dt * rhs_data[cell];
-                            max_delta = std::max(max_delta, std::abs(candidate - new_data[cell]));
+                            if (!std::isfinite(candidate) || !std::isfinite(new_data[cell]))
+                        throw std::runtime_error("advance_time: non-finite implicit iterate");
+                    max_delta = std::max(max_delta, std::abs(candidate - new_data[cell]));
                             new_data[cell] = candidate;
                         }
                     }
--- a/src/cfdx/physics/cht_solver.h
+++ b/src/cfdx/physics/cht_solver.h
@@ -15,6 +15,7 @@
 struct ChtInterfaceControls {
     std::string region1_patch;
     std::string region2_patch;
+    // Compatibility assertions only: energy controls own material properties.
     double conductivity1 = 1.0;
     double conductivity2 = 1.0;
     double tolerance = 1e-8;
@@ -23,6 +24,7 @@
     double matching_tolerance = 1e-8;
     double temperature_tolerance = 1e-8;
     double area_relative_tolerance = 1e-6;
+    double normal_alignment_tolerance = 1e-8;
 };
 
 struct ChtInterfaceFacePair {
@@ -37,13 +39,37 @@
     const cfdx::core::Mesh& mesh2, const FvGeometry& g2,
     const ChtInterfaceControls& c)
 {
+    const auto positive_finite=[](double x) { return x>0.0 && std::isfinite(x); };
     if(c.region1_patch.empty() || c.region2_patch.empty() ||
-       c.conductivity1<=0.0 || c.conductivity2<=0.0 ||
-       c.tolerance<=0.0 || c.max_iterations==0 ||
-       !(c.relaxation>0.0 && c.relaxation<=1.0) ||
-       c.matching_tolerance<=0.0 || c.temperature_tolerance<=0.0 ||
-       c.area_relative_tolerance<=0.0)
+       !positive_finite(c.conductivity1) || !positive_finite(c.conductivity2) ||
+       !positive_finite(c.tolerance) || c.max_iterations==0 ||
+       !positive_finite(c.relaxation) || c.relaxation>1.0 ||
+       !positive_finite(c.matching_tolerance) ||
+       !positive_finite(c.temperature_tolerance) ||
+       !positive_finite(c.area_relative_tolerance) ||
+       !positive_finite(c.normal_alignment_tolerance) ||
+       c.normal_alignment_tolerance>=1.0)
         throw std::invalid_argument("invalid CHT interface controls");
+    const auto finite_vec=[](const cfdx::core::Vec3& x) {
+        return std::isfinite(x.x) && std::isfinite(x.y) && std::isfinite(x.z);
+    };
+    const auto validate_geometry=[&](const cfdx::core::Mesh& mesh, const FvGeometry& g) {
+        if(g.face_centres.size()!=mesh.n_faces() ||
+           g.face_area_vectors.size()!=mesh.n_faces() ||
+           g.face_patch.size()!=mesh.n_faces() ||
+           g.cell_centres.size()!=mesh.n_cells() ||
+           g.cell_volumes.size()!=mesh.n_cells() ||
+           !mesh.faces().is_consistent() || !mesh.faces().indices_valid(mesh.n_points()))
+            throw std::invalid_argument("CHT geometry/connectivity size mismatch");
+        for(std::size_t cell=0;cell<mesh.n_cells();++cell)
+            if(!finite_vec(g.cell_centres[cell]) || !positive_finite(g.cell_volumes[cell]))
+                throw std::invalid_argument("CHT invalid cell geometry");
+        for(std::size_t f=0;f<mesh.n_faces();++f)
+            if(!finite_vec(g.face_centres[f]) || !finite_vec(g.face_area_vectors[f]))
+                throw std::invalid_argument("CHT nonfinite face geometry");
+    };
+    validate_geometry(mesh1,g1);
+    validate_geometry(mesh2,g2);
 
     std::size_t p1=mesh1.boundary().n_patches();
     std::size_t p2=mesh2.boundary().n_patches();
@@ -56,8 +82,28 @@
 
     const auto& f1=mesh1.boundary().patch(p1).face_ids;
     const auto& f2=mesh2.boundary().patch(p2).face_ids;
-    if(f1.size()!=f2.size())
-        throw std::invalid_argument("CHT interface face count mismatch");
+    if(f1.empty() || f1.size()!=f2.size())
+        throw std::invalid_argument("CHT interface must have matching nonempty face sets");
+    const auto validate_patch=[&](const cfdx::core::Mesh& mesh,
+                                  const FvGeometry& g, std::size_t patch) {
+        if(mesh.boundary().patch(patch).type==cfdx::core::PatchType::EMPTY)
+            throw std::invalid_argument("CHT interface cannot be an EMPTY patch");
+        std::vector<bool> seen(mesh.n_faces(),false);
+        for(const auto f:mesh.boundary().patch(patch).face_ids) {
+            if(f>=mesh.n_faces() || seen[f] || g.face_patch[f]!=patch ||
+               mesh.ownership().neighbour(f)>=0 || mesh.ownership().owner(f)>=mesh.n_cells())
+                throw std::invalid_argument("CHT invalid boundary face ownership or membership");
+            seen[f]=true;
+            const auto d=g.face_centres[f]-g.cell_centres[mesh.ownership().owner(f)];
+            const double area=g.face_area_vectors[f].mag();
+            const double projection=d.dot(g.face_area_vectors[f]);
+            // Do not use the transport kernel's degenerate-2D fallback for CHT.
+            if(!positive_finite(area) || !positive_finite(projection))
+                throw std::invalid_argument("CHT requires a positive outward normal distance");
+        }
+    };
+    validate_patch(mesh1,g1,p1);
+    validate_patch(mesh2,g2,p2);
 
     std::vector<bool> used(f2.size(),false);
     std::vector<ChtInterfaceFacePair> pairs;
@@ -78,6 +124,35 @@
         const double area_scale=std::max({a1,a2,1e-30});
         if(std::abs(a1-a2)/area_scale>c.area_relative_tolerance)
             throw std::runtime_error("CHT interface face area mismatch");
+        const auto& s1=g1.face_area_vectors[face1];
+        const cfdx::core::Vec3 n1(s1.x/a1,s1.y/a1,s1.z/a1);
+        const auto& s2=g2.face_area_vectors[f2[bestj]];
+        const cfdx::core::Vec3 n2(s2.x/a2,s2.y/a2,s2.z/a2);
+        if((n1+n2).mag()>c.normal_alignment_tolerance)
+            throw std::invalid_argument("CHT interface outward normals must oppose");
+        // A centre and an area do not establish a conformal interface. Check
+        // its vertex set as well; nonconformal interpolation is not implemented.
+        const auto face2=f2[bestj];
+        const auto count=mesh1.faces().face_size(face1);
+        if(count!=mesh2.faces().face_size(face2))
+            throw std::invalid_argument("CHT nonconformal interface vertex count");
+        std::vector<bool> vertex_used(count,false);
+        const auto off1=mesh1.faces().face_offset(face1);
+        const auto off2=mesh2.faces().face_offset(face2);
+        for(std::size_t v=0;v<count;++v) {
+            const auto v1=mesh1.faces().vertices_data()[off1+v];
+            const cfdx::core::Vec3 x1(mesh1.points().x(v1),mesh1.points().y(v1),mesh1.points().z(v1));
+            bool found=false;
+            for(std::size_t w=0;w<count;++w) {
+                if(vertex_used[w]) continue;
+                const auto v2=mesh2.faces().vertices_data()[off2+w];
+                const cfdx::core::Vec3 x2(mesh2.points().x(v2),mesh2.points().y(v2),mesh2.points().z(v2));
+                if(finite_vec(x1) && finite_vec(x2) && (x1-x2).mag()<=c.matching_tolerance) {
+                    vertex_used[w]=true;found=true;break;
+                }
+            }
+            if(!found) throw std::invalid_argument("CHT nonconformal interface vertices");
+        }
         used[bestj]=true;
         pairs.push_back({
             face1,f2[bestj],
@@ -93,6 +168,8 @@
     std::size_t iterations=0;
     double interface_imbalance=0.0;
     double interface_temperature_change=0.0;
+    // Boundary values actually used in the accepted region solves, pair order.
+    std::vector<double> interface_temperatures;
 };
 
 inline ChtSolveResult solve_two_region_cht(
@@ -110,8 +187,52 @@
     const ScalarBoundaryConditions& bcs1 = {},
     const ScalarBoundaryConditions& bcs2 = {})
 {
+    validate_energy_controls(energy1);
+    validate_energy_controls(energy2);
+    const auto same_property=[](double a, double b) {
+        return std::isfinite(a) && std::isfinite(b) &&
+            std::abs(a-b)<=64.0*std::numeric_limits<double>::epsilon()*
+                            std::max(std::abs(a),std::abs(b));
+    };
+    if(!(energy1.conductivity>0.0) || !(energy2.conductivity>0.0) ||
+       !same_property(controls.conductivity1,energy1.conductivity) ||
+       !same_property(controls.conductivity2,energy2.conductivity))
+        throw std::invalid_argument("CHT conductivity assertions disagree with energy materials");
+    if(!same_property(energy1.dt,energy2.dt))
+        throw std::invalid_argument("CHT requires the same physical time step in both regions");
+    const auto validate_fields=[](const cfdx::core::Mesh& mesh, const auto& flux,
+                                  const auto& temperature, const auto& source) {
+        if(flux.size()!=mesh.n_faces() || flux.dimension()!=1 ||
+           temperature.size()!=mesh.n_cells() || temperature.dimension()!=1 ||
+           source.size()!=mesh.n_cells() || source.dimension()!=1)
+            throw std::invalid_argument("CHT requires scalar fields matching each region");
+        for(std::size_t f=0;f<flux.size();++f)
+            if(!std::isfinite(flux(f))) throw std::invalid_argument("CHT nonfinite mass flux");
+        for(std::size_t cell=0;cell<temperature.size();++cell)
+            if(!std::isfinite(temperature(cell)) || !std::isfinite(source(cell)))
+                throw std::invalid_argument("CHT nonfinite temperature/source");
+    };
+    validate_fields(mesh1,flux1,T1,source1);
+    validate_fields(mesh2,flux2,T2,source2);
     const auto pairs=match_cht_interface(mesh1,g1,mesh2,g2,controls);
+    std::vector<double> conductance1(pairs.size()),conductance2(pairs.size());
+    for(std::size_t i=0;i<pairs.size();++i) {
+        const auto& p=pairs[i];
+        if(flux1(p.face1)!=0.0 || flux2(p.face2)!=0.0)
+            throw std::invalid_argument("CHT thermal interface must be impermeable");
+        // Same material, normal distance and area as assemble_energy_equation.
+        conductance1[i]=energy1.conductivity*g1.face_area_vectors[p.face1].mag()/
+            boundary_normal_distance(g1,p.face1,p.cell1);
+        conductance2[i]=energy2.conductivity*g2.face_area_vectors[p.face2].mag()/
+            boundary_normal_distance(g2,p.face2,p.cell2);
+        if(!(conductance1[i]>0.0) || !std::isfinite(conductance1[i]) ||
+           !(conductance2[i]>0.0) || !std::isfinite(conductance2[i]))
+            throw std::invalid_argument("CHT invalid interface conductance");
+    }
     ChtSolveResult result;
+    result.interface_temperatures.resize(pairs.size());
+    const auto T1_previous_time=T1;
+    const auto T2_previous_time=T2;
     std::vector<double> previous_interface_temperature(pairs.size(), std::numeric_limits<double>::quiet_NaN());
 
     for(std::size_t iter=1;iter<=controls.max_iterations;++iter) {
@@ -122,17 +243,16 @@
         // sides in the fixed-point limit.
         for(std::size_t i=0;i<pairs.size();++i) {
             const auto& p=pairs[i];
-            const double d1=(g1.face_centres[p.face1]-g1.cell_centres[p.cell1]).mag();
-            const double d2=(g2.face_centres[p.face2]-g2.cell_centres[p.cell2]).mag();
-            if(!(d1>0.0) || !(d2>0.0) || !std::isfinite(d1) || !std::isfinite(d2))
-                throw std::runtime_error("CHT interface face-to-cell distance is invalid");
-            const double h1=controls.conductivity1/d1;
-            const double h2=controls.conductivity2/d2;
-            const double Tint_new=(h1*T1(p.cell1)+h2*T2(p.cell2))/(h1+h2);
+            const double scale=std::max(conductance1[i],conductance2[i]);
+            const double h1=conductance1[i]/scale, h2=conductance2[i]/scale;
+            const double weight=h1/(h1+h2);
+            const double Tint_new=weight*T1(p.cell1)+(1.0-weight)*T2(p.cell2);
             double Tint=Tint_new;
             if(std::isfinite(previous_interface_temperature[i]))
                 Tint=controls.relaxation*Tint_new +
                      (1.0-controls.relaxation)*previous_interface_temperature[i];
+            if(!std::isfinite(Tint)) throw std::runtime_error("CHT nonfinite interface temperature");
+            result.interface_temperatures[i]=Tint;
             auto& values1=fv1.values[controls.region1_patch];
             auto& values2=fv2.values[controls.region2_patch];
             if(values1.size()!=mesh1.n_faces()) values1.resize(mesh1.n_faces(),0.0);
@@ -141,21 +261,20 @@
             values2[p.face2]=Tint;
         }
 
-        auto r1=solve_energy(mesh1,g1,flux1,T1,source1,energy1,bcs1,&fv1);
-        auto r2=solve_energy(mesh2,g2,flux2,T2,source2,energy2,bcs2,&fv2);
+        auto r1=solve_energy(mesh1,g1,flux1,T1,source1,energy1,bcs1,&fv1,nullptr,&T1_previous_time);
+        auto r2=solve_energy(mesh2,g2,flux2,T2,source2,energy2,bcs2,&fv2,nullptr,&T2_previous_time);
         if(!r1.converged || !r2.converged)
             throw std::runtime_error("CHT region energy solve did not converge");
 
         double qimb=0.0, qscale=1e-30, dtint=0.0;
         for(std::size_t i=0;i<pairs.size();++i) {
             const auto& p=pairs[i];
-            const double d1=(g1.face_centres[p.face1]-g1.cell_centres[p.cell1]).mag();
-            const double d2=(g2.face_centres[p.face2]-g2.cell_centres[p.cell2]).mag();
             const double Tint=fv1.values.at(controls.region1_patch)[p.face1];
-            const double qflux1=controls.conductivity1*
-                (T1(p.cell1)-Tint)/d1;
-            const double qflux2=controls.conductivity2*
-                (Tint-T2(p.cell2))/d2;
+            // Compare actual heat rates (W), not a second fictitious k/d model.
+            const double qflux1=conductance1[i]*(T1(p.cell1)-Tint);
+            const double qflux2=conductance2[i]*(Tint-T2(p.cell2));
+            if(!std::isfinite(qflux1) || !std::isfinite(qflux2))
+                throw std::runtime_error("CHT nonfinite interface heat rate");
             qimb=std::max(qimb,std::abs(qflux1-qflux2));
             qscale=std::max(qscale,std::abs(qflux1));
             qscale=std::max(qscale,std::abs(qflux2));
--- a/src/cfdx/physics/dual_time_driver.h
+++ b/src/cfdx/physics/dual_time_driver.h
@@ -20,7 +20,8 @@
     double dt,
     const DualTimeStepControls& controls,
     const std::function<std::vector<double>(const std::vector<double>&)>& rhs,
-    const DualTimeUpdate& update)
+    const DualTimeUpdate& update,
+    double dt_previous = 0.0)
 {
     validate_dual_time_controls(controls);
     if (initial_state.size() != state_n.size())
@@ -36,7 +37,7 @@
 
     for (std::size_t iteration = 1; iteration <= controls.max_pseudo_iterations; ++iteration) {
         const auto residual = dual_time_physical_residual(
-            state, state_n, state_nm1, dt, controls.physical_scheme, rhs);
+            state, state_n, state_nm1, dt, controls.physical_scheme, rhs, dt_previous);
         const double norm = dual_time_norm(residual);
         if (iteration == 1) report.initial_residual = norm;
         const double relative = report.initial_residual > 0.0
@@ -59,7 +60,7 @@
 
         const double next_residual = dual_time_norm(
             dual_time_physical_residual(candidate, state_n, state_nm1, dt,
-                                         controls.physical_scheme, rhs));
+                                         controls.physical_scheme, rhs, dt_previous));
         state = candidate;
         if (std::isfinite(previous_residual)) {
             if (next_residual < previous_residual)
--- a/src/cfdx/physics/dual_time_field_integrator.h
+++ b/src/cfdx/physics/dual_time_field_integrator.h
@@ -150,7 +150,7 @@
 
         auto solved = run_dual_time_step(
             initial_values, n_values, nm1_values, dt,
-            controls.pseudo_time, rhs_vector, update);
+            controls.pseudo_time, rhs_vector, update, history.dt_previous);
         return {
             unflatten(initial, solved.first, "dual_time_state"),
             std::move(solved.second)
--- a/src/cfdx/physics/dual_time_stepping.h
+++ b/src/cfdx/physics/dual_time_stepping.h
@@ -70,13 +70,13 @@
 }
 
 inline double dual_time_norm(const std::vector<double>& r) {
-    long double sum = 0.0L;
+    double sum = 0.0;
     for (double value : r) {
         if (!std::isfinite(value))
             throw std::runtime_error("dual-time residual is non-finite");
-        sum += static_cast<long double>(value) * static_cast<long double>(value);
+        sum=std::hypot(sum,value);
     }
-    return std::sqrt(static_cast<double>(sum));
+    return sum;
 }
 
 inline std::vector<double> dual_time_physical_residual(
@@ -85,7 +85,8 @@
     const std::vector<double>& state_nm1,
     double dt,
     DualTimePhysicalScheme scheme,
-    const std::function<std::vector<double>(const std::vector<double>&)>& rhs)
+    const std::function<std::vector<double>(const std::vector<double>&)>& rhs,
+    double dt_previous = 0.0)
 {
     if (!(dt > 0.0) || !std::isfinite(dt))
         throw std::invalid_argument("dual-time physical dt must be finite and positive");
@@ -104,9 +105,14 @@
         for (std::size_t i = 0; i < state.size(); ++i)
             residual[i] = state[i] - state_n[i] - dt * f[i];
     } else {
+        if (!std::isfinite(dt_previous) || dt_previous<0.0)
+            throw std::invalid_argument("invalid previous physical dt");
+        const double ratio=dt_previous>0.0 ? dt/dt_previous : 1.0;
+        const double a0=(1.0+2.0*ratio)/(1.0+ratio);
+        const double a1=-(1.0+ratio);
+        const double a2=ratio*ratio/(1.0+ratio);
         for (std::size_t i = 0; i < state.size(); ++i)
-            residual[i] = 1.5 * state[i] - 2.0 * state_n[i]
-                        + 0.5 * state_nm1[i] - dt * f[i];
+            residual[i] = a0*state[i]+a1*state_n[i]+a2*state_nm1[i]-dt*f[i];
     }
     return residual;
 }
--- a/src/cfdx/physics/energy_solver.h
+++ b/src/cfdx/physics/energy_solver.h
@@ -205,7 +205,7 @@
         const double sign=(Sf0.dot(dvec)>=0.0)?1.0:-1.0;
         const auto Sf=Sf0*sign;
         const double area=Sf.mag();
-        const double d=dvec.mag();
+        const double d=boundary_normal_distance(geometry,f,o);
         const double F=mass_flux(f);
         double Tf=temperature(o);
         bool face_override_valid=false;
--- a/src/cfdx/physics/finite_volume_transport.h
+++ b/src/cfdx/physics/finite_volume_transport.h
@@ -320,9 +320,13 @@
                 const double gn = (*cell_diffusion)[n];
                 if (go < 0.0 || gn < 0.0)
                     throw std::invalid_argument("assemble_scalar_equation: negative cell diffusion");
+                const auto delta=geometry.cell_centres[n]-geometry.cell_centres[o];
+                const double dp=(geometry.face_centres[f]-geometry.cell_centres[o]).dot(delta)/distance;
+                const double dn=distance-dp;
+                if (!(dp>0.0) || !(dn>0.0) || !std::isfinite(go) || !std::isfinite(gn))
+                    throw std::invalid_argument("invalid material face distances or coefficients");
                 gamma_face = (go > 0.0 && gn > 0.0)
-                    ? 2.0 * go * gn / (go + gn)
-                    : 0.0;
+                    ? distance/(dp/go+dn/gn) : 0.0;
             }
             const double D = gamma_face * area / distance;
 
@@ -360,7 +364,7 @@
                         gz[upwind] * (Cf.z - C_up.z);
                     phi_high = phi_up + high_increment;
                     // Keep the historical bounded SOU reconstruction.
-                    phi_high = std::clamp(
+                    if (bounded_convection) phi_high = std::clamp(
                         phi_high, std::min(phi_up, phi_other), std::max(phi_up, phi_other));
                 } else if (convection_scheme == ConvectionScheme::TVD) {
                     // Darwish-Moukalled successive-slope ratio for an
@@ -390,7 +394,7 @@
                     const double limiter =
                         cfdx::core::limiter_psi(ratio, cfdx::core::LimiterType::MINMOD);
                     phi_high = phi_up + limiter * high_increment;
-                    phi_high = std::clamp(
+                    if (bounded_convection) phi_high = std::clamp(
                         phi_high, std::min(phi_up, phi_other), std::max(phi_up, phi_other));
                 } else {
                     // Central, blended and quadratic schemes are explicit
@@ -532,9 +536,9 @@
     eq.min_face_area = min_face_area;
     eq.max_face_area = max_face_area;
 
-    double total_flux = 0.0;
-    for (std::size_t c = 0; c < nc; ++c) total_flux += div_phi[c];
-    eq.max_imbalance = std::abs(total_flux);
+    eq.max_imbalance = 0.0;
+    for (std::size_t c = 0; c < nc; ++c)
+        eq.max_imbalance=std::max(eq.max_imbalance,std::abs(div_phi[c]));
     return eq;
 }
 
@@ -614,16 +618,17 @@
                 const double go=(*cell_diffusion)[o], gn=(*cell_diffusion)[n];
                 if (!(go >= 0.0) || !(gn >= 0.0) || !std::isfinite(go) || !std::isfinite(gn))
                     throw std::invalid_argument("reconstruct_scalar_transport_flux: invalid cell diffusion");
-                gamma = (go > 0.0 && gn > 0.0) ? 2.0*go*gn/(go+gn) : 0.0;
+                const auto delta=geometry.cell_centres[n]-geometry.cell_centres[o];
+                const double dp=(geometry.face_centres[f]-geometry.cell_centres[o]).dot(delta)/d;
+                const double dn=d-dp;
+                if (!(dp>0.0) || !(dn>0.0))
+                    throw std::invalid_argument("invalid material face distances");
+                gamma = (go > 0.0 && gn > 0.0) ? d/(dp/go+dn/gn) : 0.0;
             }
             const double up = F >= 0.0 ? field(o) : field(n);
             psi_face = up;
             if (convection_scheme == ConvectionScheme::SECOND_ORDER_UPWIND ||
-                convection_scheme == ConvectionScheme::TVD ||
-                convection_scheme == ConvectionScheme::CENTRAL ||
-                convection_scheme == ConvectionScheme::BLENDED ||
-                convection_scheme == ConvectionScheme::QUICK ||
-                convection_scheme == ConvectionScheme::QUICK_BOUNDED) {
+                convection_scheme == ConvectionScheme::TVD) {
                 const std::size_t u = F >= 0.0 ? o : n;
                 const std::size_t dcell = F >= 0.0 ? n : o;
                 const double gx=grad.component_data(0)[u];
@@ -648,6 +653,11 @@
                 } else {
                     psi_face = reconstructed_face(f);
                 }
+            } else if (convection_scheme == ConvectionScheme::CENTRAL ||
+                       convection_scheme == ConvectionScheme::BLENDED ||
+                       convection_scheme == ConvectionScheme::QUICK ||
+                       convection_scheme == ConvectionScheme::QUICK_BOUNDED) {
+                psi_face = reconstructed_face(f);
             }
             result(f)=F*psi_face - gamma*area/d*(field(n)-field(o));
         } else {
@@ -685,9 +695,11 @@
                 throw std::invalid_argument("reconstruct_scalar_transport_flux: invalid boundary diffusion");
             if (bc.type == ScalarBoundaryType::FIXED_VALUE) {
                 psi_face=bc.value;
-                result(f)=F*psi_face - gamma*area/d*(psi_face-field(o));
+                result(f)=F*(F>=0.0 ? field(o) : psi_face) - gamma*area/d*(psi_face-field(o));
             } else if (bc.type == ScalarBoundaryType::FIXED_GRADIENT) {
                 psi_face=field(o)+bc.gradient*d;
+                // FIXED_GRADIENT assembly uses the extrapolated boundary
+                // state for both flux signs (unlike Dirichlet upwinding).
                 result(f)=F*psi_face - gamma*area*bc.gradient;
             } else {
                 // Zero-gradient: no diffusive boundary flux and psi_f=psi_P.
@@ -721,24 +733,27 @@
         const double candidate_value = equation.rhs(0) / diagonal;
         const double candidate_residual =
             std::abs(diagonal * candidate_value - equation.rhs(0));
-        if (candidate_residual <= controls.tolerance) {
-            // A scalar 1x1 system has an exact closed-form solution. Do not
-            // under-relax an already converged linear solve: doing so makes
-            // the returned state fail the requested residual even though the
-            // reported linear residual is zero.
+        // This tolerance has the same relative meaning as the Krylov paths.
+        // An absolute comparison rejected correctly rounded large systems.
+        const double rhs_scale = std::abs(equation.rhs(0));
+        const double relative = rhs_scale > 0.0
+            ? candidate_residual / rhs_scale : candidate_residual;
+        if (!std::isfinite(candidate_value) || !std::isfinite(candidate_residual) ||
+            !std::isfinite(relative)) {
+            return {cfdx::core::SolverStatus::DIVERGED, 1,
+                    candidate_residual, relative};
+        }
+        if (relative <= controls.tolerance) {
+            // Preserve the established exact 1x1 predictor contract.
             solution(0) = candidate_value;
-            return {
-                cfdx::core::SolverStatus::CONVERGED,
-                1, candidate_residual, candidate_residual
-            };
-        }
-        solution(0) += controls.relaxation * (candidate_value - solution(0));
+            return {cfdx::core::SolverStatus::CONVERGED, 1,
+                    candidate_residual, relative};
+        }
+        // A rejected predictor must not replace the accepted iterate.
         const double residual =
             std::abs(diagonal * solution(0) - equation.rhs(0));
-        return {
-            cfdx::core::SolverStatus::MAX_ITER_REACHED,
-            controls.max_iterations, residual, residual
-        };
+        return {cfdx::core::SolverStatus::MAX_ITER_REACHED, 1,
+                residual, rhs_scale > 0.0 ? residual / rhs_scale : residual};
     }
 
     cfdx::core::Vector candidate = solution;
--- a/src/cfdx/physics/low_storage_time_integration.h
+++ b/src/cfdx/physics/low_storage_time_integration.h
@@ -6,6 +6,17 @@
 
 namespace cfdx::physics {
 
+inline void validate_rk_rhs(
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>& u,
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k) {
+    if(u.size()!=k.size() || u.dimension()!=k.dimension())
+        throw std::invalid_argument("RK RHS shape mismatch");
+    for(std::size_t d=0;d<u.dimension();++d)
+        for(std::size_t i=0;i<u.size();++i)
+            if(!std::isfinite(u(i,d)) || !std::isfinite(k(i,d)))
+                throw std::runtime_error("RK non-finite input or RHS");
+}
+
 template<class RHS>
 inline void ssprk3_step(
     cfdx::core::Field<double,cfdx::core::Location::CELL>& u,
@@ -13,19 +24,25 @@
     if (!(dt > 0.0) || !std::isfinite(dt))
         throw std::invalid_argument("invalid dt");
     const std::size_t n = u.size();
-    cfdx::core::Field<double,cfdx::core::Location::CELL> k(n, "rhs");
+    cfdx::core::Field<double,cfdx::core::Location::CELL> k(n, "rhs", "", u.dimension());
     cfdx::core::Field<double,cfdx::core::Location::CELL> u0 = u;
 
     rhs(u, k);
-    for (std::size_t i = 0; i < n; ++i) u(i) = u0(i) + dt * k(i);
+    validate_rk_rhs(u,k);
+    for (std::size_t d = 0; d < u.dimension(); ++d)
+        for (std::size_t i = 0; i < n; ++i) u(i,d) = u0(i,d) + dt * k(i,d);
 
     rhs(u, k);
-    for (std::size_t i = 0; i < n; ++i)
-        u(i) = 0.75 * u0(i) + 0.25 * (u(i) + dt * k(i));
+    validate_rk_rhs(u,k);
+    for (std::size_t d = 0; d < u.dimension(); ++d)
+        for (std::size_t i = 0; i < n; ++i)
+        u(i,d) = 0.75 * u0(i,d) + 0.25 * (u(i,d) + dt * k(i,d));
 
     rhs(u, k);
-    for (std::size_t i = 0; i < n; ++i)
-        u(i) = (1.0 / 3.0) * u0(i) + (2.0 / 3.0) * (u(i) + dt * k(i));
+    validate_rk_rhs(u,k);
+    for (std::size_t d = 0; d < u.dimension(); ++d)
+        for (std::size_t i = 0; i < n; ++i)
+        u(i,d) = (1.0 / 3.0) * u0(i,d) + (2.0 / 3.0) * (u(i,d) + dt * k(i,d));
 }
 
 // Williamson 3-stage third-order low-storage RK. This is genuinely 2-register:
@@ -37,18 +54,20 @@
     if (!(dt > 0.0) || !std::isfinite(dt))
         throw std::invalid_argument("invalid dt");
     const std::size_t n = u.size();
-    cfdx::core::Field<double,cfdx::core::Location::CELL> r(n, "rk3_register");
+    cfdx::core::Field<double,cfdx::core::Location::CELL> r(n, "rk3_register", "", u.dimension());
     r.fill(0.0);
 
     constexpr double alpha[3] = {0.0, -5.0 / 9.0, -153.0 / 128.0};
     constexpr double beta[3]  = {1.0 / 3.0, 15.0 / 16.0, 8.0 / 15.0};
 
-    cfdx::core::Field<double,cfdx::core::Location::CELL> k(n, "rhs");
+    cfdx::core::Field<double,cfdx::core::Location::CELL> k(n, "rhs", "", u.dimension());
     for (int s = 0; s < 3; ++s) {
         rhs(u, k);
+    validate_rk_rhs(u,k);
+        for (std::size_t d = 0; d < u.dimension(); ++d)
         for (std::size_t i = 0; i < n; ++i) {
-            r(i) = alpha[s] * r(i) + dt * k(i);
-            u(i) += beta[s] * r(i);
+            r(i,d) = alpha[s] * r(i,d) + dt * k(i,d);
+            u(i,d) += beta[s] * r(i,d);
         }
     }
 }
@@ -66,11 +85,15 @@
         throw std::invalid_argument("invalid dt");
     const std::size_t n = u.size();
     cfdx::core::Field<double,cfdx::core::Location::CELL> u0 = u;
-    cfdx::core::Field<double,cfdx::core::Location::CELL> k(n, "rhs");
+    cfdx::core::Field<double,cfdx::core::Location::CELL> k(n, "rhs", "", u.dimension());
     rhs(u, k);
-    for (std::size_t i = 0; i < n; ++i) u(i) = u0(i) + 0.5 * dt * k(i);
+    validate_rk_rhs(u,k);
+    for (std::size_t d = 0; d < u.dimension(); ++d)
+        for (std::size_t i = 0; i < n; ++i) u(i,d) = u0(i,d) + 0.5 * dt * k(i,d);
     rhs(u, k);
-    for (std::size_t i = 0; i < n; ++i) u(i) = u0(i) + dt * k(i);
+    validate_rk_rhs(u,k);
+    for (std::size_t d = 0; d < u.dimension(); ++d)
+        for (std::size_t i = 0; i < n; ++i) u(i,d) = u0(i,d) + dt * k(i,d);
 }
 
 } // namespace cfdx::physics
--- a/src/cfdx/physics/m1_m4_models.h
+++ b/src/cfdx/physics/m1_m4_models.h
@@ -277,8 +277,8 @@
 {
     if (absorption < 0.0 || mean_intensity < 0.0 || T < 0.0)
         throw std::invalid_argument("invalid P1 source input");
-    // P1 source convention: S_r = 4*kappa*(I_b - J).
-    return 4.0 * absorption * (p1_blackbody_intensity(T) - mean_intensity);
+    // J is mean intensity; G=4*pi*J, so emission-absorption is 4*pi*kappa*(Ib-J).
+    return 4.0 * M_PI * absorption * (p1_blackbody_intensity(T) - mean_intensity);
 }
 
 struct Direction {
--- a/src/cfdx/physics/sst_solver.h
+++ b/src/cfdx/physics/sst_solver.h
@@ -58,7 +58,8 @@
     s.sk=Pk;
     s.spk=-c.density*c.beta_star*wi;
     const double cross=2.0*(1.0-f1)*c.density*c.sst_sigma_w2/wi*grad_k_dot_grad_omega;
-    s.sw=blend(c.gamma1,c.gamma2)*c.density*S*S+std::max(cross,0.0);
+    // SST-2003 uses the SAME limited production in both equations.
+    s.sw=blend(c.gamma1,c.gamma2)*(s.nut>0.0 ? Pk/s.nut : 0.0)+std::max(cross,0.0);
     s.spw=-c.density*blend(c.beta1,c.beta2)*wi+std::min(cross,0.0)/wi;
     s.gamma_k=c.density*(c.molecular_viscosity+
         blend(c.sst_sigma_k1,c.sst_sigma_k2)*s.nut);
--- a/src/cfdx/physics/thermal_zones.h
+++ b/src/cfdx/physics/thermal_zones.h
@@ -19,7 +19,7 @@
 //
 // Discretization: standard two-point finite volumes. Internal faces use the
 // exact series-resistance conductance
-//     G = A / ( d/(2*k_owner) + R_c + d/(2*k_neighbour) )
+//     G = A / ( d_owner/k_owner + R_c + d_neighbour/k_neighbour )
 // which reduces to the harmonic mean when R_c = 0 and is exact for linear
 // temperature profiles in uniform cells. Boundary Dirichlet faces use
 // G_b = k_owner * A / d_bc. Boundary fixed-gradient faces inject
@@ -236,8 +236,12 @@
             const double k_o = zones.cell_conductivity[o];
             const double k_n = zones.cell_conductivity[n];
             const double r_c = zones.face_contact_resistance[f];
-            const double g = area /
-                (distance / (2.0 * k_o) + r_c + distance / (2.0 * k_n));
+            const auto delta=geometry.cell_centres[n]-geometry.cell_centres[o];
+            const double dp=(geometry.face_centres[f]-geometry.cell_centres[o]).dot(delta)/distance;
+            const double dn=distance-dp;
+            if (!(dp>0.0) || !(dn>0.0))
+                throw std::runtime_error("thermal conduction: invalid face distances");
+            const double g = area / (dp/k_o + r_c + dn/k_n);
             if (!(g > 0.0) || !std::isfinite(g))
                 throw std::runtime_error(
                     "assemble_thermal_conduction: non-positive face conductance");
--- a/src/cfdx/physics/thermophysical_models.h
+++ b/src/cfdx/physics/thermophysical_models.h
@@ -2,6 +2,7 @@
 #include <algorithm>
 #include <cmath>
 #include <cstddef>
+#include <functional>
 #include <stdexcept>
 #include <string>
 #include <vector>
@@ -114,80 +115,76 @@
     double enthalpy(double T,double Tref=300,double intervals=64) const {
         if(!std::isfinite(T)||!std::isfinite(Tref)||T<=0||Tref<=0||intervals<=0)
             throw std::invalid_argument("enthalpy: invalid inputs");
+        if (!std::isfinite(intervals) || intervals<1.0 || intervals>1000000.0)
+            throw std::invalid_argument("enthalpy: invalid quadrature interval count");
         const double a=std::min(T,Tref), b=std::max(T,Tref);
+        (void)cp(a); (void)cp(b); // Validate property and table extrapolation.
+        if (a==b) return 0.0;
         double integral=0.0;
-        if (heat_capacity.model == ScalarPropertyModel::CONSTANT) {
-            integral=heat_capacity.reference_value*(b-a);
-        } else if (heat_capacity.model == ScalarPropertyModel::LINEAR) {
-            const double c0=heat_capacity.reference_value;
-            const double c1=heat_capacity.coefficients.empty()?0.0:heat_capacity.coefficients[0];
-            const double xa=a-heat_capacity.reference_temperature, xb=b-heat_capacity.reference_temperature;
-            integral=c0*(b-a)+0.5*c1*(xb*xb-xa*xa);
-        } else if (heat_capacity.model == ScalarPropertyModel::POLYNOMIAL) {
-            if (heat_capacity.coefficients.empty())
-                throw std::invalid_argument("enthalpy: polynomial heat capacity requires coefficients");
-                integral=0.0;
-            for (std::size_t i=0;i<heat_capacity.coefficients.size();++i) {
-                const double ai=heat_capacity.coefficients[i];
-                if(!std::isfinite(ai)) throw std::invalid_argument("enthalpy: non-finite heat-capacity coefficient");
-                const double power=static_cast<double>(i+1);
-                integral += ai*(std::pow(b-heat_capacity.reference_temperature,power)
-                              - std::pow(a-heat_capacity.reference_temperature,power))/power;
+        const auto model=heat_capacity.model;
+        if (model==ScalarPropertyModel::CONSTANT) {
+            integral=cp(a)*(b-a);
+        } else if (model==ScalarPropertyModel::TABLE || model==ScalarPropertyModel::LINEAR) {
+            // Every segment is affine. Insert table knots and each bound
+            // crossing; trapezoids then integrate the ACTUAL cp exactly.
+            std::vector<double> knots{a,b};
+            if(model==ScalarPropertyModel::TABLE) {
+                heat_capacity.table.validate("heat capacity table");
+                for(double x:heat_capacity.table.temperature)
+                    if(x>a && x<b) knots.push_back(x);
             }
-        } else if (heat_capacity.model == ScalarPropertyModel::TABLE) {
-            heat_capacity.table.validate("heat capacity table");
-            const auto& tab=heat_capacity.table;
-            auto integrate_in_range = [&](double lo,double hi) {
-                if (hi<=lo) return 0.0;
-                double result=0.0;
-                double x0=lo, y0=tab.evaluate(lo,"heat capacity table");
-                for (std::size_t i=0;i+1<tab.temperature.size();++i) {
-                    const double l=std::max(lo,tab.temperature[i]);
-                    const double r=std::min(hi,tab.temperature[i+1]);
-                    if (r<=l) continue;
-                    const double yl=tab.value[i] + (l-tab.temperature[i]) *
-                        (tab.value[i+1]-tab.value[i])/(tab.temperature[i+1]-tab.temperature[i]);
-                    const double yr=tab.value[i] + (r-tab.temperature[i]) *
-                        (tab.value[i+1]-tab.value[i])/(tab.temperature[i+1]-tab.temperature[i]);
-                    result += 0.5*(yl+yr)*(r-l);
-                    x0=r; y0=yr;
-                }
-                return result;
+            std::sort(knots.begin(),knots.end());
+            const auto raw=[&](double x) {
+                if(model==ScalarPropertyModel::TABLE)
+                    return heat_capacity.table.evaluate(x,"heat capacity table");
+                return heat_capacity.reference_value+
+                    (heat_capacity.coefficients.empty()?0.0:heat_capacity.coefficients[0])*
+                    (x-heat_capacity.reference_temperature);
             };
-            if (a < tab.temperature.front() || b > tab.temperature.back()) {
-                if (tab.extrapolation == ExtrapolationPolicy::REJECT)
-                    throw std::out_of_range("enthalpy: heat-capacity table outside range");
-                const double lo=tab.temperature.front(), hi=tab.temperature.back();
-                integral=integrate_in_range(std::max(a,lo),std::min(b,hi));
-                if (a<lo) {
-                    if (tab.extrapolation == ExtrapolationPolicy::CLAMP)
-                        integral += tab.value.front()*(lo-a);
-                    else {
-                        const double slope=(tab.value[1]-tab.value[0])/(tab.temperature[1]-tab.temperature[0]);
-                        const double d=lo-a;
-                        integral += tab.value.front()*d - 0.5*slope*d*d;
+            for(std::size_t i=1;i<knots.size();++i) {
+                const double l=knots[i-1],r=knots[i],vl=raw(l),vr=raw(r);
+                std::vector<double> cuts{l,r};
+                if(heat_capacity.enforce_bounds && vr!=vl) {
+                    for(double bound:{heat_capacity.minimum,heat_capacity.maximum}) {
+                        const double w=(bound-vl)/(vr-vl);
+                        if(w>0.0 && w<1.0) cuts.push_back(l+w*(r-l));
                     }
                 }
-                if (b>hi) {
-                    if (tab.extrapolation == ExtrapolationPolicy::CLAMP)
-                        integral += tab.value.back()*(b-hi);
-                    else {
-                        const double slope=(tab.value.back()-tab.value[tab.value.size()-2])/
-                                           (tab.temperature.back()-tab.temperature[tab.temperature.size()-2]);
-                        const double d=b-hi;
-                        integral += tab.value.back()*d + 0.5*slope*d*d;
-                    }
-                }
-            } else {
-                integral=integrate_in_range(a,b);
+                std::sort(cuts.begin(),cuts.end());
+                for(std::size_t j=1;j<cuts.size();++j)
+                    integral+=0.5*(cp(cuts[j-1])+cp(cuts[j]))*(cuts[j]-cuts[j-1]);
+            }
+        } else if(model==ScalarPropertyModel::POLYNOMIAL && !heat_capacity.enforce_bounds) {
+            for(std::size_t i=0;i<heat_capacity.coefficients.size();++i) {
+                const double power=static_cast<double>(i+1);
+                integral+=heat_capacity.coefficients[i]*
+                    (std::pow(b-heat_capacity.reference_temperature,power)-
+                     std::pow(a-heat_capacity.reference_temperature,power))/power;
             }
         } else {
+            // Nonlinear/bounded laws: adaptive Simpson integration of cp(),
+            // not of the unbounded raw polynomial. Exhaustion is an error.
+            std::function<double(double,double,double,double,double,double,double,int)> adapt;
+            adapt=[&](double l,double r,double fl,double fm,double fr,double whole,double tol,int depth) {
+                const double m=(l+r)*0.5;
+                const double f1=cp((l+m)*0.5), f2=cp((m+r)*0.5);
+                const double left=(m-l)*(fl+4*f1+fm)/6;
+                const double right=(r-m)*(fm+4*f2+fr)/6;
+                const double delta=left+right-whole;
+                if(std::abs(delta)<=15*tol) return left+right+delta/15;
+                if(depth==0) throw std::runtime_error("enthalpy quadrature did not converge");
+                return adapt(l,m,fl,f1,fm,left,tol/2,depth-1)+
+                       adapt(m,r,fm,f2,fr,right,tol/2,depth-1);
+            };
             const std::size_t n=static_cast<std::size_t>(intervals);
-            const double h=(b-a)/n;
-            double s=cp(a)+cp(b);
-            for(std::size_t i=1;i<n;++i) s+=2*cp(a+h*i);
-            integral=0.5*h*s;
+            for(std::size_t i=0;i<n;++i) {
+                const double l=a+(b-a)*i/n,r=a+(b-a)*(i+1)/n;
+                const double fl=cp(l),fm=cp((l+r)*0.5),fr=cp(r);
+                const double whole=(r-l)*(fl+4*fm+fr)/6;
+                integral+=adapt(l,r,fl,fm,fr,whole,1e-11*std::max(1.0,std::abs(whole)),20);
+            }
         }
+        if (!std::isfinite(integral)) throw std::runtime_error("enthalpy integral non-finite");
         return T>=Tref?integral:-integral;
     }
 };
--- a/src/cfdx/physics/transport_models.h
+++ b/src/cfdx/physics/transport_models.h
@@ -76,7 +76,7 @@
                                              double S = SUTHERTY_S,
                                              double Pr = 0.71,
                                              double Sc = 0.7,
-                                             double Cp = AIR_R * 1005.0 / 1.4) {
+                                             double Cp = 1.4 * AIR_R / (1.4 - 1.0)) {
     TransportProperties tp;
     tp.mu = sutherland_viscosity(T, mu0_ref, T_ref, S);
     tp.k = prandtl_conductivity(tp.mu, Cp, Pr);
--- a/src/cfdx/physics/turbulence_solver.h
+++ b/src/cfdx/physics/turbulence_solver.h
@@ -428,9 +428,12 @@
     const double ki=std::max(k,c.k_min), wi=std::max(omega,c.omega_min);
     const double S=std::max(strain,0.0);
     KOmegaCellSources q;
-    // Baseline Wilcox k-omega uses nu_t = k/omega. The SST-specific
-    // strain-rate limiter belongs only to the SST closure.
-    q.nut=ki/wi;
+    // Incompressible Wilcox2006-noPope: strain is sqrt(2*Sij*Sij).
+    // The diffusivities below deliberately retain k/omega (not limited nut).
+    if (!(c.beta_star>0.0) || !(c.komega_clim>=0.0))
+        throw std::invalid_argument("invalid Wilcox limiter coefficients");
+    const double omega_hat=std::max(wi,c.komega_clim*S/std::sqrt(c.beta_star));
+    q.nut=ki/omega_hat;
     const double Pk=c.density*q.nut*S*S;
     q.sk=Pk; q.spk=-c.density*c.beta_star*wi;
     const double sigma_d=grad_k_dot_grad_omega>0.0?c.komega_sigma_d0:0.0;
--- a/src/cfdx/thermodynamics/equation_of_state.h
+++ b/src/cfdx/thermodynamics/equation_of_state.h
@@ -187,8 +187,8 @@
         return params_.gamma * p / density(p, T);
     }
 
-    double dp_dT_rho(double /*p*/, double T) const override {
-        return density(0, T) * R;
+    double dp_dT_rho(double p, double T) const override {
+        return density(p, T) * R;
     }
 
     double drho_dp_T(double p, double T) const override {
--- /dev/null
+++ b/tests/audit/affine_mesh_fixture.h
@@ -0,0 +1,174 @@
+// Issue #461 — N2 gradient/reconstruction V&V campaign.
+//
+// Closes the "complete skew/non-orthogonal/high-aspect-ratio/polyhedral V&V"
+// gap for the gradient family by running the same gradient operators over an
+// orthogonal, an affine-sheared (non-orthogonal) and a stretched
+// (high-aspect-ratio) mesh family:
+//
+//   - constant field  -> exact zero gradient (Green-Gauss and least-squares);
+//   - linear field    -> exact interior gradient on every affine family;
+//   - quadratic/smooth fields on a refinement sequence -> measured observed
+//     order on interior cells with L1/L2/Linf reporting.
+//
+// Polyhedral remains out of scope here (see the #461 audit); this test makes
+// the skew/non-orthogonal/high-aspect-ratio part explicit rather than relying
+// on the Laplacian-only skew campaign that was previously cited as N2 evidence.
+
+#include "cfdx/core/numerics/gradient.h"
+#include "cfdx/core/mesh/mesh.h"
+#include "cfdx/core/geometry/geometry_cache.h"
+#include "cfdx/core/field/field.h"
+#include "verification_metrics.h"
+
+#include <array>
+#include <cmath>
+#include <cstddef>
+#include <cstdint>
+#include <iomanip>
+#include <iostream>
+#include <stdexcept>
+#include <string>
+#include <unordered_map>
+#include <vector>
+
+using namespace cfdx::core;
+using namespace cfdx::verification;
+
+namespace {
+
+struct Grid {
+    Mesh mesh;
+    GeometryCache geometry;
+    std::vector<std::size_t> interior;
+};
+
+// Unit-cube Cartesian topology, then an affine map applied to every point:
+//   x' = x + shear * y,  y' = stretch * y,  z' = z
+// Affine transforms preserve the closed-cell geometric identity that makes
+// Green-Gauss linear-exact, while shear introduces non-orthogonality and
+// stretch introduces high aspect ratio.
+Grid make_affine_cube(std::size_t n, double shear, double stretch)
+{
+    if (n < 2) throw std::invalid_argument("make_affine_cube: n must be >= 2");
+    if (!(stretch > 0.0)) throw std::invalid_argument("make_affine_cube: stretch must be positive");
+
+    const double h = 1.0 / static_cast<double>(n);
+    const std::size_t nc = n * n * n;
+    Mesh m;
+    m.points().resize(8 * nc);
+
+    const auto cell_id = [n](std::size_t i, std::size_t j, std::size_t k) {
+        return (k * n + j) * n + i;
+    };
+    const auto transform = [shear, stretch](double x, double y, double z) {
+        return std::array<double, 3>{x + shear * y, stretch * y, z};
+    };
+
+    for (std::size_t k = 0; k < n; ++k) {
+        for (std::size_t j = 0; j < n; ++j) {
+            for (std::size_t i = 0; i < n; ++i) {
+                const std::size_t c = cell_id(i, j, k);
+                const std::size_t b = 8 * c;
+                const double x0 = i * h, x1 = (i + 1) * h;
+                const double y0 = j * h, y1 = (j + 1) * h;
+                const double z0 = k * h, z1 = (k + 1) * h;
+                const double raw[8][3] = {
+                    {x0,y0,z0}, {x1,y0,z0}, {x1,y1,z0}, {x0,y1,z0},
+                    {x0,y0,z1}, {x1,y0,z1}, {x1,y1,z1}, {x0,y1,z1}
+                };
+                for (std::size_t q = 0; q < 8; ++q) {
+                    const auto p = transform(raw[q][0], raw[q][1], raw[q][2]);
+                    m.points().set(b + q, p[0], p[1], p[2]);
+                }
+            }
+        }
+    }
+
+    std::vector<std::vector<std::size_t>> faces_per_cell(nc);
+
+    struct FaceKey {
+        std::size_t i, j, k, axis;
+        bool operator==(const FaceKey& other) const {
+            return i == other.i && j == other.j && k == other.k && axis == other.axis;
+        }
+    };
+    struct FaceKeyHash {
+        std::size_t operator()(const FaceKey& f) const {
+            std::size_t h = f.i;
+            h = h * 1315423911u + f.j;
+            h = h * 1315423911u + f.k;
+            h = h * 1315423911u + f.axis;
+            return h;
+        }
+    };
+    std::unordered_map<FaceKey, std::size_t, FaceKeyHash> face_map;
+
+    auto add_or_get_face = [&](FaceKey key, std::initializer_list<std::size_t> vertices) {
+        auto it = face_map.find(key);
+        if (it != face_map.end()) return it->second;
+        const std::size_t id = m.faces().n_faces();
+        m.faces().push_face(std::vector<FaceIndex>(vertices.begin(), vertices.end()));
+        face_map.emplace(key, id);
+        return id;
+    };
+
+    const auto point = [](std::size_t c, std::size_t q) { return 8 * c + q; };
+
+    for (std::size_t k = 0; k < n; ++k) {
+        for (std::size_t j = 0; j < n; ++j) {
+            for (std::size_t i = 0; i < n; ++i) {
+                const std::size_t c = cell_id(i, j, k);
+                auto face = [&](std::size_t axis, std::size_t fi, std::size_t fj, std::size_t fk,
+                                std::initializer_list<std::size_t> v) {
+                    return add_or_get_face({fi, fj, fk, axis}, v);
+                };
+                const std::size_t xm = face(0, i, j, k, {point(c,0),point(c,4),point(c,7),point(c,3)});
+                const std::size_t xp = face(0, i + 1, j, k, {point(c,1),point(c,2),point(c,6),point(c,5)});
+                const std::size_t ym = face(1, i, j, k, {point(c,0),point(c,1),point(c,5),point(c,4)});
+                const std::size_t yp = face(1, i, j + 1, k, {point(c,3),point(c,7),point(c,6),point(c,2)});
+                const std::size_t zm = face(2, i, j, k, {point(c,0),point(c,3),point(c,2),point(c,1)});
+                const std::size_t zp = face(2, i, j, k + 1, {point(c,4),point(c,5),point(c,6),point(c,7)});
+                faces_per_cell[c] = {xm, xp, ym, yp, zm, zp};
+            }
+        }
+    }
+
+    m.ownership().resize(m.n_faces());
+    for (const auto& kv : face_map) {
+        const FaceKey key = kv.first;
+        const std::size_t f = kv.second;
+        const bool lower_boundary = (key.axis == 0 ? key.i == 0 :
+                                     key.axis == 1 ? key.j == 0 : key.k == 0);
+        const bool upper_boundary = (key.axis == 0 ? key.i == n :
+                                     key.axis == 1 ? key.j == n : key.k == n);
+        const std::size_t owner_i = key.axis == 0 ? (lower_boundary ? 0 : key.i - 1) : key.i;
+        const std::size_t owner_j = key.axis == 1 ? (lower_boundary ? 0 : key.j - 1) : key.j;
+        const std::size_t owner_k = key.axis == 2 ? (lower_boundary ? 0 : key.k - 1) : key.k;
+        const std::size_t owner = cell_id(owner_i, owner_j, owner_k);
+        m.ownership().set_owner(f, owner);
+        if (lower_boundary || upper_boundary) {
+            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
+        } else {
+            std::size_t ni = owner_i, nj = owner_j, nk = owner_k;
+            if (key.axis == 0) ++ni;
+            else if (key.axis == 1) ++nj;
+            else ++nk;
+            m.ownership().set_neighbour(f, static_cast<std::int64_t>(cell_id(ni, nj, nk)));
+        }
+    }
+
+    for (std::size_t c = 0; c < nc; ++c) m.cells().push_cell(faces_per_cell[c]);
+
+    Grid grid;
+    grid.mesh = std::move(m);
+    grid.geometry = make_geometry_cache(grid.mesh);
+    for (std::size_t k = 0; k < n; ++k)
+        for (std::size_t j = 0; j < n; ++j)
+            for (std::size_t i = 0; i < n; ++i)
+                if (i >= 1 && j >= 1 && k >= 1 && i + 1 < n && j + 1 < n && k + 1 < n)
+                    grid.interior.push_back(cell_id(i, j, k));
+    return grid;
+}
+
+
+} // namespace
--- /dev/null
+++ b/tests/audit/chain_fixture.h
@@ -0,0 +1,74 @@
+#pragma once
+#include "cfdx/physics/finite_volume_transport.h"
+using namespace cfdx::core;
+static Mesh make_1d_chain(std::size_t n)
+{
+    if (n == 0) throw std::invalid_argument("make_1d_chain: n must be positive");
+    Mesh m;
+    m.points().resize(8 * n);
+    const double dx = 1.0 / static_cast<double>(n);
+    for (std::size_t i = 0; i < n; ++i) {
+        const double x0 = dx * static_cast<double>(i);
+        const double x1 = dx * static_cast<double>(i + 1);
+        const std::size_t b = 8 * i;
+        const double p[8][3] = {
+            {x0,0,0}, {x1,0,0}, {x1,1,0}, {x0,1,0},
+            {x0,0,1}, {x1,0,1}, {x1,1,1}, {x0,1,1}
+        };
+        for (std::size_t q = 0; q < 8; ++q)
+            m.points().set(b + q, p[q][0], p[q][1], p[q][2]);
+    }
+    std::vector<std::size_t> left, right, y0, y1, z0, z1, internal;
+    auto add_face = [&](std::initializer_list<std::size_t> v) {
+        const std::size_t id = m.faces().n_faces();
+        m.faces().push_face(std::vector<FaceIndex>(v.begin(), v.end()));
+        return id;
+    };
+    left.push_back(add_face({0,4,7,3}));
+    right.push_back(add_face({8*(n-1)+1,8*(n-1)+2,8*(n-1)+6,8*(n-1)+5}));
+    for (std::size_t i = 0; i < n; ++i) {
+        const std::size_t b = 8 * i;
+        y0.push_back(add_face({b,b+1,b+5,b+4}));
+        y1.push_back(add_face({b+3,b+7,b+6,b+2}));
+        z0.push_back(add_face({b,b+3,b+2,b+1}));
+        z1.push_back(add_face({b+4,b+5,b+6,b+7}));
+    }
+    for (std::size_t i = 0; i + 1 < n; ++i) {
+        const std::size_t b = 8 * i;
+        internal.push_back(add_face({b+1,b+2,b+6,b+5}));
+    }
+    m.ownership().resize(m.n_faces());
+    m.ownership().set_owner(left[0], 0);
+    m.ownership().set_neighbour(left[0], FaceOwnership::BOUNDARY);
+    m.ownership().set_owner(right[0], n - 1);
+    m.ownership().set_neighbour(right[0], FaceOwnership::BOUNDARY);
+    for (std::size_t i = 0; i < n; ++i) {
+        for (const auto f : {y0[i], y1[i], z0[i], z1[i]}) {
+            m.ownership().set_owner(f, i);
+            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
+        }
+    }
+    for (std::size_t i = 0; i + 1 < n; ++i) {
+        m.ownership().set_owner(internal[i], i);
+        m.ownership().set_neighbour(internal[i], static_cast<std::int64_t>(i + 1));
+    }
+    for (std::size_t i = 0; i < n; ++i)
+        m.cells().push_cell({i == 0 ? left[0] : internal[i-1],
+                             i + 1 == n ? right[0] : internal[i],
+                             y0[i], y1[i], z0[i], z1[i]});
+
+    Patch inlet; inlet.name = "inlet"; inlet.type = PatchType::WALL;
+    inlet.face_ids = left;
+    Patch outlet; outlet.name = "outlet"; outlet.type = PatchType::WALL;
+    outlet.face_ids = right;
+    Patch wall; wall.name = "wall"; wall.type = PatchType::WALL;
+    for (const auto f : y0) wall.face_ids.push_back(f);
+    for (const auto f : y1) wall.face_ids.push_back(f);
+    for (const auto f : z0) wall.face_ids.push_back(f);
+    for (const auto f : z1) wall.face_ids.push_back(f);
+    m.boundary().add_patch(inlet);
+    m.boundary().add_patch(outlet);
+    m.boundary().add_patch(wall);
+    return m;
+}
+
--- /dev/null
+++ b/tests/audit/test_audit_cht_contract.cpp
@@ -0,0 +1,132 @@
+// F28: independent two-capacity oracle and invalid-interface rejection.
+// No assertion disappears with NDEBUG. Validity checks must precede mutation.
+#include "chain_fixture.h"
+#include "cfdx/physics/cht_solver.h"
+#include <functional>
+#include <iomanip>
+#include <iostream>
+#include <limits>
+#include <sstream>
+#include <string>
+#include <vector>
+using namespace cfdx::core;
+using namespace cfdx::physics;
+struct Quiet {
+    std::ostringstream stream;
+    std::streambuf* previous;
+    std::streambuf* previous_error;
+    Quiet():previous(std::cout.rdbuf(stream.rdbuf())),previous_error(std::cerr.rdbuf(stream.rdbuf())){}
+    ~Quiet(){std::cout.rdbuf(previous);std::cerr.rdbuf(previous_error);}
+};
+struct Fixture {
+    Mesh a=make_1d_chain(1), b=make_1d_chain(1);
+    FvGeometry ga,gb;
+    Field<double,Location::FACE> fa{a.n_faces(),"F","kg/s",1},fb{b.n_faces(),"F","kg/s",1};
+    Field<double,Location::CELL> ta{1,"T","K",1},tb{1,"T","K",1},qa{1,"Q","W/m3",1},qb{1,"Q","W/m3",1};
+    EnergySolverControls ea,eb;
+    ChtInterfaceControls c;
+    Fixture(double l1=1,double l2=1,double shear1=0,double shear2=0,
+            double dt=.1,double relaxation=1,double k1=1,double k2=4) {
+        for(std::size_t i=0;i<a.n_points();++i) {
+            const double x=l1*a.points().x(i);
+            a.points().set(i,x,a.points().y(i)+shear1*(x-l1),a.points().z(i));
+        }
+        for(std::size_t i=0;i<b.n_points();++i) {
+            const double x=l1+l2*b.points().x(i);
+            b.points().set(i,x,b.points().y(i)+shear2*(x-l1),b.points().z(i));
+        }
+        ga=build_fv_geometry(a);gb=build_fv_geometry(b);
+        fa.fill(0);fb.fill(0);ta.fill(300);tb.fill(400);qa.fill(10);qb.fill(-5);
+        ea.density=2;eb.density=3;ea.cp=eb.cp=1;ea.dt=eb.dt=dt;
+        ea.conductivity=k1;eb.conductivity=k2;ea.relaxation=eb.relaxation=1;
+        ea.tolerance=eb.tolerance=1e-12;ea.max_iterations=eb.max_iterations=30;
+        c.region1_patch="outlet";c.region2_patch="inlet";
+        c.conductivity1=k1;c.conductivity2=k2;c.relaxation=relaxation;
+        c.tolerance=1e-10;c.temperature_tolerance=1e-10;c.max_iterations=1500;
+    }
+    ChtSolveResult solve() {
+        Quiet quiet;
+        return solve_two_region_cht(a,ga,b,gb,fa,fb,ta,tb,qa,qb,ea,eb,c);
+    }
+};
+int main() {
+    int checks=0,failures=0;
+    double max_error=0,max_normalized_balance=0;
+    auto check=[&](const std::string& name,bool ok,double value=0.) {
+        ++checks;failures+=!ok;
+        std::cout<<(ok?"PASS":"FAIL")<<'\t'<<name<<'\t'<<std::setprecision(17)<<value<<'\n';
+    };
+    int sample=0;
+    auto positive=[&](double l1,double l2,double s1,double s2,double dt,double relaxation,double k1,double k2) {
+        Fixture f(l1,l2,s1,s2,dt,relaxation,k1,k2);
+        const long double m1=2*l1/dt,m2=3*l2/dt;
+        const long double g=1.L/(l1/(2.L*k1)+l2/(2.L*k2));
+        const long double rhs1=m1*300+10*l1,rhs2=m2*400-5*l2;
+        const long double det=m1*m2+g*(m1+m2);
+        const double exact1=static_cast<double>(((m2+g)*rhs1+g*rhs2)/det);
+        const double exact2=static_cast<double>(((m1+g)*rhs2+g*rhs1)/det);
+        const std::string name="case"+std::to_string(++sample)+"/";
+        try {
+            const auto r=f.solve();
+            const double e1=f.ta(0)-exact1,e2=f.tb(0)-exact2;
+            const double balance=static_cast<double>(m1)*(f.ta(0)-300)+
+                static_cast<double>(m2)*(f.tb(0)-400)-10*l1+5*l2;
+            const double scale=std::max({1.,std::abs(static_cast<double>(g)*(f.ta(0)-f.tb(0))),10*l1,5*l2});
+            max_error=std::max({max_error,std::abs(e1),std::abs(e2)});
+            max_normalized_balance=std::max(max_normalized_balance,std::abs(balance)/scale);
+            check(name+"converged",r.converged,static_cast<double>(r.iterations));
+            check(name+"T1",std::isfinite(e1)&&std::abs(e1)<1e-8,e1);
+            check(name+"T2",std::isfinite(e2)&&std::abs(e2)<1e-8,e2);
+            check(name+"global_energy",std::isfinite(balance)&&std::abs(balance)/scale<1e-9,balance/scale);
+            const double q1=10*l1-static_cast<double>(m1)*(f.ta(0)-300);
+            const double qseries=static_cast<double>(g)*(f.ta(0)-f.tb(0));
+            check(name+"series_flux",std::abs(q1-qseries)/scale<1e-9,(q1-qseries)/scale);
+            std::cout<<"MEASURE\t"<<sample<<"\tl1="<<l1<<"\tl2="<<l2<<"\tshear1="<<s1
+                     <<"\tshear2="<<s2<<"\tdt="<<dt<<"\trelax="<<relaxation<<"\tT1="<<f.ta(0)
+                     <<"\tT2="<<f.tb(0)<<"\texact1="<<exact1<<"\texact2="<<exact2<<'\n';
+        } catch(const std::exception& e) { check(name+"exception",false);std::cout<<e.what()<<'\n'; }
+    };
+    for(const double dt:{.01,.1}) for(const double relaxation:{.5,1.}) {
+        positive(1,1,0,0,dt,relaxation,1,4);
+        positive(.4,1.6,2,-.5,dt,relaxation,1,4);
+        positive(1,2,3,0,dt,relaxation,1,4);
+    }
+    positive(1,1,2,-.5,.001,1,.01,10);
+    using Change=std::pair<std::string,std::function<void(Fixture&)>>;
+    const double nan=std::numeric_limits<double>::quiet_NaN();
+    const double inf=std::numeric_limits<double>::infinity();
+    const std::vector<Change> bad={
+        {"material1_mismatch",[](auto& f){f.c.conductivity1=2;}},
+        {"material2_mismatch",[](auto& f){f.c.conductivity2=1;}},
+        {"different_dt",[](auto& f){f.eb.dt*=.5;}},
+        {"mixed_steady_transient",[](auto& f){f.eb.dt=0;}},
+        {"nan_conductivity",[&](auto& f){f.c.conductivity1=nan;}},
+        {"inf_conductivity",[&](auto& f){f.c.conductivity2=inf;}},
+        {"nan_balance_tolerance",[&](auto& f){f.c.tolerance=nan;}},
+        {"nan_temperature_tolerance",[&](auto& f){f.c.temperature_tolerance=nan;}},
+        {"inf_matching_tolerance",[&](auto& f){f.c.matching_tolerance=inf;}},
+        {"nan_area_tolerance",[&](auto& f){f.c.area_relative_tolerance=nan;}},
+        {"permeable_interface",[](auto& f){f.fa(1)=.1;}},
+        {"inward_normal",[](auto& f){f.gb.face_area_vectors[0]=f.gb.face_area_vectors[0]*(-1);}},
+        {"overlapping_regions",[](auto& f){f.b=f.a;f.gb=f.ga;f.c.region2_patch="outlet";}},
+        {"nonconformal_vertices",[](auto& f){
+            for(std::size_t i=0;i<f.b.n_points();++i)
+                f.b.points().set(i,f.b.points().x(i),.5+2*(f.b.points().y(i)-.5),.5+.5*(f.b.points().z(i)-.5));
+            f.gb=build_fv_geometry(f.b);
+        }},
+        {"empty_interface",[](auto& f){f.a.boundary().patch(1).face_ids.clear();f.b.boundary().patch(0).face_ids.clear();}},
+        {"duplicate_interface",[](auto& f){f.a.boundary().patch(1).face_ids.push_back(1);f.b.boundary().patch(0).face_ids.push_back(0);}},
+        {"nan_cell_centre",[&](auto& f){f.ga.cell_centres[0].x=nan;}},
+        {"nan_source",[&](auto& f){f.qa(0)=nan;}}
+    };
+    for(const auto& item:bad) {
+        Fixture f;item.second(f);bool rejected=false;
+        try{f.solve();}catch(const std::invalid_argument&){rejected=true;}
+        catch(const std::exception& e){std::cout<<"OTHER_EXCEPTION\t"<<item.first<<'\t'<<e.what()<<'\n';}
+        check("reject/"+item.first,rejected);
+        check("unchanged/"+item.first,f.ta(0)==300 && f.tb(0)==400);
+    }
+    std::cout<<"SUMMARY checks="<<checks<<" failures="<<failures<<" max_T_error="<<max_error
+             <<" max_energy_normalized="<<max_normalized_balance<<'\n';
+    return failures?1:0;
+}
--- /dev/null
+++ b/tests/audit/test_audit_cht_time.cpp
@@ -0,0 +1,54 @@
+// Two independent thermal capacities coupled by an analytic series resistance.
+// Vary the physical step and interface relaxation without changing the oracle.
+#include "chain_fixture.h"
+#include "cfdx/physics/cht_solver.h"
+#include <cmath>
+#include <iostream>
+#include <iomanip>
+using namespace cfdx::core;
+using namespace cfdx::physics;
+int main() {
+    int checks=0, failures=0;
+    for (const double dt : {0.01, 0.1, 1.0}) {
+        for (const double relaxation : {0.5, 1.0}) {
+            auto a=make_1d_chain(1), b=make_1d_chain(1);
+            for(std::size_t i=0;i<b.n_points();++i)
+                b.points().set(i,b.points().x(i)+1,b.points().y(i),b.points().z(i));
+            auto ga=build_fv_geometry(a), gb=build_fv_geometry(b);
+            Field<double,Location::FACE> fa(a.n_faces(),"F","",1),fb(b.n_faces(),"F","",1);
+            fa.fill(0); fb.fill(0);
+            Field<double,Location::CELL> ta(1,"T","K",1),tb(1,"T","K",1),qa(1,"Q","",1),qb(1,"Q","",1);
+            ta.fill(300);tb.fill(400);qa.fill(10);qb.fill(-5);
+            EnergySolverControls ea,eb;
+            ea.density=2;eb.density=3;ea.cp=eb.cp=1;
+            ea.conductivity=1;eb.conductivity=4;ea.dt=eb.dt=dt;
+            ea.relaxation=eb.relaxation=1;ea.tolerance=eb.tolerance=1e-12;
+            ChtInterfaceControls c;
+            c.region1_patch="outlet";c.region2_patch="inlet";
+            c.conductivity1=1;c.conductivity2=4;c.relaxation=relaxation;
+            c.tolerance=1e-11;c.temperature_tolerance=1e-11;c.max_iterations=500;
+            ScalarBoundaryConditions ba,bb;
+            for(const auto* n : {"inlet","outlet","wall"}) {
+                ba[n]={ScalarBoundaryType::ZERO_GRADIENT,0,0};bb[n]=ba[n];
+            }
+            ba["outlet"]={ScalarBoundaryType::FIXED_VALUE,350,0};bb["inlet"]=ba["outlet"];
+            const double g=1.0/(0.5/1.0+0.5/4.0), m1=2/dt, m2=3/dt;
+            const double d=(m1+g)*(m2+g)-g*g, b1=m1*300+10, b2=m2*400-5;
+            const double exact1=((m2+g)*b1+g*b2)/d, exact2=(g*b1+(m1+g)*b2)/d;
+            auto check=[&](const char* name,bool ok,double value){++checks;failures+=!ok;
+                std::cout<<name<<"\tdt="<<dt<<"\trelax="<<relaxation<<'\t'<<(ok?"PASS":"FAIL")<<'\t'<<std::setprecision(17)<<value<<'\n';};
+            try {
+                auto r=solve_two_region_cht(a,ga,b,gb,fa,fb,ta,tb,qa,qb,ea,eb,c,ba,bb);
+                check("cht/status",r.converged,r.iterations);
+                check("cht/temperature1",std::isfinite(ta(0))&&std::abs(ta(0)-exact1)<1e-8,ta(0)-exact1);
+                check("cht/temperature2",std::isfinite(tb(0))&&std::abs(tb(0)-exact2)<1e-8,tb(0)-exact2);
+                const double balance=2*(ta(0)-300)+3*(tb(0)-400)-dt*5;
+                check("cht/energy",std::isfinite(balance)&&std::abs(balance)<1e-7,balance);
+            } catch(const std::exception& e) {
+                ++checks;++failures;std::cerr<<"cht/exception "<<e.what()<<'\n';
+            }
+        }
+    }
+    std::cout<<"SUMMARY\tchecks="<<checks<<"\tfailures="<<failures<<'\n';
+    return failures?1:0;
+}
--- /dev/null
+++ b/tests/audit/test_audit_energy_boundary_metric.cpp
@@ -0,0 +1,38 @@
+// H05: a one-cell transient balance with an independent normal-distance oracle.
+#include "chain_fixture.h"
+#include "cfdx/physics/energy_solver.h"
+#include <cmath>
+#include <iostream>
+#include <iomanip>
+using namespace cfdx::core;
+using namespace cfdx::physics;
+int main() {
+    int checks=0,failures=0;
+    for(const double shear:{0.,1.,3.}) for(const bool override_value:{false,true}) {
+        auto m=make_1d_chain(1);
+        for(std::size_t i=0;i<m.n_points();++i)
+            m.points().set(i,m.points().x(i),m.points().y(i)+shear*(m.points().x(i)-1),m.points().z(i));
+        const auto geometry=build_fv_geometry(m);
+        Field<double,Location::FACE> flux(m.n_faces(),"F","",1);flux.fill(0);
+        Field<double,Location::CELL> t(1,"T","K",1),old(1,"Told","K",1),source(1,"Q","",1);
+        old.fill(300);source.fill(0);t.fill((10.*300+4.*400)/14.);
+        EnergySolverControls c;c.cp=1;c.density=1;c.conductivity=2;c.dt=.1;
+        ScalarBoundaryConditions bc;
+        ScalarBoundaryFaceValues fv;
+        if(override_value) fv.values["outlet"].assign(m.n_faces(),400.);
+        else bc["outlet"]={ScalarBoundaryType::FIXED_VALUE,400,0};
+        const auto* values=override_value?&fv:nullptr;
+        const auto eq=assemble_energy_equation(m,geometry,flux,source,old,c,bc,values);
+        Vector x(1,t(0));
+        const double matrix=energy_balance_relative_from_equation(eq,x);
+        const double diagnostic=energy_balance_relative(m,geometry,flux,t,old,source,c,bc,values);
+        for(const auto& pair: {std::make_pair("matrix",matrix),std::make_pair("diagnostic",diagnostic)}) {
+            const bool ok=std::isfinite(pair.second)&&pair.second<1e-12;
+            ++checks;failures+=!ok;
+            std::cout<<(ok?"PASS":"FAIL")<<'\t'<<pair.first<<"\tshear="<<shear<<"\toverride="<<override_value
+                     <<"\tvalue="<<std::setprecision(17)<<pair.second<<'\n';
+        }
+    }
+    std::cout<<"SUMMARY checks="<<checks<<" failures="<<failures<<'\n';
+    return failures?1:0;
+}
--- /dev/null
+++ b/tests/audit/test_audit_flux_contract.cpp
@@ -0,0 +1,35 @@
+#include "chain_fixture.h"
+#include <iostream>
+#include <iomanip>
+#include <cmath>
+#include <limits>
+using namespace cfdx::core;
+using namespace cfdx::physics;
+int main(){
+ auto m=make_1d_chain(5);
+ for(std::size_t i=0;i<m.n_points();++i){double x=m.points().x(i);m.points().set(i,x*x,m.points().y(i),m.points().z(i));}
+ auto g=build_fv_geometry(m);
+ Field<double,Location::CELL> p(5,"p","",1),su(5,"Su","",1),sp(5,"Sp","",1);
+ for(std::size_t c=0;c<5;++c){p(c)=std::sin(11*g.cell_centres[c].x)+.1*c;su(c)=.2*(c+1);sp(c)=-.1;}
+ Field<double,Location::FACE> f(m.n_faces(),"F","",1);
+ for(std::size_t i=0;i<f.size();++i)f(i)=g.face_area_vectors[i].x;
+ std::vector<double> k{1,4,.2,2,3};
+ ScalarBoundaryConditions bc;bc["inlet"]={ScalarBoundaryType::FIXED_VALUE,1.3,0};bc["outlet"]={ScalarBoundaryType::FIXED_VALUE,-.4,0};bc["wall"]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
+ int total=0,failed=0;
+ for(bool neumann:{false,true})for(double direction:{-1.0,1.0}){
+ if(neumann){bc["inlet"]={ScalarBoundaryType::FIXED_GRADIENT,0,.7};bc["outlet"]={ScalarBoundaryType::FIXED_GRADIENT,0,-.4};}
+ else{bc["inlet"]={ScalarBoundaryType::FIXED_VALUE,1.3,0};bc["outlet"]={ScalarBoundaryType::FIXED_VALUE,-.4,0};}
+ for(std::size_t i=0;i<f.size();++i)f(i)=direction*g.face_area_vectors[i].x;
+ for(bool bounded:{false,true})for(auto scheme:{ConvectionScheme::UPWIND,ConvectionScheme::CENTRAL,ConvectionScheme::BLENDED,ConvectionScheme::QUICK,ConvectionScheme::QUICK_BOUNDED,ConvectionScheme::SECOND_ORDER_UPWIND,ConvectionScheme::TVD}){
+  auto eq=assemble_scalar_equation(m,g,f,.7,su,sp,bc,bounded,nullptr,nullptr,nullptr,&k,scheme,&p);
+  auto flux=reconstruct_scalar_transport_flux(m,g,f,p,.7,bc,bounded,scheme,nullptr,&k);
+  Vector x(5,0);for(size_t c=0;c<5;++c)x(c)=p(c);auto ax=eq.matrix.matvec(x);
+  std::vector<double> balance(5,0),divf(5,0);
+  for(size_t i=0;i<m.n_faces();++i){auto o=m.ownership().owner(i);auto n=m.ownership().neighbour(i);balance[o]+=flux(i);divf[o]+=f(i);if(n>=0){balance[n]-=flux(i);divf[n]-=f(i);}}
+  double e=0;for(size_t c=0;c<5;++c){double r=balance[c]-g.cell_volumes[c]*(su(c)+sp(c)*p(c))-(bounded?p(c)*divf[c]:0.0);double error=std::abs(r-(ax[c]-eq.rhs(c)));e=std::isfinite(error)?std::max(e,error):std::numeric_limits<double>::infinity();}
+  bool ok=std::isfinite(e)&&e<1e-11;++total;failed+=!ok;
+  std::cout<<"flux-contract\t"<<int(scheme)<<"\tneumann="<<neumann<<"\tdirection="<<direction<<"\tbounded="<<bounded<<"\t"<<(ok?"PASS":"FAIL")<<"\t"<<std::setprecision(17)<<e<<'\n';
+ }
+ }
+ std::cout<<"SUMMARY\tchecks="<<total<<"\tfailures="<<failed<<'\n';return failed?1:0;
+}
--- /dev/null
+++ b/tests/audit/test_audit_numerics.cpp
@@ -0,0 +1,94 @@
+// Independent regression cases for the 2026-10-24 audit.
+// Checks remain active with NDEBUG. Reference constants are analytic, not
+// extracted from the assembled matrix or reconstructed by the tested code.
+#include "affine_mesh_fixture.h"
+#include "chain_fixture.h"
+#include "cfdx/core/linalg/exact_schur.h"
+#include "cfdx/core/linalg/gmres_solver.h"
+#include "cfdx/core/linalg/bicgstab_solver.h"
+#include "cfdx/core/linalg/jacobi_solver.h"
+#include "cfdx/core/linalg/gauss_seidel_solver.h"
+#include "cfdx/core/linalg/matrix_free_fv_operator.h"
+#include "cfdx/core/numerics/newton.h"
+#include "cfdx/core/numerics/interpolation.h"
+#include "cfdx/physics/sst_solver.h"
+#include "cfdx/physics/thermal_zones.h"
+#include "cfdx/physics/thermophysical_models.h"
+#include "cfdx/physics/transport_models.h"
+#include "cfdx/physics/m1_m4_models.h"
+#include "cfdx/physics/low_storage_time_integration.h"
+#include "cfdx/physics/dual_time_field_integrator.h"
+#include "cfdx/thermodynamics/equation_of_state.h"
+#include <limits>
+using namespace cfdx::core;
+using namespace cfdx::physics;
+using Cell=Field<double,Location::CELL>;
+using Face=Field<double,Location::FACE>;
+static int checks=0, failures=0;
+static double max_error(double a,double b){
+    // Fail closed: std::max alone can hide a NaN comparison.
+    if(!std::isfinite(a)||!std::isfinite(b))return std::numeric_limits<double>::infinity();
+    return std::max(a,b);
+}
+static void check(const std::string& name,bool ok,double value=0){++checks;failures+=!ok;std::cout<<name<<'\t'<<(ok?"PASS":"FAIL")<<'\t'<<std::setprecision(17)<<value<<'\n';}
+static void near(const std::string& name,double a,double b,double tol=1e-11){check(name,std::isfinite(a)&&std::abs(a-b)<=tol*std::max(1.0,std::abs(b)),a);}
+template<class F>void rejects(const std::string& name,F f){bool yes=false;try{f();}catch(const std::exception&){yes=true;}check(name,yes);}
+static SparseMatrix scalar(double v){SparseMatrix a(1,1);a.push_back(0,0,v);a.finalize();return a;}
+int main(){
+ try {
+ for(double s:{1.0,1e200,1e-200}){Vector v(1,s);near("norm/vector",v.norm2()/s,1);near("norm/krylov",krylov_norm2(v)/s,1);}
+ {Vector v(1,std::numeric_limits<double>::quiet_NaN());check("norm/NaN",std::isnan(v.norm_inf()));}
+ for(int j:{0,1}){Vector b(1,1e200),x(1,0);auto r=j?solve_gauss_seidel(scalar(1),b,x,10,1e-10):solve_jacobi(scalar(1),b,x,10,1e-10);check("stationary/huge/status",r.status==SolverStatus::CONVERGED);near("stationary/huge/solution",x(0)/b(0),1);}
+ for(double s:{1.,1e-200}){Cell u(1,"u","",1);u(0)=0;NewtonFieldControls ctl;ctl.max_iterations=1;auto r=solve_newton_field(u,[s](const Cell& q,Cell& res,SparseMatrix& J){res(0)=s*(q(0)-1);J=scalar(s);},ctl);check("newton/last-iteration/status",r.status==SolverStatus::CONVERGED);near("newton/solution",u(0),1);}
+ for(double s:{1.,1e12,1e24,1e-100,1e-200}){Vector b(1,s),x(1,0);auto r=solve_cg(scalar(s),b,x,20,1e-10);check("CG/scale/status",r.status==SolverStatus::CONVERGED);near("CG/scale/solution",x(0),1);}
+ for(int which:{0,1,2,3}){PrecisionPolicy p;p.enabled=true;Vector b(1,1),x(1,0);auto a=scalar(1.00000001);SolverResult r;
+  if(which==0)r=solve_cg(a,b,x,30,1e-12,p);
+  if(which==1)r=solve_gmres(a,b,x,1,30,1e-12,nullptr,p);
+  if(which==2)r=solve_fgmres(a,b,x,1,30,1e-12,nullptr,p);
+  if(which==3)r=solve_bicgstab(a,b,x,30,1e-12,nullptr,p);
+  check("mixed/status",r.status==SolverStatus::CONVERGED);near("mixed/FP64-residual",1-a.matvec(x)[0],0,1e-12);
+ }
+ {Vector b(1,0),x(1,1);auto r=solve_gmres(scalar(1),b,x,1,5,1e-12);check("GMRES/zero-rhs",r.status==SolverStatus::CONVERGED);near("GMRES/zero-rhs/x",x(0),0);}
+ for(double s:{1.,1e-16}){Vector b(1,s),x(1,0);auto r=solve_bicgstab(scalar(s),b,x,20,1e-12);check("BiCG/scale",r.status==SolverStatus::CONVERGED);near("BiCG/scale/x",x(0),1);}
+ {auto a=scalar(1),g=scalar(1),d=scalar(-1),c=scalar(0);BlockOperator op(a,g,d,c);ExactSchurApproximation sh([](const Vector& b,Vector& x){x=b;return true;});check("Schur/setup",sh.setup(op));Vector b(1,std::numeric_limits<double>::quiet_NaN()),x(1,0);check("Schur/reject-NaN",!sh.apply(b,x));}
+ {auto grid=make_affine_cube(3,0,1);auto& m=grid.mesh;Cell p(m.n_cells(),"p","",1);for(size_t c=0;c<m.n_cells();++c)p(c)=grid.geometry.cell_centres[c].x;
+  auto g=compute_gradient_weighted_least_squares(p,m,GradientWeighting::INVERSE_DISTANCE_SQUARED,std::numeric_limits<double>::infinity(),BoundaryGradientPolicy::EXCLUDE_BOUNDARY);double e=0;for(size_t c=0;c<m.n_cells();++c)e=max_error(e,std::abs(g(c,0)-1));near("WLS/exclude/all-cells",e,0);
+  auto ghost=compute_gradient_least_squares(p,m,BoundaryGradientPolicy::ZERO_GRADIENT_GHOST);near("LS/ghost/corner",ghost(0,0),0.5);
+  FvDiffusionOperator zero(m,grid.geometry,0);Vector diag;zero.diagonal(diag);near("matrixfree/zero-diagonal",diag(0),0);
+  rejects("matrixfree/reject-NaN",[&]{FvDiffusionOperator bad(m,grid.geometry,std::numeric_limits<double>::quiet_NaN());});
+ }
+ for(size_t n:{2u,3u}){auto grid=make_affine_cube(n,0,1);Cell p(grid.mesh.n_cells(),"p","",1);for(size_t c=0;c<p.size();++c)p(c)=grid.geometry.cell_centres[c].x;auto q=compute_quadratic_fit_least_squares(p,grid.mesh);double e=0;for(size_t c=0;c<p.size();++c)e=max_error(e,std::abs(q(c,0)-1));near("QUICK/fallback-affine",e,0);}
+ for(double l:{1.,1e-6,1e-8}){auto grid=make_affine_cube(4,0,1);auto& m=grid.mesh;for(size_t i=0;i<m.n_points();++i)m.points().set(i,l*m.points().x(i),l*m.points().y(i),l*m.points().z(i));auto geo=make_geometry_cache(m);Cell p(m.n_cells(),"p","",1);for(size_t c=0;c<p.size();++c)p(c)=std::pow(geo.cell_centres[c].x/l,2);auto q=compute_quadratic_fit_least_squares(p,m);double e=0;for(size_t c=0;c<p.size();++c)e=max_error(e,std::abs(l*l*q(c,3)-1));near("QUICK/scale/curvature",e,0);}
+ {TurbulenceTransportControls c;auto s=sst_cell_sources(1,1,100,1,1,0,c);near("SST/high-strain",s.sw,161.29032258064516);auto lo=sst_cell_sources(1,1,0.01,1,1,0,c);near("SST/low-strain",lo.sw,(5.0/9.0)*1e-4);
+  auto w=komega2006_cell_sources(1,1,100,0,c);near("Wilcox/limited-viscosity",w.nut,0.0034285714285714286);near("Wilcox/unlimited-diffusivity",w.gamma_k,c.density*(c.molecular_viscosity+c.komega_sigma_k));
+ }
+ {ThermophysicalProperties p;auto& c=p.heat_capacity;c.model=ScalarPropertyModel::TABLE;c.table.temperature={300,400};c.table.value={1000,1200};
+  near("enthalpy/clamp/below",p.enthalpy(200,100),100000);near("enthalpy/clamp/above",p.enthalpy(600,500),120000);near("enthalpy/zero",p.enthalpy(200,200),0);
+  c.table.extrapolation=ExtrapolationPolicy::LINEAR;near("enthalpy/linear/below",p.enthalpy(200,100),70000);near("enthalpy/linear/above",p.enthalpy(600,500),150000);near("enthalpy/antisymmetry",p.enthalpy(100,600)+p.enthalpy(600,100),0,1e-8);
+  near("enthalpy/additivity",p.enthalpy(600,100)-p.enthalpy(350,100)-p.enthalpy(600,350),0,1e-8);
+  c.enforce_bounds=true;c.minimum=700;c.maximum=1400;c.model=ScalarPropertyModel::LINEAR;c.reference_temperature=300;c.reference_value=1000;c.coefficients={5};near("enthalpy/clipped-linear",p.enthalpy(500,200),343000);
+  c.model=ScalarPropertyModel::POLYNOMIAL;c.coefficients={1000};c.maximum=800;near("enthalpy/clipped-polynomial",p.enthalpy(400,300),80000);
+  c.model=ScalarPropertyModel::CONSTANT;c.reference_value=1000;c.maximum=800;near("enthalpy/clipped-constant",p.enthalpy(400,300),80000);
+  rejects("enthalpy/NaN-interval-count",[&]{p.enthalpy(400,300,std::numeric_limits<double>::quiet_NaN());});
+ }
+ {auto t=compute_transport(300,1.2);near("transport/air-Cp",t.Cp,1.4*AIR_R/0.4);cfdx::thermodynamics::IdealGasEOS eos;near("EOS/dp-dT",eos.dp_dT_rho(101325,300),337.75);near("P1/full-solid-angle",m1m4::p1_source(1,0,1000),226814.97676,1e-10);}
+ for(int method:{0,1,2}){
+  auto step=[method](Cell& u,double dt,auto&& rhs){if(method==0)low_storage_rk2_step(u,dt,rhs);if(method==1)low_storage_rk3_step(u,dt,rhs);if(method==2)ssprk3_step(u,dt,rhs);};
+  Cell u(1,"u","",3);for(size_t d=0;d<3;++d)u(0,d)=d+1;
+  step(u,.1,[](const Cell& x,Cell& r){r=Cell(x.size(),"r","",x.dimension());r.fill(1);});
+  double e=0;for(size_t d=0;d<3;++d)e=max_error(e,std::abs(u(0,d)-(d+1.1)));near("RK/all-components",e,0);
+  double old=0;for(int n:{10,20,40}){Cell q(1,"q","",3);q.fill(1);for(int i=0;i<n;++i)step(q,1.0/n,[](const Cell& x,Cell& r){r=x;});double err=std::abs(q(0,2)-std::exp(1.0));if(old>0){double order=std::log(old/err)/std::log(2.0);check("RK/observed-order",order>(method==0?1.9:2.9),order);}old=err;}
+ }
+ {Cell u(1,"u","",1);u(0)=1;rejects("implicit/reject-nonfinite",[&]{advance_time(u,1,[](const Cell& q,Cell& r){r(0)=-2*std::sqrt(q(0));},TimeScheme::EULER_IMPLICIT);});}
+ {Cell u(1,"u","",1);u.fill(0);DualTimeFieldIntegrator it(u);DualTimeFieldSolveControls c;c.use_embedded_be_estimator=false;c.pseudo_time.pseudo_dt_initial=.5;c.pseudo_time.pseudo_dt_max=.5;c.pseudo_time.absolute_tolerance=1e-13;c.pseudo_time.relative_tolerance=1e-13;
+  double t=0;for(double dt:{.1,.2,.05,.08}){auto r=it.advance(dt,c,[](const Cell&,Cell& rhs){rhs(0)=1;});t+=dt;near("dual-time/variable-step",r.state(0),t,1e-11);}
+ }
+ {auto m=make_1d_chain(2);for(size_t i=0;i<m.n_points();++i){double x=m.points().x(i);x=x<=.5?.5*x:.25+1.5*(x-.5);m.points().set(i,x,m.points().y(i),m.points().z(i));}auto g=build_fv_geometry(m);Cell p(2,"p","",1),z(2,"z","",1);z.fill(0);for(size_t c=0;c<2;++c)p(c)=g.cell_centres[c].x;Face f(m.n_faces(),"f","",1);f.fill(0);auto q=interpolate_cell_to_face(p,m,InterpScheme::LINEAR);size_t inner=0;for(size_t i=0;i<m.n_faces();++i)if(m.ownership().neighbour(i)>=0)inner=i;near("linear/nonuniform",q(inner),.25);
+  std::vector<double> k{1,4};auto a=assemble_scalar_equation(m,g,f,1,z,z,{},false,nullptr,nullptr,nullptr,&k);near("diffusion/material-series",-a.matrix(0,1),32.0/7.0);
+  auto model=build_thermal_zones(m,{0,1},{{"a",0,1,0,.3,true},{"b",1,4,0,.3,true}});ScalarBoundaryConditions bc;bc["inlet"]={ScalarBoundaryType::FIXED_VALUE,300,0};auto at=assemble_thermal_conduction(m,g,model,bc);near("thermal/contact-series",-at.matrix(0,1),1/(.125+.375/4+.3));
+  f(inner)=1;auto b=assemble_scalar_equation(m,g,f,1,z,z,{},false);near("transport/local-imbalance",b.max_imbalance,1);
+  k[0]=std::numeric_limits<double>::quiet_NaN();rejects("diffusion/reject-NaN-material",[&]{assemble_scalar_equation(m,g,f,1,z,z,{},false,nullptr,nullptr,nullptr,&k);});
+ }
+ } catch(const std::exception& e) { ++failures;std::cerr<<"UNEXPECTED_EXCEPTION "<<e.what()<<'\n'; }
+ std::cout<<"SUMMARY\tchecks="<<checks<<"\tfailures="<<failures<<'\n';return failures?1:0;
+}
--- a/tests/unit/test_interpolation.cpp
+++ b/tests/unit/test_interpolation.cpp
@@ -139,20 +139,19 @@
     });
 
     run_case("interpolate_linear_scalar", []() {
-        Mesh m = make_two_cell_mesh();
-        ScalarCellField f(2, "p", "Pa", 1);
-        f(0) = 10.0;
-        f(1) = 30.0;
-
-        auto face_field = interpolate_cell_to_face(f, m, InterpScheme::LINEAR);
-        // face 0 : owner=0, boundary → 10.0
-        EXPECT_TRUE(face_field(0) == 10.0);
-        // face 1 : 0.5*(10+30) = 20.0
-        EXPECT_TRUE(face_field(1) == 20.0);
-        // face 2 : 0.5*(30+10) = 20.0
-        EXPECT_TRUE(face_field(2) == 20.0);
-        // face 3 : owner=1, boundary → 30.0
-        EXPECT_TRUE(face_field(3) == 30.0);
+        // The old fixture consisted of open coplanar triangles. Its face was
+        // not midway between its computed centres, so 1/2 was not an oracle.
+        Mesh m = make_two_cell_unit_cubes();
+        ScalarCellField f(2, "p", "Pa", 1);
+        f(0)=10.0; f(1)=30.0;
+        auto face_field=interpolate_cell_to_face(f,m,InterpScheme::LINEAR);
+        EXPECT_TRUE(face_field(0)==10.0);
+        EXPECT_TRUE(face_field(5)==20.0);
+        EXPECT_TRUE(face_field(6)==30.0);
+        // Reverse owner/neighbour ordering; physical interpolation is invariant.
+        m.ownership().set_owner(5,1);m.ownership().set_neighbour(5,0);
+        face_field=interpolate_cell_to_face(f,m,InterpScheme::LINEAR);
+        EXPECT_TRUE(face_field(5)==20.0);
     });
 
     run_case("tvd_mc_limiter_contract", []() {
@@ -263,7 +262,7 @@
     });
 
     run_case("interpolate_vector_dim3", []() {
-        Mesh m = make_two_cell_mesh();
+        Mesh m = make_two_cell_unit_cubes();
         Field<double, Location::CELL> f(2, "U", "m/s", 3);
         f.set(0, 1.0, 2.0, 3.0);
         f.set(1, 7.0, 8.0, 9.0);
@@ -271,7 +270,7 @@
         auto face_field = interpolate_cell_to_face(f, m, InterpScheme::LINEAR);
         EXPECT_TRUE(face_field.dimension() == 3);
         double x, y, z;
-        face_field.get(1, x, y, z);
+        face_field.get(5, x, y, z);
         EXPECT_TRUE(x == 4.0);
         EXPECT_TRUE(y == 5.0);
         EXPECT_TRUE(z == 6.0);
--- a/tests/unit/test_transport_conservation.cpp
+++ b/tests/unit/test_transport_conservation.cpp
@@ -85,10 +85,11 @@
         const auto reconstructed =
             reconstruct_scalar_transport_flux(
                 m,g,phi,q,0.0,bc,true,ConvectionScheme::UPWIND);
-        EXPECT_NEAR(reconstructed(0),3.0,1e-15);
+        // Outflow uses owner q=2, so 3*2=6. Dirichlet enters diffusion only.
+        EXPECT_NEAR(reconstructed(0),6.0,1e-15);
         const auto audit=audit_face_flux_conservation(m,reconstructed);
-        EXPECT_NEAR(audit.global_boundary_flux,3.0,1e-15);
-        EXPECT_NEAR(audit.global_cell_balance,3.0,1e-15);
+        EXPECT_NEAR(audit.global_boundary_flux,6.0,1e-15);
+        EXPECT_NEAR(audit.global_cell_balance,6.0,1e-15);
     });
 
 
--- a/tests/unit/test_turbulence_qualification_equations.cpp
+++ b/tests/unit/test_turbulence_qualification_equations.cpp
@@ -55,8 +55,10 @@
     {
         const double grad_cross=-0.7;
         const auto q=komega2006_cell_sources(k,omega,S,grad_cross,c);
-        const double nut=k/omega, P=c.density*nut*S*S;
-        EXPECT_CLOSE(q.nut,nut); EXPECT_CLOSE(q.sk,P);        EXPECT_CLOSE(q.spk,-c.density*c.beta_star*omega);        EXPECT_CLOSE(q.sw,c.komega_alpha*(omega/k)*P);        EXPECT_CLOSE(q.spw,-c.density*c.komega_beta0*omega);        EXPECT_CLOSE(q.gamma_k,c.density*(c.molecular_viscosity+c.komega_sigma_k*nut));        EXPECT_CLOSE(q.gamma_w,c.density*(c.molecular_viscosity+c.komega_sigma_w*nut));
+        // Wilcox2006 TMR: nu_t is strain-limited, diffusivity uses k/omega.
+        const double nut=k/std::max(omega,c.komega_clim*S/std::sqrt(c.beta_star));
+        const double P=c.density*nut*S*S;
+        EXPECT_CLOSE(q.nut,nut); EXPECT_CLOSE(q.sk,P);        EXPECT_CLOSE(q.spk,-c.density*c.beta_star*omega);        EXPECT_CLOSE(q.sw,c.komega_alpha*(omega/k)*P);        EXPECT_CLOSE(q.spw,-c.density*c.komega_beta0*omega);        EXPECT_CLOSE(q.gamma_k,c.density*(c.molecular_viscosity+c.komega_sigma_k*k/omega));        EXPECT_CLOSE(q.gamma_w,c.density*(c.molecular_viscosity+c.komega_sigma_w*k/omega));
         const auto qp=komega2006_cell_sources(k,omega,S,0.7,c);
         EXPECT_CLOSE(qp.sw,c.komega_alpha*(omega/k)*P+c.density*c.komega_sigma_d0/omega*0.7);
     }
@@ -77,7 +79,7 @@
         const double nut=c.a1*k/std::max(c.a1*omega,S*F2);
         const double P=std::min(c.density*nut*S*S,c.sst_production_limiter*c.beta_star*c.density*k*omega);
         const double cross_term=2.0*(1.0-F1)*c.density*c.sst_sigma_w2/omega*cross;
-        EXPECT_CLOSE(q.nut,nut); EXPECT_CLOSE(q.sk,P);        EXPECT_CLOSE(q.spk,-c.density*c.beta_star*omega);        EXPECT_CLOSE(q.sw,blend(c.gamma1,c.gamma2)*c.density*S*S+std::max(cross_term,0.0));        EXPECT_CLOSE(q.spw,-c.density*blend(c.beta1,c.beta2)*omega+std::min(cross_term,0.0)/omega);        EXPECT_CLOSE(q.gamma_k,c.density*(c.molecular_viscosity+blend(c.sst_sigma_k1,c.sst_sigma_k2)*nut));        EXPECT_CLOSE(q.gamma_w,c.density*(c.molecular_viscosity+blend(c.sst_sigma_w1,c.sst_sigma_w2)*nut));
+        EXPECT_CLOSE(q.nut,nut); EXPECT_CLOSE(q.sk,P);        EXPECT_CLOSE(q.spk,-c.density*c.beta_star*omega);        EXPECT_CLOSE(q.sw,blend(c.gamma1,c.gamma2)*P/nut+std::max(cross_term,0.0));        EXPECT_CLOSE(q.spw,-c.density*blend(c.beta1,c.beta2)*omega+std::min(cross_term,0.0)/omega);        EXPECT_CLOSE(q.gamma_k,c.density*(c.molecular_viscosity+blend(c.sst_sigma_k1,c.sst_sigma_k2)*nut));        EXPECT_CLOSE(q.gamma_w,c.density*(c.molecular_viscosity+blend(c.sst_sigma_w1,c.sst_sigma_w2)*nut));
     }
 
     c.model=TurbulenceModel::SPALART_ALLMARAS;
--- a/tests/validation/test_cht_validation.cpp
+++ b/tests/validation/test_cht_validation.cpp
@@ -6,13 +6,13 @@
 using namespace cfdx::physics;
 using namespace cfdx::testing;
 
-static Mesh cube(const char* patch_name)
+static Mesh cube(const char* patch_name, bool adjacent_interface=false, bool right_region=false)
 {
     Mesh m;
     m.points().resize(8);
     const double p[8][3]={{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                           {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
-    for(std::size_t i=0;i<8;++i)m.points().set(i,p[i][0],p[i][1],p[i][2]);
+    for(std::size_t i=0;i<8;++i)m.points().set(i,p[i][0]+(right_region?1.0:0.0),p[i][1],p[i][2]);
     m.faces().push_face({0,3,2,1});m.faces().push_face({4,5,6,7});
     m.faces().push_face({0,1,5,4});m.faces().push_face({3,7,6,2});
     m.faces().push_face({0,4,7,3});m.faces().push_face({1,2,6,5});
@@ -20,14 +20,22 @@
     for(std::size_t f=0;f<6;++f){m.ownership().set_owner(f,0);m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);}
     m.cells().push_cell({0,1,2,3,4,5});
     Patch patch;patch.name=patch_name;patch.type=PatchType::WALL;patch.face_ids={0,1,2,3,4,5};
+    if(adjacent_interface) patch.face_ids={right_region?4u:5u};
     m.boundary().add_patch(patch);
+    if(adjacent_interface) {
+        Patch walls;walls.name="walls";walls.type=PatchType::WALL;
+        for(std::size_t f=0;f<6;++f) if(f!=patch.face_ids[0]) walls.face_ids.push_back(f);
+        m.boundary().add_patch(walls);
+    }
     return m;
 }
 
 int main()
 {
     run_case("two_region_cht_independent_energy_conservation",[] {
-        Mesh m1=cube("interface1"),m2=cube("interface2");
+        // Adjacent regions, one common face and opposite outward normals.
+        // The old fixture overlapped both cubes and coupled all six faces.
+        Mesh m1=cube("interface1",true,false),m2=cube("interface2",true,true);
         auto g1=build_fv_geometry(m1),g2=build_fv_geometry(m2);
         Field<double,Location::FACE> phi1(m1.n_faces(),"phi1","kg/s",1);
         Field<double,Location::FACE> phi2(m2.n_faces(),"phi2","kg/s",1);
@@ -49,8 +57,10 @@
 
         // Reconstruct the interface state from the accepted temperatures and
         // independently rebuild each region's physical energy flux.
-        const double d1=(g1.face_centres[0]-g1.cell_centres[0]).mag();
-        const double d2=(g2.face_centres[0]-g2.cell_centres[0]).mag();
+        const auto face1=m1.boundary().patch(0).face_ids[0];
+        const auto face2=m2.boundary().patch(0).face_ids[0];
+        const double d1=boundary_normal_distance(g1,face1,0);
+        const double d2=boundary_normal_distance(g2,face2,0);
         const double h1=c.conductivity1/d1;
         const double h2=c.conductivity2/d2;
         const double tint=(h1*T1(0)+h2*T2(0))/(h1+h2);
@@ -75,7 +85,7 @@
 
         // Independent interface heat-flow equality: each side must carry the
         // same physical heat rate, with opposite outward signs.
-        const double interface_area=g1.face_area_vectors[0].mag()*static_cast<double>(m1.boundary().patch(0).face_ids.size());
+        const double interface_area=g1.face_area_vectors[face1].mag();
         const double qflux1=c.conductivity1*(T1(0)-tint)/d1;
         const double qflux2=c.conductivity2*(tint-T2(0))/d2;
         EXPECT_TRUE(std::isfinite(qflux1) && std::isfinite(qflux2));
--- a/tests/validation/test_level_c_coupled_verification.cpp
+++ b/tests/validation/test_level_c_coupled_verification.cpp
@@ -168,6 +168,12 @@
         {
             Mesh m1=cube_with_patches();
             Mesh m2=cube_with_patches();
+            // Cold region lies above the hot one, not on top of its volume.
+            // Its lower y face is the common interface.
+            for(std::size_t i=0;i<m2.n_points();++i)
+                m2.points().set(i,m2.points().x(i),m2.points().y(i)+1.0,m2.points().z(i));
+            m2.boundary().patch(2).name="interface";
+            m2.boundary().patch(3).name="y_max";
             auto g1=build_fv_geometry(m1);
             auto g2=build_fv_geometry(m2);
 
@@ -195,6 +201,8 @@
 
             auto bc1=fixed_wall_conditions(500.0,"interface");
             auto bc2=fixed_wall_conditions(300.0,"interface");
+            bc2.erase("y_min");
+            bc2["y_max"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
 
             const auto r=solve_two_region_cht(
                 m1,g1,m2,g2,f1,f2,T1,T2,s1,s2,ec,ec,cc,bc1,bc2);
--- a/tests/validation/test_polyhedral_gradient_campaign.cpp
+++ b/tests/validation/test_polyhedral_gradient_campaign.cpp
@@ -647,6 +647,9 @@
 
     for (const auto& mesh_case : cases) {
         for (const auto weighting : weightings) {
+            // Retain the refinement history instead of recreating a size-one
+            // vector on every mesh (which made the decay check vacuous).
+            std::vector<double> smooth_errors;
             for (const std::size_t n : {4u, 8u, 16u}) {
                 const Grid grid = make_tet_grid(
                     n, mesh_case.sx, mesh_case.sy, mesh_case.sz,
@@ -656,7 +659,6 @@
                 double min_condition = std::numeric_limits<double>::infinity();
                 std::size_t checked = 0;
                 const auto smooth = sample_field(grid, kSmooth);
-                std::vector<double> smooth_errors;
 
                 const auto* cell_faces = grid.mesh.cells().faces_data();
                 const auto* cell_offsets = grid.mesh.cells().offsets_data();
@@ -707,8 +709,27 @@
                 require(checked == grid.interior.size(),
                         "mesh WLS campaign: not all interior stencils were checked");
 
+                // Only interior errors are measured below. The field API
+                // nevertheless evaluates boundary cells as well, so request
+                // the ghost policy explicitly rather than relying on the old
+                // bug that ignored EXCLUDE_BOUNDARY.
                 const auto grad = compute_gradient_weighted_least_squares(
-                    smooth, grid.mesh, weighting);
+                    smooth, grid.mesh, weighting,
+                    std::numeric_limits<double>::infinity(),
+                    BoundaryGradientPolicy::ZERO_GRADIENT_GHOST);
+                if (n == 4u) {
+                    bool rejected = false;
+                    try {
+                        (void)compute_gradient_weighted_least_squares(
+                            smooth, grid.mesh, weighting,
+                            std::numeric_limits<double>::infinity(),
+                            BoundaryGradientPolicy::EXCLUDE_BOUNDARY);
+                    } catch (const std::runtime_error&) {
+                        rejected = true;
+                    }
+                    require(rejected,
+                        "boundary tetrahedra without ghost samples must expose rank loss");
+                }
                 const auto smooth_error = gradient_error(
                     grid, grad, kSmooth, true);
                 require(std::isfinite(smooth_error.l2),
--- a/tests/validation/test_thermal_challenge_vv.cpp
+++ b/tests/validation/test_thermal_challenge_vv.cpp
@@ -327,6 +327,7 @@
         Field<double,Location::CELL> s1(1,"s1","W/m3",1), s2(1,"s2","W/m3",1); s1.fill(0); s2.fill(0);
         EnergySolverControls e1; e1.conductivity=2; e1.relaxation=1; e1.max_iterations=20; e1.tolerance=1e-10;
         EnergySolverControls e2=e1;
+        e2.conductivity=1; // Same cold-side material as the analytical oracle.
         ChtInterfaceControls c; c.region1_patch="interface"; c.region2_patch="interface";
         c.conductivity1=2; c.conductivity2=1; c.tolerance=1e-10; c.temperature_tolerance=1e-10;
         ScalarBoundaryConditions b1{{"outer_hot",{ScalarBoundaryType::FIXED_VALUE,400,0}},
Annexe E. Manifeste machine-readable de la revision R2
{
  "revision": "R2",
  "date": "2026-10-24",
  "base_commit_from_archive_comment": "3f7a764799ba4db5fbb5e1495b85e8fc5c31c6bb",
  "archive_name": "cfdx_dev-master.zip",
  "archive_sha256": "11f92e9f842d16b184432051208c40262a0961085d86ac933946b54613b955bb",
  "patch_sha256": "8df8479192656cf47155ace39666ea245b4dcfb6d2bccb03fecb69697a12be59",
  "source_tree_sha256": "73f6e16afed4e13065abf3c2e60225f3eca2c588b739e6c22b76b9d8d8435628",
  "source_tree_files": 1181,
  "changed_or_added_files": 43,
  "ctest_targets": [
    "test_audit_numerics",
    "test_audit_flux_contract",
    "test_audit_cht_time",
    "test_audit_cht_contract",
    "test_audit_energy_boundary_metric",
    "test_level_c_coupled_verification",
    "test_cht_validation",
    "test_thermal_vv",
    "test_thermal_zones_vv",
    "test_thermal_challenge_vv"
  ],
  "checks": {
    "numerics": 87,
    "flux_contract": 56,
    "cht_time": 24,
    "cht_contract": 101,
    "energy_boundary_metric": 12,
    "total": 280
  },
  "scope": "CPU serial; selected tests only; no live PR verification; inline new paths sanitizer-instrumented, core library Release",
  "files": [
    {
      "path": "CMakeLists.txt",
      "new": false,
      "sha256_before": "a1b238fa789120dbf3ba950901d05a2c344d9156c8f56207fa911055e47503a1",
      "sha256_after": "3faed951a2f293b4b018db7ff27b18e0f2f8e60984572f13deb9af570a54ee2e"
    },
    {
      "path": "src/cfdx/core/linalg/bicgstab_solver.h",
      "new": false,
      "sha256_before": "0346943c2f52436f4f1ca2328ce6ccfde68371068c09de4b5cab175e4b6b49ce",
      "sha256_after": "e2c73692b3e17cee5d4d7c6e041af937472cdfce73dd1de4a1a3f717e9137b92"
    },
    {
      "path": "src/cfdx/core/linalg/cg_solver.h",
      "new": false,
      "sha256_before": "4a4992c67c289b44da27baaceb866885602e240dfe6dbe8c2e5b41e8dc2d9d23",
      "sha256_after": "f8f10b1a217825a33e2feaf5b3ac9777e075339e322f729d7369ce72ec7adda3"
    },
    {
      "path": "src/cfdx/core/linalg/communication_avoiding.h",
      "new": false,
      "sha256_before": "a0d50c345f0eb12999057f6ce800511b4e8664e710c8bbc3831ee6eeef945fef",
      "sha256_after": "6d674e4283cd6766e3df737a2d87b4693657bf96e917b24561257904d45f0db1"
    },
    {
      "path": "src/cfdx/core/linalg/exact_schur.h",
      "new": false,
      "sha256_before": "4c81f9ed8738a6d2e812b4b9571ccaeed14d02c71c079a672e3ba719c9ed4878",
      "sha256_after": "12e3070a773dca3be038393c4cee8419714fda5deeab1e608f225e215fe1d24b"
    },
    {
      "path": "src/cfdx/core/linalg/gauss_seidel_solver.h",
      "new": false,
      "sha256_before": "dc300ae5b6d8a15e5370e12665e3e1dabff7157df376c439eef6291cb5eeda82",
      "sha256_after": "5f91a77eb1192a267171a291f56243f8662af343863e401ac7315756f9906eb5"
    },
    {
      "path": "src/cfdx/core/linalg/gmres_solver.h",
      "new": false,
      "sha256_before": "e5696421ec7f1ee11c05d879efbbf4497ef92434589f2104f4910532034f9fb4",
      "sha256_after": "58566d6a821214589d19cdeec554c6aa7e8c5c0dce19740cd4cbc719e9a00714"
    },
    {
      "path": "src/cfdx/core/linalg/jacobi_solver.h",
      "new": false,
      "sha256_before": "316bf787d1a37994a82f4a666119cf0cdf01b74dce9b837f69076aa5c2d53a62",
      "sha256_after": "f1c6f8ff8755d73d851d3e180542d41f66faa6b8a065ccc9bfaca775682d5e99"
    },
    {
      "path": "src/cfdx/core/linalg/krylov_reductions.h",
      "new": false,
      "sha256_before": "02ded50087d9903924b684c2c2ede949d660ef75057fe18b59756f8749512346",
      "sha256_after": "8210d38a7e1b980a7580cb52d619b76cae2613101bf47487369af497616637bb"
    },
    {
      "path": "src/cfdx/core/linalg/matrix_free_fv_operator.h",
      "new": false,
      "sha256_before": "8c7427d26f7c57494f4e40edccc6f448c5ebfae8cef2f998acb4522bd438d322",
      "sha256_after": "70a14926a122d4629b272ba3019d316fab7c75b42312d0fea2c92e589ec99286"
    },
    {
      "path": "src/cfdx/core/linalg/vector.h",
      "new": false,
      "sha256_before": "8f566fccbbedeed7a226a578054594237fa7f68e96af9154aaa155baf3da31c3",
      "sha256_after": "1605f067f717c82e1bf986e1c51a727796277dc200ece17bc9871869b33d8467"
    },
    {
      "path": "src/cfdx/core/numerics/gradient.h",
      "new": false,
      "sha256_before": "3beb28b6a7332bc5fb257e775bd3bdbe1ee60bef90b0df9dbc65ec6c15e20d72",
      "sha256_after": "f44d5f199036af258700036bb74d0bbd374335ef2b0fd5d8339c8b7d7a82350c"
    },
    {
      "path": "src/cfdx/core/numerics/interpolation.h",
      "new": false,
      "sha256_before": "2140b175c35b862f50892ce6d28f76b0e39ef0324d74f0cca1940170474007ff",
      "sha256_after": "30acbfed14778a67c83904c3096fe85b8065fe0721781db7b423ae08a74ce22e"
    },
    {
      "path": "src/cfdx/core/numerics/newton.h",
      "new": false,
      "sha256_before": "7a08072c8c6348fb432f32ad87a2cc8531767537074afafcbe8329c5c0c95fa1",
      "sha256_after": "103f0e3651ce5242b3f2a543957a615c8f17fa36bbdd236b171c6d2c0cbc9aee"
    },
    {
      "path": "src/cfdx/core/numerics/temporal.h",
      "new": false,
      "sha256_before": "b396b77d97ce8e37dd493bc73f29329c83330a793593b5e274edae59faaff27f",
      "sha256_after": "6fbbb7c0fa21a216a474d6562b99ae665dac8350e53177d342ce02a85b98e798"
    },
    {
      "path": "src/cfdx/physics/cht_solver.h",
      "new": false,
      "sha256_before": "228439fb93403c4f60c68b2c08a2bff3af5bb51fe678b95e6d59dbd361e9ab22",
      "sha256_after": "e95b63a64b05a219039fcfa2a90e171bfdf33720e930e726b010af0be92a9db3"
    },
    {
      "path": "src/cfdx/physics/dual_time_driver.h",
      "new": false,
      "sha256_before": "a464f5ae6923adb0667c07af85f5f5584d9707cf4ec3c507bf95c7078bea8320",
      "sha256_after": "0477962a27c34344161123386717d1226652fe211072a69a222baa448922e4ce"
    },
    {
      "path": "src/cfdx/physics/dual_time_field_integrator.h",
      "new": false,
      "sha256_before": "68e444e8b49a9cfdb3c909816d9f671ac79a9b3d5732815945c110c76568153a",
      "sha256_after": "b737bd72ab967ec02afc62aac829bb59d0e5c0f5a77dd70384daefe7a6ca2711"
    },
    {
      "path": "src/cfdx/physics/dual_time_stepping.h",
      "new": false,
      "sha256_before": "5b90e87697c7e59e0c0357b843d884405e0ba324e16e1409073999770eb1a3f7",
      "sha256_after": "ae2f27204633041a3049707193126206b5f3dbeb3b9039ac47b3c064535fca0f"
    },
    {
      "path": "src/cfdx/physics/energy_solver.h",
      "new": false,
      "sha256_before": "95c9564b2c8d3383f8f39b70163da55e93bcb6b517e5bf360a614c53a731a25d",
      "sha256_after": "aebbd340da73cd090f0559d2cd98393dc8ca0594a6512d913c4d540563ecd1bd"
    },
    {
      "path": "src/cfdx/physics/finite_volume_transport.h",
      "new": false,
      "sha256_before": "aa8ca7a92e34a1eb9471b9f353d1b966f5951a8b0d38210f9a21ea37247b3e46",
      "sha256_after": "a73e1f0ac155d108bbd0b90a3103a67e799c00e7533f3c103dfcf29a7707bc52"
    },
    {
      "path": "src/cfdx/physics/low_storage_time_integration.h",
      "new": false,
      "sha256_before": "c90385391f9c006c3298dd9cf79fb4eb9c9ae36a8dea9674a0b142b83606eaec",
      "sha256_after": "3817aac2660a6aa7eef67bb213eb05126611977326aabc729edc2d904d7b9b7b"
    },
    {
      "path": "src/cfdx/physics/m1_m4_models.h",
      "new": false,
      "sha256_before": "07c36441ec4c359a46770650dc54dab29565c33a388a9699e2c435925020c9ab",
      "sha256_after": "fb7e36c431174c9651ee6dfbd2f00fb8cefb5f51c3f90c06026d3c240a2a2d90"
    },
    {
      "path": "src/cfdx/physics/sst_solver.h",
      "new": false,
      "sha256_before": "aae297132fe91e946baaf15c3ba20b783a990af6c5ce7a5043ce05b24cb43ba5",
      "sha256_after": "03cf25f49262f113592c9e9c5126f7fadc27b272e21080a31e1f4b00c4f6360f"
    },
    {
      "path": "src/cfdx/physics/thermal_zones.h",
      "new": false,
      "sha256_before": "fcbac3e06dd56e23fad83d91ed21387b46fdd6b6d1f0833497f42d7c51e3b999",
      "sha256_after": "1054f91d3b3b027be989398e5afb8f13446cf6d329a0e24746d1c1a9bc7d3d36"
    },
    {
      "path": "src/cfdx/physics/thermophysical_models.h",
      "new": false,
      "sha256_before": "9437a35ab092133e971e1a38fb0c1aba00f24cb43d70cc3a7141a677aed62b03",
      "sha256_after": "709aa89d250a6a3f63fafbd7e26f84a533ef064c8251e19cbe7ed1f423da352a"
    },
    {
      "path": "src/cfdx/physics/transport_models.h",
      "new": false,
      "sha256_before": "db89aca349ec64465d0d8b6325d4dc1981996241310d8773b8ce0ece52d71193",
      "sha256_after": "ab77bc31b27d1fea2b33319b663409c355fb5e524f2c2b3bda82eb0853d94eed"
    },
    {
      "path": "src/cfdx/physics/turbulence_solver.h",
      "new": false,
      "sha256_before": "63b01d1857bf0ffab4592850a7038316e1665483ce6b74813eac59709167f7e1",
      "sha256_after": "a9d576d8b80e791c522d16244c68889e3ec2b3f552a13ad441e44e1c63690742"
    },
    {
      "path": "src/cfdx/thermodynamics/equation_of_state.h",
      "new": false,
      "sha256_before": "d32ad501d57696e279521cb53f735ede77f051af6ce56115e649cdd34160925d",
      "sha256_after": "dd2b80e5bbdcc7d1a8e4718e2cdf2ae6fdc37f7ae21950150ac46b687e307050"
    },
    {
      "path": "tests/audit/affine_mesh_fixture.h",
      "new": true,
      "sha256_before": null,
      "sha256_after": "861abf16346fb172d665c728064d90534a618bb3347f252d110767819981d44e"
    },
    {
      "path": "tests/audit/chain_fixture.h",
      "new": true,
      "sha256_before": null,
      "sha256_after": "e1f8fbcd7ffe57b53dafe69db7dc08481cc3d404a162127d6ca13fa92af1ddf5"
    },
    {
      "path": "tests/audit/test_audit_cht_contract.cpp",
      "new": true,
      "sha256_before": null,
      "sha256_after": "436f4426f105ef5153bc6ce2173947a50fd2d2b8806865f5912773a568ff6eb9"
    },
    {
      "path": "tests/audit/test_audit_cht_time.cpp",
      "new": true,
      "sha256_before": null,
      "sha256_after": "9649894946df7699936aee3e8faa4ac247c4e1389705cc34d3b0da8bc156388a"
    },
    {
      "path": "tests/audit/test_audit_energy_boundary_metric.cpp",
      "new": true,
      "sha256_before": null,
      "sha256_after": "512750e05bc74d0c9e16f363c04caa555c3561f4a1185ea649dae359f2b30072"
    },
    {
      "path": "tests/audit/test_audit_flux_contract.cpp",
      "new": true,
      "sha256_before": null,
      "sha256_after": "f7f315776c7e049923eb0fbb32220082230d165be9309362a63331efb9221b9f"
    },
    {
      "path": "tests/audit/test_audit_numerics.cpp",
      "new": true,
      "sha256_before": null,
      "sha256_after": "8293a41b6fb949247683e0c74077e647f84ccfe8f27ccadf82e959e40ce1c9e7"
    },
    {
      "path": "tests/unit/test_interpolation.cpp",
      "new": false,
      "sha256_before": "b6a49b4dde7d012c519016a360bb346d5c9303535c468403cd388def14565179",
      "sha256_after": "a8baae031103861503fdbb209369c26538ca47856f54397851bc956a94f2631b"
    },
    {
      "path": "tests/unit/test_transport_conservation.cpp",
      "new": false,
      "sha256_before": "980bf545c96d0d47c477e8527c4272be145fe6fad6387339e603f10e652dc2e7",
      "sha256_after": "83998df5d7306aea7dc63740c717f05c9102ba2239ea90146c6a52aa9d1a8436"
    },
    {
      "path": "tests/unit/test_turbulence_qualification_equations.cpp",
      "new": false,
      "sha256_before": "4d9c88095df6afa6476daf5ce57c37564f7b5431ab7c7fd62faea9ea315604e9",
      "sha256_after": "398f5b6f40cad044da427609e1a318be8d4fc1223eea4da8976214b597f558eb"
    },
    {
      "path": "tests/validation/test_cht_validation.cpp",
      "new": false,
      "sha256_before": "0d98e94b9905ebbd1cf475dbdad1b65cc75e9aaeb729b5d37eb6dc78cdb30508",
      "sha256_after": "c92233d7656e36dad9ad65181e15a782163eebf226701e649d996e596e552100"
    },
    {
      "path": "tests/validation/test_level_c_coupled_verification.cpp",
      "new": false,
      "sha256_before": "0eb0f7507a49af049fe41d24140dcfcf154ef68881f7d5f93cb66c4f81121c9d",
      "sha256_after": "47991f454a7e1a6936f11978f6b7bb038676d831f7e833d5c33ac4a05cc25469"
    },
    {
      "path": "tests/validation/test_polyhedral_gradient_campaign.cpp",
      "new": false,
      "sha256_before": "403977700a366ac185c9d257a7331f8deff711b5f12f8c14968c2687580b2b8c",
      "sha256_after": "3341af70b2d88f231d66a6df1d8f75d687a1a655f3b6dd7444903a8b96260804"
    },
    {
      "path": "tests/validation/test_thermal_challenge_vv.cpp",
      "new": false,
      "sha256_before": "296f0611690332f6be72e4878e1fc279c0bf9b3248e7ffa97991f6054f917526",
      "sha256_after": "b3877f9169e051a7cd549a2914d4030dcc2a49f7066d3e85f6c3c77ad6746b4e"
    }
  ]
}
Annexe F. Journaux R2 et contre-exemples embarques
Les journaux des nouveaux tests ci-dessous sont conserves integralement, hors sorties de compilation. Les sorties de solveur internes ont ete capturees par le test pour laisser apparaitre chaque verdict individuel. Les sources exactes des tests sont dans le patch.
F.1 - Nouveau test CHT avant correction (V1)
PASS	case1/converged	34
PASS	case1/T1	6.8212102632969618e-13
PASS	case1/T2	1.6484591469634324e-12
PASS	case1/global_energy	4.0059631544067363e-12
PASS	case1/series_flux	-8.276708945368234e-13
MEASURE	1	l1=1	l2=1	shear1=0	shear2=0	dt=0.01	relax=0.5	T1=300.83894736842171	T2=399.4573684210543	exact1=300.83894736842103	exact2=399.45736842105265
PASS	case2/converged	36
FAIL	case2/T1	0.9709878295082035
FAIL	case2/T2	0.17017312475923063
FAIL	case2/global_energy	0.66666666666341678
FAIL	case2/series_flux	-0.33333333333167731
MEASURE	2	l1=0.40000000000000002	l2=1.6000000000000001	shear1=2	shear2=-0.5	dt=0.01	relax=0.5	T1=304.0340531561414	T2=399.65132890365368	exact1=303.06306532663319	exact2=399.48115577889445
PASS	case3/converged	30
FAIL	case3/T1	0.19438646339523302
FAIL	case3/T2	0.13002008478758853
FAIL	case3/global_energy	0.8856282322657526
FAIL	case3/series_flux	-0.2952094107551867
MEASURE	3	l1=1	l2=2	shear1=3	shear2=0	dt=0.01	relax=0.5	T1=300.90473888630271	T2=399.89323594381841	exact1=300.71035242290748	exact2=399.76321585903082
PASS	case4/converged	8
PASS	case4/T1	5.6843418860808015e-14
PASS	case4/T2	0
PASS	case4/global_energy	7.9254666723875083e-14
PASS	case4/series_flux	-4.5391309123673915e-14
MEASURE	4	l1=1	l2=1	shear1=0	shear2=0	dt=0.01	relax=1	T1=300.83894736842109	T2=399.45736842105265	exact1=300.83894736842103	exact2=399.45736842105265
PASS	case5/converged	9
FAIL	case5/T1	0.97098782951297835
FAIL	case5/T2	0.17017312476002644
FAIL	case5/global_energy	0.66666666666664054
FAIL	case5/series_flux	-0.33333333333333082
MEASURE	5	l1=0.40000000000000002	l2=1.6000000000000001	shear1=2	shear2=-0.5	dt=0.01	relax=1	T1=304.03405315614617	T2=399.65132890365447	exact1=303.06306532663319	exact2=399.48115577889445
PASS	case6/converged	6
FAIL	case6/T1	0.19438646339597199
FAIL	case6/T2	0.13002008478798643
FAIL	case6/global_energy	0.88562823226868426
FAIL	case6/series_flux	-0.29520941075631091
MEASURE	6	l1=1	l2=2	shear1=3	shear2=0	dt=0.01	relax=1	T1=300.90473888630345	T2=399.89323594381881	exact1=300.71035242290748	exact2=399.76321585903082
PASS	case7/converged	47
PASS	case7/T1	8.0149220593739301e-12
PASS	case7/T2	1.8530954548623413e-11
PASS	case7/global_energy	5.1113773914387357e-12
PASS	case7/series_flux	-1.021059448773527e-12
MEASURE	7	l1=1	l2=1	shear1=0	shear2=0	dt=0.10000000000000001	relax=0.5	T1=307.51176470589036	T2=395.15882352943032	exact1=307.51176470588234	exact2=395.15882352941179
PASS	case8/converged	49
FAIL	case8/T1	5.7686548083815978
FAIL	case8/T2	1.4149530662068059
FAIL	case8/global_energy	0.66666666666073104
FAIL	case8/series_flux	-0.33333333333036858
MEASURE	8	l1=0.40000000000000002	l2=1.6000000000000001	shear1=2	shear2=-0.5	dt=0.10000000000000001	relax=0.5	T1=329.01674641143501	T2=397.45693779903121	exact1=323.24809160305341	exact2=396.0419847328244
PASS	case9/converged	36
FAIL	case9/T1	1.7481830434069252
FAIL	case9/T2	1.2018758423423037
FAIL	case9/global_energy	0.88562823226525889
FAIL	case9/series_flux	-0.29520941075507878
MEASURE	9	l1=1	l2=2	shear1=3	shear2=0	dt=0.10000000000000001	relax=0.5	T1=308.32981569646813	T2=399.00799829132188	exact1=306.58163265306121	exact2=397.80612244897958
PASS	case10/converged	17
PASS	case10/T1	5.1159076974727213e-13
PASS	case10/T2	1.0231815394945443e-12
PASS	case10/global_energy	2.9590051512693915e-13
PASS	case10/series_flux	-6.505757901078594e-14
MEASURE	10	l1=1	l2=1	shear1=0	shear2=0	dt=0.10000000000000001	relax=1	T1=307.51176470588285	T2=395.15882352941281	exact1=307.51176470588234	exact2=395.15882352941179
PASS	case11/converged	18
FAIL	case11/T1	5.7686548084279252
FAIL	case11/T2	1.4149530662182315
FAIL	case11/global_energy	0.66666666666644248
FAIL	case11/series_flux	-0.33333333333321463
MEASURE	11	l1=0.40000000000000002	l2=1.6000000000000001	shear1=2	shear2=-0.5	dt=0.10000000000000001	relax=1	T1=329.01674641148134	T2=397.45693779904263	exact1=323.24809160305341	exact2=396.0419847328244
PASS	case12/converged	10
FAIL	case12/T1	1.7481830434133485
FAIL	case12/T2	1.2018758423466807
FAIL	case12/global_energy	0.8856282322685135
FAIL	case12/series_flux	-0.29520941075617052
MEASURE	12	l1=1	l2=2	shear1=3	shear2=0	dt=0.10000000000000001	relax=1	T1=308.32981569647455	T2=399.00799829132626	exact1=306.58163265306121	exact2=397.80612244897958
PASS	case13/converged	6
FAIL	case13/T1	5.0253066774530453e-07
FAIL	case13/T2	0.00033280507875588228
FAIL	case13/global_energy	0.099942029759858997
FAIL	case13/series_flux	-9.9842191173471487e-05
MEASURE	13	l1=1	l2=1	shear1=2	shear2=-0.5	dt=0.001	relax=1	T1=300.0059994202976	T2=399.99800019323413	exact1=300.00599891776693	exact2=399.99766738815538
FAIL	reject/material1_mismatch	0
FAIL	unchanged/material1_mismatch	0
FAIL	reject/material2_mismatch	0
FAIL	unchanged/material2_mismatch	0
FAIL	reject/different_dt	0
FAIL	unchanged/different_dt	0
FAIL	reject/mixed_steady_transient	0
FAIL	unchanged/mixed_steady_transient	0
FAIL	reject/nan_conductivity	0
FAIL	unchanged/nan_conductivity	0
FAIL	reject/inf_conductivity	0
FAIL	unchanged/inf_conductivity	0
FAIL	reject/nan_balance_tolerance	0
FAIL	unchanged/nan_balance_tolerance	0
FAIL	reject/nan_temperature_tolerance	0
FAIL	unchanged/nan_temperature_tolerance	0
FAIL	reject/inf_matching_tolerance	0
FAIL	unchanged/inf_matching_tolerance	0
FAIL	reject/nan_area_tolerance	0
FAIL	unchanged/nan_area_tolerance	0
FAIL	reject/permeable_interface	0
FAIL	unchanged/permeable_interface	0
FAIL	reject/inward_normal	0
FAIL	unchanged/inward_normal	0
FAIL	reject/overlapping_regions	0
FAIL	unchanged/overlapping_regions	0
FAIL	reject/nonconformal_vertices	0
FAIL	unchanged/nonconformal_vertices	0
FAIL	reject/empty_interface	0
FAIL	unchanged/empty_interface	0
FAIL	reject/duplicate_interface	0
FAIL	unchanged/duplicate_interface	0
OTHER_EXCEPTION	nan_cell_centre	CHT interface face-to-cell distance is invalid
FAIL	reject/nan_cell_centre	0
PASS	unchanged/nan_cell_centre	0
OTHER_EXCEPTION	nan_source	CHT region energy solve did not converge
FAIL	reject/nan_source	0
FAIL	unchanged/nan_source	0
SUMMARY checks=101 failures=71 max_T_error=5.7686548084279252 max_energy_normalized=0.88562823226868426
F.2 - Nouveau test CHT apres correction (R2)
PASS	case1/converged	34
PASS	case1/T1	6.8212102632969618e-13
PASS	case1/T2	1.6484591469634324e-12
PASS	case1/global_energy	4.0059631544067363e-12
PASS	case1/series_flux	-8.276708945368234e-13
MEASURE	1	l1=1	l2=1	shear1=0	shear2=0	dt=0.01	relax=0.5	T1=300.83894736842171	T2=399.4573684210543	exact1=300.83894736842103	exact2=399.45736842105265
PASS	case2/converged	37
PASS	case2/T1	-3.694822225952521e-12
PASS	case2/T2	-6.2527760746888816e-13
PASS	case2/global_energy	-2.5130217379283304e-12
PASS	case2/series_flux	1.2484340149760786e-12
MEASURE	2	l1=0.40000000000000002	l2=1.6000000000000001	shear1=2	shear2=-0.5	dt=0.01	relax=0.5	T1=303.0630653266295	T2=399.48115577889382	exact1=303.06306532663319	exact2=399.48115577889445
PASS	case3/converged	31
PASS	case3/T1	-9.0949470177292824e-13
PASS	case3/T2	-6.2527760746888816e-13
PASS	case3/global_energy	-4.3040213747176063e-12
PASS	case3/series_flux	1.3990221478519578e-12
MEASURE	3	l1=1	l2=2	shear1=3	shear2=0	dt=0.01	relax=0.5	T1=300.71035242290657	T2=399.7632158590302	exact1=300.71035242290748	exact2=399.76321585903082
PASS	case4/converged	8
PASS	case4/T1	5.6843418860808015e-14
PASS	case4/T2	0
PASS	case4/global_energy	7.9254666723875083e-14
PASS	case4/series_flux	-4.5391309123673915e-14
MEASURE	4	l1=1	l2=1	shear1=0	shear2=0	dt=0.01	relax=1	T1=300.83894736842109	T2=399.45736842105265	exact1=300.83894736842103	exact2=399.45736842105265
PASS	case5/converged	9
PASS	case5/T1	-1.7053025658242404e-13
PASS	case5/T2	0
PASS	case5/global_energy	-9.8219260906224467e-14
PASS	case5/series_flux	4.8696944482918013e-14
MEASURE	5	l1=0.40000000000000002	l2=1.6000000000000001	shear1=2	shear2=-0.5	dt=0.01	relax=1	T1=303.06306532663302	T2=399.48115577889445	exact1=303.06306532663319	exact2=399.48115577889445
PASS	case6/converged	7
PASS	case6/T1	0
PASS	case6/T2	0
PASS	case6/global_energy	-8.6080427494352361e-14
PASS	case6/series_flux	1.8937694048757522e-14
MEASURE	6	l1=1	l2=2	shear1=3	shear2=0	dt=0.01	relax=1	T1=300.71035242290748	T2=399.76321585903082	exact1=300.71035242290748	exact2=399.76321585903082
PASS	case7/converged	47
PASS	case7/T1	8.0149220593739301e-12
PASS	case7/T2	1.8530954548623413e-11
PASS	case7/global_energy	5.1113773914387357e-12
PASS	case7/series_flux	-1.021059448773527e-12
MEASURE	7	l1=1	l2=1	shear1=0	shear2=0	dt=0.10000000000000001	relax=0.5	T1=307.51176470589036	T2=395.15882352943032	exact1=307.51176470588234	exact2=395.15882352941179
PASS	case8/converged	53
PASS	case8/T1	-5.8093974075745791e-11
PASS	case8/T2	-1.4210854715202004e-11
PASS	case8/global_energy	-6.3096802950215965e-12
PASS	case8/series_flux	3.1577294097928683e-12
MEASURE	8	l1=0.40000000000000002	l2=1.6000000000000001	shear1=2	shear2=-0.5	dt=0.10000000000000001	relax=0.5	T1=323.24809160299532	T2=396.04198473281019	exact1=323.24809160305341	exact2=396.0419847328244
PASS	case9/converged	38
PASS	case9/T1	-6.3096194935496897e-12
PASS	case9/T2	-4.3200998334214091e-12
PASS	case9/global_energy	-3.1778904638289963e-12
PASS	case9/series_flux	1.0623734554263626e-12
MEASURE	9	l1=1	l2=2	shear1=3	shear2=0	dt=0.10000000000000001	relax=0.5	T1=306.5816326530549	T2=397.80612244897526	exact1=306.58163265306121	exact2=397.80612244897958
PASS	case10/converged	17
PASS	case10/T1	5.1159076974727213e-13
PASS	case10/T2	1.0231815394945443e-12
PASS	case10/global_energy	2.9590051512693915e-13
PASS	case10/series_flux	-6.505757901078594e-14
MEASURE	10	l1=1	l2=1	shear1=0	shear2=0	dt=0.10000000000000001	relax=1	T1=307.51176470588285	T2=395.15882352941281	exact1=307.51176470588234	exact2=395.15882352941179
PASS	case11/converged	20
PASS	case11/T1	-5.9685589803848416e-12
PASS	case11/T2	-1.4210854715202004e-12
PASS	case11/global_energy	-6.4485210609040001e-13
PASS	case11/series_flux	3.2594002068555703e-13
MEASURE	11	l1=0.40000000000000002	l2=1.6000000000000001	shear1=2	shear2=-0.5	dt=0.10000000000000001	relax=1	T1=323.24809160304744	T2=396.04198473282298	exact1=323.24809160305341	exact2=396.0419847328244
PASS	case12/converged	11
PASS	case12/T1	-2.2737367544323206e-13
PASS	case12/T2	-1.7053025658242404e-13
PASS	case12/global_energy	-1.3085431321649084e-13
PASS	case12/series_flux	4.112564129661141e-14
MEASURE	12	l1=1	l2=2	shear1=3	shear2=0	dt=0.10000000000000001	relax=1	T1=306.58163265306098	T2=397.80612244897941	exact1=306.58163265306121	exact2=397.80612244897958
PASS	case13/converged	6
PASS	case13/T1	0
PASS	case13/T2	-5.6843418860808015e-14
PASS	case13/global_energy	-1.7507773009128869e-11
PASS	case13/series_flux	1.2204681709704346e-12
MEASURE	13	l1=1	l2=1	shear1=2	shear2=-0.5	dt=0.001	relax=1	T1=300.00599891776693	T2=399.99766738815532	exact1=300.00599891776693	exact2=399.99766738815538
PASS	reject/material1_mismatch	0
PASS	unchanged/material1_mismatch	0
PASS	reject/material2_mismatch	0
PASS	unchanged/material2_mismatch	0
PASS	reject/different_dt	0
PASS	unchanged/different_dt	0
PASS	reject/mixed_steady_transient	0
PASS	unchanged/mixed_steady_transient	0
PASS	reject/nan_conductivity	0
PASS	unchanged/nan_conductivity	0
PASS	reject/inf_conductivity	0
PASS	unchanged/inf_conductivity	0
PASS	reject/nan_balance_tolerance	0
PASS	unchanged/nan_balance_tolerance	0
PASS	reject/nan_temperature_tolerance	0
PASS	unchanged/nan_temperature_tolerance	0
PASS	reject/inf_matching_tolerance	0
PASS	unchanged/inf_matching_tolerance	0
PASS	reject/nan_area_tolerance	0
PASS	unchanged/nan_area_tolerance	0
PASS	reject/permeable_interface	0
PASS	unchanged/permeable_interface	0
PASS	reject/inward_normal	0
PASS	unchanged/inward_normal	0
PASS	reject/overlapping_regions	0
PASS	unchanged/overlapping_regions	0
PASS	reject/nonconformal_vertices	0
PASS	unchanged/nonconformal_vertices	0
PASS	reject/empty_interface	0
PASS	unchanged/empty_interface	0
PASS	reject/duplicate_interface	0
PASS	unchanged/duplicate_interface	0
PASS	reject/nan_cell_centre	0
PASS	unchanged/nan_cell_centre	0
PASS	reject/nan_source	0
PASS	unchanged/nan_source	0
SUMMARY checks=101 failures=0 max_T_error=5.8093974075745791e-11 max_energy_normalized=1.7507773009128869e-11
F.3 - Diagnostic energie avant correction (V1)
PASS	matrix	shear=0	override=0	value=0
PASS	diagnostic	shear=0	override=0	value=7.9580786405131203e-16
PASS	matrix	shear=0	override=1	value=0
PASS	diagnostic	shear=0	override=1	value=7.9580786405131203e-16
PASS	matrix	shear=1	override=0	value=0
FAIL	diagnostic	shear=1	override=0	value=0.29289321881345193
PASS	matrix	shear=1	override=1	value=0
FAIL	diagnostic	shear=1	override=1	value=0.29289321881345193
PASS	matrix	shear=3	override=0	value=0
FAIL	diagnostic	shear=3	override=0	value=0.68377223398316178
PASS	matrix	shear=3	override=1	value=0
FAIL	diagnostic	shear=3	override=1	value=0.68377223398316178
SUMMARY checks=12 failures=4
F.4 - Diagnostic energie apres correction (R2)
PASS	matrix	shear=0	override=0	value=0
PASS	diagnostic	shear=0	override=0	value=7.9580786405131203e-16
PASS	matrix	shear=0	override=1	value=0
PASS	diagnostic	shear=0	override=1	value=7.9580786405131203e-16
PASS	matrix	shear=1	override=0	value=0
PASS	diagnostic	shear=1	override=0	value=7.9580786405131203e-16
PASS	matrix	shear=1	override=1	value=0
PASS	diagnostic	shear=1	override=1	value=7.9580786405131203e-16
PASS	matrix	shear=3	override=0	value=0
PASS	diagnostic	shear=3	override=0	value=7.9580786405131203e-16
PASS	matrix	shear=3	override=1	value=0
PASS	diagnostic	shear=3	override=1	value=7.9580786405131203e-16
SUMMARY checks=12 failures=0
F.5 - Premier passage des trois anciennes fixtures
Internal ctest changing into directory: /mnt/data/cfdx_md_r2/build
Test project /mnt/data/cfdx_md_r2/build
    Start 161: test_level_c_coupled_verification
1/3 Test #161: test_level_c_coupled_verification ...***Failed    0.00 sec
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=0.00725088 source_relative=0.00725088
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.00362544 backward_error=1.81788e-07 linear_relative=0 dT_relative=3.63576e-07 energy_balance=1.81788e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.00181272 backward_error=9.08939e-08 linear_relative=0 dT_relative=1.81788e-07 energy_balance=9.08939e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.000906359 backward_error=4.54469e-08 linear_relative=0 dT_relative=9.08939e-08 energy_balance=4.54469e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.00045318 backward_error=2.27235e-08 linear_relative=0 dT_relative=4.54469e-08 energy_balance=2.27235e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.00022659 backward_error=1.13617e-08 linear_relative=0 dT_relative=2.27235e-08 energy_balance=1.13617e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.000113295 backward_error=5.68087e-09 linear_relative=0 dT_relative=1.13617e-08 energy_balance=5.68087e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=5.66475e-05 backward_error=2.84043e-09 linear_relative=0 dT_relative=5.68087e-09 energy_balance=2.84043e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=2.83237e-05 backward_error=1.42022e-09 linear_relative=0 dT_relative=2.84043e-09 energy_balance=1.42022e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=1.41619e-05 backward_error=7.10109e-10 linear_relative=0 dT_relative=1.42022e-09 energy_balance=7.10109e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=7.08093e-06 backward_error=3.55054e-10 linear_relative=0 dT_relative=7.10109e-10 energy_balance=3.55054e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=3.54047e-06 backward_error=1.77527e-10 linear_relative=0 dT_relative=3.55054e-10 energy_balance=1.77527e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=1.77023e-06 backward_error=8.87636e-11 linear_relative=0 dT_relative=1.77527e-10 energy_balance=8.87636e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=8.85117e-07 backward_error=4.43818e-11 linear_relative=0 dT_relative=8.87636e-11 energy_balance=4.43818e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=1 dT=0.000290825 qrad_delta=0.00725088 qrad_relative=0.00725088 energy_balance=1.95884e-07
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=0.000134958 source_relative=0.000134958
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.00181301 backward_error=9.09084e-08 linear_relative=0 dT_relative=1.81817e-07 energy_balance=9.09084e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.000906504 backward_error=4.54542e-08 linear_relative=0 dT_relative=9.09084e-08 energy_balance=4.54542e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.000453252 backward_error=2.27271e-08 linear_relative=0 dT_relative=4.54542e-08 energy_balance=2.27271e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.000226626 backward_error=1.13635e-08 linear_relative=0 dT_relative=2.27271e-08 energy_balance=1.13635e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.000113313 backward_error=5.68177e-09 linear_relative=0 dT_relative=1.13635e-08 energy_balance=5.68177e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=5.66565e-05 backward_error=2.84089e-09 linear_relative=0 dT_relative=5.68177e-09 energy_balance=2.84089e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=2.83282e-05 backward_error=1.42044e-09 linear_relative=0 dT_relative=2.84089e-09 energy_balance=1.42044e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=1.41641e-05 backward_error=7.10222e-10 linear_relative=0 dT_relative=1.42044e-09 energy_balance=7.10222e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=7.08206e-06 backward_error=3.55111e-10 linear_relative=0 dT_relative=7.10222e-10 energy_balance=3.55111e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=3.54103e-06 backward_error=1.77555e-10 linear_relative=0 dT_relative=3.55111e-10 energy_balance=1.77555e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=1.77051e-06 backward_error=8.87776e-11 linear_relative=0 dT_relative=1.77555e-10 energy_balance=8.87776e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=8.85257e-07 backward_error=4.43888e-11 linear_relative=0 dT_relative=8.87777e-11 energy_balance=4.43888e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=2 dT=0.000145418 qrad_delta=0.000134958 qrad_relative=0.000134958 energy_balance=9.79687e-08
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=6.74817e-05 source_relative=6.74817e-05
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.000906759 backward_error=4.5467e-08 linear_relative=0 dT_relative=9.0934e-08 energy_balance=4.5467e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.000453379 backward_error=2.27335e-08 linear_relative=0 dT_relative=4.5467e-08 energy_balance=2.27335e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.00022669 backward_error=1.13667e-08 linear_relative=0 dT_relative=2.27335e-08 energy_balance=1.13667e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.000113345 backward_error=5.68337e-09 linear_relative=0 dT_relative=1.13667e-08 energy_balance=5.68337e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=5.66724e-05 backward_error=2.84169e-09 linear_relative=0 dT_relative=5.68337e-09 energy_balance=2.84169e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=2.83362e-05 backward_error=1.42084e-09 linear_relative=0 dT_relative=2.84169e-09 energy_balance=1.42084e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=1.41681e-05 backward_error=7.10421e-10 linear_relative=0 dT_relative=1.42084e-09 energy_balance=7.10421e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=7.08405e-06 backward_error=3.55211e-10 linear_relative=0 dT_relative=7.10421e-10 energy_balance=3.55211e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=3.54203e-06 backward_error=1.77605e-10 linear_relative=0 dT_relative=3.55211e-10 energy_balance=1.77605e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=1.77101e-06 backward_error=8.88027e-11 linear_relative=0 dT_relative=1.77605e-10 energy_balance=8.88027e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=8.85506e-07 backward_error=4.44013e-11 linear_relative=0 dT_relative=8.88027e-11 energy_balance=4.44013e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=3 dT=7.27117e-05 qrad_delta=6.74817e-05 qrad_relative=6.74817e-05 energy_balance=4.90092e-08
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=3.37421e-05 source_relative=3.37421e-05
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.000453618 backward_error=2.27454e-08 linear_relative=0 dT_relative=4.54909e-08 energy_balance=2.27454e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.000226809 backward_error=1.13727e-08 linear_relative=0 dT_relative=2.27454e-08 energy_balance=1.13727e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.000113404 backward_error=5.68636e-09 linear_relative=0 dT_relative=1.13727e-08 energy_balance=5.68636e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=5.67022e-05 backward_error=2.84318e-09 linear_relative=0 dT_relative=5.68636e-09 energy_balance=2.84318e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=2.83511e-05 backward_error=1.42159e-09 linear_relative=0 dT_relative=2.84318e-09 energy_balance=1.42159e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=1.41755e-05 backward_error=7.10795e-10 linear_relative=0 dT_relative=1.42159e-09 energy_balance=7.10795e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=7.08778e-06 backward_error=3.55397e-10 linear_relative=0 dT_relative=7.10795e-10 energy_balance=3.55397e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=3.54389e-06 backward_error=1.77699e-10 linear_relative=0 dT_relative=3.55397e-10 energy_balance=1.77699e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=1.77194e-06 backward_error=8.88494e-11 linear_relative=0 dT_relative=1.77699e-10 energy_balance=8.88494e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=8.85973e-07 backward_error=4.44247e-11 linear_relative=0 dT_relative=8.88493e-11 energy_balance=4.44247e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=4 dT=3.63572e-05 qrad_delta=3.37421e-05 qrad_relative=3.37421e-05 energy_balance=2.45286e-08
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=1.68717e-05 source_relative=1.68717e-05
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.000227039 backward_error=1.13842e-08 linear_relative=0 dT_relative=2.27685e-08 energy_balance=1.13842e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.000113519 backward_error=5.69212e-09 linear_relative=0 dT_relative=1.13842e-08 energy_balance=5.69212e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=5.67597e-05 backward_error=2.84606e-09 linear_relative=0 dT_relative=5.69212e-09 energy_balance=2.84606e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=2.83798e-05 backward_error=1.42303e-09 linear_relative=0 dT_relative=2.84606e-09 energy_balance=1.42303e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=1.41899e-05 backward_error=7.11515e-10 linear_relative=0 dT_relative=1.42303e-09 energy_balance=7.11515e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=7.09496e-06 backward_error=3.55758e-10 linear_relative=0 dT_relative=7.11515e-10 energy_balance=3.55758e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=3.54748e-06 backward_error=1.77879e-10 linear_relative=0 dT_relative=3.55758e-10 energy_balance=1.77879e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=1.77374e-06 backward_error=8.89393e-11 linear_relative=0 dT_relative=1.77879e-10 energy_balance=8.89393e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=8.86868e-07 backward_error=4.44696e-11 linear_relative=0 dT_relative=8.89394e-11 energy_balance=4.44696e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=5 dT=1.81792e-05 qrad_delta=1.68717e-05 qrad_relative=1.68717e-05 energy_balance=1.22878e-08
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=8.43614e-06 source_relative=8.43614e-06
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.000113745 backward_error=5.70345e-09 linear_relative=0 dT_relative=1.14069e-08 energy_balance=5.70345e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=5.68727e-05 backward_error=2.85173e-09 linear_relative=0 dT_relative=5.70345e-09 energy_balance=2.85173e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=2.84363e-05 backward_error=1.42586e-09 linear_relative=0 dT_relative=2.85173e-09 energy_balance=1.42586e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=1.42182e-05 backward_error=7.12931e-10 linear_relative=0 dT_relative=1.42586e-09 energy_balance=7.12931e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=7.10908e-06 backward_error=3.56466e-10 linear_relative=0 dT_relative=7.12931e-10 energy_balance=3.56466e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=3.55454e-06 backward_error=1.78233e-10 linear_relative=0 dT_relative=3.56466e-10 energy_balance=1.78233e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=1.77727e-06 backward_error=8.91164e-11 linear_relative=0 dT_relative=1.78233e-10 energy_balance=8.91164e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=8.88636e-07 backward_error=4.45583e-11 linear_relative=0 dT_relative=8.91163e-11 energy_balance=4.45583e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=6 dT=9.08988e-06 qrad_delta=8.43614e-06 qrad_relative=8.43614e-06 energy_balance=6.16729e-09
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=4.2182e-06 source_relative=4.2182e-06
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=5.70969e-05 backward_error=2.86297e-09 linear_relative=0 dT_relative=5.72594e-09 energy_balance=2.86297e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=2.85485e-05 backward_error=1.43149e-09 linear_relative=0 dT_relative=2.86297e-09 energy_balance=1.43149e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=1.42742e-05 backward_error=7.15743e-10 linear_relative=0 dT_relative=1.43149e-09 energy_balance=7.15743e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=7.13712e-06 backward_error=3.57871e-10 linear_relative=0 dT_relative=7.15743e-10 energy_balance=3.57871e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=3.56856e-06 backward_error=1.78936e-10 linear_relative=0 dT_relative=3.57871e-10 energy_balance=1.78936e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=1.78428e-06 backward_error=8.94679e-11 linear_relative=0 dT_relative=1.78936e-10 energy_balance=8.94679e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=8.92142e-07 backward_error=4.4734e-11 linear_relative=0 dT_relative=8.94678e-11 energy_balance=4.4734e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=7 dT=4.54497e-06 qrad_delta=4.2182e-06 qrad_relative=4.2182e-06 energy_balance=3.10699e-09
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=2.10911e-06 source_relative=2.10911e-06
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=2.87726e-05 backward_error=1.44272e-09 linear_relative=0 dT_relative=2.88544e-09 energy_balance=1.44272e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=1.43863e-05 backward_error=7.21361e-10 linear_relative=0 dT_relative=1.44272e-09 energy_balance=7.21361e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=7.19314e-06 backward_error=3.60681e-10 linear_relative=0 dT_relative=7.21361e-10 energy_balance=3.60681e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=3.59657e-06 backward_error=1.8034e-10 linear_relative=0 dT_relative=3.6068e-10 energy_balance=1.8034e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=1.79828e-06 backward_error=9.01701e-11 linear_relative=0 dT_relative=1.8034e-10 energy_balance=9.01701e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=8.99143e-07 backward_error=4.50851e-11 linear_relative=0 dT_relative=9.01701e-11 energy_balance=4.50851e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=8 dT=2.27229e-06 qrad_delta=2.10911e-06 qrad_relative=2.10911e-06 energy_balance=1.57696e-09
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=1.05446e-06 source_relative=1.05446e-06
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=1.46116e-05 backward_error=7.32659e-10 linear_relative=0 dT_relative=1.46532e-09 energy_balance=7.32659e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=7.30579e-06 backward_error=3.66329e-10 linear_relative=0 dT_relative=7.32659e-10 energy_balance=3.66329e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=3.6529e-06 backward_error=1.83165e-10 linear_relative=0 dT_relative=3.66329e-10 energy_balance=1.83165e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=1.82645e-06 backward_error=9.15823e-11 linear_relative=0 dT_relative=1.83165e-10 energy_balance=9.15823e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=9.13225e-07 backward_error=4.57912e-11 linear_relative=0 dT_relative=9.15824e-11 energy_balance=4.57912e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=9 dT=1.13562e-06 qrad_delta=1.05446e-06 qrad_relative=1.05446e-06 energy_balance=8.12276e-10
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=5.2699e-07 source_relative=5.2699e-07
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=7.53437e-06 backward_error=3.7779e-10 linear_relative=0 dT_relative=7.55581e-10 energy_balance=3.7779e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=3.76718e-06 backward_error=1.88895e-10 linear_relative=0 dT_relative=3.7779e-10 energy_balance=1.88895e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=1.88359e-06 backward_error=9.44476e-11 linear_relative=0 dT_relative=1.88895e-10 energy_balance=9.44476e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=9.41794e-07 backward_error=4.72237e-11 linear_relative=0 dT_relative=9.44475e-11 energy_balance=4.72237e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=10 dT=5.66686e-07 qrad_delta=5.2699e-07 qrad_relative=5.2699e-07 energy_balance=4.30651e-10
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=2.62973e-07 source_relative=2.62973e-07
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=4.00276e-06 backward_error=2.00708e-10 linear_relative=0 dT_relative=4.01415e-10 energy_balance=2.00708e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=2.00138e-06 backward_error=1.00354e-10 linear_relative=0 dT_relative=2.00708e-10 energy_balance=1.00354e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=1.00069e-06 backward_error=5.01769e-11 linear_relative=0 dT_relative=1.00354e-10 energy_balance=5.01769e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=5.00346e-07 backward_error=2.50885e-11 linear_relative=0 dT_relative=5.01769e-11 energy_balance=2.50885e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=11 dT=3.01062e-07 qrad_delta=2.62973e-07 qrad_relative=2.62973e-07 energy_balance=2.28791e-10
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=1.39709e-07 source_relative=1.39709e-07
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=2.12654e-06 backward_error=1.06629e-10 linear_relative=0 dT_relative=2.13259e-10 energy_balance=1.06629e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=1.06327e-06 backward_error=5.33147e-11 linear_relative=0 dT_relative=1.06629e-10 energy_balance=5.33147e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=5.31634e-07 backward_error=2.66574e-11 linear_relative=0 dT_relative=5.33148e-11 energy_balance=2.66574e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=12 dT=1.49281e-07 qrad_delta=1.39709e-07 qrad_relative=1.39709e-07 energy_balance=1.28213e-10
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=6.92746e-08 source_relative=6.92746e-08
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=1.19621e-06 backward_error=5.99808e-11 linear_relative=0 dT_relative=1.19962e-10 energy_balance=5.99808e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=5.98106e-07 backward_error=2.99904e-11 linear_relative=0 dT_relative=5.99808e-11 energy_balance=2.99904e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=13 dT=7.19771e-08 qrad_delta=6.92746e-08 qrad_relative=6.92746e-08 energy_balance=7.96196e-11
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=3.34013e-08 source_relative=3.34013e-08
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=7.47648e-07 backward_error=3.74888e-11 linear_relative=0 dT_relative=7.49776e-11 energy_balance=3.74888e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=14 dT=2.9991e-08 qrad_delta=3.34013e-08 qrad_relative=3.34013e-08 energy_balance=5.91357e-11
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=1.39175e-08 source_relative=1.39175e-08
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=5.60744e-07 backward_error=2.8117e-11 linear_relative=0 dT_relative=5.6234e-11 energy_balance=2.8117e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=15 dT=2.24936e-08 qrad_delta=1.39175e-08 qrad_relative=1.39175e-08 energy_balance=4.43522e-11
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=1.04383e-08 source_relative=1.04383e-08
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=4.20563e-07 backward_error=2.1088e-11 linear_relative=0 dT_relative=4.21758e-11 energy_balance=2.1088e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=16 dT=1.68703e-08 qrad_delta=1.04383e-08 qrad_relative=1.04383e-08 energy_balance=3.32646e-11
RADIATION_RESIDUAL: iteration=1 intensity_delta=7579.4 intensity_relative=0.0815773 source_delta=7.82874e-09 source_relative=7.82874e-09
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=3.15427e-07 backward_error=1.58162e-11 linear_relative=0 dT_relative=3.16323e-11 energy_balance=1.58162e-11 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=17 dT=1.2653e-08 qrad_delta=7.82874e-09 qrad_relative=7.82874e-09 energy_balance=2.49487e-11
LEVEL_C_COUPLED_VERIFICATION: FAIL: CHT interface outward normals must oppose

    Start 162: test_cht_validation
2/3 Test #162: test_cht_validation .................***Failed    0.00 sec
Running 3 test cases...
  FAIL: two_region_cht_independent_energy_conservation - CHT interface outward normals must oppose
  PASS: cht_independent_energy_balance_detects_nonconservative_state
  PASS: cht_face_override_is_used_by_independent_energy_balance

Results: 2 passed, 1 failed

    Start 173: test_thermal_challenge_vv
3/3 Test #173: test_thermal_challenge_vv ...........***Failed    0.01 sec
Running 8 test cases...
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0 backward_error=0 linear_relative=0 dT_relative=0.121569 energy_balance=0 accepted_state=1
THERMAL_CHALLENGE: case=multicell_conduction iterations=1 residual=0 backward_error=0 energy_imbalance=0
  PASS: thermal_multicell_analytical_conduction
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=8.88178e-16 backward_error=8.88178e-16 linear_relative=1.91874e-16 dT_relative=0.5 energy_balance=4.93432e-17 accepted_state=1
THERMAL_CHALLENGE: case=mms_quadratic iterations=1 residual=8.88178e-16 backward_error=8.88178e-16 energy_imbalance=4.93432e-17
THERMAL_CHALLENGE: case=mms_quadratic n=4 L_inf=0.015625
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=7.9492e-14 backward_error=3.13083e-14 linear_relative=8.96874e-15 dT_relative=0.625 energy_balance=8.77605e-15 accepted_state=1
THERMAL_CHALLENGE: case=mms_quadratic iterations=1 residual=7.9492e-14 backward_error=3.13083e-14 energy_imbalance=8.77605e-15
THERMAL_CHALLENGE: case=mms_quadratic n=8 L_inf=0.00390625 observed_order=2
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=2.1922e-11 backward_error=6.77902e-12 linear_relative=1.73274e-12 dT_relative=0.6875 energy_balance=2.23802e-12 accepted_state=1
THERMAL_CHALLENGE: case=mms_quadratic iterations=1 residual=2.1922e-11 backward_error=6.77902e-12 energy_imbalance=2.23802e-12
THERMAL_CHALLENGE: case=mms_quadratic n=16 L_inf=0.000976563 observed_order=2
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=3.55271e-15 backward_error=3.55271e-15 linear_relative=1.74609e-16 dT_relative=0.71875 energy_balance=3.82012e-17 accepted_state=1
THERMAL_CHALLENGE: case=mms_quadratic iterations=1 residual=3.55271e-15 backward_error=3.55271e-15 energy_imbalance=3.82012e-17
THERMAL_CHALLENGE: case=mms_quadratic n=32 L_inf=0.000244141 observed_order=2
  PASS: thermal_mms_quadratic_refinement
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0 backward_error=0 linear_relative=0 dT_relative=0.0204082 energy_balance=0 accepted_state=1
THERMAL_CHALLENGE: case=volumetric_generation iterations=1 residual=0 backward_error=0 energy_imbalance=0
THERMAL_CHALLENGE: case=volumetric_generation n=4 L_inf=0.390625 energy_imbalance=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=9.09495e-13 backward_error=4.67607e-17 linear_relative=1.63879e-16 dT_relative=0.0204082 energy_balance=3.14161e-17 accepted_state=1
THERMAL_CHALLENGE: case=volumetric_generation iterations=1 residual=9.09495e-13 backward_error=4.67607e-17 energy_imbalance=3.14161e-17
THERMAL_CHALLENGE: case=volumetric_generation n=8 L_inf=0.0976563 observed_order=2 energy_imbalance=3.14161e-17
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=1.09139e-11 backward_error=2.78417e-16 linear_relative=9.47082e-16 dT_relative=0.0204082 energy_balance=0 accepted_state=1
THERMAL_CHALLENGE: case=volumetric_generation iterations=1 residual=1.09139e-11 backward_error=2.78417e-16 energy_imbalance=0
THERMAL_CHALLENGE: case=volumetric_generation n=16 L_inf=0.0244141 observed_order=2 energy_imbalance=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0 backward_error=0 linear_relative=0 dT_relative=0.0204082 energy_balance=0 accepted_state=1
THERMAL_CHALLENGE: case=volumetric_generation iterations=1 residual=0 backward_error=0 energy_imbalance=0
THERMAL_CHALLENGE: case=volumetric_generation n=32 L_inf=0.00610352 observed_order=2 energy_imbalance=0
  PASS: thermal_volumetric_generation_refinement
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0 backward_error=0 linear_relative=0 dT_relative=0.0204082 energy_balance=0 accepted_state=1
THERMAL_CHALLENGE: case=volumetric_power_sweep iterations=1 residual=0 backward_error=0 energy_imbalance=0
THERMAL_POWER_STUDY: qvol=100 W/m3 Tmax=306.25 dTmax=6.25 exact_dTmax_discrete=6.25 exact_dTmax_continuous=6.25 expected_discretization_error=0.00610352 observed_continuous_cell_error=0.00610352 generated_power=100 expected_power=100 Tmin=300.391
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=3.63798e-12 backward_error=4.03435e-17 linear_relative=1.33872e-16 dT_relative=0.172414 energy_balance=0 accepted_state=1
THERMAL_CHALLENGE: case=volumetric_power_sweep iterations=1 residual=3.63798e-12 backward_error=4.03435e-17 energy_imbalance=0
THERMAL_POWER_STUDY: qvol=1000 W/m3 Tmax=362.5 dTmax=62.5 exact_dTmax_discrete=62.5 exact_dTmax_continuous=62.5 expected_discretization_error=0.0610352 observed_continuous_cell_error=0.0610352 generated_power=1000 expected_power=1000 Tmin=303.906
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=1.45519e-11 backward_error=9.43092e-17 linear_relative=7.60182e-16 dT_relative=0.675676 energy_balance=1.53631e-17 accepted_state=1
THERMAL_CHALLENGE: case=volumetric_power_sweep iterations=1 residual=1.45519e-11 backward_error=9.43092e-17 energy_imbalance=1.53631e-17
THERMAL_POWER_STUDY: qvol=10000 W/m3 Tmax=925 dTmax=625 exact_dTmax_discrete=625 exact_dTmax_continuous=625 expected_discretization_error=0.610352 observed_continuous_cell_error=0.610352 generated_power=10000 expected_power=10000 Tmin=339.062
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=1.16415e-10 backward_error=7.83555e-17 linear_relative=4.81118e-15 dT_relative=0.954198 energy_balance=8.67838e-18 accepted_state=1
THERMAL_CHALLENGE: case=volumetric_power_sweep iterations=1 residual=1.16415e-10 backward_error=7.83555e-17 energy_imbalance=8.67838e-18
THERMAL_POWER_STUDY: qvol=100000 W/m3 Tmax=6550 dTmax=6250 exact_dTmax_discrete=6250 exact_dTmax_continuous=6250 expected_discretization_error=6.10352 observed_continuous_cell_error=6.10352 generated_power=100000 expected_power=100000 Tmin=690.625
  PASS: thermal_volumetric_power_sweep
  FAIL: thermal_cht_interface_conservation - CHT conductivity assertions disagree with energy materials
RADIATION_RESIDUAL: iteration=1 intensity_delta=146.2 intensity_relative=0.0795775 source_delta=1.13687e-13 source_relative=1.13687e-13
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
RADIATION_EQUILIBRIUM: iterations=2 G=1837.2 expected_G=1837.2 qrad=1.13687e-13 weight_sum=12.5664
  PASS: radiation_uniform_thermal_equilibrium
RADIATION_RESIDUAL: iteration=1 intensity_delta=146.2 intensity_relative=0.0795775 source_delta=1.13687e-13 source_relative=1.13687e-13
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=50 backward_error=0.0124386 linear_relative=1.11738e-16 dT_relative=0.0251904 energy_balance=0.0124386 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=25 backward_error=0.00618083 linear_relative=1.11738e-16 dT_relative=0.0124386 energy_balance=0.00618083 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=12.5 backward_error=0.0030809 linear_relative=1.11738e-16 dT_relative=0.00618083 energy_balance=0.0030809 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=6.25 backward_error=0.00153808 linear_relative=1.11738e-16 dT_relative=0.0030809 energy_balance=0.00153808 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=3.125 backward_error=0.000768448 linear_relative=1.11738e-16 dT_relative=0.00153808 energy_balance=0.000768448 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=1.5625 backward_error=0.000384077 linear_relative=1.11738e-16 dT_relative=0.000768448 energy_balance=0.000384077 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=0.78125 backward_error=0.000192001 linear_relative=1.11738e-16 dT_relative=0.000384077 energy_balance=0.000192001 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=0.390625 backward_error=9.59915e-05 linear_relative=1.11738e-16 dT_relative=0.000192001 energy_balance=9.59915e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=0.195312 backward_error=4.79935e-05 linear_relative=1.11738e-16 dT_relative=9.59915e-05 energy_balance=4.79935e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=0.0976562 backward_error=2.39961e-05 linear_relative=1.11738e-16 dT_relative=4.79935e-05 energy_balance=2.39961e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=0.0488281 backward_error=1.19979e-05 linear_relative=1.11738e-16 dT_relative=2.39961e-05 energy_balance=1.19979e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=0.0244141 backward_error=5.99893e-06 linear_relative=1.11738e-16 dT_relative=1.19979e-05 energy_balance=5.99893e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=0.012207 backward_error=2.99946e-06 linear_relative=1.11738e-16 dT_relative=5.99893e-06 energy_balance=2.99946e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=0.00610352 backward_error=1.49973e-06 linear_relative=1.11738e-16 dT_relative=2.99946e-06 energy_balance=1.49973e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=0.00305176 backward_error=7.49862e-07 linear_relative=1.11738e-16 dT_relative=1.49973e-06 energy_balance=7.49862e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=0.00152588 backward_error=3.74931e-07 linear_relative=1.11738e-16 dT_relative=7.49862e-07 energy_balance=3.74931e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=0.000762939 backward_error=1.87465e-07 linear_relative=1.11738e-16 dT_relative=3.74931e-07 energy_balance=1.87465e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=0.00038147 backward_error=9.37327e-08 linear_relative=1.11738e-16 dT_relative=1.87465e-07 energy_balance=9.37327e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=0.000190735 backward_error=4.68664e-08 linear_relative=1.11738e-16 dT_relative=9.37327e-08 energy_balance=4.68664e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=9.53674e-05 backward_error=2.34332e-08 linear_relative=1.11738e-16 dT_relative=4.68664e-08 energy_balance=2.34332e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=4.76837e-05 backward_error=1.17166e-08 linear_relative=1.11738e-16 dT_relative=2.34332e-08 energy_balance=1.17166e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=2.38419e-05 backward_error=5.85829e-09 linear_relative=1.11738e-16 dT_relative=1.17166e-08 energy_balance=5.85829e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=23 linear_status=0 residual=1.19209e-05 backward_error=2.92915e-09 linear_relative=1.11738e-16 dT_relative=5.85829e-09 energy_balance=2.92915e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=24 linear_status=0 residual=5.96046e-06 backward_error=1.46457e-09 linear_relative=1.11738e-16 dT_relative=2.92915e-09 energy_balance=1.46457e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=25 linear_status=0 residual=2.98023e-06 backward_error=7.32287e-10 linear_relative=1.11738e-16 dT_relative=1.46457e-09 energy_balance=7.32287e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=26 linear_status=0 residual=1.49012e-06 backward_error=3.66143e-10 linear_relative=1.11738e-16 dT_relative=7.32287e-10 energy_balance=3.66143e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=27 linear_status=0 residual=7.45058e-07 backward_error=1.83072e-10 linear_relative=1.11738e-16 dT_relative=3.66143e-10 energy_balance=1.83072e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=28 linear_status=0 residual=3.72529e-07 backward_error=9.15359e-11 linear_relative=1.11738e-16 dT_relative=1.83072e-10 energy_balance=9.15359e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=29 linear_status=0 residual=1.86265e-07 backward_error=4.57679e-11 linear_relative=1.11738e-16 dT_relative=9.1536e-11 energy_balance=4.57679e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=30 linear_status=0 residual=9.31323e-08 backward_error=2.2884e-11 linear_relative=1.11738e-16 dT_relative=4.57679e-11 energy_balance=2.2884e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=31 linear_status=0 residual=4.65664e-08 backward_error=1.1442e-11 linear_relative=1.11738e-16 dT_relative=2.28839e-11 energy_balance=1.1442e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=32 linear_status=0 residual=2.32833e-08 backward_error=5.72105e-12 linear_relative=1.11738e-16 dT_relative=1.1442e-11 energy_balance=5.72105e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=33 linear_status=0 residual=1.16418e-08 backward_error=2.86055e-12 linear_relative=1.11738e-16 dT_relative=5.72101e-12 energy_balance=2.86055e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=34 linear_status=0 residual=5.82099e-09 backward_error=1.4303e-12 linear_relative=1.11738e-16 dT_relative=2.8605e-12 energy_balance=1.4303e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=35 linear_status=0 residual=2.91038e-09 backward_error=7.15124e-13 linear_relative=1.11738e-16 dT_relative=1.43034e-12 energy_balance=7.15124e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=36 linear_status=0 residual=1.45519e-09 backward_error=3.57562e-13 linear_relative=1.11738e-16 dT_relative=7.15081e-13 energy_balance=3.57562e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=1 dT=7.75242 qrad_delta=1.13687e-13 qrad_relative=1.13687e-13 energy_balance=0.027258
RADIATION_RESIDUAL: iteration=1 intensity_delta=147.166 intensity_relative=0.0795863 source_delta=18.5453 source_relative=1
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=25.2225 backward_error=0.0060553 linear_relative=0 dT_relative=0.0121844 energy_balance=0.0060553 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=12.6113 backward_error=0.00301851 linear_relative=0 dT_relative=0.0060553 energy_balance=0.00301851 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=6.30563 backward_error=0.00150698 linear_relative=0 dT_relative=0.00301851 energy_balance=0.00150698 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=3.15281 backward_error=0.000752923 linear_relative=0 dT_relative=0.00150698 energy_balance=0.000752923 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=1.57641 backward_error=0.00037632 linear_relative=0 dT_relative=0.000752923 energy_balance=0.00037632 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.788203 backward_error=0.000188124 linear_relative=0 dT_relative=0.00037632 energy_balance=0.000188124 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=0.394102 backward_error=9.40534e-05 linear_relative=0 dT_relative=0.000188124 energy_balance=9.40534e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=0.197051 backward_error=4.70245e-05 linear_relative=0 dT_relative=9.40534e-05 energy_balance=4.70245e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=0.0985254 backward_error=2.35117e-05 linear_relative=0 dT_relative=4.70245e-05 energy_balance=2.35117e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=0.0492627 backward_error=1.17557e-05 linear_relative=0 dT_relative=2.35117e-05 energy_balance=1.17557e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=0.0246314 backward_error=5.87782e-06 linear_relative=0 dT_relative=1.17557e-05 energy_balance=5.87782e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=0.0123157 backward_error=2.9389e-06 linear_relative=0 dT_relative=5.87782e-06 energy_balance=2.9389e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=0.00615784 backward_error=1.46945e-06 linear_relative=0 dT_relative=2.9389e-06 energy_balance=1.46945e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=0.00307892 backward_error=7.34724e-07 linear_relative=0 dT_relative=1.46945e-06 energy_balance=7.34724e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=0.00153946 backward_error=3.67362e-07 linear_relative=0 dT_relative=7.34724e-07 energy_balance=3.67362e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=0.00076973 backward_error=1.83681e-07 linear_relative=0 dT_relative=3.67362e-07 energy_balance=1.83681e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=0.000384865 backward_error=9.18404e-08 linear_relative=0 dT_relative=1.83681e-07 energy_balance=9.18404e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=0.000192432 backward_error=4.59202e-08 linear_relative=0 dT_relative=9.18404e-08 energy_balance=4.59202e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=9.62162e-05 backward_error=2.29601e-08 linear_relative=0 dT_relative=4.59202e-08 energy_balance=2.29601e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=4.81081e-05 backward_error=1.148e-08 linear_relative=0 dT_relative=2.29601e-08 energy_balance=1.148e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=2.40541e-05 backward_error=5.74002e-09 linear_relative=0 dT_relative=1.148e-08 energy_balance=5.74002e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=1.2027e-05 backward_error=2.87001e-09 linear_relative=0 dT_relative=5.74002e-09 energy_balance=2.87001e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=23 linear_status=0 residual=6.01351e-06 backward_error=1.43501e-09 linear_relative=0 dT_relative=2.87001e-09 energy_balance=1.43501e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=24 linear_status=0 residual=3.00676e-06 backward_error=7.17503e-10 linear_relative=0 dT_relative=1.43501e-09 energy_balance=7.17503e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=25 linear_status=0 residual=1.50338e-06 backward_error=3.58751e-10 linear_relative=0 dT_relative=7.17503e-10 energy_balance=3.58751e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=26 linear_status=0 residual=7.51689e-07 backward_error=1.79376e-10 linear_relative=0 dT_relative=3.58751e-10 energy_balance=1.79376e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=27 linear_status=0 residual=3.75844e-07 backward_error=8.96878e-11 linear_relative=0 dT_relative=1.79376e-10 energy_balance=8.96878e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=28 linear_status=0 residual=1.87922e-07 backward_error=4.48439e-11 linear_relative=0 dT_relative=8.96878e-11 energy_balance=4.48439e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=29 linear_status=0 residual=9.39608e-08 backward_error=2.24219e-11 linear_relative=0 dT_relative=4.48439e-11 energy_balance=2.24219e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=30 linear_status=0 residual=4.69804e-08 backward_error=1.12109e-11 linear_relative=0 dT_relative=2.2422e-11 energy_balance=1.12109e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=31 linear_status=0 residual=2.349e-08 backward_error=5.60542e-12 linear_relative=0 dT_relative=1.1211e-11 energy_balance=5.60542e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=32 linear_status=0 residual=1.17452e-08 backward_error=2.80276e-12 linear_relative=0 dT_relative=5.60549e-12 energy_balance=2.80276e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=33 linear_status=0 residual=5.87215e-09 backward_error=1.40127e-12 linear_relative=0 dT_relative=2.80283e-12 energy_balance=1.40127e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=34 linear_status=0 residual=2.9363e-09 backward_error=7.00691e-13 linear_relative=0 dT_relative=1.40133e-12 energy_balance=7.00691e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=35 linear_status=0 residual=1.46792e-09 backward_error=3.50291e-13 linear_relative=0 dT_relative=7.00664e-13 energy_balance=3.50291e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=2 dT=3.79602 qrad_delta=18.5453 qrad_relative=1 energy_balance=0.0139501
RADIATION_RESIDUAL: iteration=1 intensity_delta=147.667 intensity_relative=0.0795909 source_delta=9.6072 source_relative=0.341255
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=12.8269 backward_error=0.00302495 linear_relative=0 dT_relative=0.00606826 energy_balance=0.00302495 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=6.41343 backward_error=0.00151019 linear_relative=0 dT_relative=0.00302495 energy_balance=0.00151019 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=3.20671 backward_error=0.000754527 linear_relative=0 dT_relative=0.00151019 energy_balance=0.000754527 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=1.60336 backward_error=0.000377121 linear_relative=0 dT_relative=0.000754527 energy_balance=0.000377121 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.801679 backward_error=0.000188525 linear_relative=0 dT_relative=0.000377121 energy_balance=0.000188525 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.400839 backward_error=9.42536e-05 linear_relative=0 dT_relative=0.000188525 energy_balance=9.42536e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=0.20042 backward_error=4.71246e-05 linear_relative=0 dT_relative=9.42536e-05 energy_balance=4.71246e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=0.10021 backward_error=2.35617e-05 linear_relative=0 dT_relative=4.71246e-05 energy_balance=2.35617e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=0.0501049 backward_error=1.17807e-05 linear_relative=0 dT_relative=2.35617e-05 energy_balance=1.17807e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=0.0250525 backward_error=5.89033e-06 linear_relative=0 dT_relative=1.17807e-05 energy_balance=5.89033e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=0.0125262 backward_error=2.94516e-06 linear_relative=0 dT_relative=5.89033e-06 energy_balance=2.94516e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=0.00626311 backward_error=1.47258e-06 linear_relative=0 dT_relative=2.94516e-06 energy_balance=1.47258e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=0.00313156 backward_error=7.36287e-07 linear_relative=0 dT_relative=1.47258e-06 energy_balance=7.36287e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=0.00156578 backward_error=3.68144e-07 linear_relative=0 dT_relative=7.36287e-07 energy_balance=3.68144e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=0.000782889 backward_error=1.84072e-07 linear_relative=0 dT_relative=3.68144e-07 energy_balance=1.84072e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=0.000391445 backward_error=9.20359e-08 linear_relative=0 dT_relative=1.84072e-07 energy_balance=9.20359e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=0.000195722 backward_error=4.60179e-08 linear_relative=0 dT_relative=9.20359e-08 energy_balance=4.60179e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=9.78612e-05 backward_error=2.3009e-08 linear_relative=0 dT_relative=4.60179e-08 energy_balance=2.3009e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=4.89306e-05 backward_error=1.15045e-08 linear_relative=0 dT_relative=2.3009e-08 energy_balance=1.15045e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=2.44653e-05 backward_error=5.75224e-09 linear_relative=0 dT_relative=1.15045e-08 energy_balance=5.75224e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=1.22326e-05 backward_error=2.87612e-09 linear_relative=0 dT_relative=5.75224e-09 energy_balance=2.87612e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=6.11632e-06 backward_error=1.43806e-09 linear_relative=0 dT_relative=2.87612e-09 energy_balance=1.43806e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=23 linear_status=0 residual=3.05816e-06 backward_error=7.1903e-10 linear_relative=0 dT_relative=1.43806e-09 energy_balance=7.1903e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=24 linear_status=0 residual=1.52908e-06 backward_error=3.59515e-10 linear_relative=0 dT_relative=7.1903e-10 energy_balance=3.59515e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=25 linear_status=0 residual=7.6454e-07 backward_error=1.79758e-10 linear_relative=0 dT_relative=3.59515e-10 energy_balance=1.79758e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=26 linear_status=0 residual=3.8227e-07 backward_error=8.98787e-11 linear_relative=0 dT_relative=1.79758e-10 energy_balance=8.98787e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=27 linear_status=0 residual=1.91135e-07 backward_error=4.49393e-11 linear_relative=0 dT_relative=8.98787e-11 energy_balance=4.49393e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=28 linear_status=0 residual=9.55674e-08 backward_error=2.24697e-11 linear_relative=0 dT_relative=4.49394e-11 energy_balance=2.24697e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=29 linear_status=0 residual=4.77839e-08 backward_error=1.12349e-11 linear_relative=0 dT_relative=2.24697e-11 energy_balance=1.12349e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=30 linear_status=0 residual=2.3892e-08 backward_error=5.61744e-12 linear_relative=0 dT_relative=1.12348e-11 energy_balance=5.61744e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=31 linear_status=0 residual=1.19458e-08 backward_error=2.80867e-12 linear_relative=0 dT_relative=5.61742e-12 energy_balance=2.80867e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=32 linear_status=0 residual=5.97311e-09 backward_error=1.40439e-12 linear_relative=0 dT_relative=2.80871e-12 energy_balance=1.40439e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=33 linear_status=0 residual=2.98633e-09 backward_error=7.0214e-13 linear_relative=0 dT_relative=1.40436e-12 energy_balance=7.0214e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=34 linear_status=0 residual=1.49339e-09 backward_error=3.51124e-13 linear_relative=0 dT_relative=7.02087e-13 energy_balance=3.51124e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=3 dT=1.9021 qrad_delta=9.6072 qrad_relative=0.341255 energy_balance=0.00714482
RADIATION_RESIDUAL: iteration=1 intensity_delta=147.925 intensity_relative=0.0795932 source_delta=4.94791 source_relative=0.149482
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=6.5487 backward_error=0.00153048 linear_relative=0 dT_relative=0.00306565 energy_balance=0.00153048 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=3.27435 backward_error=0.000764653 linear_relative=0 dT_relative=0.00153048 energy_balance=0.000764653 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=1.63717 backward_error=0.000382181 linear_relative=0 dT_relative=0.000764653 energy_balance=0.000382181 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.818587 backward_error=0.000191054 linear_relative=0 dT_relative=0.000382181 energy_balance=0.000191054 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.409294 backward_error=9.55178e-05 linear_relative=0 dT_relative=0.000191054 energy_balance=9.55178e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.204647 backward_error=4.77566e-05 linear_relative=0 dT_relative=9.55178e-05 energy_balance=4.77566e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=0.102323 backward_error=2.38777e-05 linear_relative=0 dT_relative=4.77566e-05 energy_balance=2.38777e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=0.0511617 backward_error=1.19387e-05 linear_relative=0 dT_relative=2.38777e-05 energy_balance=1.19387e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=0.0255809 backward_error=5.96933e-06 linear_relative=0 dT_relative=1.19387e-05 energy_balance=5.96933e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=0.0127904 backward_error=2.98465e-06 linear_relative=0 dT_relative=5.96933e-06 energy_balance=2.98465e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=0.00639521 backward_error=1.49232e-06 linear_relative=0 dT_relative=2.98465e-06 energy_balance=1.49232e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=0.00319761 backward_error=7.46162e-07 linear_relative=0 dT_relative=1.49232e-06 energy_balance=7.46162e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=0.0015988 backward_error=3.73081e-07 linear_relative=0 dT_relative=7.46162e-07 energy_balance=3.73081e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=0.000799402 backward_error=1.8654e-07 linear_relative=0 dT_relative=3.73081e-07 energy_balance=1.8654e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=0.000399701 backward_error=9.32702e-08 linear_relative=0 dT_relative=1.8654e-07 energy_balance=9.32702e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=0.00019985 backward_error=4.66351e-08 linear_relative=0 dT_relative=9.32702e-08 energy_balance=4.66351e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=9.99252e-05 backward_error=2.33175e-08 linear_relative=0 dT_relative=4.66351e-08 energy_balance=2.33175e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=4.99626e-05 backward_error=1.16588e-08 linear_relative=0 dT_relative=2.33175e-08 energy_balance=1.16588e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=2.49813e-05 backward_error=5.82938e-09 linear_relative=0 dT_relative=1.16588e-08 energy_balance=5.82938e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=1.24907e-05 backward_error=2.91469e-09 linear_relative=0 dT_relative=5.82938e-09 energy_balance=2.91469e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=6.24533e-06 backward_error=1.45735e-09 linear_relative=0 dT_relative=2.91469e-09 energy_balance=1.45735e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=3.12266e-06 backward_error=7.28673e-10 linear_relative=0 dT_relative=1.45735e-09 energy_balance=7.28673e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=23 linear_status=0 residual=1.56133e-06 backward_error=3.64336e-10 linear_relative=0 dT_relative=7.28673e-10 energy_balance=3.64336e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=24 linear_status=0 residual=7.80666e-07 backward_error=1.82168e-10 linear_relative=0 dT_relative=3.64337e-10 energy_balance=1.82168e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=25 linear_status=0 residual=3.90333e-07 backward_error=9.10841e-11 linear_relative=0 dT_relative=1.82168e-10 energy_balance=9.10841e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=26 linear_status=0 residual=1.95167e-07 backward_error=4.55421e-11 linear_relative=0 dT_relative=9.10841e-11 energy_balance=4.55421e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=27 linear_status=0 residual=9.75833e-08 backward_error=2.27711e-11 linear_relative=0 dT_relative=4.5542e-11 energy_balance=2.27711e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=28 linear_status=0 residual=4.87921e-08 backward_error=1.13856e-11 linear_relative=0 dT_relative=2.2771e-11 energy_balance=1.13856e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=29 linear_status=0 residual=2.43958e-08 backward_error=5.69276e-12 linear_relative=0 dT_relative=1.13857e-11 energy_balance=5.69276e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=30 linear_status=0 residual=1.21981e-08 backward_error=2.84643e-12 linear_relative=0 dT_relative=5.69266e-12 energy_balance=2.84643e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=31 linear_status=0 residual=6.09907e-09 backward_error=1.42322e-12 linear_relative=0 dT_relative=2.84633e-12 energy_balance=1.42322e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=32 linear_status=0 residual=3.04999e-09 backward_error=7.11715e-13 linear_relative=0 dT_relative=1.42317e-12 energy_balance=7.11715e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=33 linear_status=0 residual=1.52522e-09 backward_error=3.5591e-13 linear_relative=0 dT_relative=7.11583e-13 energy_balance=3.5591e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=4 dT=0.963883 qrad_delta=4.94791 qrad_relative=0.149482 energy_balance=0.00366079
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.057 intensity_relative=0.0795944 source_delta=2.54197 source_relative=0.0713188
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=3.34995 backward_error=0.000779307 linear_relative=0 dT_relative=0.00155983 energy_balance=0.000779307 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=1.67497 backward_error=0.000389502 linear_relative=0 dT_relative=0.000779307 energy_balance=0.000389502 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.837487 backward_error=0.000194713 linear_relative=0 dT_relative=0.000389502 energy_balance=0.000194713 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.418743 backward_error=9.7347e-05 linear_relative=0 dT_relative=0.000194713 energy_balance=9.7347e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.209372 backward_error=4.86711e-05 linear_relative=0 dT_relative=9.7347e-05 energy_balance=4.86711e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.104686 backward_error=2.4335e-05 linear_relative=0 dT_relative=4.86711e-05 energy_balance=2.4335e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=0.0523429 backward_error=1.21673e-05 linear_relative=0 dT_relative=2.4335e-05 energy_balance=1.21673e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=0.0261715 backward_error=6.08363e-06 linear_relative=0 dT_relative=1.21673e-05 energy_balance=6.08363e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=0.0130857 backward_error=3.04181e-06 linear_relative=0 dT_relative=6.08363e-06 energy_balance=3.04181e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=0.00654286 backward_error=1.5209e-06 linear_relative=0 dT_relative=3.04181e-06 energy_balance=1.5209e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=0.00327143 backward_error=7.6045e-07 linear_relative=0 dT_relative=1.5209e-06 energy_balance=7.6045e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=0.00163572 backward_error=3.80225e-07 linear_relative=0 dT_relative=7.6045e-07 energy_balance=3.80225e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=0.000817858 backward_error=1.90112e-07 linear_relative=0 dT_relative=3.80225e-07 energy_balance=1.90112e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=0.000408929 backward_error=9.50562e-08 linear_relative=0 dT_relative=1.90112e-07 energy_balance=9.50562e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=0.000204465 backward_error=4.75281e-08 linear_relative=0 dT_relative=9.50562e-08 energy_balance=4.75281e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=0.000102232 backward_error=2.3764e-08 linear_relative=0 dT_relative=4.75281e-08 energy_balance=2.3764e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=5.11161e-05 backward_error=1.1882e-08 linear_relative=0 dT_relative=2.3764e-08 energy_balance=1.1882e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=2.55581e-05 backward_error=5.94101e-09 linear_relative=0 dT_relative=1.1882e-08 energy_balance=5.94101e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=1.2779e-05 backward_error=2.97051e-09 linear_relative=0 dT_relative=5.94101e-09 energy_balance=2.97051e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=6.38952e-06 backward_error=1.48525e-09 linear_relative=0 dT_relative=2.97051e-09 energy_balance=1.48525e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=3.19476e-06 backward_error=7.42627e-10 linear_relative=0 dT_relative=1.48525e-09 energy_balance=7.42627e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=1.59738e-06 backward_error=3.71313e-10 linear_relative=0 dT_relative=7.42626e-10 energy_balance=3.71313e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=23 linear_status=0 residual=7.9869e-07 backward_error=1.85657e-10 linear_relative=0 dT_relative=3.71313e-10 energy_balance=1.85657e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=24 linear_status=0 residual=3.99345e-07 backward_error=9.28284e-11 linear_relative=0 dT_relative=1.85657e-10 energy_balance=9.28284e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=25 linear_status=0 residual=1.99673e-07 backward_error=4.64142e-11 linear_relative=0 dT_relative=9.28284e-11 energy_balance=4.64142e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=26 linear_status=0 residual=9.98366e-08 backward_error=2.32072e-11 linear_relative=0 dT_relative=4.64142e-11 energy_balance=2.32072e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=27 linear_status=0 residual=4.99181e-08 backward_error=1.16035e-11 linear_relative=0 dT_relative=2.32071e-11 energy_balance=1.16035e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=28 linear_status=0 residual=2.49588e-08 backward_error=5.80171e-12 linear_relative=0 dT_relative=1.16036e-11 energy_balance=5.80171e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=29 linear_status=0 residual=1.24796e-08 backward_error=2.90091e-12 linear_relative=0 dT_relative=5.80173e-12 energy_balance=2.90091e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=30 linear_status=0 residual=6.24004e-09 backward_error=1.45051e-12 linear_relative=0 dT_relative=2.90077e-12 energy_balance=1.45051e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=31 linear_status=0 residual=3.12002e-09 backward_error=7.25254e-13 linear_relative=0 dT_relative=1.45048e-12 energy_balance=7.25254e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=32 linear_status=0 residual=1.56024e-09 backward_error=3.6268e-13 linear_relative=0 dT_relative=7.25238e-13 energy_balance=3.6268e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=5 dT=0.491199 qrad_delta=2.54197 qrad_relative=0.0713188 energy_balance=0.00187604
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.125 intensity_relative=0.079595 source_delta=1.30443 source_relative=0.0353055
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=1.71534 backward_error=0.000398105 linear_relative=0 dT_relative=0.000796527 energy_balance=0.000398105 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.857668 backward_error=0.000199013 linear_relative=0 dT_relative=0.000398105 energy_balance=0.000199013 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.428834 backward_error=9.94965e-05 linear_relative=0 dT_relative=0.000199013 energy_balance=9.94965e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.214417 backward_error=4.97458e-05 linear_relative=0 dT_relative=9.94965e-05 energy_balance=4.97458e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.107208 backward_error=2.48723e-05 linear_relative=0 dT_relative=4.97458e-05 energy_balance=2.48723e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.0536042 backward_error=1.2436e-05 linear_relative=0 dT_relative=2.48723e-05 energy_balance=1.2436e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=0.0268021 backward_error=6.21795e-06 linear_relative=0 dT_relative=1.2436e-05 energy_balance=6.21795e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=0.0134011 backward_error=3.10897e-06 linear_relative=0 dT_relative=6.21795e-06 energy_balance=3.10897e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=0.00670053 backward_error=1.55448e-06 linear_relative=0 dT_relative=3.10897e-06 energy_balance=1.55448e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=0.00335026 backward_error=7.77239e-07 linear_relative=0 dT_relative=1.55448e-06 energy_balance=7.77239e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=0.00167513 backward_error=3.8862e-07 linear_relative=0 dT_relative=7.77239e-07 energy_balance=3.8862e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=0.000837566 backward_error=1.9431e-07 linear_relative=0 dT_relative=3.8862e-07 energy_balance=1.9431e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=0.000418783 backward_error=9.71549e-08 linear_relative=0 dT_relative=1.9431e-07 energy_balance=9.71549e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=0.000209392 backward_error=4.85774e-08 linear_relative=0 dT_relative=9.71549e-08 energy_balance=4.85774e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=0.000104696 backward_error=2.42887e-08 linear_relative=0 dT_relative=4.85774e-08 energy_balance=2.42887e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=5.23479e-05 backward_error=1.21444e-08 linear_relative=0 dT_relative=2.42887e-08 energy_balance=1.21444e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=2.61739e-05 backward_error=6.07218e-09 linear_relative=0 dT_relative=1.21444e-08 energy_balance=6.07218e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=1.3087e-05 backward_error=3.03609e-09 linear_relative=0 dT_relative=6.07218e-09 energy_balance=3.03609e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=6.54349e-06 backward_error=1.51804e-09 linear_relative=0 dT_relative=3.03609e-09 energy_balance=1.51804e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=3.27174e-06 backward_error=7.59022e-10 linear_relative=0 dT_relative=1.51804e-09 energy_balance=7.59022e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=1.63587e-06 backward_error=3.79511e-10 linear_relative=0 dT_relative=7.59022e-10 energy_balance=3.79511e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=8.17936e-07 backward_error=1.89756e-10 linear_relative=0 dT_relative=3.79511e-10 energy_balance=1.89756e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=23 linear_status=0 residual=4.08968e-07 backward_error=9.48778e-11 linear_relative=0 dT_relative=1.89756e-10 energy_balance=9.48778e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=24 linear_status=0 residual=2.04484e-07 backward_error=4.7439e-11 linear_relative=0 dT_relative=9.48777e-11 energy_balance=4.7439e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=25 linear_status=0 residual=1.02242e-07 backward_error=2.37195e-11 linear_relative=0 dT_relative=4.74389e-11 energy_balance=2.37195e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=26 linear_status=0 residual=5.11209e-08 backward_error=1.18597e-11 linear_relative=0 dT_relative=2.37196e-11 energy_balance=1.18597e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=27 linear_status=0 residual=2.55604e-08 backward_error=5.92985e-12 linear_relative=0 dT_relative=1.18597e-11 energy_balance=5.92985e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=28 linear_status=0 residual=1.27802e-08 backward_error=2.96492e-12 linear_relative=0 dT_relative=5.92985e-12 energy_balance=2.96492e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=29 linear_status=0 residual=6.39056e-09 backward_error=1.48257e-12 linear_relative=0 dT_relative=2.96483e-12 energy_balance=1.48257e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=30 linear_status=0 residual=3.19551e-09 backward_error=7.41337e-13 linear_relative=0 dT_relative=1.48251e-12 energy_balance=7.41337e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=31 linear_status=0 residual=1.59753e-09 backward_error=3.70616e-13 linear_relative=0 dT_relative=7.41253e-13 energy_balance=3.70616e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=6 dT=0.251031 qrad_delta=1.30443 qrad_relative=0.0353055 energy_balance=0.000961512
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.16 intensity_relative=0.0795953 source_delta=0.668997 source_relative=0.017785
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.878776 backward_error=0.000203706 linear_relative=0 dT_relative=0.000407494 energy_balance=0.000203706 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.439388 backward_error=0.000101842 linear_relative=0 dT_relative=0.000203706 energy_balance=0.000101842 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.219694 backward_error=5.09186e-05 linear_relative=0 dT_relative=0.000101842 energy_balance=5.09186e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.109847 backward_error=2.54587e-05 linear_relative=0 dT_relative=5.09186e-05 energy_balance=2.54587e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.0549235 backward_error=1.27292e-05 linear_relative=0 dT_relative=2.54587e-05 energy_balance=1.27292e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.0274617 backward_error=6.36455e-06 linear_relative=0 dT_relative=1.27292e-05 energy_balance=6.36455e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=0.0137309 backward_error=3.18226e-06 linear_relative=0 dT_relative=6.36455e-06 energy_balance=3.18226e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=0.00686544 backward_error=1.59113e-06 linear_relative=0 dT_relative=3.18226e-06 energy_balance=1.59113e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=0.00343272 backward_error=7.95564e-07 linear_relative=0 dT_relative=1.59113e-06 energy_balance=7.95564e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=0.00171636 backward_error=3.97782e-07 linear_relative=0 dT_relative=7.95564e-07 energy_balance=3.97782e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=0.00085818 backward_error=1.98891e-07 linear_relative=0 dT_relative=3.97782e-07 energy_balance=1.98891e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=0.00042909 backward_error=9.94454e-08 linear_relative=0 dT_relative=1.98891e-07 energy_balance=9.94454e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=0.000214545 backward_error=4.97227e-08 linear_relative=0 dT_relative=9.94454e-08 energy_balance=4.97227e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=0.000107272 backward_error=2.48613e-08 linear_relative=0 dT_relative=4.97227e-08 energy_balance=2.48613e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=5.36362e-05 backward_error=1.24307e-08 linear_relative=0 dT_relative=2.48614e-08 energy_balance=1.24307e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=2.68181e-05 backward_error=6.21534e-09 linear_relative=0 dT_relative=1.24307e-08 energy_balance=6.21534e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=1.34091e-05 backward_error=3.10767e-09 linear_relative=0 dT_relative=6.21534e-09 energy_balance=3.10767e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=6.70453e-06 backward_error=1.55383e-09 linear_relative=0 dT_relative=3.10767e-09 energy_balance=1.55383e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=3.35226e-06 backward_error=7.76917e-10 linear_relative=0 dT_relative=1.55383e-09 energy_balance=7.76917e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=1.67613e-06 backward_error=3.88459e-10 linear_relative=0 dT_relative=7.76917e-10 energy_balance=3.88459e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=8.38066e-07 backward_error=1.94229e-10 linear_relative=0 dT_relative=3.88459e-10 energy_balance=1.94229e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=4.19033e-07 backward_error=9.71146e-11 linear_relative=0 dT_relative=1.94229e-10 energy_balance=9.71146e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=23 linear_status=0 residual=2.09516e-07 backward_error=4.85572e-11 linear_relative=0 dT_relative=9.71146e-11 energy_balance=4.85572e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=24 linear_status=0 residual=1.04758e-07 backward_error=2.42787e-11 linear_relative=0 dT_relative=4.85572e-11 energy_balance=2.42787e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=25 linear_status=0 residual=5.23792e-08 backward_error=1.21393e-11 linear_relative=0 dT_relative=2.42787e-11 energy_balance=1.21393e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=26 linear_status=0 residual=2.61894e-08 backward_error=6.06962e-12 linear_relative=0 dT_relative=1.21393e-11 energy_balance=6.06962e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=27 linear_status=0 residual=1.30949e-08 backward_error=3.03486e-12 linear_relative=0 dT_relative=6.06958e-12 energy_balance=3.03486e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=28 linear_status=0 residual=6.54745e-09 backward_error=1.51743e-12 linear_relative=0 dT_relative=3.03488e-12 energy_balance=1.51743e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=29 linear_status=0 residual=3.27373e-09 backward_error=7.58715e-13 linear_relative=0 dT_relative=1.51744e-12 energy_balance=7.58715e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=30 linear_status=0 residual=1.63664e-09 backward_error=3.79305e-13 linear_relative=0 dT_relative=7.58721e-13 energy_balance=3.79305e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=7 dT=0.128477 qrad_delta=0.668997 qrad_relative=0.017785 energy_balance=0.00049282
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.178 intensity_relative=0.0795955 source_delta=0.34301 source_relative=0.00903636
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.450317 backward_error=0.000104322 linear_relative=0 dT_relative=0.000208666 energy_balance=0.000104322 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.225159 backward_error=5.21582e-05 linear_relative=0 dT_relative=0.000104322 energy_balance=5.21582e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.112579 backward_error=2.60784e-05 linear_relative=0 dT_relative=5.21582e-05 energy_balance=2.60784e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.0562897 backward_error=1.3039e-05 linear_relative=0 dT_relative=2.60784e-05 energy_balance=1.3039e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.0281448 backward_error=6.51948e-06 linear_relative=0 dT_relative=1.3039e-05 energy_balance=6.51948e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.0140724 backward_error=3.25973e-06 linear_relative=0 dT_relative=6.51948e-06 energy_balance=3.25973e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=0.00703621 backward_error=1.62986e-06 linear_relative=0 dT_relative=3.25973e-06 energy_balance=1.62986e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=0.0035181 backward_error=8.1493e-07 linear_relative=0 dT_relative=1.62986e-06 energy_balance=8.1493e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=0.00175905 backward_error=4.07465e-07 linear_relative=0 dT_relative=8.1493e-07 energy_balance=4.07465e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=0.000879526 backward_error=2.03732e-07 linear_relative=0 dT_relative=4.07465e-07 energy_balance=2.03732e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=0.000439763 backward_error=1.01866e-07 linear_relative=0 dT_relative=2.03732e-07 energy_balance=1.01866e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=0.000219881 backward_error=5.09331e-08 linear_relative=0 dT_relative=1.01866e-07 energy_balance=5.09331e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=0.000109941 backward_error=2.54666e-08 linear_relative=0 dT_relative=5.09331e-08 energy_balance=2.54666e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=5.49704e-05 backward_error=1.27333e-08 linear_relative=0 dT_relative=2.54666e-08 energy_balance=1.27333e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=2.74852e-05 backward_error=6.36664e-09 linear_relative=0 dT_relative=1.27333e-08 energy_balance=6.36664e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=1.37426e-05 backward_error=3.18332e-09 linear_relative=0 dT_relative=6.36664e-09 energy_balance=3.18332e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=6.8713e-06 backward_error=1.59166e-09 linear_relative=0 dT_relative=3.18332e-09 energy_balance=1.59166e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=3.43565e-06 backward_error=7.9583e-10 linear_relative=0 dT_relative=1.59166e-09 energy_balance=7.9583e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=1.71782e-06 backward_error=3.97915e-10 linear_relative=0 dT_relative=7.9583e-10 energy_balance=3.97915e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=8.58912e-07 backward_error=1.98958e-10 linear_relative=0 dT_relative=3.97915e-10 energy_balance=1.98958e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=4.29456e-07 backward_error=9.94788e-11 linear_relative=0 dT_relative=1.98958e-10 energy_balance=9.94788e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=2.14728e-07 backward_error=4.97394e-11 linear_relative=0 dT_relative=9.94786e-11 energy_balance=4.97394e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=23 linear_status=0 residual=1.07364e-07 backward_error=2.48697e-11 linear_relative=0 dT_relative=4.97395e-11 energy_balance=2.48697e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=24 linear_status=0 residual=5.3682e-08 backward_error=1.24348e-11 linear_relative=0 dT_relative=2.48696e-11 energy_balance=1.24348e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=25 linear_status=0 residual=2.6841e-08 backward_error=6.21742e-12 linear_relative=0 dT_relative=1.2435e-11 energy_balance=6.21742e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=26 linear_status=0 residual=1.342e-08 backward_error=3.10861e-12 linear_relative=0 dT_relative=6.21748e-12 energy_balance=3.10861e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=27 linear_status=0 residual=6.7098e-09 backward_error=1.55425e-12 linear_relative=0 dT_relative=3.10874e-12 energy_balance=1.55425e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=28 linear_status=0 residual=3.35513e-09 backward_error=7.77178e-13 linear_relative=0 dT_relative=1.55419e-12 energy_balance=7.77178e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=29 linear_status=0 residual=1.67756e-09 backward_error=3.88589e-13 linear_relative=0 dT_relative=7.77095e-13 energy_balance=3.88589e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=8 dT=0.0658029 qrad_delta=0.34301 qrad_relative=0.00903636 energy_balance=0.0002526
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.187 intensity_relative=0.0795956 source_delta=0.175844 source_relative=0.00461113
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.230789 backward_error=5.34485e-05 linear_relative=0 dT_relative=0.000106903 energy_balance=5.34485e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.115395 backward_error=2.67235e-05 linear_relative=0 dT_relative=5.34485e-05 energy_balance=2.67235e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.0576973 backward_error=1.33616e-05 linear_relative=0 dT_relative=2.67235e-05 energy_balance=1.33616e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.0288487 backward_error=6.68075e-06 linear_relative=0 dT_relative=1.33616e-05 energy_balance=6.68075e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.0144243 backward_error=3.34036e-06 linear_relative=0 dT_relative=6.68075e-06 energy_balance=3.34036e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.00721217 backward_error=1.67018e-06 linear_relative=0 dT_relative=3.34036e-06 energy_balance=1.67018e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=0.00360608 backward_error=8.35089e-07 linear_relative=0 dT_relative=1.67018e-06 energy_balance=8.35089e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=0.00180304 backward_error=4.17544e-07 linear_relative=0 dT_relative=8.35089e-07 energy_balance=4.17544e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=0.000901521 backward_error=2.08772e-07 linear_relative=0 dT_relative=4.17544e-07 energy_balance=2.08772e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=0.000450761 backward_error=1.04386e-07 linear_relative=0 dT_relative=2.08772e-07 energy_balance=1.04386e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=0.00022538 backward_error=5.2193e-08 linear_relative=0 dT_relative=1.04386e-07 energy_balance=5.2193e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=0.00011269 backward_error=2.60965e-08 linear_relative=0 dT_relative=5.2193e-08 energy_balance=2.60965e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=5.63451e-05 backward_error=1.30482e-08 linear_relative=0 dT_relative=2.60965e-08 energy_balance=1.30482e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=2.81725e-05 backward_error=6.52412e-09 linear_relative=0 dT_relative=1.30482e-08 energy_balance=6.52412e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=1.40863e-05 backward_error=3.26206e-09 linear_relative=0 dT_relative=6.52412e-09 energy_balance=3.26206e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=7.04313e-06 backward_error=1.63103e-09 linear_relative=0 dT_relative=3.26206e-09 energy_balance=1.63103e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=3.52157e-06 backward_error=8.15516e-10 linear_relative=0 dT_relative=1.63103e-09 energy_balance=8.15516e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=1.76078e-06 backward_error=4.07758e-10 linear_relative=0 dT_relative=8.15516e-10 energy_balance=4.07758e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=8.80391e-07 backward_error=2.03879e-10 linear_relative=0 dT_relative=4.07758e-10 energy_balance=2.03879e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=4.40196e-07 backward_error=1.01939e-10 linear_relative=0 dT_relative=2.03879e-10 energy_balance=1.01939e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=2.20098e-07 backward_error=5.09698e-11 linear_relative=0 dT_relative=1.01939e-10 energy_balance=5.09698e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=1.10049e-07 backward_error=2.54848e-11 linear_relative=0 dT_relative=5.09698e-11 energy_balance=2.54848e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=23 linear_status=0 residual=5.50244e-08 backward_error=1.27424e-11 linear_relative=0 dT_relative=2.54849e-11 energy_balance=1.27424e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=24 linear_status=0 residual=2.75122e-08 backward_error=6.37121e-12 linear_relative=0 dT_relative=1.27425e-11 energy_balance=6.37121e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=25 linear_status=0 residual=1.37561e-08 backward_error=3.1856e-12 linear_relative=0 dT_relative=6.37117e-12 energy_balance=3.1856e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=26 linear_status=0 residual=6.87805e-09 backward_error=1.5928e-12 linear_relative=0 dT_relative=3.18568e-12 energy_balance=1.5928e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=27 linear_status=0 residual=3.4388e-09 backward_error=7.96349e-13 linear_relative=0 dT_relative=1.59275e-12 energy_balance=7.96349e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=28 linear_status=0 residual=1.7194e-09 backward_error=3.98174e-13 linear_relative=0 dT_relative=7.96374e-13 energy_balance=3.98174e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=9 dT=0.0337155 qrad_delta=0.175844 qrad_relative=0.00461113 energy_balance=0.000129475
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.192 intensity_relative=0.0795956 source_delta=0.0901401 source_relative=0.00235815
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.118288 backward_error=2.739e-05 linear_relative=0 dT_relative=5.47814e-05 energy_balance=2.739e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.0591442 backward_error=1.36948e-05 linear_relative=0 dT_relative=2.739e-05 energy_balance=1.36948e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.0295721 backward_error=6.84735e-06 linear_relative=0 dT_relative=1.36948e-05 energy_balance=6.84735e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.0147861 backward_error=3.42366e-06 linear_relative=0 dT_relative=6.84735e-06 energy_balance=3.42366e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.00739303 backward_error=1.71183e-06 linear_relative=0 dT_relative=3.42366e-06 energy_balance=1.71183e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.00369651 backward_error=8.55914e-07 linear_relative=0 dT_relative=1.71183e-06 energy_balance=8.55914e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=0.00184826 backward_error=4.27957e-07 linear_relative=0 dT_relative=8.55914e-07 energy_balance=4.27957e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=0.000924128 backward_error=2.13978e-07 linear_relative=0 dT_relative=4.27957e-07 energy_balance=2.13978e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=0.000462064 backward_error=1.06989e-07 linear_relative=0 dT_relative=2.13978e-07 energy_balance=1.06989e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=0.000231032 backward_error=5.34946e-08 linear_relative=0 dT_relative=1.06989e-07 energy_balance=5.34946e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=0.000115516 backward_error=2.67473e-08 linear_relative=0 dT_relative=5.34946e-08 energy_balance=2.67473e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=5.7758e-05 backward_error=1.33736e-08 linear_relative=0 dT_relative=2.67473e-08 energy_balance=1.33736e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=2.8879e-05 backward_error=6.68682e-09 linear_relative=0 dT_relative=1.33736e-08 energy_balance=6.68682e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=1.44395e-05 backward_error=3.34341e-09 linear_relative=0 dT_relative=6.68682e-09 energy_balance=3.34341e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=7.21975e-06 backward_error=1.6717e-09 linear_relative=0 dT_relative=3.34341e-09 energy_balance=1.6717e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=3.60988e-06 backward_error=8.35852e-10 linear_relative=0 dT_relative=1.6717e-09 energy_balance=8.35852e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=1.80494e-06 backward_error=4.17926e-10 linear_relative=0 dT_relative=8.35852e-10 energy_balance=4.17926e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=9.02469e-07 backward_error=2.08963e-10 linear_relative=0 dT_relative=4.17926e-10 energy_balance=2.08963e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=4.51235e-07 backward_error=1.04482e-10 linear_relative=0 dT_relative=2.08963e-10 energy_balance=1.04482e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=2.25617e-07 backward_error=5.22407e-11 linear_relative=0 dT_relative=1.04482e-10 energy_balance=5.22407e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=1.12808e-07 backward_error=2.61203e-11 linear_relative=0 dT_relative=5.22409e-11 energy_balance=2.61203e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=5.64041e-08 backward_error=1.30602e-11 linear_relative=0 dT_relative=2.61204e-11 energy_balance=1.30602e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=23 linear_status=0 residual=2.82016e-08 backward_error=6.52997e-12 linear_relative=0 dT_relative=1.30602e-11 energy_balance=6.52997e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=24 linear_status=0 residual=1.41013e-08 backward_error=3.26509e-12 linear_relative=0 dT_relative=6.52993e-12 energy_balance=3.26509e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=25 linear_status=0 residual=7.0504e-09 backward_error=1.63249e-12 linear_relative=0 dT_relative=3.26514e-12 energy_balance=1.63249e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=26 linear_status=0 residual=3.5252e-09 backward_error=8.16246e-13 linear_relative=0 dT_relative=1.63239e-12 energy_balance=8.16246e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=27 linear_status=0 residual=1.7626e-09 backward_error=4.08123e-13 linear_relative=0 dT_relative=8.16376e-13 energy_balance=4.08123e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=10 dT=0.0172782 qrad_delta=0.0901401 qrad_relative=0.00235815 energy_balance=6.6365e-05
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.194 intensity_relative=0.0795956 source_delta=0.0462053 source_relative=0.00120732
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.0606295 backward_error=1.40377e-05 linear_relative=0 dT_relative=2.80758e-05 energy_balance=1.40377e-05 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.0303147 backward_error=7.01881e-06 linear_relative=0 dT_relative=1.40377e-05 energy_balance=7.01881e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.0151574 backward_error=3.50939e-06 linear_relative=0 dT_relative=7.01881e-06 energy_balance=3.50939e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.00757868 backward_error=1.75469e-06 linear_relative=0 dT_relative=3.50939e-06 energy_balance=1.75469e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.00378934 backward_error=8.77346e-07 linear_relative=0 dT_relative=1.75469e-06 energy_balance=8.77346e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.00189467 backward_error=4.38673e-07 linear_relative=0 dT_relative=8.77346e-07 energy_balance=4.38673e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=0.000947335 backward_error=2.19336e-07 linear_relative=0 dT_relative=4.38673e-07 energy_balance=2.19336e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=0.000473668 backward_error=1.09668e-07 linear_relative=0 dT_relative=2.19336e-07 energy_balance=1.09668e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=0.000236834 backward_error=5.48341e-08 linear_relative=0 dT_relative=1.09668e-07 energy_balance=5.48341e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=0.000118417 backward_error=2.7417e-08 linear_relative=0 dT_relative=5.48341e-08 energy_balance=2.7417e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=5.92084e-05 backward_error=1.37085e-08 linear_relative=0 dT_relative=2.7417e-08 energy_balance=1.37085e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=2.96042e-05 backward_error=6.85426e-09 linear_relative=0 dT_relative=1.37085e-08 energy_balance=6.85426e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=1.48021e-05 backward_error=3.42713e-09 linear_relative=0 dT_relative=6.85426e-09 energy_balance=3.42713e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=7.40106e-06 backward_error=1.71357e-09 linear_relative=0 dT_relative=3.42713e-09 energy_balance=1.71357e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=3.70053e-06 backward_error=8.56783e-10 linear_relative=0 dT_relative=1.71357e-09 energy_balance=8.56783e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=1.85026e-06 backward_error=4.28391e-10 linear_relative=0 dT_relative=8.56783e-10 energy_balance=4.28391e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=9.25132e-07 backward_error=2.14196e-10 linear_relative=0 dT_relative=4.28391e-10 energy_balance=2.14196e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=4.62566e-07 backward_error=1.07098e-10 linear_relative=0 dT_relative=2.14196e-10 energy_balance=1.07098e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=2.31283e-07 backward_error=5.3549e-11 linear_relative=0 dT_relative=1.07098e-10 energy_balance=5.3549e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=1.15642e-07 backward_error=2.67745e-11 linear_relative=0 dT_relative=5.3549e-11 energy_balance=2.67745e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=5.78211e-08 backward_error=1.33873e-11 linear_relative=0 dT_relative=2.67744e-11 energy_balance=1.33873e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=2.89106e-08 backward_error=6.69366e-12 linear_relative=0 dT_relative=1.33873e-11 energy_balance=6.69366e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=23 linear_status=0 residual=1.44551e-08 backward_error=3.34678e-12 linear_relative=0 dT_relative=6.69373e-12 energy_balance=3.34678e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=24 linear_status=0 residual=7.22775e-09 backward_error=1.67344e-12 linear_relative=0 dT_relative=3.34678e-12 energy_balance=1.67344e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=25 linear_status=0 residual=3.61388e-09 backward_error=8.36721e-13 linear_relative=0 dT_relative=1.6733e-12 energy_balance=8.36721e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=26 linear_status=0 residual=1.80717e-09 backward_error=4.18413e-13 linear_relative=0 dT_relative=8.36739e-13 energy_balance=4.18413e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=11 dT=0.00885542 qrad_delta=0.0462053 qrad_relative=0.00120732 energy_balance=3.40169e-05
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.195 intensity_relative=0.0795956 source_delta=0.0236841 source_relative=0.00061847
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.0310765 backward_error=7.19494e-06 linear_relative=0 dT_relative=1.439e-05 energy_balance=7.19494e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.0155383 backward_error=3.59746e-06 linear_relative=0 dT_relative=7.19494e-06 energy_balance=3.59746e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.00776914 backward_error=1.79873e-06 linear_relative=0 dT_relative=3.59746e-06 energy_balance=1.79873e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.00388457 backward_error=8.99362e-07 linear_relative=0 dT_relative=1.79873e-06 energy_balance=8.99362e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.00194228 backward_error=4.49681e-07 linear_relative=0 dT_relative=8.99362e-07 energy_balance=4.49681e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.000971142 backward_error=2.2484e-07 linear_relative=0 dT_relative=4.49681e-07 energy_balance=2.2484e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=0.000485571 backward_error=1.1242e-07 linear_relative=0 dT_relative=2.2484e-07 energy_balance=1.1242e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=0.000242785 backward_error=5.62101e-08 linear_relative=0 dT_relative=1.1242e-07 energy_balance=5.62101e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=0.000121393 backward_error=2.8105e-08 linear_relative=0 dT_relative=5.62101e-08 energy_balance=2.8105e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=6.06964e-05 backward_error=1.40525e-08 linear_relative=0 dT_relative=2.8105e-08 energy_balance=1.40525e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=3.03482e-05 backward_error=7.02626e-09 linear_relative=0 dT_relative=1.40525e-08 energy_balance=7.02626e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=1.51741e-05 backward_error=3.51313e-09 linear_relative=0 dT_relative=7.02626e-09 energy_balance=3.51313e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=7.58705e-06 backward_error=1.75657e-09 linear_relative=0 dT_relative=3.51313e-09 energy_balance=1.75657e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=3.79352e-06 backward_error=8.78282e-10 linear_relative=0 dT_relative=1.75656e-09 energy_balance=8.78282e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=1.89676e-06 backward_error=4.39141e-10 linear_relative=0 dT_relative=8.78282e-10 energy_balance=4.39141e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=9.48381e-07 backward_error=2.19571e-10 linear_relative=0 dT_relative=4.39141e-10 energy_balance=2.19571e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=4.74191e-07 backward_error=1.09785e-10 linear_relative=0 dT_relative=2.19571e-10 energy_balance=1.09785e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=2.37095e-07 backward_error=5.48927e-11 linear_relative=0 dT_relative=1.09785e-10 energy_balance=5.48927e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=1.18548e-07 backward_error=2.74463e-11 linear_relative=0 dT_relative=5.48926e-11 energy_balance=2.74463e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=5.9274e-08 backward_error=1.37232e-11 linear_relative=0 dT_relative=2.74462e-11 energy_balance=1.37232e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=2.96372e-08 backward_error=6.86166e-12 linear_relative=0 dT_relative=1.37232e-11 energy_balance=6.86166e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=1.48184e-08 backward_error=3.43078e-12 linear_relative=0 dT_relative=6.86169e-12 energy_balance=3.43078e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=23 linear_status=0 residual=7.4092e-09 backward_error=1.71539e-12 linear_relative=0 dT_relative=3.43075e-12 energy_balance=1.71539e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=24 linear_status=0 residual=3.70437e-09 backward_error=8.57642e-13 linear_relative=0 dT_relative=1.71547e-12 energy_balance=8.57642e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=25 linear_status=0 residual=1.85219e-09 backward_error=4.28821e-13 linear_relative=0 dT_relative=8.57644e-13 energy_balance=4.28821e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=12 dT=0.00453882 qrad_delta=0.0236841 qrad_relative=0.00061847 energy_balance=1.74361e-05
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.196 intensity_relative=0.0795956 source_delta=0.01214 source_relative=0.000316915
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.0159289 backward_error=3.68783e-06 linear_relative=0 dT_relative=7.37569e-06 energy_balance=3.68783e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.00796445 backward_error=1.84391e-06 linear_relative=0 dT_relative=3.68783e-06 energy_balance=1.84391e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.00398222 backward_error=9.21955e-07 linear_relative=0 dT_relative=1.84391e-06 energy_balance=9.21955e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.00199111 backward_error=4.60977e-07 linear_relative=0 dT_relative=9.21955e-07 energy_balance=4.60977e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.000995556 backward_error=2.30489e-07 linear_relative=0 dT_relative=4.60977e-07 energy_balance=2.30489e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.000497778 backward_error=1.15244e-07 linear_relative=0 dT_relative=2.30489e-07 energy_balance=1.15244e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=0.000248889 backward_error=5.76221e-08 linear_relative=0 dT_relative=1.15244e-07 energy_balance=5.76221e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=0.000124444 backward_error=2.88111e-08 linear_relative=0 dT_relative=5.76221e-08 energy_balance=2.88111e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=6.22222e-05 backward_error=1.44055e-08 linear_relative=0 dT_relative=2.88111e-08 energy_balance=1.44055e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=3.11111e-05 backward_error=7.20276e-09 linear_relative=0 dT_relative=1.44055e-08 energy_balance=7.20276e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=1.55556e-05 backward_error=3.60138e-09 linear_relative=0 dT_relative=7.20276e-09 energy_balance=3.60138e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=7.77778e-06 backward_error=1.80069e-09 linear_relative=0 dT_relative=3.60138e-09 energy_balance=1.80069e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=3.88889e-06 backward_error=9.00345e-10 linear_relative=0 dT_relative=1.80069e-09 energy_balance=9.00345e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=1.94444e-06 backward_error=4.50173e-10 linear_relative=0 dT_relative=9.00346e-10 energy_balance=4.50173e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=9.72222e-07 backward_error=2.25086e-10 linear_relative=0 dT_relative=4.50173e-10 energy_balance=2.25086e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=4.86111e-07 backward_error=1.12543e-10 linear_relative=0 dT_relative=2.25086e-10 energy_balance=1.12543e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=2.43055e-07 backward_error=5.62715e-11 linear_relative=0 dT_relative=1.12543e-10 energy_balance=5.62715e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=1.21528e-07 backward_error=2.81357e-11 linear_relative=0 dT_relative=5.62714e-11 energy_balance=2.81357e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=6.07633e-08 backward_error=1.40678e-11 linear_relative=0 dT_relative=2.81359e-11 energy_balance=1.40678e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=3.03817e-08 backward_error=7.03388e-12 linear_relative=0 dT_relative=1.40679e-11 energy_balance=7.03388e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=1.51908e-08 backward_error=3.51694e-12 linear_relative=0 dT_relative=7.03379e-12 energy_balance=3.51694e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=7.59519e-09 backward_error=1.75842e-12 linear_relative=0 dT_relative=3.51708e-12 energy_balance=1.75842e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=23 linear_status=0 residual=3.79714e-09 backward_error=8.79104e-13 linear_relative=0 dT_relative=1.75854e-12 energy_balance=8.79104e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=24 linear_status=0 residual=1.89857e-09 backward_error=4.39552e-13 linear_relative=0 dT_relative=8.79089e-13 energy_balance=4.39552e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=13 dT=0.00232642 qrad_delta=0.01214 qrad_relative=0.000316915 energy_balance=8.93732e-06
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.196 intensity_relative=0.0795956 source_delta=0.00622269 source_relative=0.000162417
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.00816471 backward_error=1.89026e-06 linear_relative=0 dT_relative=3.78052e-06 energy_balance=1.89026e-06 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.00408235 backward_error=9.45128e-07 linear_relative=0 dT_relative=1.89026e-06 energy_balance=9.45128e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.00204118 backward_error=4.72564e-07 linear_relative=0 dT_relative=9.45128e-07 energy_balance=4.72564e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.00102059 backward_error=2.36282e-07 linear_relative=0 dT_relative=4.72564e-07 energy_balance=2.36282e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.000510294 backward_error=1.18141e-07 linear_relative=0 dT_relative=2.36282e-07 energy_balance=1.18141e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.000255147 backward_error=5.90704e-08 linear_relative=0 dT_relative=1.18141e-07 energy_balance=5.90704e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=0.000127574 backward_error=2.95352e-08 linear_relative=0 dT_relative=5.90704e-08 energy_balance=2.95352e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=6.37868e-05 backward_error=1.47676e-08 linear_relative=0 dT_relative=2.95352e-08 energy_balance=1.47676e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=3.18934e-05 backward_error=7.3838e-09 linear_relative=0 dT_relative=1.47676e-08 energy_balance=7.3838e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=1.59467e-05 backward_error=3.6919e-09 linear_relative=0 dT_relative=7.3838e-09 energy_balance=3.6919e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=7.97335e-06 backward_error=1.84595e-09 linear_relative=0 dT_relative=3.6919e-09 energy_balance=1.84595e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=3.98667e-06 backward_error=9.22975e-10 linear_relative=0 dT_relative=1.84595e-09 energy_balance=9.22975e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=1.99334e-06 backward_error=4.61488e-10 linear_relative=0 dT_relative=9.22975e-10 energy_balance=4.61488e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=9.96668e-07 backward_error=2.30744e-10 linear_relative=0 dT_relative=4.61488e-10 energy_balance=2.30744e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=4.98334e-07 backward_error=1.15372e-10 linear_relative=0 dT_relative=2.30744e-10 energy_balance=1.15372e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=2.49167e-07 backward_error=5.76859e-11 linear_relative=0 dT_relative=1.15372e-10 energy_balance=5.76859e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=1.24583e-07 backward_error=2.8843e-11 linear_relative=0 dT_relative=5.76859e-11 energy_balance=2.8843e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=6.22917e-08 backward_error=1.44215e-11 linear_relative=0 dT_relative=2.88429e-11 energy_balance=1.44215e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=3.11461e-08 backward_error=7.2108e-12 linear_relative=0 dT_relative=1.44215e-11 energy_balance=7.2108e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=1.55728e-08 backward_error=3.60535e-12 linear_relative=0 dT_relative=7.21076e-12 energy_balance=3.60535e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=7.78664e-09 backward_error=1.80273e-12 linear_relative=0 dT_relative=3.60538e-12 energy_balance=1.80273e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=3.89309e-09 backward_error=9.0131e-13 linear_relative=0 dT_relative=1.80269e-12 energy_balance=9.0131e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=23 linear_status=0 residual=1.94677e-09 backward_error=4.50708e-13 linear_relative=0 dT_relative=9.01255e-13 energy_balance=4.50708e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=14 dT=0.00119245 qrad_delta=0.00622269 qrad_relative=0.000162417 energy_balance=4.58104e-06
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=0.0031896 source_relative=8.3244e-05
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.00418501 backward_error=9.6889e-07 linear_relative=0 dT_relative=1.93778e-06 energy_balance=9.6889e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.00209251 backward_error=4.84445e-07 linear_relative=0 dT_relative=9.6889e-07 energy_balance=4.84445e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.00104625 backward_error=2.42222e-07 linear_relative=0 dT_relative=4.84445e-07 energy_balance=2.42222e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.000523126 backward_error=1.21111e-07 linear_relative=0 dT_relative=2.42222e-07 energy_balance=1.21111e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.000261563 backward_error=6.05556e-08 linear_relative=0 dT_relative=1.21111e-07 energy_balance=6.05556e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=0.000130782 backward_error=3.02778e-08 linear_relative=0 dT_relative=6.05556e-08 energy_balance=3.02778e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=6.53908e-05 backward_error=1.51389e-08 linear_relative=0 dT_relative=3.02778e-08 energy_balance=1.51389e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=3.26954e-05 backward_error=7.56945e-09 linear_relative=0 dT_relative=1.51389e-08 energy_balance=7.56945e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=1.63477e-05 backward_error=3.78472e-09 linear_relative=0 dT_relative=7.56944e-09 energy_balance=3.78472e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=8.17385e-06 backward_error=1.89236e-09 linear_relative=0 dT_relative=3.78472e-09 energy_balance=1.89236e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=4.08692e-06 backward_error=9.46181e-10 linear_relative=0 dT_relative=1.89236e-09 energy_balance=9.46181e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=2.04346e-06 backward_error=4.7309e-10 linear_relative=0 dT_relative=9.46181e-10 energy_balance=4.7309e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=1.02173e-06 backward_error=2.36545e-10 linear_relative=0 dT_relative=4.7309e-10 energy_balance=2.36545e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=5.10865e-07 backward_error=1.18273e-10 linear_relative=0 dT_relative=2.36545e-10 energy_balance=1.18273e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=2.55433e-07 backward_error=5.91363e-11 linear_relative=0 dT_relative=1.18273e-10 energy_balance=5.91363e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=1.27716e-07 backward_error=2.95681e-11 linear_relative=0 dT_relative=5.91363e-11 energy_balance=2.95681e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=6.38584e-08 backward_error=1.47841e-11 linear_relative=0 dT_relative=2.95681e-11 energy_balance=1.47841e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=3.19292e-08 backward_error=7.39205e-12 linear_relative=0 dT_relative=1.47841e-11 energy_balance=7.39205e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=1.59648e-08 backward_error=3.69608e-12 linear_relative=0 dT_relative=7.39206e-12 energy_balance=3.69608e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=7.98218e-09 backward_error=1.84799e-12 linear_relative=0 dT_relative=3.69603e-12 energy_balance=1.84799e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=3.99086e-09 backward_error=9.23941e-13 linear_relative=0 dT_relative=1.8481e-12 energy_balance=9.23941e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=22 linear_status=0 residual=1.99543e-09 backward_error=4.6197e-13 linear_relative=0 dT_relative=9.23962e-13 energy_balance=4.6197e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=15 dT=0.000611214 qrad_delta=0.0031896 qrad_relative=8.3244e-05 energy_balance=2.34812e-06
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=0.00163491 source_relative=4.2667e-05
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.00214513 backward_error=4.96626e-07 linear_relative=0 dT_relative=9.93253e-07 energy_balance=4.96626e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.00107256 backward_error=2.48313e-07 linear_relative=0 dT_relative=4.96626e-07 energy_balance=2.48313e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.000536282 backward_error=1.24157e-07 linear_relative=0 dT_relative=2.48313e-07 energy_balance=1.24157e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.000268141 backward_error=6.20783e-08 linear_relative=0 dT_relative=1.24157e-07 energy_balance=6.20783e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=0.000134071 backward_error=3.10391e-08 linear_relative=0 dT_relative=6.20783e-08 energy_balance=3.10391e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=6.70353e-05 backward_error=1.55196e-08 linear_relative=0 dT_relative=3.10391e-08 energy_balance=1.55196e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=3.35176e-05 backward_error=7.75978e-09 linear_relative=0 dT_relative=1.55196e-08 energy_balance=7.75978e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=1.67588e-05 backward_error=3.87989e-09 linear_relative=0 dT_relative=7.75978e-09 energy_balance=3.87989e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=8.37941e-06 backward_error=1.93995e-09 linear_relative=0 dT_relative=3.87989e-09 energy_balance=1.93995e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=4.1897e-06 backward_error=9.69973e-10 linear_relative=0 dT_relative=1.93995e-09 energy_balance=9.69973e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=2.09485e-06 backward_error=4.84986e-10 linear_relative=0 dT_relative=9.69973e-10 energy_balance=4.84986e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=1.04743e-06 backward_error=2.42493e-10 linear_relative=0 dT_relative=4.84987e-10 energy_balance=2.42493e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=5.23713e-07 backward_error=1.21247e-10 linear_relative=0 dT_relative=2.42493e-10 energy_balance=1.21247e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=2.61857e-07 backward_error=6.06234e-11 linear_relative=0 dT_relative=1.21246e-10 energy_balance=6.06234e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=1.30929e-07 backward_error=3.03117e-11 linear_relative=0 dT_relative=6.06232e-11 energy_balance=3.03117e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=6.54641e-08 backward_error=1.51558e-11 linear_relative=0 dT_relative=3.03118e-11 energy_balance=1.51558e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=3.27323e-08 backward_error=7.57796e-12 linear_relative=0 dT_relative=1.51557e-11 energy_balance=7.57796e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=1.63659e-08 backward_error=3.78893e-12 linear_relative=0 dT_relative=7.57804e-12 energy_balance=3.78893e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=8.18318e-09 backward_error=1.89452e-12 linear_relative=0 dT_relative=3.78884e-12 energy_balance=1.89452e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=4.09182e-09 backward_error=9.47311e-13 linear_relative=0 dT_relative=1.89442e-12 energy_balance=9.47311e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=21 linear_status=0 residual=2.04591e-09 backward_error=4.73655e-13 linear_relative=0 dT_relative=9.4721e-13 energy_balance=4.73655e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=16 dT=0.000313292 qrad_delta=0.00163491 qrad_relative=4.2667e-05 energy_balance=1.20359e-06
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=0.000838014 source_relative=2.18696e-05
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.00109954 backward_error=2.54558e-07 linear_relative=0 dT_relative=5.09115e-07 energy_balance=2.54558e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.000549769 backward_error=1.27279e-07 linear_relative=0 dT_relative=2.54558e-07 energy_balance=1.27279e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.000274884 backward_error=6.36394e-08 linear_relative=0 dT_relative=1.27279e-07 energy_balance=6.36394e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=0.000137442 backward_error=3.18197e-08 linear_relative=0 dT_relative=6.36394e-08 energy_balance=3.18197e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=6.87211e-05 backward_error=1.59098e-08 linear_relative=0 dT_relative=3.18197e-08 energy_balance=1.59098e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=3.43606e-05 backward_error=7.95492e-09 linear_relative=0 dT_relative=1.59098e-08 energy_balance=7.95492e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=1.71803e-05 backward_error=3.97746e-09 linear_relative=0 dT_relative=7.95492e-09 energy_balance=3.97746e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=8.59014e-06 backward_error=1.98873e-09 linear_relative=0 dT_relative=3.97746e-09 energy_balance=1.98873e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=4.29507e-06 backward_error=9.94365e-10 linear_relative=0 dT_relative=1.98873e-09 energy_balance=9.94365e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=2.14753e-06 backward_error=4.97183e-10 linear_relative=0 dT_relative=9.94365e-10 energy_balance=4.97183e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=1.07377e-06 backward_error=2.48591e-10 linear_relative=0 dT_relative=4.97183e-10 energy_balance=2.48591e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=5.36883e-07 backward_error=1.24296e-10 linear_relative=0 dT_relative=2.48591e-10 energy_balance=1.24296e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=2.68441e-07 backward_error=6.21477e-11 linear_relative=0 dT_relative=1.24296e-10 energy_balance=6.21477e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=1.3422e-07 backward_error=3.10738e-11 linear_relative=0 dT_relative=6.21479e-11 energy_balance=3.10738e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=6.71107e-08 backward_error=1.5537e-11 linear_relative=0 dT_relative=3.10737e-11 energy_balance=1.5537e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=3.35549e-08 backward_error=7.7684e-12 linear_relative=0 dT_relative=1.55371e-11 energy_balance=7.7684e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=1.67779e-08 backward_error=3.88431e-12 linear_relative=0 dT_relative=7.76835e-12 energy_balance=3.88431e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=8.38872e-09 backward_error=1.9421e-12 linear_relative=0 dT_relative=3.88435e-12 energy_balance=1.9421e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=4.19459e-09 backward_error=9.71103e-13 linear_relative=0 dT_relative=1.942e-12 energy_balance=9.71103e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=20 linear_status=0 residual=2.09729e-09 backward_error=4.85551e-13 linear_relative=0 dT_relative=9.70998e-13 energy_balance=4.85551e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=17 dT=0.000160585 qrad_delta=0.000838014 qrad_relative=2.18696e-05 energy_balance=6.1693e-07
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=0.000429545 source_relative=1.12097e-05
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.000563595 backward_error=1.3048e-07 linear_relative=0 dT_relative=2.60959e-07 energy_balance=1.3048e-07 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.000281798 backward_error=6.52399e-08 linear_relative=0 dT_relative=1.3048e-07 energy_balance=6.52399e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=0.000140899 backward_error=3.26199e-08 linear_relative=0 dT_relative=6.52399e-08 energy_balance=3.26199e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=7.04494e-05 backward_error=1.631e-08 linear_relative=0 dT_relative=3.26199e-08 energy_balance=1.631e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=3.52247e-05 backward_error=8.15498e-09 linear_relative=0 dT_relative=1.631e-08 energy_balance=8.15498e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=1.76124e-05 backward_error=4.07749e-09 linear_relative=0 dT_relative=8.15498e-09 energy_balance=4.07749e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=8.80618e-06 backward_error=2.03875e-09 linear_relative=0 dT_relative=4.07749e-09 energy_balance=2.03875e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=4.40309e-06 backward_error=1.01937e-09 linear_relative=0 dT_relative=2.03875e-09 energy_balance=1.01937e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=2.20154e-06 backward_error=5.09686e-10 linear_relative=0 dT_relative=1.01937e-09 energy_balance=5.09686e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=1.10077e-06 backward_error=2.54843e-10 linear_relative=0 dT_relative=5.09686e-10 energy_balance=2.54843e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=5.50386e-07 backward_error=1.27422e-10 linear_relative=0 dT_relative=2.54843e-10 energy_balance=1.27422e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=2.75193e-07 backward_error=6.37107e-11 linear_relative=0 dT_relative=1.27422e-10 energy_balance=6.37107e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=1.37597e-07 backward_error=3.18554e-11 linear_relative=0 dT_relative=6.37107e-11 energy_balance=3.18554e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=6.87978e-08 backward_error=1.59276e-11 linear_relative=0 dT_relative=3.18555e-11 energy_balance=1.59276e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=3.43989e-08 backward_error=7.9638e-12 linear_relative=0 dT_relative=1.59278e-11 energy_balance=7.9638e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=1.71995e-08 backward_error=3.9819e-12 linear_relative=0 dT_relative=7.9637e-12 energy_balance=3.9819e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=8.59927e-09 backward_error=1.99084e-12 linear_relative=0 dT_relative=3.98203e-12 energy_balance=1.99084e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=4.29964e-09 backward_error=9.95422e-13 linear_relative=0 dT_relative=1.99101e-12 energy_balance=9.95422e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=2.14959e-09 backward_error=4.97658e-13 linear_relative=0 dT_relative=9.95507e-13 energy_balance=4.97658e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=18 dT=8.23117e-05 qrad_delta=0.000429545 qrad_relative=1.12097e-05 energy_balance=3.16223e-07
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=0.000220174 source_relative=5.74575e-06
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.000288885 backward_error=6.68807e-08 linear_relative=0 dT_relative=1.33761e-07 energy_balance=6.68807e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=0.000144443 backward_error=3.34403e-08 linear_relative=0 dT_relative=6.68807e-08 energy_balance=3.34403e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=7.22213e-05 backward_error=1.67202e-08 linear_relative=0 dT_relative=3.34403e-08 energy_balance=1.67202e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=3.61106e-05 backward_error=8.36008e-09 linear_relative=0 dT_relative=1.67202e-08 energy_balance=8.36008e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=1.80553e-05 backward_error=4.18004e-09 linear_relative=0 dT_relative=8.36008e-09 energy_balance=4.18004e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=9.02766e-06 backward_error=2.09002e-09 linear_relative=0 dT_relative=4.18004e-09 energy_balance=2.09002e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=4.51383e-06 backward_error=1.04501e-09 linear_relative=0 dT_relative=2.09002e-09 energy_balance=1.04501e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=2.25692e-06 backward_error=5.22505e-10 linear_relative=0 dT_relative=1.04501e-09 energy_balance=5.22505e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=1.12846e-06 backward_error=2.61253e-10 linear_relative=0 dT_relative=5.22505e-10 energy_balance=2.61253e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=5.64229e-07 backward_error=1.30626e-10 linear_relative=0 dT_relative=2.61253e-10 energy_balance=1.30626e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=2.82115e-07 backward_error=6.53132e-11 linear_relative=0 dT_relative=1.30626e-10 energy_balance=6.53132e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=1.41058e-07 backward_error=3.26567e-11 linear_relative=0 dT_relative=6.53131e-11 energy_balance=3.26567e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=7.05286e-08 backward_error=1.63283e-11 linear_relative=0 dT_relative=3.26568e-11 energy_balance=1.63283e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=3.52643e-08 backward_error=8.16414e-12 linear_relative=0 dT_relative=1.63282e-11 energy_balance=8.16414e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=1.76324e-08 backward_error=4.08212e-12 linear_relative=0 dT_relative=8.1641e-12 energy_balance=4.08212e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=8.81619e-09 backward_error=2.04106e-12 linear_relative=0 dT_relative=4.08223e-12 energy_balance=2.04106e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=4.40787e-09 backward_error=1.02048e-12 linear_relative=0 dT_relative=2.04111e-12 energy_balance=1.02048e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=2.20416e-09 backward_error=5.10292e-13 linear_relative=0 dT_relative=1.02038e-12 energy_balance=5.10292e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=19 linear_status=0 residual=1.10185e-09 backward_error=2.55093e-13 linear_relative=0 dT_relative=5.10369e-13 energy_balance=2.55093e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=19 dT=4.21909e-05 qrad_delta=0.000220174 qrad_relative=5.74575e-06 energy_balance=1.62088e-07
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=0.000112856 source_relative=2.94512e-06
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0.000148075 backward_error=3.42814e-08 linear_relative=0 dT_relative=6.85628e-08 energy_balance=3.42814e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=7.40377e-05 backward_error=1.71407e-08 linear_relative=0 dT_relative=3.42814e-08 energy_balance=1.71407e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=3.70189e-05 backward_error=8.57035e-09 linear_relative=0 dT_relative=1.71407e-08 energy_balance=8.57035e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=1.85094e-05 backward_error=4.28517e-09 linear_relative=0 dT_relative=8.57035e-09 energy_balance=4.28517e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=9.25472e-06 backward_error=2.14259e-09 linear_relative=0 dT_relative=4.28517e-09 energy_balance=2.14259e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=4.62736e-06 backward_error=1.07129e-09 linear_relative=0 dT_relative=2.14259e-09 energy_balance=1.07129e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=2.31368e-06 backward_error=5.35647e-10 linear_relative=0 dT_relative=1.07129e-09 energy_balance=5.35647e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=1.15684e-06 backward_error=2.67823e-10 linear_relative=0 dT_relative=5.35647e-10 energy_balance=2.67823e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=5.7842e-07 backward_error=1.33912e-10 linear_relative=0 dT_relative=2.67823e-10 energy_balance=1.33912e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=2.8921e-07 backward_error=6.69558e-11 linear_relative=0 dT_relative=1.33912e-10 energy_balance=6.69558e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=1.44605e-07 backward_error=3.34778e-11 linear_relative=0 dT_relative=6.69558e-11 energy_balance=3.34778e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=7.23021e-08 backward_error=1.67389e-11 linear_relative=0 dT_relative=3.3478e-11 energy_balance=1.67389e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=3.61511e-08 backward_error=8.36944e-12 linear_relative=0 dT_relative=1.67389e-11 energy_balance=8.36944e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=1.80758e-08 backward_error=4.18477e-12 linear_relative=0 dT_relative=8.36936e-12 energy_balance=4.18477e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=9.03765e-09 backward_error=2.09233e-12 linear_relative=0 dT_relative=4.18477e-12 energy_balance=2.09233e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=4.51882e-09 backward_error=1.04617e-12 linear_relative=0 dT_relative=2.09248e-12 energy_balance=1.04617e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=2.25918e-09 backward_error=5.23031e-13 linear_relative=0 dT_relative=1.04615e-12 energy_balance=5.23031e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=18 linear_status=0 residual=1.12959e-09 backward_error=2.61515e-13 linear_relative=0 dT_relative=5.23164e-13 energy_balance=2.61515e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=20 dT=2.1626e-05 qrad_delta=0.000112856 qrad_relative=2.94512e-06 energy_balance=8.30824e-08
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=5.7847e-05 source_relative=1.50959e-06
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=7.59e-05 backward_error=1.75718e-08 linear_relative=0 dT_relative=3.51437e-08 energy_balance=1.75718e-08 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=3.795e-05 backward_error=8.78591e-09 linear_relative=0 dT_relative=1.75718e-08 energy_balance=8.78591e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=1.8975e-05 backward_error=4.39296e-09 linear_relative=0 dT_relative=8.78591e-09 energy_balance=4.39296e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=9.4875e-06 backward_error=2.19648e-09 linear_relative=0 dT_relative=4.39296e-09 energy_balance=2.19648e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=4.74375e-06 backward_error=1.09824e-09 linear_relative=0 dT_relative=2.19648e-09 energy_balance=1.09824e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=2.37187e-06 backward_error=5.4912e-10 linear_relative=0 dT_relative=1.09824e-09 energy_balance=5.4912e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=1.18594e-06 backward_error=2.7456e-10 linear_relative=0 dT_relative=5.4912e-10 energy_balance=2.7456e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=5.92969e-07 backward_error=1.3728e-10 linear_relative=0 dT_relative=2.7456e-10 energy_balance=1.3728e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=2.96485e-07 backward_error=6.86401e-11 linear_relative=0 dT_relative=1.3728e-10 energy_balance=6.86401e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=1.48243e-07 backward_error=3.43201e-11 linear_relative=0 dT_relative=6.86399e-11 energy_balance=3.43201e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=7.41211e-08 backward_error=1.716e-11 linear_relative=0 dT_relative=3.43201e-11 energy_balance=1.716e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=3.70605e-08 backward_error=8.58e-12 linear_relative=0 dT_relative=1.71601e-11 energy_balance=8.58e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=1.853e-08 backward_error=4.28994e-12 linear_relative=0 dT_relative=8.58003e-12 energy_balance=4.28994e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=9.26548e-09 backward_error=2.14508e-12 linear_relative=0 dT_relative=4.28984e-12 energy_balance=2.14508e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=4.63297e-09 backward_error=1.07259e-12 linear_relative=0 dT_relative=2.14492e-12 energy_balance=1.07259e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=2.31648e-09 backward_error=5.36296e-13 linear_relative=0 dT_relative=1.07264e-12 energy_balance=5.36296e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=17 linear_status=0 residual=1.15824e-09 backward_error=2.68148e-13 linear_relative=0 dT_relative=5.3632e-13 energy_balance=2.68148e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=21 dT=1.10849e-05 qrad_delta=5.7847e-05 qrad_relative=1.50959e-06 energy_balance=4.25862e-08
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=2.96509e-05 source_relative=7.73779e-07
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=3.89047e-05 backward_error=9.00694e-09 linear_relative=0 dT_relative=1.80139e-08 energy_balance=9.00694e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=1.94523e-05 backward_error=4.50347e-09 linear_relative=0 dT_relative=9.00694e-09 energy_balance=4.50347e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=9.72617e-06 backward_error=2.25173e-09 linear_relative=0 dT_relative=4.50347e-09 energy_balance=2.25173e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=4.86309e-06 backward_error=1.12587e-09 linear_relative=0 dT_relative=2.25173e-09 energy_balance=1.12587e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=2.43154e-06 backward_error=5.62933e-10 linear_relative=0 dT_relative=1.12587e-09 energy_balance=5.62933e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=1.21577e-06 backward_error=2.81467e-10 linear_relative=0 dT_relative=5.62933e-10 energy_balance=2.81467e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=6.07885e-07 backward_error=1.40733e-10 linear_relative=0 dT_relative=2.81467e-10 energy_balance=1.40733e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=3.03943e-07 backward_error=7.03667e-11 linear_relative=0 dT_relative=1.40733e-10 energy_balance=7.03667e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=1.51972e-07 backward_error=3.51834e-11 linear_relative=0 dT_relative=7.03667e-11 energy_balance=3.51834e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=7.59856e-08 backward_error=1.75916e-11 linear_relative=0 dT_relative=3.51834e-11 energy_balance=1.75916e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=3.79928e-08 backward_error=8.79582e-12 linear_relative=0 dT_relative=1.75917e-11 energy_balance=8.79582e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=1.89966e-08 backward_error=4.39796e-12 linear_relative=0 dT_relative=8.79575e-12 energy_balance=4.39796e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=9.49831e-09 backward_error=2.19898e-12 linear_relative=0 dT_relative=4.39797e-12 energy_balance=2.19898e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=4.74893e-09 backward_error=1.09944e-12 linear_relative=0 dT_relative=2.19898e-12 energy_balance=1.09944e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=2.37424e-09 backward_error=5.49666e-13 linear_relative=0 dT_relative=1.09949e-12 energy_balance=5.49666e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=16 linear_status=0 residual=1.18689e-09 backward_error=2.7478e-13 linear_relative=0 dT_relative=5.49836e-13 energy_balance=2.7478e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=22 dT=5.68185e-06 qrad_delta=2.96509e-05 qrad_relative=7.73779e-07 energy_balance=2.18288e-08
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=1.51983e-05 source_relative=3.96619e-07
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=1.99418e-05 backward_error=4.61679e-09 linear_relative=0 dT_relative=9.23359e-09 energy_balance=4.61679e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=9.97092e-06 backward_error=2.3084e-09 linear_relative=0 dT_relative=4.61679e-09 energy_balance=2.3084e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=4.98546e-06 backward_error=1.1542e-09 linear_relative=0 dT_relative=2.3084e-09 energy_balance=1.1542e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=2.49273e-06 backward_error=5.77099e-10 linear_relative=0 dT_relative=1.1542e-09 energy_balance=5.77099e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=1.24636e-06 backward_error=2.8855e-10 linear_relative=0 dT_relative=5.77099e-10 energy_balance=2.8855e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=6.23183e-07 backward_error=1.44275e-10 linear_relative=0 dT_relative=2.8855e-10 energy_balance=1.44275e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=3.11591e-07 backward_error=7.21374e-11 linear_relative=0 dT_relative=1.44275e-10 energy_balance=7.21374e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=1.55796e-07 backward_error=3.60687e-11 linear_relative=0 dT_relative=7.21373e-11 energy_balance=3.60687e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=7.78978e-08 backward_error=1.80343e-11 linear_relative=0 dT_relative=3.60688e-11 energy_balance=1.80343e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=3.89487e-08 backward_error=9.01712e-12 linear_relative=0 dT_relative=1.80343e-11 energy_balance=9.01712e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=1.94746e-08 backward_error=4.50861e-12 linear_relative=0 dT_relative=9.01706e-12 energy_balance=4.50861e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=9.73751e-09 backward_error=2.25436e-12 linear_relative=0 dT_relative=4.50862e-12 energy_balance=2.25436e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=4.86853e-09 backward_error=1.12713e-12 linear_relative=0 dT_relative=2.25431e-12 energy_balance=1.12713e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=2.43426e-09 backward_error=5.63563e-13 linear_relative=0 dT_relative=1.12724e-12 energy_balance=5.63563e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=15 linear_status=0 residual=1.2169e-09 backward_error=2.81729e-13 linear_relative=0 dT_relative=5.63532e-13 energy_balance=2.81729e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=23 dT=2.91237e-06 qrad_delta=1.51983e-05 qrad_relative=3.96619e-07 energy_balance=1.11891e-08
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=7.79025e-06 source_relative=2.03297e-07
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=1.0222e-05 backward_error=2.36652e-09 linear_relative=0 dT_relative=4.73304e-09 energy_balance=2.36652e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=5.11099e-06 backward_error=1.18326e-09 linear_relative=0 dT_relative=2.36652e-09 energy_balance=1.18326e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=2.55549e-06 backward_error=5.9163e-10 linear_relative=0 dT_relative=1.18326e-09 energy_balance=5.9163e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=1.27775e-06 backward_error=2.95815e-10 linear_relative=0 dT_relative=5.9163e-10 energy_balance=2.95815e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=6.38874e-07 backward_error=1.47907e-10 linear_relative=0 dT_relative=2.95815e-10 energy_balance=1.47907e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=3.19437e-07 backward_error=7.39537e-11 linear_relative=0 dT_relative=1.47907e-10 energy_balance=7.39537e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=1.59719e-07 backward_error=3.69769e-11 linear_relative=0 dT_relative=7.39537e-11 energy_balance=3.69769e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=7.98595e-08 backward_error=1.84885e-11 linear_relative=0 dT_relative=3.69769e-11 energy_balance=1.84885e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=3.993e-08 backward_error=9.24431e-12 linear_relative=0 dT_relative=1.84884e-11 energy_balance=9.24431e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=1.99652e-08 backward_error=4.62221e-12 linear_relative=0 dT_relative=9.24413e-12 energy_balance=4.62221e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=9.98261e-09 backward_error=2.3111e-12 linear_relative=0 dT_relative=4.62215e-12 energy_balance=2.3111e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=4.99131e-09 backward_error=1.15555e-12 linear_relative=0 dT_relative=2.31108e-12 energy_balance=1.15555e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=2.49565e-09 backward_error=5.77776e-13 linear_relative=0 dT_relative=1.15554e-12 energy_balance=5.77776e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=14 linear_status=0 residual=1.24783e-09 backward_error=2.88888e-13 linear_relative=0 dT_relative=5.77769e-13 energy_balance=2.88888e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=24 dT=1.4928e-06 qrad_delta=7.79025e-06 qrad_relative=2.03297e-07 energy_balance=5.7355e-09
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=3.99308e-06 source_relative=1.04204e-07
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=5.23983e-06 backward_error=1.21309e-09 linear_relative=0 dT_relative=2.42618e-09 energy_balance=1.21309e-09 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=2.61992e-06 backward_error=6.06544e-10 linear_relative=0 dT_relative=1.21309e-09 energy_balance=6.06544e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=1.30996e-06 backward_error=3.03272e-10 linear_relative=0 dT_relative=6.06544e-10 energy_balance=3.03272e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=6.54979e-07 backward_error=1.51636e-10 linear_relative=0 dT_relative=3.03272e-10 energy_balance=1.51636e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=3.27489e-07 backward_error=7.58179e-11 linear_relative=0 dT_relative=1.51636e-10 energy_balance=7.58179e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=1.63745e-07 backward_error=3.79091e-11 linear_relative=0 dT_relative=7.58179e-11 energy_balance=3.79091e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=8.18727e-08 backward_error=1.89546e-11 linear_relative=0 dT_relative=3.79089e-11 energy_balance=1.89546e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=4.09364e-08 backward_error=9.47729e-12 linear_relative=0 dT_relative=1.89546e-11 energy_balance=9.47729e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=2.04682e-08 backward_error=4.73865e-12 linear_relative=0 dT_relative=9.47714e-12 energy_balance=4.73865e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=1.02345e-08 backward_error=2.36943e-12 linear_relative=0 dT_relative=4.73857e-12 energy_balance=2.36943e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=5.11727e-09 backward_error=1.18471e-12 linear_relative=0 dT_relative=2.36947e-12 energy_balance=1.18471e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=2.55841e-09 backward_error=5.92305e-13 linear_relative=0 dT_relative=1.18473e-12 energy_balance=5.92305e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=13 linear_status=0 residual=1.27966e-09 backward_error=2.96258e-13 linear_relative=0 dT_relative=5.92186e-13 energy_balance=2.96258e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=25 dT=7.65171e-07 qrad_delta=3.99308e-06 qrad_relative=1.04204e-07 energy_balance=2.94012e-09
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=2.04674e-06 source_relative=5.34124e-08
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=2.68612e-06 backward_error=6.2187e-10 linear_relative=0 dT_relative=1.24374e-09 energy_balance=6.2187e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=1.34306e-06 backward_error=3.10935e-10 linear_relative=0 dT_relative=6.2187e-10 energy_balance=3.10935e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=6.71529e-07 backward_error=1.55468e-10 linear_relative=0 dT_relative=3.10935e-10 energy_balance=1.55468e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=3.35765e-07 backward_error=7.77338e-11 linear_relative=0 dT_relative=1.55467e-10 energy_balance=7.77338e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=1.67882e-07 backward_error=3.88669e-11 linear_relative=0 dT_relative=7.77339e-11 energy_balance=3.88669e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=8.39409e-08 backward_error=1.94334e-11 linear_relative=0 dT_relative=3.8867e-11 energy_balance=1.94334e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=4.19705e-08 backward_error=9.7167e-12 linear_relative=0 dT_relative=1.94333e-11 energy_balance=9.7167e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=2.09852e-08 backward_error=4.85835e-12 linear_relative=0 dT_relative=9.71683e-12 energy_balance=4.85835e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=1.04928e-08 backward_error=2.42923e-12 linear_relative=0 dT_relative=4.85823e-12 energy_balance=2.42923e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=5.24642e-09 backward_error=1.21461e-12 linear_relative=0 dT_relative=2.4293e-12 energy_balance=1.21461e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=2.62298e-09 backward_error=6.07254e-13 linear_relative=0 dT_relative=1.21465e-12 energy_balance=6.07254e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=12 linear_status=0 residual=1.31149e-09 backward_error=3.03627e-13 linear_relative=0 dT_relative=6.07324e-13 energy_balance=3.03627e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=26 dT=3.92205e-07 qrad_delta=2.04674e-06 qrad_relative=5.34124e-08 energy_balance=1.50728e-09
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=1.0491e-06 source_relative=2.73777e-08
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=1.37715e-06 backward_error=3.18829e-10 linear_relative=0 dT_relative=6.37658e-10 energy_balance=3.18829e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=6.88577e-07 backward_error=1.59414e-10 linear_relative=0 dT_relative=3.18829e-10 energy_balance=1.59414e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=3.44289e-07 backward_error=7.97073e-11 linear_relative=0 dT_relative=1.59414e-10 energy_balance=7.97073e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=1.72144e-07 backward_error=3.98536e-11 linear_relative=0 dT_relative=7.97073e-11 energy_balance=3.98536e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=8.60723e-08 backward_error=1.99268e-11 linear_relative=0 dT_relative=3.98536e-11 energy_balance=1.99268e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=4.30359e-08 backward_error=9.96337e-12 linear_relative=0 dT_relative=1.99269e-11 energy_balance=9.96337e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=2.15177e-08 backward_error=4.98163e-12 linear_relative=0 dT_relative=9.96336e-12 energy_balance=4.98163e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=1.07589e-08 backward_error=2.49082e-12 linear_relative=0 dT_relative=4.98168e-12 energy_balance=2.49082e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=5.37966e-09 backward_error=1.24546e-12 linear_relative=0 dT_relative=2.49075e-12 energy_balance=1.24546e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=2.68983e-09 backward_error=6.2273e-13 linear_relative=0 dT_relative=1.24547e-12 energy_balance=6.2273e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=11 linear_status=0 residual=1.34514e-09 backward_error=3.11418e-13 linear_relative=0 dT_relative=6.22643e-13 energy_balance=3.11418e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=27 dT=2.01032e-07 qrad_delta=1.0491e-06 qrad_relative=2.73777e-08 energy_balance=7.7285e-10
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=5.37737e-07 source_relative=1.40329e-08
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=7.06222e-07 backward_error=1.635e-10 linear_relative=0 dT_relative=3.26999e-10 energy_balance=1.635e-10 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=3.53111e-07 backward_error=8.17498e-11 linear_relative=0 dT_relative=1.63499e-10 energy_balance=8.17498e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=1.76556e-07 backward_error=4.08749e-11 linear_relative=0 dT_relative=8.17496e-11 energy_balance=4.08749e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=8.82778e-08 backward_error=2.04375e-11 linear_relative=0 dT_relative=4.08749e-11 energy_balance=2.04375e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=4.41391e-08 backward_error=1.02188e-11 linear_relative=0 dT_relative=2.04375e-11 energy_balance=1.02188e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=2.20698e-08 backward_error=5.10944e-12 linear_relative=0 dT_relative=1.02187e-11 energy_balance=5.10944e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=1.10349e-08 backward_error=2.55472e-12 linear_relative=0 dT_relative=5.10927e-12 energy_balance=2.55472e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=5.51745e-09 backward_error=1.27736e-12 linear_relative=0 dT_relative=2.55473e-12 energy_balance=1.27736e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=2.75895e-09 backward_error=6.38733e-13 linear_relative=0 dT_relative=1.27736e-12 energy_balance=6.38733e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=10 linear_status=0 residual=1.3797e-09 backward_error=3.19419e-13 linear_relative=0 dT_relative=6.38682e-13 energy_balance=3.19419e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=28 dT=1.03041e-07 qrad_delta=5.37737e-07 qrad_relative=1.40329e-08 energy_balance=3.96407e-10
RADIATION_RESIDUAL: iteration=1 intensity_delta=148.197 intensity_relative=0.0795956 source_delta=2.75624e-07 source_relative=7.19275e-09
RADIATION_RESIDUAL: iteration=2 intensity_delta=0 intensity_relative=0 source_delta=0 source_relative=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=3.62327e-07 backward_error=8.38834e-11 linear_relative=0 dT_relative=1.67767e-10 energy_balance=8.38834e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=2 linear_status=0 residual=1.81164e-07 backward_error=4.19417e-11 linear_relative=0 dT_relative=8.38836e-11 energy_balance=4.19417e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=3 linear_status=0 residual=9.05816e-08 backward_error=2.09708e-11 linear_relative=0 dT_relative=4.19418e-11 energy_balance=2.09708e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=4 linear_status=0 residual=4.52906e-08 backward_error=1.04853e-11 linear_relative=0 dT_relative=2.09709e-11 energy_balance=1.04853e-11 accepted_state=0
THERMAL_RESIDUAL: iteration=5 linear_status=0 residual=2.26455e-08 backward_error=5.24273e-12 linear_relative=0 dT_relative=1.04853e-11 energy_balance=5.24273e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=6 linear_status=0 residual=1.13223e-08 backward_error=2.62126e-12 linear_relative=0 dT_relative=5.24281e-12 energy_balance=2.62126e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=7 linear_status=0 residual=5.66115e-09 backward_error=1.31063e-12 linear_relative=0 dT_relative=2.62141e-12 energy_balance=1.31063e-12 accepted_state=0
THERMAL_RESIDUAL: iteration=8 linear_status=0 residual=2.8308e-09 backward_error=6.55367e-13 linear_relative=0 dT_relative=1.31052e-12 energy_balance=6.55367e-13 accepted_state=0
THERMAL_RESIDUAL: iteration=9 linear_status=0 residual=1.41563e-09 backward_error=3.27736e-13 linear_relative=0 dT_relative=6.55262e-13 energy_balance=3.27736e-13 accepted_state=1
THERMAL_RADIATION_RESIDUAL: iteration=29 dT=5.28137e-08 qrad_delta=2.75624e-07 qrad_relative=7.19275e-09 energy_balance=2.03459e-10
THERMAL_CHALLENGE: case=radiation_energy_coupled iterations=29 source_residual=2.75624e-07 energy_balance=2.03459e-10
  PASS: thermal_radiation_coupled_conservation
RADIATION_RESIDUAL: angular=S2 directions=8 quadrature=PASS
RADIATION_RESIDUAL: angular=S4 directions=32 quadrature=PASS
RADIATION_RESIDUAL: angular=S6 directions=72 quadrature=PASS
RADIATION_RESIDUAL: angular=S8 directions=128 quadrature=PASS
  PASS: radiation_SN_angular_moment_validation

Results: 7 passed, 1 failed


0% tests passed, 3 tests failed out of 3

Label Time Summary:
phase13       =   0.02 sec*proc (3 tests)
validation    =   0.02 sec*proc (3 tests)

Total Test time (real) =   0.03 sec

The following tests FAILED:
    161 - test_level_c_coupled_verification (Failed)        phase13 validation
    162 - test_cht_validation (Failed)                      phase13 validation
    173 - test_thermal_challenge_vv (Failed)                phase13 validation
Errors while running CTest
F.6 - Selection CTest finale R2
Internal ctest changing into directory: /mnt/data/cfdx_md_r2/build
Test project /mnt/data/cfdx_md_r2/build
      Start   1: test_audit_numerics
 1/10 Test   #1: test_audit_numerics .................   Passed    0.00 sec
      Start   2: test_audit_flux_contract
 2/10 Test   #2: test_audit_flux_contract ............   Passed    0.00 sec
      Start   3: test_audit_cht_time
 3/10 Test   #3: test_audit_cht_time .................   Passed    0.01 sec
      Start   4: test_audit_cht_contract
 4/10 Test   #4: test_audit_cht_contract .............   Passed    0.01 sec
      Start   5: test_audit_energy_boundary_metric
 5/10 Test   #5: test_audit_energy_boundary_metric ...   Passed    0.00 sec
      Start 161: test_level_c_coupled_verification
 6/10 Test #161: test_level_c_coupled_verification ...   Passed    0.01 sec
      Start 162: test_cht_validation
 7/10 Test #162: test_cht_validation .................   Passed    0.00 sec
      Start 163: test_thermal_vv
 8/10 Test #163: test_thermal_vv .....................   Passed    0.00 sec
      Start 164: test_thermal_zones_vv
 9/10 Test #164: test_thermal_zones_vv ...............   Passed    0.00 sec
      Start 173: test_thermal_challenge_vv
10/10 Test #173: test_thermal_challenge_vv ...........   Passed    0.02 sec

100% tests passed, 0 tests failed out of 10

Label Time Summary:
audit            =   0.03 sec*proc (5 tests)
phase13          =   0.04 sec*proc (5 tests)
qualification    =   0.00 sec*proc (1 test)
regression       =   0.03 sec*proc (5 tests)
thermal          =   0.01 sec*proc (2 tests)
validation       =   0.04 sec*proc (5 tests)

Total Test time (real) =   0.07 sec
F.7 - Reexecution des suites anterieures et sanitizers
test_audit_numerics: exit=0 ; SUMMARY	checks=87	failures=0
test_audit_flux_contract: exit=0 ; SUMMARY	checks=56	failures=0
test_audit_cht_time: exit=0 ; SUMMARY	checks=24	failures=0
test_audit_cht_contract (ASan/UBSan): exit=0 ; SUMMARY checks=101 failures=0 max_T_error=5.8093974075745791e-11 max_energy_normalized=1.7507773009128869e-11
test_audit_energy_boundary_metric (ASan/UBSan): exit=0 ; SUMMARY checks=12 failures=0
F.8 - Verification de livraison
Le script embarque a ete extrait et execute sur le ZIP original : application du patch reussie (code 0). Les 1181 fichiers reconstruits a partir de ce Markdown ont ete compares par SHA-256 a la copie effectivement testee : 1181/1181 identiques. Les 1174 fichiers de la copie originale restent identiques au contenu du ZIP.
Le script refuse effectivement une destination existante et une archive source differente. Le code Python embarque correspond octet pour octet au script verifie, les blocs Markdown sont fermes, et le registre a ete recompte : 53 entrees, dont 34 corrigees/testees, 2 partielles et 17 ouvertes.
La reconstruction de livraison est une verification d'application et d'identite des sources ; elle n'est pas comptee comme une seconde execution de la campagne CTest. Un avertissement git apply relatif a une ligne blanche finale heritee de V1 a ete conserve ; il ne bloque pas l'application ni l'identite des fichiers.