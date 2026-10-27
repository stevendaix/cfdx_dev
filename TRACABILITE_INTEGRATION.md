# CFDX — Matrice de traçabilité de l’intégration combinée

Date de préparation : 2026-10-27  
Branche d’intégration : `integration/combined-tree-2026-10-27`  
Base déclarée : `master` au démarrage de l’étape 1. Le SHA exact de la base doit être consigné dans le rapport d’exécution et vérifié avant toute fusion.

## Règles de preuve

- Une PR ouverte ou en brouillon n’est pas considérée comme intégrée.
- Un SHA de merge affiché pour une PR non fusionnée ne constitue pas une preuve d’intégration dans `master`.
- La compilation et les tests doivent être exécutés sur le SHA combiné exact.
- Aucun lot ne peut être marqué « clôturé » avant validation de ses critères numériques et logiciels.
- Les branches `patch-r4` et `lot-thermique-v2` sont des hypothèses du plan : leur existence, leur propriétaire et leur contenu doivent être vérifiés avant l’étape 3.
- Les corrections indépendantes (Krylov, F06, viscosité variable, SST) doivent conserver leurs tests et contrats sans résolution de conflit par simple préférence de version.

## Matrice initiale

| ID / lot | Source prévue | Sujet | État constaté au démarrage | Gate avant intégration |
|---|---|---|---|---|
| BASE | `master` | Base de travail | SHA à confirmer dans le relevé d’exécution ; le commit récent annoncé porte sur la certification GMRES et la qualification MPI | Vérifier HEAD et absence de modifications locales avant tout reset local |
| F03 / #1078 | PR #1078 | Couplage pression-vitesse Rhie–Chow / Schur | Non qualifié par cette matrice ; état et HEAD à relire | Contrat `H/A = U + R*G*p`, continuité cohérente, anti-damier et tests physiques sur le SHA combiné |
| F06 / #1129 | `fix/f06-nonorthogonal-transport-contract` | Transport non orthogonal, métriques de résidu | PR ouverte, brouillon, HEAD `dc319f916e5297fcacb259e0f9a5e9ed617cc532`; mergeabilité déclarée false au relevé | Résoudre la non-mergeabilité, puis CI sur le HEAD exact ; tests transport, énergie, PCD/N9 et maillages skewés |
| F11 / #1139 | `fix/skewed-variable-viscosity-face-reconstruction` | Viscosité variable sur faces skewées | PR ouverte, brouillon, HEAD `8a419c90484b555631c885ae9bdbbd07b3fd2b25`; CI complète à confirmer | Stress affine, MMS, Release, DebugSanitizers et maturité numérique sur HEAD exact |
| SST / #1140 | `fix/sst-persisted-cli-restart` | Reprise SST depuis cas persistant | PR ouverte, brouillon, HEAD `8aec53f166c7df057ba1e7c5d2143d074d05496f`; CI exacte à confirmer | Scénarios CLI continus/interrompus/repris, champs U/p/k/omega, CTest exécutant réellement pytest, Release et sanitizers |
| R4 | Branche `patch-r4` (à vérifier) | Infrastructure solveur, Krylov et gardes robustesse | Existence/HEAD non encore établis | Inventorier les commits et diffs ; préserver les correctifs GMRES récents de master ; tests indépendants et certification FP64 |
| Thermique v2 | Branche `lot-thermique-v2` (à vérifier) | Physique thermique, rayonnement, turbulence et raccordements | Existence/HEAD non encore établis | Vérifier provenance et base ; audit des fichiers communs ; ne pas écraser les gardes et diagnostics existants |
| Lot 0 | À déterminer après inventaire | Provenance et harnais | Non clôturé | Matrice des sources, dépendances, SHA et tests |
| Lot 1 | À déterminer après inventaire | Correctifs transversaux / robustesse | Non clôturé | Tests unitaires, intégration, Release et DebugSanitizers |
| Lot 2 | #1078 et dépendances associées | Contrat pression-vitesse | Non clôturé | Gates F03 et qualification physique sur l’arbre combiné |

## Journal d’exécution

| Étape | Résultat | Preuve |
|---|---|---|
| 1 — création de branche et matrice initiale | En cours ; le commit du fichier doit être vérifié | Branche GitHub et commit de création à renseigner après écriture |
| 2 — correctifs critiques | Non commencé | Ne pas fusionner de PR en brouillon sans revue des checks et de la mergeabilité |
| 3 — R4 / v2 | Bloqué jusqu’à vérification des branches sources et des bases | Diff d’inventaire à produire |
| 4 — F03 | Non commencé | Logs CTest et résultats numériques à joindre |
| 5 — validation finale et audit | Non commencé | Rapport CTest exact-SHA, compilation et hashes réels requis |

## Décision de clôture

À ce stade, aucun lot n’est déclaré « fusionné et testé ». Cette matrice est un registre de suivi, pas une preuve de compilation ou de validation numérique.
