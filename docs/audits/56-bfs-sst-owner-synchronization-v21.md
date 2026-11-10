# CFDX V21 — synchronisation de la contrainte omegaWallFunction et diagnostic BFS

**Date** : 10 novembre 2026. **Base de code** : sources complètes CFDX V20 conservées dans `CFDX_V20_sources.zip` (elles-mêmes obtenues à partir de l'archive d'origine du projet). **Périmètre** : solveur SST de production, BFS diagnostic 712 cellules, sans modification des seuils et sans push distant.

## Résultat établi, sans ambiguïté

1. Le mode expérimental `sst_owner_consistent_projection_v21_experimental=true` corrige le défaut entre la valeur de omega imposée dans les 69 cellules propriétaires et la cible `omegaWallFunction` recalculée après la modification de `k`.
2. Dans le BFS de 712 cellules, le défaut maximal de la contrainte propriétaire passe de **396,4000173842** sur V19 à **0** sur le deuxième état extérieur V21. Le nombre de cellules propriétaires est inchangé (69).
3. **La convergence du deuxième état reste REFUSÉE à 800 itérations SST** : `k_returned_nonlinear_raw=0.008735979683`, `omega_returned_nonlinear_raw=24.41970936`, `k_returned_row_scaled=0.548989493`, `omega_returned_row_scaled=1451.800745`. La cellule 357 porte désormais un défaut de `k` plutôt qu'un défaut de ligne propriétaire `omega`; la cellule 363 apparaît au maximum du défaut local relatif de `omega`.
4. Le mode Anderson V20 ajouté au correctif V21 n'est **pas** une solution suffisante : 18 itérations SST, puis arrêt fail-closed au premier état extérieur, sans état extérieur accepté dans le journal observé.
5. Le mode **désactivé** reproduit le premier état du journal V19 à l'identique sur le binaire de production V21, avec 217 itérations SST, `U=0.5251832153953`, continuité `0.001771416233814`, `k_local=1.805716607269e-13` et `omega_local=1.994937388111e-13`.

Ces faits **ne** constituent **ni** une convergence globale BFS, **ni** une validation NASA, **ni** une qualification des six Schur/PCD/AMG/MGR.

## Cause physique et mathématique traitée

Dans les cellules propriétaires de la condition `omegaWallFunction`, l'assemblage remplace la ligne FV ordinaire par la contrainte :

`F_wall,i(k,omega,U) = omega_i - omega_target,i(k,U,nu_t,...) = 0`.

La valeur cible dépend du champ `k` et de l'état de turbulence. Dans le chemin antérieur, la ligne était assemblée avant l'évolution du champ `k`, mais sa fermeture résiduelle était testée contre la cible actualisée. Le résidu algébrique du système gelé ne mesurait donc pas le défaut de la contrainte sur l'état accepté.

L'option V21 procède transactionnellement pour chaque candidat de transport SST :

1. Évalue le triplet `(omega_face, G_owner, omega_owner_target)` avec le champ candidat `k,omega`.
2. Projette `omega_i = omega_owner_target_i` **uniquement sur les lignes propriétaires**.
3. Réévalue les valeurs de face, la production `G` et les cibles propriétaires sur **le même état projeté**.
4. Si la propriété d'owner a changé, si une valeur est non finie/non positive, ou si la cible dépend encore de la valeur projetée au-delà de `64 * epsilon_machine`, rejette sans publier de champ, `G` ni face.
5. Après l'étape Picard ou la recherche linéaire Anderson, réévalue le vrai opérateur SST non linéaire avant de tester les conditions de convergence déjà présentes. Les modes V20 Anderson et line-search restent facultatifs ; aucun critère n'est atténué.

La projection n'annule **pas** le défaut des équations de transport des cellules non contraintes. Une projection compatible avec la contrainte propriétaire n'est pas un pas Newton complet, ni une preuve de diminution du mérite global si elle est utilisée avec Picard non globalisé.

## Bilan réellement réassemblé à la 100e itération SST du deuxième état

Le second lancement V21 avec `CFDX_SST_ASSEMBLY_TRACE_V20=1` retrouve pour la ligne `k`, cellule 357 (unités du bilan discret en cohérence avec l'assemblage FV) :

| Contribution signée | Valeur |
| --- | ---: |
| Convection | `+0.000365748847211` |
| Diffusion | `+0.006765443061181` |
| Production | `-0.1412872344173` |
| Destruction | `+0.1254200821097` |
| Défaut de ligne | `-0.008735960399175` |
| Erreur de recomposition | `-1.387778780781e-17` |

La ligne `omega` propriétaire 357 a `constrained=1`, convection/diffusion/production/destruction `0`, `balance=0`. Pour la ligne **non propriétaire** `omega`, cellule 363, le même instant donne :

| Contribution signée | Valeur |
| --- | ---: |
| Convection | `-0.5500659747829` |
| Diffusion | `-79.86664537065` |
| Production | `-208.2781357151` |
| Destruction | `+313.1145644987` |
| Diffusion croisée | `-0.00007704471761874` |
| Défaut de ligne | `+24.41964039343` |
| Erreur de recomposition | `+2.486899575160e-14` |

**Interprétation** : le défaut restant est bien porté par les **équations de transport et leur couplage production/destruction/diffusion**, et non uniquement par la contrainte `omegaWallFunction` gelée. La correction locale ne peut suffire ; il faut un opérateur résiduel joint `R_k, R_omega` et une direction Newton–Krylov/préconditionnement ou une globalisation montrant une réelle contraction. Les valeurs ci-dessus sont des résidus signés, pas des indicateurs d'erreur physique indépendante.

## Options et diagnostics

- Mode de production historique inchangé : `sst_owner_consistent_projection_v21_experimental=false` (par défaut).
- Synchronisation opt-in : `sst_owner_consistent_projection_v21_experimental=true` dans `numerics.raw_settings`, avec callback actif V17.
- Anderson reste **distinct** et opt-in : `sst_anderson_m2_v20_experimental=true` ; l'essai conjoint V21 n'a pas convergé.
- `CFDX_SST_DIAGNOSTICS=1` émet les comptes de propriétaires, écarts de cible, résidus réassemblés et tests d'acceptation.
- `CFDX_SST_ASSEMBLY_TRACE_V20=1` (maillages ≤4096 cellules) audite en plus les lignes 77, 357, **363**, 369, 426 et les deux pires cellules selon le défaut local relatif.
- `scripts/diagnose_bfs_sst_v21.py --log journal.log --output diagnostic.json` génère une synthèse **non qualifiante**. La stabilité de la projection rapportée porte seulement sur les échantillons observés.

## Contrôles exécutés

- CTest ciblé : 4/4 OK (`test_sst_wall_owner_sync_v21`, `test_sst_wall_owner_v18`, `test_sst_equilibrated_v19`, `test_sst_anderson_v20`). Le premier contient quatre cas de transaction et de cohérence.
- Pytest ciblé : 18/18 OK, dont trois tests V21 pour la fermeture/validation fail-closed des journaux.
- Compilation et lien du **solveur de production** : OK, build diagnostique `-O0`, non représentatif d'un benchmark Release.
- Rejeu réel HDF5 712 cellules, mode désactivé : 217 itérations SST jusqu'au premier état, identique à V19 ; retour du programme `1` car l'outer BFS ne converge pas en une itération.
- Rejeu réel HDF5 712 cellules, projection activée : premier état SST convergé en 217 itérations ; deuxième état SST rejeté à 800 itérations, rollback, seul un état externe accepté. Une campagne à trois états convergés n'existe pas.
- Rejeu projection + Anderson : pas d'état externe accepté dans le journal ; stagnation fail-closed après 18 itérations SST au premier essai.
- Patch cumulatif depuis l'archive originale : `git apply --check` et application sur la reconstruction indépendante de l'original ; empreintes SHA-256 identiques des **1777 fichiers sources** (les fichiers temporaires `.pytest_cache` sont exclus).

## Suites de qualification impératives

1. Diagnostiquer localement les cellules 357/363 et les sources `Pk`, `Dk`, `Pomega`, `Domega`, la diffusion croisée, `G`, `nu_t`, F1/F2 **sur deux itérations SST successives**, pour détecter la dynamique oscillante.
2. Raccorder à la résolution de production un `F(k,omega)` joint, différencié directionnellement à **conditions de vitesse/pression fixées**, y compris les dépendances des cibles propriétaires de `k` ; comparer JFNK et Picard stabilisé avec évaluation réassemblée de `F`, sécurité positivity et rollback.
3. Après fermeture SST seulement, valider le résidu global indépendant `U,p,k,omega`, le bilan masse et la fermeture de `nu_t`, sans modifier les critères. Rejouer la comparaison avec **les six Schur/PCD/AMG/MGR** à budget et matériel égaux.
4. Répéter sur le BFS 20 540 cellules et sur les cas NACA et de cavité avant promotion du changement vers les paramètres par défaut. La validation Driver–Seegmiller requiert encore une solution globalement convergée et les données physiques expérimentales comparables.

**Limites de la présente livraison** : aucun test CI complet, pas de test GPU/MPI, pas de run BFS 20540 finalisé, pas de push GitHub. L'équivalence de `updateCoeffs()` et `manipulateMatrix()` OpenFOAM n'est pas démontrée par la seule fermeture de la ligne propriétaire.
