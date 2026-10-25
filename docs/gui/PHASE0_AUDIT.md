# Phase 0 — Confrontation du code GUI à la feuille de route

**Référentiel** : feuille de route « Scientific Workbench nouvelle génération » adoptée comme référence produit (8 phases, 5 lots).
**Base auditée** : `master` `e29c9a64` (2026-10-23) — `workbench.py` (959 lignes), `gui_3d.py` (253), `run_center_panel.py` (80), `visualization_pipeline.py`, `setup_panel.py`, `probes_panel.py`, `mesh_browser_panel.py`, `results_panel.py`, `theme.py`.
**Objectif** : chiffrer l'écart réel avant tout prototypage, décider quoi réutiliser, cadrer le prototype de référence du Lot 1.

---

## 1. Conformités déjà réelles (à réutiliser, pas à reconstruire)

| Acquis | Localisation | Phase couverte | Décision |
|---|---|---|---|
| Docks masquables + persistance QSettings | `workbench.py` | Phase 0 / Lot 1 | Réutiliser tel quel |
| Workspaces setup/mesh/run/analyze (#1036) | `workbench.py` | Phase 0 / Lot 1 | Réutiliser, étendre |
| Palette de commandes Ctrl+K | `workbench.py` (`findChildren(QAction)`) | Phase 6 | Réutiliser le geste, remplacer le mécanisme par un registre de capacités (Lot 2) |
| Pipeline non destructif append/remove/clear/replay (#1041/#1044/#1046) | `visualization_pipeline.py` | Phase 3 / Lot 3 | Réutiliser le cœur ; rendre les filtres éditables |
| Sélection de patch 3D (opacité) | `gui_3d.py` | Phase 2 / Lot 2 | Étendre vers sélection faces/régions |
| Arbre workflow, états [OK]/[!]/[X]/[~], `WorkflowStatus.CONDITIONAL` (#959) | `workbench.py` | **Phase 7 — partiellement déjà là** | Assumer l'existant, pas de refonte |
| État applicatif central Application/events (signaux queued) | `workbench.py` | Socle transverse | Conserver comme ancre du prototype |

**Conséquence** : la Phase 0 ne part pas de zéro et ne doit **pas** créer d'architecture parallèle. Le prototype s'ancre sur `workbench.py` / `gui_3d.py` existants.

## 2. Gaps mesurés, mappés aux lots de la feuille de route

| # | Gap (mesuré dans le code) | Lot cible | Ampleur |
|---|---|---|---|
| G1 | Pipeline modal : contour/slice/glyph via `QInputDialog` à la création ; pas de re-paramétrage, pas de réordonnancement, pas de threshold/clip ; normale de slice fixe | Lot 3 | Majeur |
| G2 | Inspecteur Properties générique (texte + Apply par propriété) ; sans unités, contraintes ni validation typée | Lot 2 | Majeur |
| G3 | Actions désactivées sans explication de prérequis ; palette = `findChildren(QAction)`, pas de registre de capacités formel | Lot 2 | Moyen |
| G4 | Aucun undo/redo (0 occurrence) | Lot 2 | Majeur |
| G5 | Run Center minimal (~80 lignes) : texte brut, aucune courbe, aucun parsing structuré — alors que le solveur produit `to_string(convergence_status)` (STAGNATED/DIVERGED/MAX_ITERATIONS) et la trace INCOMPRESSIBLE_ITER à 12 champs (#910) | Lot 4 | Moyen (données solveur déjà prêtes) |
| G6 | Pas de comparaison multi-vues synchronisée | Lot 4 | Majeur |
| G7 | Pas de sélection fine faces/régions | Lot 2 | Moyen |

**Rejet explicite** (décision produit 2026-10-23) : G5/Run Center ne remonte pas dans le Lot 1 malgré sa rentabilité — l'ordre des lots est non négociable. Il faut d'abord une expérience cohérente centrée viewport.

**Mise à jour G3 (Lot 2 — palette pilotée par le registre, 2026-10-23)** : le registre de capacités couvre désormais toutes les commandes du Workbench (actions fichier, Check/Run/Stop, contour/slice/glyph, workspaces — 16 entrées, labels « X workspace » pour désambiguïser). La palette Ctrl+K énumère le registre — plus aucun `findChildren(QAction)` pour l'énumération — et évalue les prérequis en direct à chaque rafraîchissement ; `_refresh_action_availability` est piloté par la même évaluation. G3 est clos (Lot 1 + Lot 2). Les actions DAT gardent leur garde au déclenchement ; un prérequis « case ouvert » reste envisageable en suivi.

**Mise à jour G7 (Lot 2 — sélection fine faces/régions, 2026-10-23)** : le browser de maillage devient un arbre régions → patches → faces (les faces sont créées à l'expansion du nœud patch, jamais au chargement) et émet des `MeshSelection` dont les `stable_id` sont les ids canoniques partagés — `region:<nom>`, `[region:<r>/]patch:<nom>`, `<patch>/face:<k>` — produits et parsés strictement par le nouveau module `cfdx.selection.py`. La vue 3D interprète les mêmes ids : patch = opacité (inchangé), région = opacité par région, face = acteur dédié `workbench.selection.face` avec atténuation 0.2 des autres acteurs ; la barre de statut du Workbench affiche le libellé lisible de la sélection courante (`selection_summary()`). G7 est clos (browser + vue 3D). Hors périmètre : chargement 3D multi-régions (la vue 3D charge toujours un seul jeu de points ; une région sans acteur correspondant est ignorée silencieusement) et picking souris direct dans le viewport (P1).

**Mise à jour G2 (Lot 2 — inspecteur contextuel typé, 2026-10-24)** : le panneau Properties devient un inspecteur typé piloté par les contrats `PropertyState` existants — chaque ligne déclare son `ParameterType` et reçoit l'éditeur correspondant (réel → `QDoubleSpinBox`, entier → `QSpinBox`, booléen → `QCheckBox`, texte → `QLineEdit`), le libellé porte l'unité `[unit]` quand le contrat en fournit une, et l'impact `REQUIRES_RESTART` est annoncé en infobulle. Le nouveau module `cfdx/properties_panel.py` est un adaptateur Qt fin ; le cast ad hoc `float(...)` de `numerics.cfl` disparaît (l'éditeur fournit la valeur déjà typée), la validation passe par le contrat `set_property` et les erreurs s'affichent en clair dans la barre de statut. L'inspecteur est contextuel sur **toute** sélection : workflow et maillage convergent vers `Application.select` + événement `SelectionChanged`, et les sélections maillage résolvent des propriétés typées read-only (région : type objet/cells/faces/points/patches ; patch : type objet/type de patch/faces/région ; face : type objet/index/patch/région) depuis le `MeshCatalog` conservé par `Application.load_mesh` — la topologie reste possédée par C++/HDF5, jamais éditable. Hors périmètre : bornes min/max des réels (le contrat de domaine ne fournit que la validation « réel fini ») et édition des propriétés maillage. G2 est clos.

**Mise à jour G4 (Lot 2 — undo/redo des éditions applicatives, 2026-10-24)** : la couche application enregistre elle-même ses éditions réversibles — le nouveau module `cfdx/application/history.py` porte un historique borné (`EditHistory`, deux piles, profondeur maximale `MAX_UNDO_DEPTH` = 100 ; un nouvel enregistrement purge la branche redo) d'entrées `HistoryEntry` (libellé en clair + callables inverse/rejeu). `Application.set_property` capture la valeur précédente via le nouveau contrat symétrique `get_property` avant d'exécuter, n'enregistre que les changements réels, et republie l'état après enregistrement pour que l'activation des actions Undo/Redo passe par le chemin habituel `ApplicationStateChanged` → `_refresh_action_availability`. Les actions Undo (Ctrl+Z) et Redo (Ctrl+Shift+Z, Ctrl+Y) entrent au registre des capacités avec les prérequis « an edit to undo/redo » — indisponibles donc expliquées en clair dans la palette et l'infobulle, jamais silencieusement désactivées (principe G3) ; la barre de statut rapporte `Undid: …`/`Redid: …` ou l'échec en texte clair. Une nouvelle session (`replace_session`, donc `new_project`/`open_project`) vide l'historique : des entrées d'un case précédent seraient incohérentes. Hors périmètre : éditions du panneau setup spécialisé et sondes (non réversibles) ; l'undo restaure les valeurs typées, pas le marqueur restart-required (comportement conservateur). G4 est clos — le Lot 2 (G2/G3/G4/G7) est complet.

## 3. Décisions de réutilisation vs prototypage

**Réutiliser sans réécrire** : état applicatif central, docks/workspaces, cœur du pipeline non destructif, arbre workflow, sélecteur de patch.

**Prototyper dans le Lot 1 (nouveau)** : shell de référence et viewport central — agencement, densité d'information scientifique, navigation, comportement des docks, valeurs de thème déjà validées (palette, contrastes — cf. SCIENTIFIC_DESIGN_SYSTEM.md).

**Interdit pendant la Phase 0** : architecture parallèle, réécriture totale du GUI d'un coup, nouvelle passe de thème « opportuniste », tout changement de tolérance ou de gate depuis le GUI.

## 4. Périmètre du prototype de référence (Lot 1)

1. Audit chiffré du shell actuel (ce document ; à maintenir à jour à chaque PR fusionnée).
2. Maquette interactive validée par Steven avant tout code.
3. Prototype incrémental sur l'existant : un seul `workbench` cible, viewport 3D au centre, docks repositionnés, palette et arbre conservés.
4. Validation visuelle : captures avant/après, tests Qt headless sur le comportement (activation/désactivation, persistance, navigation), preuve par PR.

## 5. Critères de sortie de la Phase 0

- Document d'audit figé dans `docs/gui/` (ce document).
- Maquette interactive validée (go/no-go explicite de Steven).
- Prototype du shell + viewport manipulable : viewport central, docks masquables persistés, palette Ctrl+K fonctionnelle, aucune perte de fonctionnalité existante.
- Tests Qt headless verts sur le comportement du prototype ; captures avant/après produites.
- Décision explicite de Steven (go/no-go) avant d'ouvrir le Lot 2.

## 6. Gardes de validation (rappel des conventions repo)

- Une validation GUI réussie n'est **pas** une validation numérique ; INCONCLUSIVE ≠ PASS.
- CI verte par exécution sur le head exact ; pas de merge sans feu vert explicite de Steven (squash, `expectedHeadSha` épinglé).
- Les workflows validation-total ne tournent sur une PR que via le label `validation-total` (les skipped sont normaux).
- Jamais de modification de tolérances/gates depuis le GUI ; sources jamais détruites ; une action = une implémentation.

---

*Document de référence Lot 1 — étape 1 de la Phase 0. Toute PR GUI doit se rattacher à un lot et à un gap identifié ici.*
