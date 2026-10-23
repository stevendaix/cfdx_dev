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
