# FICHE MÉTHODE — <Nom de la méthode>

> À copier vers `docs/audits/NN-nom-de-la-methode.md’ et à remplir **intégralement** avant la première ligne de code.
> Les sections 1–6 sont le contrat mathématique : elles ne peuvent pas rester PARTIAL.

**En-tête d'évidence (gouvernance §16.5)**
- Scope : <méthode couverte / non couverte>
- Statut : DRAFT | APPROVED | SUPERSEDED (date)
- Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json (cas) + NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json (exigences)
- Issue / PR liées : <#N>
- Dernier commit vérifié : <SHA> (<date>)
- Limitations connues : <…>

## 1. Références sources (obligatoire : DOI ou ISBN)
- <Auteur, « Titre », revue, vol., pp. année — DOI: …>
- <Livre de référence — ISBN/éditeur>

## 2. Formulation continue
EDP exacte + hypothèses (incompressible, newtonien, régime, …).
$$ <equation> $$

## 3. Discrétisation Volume Fini
Forme discrète exacte : somme sur les faces $S_f \\cdot \\phi_f$, schéma d'interpolation aux faces, traitement du terme source (explicite / linéarisé / implicite), coefficients $a_P$, $a_{nb}$.

## 4. Traitement des frontières
Dirichlet / Neumann / Robin : contribution explicite de chaque CL à la matrice $A$ et au second membre $b$ (coefficient diagonal, flux imposé, valeur limite).
Frontières spéciales : WALL, EMPTY (2D/axisym.), périodiques, sortie de pression.

## 5. Pièges connus (Gotchas)
- <ex. checkerboarding → Rhie–Chow ; jauge de pression (null-space) ; ordre 1 sur maillage non orthogonal ; …>

## 6. Ordre de précision attendu et validation exigée
- Ordre attendu : <ex. $O(h^2)$ espace, $O(\\Delta t)$ temps> — à valider par MMS sur ≥ 3 niveaux de maillage (pente log–log).
- Oracles quantitatifs : <analytique | Ghia | OpenFOAM/SU2 | …> avec tolérances L2/L∞ **pré-déclarées**.
- GCI : <rapports de raffinement, facteur de sécurité, grille applicable>.
- Conservation : <bilans à vérifier, tolérance>.
- Parallèle : <CPU vs MPI vs GPU, même maillage, tolérance 1e-12>.

## 7. Statut d'implémentation CFDX
- Fichiers source concernés : <src/…>
- Tests exécutables : <tests/…>
- Écarts connus code ↔ fiche : <… ou «à auditer»>
- Gaps et prochaines actions : <…>

## 8. Journal des modifications de la fiche
- <date> — <auteur/PR> — <création / révision>
