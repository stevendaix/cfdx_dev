# FICHE MÉTHODE 17 — Validation indépendante de puissance momentum

**Statut :** DRAFT — identité de comptabilité discrète ; aucune preuve de stabilité ni qualification physique.
**Périmètre :** oracle de puissance sur forces intégrées, tests déterministes et raccordement futur aux opérateurs de production.

## 1. Objectif et domaine de validité

L'objectif est de contrôler indépendamment le bilan de puissance mécanique associé à un ensemble de forces de faces et de sources cellulaires fournis par l'appelant. Le contrôle compare deux assemblages algébriquement équivalents mais calculés par des chemins distincts : un oracle basé sur les sauts de vitesse aux faces et un scatter signé des forces vers les cellules.

Cette fiche établit une identité de comptabilité finie-dimensionnelle. Elle ne constitue pas une preuve de positivité de la diffusion complète, de stabilité du solveur, de convergence non linéaire, de réalisabilité turbulente ou de précision physique.

## 2. Formulation continue

Pour un champ de vitesse incompressible, le produit de l'équation de quantité de mouvement par la vitesse donne un bilan d'énergie cinétique. Le travail des contraintes apparaît comme un flux aux frontières du volume de contrôle et comme un échange volumique associé au gradient de vitesse. La convection et la pression transportent également de l'énergie ; leur contribution locale ou globale dépend des conditions aux limites et de la discrétisation.

Il est donc incorrect d'imposer sans hypothèses que chaque canal de puissance soit dissipatif. En particulier, la convection, la pression, les efforts imposés aux frontières et les forces volumiques peuvent fournir ou retirer de l'énergie au domaine.

Référence de contexte sur les propriétés énergétiques des discrétisations incompressibles : Morinishi, Lund, Vasilyev et Moin, *Fully conservative higher order finite difference schemes for incompressible flow*, Journal of Computational Physics (1998), DOI: 10.1006/jcph.1997.5882.

## 3. Discrétisation Volume Fini

Pour une face interne orientée de son propriétaire P vers son voisin N, une force intégrée unique \(\mathbf F_f\) [N] contribue à la puissance par :
\[
P_f=(\mathbf U_P-\mathbf U_N)\cdot\mathbf F_f.
\]
Pour une face frontière :
\[
P_f=\mathbf U_P\cdot\mathbf F_f.
\]
Une source cellulaire intégrée \(\mathbf S_c\) [N] contribue par \(\mathbf U_c\cdot\mathbf S_c\). Avec \(\mathbf U\) en [m/s], la puissance est en [W].

Le scatter signé indépendant applique \(+\mathbf F_f\) au propriétaire et \(-\mathbf F_f\) au voisin pour une face interne, et seulement \(+\mathbf F_f\) au propriétaire pour une face frontière. Le produit scalaire des forces cellulaires résultantes avec les vitesses doit retrouver l'oracle face, à l'arrondi flottant près.

Chaque face interne doit apparaître une seule fois dans la liste des faces orientées. Le voisin d'une face frontière vaut exclusivement `FaceOwnership::BOUNDARY` (-1) ; toute autre sentinelle négative est invalide. Les forces et sources sont déjà intégrées [N] : aucun facteur de volume ou d'aire ne doit être appliqué une seconde fois.

Les canaux convection, diffusion moléculaire, traction visqueuse transposée, pression, force extérieure et accumulation sont séparés. Les puissances internes et frontières sont rapportées distinctement. La comparaison oracle/scatter est effectuée opérateur par opérateur afin qu'une compensation entre canaux ne masque pas une erreur de signe.

## 4. Traitement des frontières

Les faces frontières ne possèdent pas de cellule voisine dans le graphe de puissance. Leur contribution utilise la vitesse de la cellule propriétaire et doit rester explicite dans le rapport. Le travail frontière n'est pas supposé nul : il dépend du flux de quantité de mouvement et de la condition physique imposée.

Le test de positivité actuellement présent est volontairement plus restreint : graphe à conductances positives et frontière de Dirichlet homogène. Dans ce cas seulement, le travail du Laplacien discret s'écrit comme une somme de termes \(g_f\lVert\Delta\mathbf U_f\rVert^2\), non négative. Cette propriété ne se transpose pas automatiquement à des conditions frontières générales ni au flux visqueux de production.

## 5. Pièges connus

- Confondre force intégrée [N] et densité de force ; appliquer un volume ou une aire une seconde fois.
- Compter deux fois une face interne ou inverser le signe sans inverser simultanément l'orientation et la force.
- Laisser une sentinelle négative inconnue être interprétée comme une frontière.
- Comparer seulement la puissance totale : une erreur de signe sur un canal peut être masquée par un autre canal.
- Déduire une dissipation positive de la seule identité de comptabilité.
- Extrapoler le résultat du Laplacien de graphe à viscosité constante au schéma complet non orthogonal, à viscosité variable ou avec correction différée.
- Assimiler une identité exacte en arithmétique réelle à une preuve sur toutes les exécutions flottantes.
- Déduire stabilité temporelle, convergence du couplage pression-vitesse, réalisabilité des contraintes de Reynolds ou précision d'une fermeture turbulente d'un bilan de puissance fermé.

## 6. Protocole de vérification et essais

Les tests déterministes couvrent un maillage à deux cellules avec faces internes et frontières :
- comparaison oracle face / scatter signé pour convection, viscosité, viscosité transposée, pression, force extérieure et accumulation séparément ;
- séparation des contributions internes et frontières pour chaque opérateur face ;
- invariance au renversement de l'orientation d'une face interne si la force est inversée simultanément ;
- rejet des tailles incohérentes, propriétaire/voisin invalide, sentinelle négative inconnue, NaN et overflow ;
- test de signe limité au Laplacien de graphe à conductances positives avec Dirichlet homogène.

Critère : égalité des deux chemins pour chaque canal et chaque composante, dans une tolérance numérique justifiée par les opérations flottantes du test. Aucun seuil ne doit être élargi pour faire passer un résultat. Les erreurs de topologie et les valeurs non finies doivent échouer explicitement.

## 7. Statut d'implémentation

L'oracle est un vérificateur de contrat pour des forces intégrées fournies par son appelant. La capture des forces brutes directement depuis les opérateurs de production à l'état accepté n'est pas encore raccordée. Cette PR ne valide donc pas encore les opérateurs réellement assemblés sur les cas de production.

Audit du point de raccordement actuel (lecture de `steady_incompressible_solver.h` et `diagnose_momentum_face_budget`) :
- le solveur reconstruit `full_flux` et le flux convectif avec le même `mass_flux` autoritaire ;
- le flux diffusif est candidat à être isolé par différence `full_flux - convection_flux`, sous réserve de confirmer le périmètre exact des autres contributions ;
- la traction visqueuse transposée est fournie séparément par `audit_transpose_flux`, avec contribution résiduelle de signe opposé ;
- les sources pression, force extérieure et accumulation sont intégrées par cellule dans le diagnostic existant.

Le contrat de raccordement envisagé est de collecter, pour chaque composante, les forces face réellement utilisées par l'assemblage et les sources cellulaires correspondantes, toutes issues du même état accepté et du même flux massique. Les trois composantes doivent être réunies avant l'appel à l'oracle. Chaque canal devra ensuite être réconcilié avec le budget cellulaire et le résidu CSR signé, cellule par cellule. Les signes candidats devront être confirmés sur les opérateurs effectifs avant toute implémentation.

La suite est donc : (1) raccorder les forces de production sans reconstruire une seconde physique, (2) tester l'équivalence canal par canal sur un cas manufacturé, (3) comparer la somme cellulaire avec le résidu CSR signé, puis (4) mener des études de convergence maillage/temps et des validations a-posteriori sur des références indépendantes de la calibration.

La relation à Lean reste limitée aux contrats abstraits owner/neighbour, incidence et travail. Les tests C++ vérifient des maillages runtime déterministes ; la correspondance générale entre topologie Lean, `Mesh::ownership()` et sorties des opérateurs de production reste à formaliser.
