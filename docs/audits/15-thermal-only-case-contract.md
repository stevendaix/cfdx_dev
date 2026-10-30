# FICHE MÉTHODE 15 — Contrat de données thermique autonome (A1)

**En-tête d'évidence (gouvernance §16.5)**  
- Scope : contrat typé de cas thermique autonome, séparation matériaux / régions solides et fluides / nœuds fluides 0D / interfaces / liens conductifs. Cette fiche spécifie les sémantiques et invariants de données ; elle ne qualifie aucun opérateur numérique ni solveur.
- Statut : IMPLÉMENTÉ — A1 (contrat en mémoire + validation fail-closed + tests unitaires enregistrés). CI exacte HEAD en cours ; aucune qualification physique revendiquée.
- Issue / PR : #1150 / #1151.
- Spécification détaillée : \`docs/thermal/thermal-case-data-contract.md\`.
- Sources code : \`src/cfdx/physics/thermal_case.h\`, \`src/cfdx/physics/thermophysical_models.h\`.
- Limites : pas de sérialisation HDF5, pas de résolution mesh/patch/champ, pas d'assemblage, pas de driver thermique, pas de qualification physique.
- Convention : SI partout ; température absolue en K ; pression absolue en Pa ; aucun changement d'unité implicite.

## 1. Références sources

- S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980 — DOI: 10.1201/9781482234213. Référence de base pour les bilans thermiques et volumes finis.
- F. Moukalled, L. Mangani, M. Darwish, *The Finite Volume Method in Computational Fluid Dynamics*, Springer, 2016 — ISBN: 978-3-319-16873-9.
- CFDX, fiches méthodes existantes : \`docs/audits/10-thermique-cht.md\`, \`docs/audits/13-conduction-solide-multimateriaux.md\`, \`docs/audits/14-energie-fluide-production.md\`.

## 2. Formulation continue

Le contrat de données prépare les formulations thermiques suivantes, sans les assembler dans A1.

Pour une région maillée, la forme conservative de l'énergie est :
$$
\frac{\partial(\rho h)}{\partial t}
+ \nabla\cdot(\dot{m}_v h)
= \nabla\cdot(k\nabla T)+q'''
$$
où $h(T)=\int c_p(T)\,dT$, $\dot{m}_v=\rho\mathbf{u}$ est le flux massique volumique et $q'''$ une source volumique.

- Solide immobile : $\dot{m}_v=0$.
- Fluide stagnant : $\dot{m}_v=0$, sans équation de quantité de mouvement.
- Fluide à advection prescrite : flux massique de faces fourni/importé explicitement ; aucune reconstruction implicite depuis une vitesse cellule.
- Nœud fluide fermé 0D :
$$
\frac{d(mh)}{dt}=\sum_j G_j(T_j-T_i)+Q_i
$$
avec $G_j\ge 0$, $Q_i$ source en W et $m=\rho(T,p)V$. Pour $c_p$ constant : $mc_p\,dT/dt=\sum_jG_j(T_j-T_i)+Q_i$.

Une interface thermique entre deux domaines doit échanger des puissances égales et opposées. La convention du lien conductif est $Q_{source\to target}=G(T_{source}-T_{target})$.

## 3. Discrétisation Volume Fini

A1 ne réalise aucune discrétisation. Il fixe seulement les données nécessaires pour que les futurs opérateurs puissent conserver les bilans :

- propriétés $k(T)$, $c_p(T)$, $\rho(T,p)$ réutilisées depuis \`ThermophysicalProperties\` ;
- maillage référencé par \`mesh_id\`, sans dupliquer la topologie ;
- transport prescrit référencé par un champ massique de faces, dont l'orientation devra être validée contre la convention owner/neighbour du maillage ;
- interfaces identifiées par paires $(region, patch)$ ;
- liens 0D utilisant une conductance $G$ en W/K.

Le futur assemblage devra intégrer chaque flux interne une fois avec signes opposés, et vérifier un bilan énergétique indépendant de la matrice. Aucune précision d'ordre, stabilité ou convergence numérique n'est attribuée à A1.

## 4. Traitement des frontières

Le contrat A1 identifie les patches d'interface, mais ne définit pas encore la collection générale des conditions limites thermiques. Les conditions Dirichlet, flux imposé/Neumann, convection/Robin, puissance et adiabatique seront définies et validées avec le driver A2.

Les interfaces CHT A1 exigent deux régions distinctes et deux patches nommés. Le contrôle de l'existence des patches, de leur topologie et de l'appariement appartient à l'adaptateur mesh-aware A3 ; le validateur A1 ne prétend pas les vérifier.

## 5. Pièges connus

- Ne pas confondre type physique de région et matériau thermophysique.
- Un fluide stagnant n'est pas un fluide en écoulement gelé : le premier ne possède pas de flux massique ; le second consomme un flux fourni par un calcul d'écoulement préalable.
- La convection naturelle ne découle pas de la diffusion seule ; aucune conductivité artificiellement augmentée ne doit la remplacer.
- $R_c$ en m².K/W est une résistance surfacique ; $G$ en W/K est une conductance globale. Elles ne sont pas interchangeables sans aire.
- Les nœuds 0D A1 sont fermés : débits massiques et transport d'enthalpie de réseau sont reportés en A5.
- La masse d'un nœud dépend de $\rho(T,p)V$ ; la pression initiale est explicite pour les modèles de densité dépendant de la pression.
- Les propriétés doivent rester finies et positives dans leur plage déclarée ; une valeur par défaut de viscosité ne doit pas être exigée par un cas thermique autonome.
- Les valeurs numériques sont en SI et la température est absolue ; aucune conversion automatique ni interprétation °C/K implicite.

## 6. Ordre de précision attendu et validation exigée

A1 ne possède pas d'ordre de précision numérique : aucun solveur ni opérateur n'est livré.

Tests de contrat enregistrés dans \`test_thermal_case_contracts\` :
- cas combinant régions solides, fluide stagnant et nœud 0D ;
- indépendance vis-à-vis de la viscosité ;
- cohérence mode temporel / pas de temps ;
- flux massique de faces obligatoire pour advection prescrite ;
- rejet du transport fluide sur solide ;
- références matériau invalides ;
- volume 0D invalide ;
- appariement solide-fluide des interfaces CHT ;
- validation des liens, conductances et références ;
- rejet des énumérations inconnues ;
- pression explicite pour densité fluide dépendante de $p$.

Ces tests vérifient des invariants de données, pas la physique du solveur. Les oracles analytiques, bilans globaux, tests HDF5 round-trip et qualification solver-level restent requis en A2–A8. CI exacte HEAD doit être verte avant fusion ; aucun PASS de qualification n'est déduit de cette fiche.

## 7. Statut d'implémentation

- Modèle : \`src/cfdx/physics/thermal_case.h\`.
- Propriétés réutilisées : \`src/cfdx/physics/thermophysical_models.h\`.
- Tests : \`tests/unit/test_thermal_case_contracts.cpp\`.
- Enregistrement : \`CMakeLists.txt\`.
- Contrat de persistance cible : \`docs/thermal/thermal-case-data-contract.md\`.
- A1 implémente les types et validations de données. Les groupes HDF5 décrits sont une cible de schéma, pas une sérialisation déjà disponible.
- Aucun changement de tolérance, de solveur existant, de registre de qualification ou de résultat numérique.
- Prochaine étape : A2, adaptateur de cas et driver autonome stationnaire/transitoire pour conduction solide et diffusion fluide stagnante.
