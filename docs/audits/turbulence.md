# Suivi actif — intégration turbulence / thermique / rayonnement CFDX

> **Fichier de suivi de référence.** Mettre à jour ce registre à chaque session de travail, commit, PR, résultat CI ou changement de verdict. Ne pas recopier aveuglément les annexes historiques : vérifier le code sur `master`, les PR ouvertes et le HEAD exact. Une case « présent » signifie code trouvé ; elle ne signifie pas qualification physique.

- **Dernière revue du registre :** 2026-10-24.
- **Dépôt :** [stevendaix/cfdx_dev](https://github.com/stevendaix/cfdx_dev).
- **Source de vérité technique :** ce fichier, y compris sa matrice de suivi et les annexes d'évidence.
- **Règles de travail :** pas de tolérance élargie, pas de tests désactivés, pas de fallback silencieux, pas d'activation générale de RANS sur la seule base de tests unitaires. Avant tout changement, inspecter `master`, le diff et les PR qui se chevauchent. Après tout changement, exécuter les contrôles applicables sur le HEAD exact, documenter les échecs et distinguer vérification logicielle, convergence numérique et validation physique.
- **État de ce commit :** mise à jour documentaire de suivi ; ne prétend pas que les chantiers ouverts sont implémentés. Les résultats de la révision 2 plus bas restent des résultats historiques du checkout et des campagnes qui y sont décrites, pas une preuve automatique du HEAD actuel.

## 1. Objet et domaine de validité

Cette fiche suit l'intégration du transport SST, de l'énergie turbulente, du rayonnement de surface à surface (S2S) et du transfert conjugué (CHT) dans CFDX. Elle distingue explicitement : (a) une équation ou un opérateur présent, (b) son raccordement au solveur de production, et (c) sa qualification par vérification numérique et validation physique. Les résultats d'un test unitaire ne suffisent pas à affirmer qu'une simulation industrielle est qualifiée.

**Statut :** intégration partielle. Les sections détaillées et la matrice TR-00 à TR-16 ci-dessous font foi pour les travaux restants ; les résultats du dossier v2 plus bas sont historiques et doivent être rattachés à leur SHA exact avant réutilisation.

## 2. Formulation continue

Pour un fluide incompressible à propriétés constantes, le système porteur comprend la continuité et la quantité de mouvement :

\[
\nabla\cdot\mathbf{U}=0,\qquad
\rho\left(\frac{\partial\mathbf{U}}{\partial t}
+\nabla\cdot(\mathbf{U}\otimes\mathbf{U})\right)
=-\nabla p^*+\nabla\cdot\boldsymbol{\tau}_{eff}+\mathbf{S}.
\]

Pour la variante SST-2003, les transports de l'énergie cinétique turbulente \(k\) et de la fréquence spécifique \(\omega\) utilisent une production limitée cohérente :

\[
\widetilde P_k=\min\!\left(\rho\nu_t S^2,\;10\beta^*\rho k\omega\right),\qquad
P_\omega=\begin{cases}\gamma\widetilde P_k/\nu_t,&\nu_t>0,\\0,&\nu_t=0.\end{cases}
\]

Le flux de chaleur turbulent est modélisé par \(\lambda_{eff}=\lambda+\rho c_p\nu_t/Pr_t\), avec \(\nu_t\ge0\) et \(Pr_t>0\). Le couplage doit consommer l'état SST accepté et le flux massique autoritaire, et ne doit pas annoncer de convergence globale lorsque seule une sous-étape a convergé.

Pour deux surfaces grises diffuses, le S2S utilise les radiosités et facteurs de forme qui doivent respecter la fermeture des lignes et la réciprocité pondérée par l'aire. En CHT, les puissances transférées de part et d'autre d'une interface imperméable doivent être égales et opposées ; une résistance de contact impose un saut de température explicite.

Référence de formulation SST : F. R. Menter, « Two-equation eddy-viscosity turbulence models for engineering applications », *AIAA Journal*, 32(8), 1598–1605 (1994), DOI 10.2514/3.12149. La variante SST-2003 et ses limiteurs doivent être comparés à la référence NASA TMR avant qualification : https://tmbwg.github.io/turbmodels/sst.html.

## 3. Discrétisation Volume Fini

Les flux convectifs et diffusifs sont assemblés par face sur volumes de contrôle. Le bilan de conservation est calculé séparément de la matrice et du second membre, à partir des flux effectivement appliqués et de l'état retourné. La correction non orthogonale, les gradients de face et les flux aux frontières doivent être partagés ou comparés explicitement entre l'assemblage et l'audit indépendant.

Pour la chaleur, le flux diffusif interne avec conductivité discontinue doit employer une résistance série, et non une moyenne arithmétique non justifiée. Pour les résidus SST, il faut distinguer le résidu du système linéaire candidat, le résidu non linéaire avant relaxation et celui de l'état effectivement retourné après relaxation/bornage.

Les preuves de fermeture attendues comprennent des oracles analytiques, des MMS sur plusieurs maillages avec ordre observé, des bilans locaux et globaux indépendants, ainsi que des tests négatifs qui montrent que les portes rejettent un état physiquement invalide même si le solveur linéaire annonce une convergence.

## 4. Traitement des frontières

Les conditions de paroi SST doivent valider avant mutation les patches non vides, les identifiants de faces frontières valides, l'absence de doublons, les propriétaires géométriques et la finitude des valeurs imposées ; la condition de paroi sur \(\omega\) exige des valeurs strictement positives. Les gradients de face doivent respecter les conditions scalaires utilisées par le transport.

En thermique, les flux imposés ou radiatifs doivent conserver leurs unités et leur signe en W/m², indépendamment de la conductivité de la cellule adjacente. Les interfaces CHT doivent vérifier la correspondance géométrique, l'opposition des normales, l'absence de fuite massique pour une interface imperméable et la cohérence entre la conductivité de l'opérateur énergie et celle utilisée pour l'interface.

Les validations d'entrée doivent précéder tout appel de solveur ou mutation d'état ; une entrée invalide doit être rejetée sans écrêtage silencieux ni changement partiel des champs.

## 5. Pièges connus

- Une valeur \(\Delta k\) ou \(\Delta\omega\) faible ne prouve pas que le résidu de l'équation évalué sur l'état borné soit faible.
- Utiliser une production brute dans l'équation de \(\omega\) alors que celle de \(k\) est limitée mélange deux contrats de fermeture.
- Un bilan global fermé peut cacher deux erreurs locales qui se compensent ; publier les deux métriques.
- Une conductivité de contact différente de celle de l'opérateur énergie crée un flux incohérent même si chaque sous-solveur converge.
- Une température de cellule n'est pas automatiquement la température radiative de face ; reconstruire la température physique cohérente avec la condition thermique appliquée.
- Une exécution locale avec un binaire lié à MPI ne constitue pas une preuve de parallélisme multi-rang.
- Un checkpoint sans flux massique de face autoritaire et sans historique requis ne prouve pas l'équivalence de trajectoire après reprise.
- Une campagne ciblée ne qualifie pas les maillages skew/non orthogonaux, les géométries radiatives générales, les propriétés variables ou les autres modèles RANS non testés.

## 6. Vérification, validation et reproductibilité

Pour chaque changement, consigner le commit exact, les versions d'outils, la configuration de build, les commandes complètes, les codes de retour, les nombres de tests, les journaux bruts et leurs empreintes SHA-256. Les tests de contrat et les oracles analytiques établissent la vérification de code ; les études de raffinement et d'erreur itérative établissent la vérification de solution ; la comparaison aux données DNS/expérimentales avec conditions comparables est nécessaire à la validation physique.

Les références de qualification SST doivent inclure des cas adaptés au modèle et au régime étudié, par exemple les données de canal turbulent de Moser, Kim et Mansour et les cas de la NASA Turbulence Modeling Resource. Toute différence de conditions, de convention de pression, de définition de coefficient ou de traitement de paroi doit être consignée plutôt que masquée par une tolérance élargie.

## 7. Statut d'implémentation

Le noyau SST, certains chemins transactionnels de production, les stages énergie/S2S et les contrats CHT ont des implémentations et des tests ciblés décrits dans le dossier détaillé ci-dessous. L'intégration complète du schéma de régions, du pilote multi-régions, du contact par interface, du restart global et de la qualification physique n'est pas considérée comme terminée tant que les critères de fermeture du registre TR-00 à TR-16 ne sont pas satisfaits sur le code fusionné dans `master`.

## Registre d'avancement (à maintenir en premier)

### Contrôles d'intégration à chaque reprise

- [ ] Lire ce registre et la section « Travail restant » avant de modifier le code.
- [ ] Relire le HEAD de `master`, l'arbre du dépôt, et les PR ouvertes pertinentes ; comparer chaque fichier concerné.
- [ ] Pour chaque correction candidate, vérifier si elle est déjà sur `master`, dans une PR ouverte, ou seulement dans l'annexe patch de ce document.
- [ ] Ne pas compter deux PR concurrentes comme deux corrections indépendantes. Une correction n'est « intégrée » qu'après vérification du contenu de `master`.
- [ ] Pour une branche de travail, consigner SHA de base, SHA de tête, fichiers modifiés, commandes, versions, code de retour, résultats et limites.
- [ ] Vérifier les contrôles CI après le dernier push ; des runs verts sur un HEAD antérieur ne certifient pas le nouveau HEAD.
- [ ] Ne fusionner qu'après résolution des conflits, audit du diff complet, CI verte et absence de défaut bloquant. Préférer plusieurs incréments ciblés plutôt qu'un patch géant non attribuable.

### État initial confirmé le 2026-10-24

| Sujet | État constaté | Décision / prochaine action |
|---|---|---|
| Patch v2 de l'annexe | Le dossier contient un patch historique annoncé comme 58 fichiers ; cela ne prouve pas qu'il soit encore intégralement absent de `master`. | Comparer fichier par fichier avec `master` et PR actives avant toute réapplication. |
| Intégration thermique/rayonnement/turbulence | PR [#1101](https://github.com/stevendaix/cfdx_dev/pull/1101) ouverte ; elle a un HEAD historique avec des runs verts, mais sa branche était en retard/divergence par rapport à `master` au dernier contrôle. | Ne pas fusionner sur la base de runs historiques ; rebaser/réconcilier, auditer le diff et relancer la CI sur le HEAD final. |
| Limiteur de production SST ω | PR concurrentes [#1094](https://github.com/stevendaix/cfdx_dev/pull/1094) et [#1099](https://github.com/stevendaix/cfdx_dev/pull/1099). | Choisir un seul correctif après comparaison exacte ; vérifier la source de production et les tests d'oracle avant fusion. |
| Oracles non finis | PR [#1093](https://github.com/stevendaix/cfdx_dev/pull/1093), [#1098](https://github.com/stevendaix/cfdx_dev/pull/1098) et [#1100](https://github.com/stevendaix/cfdx_dev/pull/1100) se recouvrent partiellement. | Consolider sans perdre les cas NaN/Inf ; une comparaison non finie doit échouer explicitement. |
| Conductivité CHT | PR [#1095](https://github.com/stevendaix/cfdx_dev/pull/1095) et [#1097](https://github.com/stevendaix/cfdx_dev/pull/1097) se recouvrent. | Vérifier le contrat commun de conductivité et l'ordre des validations avant de choisir un seul patch. |
| SST production E2E/restart | PR [#1045](https://github.com/stevendaix/cfdx_dev/pull/1045) et [#1068](https://github.com/stevendaix/cfdx_dev/pull/1068) restent à auditer selon leur HEAD réel. | Vérifier l'acceptation transactionnelle, les flux autoritaires et la reprise ; un titre WIP n'est pas une preuve d'intégration. |
| Transport thermique | La v2 documente l'API de conductivité turbulente et des tests ciblés ; elle note encore l'absence de contrat de cas/régions complet et de pilote général CHT. | Garder le sélecteur RANS fermé jusqu'à l'intégration de production et à la qualification. |
| S2S | Le noyau et des chemins de production existent, mais les limites de géométrie, de relaxation, de reprise et de MPI sont explicitées dans l'audit. | Conserver les limites et les cas négatifs ; ajouter les campagnes manquantes sans prétendre à une qualification générale. |

### Journal d'exécution — 2026-10-24

| SHA | Fichiers | Changement | Vérification | État |
|---|---|---|---|---|
| `04f0af55162f07e20278108a63318e6401eed77c` | `src/cfdx/physics/sst_solver.h` | La source ω SST-2003 réutilise la production P̃k limitée, avec branche explicite νt=0. | Revue du diff effectuée ; CI de cette branche à exécuter. | EN COURS |
| `6ab98658cd53c0bbaba8ed8edb6f3882fae6f377` | `tests/unit/test_turbulence_qualification_equations.cpp` | Oracle de limiteur actif et cas νt=0 ajoutés. | Valeurs de référence indépendantes codées ; exécution CTest encore requise. | EN COURS |

| `9148030077a60d4d96c66946f00f1af1e2462607` + `e4089770512f147e52b5470a5f3e3ef783800e95` | `src/cfdx/physics/cht_solver.h`, `tests/validation/test_cht_validation.cpp`, `docs/audits/10-thermique-cht.md` | Rejet avant mutation des conductivités d'interface incohérentes avec les deux solveurs énergie ; test négatif vérifiant la préservation des températures. | Revue du diff effectuée ; CI de la branche empilée à exécuter. | EN COURS |

| `8a4accbf2eb846449463a16db174c6d6f0b12330` | `tests/unit/test_turbulence_qualification_equations.cpp` | Corrigée une erreur de portée détectée par le premier Fast Guard : le test du cas νt=0 utilise maintenant sa propre valeur de strain. | Fast Guard initial en échec sur ce défaut de compilation ; correction poussée, CI relancée sur le HEAD courant. | CORRIGÉ — revalidation en cours |

| `6ef95ee36b3870498f366429f71c483dbf500aa4` | `tests/validation/test_cht_validation.cpp` | Étendu le test négatif aux conductivités NaN, infinie, nulles et contradictoires ; chaque rejet vérifie que T1/T2 restent inchangées. | Revue statique faite ; exécution CTest sur le HEAD final encore attendue. | EN COURS |

| `467255ef592dfdc82b34036222753cd2e0ad97d9` | `scripts/check_audit_gate.py` | Le gate classait à tort `audit_v2.md` comme fiche méthode à sections strictes 1–7 ; le dossier est maintenant explicitement classé comme rapport non-fiche, sans relâcher les contrôles des fiches méthode. | Le gate avait échoué sur cette classification ; revalidation CI requise sur le HEAD courant. | EN COURS |

**Décision TR-01 :** correction proposée dans une branche fondée sur le `master` courant. Les PR #1094 et #1099 sont des variantes historiques concurrentes ; leurs HEAD ne sont pas utilisés directement car ils ont divergé du `master`. Cette implémentation doit passer CI avant d'être considérée intégrée.

**État des preuves :** aucune preuve physique SST n'est revendiquée par ces deux tests algébriques. N9/qualification physique ne doivent pas être considérés comme réussis du seul fait que le test d'équation passe.

### Matrice de travaux — état à mettre à jour après chaque PR

| ID | Chantier | État de suivi | Dépendances | Preuve exigée pour fermer |
|---|---|---|---|---|
| TR-00 | Déduplication patch v2 vs `master` et PR ouvertes | EN COURS | Aucun | Inventaire des 58 chemins, diff par chemin, décision `already-present / PR-open / missing / conflicting`, zéro écrasement silencieux |
| TR-01 | SST : production limitée de ω et cohérence production de (k) | EN COURS — correction proposée sur branche `fix/sst2003-omega-production-current-master`; CI en attente | TR-00 | Référence SST-2003, test de non-finis, contre-exemple pour ν_t=0, résidu évalué sur l'état accepté après bornage |
| TR-02 | Conditions de paroi SST et données d'entrée | À AUDITER | TR-00 | Patch non vide, identifiants de faces valides et uniques, omega strictement positif, validation avant mutation de U/p/k/ω/ν_t |
| TR-03 | Résidus SST de l'état accepté | À QUALIFIER | TR-01, TR-02 | Résidu non linéaire recalculé indépendamment après projection/bornage ; un état bloqué au plancher ne doit pas être déclaré convergé par le seul petit Δchamp |
| TR-04 | Gradient de face et contrainte transposée turbulente | À QUALIFIER | TR-02 | MMS sur au moins trois maillages, ordre observé documenté, bilan indépendant et cas μ uniforme |
| TR-05 | Cas déclaré et schéma régions/interfaces | OUVERT | Contrat HDF5/C++/Python | Aller-retour identités/unités/interfaces et rejets déterministes avant allocation/solveur |
| TR-06 | Énergie de production avec λ_eff = λ + ρ c_p ν_t/Pr_t | PARTIEL — API/tests décrits, pilote général ouvert | TR-01, TR-05 | État SST accepté, flux massique autoritaire, variations de ν_t et Pr_t, limite laminaire, bilans globaux/locaux et rejets NaN/≤0 avant mutation |
| TR-07 | Pilote CHT multi-régions | OUVERT | TR-05, TR-06 | Oracle analytique deux régions, bilan par région et somme globale, perturbation locale détectée, rollback de toutes les régions sur échec |
| TR-08 | Conductivité variable aux interfaces CHT | OUVERT | TR-07 | Conductances utilisant les champs des cellules propriétaires, maillage asymétrique, λ discontinu et oracle analytique |
| TR-09 | Résistance de contact par interface | OUVERT | TR-05, TR-08 | Trois zones avec deux R'' distinctes, températures de surface distinctes et limite R''→0 |
| TR-10 | Couplage S2S et température réelle des parois | PARTIEL — chemins ciblés décrits | TR-06, TR-07 | Constitutif radiatif recalculé depuis les flux/températures appliqués, point fixe convergé, cas non physiques rejetés |
| TR-11 | Restart transactionnel multiphysique | OUVERT | TR-05, TR-07, TR-10 | Calcul continu vs interruption/reprise : U/p/T/k/ω/flux de face/historique identiques sous seuil fixé ; rollback atomique en cas de checkpoint invalide |
| TR-12 | Propriétés thermophysiques variables | OUVERT | TR-05, TR-06 | λ(T) avec MMS ; cp(T) via enthalpie h(T)=∫cp(T)dT ; unités, bornes et extrapolation explicites |
| TR-13 | Import OpenFOAM/Fluent et exécutable de production | PARTIEL | TR-05 | Cas E2E chargé par le parseur réel, propriétés/gravitée/BC préservées et physique exécutée déclarée explicitement |
| TR-14 | Qualification RANS de référence | OUVERT | TR-01 à TR-04 | Canal DNS Moser–Kim–Mansour, plaque plane, NACA0012 ; y+, frottement/Cp, GCI, bilan et conditions comparables documentés |
| TR-15 | Campagnes MPI et tailles représentatives | OUVERT | TR-07, TR-10, TR-11 | Exécutions multi-rang réelles, réduction de sommes/échanges vérifiée, parité sérielle/distribuée et budget mesuré |
| TR-16 | Revue des autres modèles RANS | OUVERT — SST seulement ciblé | TR-01 à TR-04 | Audit spécifique de bornage, résidus, BC, fermeture et benchmarks avant activation de chaque modèle |

### Règles de mise à jour de ce suivi

1. À l'ouverture d'une session, relire cette matrice puis contrôler le code courant ; ne pas supposer qu'un statut historique est toujours actuel.
2. Après chaque correction, ajouter le SHA du commit, le chemin, les tests exacts et leur résultat à la section de preuves correspondante. Si une PR concurrente couvre le même défaut, l'indiquer et ne pas dupliquer le patch.
3. Déplacer une ligne vers « fermé » uniquement lorsque le code est sur `master` et que la preuve requise existe sur le commit fusionné. « PR verte », « test unitaire vert » ou « API présente » ne suffisent pas à qualifier une capacité physique.
4. Toute campagne non exécutée, non distribuée ou sans référence mesurable doit rester explicitement ouverte.
5. Les annexes historiques et le patch reproduit ci-dessous ne doivent pas être modifiés pour faire disparaître un échec : ajouter une nouvelle entrée datée et préserver les résultats antérieurs.

---

CFDX — Intégration thermique, rayonnement et turbulence
Révision 2 — 24 octobre 2026. Code amélioré, essais exécutés et dossier de revue.
Dépôt : https://github.com/stevendaix/cfdx_dev. État GitHub actualisé à 16:46 UTC ; dossier assemblé à 2026-10-24T17:10:00.944265+00:00.
1. Résultat et portée exacte
Le lot local a été amélioré et vérifié : 43/43 cibles CTest en Release, 9/9 cibles sous AddressSanitizer et UndefinedBehaviorSanitizer, et 87 tests Python existants réussis, avec 4 tests ignorés faute de dépendance PyFluent. L’intégration générale fluide–solide–SST–rayonnement reste à terminer.
La révision 1 recensait les PR et assemblait 27 fichiers de corrections. Cette révision conserve ce socle et corrige des défauts supplémentaires dans les opérateurs, les frontières, les bilans, la reprise et les exécutables. Elle fournit le patch complet de 58 fichiers contre le master, et pas seulement une liste de recommandations. Aucun merge ou push distant n’est effectué par ce dossier.
Le point essentiel est la séparation entre trois états : un noyau mathématique présent, un chemin réellement raccordé à une application, et une méthode qualifiée sur des références physiques. La v2 améliore les deux premiers dans un périmètre explicite ; elle ne transforme pas les tests locaux en qualification générale RANS, CHT ou multiphysique.
Identité du code livré
Élément
Identifiant exact
Base d’application du patch
038e3aeff40f02af95b8845b931245ce23a1a094
Arbre Git du lot v2 testé
89adbb733827bdc932ff546ecf701eec11fe65f6
SHA-256 du patch UTF-8 en annexe A
1cef46e6ad105f0b1818bceeffeff313f5af1bacb8af4546b6a6fa9dbde1ebf3
Fichiers du patch complet
58
Arbre du précédent lot v1
407288b3726ab94acdd90f62cabdd863105bb499
SHA-256 du précédent patch v1
779064319468d3a40e94513621efcbcdd997892b2e1588d777655f7bec156b51
Un arbre Git identifie le contenu complet des sources ; il ne s’agit pas d’un commit publié. La reproduction indépendante en annexe reconstruit cet arbre à partir du master et du patch extrait du Markdown. Les chemins de fichiers cités dans les chapitres sont ceux du dépôt ; leurs modifications intégrales figurent en annexe A. Les numéros de lignes des notes sont des repères de lecture, les symboles et l’arbre étant les identifiants stables de cette livraison.
État après les corrections
Capacité
Ce qui est établi dans cette v2
Ce qui reste ouvert
Énergie
Source implicite dans les bilans, conductivité par cellule, frontières cohérentes, admissibilité absolue de T
Propriétés variables en T et formulation en enthalpie
Énergie issue de SST
lambda_eff = lambda + rho*cp*nu_t/Pr_t et test sur un vrai état accepté
Configuration Pr_t, construction de l’état SST et composition des boucles dans l’application
CHT
Contrats géométriques et physiques de l’interface renforcés, deux régions contrôlées avant mutation
Schéma et pilote multi-régions, contact par interface, reprise globale
Rayonnement S2S
Température de face reconstruite selon la condition thermique, flux appliqué traçable, conductivité effective transmise
Géométrie générale, parois non isothermes, qualification physique élargie
P1
Sources auxiliaires corrigées et contrats numériques renforcés
Raccordement volumique production, parois et cas couplés de référence
SST
Parois complètes, gradients cohérents, résidus après bornage, tenseur visqueux variable, retry et reprise ciblée
Activation depuis un cas, MMS complètes k/omega, références DNS/expérimentales
Reprise thermique
T enregistrée/rechargée et transmise au premier forçage Boussinesq
Persistance du flux de faces et des historiques nécessaires à une équivalence de trajectoire
Exécutables
Physiques non traitées et temps physique refusés ; verdict global honnête
Pilotes transitoire et multiphysique complets
Vérifié ici signifie qu’un oracle analytique, discret indépendant ou un contrat discriminant a été satisfait sur les cas décrits. Validation physique exige encore les références et études de sensibilité indiquées au chapitre 7. Aucun nom de fichier contenant « validation », résidu convergé ou test vert ne suffit à promouvoir une méthode.
Organisation du dossier
1.	Résultat et identité du code.
2.	PR, changements de versions et composition.
3.	Améliorations mesurées et limites des comparaisons.
4.	Détail thermique et CHT.
5.	Détail rayonnement.
6.	Détail SST.
7.	Travail restant, dépendances et critères de clôture.
8.	Tests, résultats négatifs conservés et reproduction.
9.	Inventaire exact des fichiers.
10.	Annexes : patch complet, preuves consolidées, historique et sondes.
2. PR examinées et choix d’intégration
Les socles du dépôt étaient déjà présents : conduction T1/énergie T2, boucle S2S, boucle de flottabilité et runtime SST transactionnel issu de #1047. Le problème d’intégration tient aux interfaces entre ces composants, aux contrats incomplets et à la preuve numérique de leur combinaison.
Le socle de corrections v1 conservé reprend #1079, #1080, #1081, #1086, #1088 corrigée, #1089, #1094, #1097 et #1100. La base contient déjà #1084. S’y ajoutaient les sous-ensembles F04/H04/H05 de #1096 et les validations précoces des deux régions de #1095. Les résolutions ont conservé simultanément les tests et les contrats utiles, sans remettre une ancienne version des solveurs sur les versions actuelles.
Actualisation des PR depuis la première révision
PR
HEAD observé le 24 octobre à 16:46 UTC
Décision actuelle
#1096
265ad02b432fcf45c63ded89926eb75d422bfb10
La fin de CMake a été restaurée : la configuration complète passe. L’ancien défaut de troncature n’est plus un constat courant. La PR complète n’a pas été compilée et qualifiée ici ; le lot continue à utiliser son sous-ensemble thermique sélectionné.
#1092
15bb69c7917518e2bacf1017a4002bbddc65d57b
Version mise à jour, avec différences dans CG, communication avoiding et tests. Non intégrée dans le présent arbre ; interactions gradients/Krylov à vérifier sur une combinaison séparée.
#1078
8c00dfd20e90ec648ff200a3acd40d918ff99a74
HEAD récupéré ; correctif du chemin couplé hors du lot testé. Aucun résultat de cette v2 ne le qualifie.
#1088
fefc7bb9b87cb8ae01ca217d3af71621f5e212b7
Version corrigée déjà représentée dans le lot v1 ; profils Diffusion des deux tests CG conservés.
#1084
d2f1c36426669623ce3c9d23706bbdfc3e593078
Fusionnée dans la base ; ne pas appliquer une seconde fois.
#1099
7ad610f9e3d189c6db8912c0098d12318982523d
HEAD inchangé par rapport à la reproduction d’échec de compilation ; exclue. Le correctif SST retenu vient de #1094.
La configuration de #1096 a été exécutée dans un checkout isolé avec CMake 4.4.4, Release, MPI et HDF5. Son succès clôt le blocage CMake de l’ancienne version 314f3603…, sans fournir une preuve sur tous les autres changements de la PR. Les métadonnées API et les références Git récupérées priment sur une page de PR mise en cache.
Ajouts de la révision 2
Les corrections de cette révision sont des changements locaux au-dessus du lot v1 : thermique et conductivité effective, SST et contraintes, température de paroi radiative, import des propriétés et de la gravité, reprise T, gardes des exécutables et tests discriminants. Le script T2 reprend sélectivement la logique de verdict de campagne de #1055 ; il conserve le solveur énergie actuel.
Les PR #1045, #1093 et #1098 ne sont pas additionnées mécaniquement à leurs successeurs ou concurrents. #1068 reste une liste de travaux : ses intitulés ne prouvent pas les capacités E2E ou de reprise. L’annexe historique conserve la table détaillée de la première inspection et ses anciennes versions ; les décisions courantes sont celles de cette section.
3. Résultats améliorés, avec oracles et limites
Défaut ou capacité
Avant correction correspondante
Après correction
Portée de la preuve
Source implicite énergie, Su=600, Sp=-2, T=300
Diagnostic indépendant à −600 W malgré l’équilibre
0 W ; à T=250 K, le défaut attendu reste −100 W
Source physique évaluée sur le champ accepté, pas seulement sur l’assemblage
Conduction à deux matériaux
Enveloppe énergie sans champ de conductivité transmis de bout en bout
T={360,310} K, flux 160 W pour lambda={1,4}
Solution analytique en résistances série, erreurs exigées ≤1e−10
Chaleur turbulente
Adaptateur nu_t → lambda_eff absent
Écart de 0,0179786 K par rapport au calcul moléculaire ; bilan 4,51195e−13 W
Véritable état SST accepté sur 3 cellules ; oracle tridiagonal indépendant ; pas un canal turbulent globalement qualifié
S2S évalué avec les températures de frontière
Erreur de loi radiative 13,9313257111 W/m²
1,4211e−14 W/m²
Même cas discret, seule l’ancienne reconstruction de paroi est isolée pour la comparaison
Flux S2S du cas comparatif
18,0171885342 W/m²
20,5158709050 W/m²
Correction du point fixe physique discret ; une hausse de flux n’est pas un critère universel d’amélioration
Gradients k/omega aux frontières
Produit grad(k)·grad(omega)=0,5 sur les cellules de bord
2, conformément au champ affine exact
Frontières réelles et omega de face utilisés
SST après plancher
converged=true, résidus annoncés nuls
Refus correct : résidus k=0,09 et omega=0,075
Les diagnostics décrivent le champ effectivement retourné
Divergence du stress transposé variable
Terme absent, erreur L2 voisine de 0,138 sans convergence
L2 : 0,00437466 → 0,00111903 → 0,000282884 ; ordres 1,96692 et 1,98397
MMS du terme visqueux sur 16/32/64 cellules ; ce n’est pas une MMS des deux équations SST
Reprise chaude SST
Arrêt linéaire relatif trop tôt, défaut omega 1,000015274e−8 pour exigence 1e−8
Tolérance interne resserrée ; défaut 2,766735478e−9
La continuation sans DAT reproduisait le défaut ; aucune limite externe n’a été relâchée
Température finale non physique
Un résultat à −100 K pouvait fermer linéaire/global/local et être accepté
Verdict refusé avec compteur et porte d’admissibilité explicites
Contrat de température absolue du stage production ; pas d’écrêtage dans le noyau mathématique
Reprise thermique applicative
Le checkpoint courant omettait T
T=350 K conservée et relue ; reprise sans T ou avec T≤0 rejetée
Exécutable et HDF5/DAT réels ; reprise stationnaire uniquement
Sélection solide/transitoire
Déclarations de physique non traitées parfois ignorées
Rejet explicite avec code de sortie et résumé d’échec
Inclut sélection canonique BDF2, inconnues temporelles et propriétés absentes
La correction S2S demande 23 itérations externes au lieu de 11 sur le cas comparé, avec la même relaxation 0,25, le même maximum de 200 itérations et la même tolérance 1e−10. La quantité améliorée est la cohérence entre loi radiative et température de frontière. L’ancien bilan global nul masquait cette incohérence ; le rapport ne présente donc pas la correction comme une accélération.
Contrats numériques préservés
Les seuils externes, références analytiques et budgets d’itérations n’ont pas été relâchés pour produire un résultat vert. Le critère SST vérifie en plus le défaut |A*phi-b|/sum(|a_ij|) après bornage. La précision du solveur linéaire interne est resserrée lorsque son résidu relatif ne peut garantir ce défaut absolu ; la formule et le cas de reprise sont détaillés au chapitre 6.
Une assertion ancienne acceptant un champ de déformation SST de dimension 2 a été remplacée par l’attente d’un rejet précoce : c’était une entrée incompatible avec le contrat scalaire. Les nouveaux essais positifs couvrent la physique raccordée. Les changements de fixtures d’application fournissent les métadonnées HDF5 et la sélection numérique explicitement exigées par le code ; ils ne suppriment aucune condition d’acceptation.
Les quatre campagnes Python nouvelles couvrent aussi le verdict T2 : un succès d’exécution ne vaut plus succès Poiseuille ou Nusselt. --require-qualification permet d’exiger le succès des portes de cette campagne. La campagne T2 physique complète n’a pas été relancée ici : six tests vérifient la logique de décision, et non une nouvelle valeur de Nusselt.
La suppression des propriétés d’air implicites est un changement volontaire de compatibilité : un ancien cas incomplet peut désormais être rejeté et doit déclarer les propriétés requises par sa physique. Les tests de schéma, d’import et de conversion couverts passent, mais toute l’interface graphique et tous les producteurs de cas n’ont pas été rejoués. Cette limite doit être prise en compte lors de la revue d’intégration.
4. Thermique, CHT et chaleur turbulente — analyse détaillée
Clôture commune v2 : arbre 89adbb733827bdc932ff546ecf701eec11fe65f6, 43/43 CTest Release et 9/9 CTest ASan/UBSan. Les essais ciblés ci-dessous conservent leur chronologie ; la campagne finale et son périmètre font foi au chapitre 8. Les noms des preuves locales servent à les identifier ; les résultats principaux sont repris en annexes B et D.
Date : 24 octobre 2026. Inspection et modifications locales dans cfdx_integration_v2, basé sur le master 038e3aeff40f02af95b8845b931245ce23a1a094. Le lot précédemment vérifié de 27 fichiers (arbre 407288b3726ab94acdd90f62cabdd863105bb499) était déjà appliqué avant cette nouvelle intervention. Cette note distingue ce master, ce lot antérieur et les améliorations v2 ; les améliorations locales ne sont pas présentées comme fusionnées.
Ce qui a été réellement amélioré dans la v2
Trois raccordements du noyau sont maintenant présents : la source implicite dans les bilans indépendants, une conductivité effective par cellule dans toute la chaîne énergie, et des contrôles explicites du contrat CHT. La température reprise peut aussi alimenter la première source de flottabilité, avec validation avant mutation. Ils réduisent des erreurs vérifiables ; ils ne constituent pas à eux seuls un pilote général de simulation fluide–solide–SST–rayonnement.
Élément
État avant v2
Modification locale v2
Point d’entrée
Source implicite
L’assemblage utilisait Su+Sp*T, mais le diagnostic ne comptait que Su
Sp*T_final est compté dans la source physique du bilan global, du bilan par cellule et de l’historique du solveur
src/cfdx/physics/energy_solver.h:reconstruct_energy_balance, energy_balance_relative, solve_energy
Conductivité variable
L’opérateur scalaire savait traiter cell_diffusion, mais l’enveloppe énergie et les audits production ne transmettaient pas ce champ
Paramètres optionnels conductivity_values et face_conditions, rétrocompatibles, transmis à l’assemblage et à toutes les reconstructions de flux
energy_solver.h:assemble_energy_equation, reconstruct_energy_face_flux, solve_energy
Chaleur turbulente
TurbulenceTransportControls::turbulent_prandtl existait, mais aucun adaptateur énergie ne calculait lambda_t
Ajout de build_effective_thermal_conductivity avec validation de nu_t, de Pr_t, des propriétés et du débordement
production_energy_stage.h:228–253
Paroi à flux imposé sur cellules de conductivités différentes
Un gradient unique q''/lambda_moleculaire était associé au patch
Le gradient effectivement utilisé devient q''/lambda_eff_owner par face ; les audits emploient la même condition
production_energy_stage.h:build_production_temperature_face_conditions:413–457
Interface CHT
Appariement par centre et aire ; interface vide, domaines superposés, pas de temps distincts et transport massique à travers la paroi n’étaient pas explicitement exclus
Rejets précoces avec diagnostics précis, avant mutation pour les contrats couverts
cht_solver.h:match_cht_interface:35–140, solve_two_region_cht:165–222
Équation et conservation
La formulation conservée est une équation en température avec masse volumique et capacité thermique constantes :
[ \rho c_p V_P\frac{T_P^{n+1}-T_P^n}{\Delta t} +\sum_f\left(\dot m_f c_p T_{upwind,f}-\lambda_{eff,f} A_f\nabla_nT_f\right) =(S_{u,P}+S_{p,P}T_P^{n+1})V_P. ]
La branche stationnaire omet le terme temporel. L’assemblage utilise encore la correction de convection bornée du noyau scalaire ; le bilan physique indépendant reste sensible au défaut du flux massique, ce qui est souhaité pour l’acceptation production.
Le nouveau reconstruct_energy_face_flux reconstruit directement des watts avec cp*mass_flux et la conductivité en W/(m·K). Il remplace la reconstruction en unités température puis multiplication par cp. Ce changement conserve la formulation physique et supprime la nécessité de diviser un nouveau champ de conductivité par cp pour les seuls diagnostics. L’assemblage et la reconstruction restent deux chemins distincts : le second ne lit ni la matrice ni le second membre.
Pour une face interne à conductivité variable, le noyau scalaire préexistant calcule la résistance série géométrique
[ D_f=\frac{A_f}{d_{Pf}/\lambda_P+d_{Nf}/\lambda_N}. ]
La v2 réutilise cette discrétisation ; elle n’ajoute pas une moyenne arithmétique de conductivité. Les distances sont les demi-distances normales déjà corrigées dans le lot antérieur.
La correction de source est évaluée sur la température acceptée, après relaxation. Un état intermédiaire incorrect doit donc toujours présenter un défaut de conservation. Le noyau solve_energy conserve ses critères de convergence existants ; le diagnostic indépendant y est enregistré, tandis que le stage production possède déjà ses portes bloquantes globales et locales. Il ne faut pas affirmer que l’ancien défaut Sp était à lui seul une porte bloquante du noyau.
API de chaleur turbulente et limites précises
L’adaptateur calcule
[ \lambda_{eff,P}=\lambda+\rho c_p\frac{\nu_{t,P}}{Pr_t}, \qquad \nu_t\ge0,\quad Pr_t>0. ]
Il reçoit une viscosité turbulente cinématique en m²/s. Utiliser une viscosité dynamique mu_t à cette entrée multiplierait à tort par une masse volumique supplémentaire. Le champ doit provenir d’un état SST accepté ; ce helper ne résout pas le transport SST et n’en contrôle pas le cycle d’acceptation.
Les signatures nouvelles sont optionnelles pour préserver les appelants existants :
build_effective_thermal_conductivity(energy_case, nu_t, turbulent_prandtl);
run_production_energy_stage(mesh, geometry, setup, mass_flux, T,
                            applied_radiative_flux, &conductivity_values);
Le solveur énergie et la reconstruction acceptent également le champ et un ScalarBoundaryFaceConditions. Le stage construit cette dernière structure lui-même, à partir des BC canoniques. Cela évite qu’un flux prescrit de 10 W/m² devienne involontairement 40 W/m² sur une cellule dont la conductivité vaut quatre fois la valeur moléculaire.
Les propriétés rho et cp restent constantes. Un champ lambda_eff fourni au solveur est gelé pendant sa résolution : la réévaluation en fonction de T, de nu_t ou de y+ appartient au futur pilote externe. La v2 n’ajoute aucune activation RANS au sélecteur de cas production.
Contrat CHT explicitement imposé
ChtInterfaceControls ajoute normal_opposition_tolerance=1e-8, appliqué à |1+n1·n2| pour les normales extérieures unitaires. Tous les paramètres numériques de l’interface sont vérifiés finis. Les patches doivent être non vides, de même cardinal, constitués de faces frontières valides et sans doublon. Les centres sont appariés à la tolérance configurée, les aires restent contrôlées par le contrat existant, et les normales doivent s’opposer.
Les deux régions doivent utiliser exactement le même pas physique, y compris le choix stationnaire dt=0. Le noyau n’implémente pas de sous-cyclage temporel : accepter deux dt différents ne correspondait à aucun schéma déclaré. Les deux champs sources, températures, flux et volumes nécessaires sont vérifiés avant la première résolution.
Le flux massique sur chaque face d’interface doit être nul. Le noyau réalise un échange conductif à interface imperméable ; sa loi de résistance ne contient aucun transfert advectif d’enthalpie entre régions. Ce contrôle s’applique aux interfaces traitées par ce noyau, et ne constitue pas une implémentation d’interfaces fluides perméables.
Les contrôles v2 ne rendent pas transactionnelle toute la résolution couplée. Si un solveur régional échoue après le début des itérations, l’ancien comportement de retour par exception demeure ; le rollback global multiphysique fait partie des travaux restants.
Preuves exécutées et portée
Compilation ciblée avec mpic++, C++17, -O2 -DNDEBUG, mêmes macros MPI/HDF5 que le core construit ; lien avec build_baseline/libcfdx_core.a. Les commandes, codes de retour, versions, nombres de cas et SHA-256 des sources/journaux sont enregistrés dans audit_evidence/integration_v2_thermal/manifest.json. Ces exécutables ont été lancés comme processus locaux uniques. Le dernier 11/11 énergie inclut la porte TH-10 ; les journaux CHT 8/8 et régression production 5/5 sont les exécutions ciblées réalisées avant cette dernière porte. La campagne consolidée du chapitre 8 confirme l’arbre final complet. La présence des dépendances MPI dans la compilation ne constitue pas une campagne de calcul distribué.
Test
Avant v2
Après modifications locales
Oracle ou exigence
Source implicite à l’équilibre et température acceptée relaxée
0/2 cas passent
2/2 cas passent
V=1, Su=600, Sp=-2. Équilibre T=300 : bilan 0 W. État relaxé T=250 : bilan -100 W, solveur non convergé
Ensemble énergie dédié
2 cas précédents seulement : nouvelle API absente
11/11 cas passent
Source perturbée, résistance série, Prandtl/nu_t, paroi multi-conductivité, état SST accepté, reprise de flottabilité, admissibilité physique et rejets avant mutation
Ensemble contrats CHT
1/8 cas passe
8/8 cas passent
Adjacence valide ; rejet de superposition, interface vide, NaN, dt différents, fuite massique, source région 2 invalide, face non frontière
Admissibilité physique TH-10
T=-100 K déclaré convergé malgré le rejet attendu (10/11 cas passent dans la v2 avant cette porte)
T=-100 et T=0 refusés ; T=1 accepté
Bilans et résolution linéaire restent verts dans les trois cas ; aucune température écrêtée
Régression production énergie préexistante
Pas relancée pour une nouvelle valeur « avant »
5/5 cas passent
Oracle Thomas, parois chauffées, flux imposé/radiatif, contrats d’entrée ; aucune tolérance modifiée
Les journaux bruts sont dans audit_evidence/integration_v2_thermal/before/energy.log, before/cht.log, after/energy.log, after/cht.log et after/production_energy.log. Le macro CFDX_THERMAL_CONTRACTS_BASELINE ne lance dans le test énergie que les deux contrats exprimables par l’ancienne API ; il ne fournit pas un faux calcul « avant » pour une capacité absente.
Les deux nouvelles cibles thermiques à conserver dans CMake correspondent à tests/validation/test_energy_integration_contracts.cpp et tests/validation/test_cht_integration_contracts.cpp. L’ajout CMake et les résultats consolidés sont consignés au chapitre 8 du présent dossier.
Résultats numériques discriminants
Dans le test de conduction à deux cellules, la longueur totale vaut 1 m, l’aire 1 m², les deux matériaux occupent chacun 0,5 m, lambda={1,4} et les frontières valent 400 et 300 K. La résistance exacte vaut 0,625 m²K/W ; le flux vaut 160 W/m² ; les centres doivent être à 360 et 310 K. Le test exige une erreur inférieure à 1e-10 K et des flux à 1e-10 W. Perturber ensuite la seconde cellule de 1 K produit un défaut local supérieur à 10 W : le diagnostic ne se contente pas de recopier le résultat d’une résolution réussie.
Dans le test de paroi à flux, chaque cellule reçoit 20 W via ses quatre faces latérales. La paroi gauche à 300 K évacue 40 W ; la droite est adiabatique. Les valeurs analytiques sont 310 et 316,25 K. Le même patch de paroi exige les gradients 10 et 2,5 K/m selon la cellule adjacente, tout en conservant 10 W/m² partout. Le stage ferme ses deux portes de conservation.
Le raccordement avec un vrai état SSTProductionState accepté est également exécuté. Le fixture tests/audit/sst_channel_fixture.h résout une étape de production sur trois cellules avec flux massique ±1,2 kg/s ; les champs k et omega sont effectivement modifiés par le transport SST. Le nu_t accepté alimente l’adaptateur de conductivité, et le flux massique autoritaire alimente l’énergie. L’oracle indépendant est un système tridiagonal de convection amont et résistances série résolu par Thomas ; l’accord exigé est 1e-9 K. L’effet maximal de lambda_t mesuré par rapport au calcul moléculaire est 0,0179786 K et le résidu énergétique indépendant est 4,51195e-13 W. Ce test ne revendique pas une convergence globale du porteur : une seule itération acceptée suffit à vérifier le passage de données, pas à qualifier un canal turbulent physique.
Le champ T repris est maintenant accepté par run_production_buoyant_outer_loop via un paramètre final optionnel. Un test avec T={310,320,330} K, Tref=300 K, beta=0,001/K et g=(0,-9,81,0) compare le premier porteur à une résolution indépendante alimentée par les accélérations manuelles {0,0981;0,1962;0,2943} m/s². Les vitesses et pressions concordent à 1e-13 ; le maximum de Uy du premier état vaut 0,104895 m/s. Les champs de taille/dimension incorrecte et les températures nulles, négatives, NaN ou infinies sont rejetés avant tout callback du porteur, sans modifier U ou p.
La régression production existante conserve des erreurs à l’oracle de l’ordre de 8.53e-13 K pour les parois chauffées et 3.02e-12 K pour le flux imposé. Le dernier cas de flux imposé présente un défaut global normalisé 2.59e-15 et un maximum local normalisé par cellule 2.41e-13, sous les tolérances inchangées. Ce sont des résultats de cette exécution locale, pas une qualification générale du transfert thermique turbulent.
Ce qui manque encore pour une intégration générale
TH-01 — Un contrat de cas multi-régions et un chargement correspondant
Code existant. src/cfdx/io/cfdx_io/case_schema.h:MaterialSpec:73–83 décrit des constantes matérielles ; BoundarySpec:104–131 décrit une frontière de maillage ; CaseSetup:182–238 porte une liste de matériaux, un unique jeu de frontières et une condition initiale. Il ne possède pas de liste typée de régions thermiques ni de table d’interfaces reliant deux patches de deux maillages. read_case_cfdx_h5 charge un Mesh. Le solveur thermique dédié exploite déjà un champ entier material_id ou zone_id via apps/thermal_solver/main.cpp:resolve_cell_zones:239.
Code absent. Aucun objet de cas n’associe de façon vérifiable une région à son maillage, son matériau, son statut fluide/solide, ses conditions initiales et ses champs. Aucun objet d’interface ne donne region1/patch1, region2/patch2, la loi de contact, les unités de résistance et l’unicité de l’appariement.
Implémentation proposée. Ajouter un schéma versionné de régions et d’interfaces, le porter symétriquement en C++, Python et HDF5, puis construire un registre de régions validé avant allocation du solveur. Garder la voie monomaillage existante explicite ; ne pas déduire la physique fluide/solide du nom d’un matériau. Les champs de zone existants peuvent représenter un cas conforme sur maillage unique, mais ils ne remplacent pas l’identité d’une interface entre deux maillages.
Critères de clôture. Aller-retour cas C++↔HDF5↔Python préservant identités, unités et appariements ; rejets déterministes de région absente, patch inconnu, matériau absent, zone non entière, interface utilisée deux fois ; exemple fluide + solide minimal chargé par le même exécutable utilisé en production. Dépendance : ce contrat précède le pilote CHT et le restart multi-régions.
TH-02 — Un pilote CHT qui relie réellement les exécutables au noyau
Code existant. cht_solver.h:solve_two_region_cht résout deux régions avec une interface conforme. production_energy_stage.h:build_production_energy_case exige exactement un matériau fluide. apps/thermal_solver/main.cpp:build_zone_specs:185 assigne toujours solid=true à la ligne 196 et résout une conduction stationnaire de zones ; ThermalZoneSpec::solid est encore marqué réservé à la CHT dans thermal_zones.h:53. Aucun appel applicatif à solve_two_region_cht n’a été trouvé dans apps/.
Code absent. Le pilote ne répartit pas encore l’écoulement sur les régions fluides, la conduction sur les solides, ni leurs flux sur un registre d’interfaces. Il n’existe pas de convergence globale liant quantité de mouvement, continuité, énergie régionale, échange CHT et éventuellement rayonnement.
Implémentation proposée. Introduire un état par région et une boucle externe qui possède explicitement l’état accepté. Geler T^n pendant toutes les itérations d’un pas physique ; propager seulement les nouveaux états non linéaires entre régions. Utiliser le flux massique autoritaire du porteur et imposer zéro à la frontière solide. Calculer et rapporter les puissances d’interface avec des signes opposés et les mêmes aires. Définir un mécanisme d’acceptation/rejet commun au lieu de laisser une région validée et l’autre partiellement modifiée après exception.
Critères de clôture. Cas CLI deux régions avec oracle de résistances, bilan individuel de chaque région et somme globale ; détection d’une perturbation d’un seul côté même si le bilan global se compense ; échec régional entraînant restauration de toutes les régions ; mêmes valeurs finales pour au moins deux relaxations à tolérance inchangée ; absence de flux massique à l’interface vérifiée dans le résultat écrit.
TH-03 — Conductivité variable et résistance de contact dans le noyau CHT
Code existant. La v2 sait résoudre l’énergie avec lambda_eff par cellule, mais solve_two_region_cht prend encore deux conductivités scalaires dans ChtInterfaceControls. Ses coefficients h1, h2 à cht_solver.h:242–244 et ses flux à 268–271 utilisent ces scalaires. Ses appels énergie à 257–258 ne transmettent pas de champ de conductivité. Le contrôle antérieur entre conductivité d’interface et conductivité régionale reste nécessaire pour ce contrat constant.
Code absent. Passer un champ lambda_eff à l’énergie ne suffit donc pas à obtenir de la CHT turbulente : l’interface doit employer les valeurs des mêmes cellules propriétaires. Elle suppose actuellement une température commune Tgamma, sans saut de contact.
Implémentation proposée. Étendre le contrat régional pour transporter les champs gelés, puis construire chaque conductance avec lambda_owner/d_normal. Pour une résistance de contact surfacique, employer q''=(TP1-TP2)/(d1/lambda1 + R'' + d2/lambda2) et deux températures de surface satisfaisant Tgamma1-Tgamma2=R'' q''. Faire passer ces deux températures aux opérateurs régionaux et à leurs audits. Ne pas appliquer la même température de surface des deux côtés lorsqu’un contact résistif est déclaré.
Critères de clôture. Oracle à deux couches de conductivités discontinues, maillage asymétrique et R''>0 ; comparaison des deux puissances avec la solution analytique ; limite R''→0 retrouvant le cas continu ; cas lambda_t=0 retrouvant le cas laminaire ; rejet d’un champ mal dimensionné avant mutation.
TH-04 — Une loi de contact par interface, plutôt qu’une seule valeur par matériau
Code existant. thermal_zones.h:ThermalZoneSpec:47–54 fournit une valeur contact_resistance par zone. build_thermal_zones:145–165 pose la résistance des faces entre zones à partir de ces deux valeurs et rejette deux valeurs positives différentes. Le calcul de flux asymétrique de #1084 est dans la base actuelle.
Code absent. Un même matériau au contact de deux voisins ne peut pas déclarer deux résistances différentes sans dupliquer artificiellement sa déclaration de zone. Le problème n’est plus le calcul de la résistance série corrigé par #1084 ; il s’agit de la représentation physique des interfaces.
Implémentation proposée et clôture. Introduire une loi identifiée par paire de régions/patches ou par ensemble de faces, avec priorité et unités explicites. Test obligatoire : trois zones, matériau central identique, deux interfaces ayant deux R'' différents, et deux sauts de température distincts conformes à l’oracle. Vérifier qu’une propriété de zone héritée ne remplace pas silencieusement une loi d’interface explicite.
TH-05 — Raccordement de l’état SST accepté et déclaration de Prandtl
Code existant. Le runtime SST possède SSTProductionState::nu_t et des étapes transactionnelles. La v2 ajoute le calcul thermique qui manquait. TurbulenceTransportControls::turbulent_prandtl existe à turbulence_transport.h:21, mais le schéma de matériau/cas ne définit pas encore une sélection thermique typée et versionnée de Pr_t raccordée au pilote énergétique. Le chemin production courant appelle le stage énergie sans champ effectif.
Code absent. Le pilote doit choisir explicitement un Pr_t, lire le nu_t de l’état accepté, construire lambda_eff, le maintenir cohérent pendant l’itération énergie/S2S et le restaurer ou le recalculer après un rollback SST. Il faut aussi décider quel état fait foi si la température modifie ensuite viscosité, flottabilité ou propriétés radiatives.
Implémentation proposée. Commencer par une chaîne monofluide SST→énergie à rho,cp,lambda,Pr_t constants avec nu_t variable. Ajouter une politique de Pr_t déclarée, contrôlée et écrite dans le rapport. Les lois de paroi thermiques et une loi Pr_t(T,y+) constituent des capacités suivantes ; la simple formule de diffusion turbulente ne les implémente pas.
Critères de clôture. Cas de bout en bout avec transport SST réellement accepté et flux massique autoritaire, changement mesurable du champ thermique quand nu_t ou Pr_t varie, limite laminaire exacte, conservation globale et locale, refus d’un Pr_t NaN/nul/négatif et absence de mutation énergétique après rejet de l’état porteur. Pour qualification physique : benchmark de canal chauffé avec maillage/paroi et Reynolds/Prandtl documentés, puis étude de raffinement et comparaison de Nusselt à une référence appropriée ; la réussite d’un raccordement logiciel ne suffit pas.
TH-06 — Propriétés dépendant de la température et équation en enthalpie
Code existant. thermophysical_models.h:ScalarPropertyControls:38–44, evaluate_scalar_property:46–68 et ThermophysicalProperties:102–114 proposent plusieurs lois de lambda, cp, rho, mu et Pr_t. Le chemin production lit pourtant les scalaires du matériau. Les recherches d’appel dans apps/ et les stages de production ne montrent pas de boucle reconstruisant ces propriétés depuis la température.
Code absent. Les lois ne sont pas déclarées dans le schéma de cas et ne sont pas itérées avec le champ résolu. La nouvelle conductivité variable est un raccordement utile pour lambda(T) mais n’autorise pas à annoncer cp(T) ou rho(T) : l’accumulation et le transport restent ceux de rho*cp*T à constantes.
Implémentation proposée. D’abord sélectionner et sérialiser les lois, leurs unités, bornes et politique d’extrapolation. Pour lambda(T) seul, geler un champ par itération et vérifier sa convergence. Pour cp(T), passer à un transport conservatif d’enthalpie h(T)=integral cp(T)dT et à une inversion contrôlée ; vérifier que l’intégrale et la dérivée respectent les mêmes bornes. La nouvelle #1096 contient des modifications de thermophysical_models.h, mais elles n’ont pas été incorporées automatiquement dans cette v2 ; les auditer et les tester séparément avant d’en dépendre.
Critères de clôture. Solution manufacturée avec lambda(T) connue, conservation d’une cellule adiabatique chauffée avec cp(T) non constant, tests table/bornes/extrapolation et inversibilité, comparaison en raffinement temporel. Aucun écrêtage silencieux de propriété ne doit masquer une sortie de domaine déclarée interdite.
TH-07 — Temps physique et restart thermique cohérents
Code existant. Le lot antérieur corrige le gel de T^n dans le couplage CHT ; la v2 impose le même dt aux deux régions. Le format restart sait déjà lire/écrire T via DatRestartFields::temperature (dat_restart.h:133–139, écriture 215, lecture et engagement 329–356). Le stage production énergétique impose encore dt=0 (build_production_energy_case, ligne 205 après ajout des commentaires).
Défaut identifié puis corrigé dans la v2. L’ancien exécutable écrivait le restart via l’API U,p sans T et la boucle de flottabilité réinitialisait T depuis le cas. Le lot global emploie maintenant DatRestartFields::temperature à la lecture avant calcul (apps/production_solver/main.cpp:886–905) et à l’écriture (1134–1138). Il refuse un checkpoint thermique dépourvu de T ou avec une température non admissible. La boucle reçoit initial_temperature_field (production_buoyancy_stage.h:367–422), valide entièrement le champ avant mutation et l’utilise pour sa première source. Le test noyau de première source et de rejet atomique est vert ; les tests CLI associés appartiennent à la campagne consolidée.
Code encore absent. Le pilote n’a pas de boucle thermique multi-régions transitoire. La persistance thermique stationnaire ne stocke pas le flux massique autoritaire, l’historique du point fixe, les états de sous-boucles et les niveaux temporels nécessaires à une trajectoire multiphysique identique après reprise.
Implémentation proposée. La reprise thermique stationnaire doit maintenant être maintenue par ses tests CLI. Pour le transitoire, déclarer le niveau temporel des champs (T^n, éventuellement T^{n-1}), l’horloge acceptée, la correspondance des régions et les données d’interface nécessaires. Une reprise de point de départ stationnaire et une reprise exacte de trajectoire temporelle sont deux contrats différents. Le pilote de reprise doit effectuer une lecture transactionnelle de l’ensemble des régions et champs couplés.
Critères de clôture. Simulation continue versus arrêt/reprise au même instant : mêmes températures, mêmes flux et même compteur de temps selon un seuil déclaré ; champ T absent/maillage différent/région absente rejetés sans mutation partielle ; pas de réinitialisation silencieuse de T lors d’un restart demandé ; perturbation du checkpoint détectée. Le lot global v2 ajoute production_time_contract.h:validate_production_steady_time_selection, appelé par les deux exécutables stationnaires : les déclarations transitoires sont explicitement refusées en attendant le pilote temporel.
TH-08 — Contrats d’import : modèle d’énergie, gravité et propriétés manquantes
Code existant et incohérence. adapters/openfoam.py:150–152 produit energy_model="energy" quand il trouve T ; adapters/fluent.py:1037 produit le même libellé. validate_production_energy_model ne l’accepte pas : il faut choisir une physique canonique. Le format de gravité Fluent était également incompatible ; le lot global v2 corrige maintenant sa sortie en trois scalaires séparés par des espaces, compatibles avec production_buoyancy_stage.h:parse_gravity_vector:132–153.
Le lot global v2 supprime les valeurs d’air implicitement inventées par le lecteur C++ HDF5 et aligne maintenant src/cfdx/python/cfdx/io/schema.py:MaterialSpec:60–72 sur des propriétés nulles/vides lorsque non déclarées. tests/python/test_imported_thermal_contracts.py vérifie les valeurs absentes, les valeurs explicites et le format de gravité. Cela ne doit pas être confondu avec un matériau « air » explicitement choisi et documenté, qui est une déclaration légitime. Le point restant porte surtout sur la classification canonique des modèles et sur les autres producteurs de valeurs par défaut dans chaque adaptateur.
Implémentation proposée. Uniformiser les modèles canoniques et la représentation de gravité à travers schéma, adaptateurs et solveurs. Ne pas convertir aveuglément toute présence de T en conjugate : un cas solide, un gaz compressible et un fluide Boussinesq ont des exigences différentes. Conserver un rapport de conversion incomplet si la source ne permet pas de décider. Distinguer propriété absente, propriété explicitement nulle et matériau de bibliothèque choisi.
Critères de clôture. Tests E2E OpenFOAM/Fluent minimal→CFDX→exécutable, avec identification explicite de la physique réellement exécutée ; cas incomplet refusé avant calcul ; un vrai cas de flottabilité importé conserve exactement le vecteur gravité ; tous les producteurs/consommateurs de CaseSetup partagent la même absence de défaut matériel implicite.
TH-09 — Limite de validité Boussinesq et composition des boucles
Code existant. production_buoyancy_stage.h:build_cell_buoyancy_acceleration:242 construit -beta*(T-Tref)*g. La boucle run_production_buoyant_outer_loop:367 résout le porteur puis énergie ou énergie–S2S, et contrôle une réponse de point fixe. thermophysical_models.h:evaluate_density:84–91 possède déjà un contrôle de |beta*DeltaT|, mais ce contrôle n’est pas employé dans la construction de l’accélération production.
Code absent. Le domaine de validité ne possède pas de seuil de cas/résultat effectivement contrôlé dans cette voie. La composition générale avec régions solides et SST n’est pas orchestrée ; les appels aux stages à production_buoyancy_stage.h:507–521 n’apportent pas encore de lambda_eff provenant d’un état SST.
Implémentation proposée et clôture. Déclarer un seuil de validité Boussinesq avec un diagnostic maximum, rejeter ou signaler selon un contrat explicite lorsque la température sort de ce domaine, et tester un cas valide ainsi qu’un cas volontairement hors domaine. À l’ajout de CHT/SST, conserver un propriétaire unique de l’état accepté, du flux autoritaire et du rollback. Tester le couplage complet avec stationnarité de U,p,T,k,omega et conservation régionale, sans déclarer convergée une seule sous-boucle.
TH-10 — Admissibilité physique de la température finale : fermé pour le stage production
Défaut reproduit. Avant cette correction, run_production_energy_stage rejetait une température initiale non positive, mais son verdict final ne requérait que des minima/maxima finis. La cellule cubique V=1 m³, six parois à 300 K, lambda=1 W/(m·K), flux massique nul et source de -4800 W/m³ a été exécutée. L’oracle 12*T=3600-4800 donne -100 K. Le résultat observé était exactement linear=1, global=1, local=1, residual=0, converged=1. Le test discriminant échouait donc sur le verdict physique, et non sur une tolérance numérique. Le journal est conservé dans integration_v2_thermal/th10_before/energy.log, avec le header d’avant correction.
Correction exécutée. ProductionEnergyStageResult (production_energy_stage.h:460–503) expose maintenant temperature_admissibility_passed et nonpositive_temperature. Les températures cellulaires finales doivent toutes être finies et strictement positives pour que cette porte soit vraie. Elle fait partie du verdict global à production_energy_stage.h:694–710, séparément de la convergence linéaire et des bilans. Les valeurs brutes restent dans le champ et les diagnostics : aucun écrêtage ne maquille l’échec. Le lot applicatif expose les deux nouveaux champs dans ses JSON et sa sortie console.
Preuve après correction. Le même oracle donne T=-100 K puis T=0 K pour les sources -4800 et -3600 W/m³ : dans les deux cas la résolution et les bilans sont verts, le résidu indépendant vaut exactement zéro, nonpositive_temperature=1, la porte physique est fausse et le verdict du stage est faux. Avec -3588 W/m³, T=1 K est accepté, avec compteur nul et porte vraie. Le test vérifie aussi que le noyau générique solve_energy résout toujours T=-100 avec succès : il reste utilisable pour les températures relatives et les solutions manufacturées.
Périmètre de clôture. Cette porte concerne le champ cellulaire en kelvins du stage production monofluide, directement ou comme sous-étape de flottabilité/S2S. L’outer S2S et la boucle Boussinesq propagent déjà son verdict négatif. Le helper S2S à deux passes retourne maintenant immédiatement le résultat structuré si la seconde énergie est refusée (production_radiation_s2s_stage.h:623–624), avant de reconstruire une température de surface à partir de cet état. Les contrôles spécifiques de surface S2S et l’admissibilité du solveur solide dédié restent des contrats distincts ; cette clôture ne revendique pas une qualification globale multiphysique.
Propagation S2S exécutée. Le cas production_s2s_second_energy_pass_preserves_temperature_rejection (tests/validation/test_production_s2s_energy_stage_vv.cpp:398–428) vérifie une première passe et une radiosité valides, puis une seconde solution énergétique négative avec convergence linéaire et bilans encore verts. Avant le retour anticipé, le calcul jetait une exception de température de paroi et perdait son résultat structuré : 5 cas passent et ce sixième échoue. Après le correctif, 6/6 cas passent, le compteur de températures non positives et le champ brut négatif sont conservés, et les verdicts seconde passe/ensemble restent faux. Preuves lues : audit_evidence/radiation_pass2_physical_gate_before.json et audit_evidence/radiation_pass2_physical_gate_after.json.
Ordre concret pour fermer les écarts
1.	Intégrer et exécuter ensemble les corrections v2 avec leurs tests, les contrôles applicatifs et les diagnostics S2S des autres volets ; conserver les journaux avant/après.
2.	Fermer les contrats d’import et de déclaration ; cela permet de construire des cas physiques sans interprétation implicite.
3.	Raccorder le runtime SST accepté au stage énergie monofluide avec Prandtl déclaré ; valider la chaleur turbulente indépendamment des régions solides.
4.	Introduire le contrat régions/interfaces puis le pilote CHT stationnaire conservatif, en utilisant exactement les mêmes champs de conductivité dans les régions et aux interfaces.
5.	Ajouter le contact par interface et les deux températures de surface nécessaires ; raccorder le rayonnement à ces surfaces physiques.
6.	Rendre l’état couplé transactionnel et reproductible en restart ; ajouter ensuite la boucle temporelle et les propriétés variables avec l’équation conservatrice appropriée.
7.	Qualifier chaque capacité avec un benchmark adapté et des critères publiés avant d’ouvrir son sélecteur applicatif. Les campagnes de Nusselt doivent distinguer exécution, conservation numérique et accord à la référence.
Actualisation de la PR #1096
La nouvelle ref refs/audit/integration_v2_pr1096, SHA 265ad02b432fcf45c63ded89926eb75d422bfb10, a été comparée au SHA 314f3603ca43d5d1eaa6151ea40e2f73f8742e3f audité précédemment : la différence est uniquement 135 lignes rétablies dans CMakeLists.txt. La fin du fichier et les blocs endif() sont restaurés. L’ancienne erreur de parsing CMake doit être conservée comme preuve historique sur l’ancien SHA et ne peut plus justifier un rejet de la nouvelle tête.
Cela ne signifie pas que les 32 fichiers et tous les sujets de la PR ont été validés ensemble : ils englobent d’autres intégrateurs temporels, solveurs linéaires et propriétés thermophysiques. Les corrections F04/H05 et le correctif scalaire H04 nécessaires à la thermique étaient déjà sélectionnés dans le lot antérieur. La nouvelle tête ne fournit toujours pas les contrôles d’interface CHT ajoutés ici ni le pilote de production multi-régions. L’intégration du reste doit être décidée par modification et preuves associées, sans reprendre une critique de compilation devenue périmée.
5. Rayonnement — analyse détaillée
Clôture commune v2 : arbre 89adbb733827bdc932ff546ecf701eec11fe65f6, 43/43 CTest Release et 9/9 CTest ASan/UBSan. Les essais ciblés ci-dessous conservent leur chronologie ; la campagne finale et son périmètre font foi au chapitre 8. Les noms des preuves locales servent à les identifier ; les résultats principaux sont repris en annexes B et D.
Date : 24 octobre 2026. Revue et modifications locales dans /workspace/scratch/cae91c5c59b2/cfdx_integration_v2, sur le master 038e3aeff40f02af95b8845b931245ce23a1a094 et le lot d'intégration déjà préparé, arbre 407288b3726ab94acdd90f62cabdd863105bb499. Aucun push ni merge réalisé par ce volet. Les identifiants ci-dessus sont les bases, pas le SHA d'un commit contenant les nouveaux changements.
Le integration_radiation_v2_manifest.json contient les SHA-256 des huit fichiers modifiés, de la bibliothèque liée et du integration_radiation_v2.patch. Ce patch est un différentiel par rapport au lot déjà stagé ; il ne remplace pas le patch complet de l'intégration.
Résultat principal et portée
Trois écarts exécutables sont corrigés : acceptation de valeurs non finies dans les helpers P1 ; géométrie S2S invalide acceptée par des entrées directes ; utilisation des températures de cellules adjacentes comme températures de paroi S2S. La dernière correction raccorde le rayonnement à la condition thermique réellement appliquée, y compris lorsque la conductivité varie par cellule.
Les quatre exécutables ciblés passent 24 cas sur 24 : 6 pour les noyaux radiatifs/P1, 8 pour S2S et sa géométrie, 6 pour le stage de production, 4 pour la boucle externe. Les nouveaux tests ont d'abord été exécutés sur les chemins défectueux afin de conserver une preuve négative. Les tolérances existantes, la relaxation externe de 0,25 et le plafond de 200 itérations sont conservés.
La preuve la plus discriminante est un cas couplé fermé : avant correction, les bilans énergétiques et radiatifs ainsi que le résidu de point fixe sont verts, mais la loi radiative évaluée avec les vraies températures de frontière présente une erreur de 13,9313257111 W/m². Après correction, cette erreur est 1,4210854715 × 10⁻¹⁴ W/m² pour le même système discret. Un bilan conservatif ne suffisait donc pas à vérifier le bon raccordement thermique. Voir les radiation_fixed_point_metrics_v2.json et la radiation_wall_fixed_point_probe.cpp.
Ces résultats sont des vérifications numériques locales de contrats précis. Ils ne raccordent pas P1 ou DOM au solveur de production et ne qualifient pas les facteurs de vue pour une géométrie industrielle quelconque.
Inventaire des changements
Les chemins ci-dessous sont relatifs au worktree cfdx_integration_v2.
Fichier
Fonctions ou contrat
Changement et raison
src/cfdx/physics/m1_m4_models.h
p1_blackbody_intensity, p1_source
Rejet de NaN/Inf et des valeurs négatives avant calcul ; rejet explicite du débordement de résultat ; conservation de la normalisation 4π, des limites transparentes et du chauffage radiatif.
src/cfdx/physics/radiation_advanced.h
validate_view_factor_patch, patch_pair_view_factor, estimate_view_factor_matrix
Validation de tous les centres, normales et aires, y compris appel par paire et ensemble singleton ; normalisation stable des normales ; rejet des séparations ou échanges non représentables ; diagonale nulle et réciprocité conservées.
Même fichier
validate_radiation_triangle_coordinates, aire/normale de triangle, intersection rayon-triangle, AABB, aire totale de la source
Validation avant calcul, commune aux entrées directes et à la construction BVH ; les rayons non finis et résultats intermédiaires non représentables ne peuvent plus être interprétés comme une intersection valide.
src/cfdx/physics/production_radiation_s2s_stage.h
production_s2s_wall_temperatures
Reconstruction face par face avec la condition thermique canonique et la distance normale ; moyenne surfacique après validation de chaque température.
Même fichier
solve_production_s2s_radiation
Paramètres optionnels pour le flux radiatif réellement appliqué et le champ de conductivité ; sortie applied_heat_flux distincte de la réponse injected_heat_flux.
Même fichier
run_production_s2s_energy_stage, run_production_s2s_energy_outer_loop
Propagation cohérente des flux appliqués et de la conductivité entre énergie, reconstruction et radiosité ; diagnostic de retard reconstruit sur le bon état de frontière.
docs/audits/11-rayonnement.md
Fiche méthode
Correction des équations S2S, unités et signes, du null-space P1/S2S, de la linéarisation partielle de l'émission ; limites explicites entre noyaux, production et qualification.
tests/validation/test_radiation_vv.cpp
6 cas
Tests P1 positifs, négatifs et de débordement.
tests/validation/test_s2s_radiation_vv.cpp
8 cas
Tests géométriques, conservation des résultats valides, rayons/triangles non finis et maintien du cas de visibilité avec écran.
tests/validation/test_production_s2s_energy_stage_vv.cpp
6 cas
Oracle de frontière, décalage tangentiel, conductivité locale, températures invalides avant moyenne, Cramer/replay et retour structuré après rejet physique de pass2.
tests/validation/test_production_s2s_energy_outer_loop_vv.cpp
4 cas
Point fixe à conductivité constante et variable, replay indépendant et cas isotherme convergeant avec une map initiale vide.
La prise en charge de cell_diffusion, le helper build_production_temperature_face_conditions et la correction du bilan S_u+S_pT sont réalisés dans le volet thermique. Ce volet les consomme ; il ne revendique pas leur implémentation séparée. Le sérialiseur de production est coordonné avec le responsable de l'application.
Formulation et conventions conservées
P1 : source et unités
Pour une température absolue T :
[ E_b=\sigma T^4\quad[\mathrm{W/m^2}],\qquad I_b=E_b/\pi\quad[\mathrm{W/(m^2,sr)}],\qquad G=4\pi J_{moy}. ]
Le helper renvoie la divergence radiative, positive pour une perte d'énergie de la matière :
[ q_{rad}=4\pi\kappa_a(I_b-J_{moy}) =\kappa_a(4\sigma T^4-G)\quad[\mathrm{W/m^3}]. ]
La source de l'énergie a le signe opposé, S_rad = −q_rad : positive, elle chauffe. Le nom J_moy évite de confondre une intensité moyenne P1 avec la radiosité surfacique J_i de S2S.
Le facteur 4π était déjà corrigé dans le premier lot. La présente révision ferme les entrées NaN/±Inf et les débordements. La validation intervient avant le raccourci κ_a = 0, si bien qu'une entrée NaN n'est pas cachée par une absorption nulle. Le calcul de l'intensité applique le préfacteur avant les puissances de T, et le calcul de la source forme d'abord la différence d'intensités : une absorption très grande et finie ne transforme plus un équilibre exact en Inf × 0. T = 0, absorption nulle, source nulle à l'équilibre et source négative en chauffage restent admis.
Pour le couplage volumique ségrégé, à G figé, S_p = −16 κ_a σ T0³ et S_u = S_rad(T0) − S_p T0. L'assemblage ajoute −S_p V à la diagonale ; le bilan indépendant doit reconstruire S_u + S_p T avec la température effectivement résolue. Il ne s'agit pas du Jacobien complet incluant la réponse de G à T. La fiche corrigée n'attribue plus cette linéarisation volumique au couplage de frontière S2S.
S2S : le système conservatif ne change pas
[ G_i=\sum_jF_{ij}J_j+G_{ext,i},\qquad J_i=\varepsilon_i\sigma T_{w,i}^4+(1-\varepsilon_i)G_i, ]
[ q_i''=J_i-G_i=\varepsilon_i(\sigma T_{w,i}^4-G_i),\qquad Q_i=A_i q_i'',\qquad q_{entrant,i}''=-q_i''. ]
Les unités de J, G et q'' sont W/m², celles de Q sont W. La formule précédente de la fiche additionnait deux expressions du même échange et était singulière pour une paroi noire ; elle a été remplacée par ces identités.
Le stage de production continue d'exiger au moins deux patches WALL participants, une émissivité dans (0,1], des facteurs fournis par le cas, une diagonale nulle, la fermeture de chaque ligne à 1e−10 et la réciprocité des échanges pondérés par les aires. Les températures imposées sur les patches participants restent rejetées par le contrat de rétroaction actuel. Les bilans de radiosité et d'énergie restent bloquants.
Les capacités ouvertes du noyau S2S ne sont pas confondues avec ce sous-ensemble fermé de production. Dans le noyau, une fraction de vue manquante peut voir l'ambiance, et une irradiation extérieure explicite remplace la contribution ambiante. Une normalisation d'une matrice ne démontre pas la visibilité physique de ses surfaces.
Reconstruction de température à la frontière
L'énergie utilise un flux positif entrant :
[ q_{entrant,f}''=q_{imposé,f}''+q_{rad,appliqué,f}'',\qquad g_f=\frac{q_{entrant,f}''}{\lambda_{eff,P}},\qquad T_{w,f}=T_P+g_f d_{n,f}. ]
La distance est celle de l'opérateur volumes finis :
[ d_{n,f}=\frac{|(\mathbf{x}_f-\mathbf{x}_P)\cdot\mathbf{S}_f|}{|\mathbf{S}_f|}. ]
Le helper réutilise boundary_normal_distance. Le test avec déplacement tangentiel du centre de face empêche de substituer la distance euclidienne. Les conditions par face provenant du volet thermique divisent le même flux imposé par la conductivité du propriétaire, et non par une moyenne de patch. Une condition de température fixée est reconstruite par sa valeur ; le contrat de sélection du stage conserve néanmoins son interdiction sur les participants radiatifs.
Chaque T_w,f doit être finie et strictement positive avant pondération. Une face négative n'est donc pas cachée par une autre face plus chaude. La température passée au système de radiosité reste :
[ T_{w,i}=\frac{\sum_{f\in i}A_fT_{w,f}}{\sum_{f\in i}A_f}. ]
Cette moyenne conserve l'approximation d'un patch isotherme. Elle ne représente pas exactement l'émission d'une paroi fortement non uniforme, puisque (moyenne T_w)⁴ diffère de moyenne(T_w⁴).
État appliqué et réponse calculée
À l'itération k, l'énergie est résolue avec F_k, puis T_w est reconstruit avec ce même F_k. La radiosité donne une nouvelle réponse R_k. Le résidu de point fixe compare R_k et F_k avant relaxation ; la mise à jour pour la prochaine énergie reste F_(k+1) = F_k + 0,25 (R_k − F_k).
Les deux maps sont conservées dans le résultat radiatif : applied_heat_flux représente l'état appliqué ; injected_heat_flux représente la réponse calculée. Tous les participants ont une entrée dans la première, y compris une valeur zéro lorsque la map d'entrée est vide. Le test isotherme converge dès la première itération et vérifie l'accès .at(name) sur chaque surface ; il protège la sérialisation et le replay.
La boucle deux-passes reste un passage retardé, sans revendication de convergence du point fixe. La boucle externe reste gouvernée par ses deux critères séparés à 1e−10, flux non relaxé et stationnarité de température.
La nouvelle porte d'admissibilité du volet thermique peut rejeter une solution linéaire à température négative malgré des bilans fermés. Après une pass2 rejetée, le stage S2S retourne désormais immédiatement le résultat énergie et le champ non écrêté ; il ne tente pas de calculer le diagnostic de retard pariétal sur ce champ. Ce diagnostic n'est pas exploitable après ce rejet. Le test dédié réduit λ à 0,01 W/(m K) dans la fixture 4×4 : pass1 et radiosité passent, puis pass2 conserve ses bilans et son statut linéaire mais échoue au critère thermique. Avant le retour immédiat, une exception de reconstruction masquait ce verdict structuré ; après correction, converged=false, temperature_admissibility_passed=false, nonpositive_temperature>0 et min_temperature<0 restent accessibles.
Preuves avant et après
Reproductions négatives
Défaut reproduit
Résultat avant correction
Résultat après correction
Preuve
P1, NaN/Inf et résultat non représentable
Deux nouveaux cas échouent ; quatre cas existants passent
6/6
radiation_v2_before.json, radiation_v2_after.json
Centroïdes/normales/aires S2S invalides ; normale 1e200
Deux nouveaux cas échouent ; cinq cas existants passent
Les huit cas S2S finaux passent
Mêmes fichiers, volet test_s2s_radiation_vv
Intersection directe avec triangle non fini
Le nouveau cas ne reçoit pas l'exception requise ; 7 autres cas passent
8/8
radiation_triangle_v2_before.json, radiation_v2_after.json
Paroi à flux imposé : T_P=350 K, q=50 W/m², d_n=0,125 m, λ=1 W/(m K)
T_w renvoyée 350 K au lieu de 356,25 K : erreur 6,25 K
Égalité à l'oracle à 1e−12 K
radiation_wall_temperature_v2_before.json, radiation_wall_temperature_v2_verified.json
Température reconstruite négative
Aucun rejet alors que la BC donne une paroi non physique
Rejet, même si la moyenne du patch serait positive
Mêmes fichiers et source du test stage
Pass2 énergie conservatrice mais physiquement rejetée
Le diagnostic pariétal jette après la pass2 et masque son verdict structuré
Retour immédiat de l'échec physique et champ brut conservé ; stage 6/6
radiation_pass2_physical_gate_before.json, radiation_pass2_physical_gate_after.json
Les échecs avant correction correspondent à de nouveaux tests discriminants ajoutés à une base dont les anciens tests passaient ; ils ne doivent pas être présentés comme les résultats initiaux de la suite historique entière.
Mesure du couplage complet
La sonde utilise le même problème discret que les tests S2S : boîte 4×4×1, conduction sans convection, parois gauche/droite 400/300 K, flux inférieur 50 W/m², λ = 1 W/(m K), ε inférieur/supérieur = 0,8/0,9 et facteurs de vue déclarés F12 = F21 = 1. Cette matrice donnée définit l'enceinte modèle ; elle n'est pas une preuve de facteurs de vue géométriquement exacts pour toutes les faces de cette boîte.
L'expérience compare l'ancien seul header S2S extrait du HEAD de base à son nouveau contenu, les autres éléments thermiques étant ceux du worktree courant. Il s'agit d'une isolation de la reconstruction de paroi, pas d'une comparaison globale entre tous les changements de master et du lot. Le flux réellement appliqué est récupéré par le bilan indépendant de l'énergie, puis la température de frontière et l'oracle de deux surfaces sont recalculés.
Grandeur
Ancienne reconstruction par cellules
Nouvelle reconstruction par BC
Itérations externes
11
23
T_w inférieure fournie à la radiosité (K)
357,0037421510
360,2597972962
T_w inférieure exigée par le flux thermique appliqué (K)
361,0015935856
360,2597972962
T_w supérieure fournie à la radiosité (K)
354,6034007061
357,5973455610
T_w supérieure exigée par le flux thermique appliqué (K)
356,8555492716
357,5973455610
Flux radiatif inférieur renvoyé (W/m²)
18,0171885342
20,5158709050
Flux attendu d'après les températures de frontière (W/m²)
31,9485142453
20,5158709050
Erreur de loi radiative (W/m²)
13,9313257111
1,4211e−14
Bilan énergie global normalisé
4,5475e−16
0
Bilan radiatif normalisé
0
0
Résidu externe de flux normalisé
6,3757e−12
8,4545e−11
Le nombre d'itérations augmente parce que la rétroaction porte désormais sur la température de frontière. Aucun critère n'a été relâché pour compenser cette modification. L'amélioration mesurée est la cohérence du problème discret et de sa loi de frontière, pas une baisse artificielle du résidu ni une preuve expérimentale.
Les valeurs exactes et les commandes de compilation sont dans radiation_fixed_point_metrics_v2.json. Le script source ne dépend d'aucun accès réseau.
Exécution et reproductibilité
Les commandes exactes sont enregistrées dans les JSON. Les tests ont été compilés avec mpic++ -std=c++17 -O2 -DNDEBUG, macros CFDX_HAS_MPI=1 et CFDX_HAS_HDF5=1, en incluant les headers modifiés du worktree et en liant la bibliothèque build_baseline/libcfdx_core.a déjà construite. Son SHA-256 figure dans le manifeste. La bibliothèque fournit les composants communs ; les changements radiatifs et thermiques concernés sont dans les headers compilés avec chaque exécutable.
Les exécutions sont locales, en un processus. Un binaire lié à MPI ne démontre pas une exécution multi-rang ni la parité MPI. Les quatre cibles CMake préexistantes restent disponibles et aucun changement CMake n'est nécessaire pour ce volet :
cmake --build build --target test_radiation_vv test_s2s_radiation_vv \
  test_production_s2s_energy_stage_vv test_production_s2s_energy_outer_loop_vv -j2
ctest --test-dir build --output-on-failure \
  -R '^(test_radiation_vv|test_s2s_radiation_vv|test_production_s2s_energy_stage_vv|test_production_s2s_energy_outer_loop_vv)$'
La reconstruction du projet complet et ses tests de non-régression appartiennent au lot global coordonné. Les résultats ciblés de cette note ne doivent pas être substitués à cette vérification d'intégration.
Ce qui reste nécessaire pour « tout intégrer »
Production P1 et DOM
Le contrat src/cfdx/physics/production_radiation_contract.h continue d'autoriser S2S et de rejeter les autres modèles radiatifs comme non raccordés à l'énergie de production. Pour ouvrir P1/DOM, il faut au minimum :
1.	Un schéma de cas pour les propriétés optiques par matériau/zone et les conditions radiatives, avec rejet des combinaisons incompatibles.
2.	La création des champs G/intensités et des quadratures, puis l'appel effectif du transport radiatif dans l'orchestration de production.
3.	L'injection énergétique selon les signes ci-dessus, avec traitement cohérent des termes implicites et reconstruction indépendante de leur bilan.
4.	Une boucle externe qui contrôle séparément la convergence du transport, la stationnarité thermique, le point fixe et les bilans globaux/locaux.
5.	Les sorties, checkpoints et reprises cohérents avec les champs nouveaux, puis une vérification MPI des sommes et des échanges.
6.	Des tests de bout en bout qui utilisent le parseur et l'exécutable de production. La présence d'un noyau et de ses tests unitaires ne clôt pas ces points.
Aucun de ces raccordements n'est annoncé comme réalisé par les garde-fous du helper P1.
Géométrie S2S et résolution des surfaces
La production utilise toujours les facteurs déclarés dans le cas. L'estimateur par centroïdes n'est pas appelé pour inventer automatiquement cette matrice. Il ne traite pas les obstacles et conserve son comportement de contribution nulle pour des centres confondus, après validation des entrées. Une aire positive et une normale finie ne prouvent pas une géométrie physique sans recouvrement.
Le noyau triangulaire dispose d'un test déterministe d'écran et d'une accélération BVH. La fermeture des entrées NaN/Inf ne qualifie pas son traitement de toutes les géométries concaves, des intersections multiples, des surfaces proches ou des tolérances absolues d'intersection. Il reste à vérifier la visibilité et la convergence en échantillons contre des oracles indépendants, ainsi que la dépendance à l'échelle géométrique et au découpage triangulaire.
Une qualification générale exige aussi une étude du découpage des patches, de l'erreur de température moyenne et de la matrice de facteurs associée. Un traitement par face ou par sous-surfaces doit conserver les identités d'aires et l'intégration du flux réinjecté. Remplacer seulement T_moy⁴ par moyenne(T⁴) ne suffit pas à définir un modèle complet de radiosités et de flux locaux non uniformes.
Stabilité, précision et domaine de validité
La relaxation fixe 0,25 passe les cas ciblés constant et variable. Elle n'est pas une preuve de robustesse pour une forte résistance conductrice, des températures beaucoup plus élevées ou un maillage arbitrairement grossier. L'épuisement des 200 itérations doit rester un échec explicite ; une méthode d'accélération éventuelle nécessiterait ses propres essais sans modifier le critère de point fixe.
Les MMS P1 à trois maillages, études de convergence angulaire DOM, GCI sur les quantités d'intérêt, références physiques du domaine visé et exécutions multi-rang restent nécessaires avant extension de qualification. La documentation corrigée sépare ces exigences de la vérification ciblée présente.
Décision d'intégration de ce volet
Les modifications ont rejoint le lot local dont la compilation et les tests ciblés consolidés sont rapportés au chapitre 8 ; la qualification globale du dépôt demeure distincte. Les interfaces optionnelles conservent les appels existants et les maps retournées sont complètes au premier état isotherme. La correction du bilan volumique S_pT et la fourniture de λ_eff doivent être intégrées avec le volet thermique ; le paramètre de conductivité doit être transmis depuis l'orchestration qui le calcule.
Les critères ouverts ci-dessus restent visibles. Le résultat final peut annoncer garde-fous radiatifs et couplage S2S de frontière vérifiés sur les cas définis, sans annoncer P1/DOM de production, occlusion générale ni qualification industrielle achevée.
6. Turbulence SST — analyse détaillée
Clôture commune v2 : arbre 89adbb733827bdc932ff546ecf701eec11fe65f6, 43/43 CTest Release et 9/9 CTest ASan/UBSan. Les essais ciblés ci-dessous conservent leur chronologie ; la campagne finale et son périmètre font foi au chapitre 8. Les noms des preuves locales servent à les identifier ; les résultats principaux sont repris en annexes B et D.
Date : 24 octobre 2026. Travail local autorisé par la demande « regarde ce qu'il manque pour tout intégrer dans cfdx, fais-le en détail et améliore le code et les résultats ». Aucun push, merge ni activation générale du sélecteur RANS.
Périmètre et états de référence
•	Worktree courant : /workspace/scratch/cae91c5c59b2/cfdx_integration_v2.
•	HEAD de départ : 038e3aeff40f02af95b8845b931245ce23a1a094.
•	Lot précédent déjà appliqué : arbre 407288b3726ab94acdd90f62cabdd863105bb499. La copie /workspace/scratch/cae91c5c59b2/cfdx_candidate expose ce même arbre dans son index, avec HEAD historique 815c000d428b74bdb7ba568853f76f725a224fff ; elle a servi au contre-exemple « avant ».
•	PR #1092 inspectée statiquement : ref refs/audit/integration_v2_pr1092, SHA 15bb69c7917518e2bacf1017a4002bbddc65d57b. Elle n'a pas été intégrée en bloc.
•	Base de règles lue : .ai/AGENTS.md. Les modifications ci-dessous sont une implémentation et des vérifications ciblées ; elles ne valent pas validation physique ou qualification de l'ensemble RANS.
Le correctif limité de production omega de #1094, le retrait de l'ancienne API ambiguë #1081 et les correctifs du lot précédent sont conservés. Le travail nouveau porte sur la couverture réelle des conditions de paroi, la cohérence des gradients turbulents, les résidus de l'état borné, le tenseur visqueux variable, le bilan de quantité de mouvement et les preuves du cycle de vie couplé.
Inventaire des fichiers et points d'intégration
Chemins relatifs au worktree ci-dessus ; lignes relevées sur la version locale du 24 octobre après corrections. Les symboles sont les références stables si un autre agent décale les lignes.
Fichier
Symboles / lignes
Intervention
src/cfdx/physics/scalar_transport_gradient.h — nouveau
validate_transport_boundary_face_values:10 ; scalar_gradient_with_transport_boundaries:43
Couverture de tous les IDs de face d'un patch ; gradient Gauss avec véritables BC scalaires et données de face prioritaires.
src/cfdx/physics/sst_solver.h
sst_boundary_gradient_dot:12 ; sst_equation_solution_defect:39 ; sst_linear_tolerance_for_solution_defect:65 ; solve_sst_transport:157 ; solve_sst_transport_dynamic_blending:252
Gradients BC utilisés pour F1 et diffusion croisée ; dimensions/itérations validées ; résidus recalculés après bornage ; gate de défaut de l'état retourné et traduction du seuil linéaire.
src/cfdx/physics/sst_production_coupling.h
restauration DAT:328 ; refresh dérivé:476 ; constructeur omega de paroi:639 ; gardes:714 et 728
Ancien contrôle sans maillage refuse les obligations wall omega ; nouvelle surcharge vérifie chaque face, la métrique géométrique et la valeur ; F2 seul ne calcule plus un gradient inutile.
src/cfdx/physics/momentum_viscous_stress.h — nouveau
transpose_viscous_stress_flux:16 ; viscous_face_flux_divergence:78
Flux conservatif complémentaire du tenseur visqueux variable.
src/cfdx/physics/steady_incompressible_solver.h
runtime:177 ; obligations wall omega:196 ; callback prédicteur:293 ; préflight:2546 ; activation stress:2561 ; prédicteur:2914 ; résidu final:3598 ; bilan:4596–4651
Préflight avant initialisation U/p ; données réellement consommées au premier prédicteur observables ; stress transposé branché aux trois chemins ; bilan corrigé en signe et densité.
tests/unit/test_sst_blending.cpp
tests dimensions et wall omega
Migration de la surcharge stricte. Une ancienne attente reposait sur l'acceptation volontaire d'un champ de strain de dimension 2 : elle exige maintenant son rejet avant mutation. Aucun seuil n'a été desserré.
tests/audit/sst_channel_fixture.h — nouveau
SSTChannelCase:10
Fixture trois cellules, vrai transport SST accepté, flux traversant non nul. Réutilisé par la vérification thermique de raccordement.
tests/audit/test_sst_integration_gaps.cpp — nouveau
9 cas à partir de :39, :79, :118, :151, :179, :228, :270, :320 et :352
Oracles décrits ci-dessous. La cible CMake est enregistrée par l'agent coordinateur.
docs/audits/09-turbulence-rans.md
section 13:111
Méthode actuelle, pression modifiée, limites, reprise et niveau des preuves. Les sections historiques sont explicitement datées.
Défaut 1 — couverture partielle WALL_OMEGA_SST
Avant
Le garde vérifiait qu'un vecteur correspondant au nom de patch comportait au moins une entrée finie. Il ne comparait pas les entrées aux face_ids du maillage. Une valeur manquante sur une autre face était ignorée par l'assemblage scalaire : le transport retombait sur sa BC de patch, et notamment sur le gradient nul implicite lorsqu'aucune BC omega de patch n'existait.
L'ancien contre-exemple est conservé dans audit_evidence/turbulence_counterexamples.cpp et son journal. Le nouvel oracle inspecte une vraie paroi de 12 faces sur un maillage de trois cellules.
Après
validate_transport_boundary_face_values exige, pour tout patch de données explicites :
1.	un patch réellement présent et non vide ;
2.	un vecteur de longueur mesh.n_faces() indexé par face globale ;
3.	chaque face_id dans l'intervalle, réellement frontière, propriétaire valide et même appartenance dans geometry.face_patch ;
4.	toutes les valeurs des faces du patch finies, et strictement positives pour omega.
Le garde spécialisé exige en plus un patch WALL, des noms uniques, une valeur pour chaque obligation déclarée, une aire positive et une distance géométrique normale positive. Il ne réemploie pas le fallback de longueur caractéristique de boundary_normal_distance pour certifier une paroi SST résolue. Une géométrie où le centre de face coïncide avec le centre de cellule est refusée même si des omega positifs ont été fournis manuellement.
Ces contrôles sont exécutés dans solve_steady_incompressible avant apply_initialization, donc avant toute modification du champ porteur due au mode Uniform ou Restart. Le test injecte NaN, +Inf, zéro et une valeur négative dans la dernière face de paroi ; U, p, k, omega, nu_t, distance et strain doivent rester bit à bit identiques. Les vecteurs tronqués, l'appartenance géométrique erronée et la métrique normale dégénérée sont aussi refusés.
Interface à respecter dans le driver
auto contract = make_sst_runtime_boundary_contract(setup);
auto geometry = build_fv_geometry(mesh);
auto values = make_sst_wall_omega_face_values(
    mesh, geometry, contract.wall_omega_sst_patches, transport_controls);
require_supported_sst_runtime_boundaries(contract, mesh, geometry, &values);

runtime.state = &accepted_sst_state;
runtime.transport = transport_controls;
runtime.k_boundary_conditions = contract.k;
runtime.omega_boundary_conditions = contract.omega;
runtime.wall_omega_sst_patches = contract.wall_omega_sst_patches;
runtime.omega_face_values = values;
// wall_patch_names reste la liste des parois servant à la distance minimale.
La surcharge sans maillage est conservée seulement pour les contrats sans obligation wall omega. Elle échoue explicitement pour les autres : aucune validation géométrique ne peut être prouvée à partir du vecteur seul. La compatibilité de l'ancien appel avec wall omega est volontairement resserrée ; le caller doit migrer vers la surcharge ci-dessus.
Défaut 2 — gradients de k et omega calculés avec une autre frontière
Formulation et discrétisation
La diffusion croisée SST utilise grad(k)·grad(omega), et cette contraction entre dans CD_komega, donc dans F1. Elle doit employer les BC du transport réel. Le nouveau gradient applique
[ \nabla\phi_P=\frac{1}{V_P}\sum_{f\in\partial P}\phi_f,\mathbf S_{P,f}. ]
Sur une face interne, la valeur est interpolée géométriquement par la projection longitudinale déjà utilisée dans le noyau Gauss. Sur une frontière, la résolution suit cet ordre : donnée explicite de face, valeur fixe de patch, extrapolation par gradient normal fixé, puis valeur propriétaire pour gradient nul. Les faces EMPTY sont exclues. La correction de skewness n'a pas été introduite.
sst_boundary_gradient_dot est appelé une fois par itération et sert aux termes sources et à F1. Les anciens calculs de gradients génériques sans BC, utilisés seulement pour recalculer F2 dans les routines de reprise/initialisation/refresh, ont été retirés : F2 dépend de k, omega, nu et distance, pas de grad(k)·grad(omega). Cela évite d'étendre inutilement les signatures DAT.
Oracle discriminante
Sur trois cellules orthogonales, k=2+x, omega=3+2x+3y, avec valeurs et gradients imposés cohérents :
•	contraction exacte : 2 dans chaque cellule ;
•	ancien gradient générique sans BC : 0,5 dans les deux cellules de bord ;
•	nouveau gradient : composantes omega (2,3,0) et contraction 2, vérifiées à 2e-14 ;
•	une valeur de patch omega artificielle 999 est surchargée par les valeurs de face : le test certifie cette priorité effective ;
•	les gradients normaux k imposés −1 à gauche et +1 à droite reproduisent le gradient x égal à 1.
Il s'agit d'une vérification de reconstruction affine sur le maillage indiqué ; ni le transport non linéaire complet ni le comportement sur maillages déformés ne sont validés par ce seul test.
Défaut 3 — faux résidu et fausse convergence après plancher
Mesure avant/après exécutée
Source reproductible : audit_evidence/integration_v2/sst_bounded_residual_probe.cpp. Même géométrie et mêmes données dans les deux compilations : une cellule de volume 1, flux et strain nuls, k=omega=k_min=omega_min=1, F1=F2=1, trois itérations au maximum, tolérance 1e-8, densité 1.
Mesure
Lot précédent, arbre 407288b3
Lot nouveau
converged
1
0
k retourné
1
1
omega retourné
1
1
résidu k rapporté
0
0,089999999999999997
résidu omega rapporté
0
0,074999999999999997
itérations
1
3
corrections de plancher k / omega
1 / 1
3 / 3
Journaux : sst_bound_before.log, sst_bound_after.log et les deux *_compile.log dans audit_evidence/integration_v2/.
L'équation homogène résolue linéairement a la solution zéro, puis le plancher replace le champ à 1. Les valeurs exactes des résidus de l'état retourné sont donc rho*beta_star=0.09 et rho*beta1=0.075. L'ancien code mesurait les vecteurs linéaires encore nuls ; le nouveau les met à jour après enforce_turbulence_bounds.
Critère conservé et gate supplémentaire
Le statut convergé des deux solveurs scalaires, le contrôle des non-finis et le critère max(|delta k|, |delta omega|) <= tolerance restent requis. Le nouveau gate ajoute
[ e_\phi=\max_i\frac{|(A\phi-b)i|}{\sum_j|a{ij}|}\leq\mathrm{tolerance}. ]
Ce défaut a les unités du champ, comme le seuil absolu de changement déjà utilisé ; il n'est pas une estimation de l'erreur de discrétisation ou une validation physique. Une ligne de norme nulle/non finie ou un résidu non fini donne un défaut infini et interdit la convergence. La normalisation reste définie avec un RHS nul, contrairement à une division pure par ||b||. Les planchers par défaut 1e-12 peuvent toujours être inférieurs à la tolérance 1e-8 ; leur résidu réellement non nul reste visible, avec leur télémétrie. Les critères et seuils externes sont conservés ; le seuil linéaire relatif est désormais traduit de façon cohérente, comme expliqué ci-dessous.
Cause de l'échec du nouveau test de reprise et correction
Le premier essai du neuvième cas DAT échouait après trois retries. Les comparaisons U/p/k/omega/nu_t restaurés et le premier prédicteur passaient ; l'échec se produisait au transport suivant. Une sonde a reproduit exactement cet échec en continuant le même état sans aucune sérialisation DAT : le format de fichier et ses valeurs n'étaient pas la cause.
Au point stationnaire approché, le solveur linéaire omega retournait CONVERGED sans changer le champ car son résidu relatif au RHS était 7.9645597175522109e-9 < 1e-8. Pourtant le défaut en unités omega valait 1.0000152740357402e-8 > 1e-8. Le nouveau gate extérieur rejetait correctement cet état ; répéter le solveur intérieur avec le même seuil relatif reproduisait le même champ inchangé (delta omega=0) et ne pouvait jamais satisfaire le gate. Augmenter le budget n'aurait pas résolu la cause.
Le correctif traduit maintenant le seuil absolu demandé en une tolérance relative interne suffisante :
[ \mathrm{tol}{lin}=\mathrm{tolerance}\min\left(1, \frac{\min_i\sum_j|a{ij}|}{\max(|b|_2,10^{-15})}\right). ]
La constante de normalisation 1e-15 correspond au plancher RHS de CG/BiCGStab existant. Si ||r||2 satisfait ce seuil, chaque résidu de ligne peut atteindre le défaut absolu demandé ; la tolérance ne peut qu'être resserrée. Les 100 itérations SST, 2000 itérations linéaires, relaxation 0,7, seuil externe 1e-8, gates de porteur et critères de tests restent identiques. RHS/normes non finis, ligne nulle, norme RHS débordée et résultat sous-débordé à zéro sont refusés explicitement. L'oracle permanent couvre notamment le sous-débordement à partir de denorm_min, le RHS infini et la ligne nulle.
Preuve négative conservée : integration_v2/sst_dat_failure.log. Sonde, sorties détaillées et commande de compilation : sst_restart_diagnosis.cpp, sst_restart_diagnosis.log, sst_restart_diagnosis_build.log. Après correction, les 9 cas dédiés passent, dont le parcours DAT précédemment en échec, avec tous les seuils externes inchangés (sst_final_agent_test.log). La même sonde recompilée donne aussi continued_without_DAT=success et un défaut omega 2.766735478334318e-9 (sst_restart_diagnosis_after.log).
Ce correctif cible SST. Les solveurs génériques k-epsilon/k-omega/SA dans turbulence_solver.h n'ont pas été promus ni globalement réécrits ; leur vérification reste requise avant toute activation des autres modèles.
Défaut 4 — tenseur visqueux variable incomplet
Équation et convention de pression
Pour l'incompressible avec viscosité tourbillonnaire :
[ \nabla\cdot(\rho U_i\mathbf U) =-\partial_i p^+\partial_j{\mu_{eff}(\partial_jU_i+\partial_iU_j)}+\rho f_i, \quad \mu_{eff}=\rho(\nu+\nu_t), \quad p^=p_{statique}+\tfrac23\rho k. ]
Les valeurs de pression et leurs BC dans le runtime SST représentent p*. Le noyau ne convertit pas automatiquement une pression statique imposée. Une sortie Cp ou pression statique doit soustraire 2*rho*k/3, avec une convention identique aux frontières ; aucun résultat aérodynamique n'est qualifié ici.
La diffusion scalaire de chaque composante était déjà dans la matrice, mais la divergence du tenseur transposé manquait lorsque mu varie. Elle n'est pas supprimable par l'incompressibilité : U=(y,0,0) et mu=1+x² donnent une contribution (0,2x,0), non nulle et de rotationnel non nul.
Discrétisation nouvelle
Le nouveau flux complémentaire est
[ F^T_{i,f}=\mu_f\sum_j(\partial_iU_j)fS{f,j}. ]
La valeur de face est unique, orientée propriétaire ; sa contribution s'ajoute à la cellule propriétaire et se retranche à la voisine. Les gradients U utilisent les BC réelles du transport. La viscosité interne suit le même coefficient harmonique pondéré par distances normales que la diffusion scalaire :
[ \mu_f=\frac{d_{Pf}+d_{Nf}}{d_{Pf}/\mu_P+d_{Nf}/\mu_N}. ]
En frontière, mu et le gradient de face restent ceux du propriétaire. Le terme scalaire est implicite et le terme transposé est explicite par itération externe. Cette option est activée pour un runtime SST ou un champ cellulaire explicite de viscosité turbulente. Le parcours historique strictement uniforme, sans l'un de ces deux modes, garde sa réduction incompressible en laplacien.
Le prédicteur reçoit div(F^T), puis le résidu non linéaire reconstruit ce terme avec U corrigé. Le bilan indépendant reconstruit le flux final convection - diffusion scalaire - F^T. Il ne réemploie pas le RHS du prédicteur. La validation de l'égalité de deux assembleurs scalaires à mu uniforme reste valide ; cela n'implique pas qu'une itération complète avec un opérateur transposé discret soit identique pour toute vitesse discrète imparfaitement solénoïdale.
MMS exécutée et limites
MMS : U=(y,0,0), mu=1+x²(1−x)², donc div(mu*grad(U)^T)=(0,2x(1−x)(1−2x),0).
N
Erreur L2 nouvelle
Erreur L2 si terme omis
Ordre observé
16
0,00437466
0,138032
—
32
0,00111903
0,138014
1,96692
64
0,000282884
0,138013
1,98397
Source et résultat exploratoires : stress_mms_probe.cpp, stress_mms_probe.log (mailles 8 à 128). L'oracle permanent exige un ordre supérieur à 1,8 sur les trois maillages 16/32/64, vérifie l'égalité du flux frontière et de la somme volumique, et vérifie une force nulle pour mu constant. L'erreur du terme omis reste d'ordre un alors que celle du nouvel opérateur décroît.
Le prédicteur de production avec viscosité variable, sans forçage transverse et avec Uy initialement nul, produit Uy_gauche=0.000596783, Uy_droite=-0.000619322 dans le test huit cellules. L'ancien chemin scalaire seul garde Uy identiquement nul. Cette oracle démontre le raccordement effectif du terme et son signe, au-delà de son existence dans un helper.
Portée : la dérivée de mu s'annule aux frontières dans cette MMS. Les valeurs propriétaire en frontière, les forts gradients de mu aux murs, les contrastes élevés, la stabilité du traitement explicite, les maillages non orthogonaux et les conditions de paroi RANS nécessitent encore des campagnes. Le résultat ne prouve pas la qualification d'un canal turbulent, même si l'ordre est proche de deux sur cette famille.
Défaut 5 — bilan indépendant : signe et densité de la force
audit_integrated_balance additionne le flux sortant et le terme source. Avec un flux convection - viscosité, la source à fournir est grad(p*) - rho*f. L'ancien code fournissait f - grad(p) sans multiplier la force massique par rho. Les deux erreurs sont corrigées dans le bilan de quantité de mouvement.
Oracle non symétrique en source : une cellule cubique, six faces fixes à vitesse nulle, rho=2,5, nu=0,4 (donc mu=1), force massique x=3, relaxation vitesse 0,7. La conductance totale vaut 12. Après le premier pas, Ux=0.7*7.5/12=0.4375 et le bilan signé exact vaut 12*Ux - 7.5 = -2.25. L'oracle passe dans le chemin viscosité cellulaire nulle et dans la limite historique uniforme. Ce test ne repose pas sur une annulation globale d'un forçage antisymétrique.
Le bilan reste un diagnostic physique stationnaire ; les contrôles de continuité et les résidus de l'équation bornée restent distincts. Une certification du bilan transitoire avec accumulation et historique appartient à une campagne séparée.
Preuves du cycle de vie : premier prédicteur, retry et DAT
Fixture réellement raccordé
Le fixture utilise trois cellules en chaîne, rho=1,2, nu=0,01, U initialisée uniformément à (1,0,0), des parois mobiles à la même vitesse, inlet k=0,01/omega=2, paroi k=0 et omega géométrique, outlet gradients nuls. Le flux autoritaire traverse le domaine à 1,2 kg/s. Les valeurs nu_t, distance et strain fournies initialement sont volontairement incohérentes (9,17,11) pour rendre le défaut de refresh détectable.
Le callback momentum_predictor_callback observe les tableaux exacts consommés juste avant assemblage. Le premier mu attendu est 1.2*(0.01+0.01/2)=0.018, la distance 0,5 et le strain nul. Le test impose ces valeurs et un flux interne non nul. L'état SST accepté après ce pas donne environ :
•	k : [0.00916065, 0.00838028, 0.00770093] ;
•	omega : [1.94838, 1.90061, 1.85891].
Ces champs non uniformes et le flux traversant rendent le test plus discriminant qu'une cellule au repos. Le fixture reste un test logiciel à parois mobiles, pas un cas physique de validation de turbulence.
Rejet après commit SST
Le callback de sortie rejette le premier essai après que le transport SST a commis son nouvel état. Le prédicteur de retry doit retrouver exactement U, phi, k, omega, nu_t, distance, strain et mu du début d'essai. Il ne peut conserver un champ issu du candidat rejeté.
Le contrôleur réduit effectivement les relaxations à alpha_u=0,35, alpha_p=0,15. L'état accepté est comparé à un calcul propre employant ces mêmes relaxations, et pas au candidat rejeté. Une seule ligne d'historique est publiée, avec les flux acceptés. Cela vérifie aussi que l'historique provisoire du premier essai ne fuit pas.
Reprise DAT bornée
Le test neuf sauvegarde U/p/k/omega provenant d'un état réellement accepté, restaure les trois cellules via restore_sst_production_state_from_dat, reconstruit les dérivés et observe le premier prédicteur. La comparaison porte sur les champs primaires et mu dérivé du k/omega restauré.
DAT v4 ne stocke pas le flux de face autoritaire. Le présent test ne certifie donc pas une équivalence de trajectoire couplée générale ni la conservation d'un historique de flux/temps à travers une interruption. La reconstruction de flux à partir d'une vitesse uniforme dans ce fixture ne suffit pas à résoudre ce manque architectural.
Exécutions et état des preuves
Toutes les exécutions mentionnées ici sont locales, avec GCC 13.3, C++17. Aucune CI distante n'a été lancée par cet agent.
Exécution
Observation
Journal
Ancien test SST compilé après premiers correctifs
39/40 au premier essai : le seul échec attendait la convergence d'un champ strain non scalaire. Attente migrée vers rejet précoce sans mutation.
integration_v2/sst_run_agent.log
Premier build commun
test_sst_blending passe ; un premier snapshot a donné CTest 41/41, observé pendant la session. Le journal commun est réécrit par les exécutions ultérieures : seule sa dernière version couvre l'arbre final.
coordination root
Cible dédiée initiale
5/5 : BC, gradient, résidu, prédicteur et retry.
integration_v2/sst_gaps_first_run.log
Cible dédiée avec stress et bilan
8/8 ; mesures MMS, Uy et bilan ci-dessus.
integration_v2/sst_focused_test.log
Premier essai du cas DAT
8/9 ; défaut de cohérence entre tolérance linéaire relative et gate SST absolu identifié par continuation directe sans DAT.
integration_v2/sst_dat_failure.log, sst_restart_diagnosis.log
Probe avant/après résidu après plancher
Ancien faux positif reproduit, nouveau rejet correct.
integration_v2/sst_bound_before.log, sst_bound_after.log
Probe MMS indépendante du solveur stationnaire
Erreurs/ordres mesurés de N8 à N128.
integration_v2/stress_mms_probe.log
Derniers sources avec DAT, métrique stricte, préflight dynamique et traduction de seuil
9/9 en compilation/exécution indépendante g++ -O0, bibliothèque core du lot précédent. La clôture commune figure au chapitre 8 : 43/43 Release et 9/9 ASan/UBSan.
integration_v2/sst_final_agent_build.log, sst_final_agent_test.log
Commandes de la cible permanente, après configuration du build commun :
/workspace/scratch/cae91c5c59b2/audit_tools/cmake/data/bin/cmake --build /workspace/scratch/cae91c5c59b2/build_integration_v2 --target test_sst_integration_gaps test_sst_blending test_steady_incompressible_solver
/workspace/scratch/cae91c5c59b2/audit_tools/cmake/data/bin/ctest --test-dir /workspace/scratch/cae91c5c59b2/build_integration_v2 --output-on-failure -R '^(test_sst_integration_gaps|test_sst_blending|test_steady_incompressible_solver)$'
Le probe MMS se compile avec g++ -std=c++17 -O1 -I src -I tests .../stress_mms_probe.cpp. Le probe de transport SST nécessite également la bibliothèque libcfdx_core.a et les mêmes bibliothèques MPI/HDF5 que le build : une première tentative sans cette bibliothèque a échoué au lien sur NativeAMGPreconditioner, puis les compilations correctement liées ont réussi. Cet échec d'outillage ne constitue pas un échec numérique.
git diff --check a été exécuté sur les fichiers SST modifiés, sans erreur d'espacement. La compilation commune et son arbre exact sont consignés aux chapitres 1 et 8.
PR #1092 : interactions et choix
Le diff de #1092 comprend 16 fichiers et modifie notamment les critères de solveurs CG/BiCGStab/GMRES, les cas de stencil déficient et les politiques de frontière des gradients LS/WLS. La valeur par défaut de BoundaryGradientPolicy devient ZERO_GRADIENT_GHOST et des paramètres de BC sont ajoutés à certains appels.
Le présent gradient SST ne dépend pas de ce changement de défaut : il résout explicitement les mêmes BC que son transport et réemploie uniquement l'interpolation géométrique existante. Il évite donc d'obtenir accidentellement un nouveau modèle de frontière au gré de l'ordre d'intégration des PR. Les améliorations de Krylov de #1092 peuvent être utiles, mais exigent leur propre lot de certification et le rejeu des transports SST. Le PCD du parcours couplé utilise encore la viscosité moléculaire dans sa préconditionneuse, alors que la matrice de quantité de mouvement consomme mu_eff ; cela demeure un sujet de robustesse à qualifier, pas une preuve d'absence de mu_eff de l'opérateur.
Ce qui reste à accomplir avant une promotion RANS complète
1.	Qualification des deux transports SST. MMS indépendantes de k et omega, avec production limitée, F1/F2 variables, diffusion croisée et convection non nulle sur au moins trois maillages. Les oracles de fermeture et le gradient affine ne remplacent pas cette campagne.
2.	Cas physiques déclarés. Canal contre DNS Moser–Kim–Mansour, plaque plane puis NACA0012 ; y+, frottement, Cp avec convention de pression explicitée, raffinement/GCI et bilans indépendants.
3.	Frontières et maillages. Étude de mu de face aux murs, gradients de face, maillages skew/non orthogonaux, frontières spécifiques non couvertes, contrastes de viscosité et robustesse du terme transposé explicite.
4.	Persistance complète. Flux de face autoritaire, historique transitoire et métadonnées de pression/model/BC cohérents dans un checkpoint ; essai interrompu contre calcul continu à flux et coefficients non uniformes, avec échecs de reprise atomiques.
5.	Algorithmes et parallélisme. Qualification des variantes SIMPLE/SIMPLEC/PISO/PIMPLE/COUPLED effectivement promises, préconditionnement variable, tailles représentatives et MPI ; pas seulement le fixture SIMPLE série.
6.	Boucle thermique couplée. L'énergie peut consommer nu_t accepté et lambda_t=rho*cp*nu_t/Pr_t via le travail thermique, mais une vérification de passage de champ ne qualifie pas une boucle RANS–énergie–rayonnement entièrement convergée.
7.	Extrêmes et autres modèles. Le comportement de production omega quand le calcul de nu_t sous-flue au niveau sous-normal reste borné par la convention du lot #1094 ; les vieux probes extrêmes sont conservés. Les autres modèles RANS demandent une revue équivalente de bornage/résidus/BC avant activation.
8.	Driver et contrats publics. Le driver futur doit conserver explicitement les obligations wall omega, porter p* et sa conversion statique, rapporter les bornages, sauvegarder tous les états nécessaires et ne jamais promouvoir un modèle au seul motif que des tests unitaires passent.
Le présent lot réduit des défauts d'intégration observables et ajoute des preuves discriminantes. La promotion physique et l'ouverture du sélecteur RANS restent bloquées par les travaux ci-dessus.
7. Ce qui manque encore pour tout intégrer dans CFDX
Cette matrice décrit les travaux qui restent après les correctifs présents dans le patch. Elle précise la dépendance et le résultat observable qui permettra de fermer chaque chantier. Une case ouverte n’est pas présentée comme une implémentation fournie.
Priorité / lot
Code ou contrat à compléter
Dépendances
Critère de clôture
P0 — Cas multiphysique
Régions, affectation des matériaux, maillages, types d’interfaces et identifiants stables dans case_schema.h, schéma Python et HDF5
Définition du format commun
Aller-retour C++/Python/HDF5 sans perte ; cas incohérents rejetés avant résolution
P0 — RANS applicatif
Fabrique de SSTProductionState, contrôles, BC k/omega, distances de paroi, faces omega, ownership et choix explicite de pression p*
Contrat du cas et parois
Un vrai cas utilisateur pilote le runtime SST, ses callbacks, son retry et ses sorties ; absence de fallback laminaire
P0 — Chaleur turbulente
Pr_t déclaré et persisté ; nu_t accepté fourni à build_effective_thermal_conductivity, puis énergie et S2S
État SST accepté
La dernière viscosité acceptée alimente le bon pas énergie ; une réponse rejetée ne modifie pas le champ thermique accepté
P0 — Pilote CHT
Chargement de plusieurs régions, résolution fluides/solides, appariement et échange conservatif des interfaces
Schéma de régions
Cas utilisateur fluide–solide réellement exécuté, avec conservation régionale et interface ; pas seulement appel unitaire du noyau
P1 — Contact et lambda variable CHT
Propriétés par cellule côté interface, résistances de contact par interface et deux températures de surface
Pilote CHT
Résistance série exacte sur matériaux asymétriques et plusieurs contacts ; pas de double comptage de Rc
P1 — Transaction couplée
Critères U/p/k/omega/T/rayonnement et acceptation commune, rollback de tous les états et caches
RANS + énergie + CHT/S2S
Rejet injecté à chaque étape puis reprise identique à un essai propre, avec bilans et état final comparés
P1 — Checkpoint complet
Flux de faces autoritaire, champs par région, temps, historiques BDF, métadonnées du modèle et version
Transaction et format de cas
Continuation interrompue comparée au calcul continu sur plusieurs pas ; la reprise actuelle de T seule ne suffit pas
P1 — Temps physique
Pilote temporel, historiques et séparation temps physique/itérations externes
Checkpoint et opérateurs cohérents
MMS transitoire sur au moins trois pas, conservation et restart ; les exécutables actuels refusent ce régime
P1 — Rayonnement général
Raccordements production P1/DOM/Rosseland, parois et sources volumétriques ; géométrie S2S et sous-surfaces non isothermes
Énergie et matériaux
Slab/enclosure de référence, fermeture source-flux, visibilité testée, convergence angulaire/spatiale selon le modèle
P1 — Import et interface utilisateur
Alias d’énergie Fluent energy, déclaration explicite des températures initiales et matrice de capacités commune
Contrat du cas
Un cas importé ne sélectionne aucune physique ignorée ; GUI/CLI et exécutables annoncent les mêmes capacités
P2 — Thermophysique
Cp(T), lambda(T), viscosité(T), enthalpie et densité selon la formulation choisie
Fiche mathématique et contrat matériaux
Conservation avec propriétés variables et références indépendantes ; ne pas étendre artificiellement la formule à cp constant
P2 — Qualification
MMS k/omega, canal turbulent, plaque, profils/Nu/Cf puis configurations plus complexes
Raccordements précédents
Études de maillage, erreurs et incertitudes tracées ; critères du registre satisfaits sur le commit final
P2 — Parallèle et matériel
Flux, parois, interfaces, réductions, reprise et stress transposé en multi-rang ; GPU si revendiqué
Contrats distribués
Parité et convergence MPI/CPU démontrées, et campagnes GPU distinctes
Ordre d’intégration proposé
1.	Intégrer le patch correctif exact, avec la CI de son commit final, tout en conservant les gardes de production.
2.	Définir le contrat de cas et de régions. Le parseur, les adaptateurs et les exécutables doivent partager le même sens des champs.
3.	Raccorder le pilote SST au cas, en conservant le cycle accepter/rejeter et en spécifiant p*. Qualifier les parois et les deux transports avant d’annoncer RANS utilisable.
4.	Raccorder la conductivité turbulente à l’énergie, puis au rayonnement. Commencer par une énergie passive sur écoulement accepté ; activer les rétroactions uniquement avec leurs critères couplés.
5.	Construire le pilote CHT, puis traiter contact par interface, propriétés variables et transactions communes.
6.	Compléter la reprise et le temps physique, puis les modèles radiatifs volumiques et les campagnes distribuées/physiques.
Ce découpage évite une PR unique impossible à attribuer : chaque lot doit contenir son contrat de cas, son raccordement, son contre-exemple et son oracle. Les changements de Krylov de #1092 et du porteur couplé de #1078 doivent être intégrés dans des étapes séparées, avec leurs propres preuves, avant d’affirmer la validité de leur combinaison avec ce lot.
Qualification minimale à ne pas confondre avec les tests présents
La MMS du stress transposé prouve un opérateur précis. Il faut encore une MMS des deux équations SST, au moins trois maillages, et des références de canal, paroi et transfert thermique. Un petit cas sur trois cellules prouve que les bons champs circulent entre composants ; il ne fournit ni un profil de loi de paroi, ni un y+, ni un Nusselt industriel validé.
Le test d’import de gravité ne qualifie pas tout le parcours Fluent : l’alias energy_model="energy" reste à router explicitement, et la température initiale Python par défaut de 288,15 K doit être distinguée d’une valeur déclarée. La nouvelle porte T>0 couvre le stage énergie monofluide et sa propagation dans les couplages ; elle ne démontre pas l’admissibilité générale de toutes les sorties du solveur solide.
La pression du runtime SST contient le terme isotrope turbulent : p*=p_statique+2*rho*k/3. Les BC et la conversion des résultats doivent respecter ce contrat avant de publier une pression statique ou un Cp. Boussinesq reste limité par l’hypothèse beta*DeltaT << 1 ; une variation de densité générale demande une autre formulation.
Enfin, un checkpoint qui contient U/p/T/k/omega mais pas le flux de face autoritaire ni les historiques n’est pas une preuve de reprise exacte d’une trajectoire couplée. Le présent lot améliore les champs sauvegardés et teste des reprises ciblées ; le critère de trajectoire reste ouvert.
8. Tests exécutés, résultats négatifs conservés et reproduction
Environnement et niveaux de preuve
Les calculs sont locaux sur CPU, dans un processus par exécutable : GCC 13.3, C++17, CMake 4.4.4, Ninja 1.13.2, HDF5 1.10.10 et OpenMPI 4.1.6. Les cibles sont compilées avec MPI disponible, mais aucun lancement distribué multi-rang n’est revendiqué. La configuration DebugSanitizers du dépôt active AddressSanitizer et UndefinedBehaviorSanitizer ; les tests ciblés n’ajoutent aucune suppression locale de diagnostic. Les tests longs et GPU ne font pas partie de cette campagne.
Le Python de validation fournit h5py 3.16.0, NumPy 2.5.4 et pytest 9.1.1, ainsi que les dépendances du schéma. Les quatre tests PyFluent ignorés sont signalés explicitement par pytest ; leur non-exécution n’est pas comptée comme une réussite.
Ensemble
Résultat
Ce que le résultat prouve
Base historique initiale
31/31 CTest
Les anciennes régressions passaient avant les nouveaux contre-exemples
Lot historique v1
33/33 CTest
Vérification du premier lot à arbre 407288b3…
Lot v2 consolidé
43/43 CTest
Non-régression ciblée et nouvelles interfaces sur l’arbre courant
ASan/UBSan ciblés v2
9/9 CTest
Aucun diagnostic mémoire/UB détecté sur ces chemins et ces cas
Python existant : schéma, Fluent, convertisseur
87 réussites, 4 ignorés
Compatibilité des contrats couverts ; PyFluent externe non exécuté
Campagne physique T2 complète
Non relancée
Les tests nouveaux portent sur le verdict, pas sur une nouvelle qualification Nusselt
CI du nouveau commit global
Non exécutée
Le lot n’a pas été poussé ; les contrôles obligatoires du dépôt restent à lancer sur le commit proposé
Les nombres de cibles CTest ne sont pas les nombres de cas internes : une cible C++ exécute plusieurs run_case, et une cible Python plusieurs tests. Les totaux ne doivent donc pas être additionnés comme des expériences indépendantes.
Échecs rencontrés et traitement
Les preuves négatives sont conservées. Le premier essai de 43 cibles avait deux erreurs de lancement dues au bit exécutable de deux binaires construits, un fixture de reprise thermique sans la sélection numérique exigée, et le défaut SST chaud décrit au chapitre 6. Les deux binaires ont été rendus exécutables, le fixture a reçu une configuration canonique complète, et la précision linéaire SST a été corrigée sans affaiblir son critère externe.
La toute première fixture solide manquait les listes de faces des patches HDF5 : cet échec de fixture a été séparé de la comparaison avant/après. La reproduction valide sur le code antérieur montrait six déclarations incompatibles silencieusement acceptées, tandis que le cas de conduction autorisé passait. Les nouveaux tests couvrent en plus la sélection temporelle canonique.
La revue croisée a reproduit un résumé converged=true après une exception énergie parce que le seul porteur avait convergé. Le catch publie maintenant converged=false et la cause réelle ; la voie normale indique aussi quelle porte énergie/rayonnement/flottabilité a échoué. Les tests d’application vérifient cette différence.
Un test nouveau a enfin reproduit une température négative avec résidus et bilans nuls. La porte physique la refuse, et le stage S2S retourne ce diagnostic structuré si sa seconde passe énergie échoue au lieu de lancer une reconstruction radiative sur un état inadmissible. Cette preuve montre pourquoi un bilan global fermé est nécessaire mais insuffisant.
La première campagne ASan/UBSan s’est terminée avec neuf erreurs de LeakSanitizer liées à l’accès /proc/ptrace de l’environnement restreint, après exécution des assertions. Elle n’est pas comptée comme un succès. La campagne finale a été relancée dans l’environnement permettant cet accès, sans désactiver LeakSanitizer ; ses journaux déterminent le verdict publié.
Les anciens échecs CMake de #1096 et compilation de #1099 restent historiques. Pour #1096, le HEAD récent a été configuré avec succès : le dossier ne maintient pas un blocage désormais corrigé. Pour #1099, le HEAD contrôlé est inchangé.
Détail des cibles Release
Cible
Résultat
Cas internes repérés dans le journal
test_audit_cht_time
PASS
—
test_audit_energy_boundary_metric
PASS
—
test_phase12_radiation_vv
PASS
3
test_dat_restart_identity
PASS
10
test_m1_m4_physics
PASS
8
test_finite_volume_transport
PASS
21
test_steady_incompressible_solver
PASS
24
test_m2_m4_solver_validation
PASS
10
test_phase10_turbulence_hardening
PASS
17
test_turbulence_qualification_equations
PASS
—
test_production_turbulence_contract
PASS
5
test_solid_conduction_selection
PASS
5
test_energy_integration_contracts
PASS
11
test_cht_integration_contracts
PASS
8
test_sst_integration_gaps
PASS
9
test_production_radiation_contract
PASS
6
test_phase10_sst_limiter
PASS
—
test_benchmark_matrix
PASS
38
test_level_b_reference_benchmarks
PASS
—
test_level_c_coupled_verification
PASS
—
test_cht_validation
PASS
5
test_thermal_vv
PASS
3
test_thermal_zones_vv
PASS
9
test_production_energy_stage_vv
PASS
5
test_production_s2s_energy_stage_vv
PASS
6
test_production_s2s_energy_outer_loop_vv
PASS
4
test_production_buoyancy_outer_loop_vv
PASS
—
test_radiation_vv
PASS
6
test_s2s_radiation_vv
PASS
8
test_thermal_radiation_model_matrix
PASS
7
test_thermal_radiation_regression
PASS
16
test_thermal_challenge_vv
PASS
—
test_thermophysical_models_vv
PASS
4
test_sst_blending
PASS
40
test_thermo_transport_contracts
PASS
6
test_thermophysical_models
PASS
5
test_hdf5_roundtrip
PASS
18
test_field_roundtrip
PASS
1
test_case_hdf5_io
PASS
22
test_solid_conduction_contract_e2e
PASS
8
test_t2_campaign_verdict
PASS
6
test_production_thermal_restart_e2e
PASS
6
test_imported_thermal_contracts
PASS
3
Cibles supplémentaires ASan/UBSan
test_solid_conduction_selection, test_energy_integration_contracts, test_cht_integration_contracts, test_sst_integration_gaps, test_production_s2s_energy_stage_vv, test_production_s2s_energy_outer_loop_vv, test_radiation_vv, test_s2s_radiation_vv, test_case_hdf5_io.
Reproduire dans un checkout propre
Le fichier contient un patch complet. Ne pas appliquer successivement le patch v1 et le patch v2 : le second inclut déjà le premier et s’applique directement au master indiqué. Utiliser un checkout propre ou un worktree isolé, puis extraire le patch sans traitement de texte qui modifierait les espaces ou les fins de ligne.
from pathlib import Path
import hashlib

report = Path("CFDX_integration_thermique_rayonnement_turbulence_2026-10-24.md").read_text()
patch = report.split("<!-- CFDX_PATCH_BEGIN -->\n````diff\n", 1)[1]
patch = patch.split("\n````\n<!-- CFDX_PATCH_END -->", 1)[0] + "\n"
assert hashlib.sha256(patch.encode()).hexdigest() == "1cef46e6ad105f0b1818bceeffeff313f5af1bacb8af4546b6a6fa9dbde1ebf3"
Path("candidate_integration_v2.patch").write_text(patch)
git worktree add --detach ../cfdx-integration-review 038e3aeff40f02af95b8845b931245ce23a1a094
cd ../cfdx-integration-review
git apply --check /chemin/candidate_integration_v2.patch
git apply /chemin/candidate_integration_v2.patch
git add -A
test "$(git write-tree)" = "89adbb733827bdc932ff546ecf701eec11fe65f6"

# Utiliser un Python où les dépendances du projet et h5py sont installées.
cmake -S . -B build-integration -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCFDX_BUILD_TESTS=ON \
  -DCFDX_ENABLE_LONG_VALIDATION=OFF -DCFDX_ENABLE_MPI=ON \
  -DCFDX_USE_HDF5_SERIAL=OFF -DPython3_EXECUTABLE=/chemin/venv/bin/python3
cmake --build build-integration --parallel 4 --target \
  test_m1_m4_physics \
  test_m2_m4_solver_validation \
  test_finite_volume_transport \
  test_turbulence_qualification_equations \
  test_production_turbulence_contract \
  test_production_radiation_contract \
  test_sst_blending \
  test_phase10_sst_limiter \
  test_phase10_turbulence_hardening \
  test_steady_incompressible_solver \
  test_cht_validation \
  test_thermal_vv \
  test_thermal_zones_vv \
  test_production_energy_stage_vv \
  test_production_s2s_energy_stage_vv \
  test_production_s2s_energy_outer_loop_vv \
  test_production_buoyancy_outer_loop_vv \
  test_radiation_vv \
  test_s2s_radiation_vv \
  test_thermal_radiation_model_matrix \
  test_thermal_radiation_regression \
  test_thermal_challenge_vv \
  test_thermophysical_models_vv \
  test_phase12_radiation_vv \
  test_thermo_transport_contracts \
  test_thermophysical_models \
  test_dat_restart_identity \
  test_case_hdf5_io \
  test_benchmark_matrix \
  test_level_b_reference_benchmarks \
  test_level_c_coupled_verification \
  test_audit_cht_time \
  test_audit_energy_boundary_metric \
  test_solid_conduction_selection \
  test_energy_integration_contracts \
  test_cht_integration_contracts \
  test_sst_integration_gaps \
  test_hdf5_roundtrip \
  test_field_roundtrip \
  cfdx_thermal_solver \
  cfdx_production_solver
ctest --test-dir build-integration --output-on-failure --no-tests=error \
  --parallel 3 --timeout 180 -R '^(test_m1_m4_physics|test_m2_m4_solver_validation|test_finite_volume_transport|test_turbulence_qualification_equations|test_production_turbulence_contract|test_production_radiation_contract|test_sst_blending|test_phase10_sst_limiter|test_phase10_turbulence_hardening|test_steady_incompressible_solver|test_cht_validation|test_thermal_vv|test_thermal_zones_vv|test_production_energy_stage_vv|test_production_s2s_energy_stage_vv|test_production_s2s_energy_outer_loop_vv|test_production_buoyancy_outer_loop_vv|test_radiation_vv|test_s2s_radiation_vv|test_thermal_radiation_model_matrix|test_thermal_radiation_regression|test_thermal_challenge_vv|test_thermophysical_models_vv|test_phase12_radiation_vv|test_thermo_transport_contracts|test_thermophysical_models|test_dat_restart_identity|test_case_hdf5_io|test_benchmark_matrix|test_level_b_reference_benchmarks|test_level_c_coupled_verification|test_audit_cht_time|test_audit_energy_boundary_metric|test_solid_conduction_selection|test_energy_integration_contracts|test_cht_integration_contracts|test_sst_integration_gaps|test_solid_conduction_contract_e2e|test_t2_campaign_verdict|test_hdf5_roundtrip|test_field_roundtrip|test_production_thermal_restart_e2e|test_imported_thermal_contracts)$'
PYTHONPATH=src/cfdx/python /chemin/venv/bin/python3 -m pytest -q -rs \
  tests/python/test_schema.py tests/python/test_fluent_adapter.py \
  tests/python/test_converter.py
Pour les sanitizers, utiliser un répertoire de build distinct, le type DebugSanitizers prévu par CMake et les neuf cibles listées ci-dessus. Les commandes absolues réellement exécutées, leurs codes de retour et leurs durées sont enregistrés dans l’annexe B.
Après toute modification de la base, résolution de conflit ou modification du patch, l’identité de l’arbre change : il faut recompiler et requalifier les chemins affectés. Le succès local ne remplace pas les checks obligatoires de la PR finale.
9. Inventaire exact du code à intégrer
Le tableau est généré depuis le diff contre la base. A signifie nouveau fichier ; M signifie fichier existant modifié. Les colonnes +/− comptent les lignes du patch, documentation et tests compris.
État
Fichier
+
−
M
CMakeLists.txt
38
0
M
apps/production_solver/main.cpp
70
22
M
apps/thermal_solver/main.cpp
2
8
M
docs/audits/09-turbulence-rans.md
119
4
M
docs/audits/10-thermique-cht.md
85
2
M
docs/audits/11-rayonnement.md
28
18
M
scripts/t2_heated_channel_calibration.py
48
12
M
src/cfdx/core/linalg/solver_status.h
8
0
M
src/cfdx/io/hdf5/case_hdf5_io.cpp
10
6
M
src/cfdx/physics/cht_solver.h
113
3
M
src/cfdx/physics/energy_solver.h
104
88
M
src/cfdx/physics/finite_volume_transport.h
86
13
M
src/cfdx/physics/m1_m4_models.h
22
6
A
src/cfdx/physics/momentum_viscous_stress.h
111
0
M
src/cfdx/physics/production_buoyancy_stage.h
22
5
M
src/cfdx/physics/production_energy_stage.h
122
16
M
src/cfdx/physics/production_radiation_s2s_stage.h
115
25
A
src/cfdx/physics/production_solid_conduction_contract.h
58
0
A
src/cfdx/physics/production_time_contract.h
34
0
M
src/cfdx/physics/radiation_advanced.h
94
37
A
src/cfdx/physics/scalar_transport_gradient.h
114
0
M
src/cfdx/physics/sst_production_coupling.h
52
14
M
src/cfdx/physics/sst_solver.h
119
9
M
src/cfdx/physics/steady_incompressible_solver.h
73
4
M
src/cfdx/physics/thermal_zones.h
16
2
M
src/cfdx/physics/turbulence.h
6
8
M
src/cfdx/python/cfdx/io/adapters/fluent.py
2
2
M
src/cfdx/python/cfdx/io/schema.py
8
6
A
tests/audit/chain_fixture.h
73
0
A
tests/audit/sst_channel_fixture.h
60
0
A
tests/audit/test_audit_cht_time.cpp
54
0
A
tests/audit/test_audit_energy_boundary_metric.cpp
38
0
A
tests/audit/test_sst_integration_gaps.cpp
382
0
M
tests/common/test_harness.h
10
2
A
tests/python/test_imported_thermal_contracts.py
36
0
A
tests/python/test_production_thermal_restart_e2e.py
157
0
A
tests/python/test_solid_conduction_contract_e2e.py
115
0
A
tests/python/test_t2_campaign_verdict.py
69
0
M
tests/unit/test_case_hdf5_io.cpp
37
0
M
tests/unit/test_finite_volume_transport.cpp
153
0
M
tests/unit/test_m1_m4_physics.cpp
20
1
A
tests/unit/test_solid_conduction_selection.cpp
104
0
M
tests/unit/test_sst_blending.cpp
11
11
M
tests/unit/test_turbulence_qualification_equations.cpp
61
1
M
tests/validation/test_benchmark_matrix.cpp
2
2
A
tests/validation/test_cht_integration_contracts.cpp
103
0
M
tests/validation/test_cht_validation.cpp
57
6
A
tests/validation/test_energy_integration_contracts.cpp
399
0
M
tests/validation/test_level_b_reference_benchmarks.cpp
1
1
M
tests/validation/test_level_c_coupled_verification.cpp
8
0
M
tests/validation/test_phase10_sst_limiter.cpp
3
1
M
tests/validation/test_production_s2s_energy_outer_loop_vv.cpp
88
8
M
tests/validation/test_production_s2s_energy_stage_vv.cpp
138
8
M
tests/validation/test_radiation_vv.cpp
35
0
M
tests/validation/test_s2s_radiation_vv.cpp
131
0
M
tests/validation/test_thermal_challenge_vv.cpp
1
1
M
tests/validation/test_thermal_radiation_regression.cpp
5
3
M
tests/validation/test_thermal_zones_vv.cpp
22
0
Les nouveaux headers sont inclus par leurs consommateurs ; CMake enregistre les nouvelles cibles de tests. Le patch met aussi à jour les fiches 09, 10 et 11. Les scripts d’audit et journaux situés hors du dépôt servent de preuve : ils ne doivent pas être ajoutés au runtime de CFDX.
Sources et traçabilité
Les preuves principales sont le code et les tests exacts de l’annexe A, les journaux de l’annexe B, les versions de PR du chapitre 2 et l’historique de l’annexe C. Les notes détaillées citent chaque symbole concerné et décrivent son oracle. Les références de formulation SST et rayonnement déjà utilisées dans la première revue sont conservées dans les fiches de méthode du dépôt ; elles ne sont pas remplacées par les seuls résultats de tests.
Annexe A — Patch complet du lot v2 vérifié
Base 038e3aeff40f02af95b8845b931245ce23a1a094 ; arbre reconstruit attendu 89adbb733827bdc932ff546ecf701eec11fe65f6 ; SHA-256 1cef46e6ad105f0b1818bceeffeff313f5af1bacb8af4546b6a6fa9dbde1ebf3. Le payload ci-dessous est intégral, sans ellipse.
diff --git a/CMakeLists.txt b/CMakeLists.txt
index 9fd63523..b0e43e18 100644
--- a/CMakeLists.txt
+++ b/CMakeLists.txt
@@ -419,6 +419,13 @@ if(CFDX_BUILD_TESTS)
         add_test(NAME ${test_name} COMMAND ${test_name})
     endmacro()
 
+    # Focused thermal audit regressions, extracted from PR #1096 without its
+    # unrelated changes or truncated build configuration.
+    add_cfdx_test(test_audit_cht_time tests/audit/test_audit_cht_time.cpp)
+    add_cfdx_test(test_audit_energy_boundary_metric tests/audit/test_audit_energy_boundary_metric.cpp)
+    set_tests_properties(test_audit_cht_time test_audit_energy_boundary_metric
+        PROPERTIES LABELS "audit;regression;thermal" TIMEOUT 60)
+
     # N8 physical/production qualification reuses the existing validation
     # executables and records a machine-readable campaign report. It is
     # intentionally a separate CTest entry so the individual physical or
@@ -655,6 +662,10 @@ set_tests_properties(test_phase12_radiation_vv PROPERTIES LABELS "phase12;valida
     add_cfdx_test(test_phase10_turbulence_hardening tests/validation/test_phase10_turbulence_hardening.cpp)
     add_cfdx_test(test_turbulence_qualification_equations tests/unit/test_turbulence_qualification_equations.cpp)
     add_cfdx_test(test_production_turbulence_contract tests/unit/test_production_turbulence_contract.cpp)
+    add_cfdx_test(test_solid_conduction_selection tests/unit/test_solid_conduction_selection.cpp)
+    add_cfdx_test(test_energy_integration_contracts tests/validation/test_energy_integration_contracts.cpp)
+    add_cfdx_test(test_cht_integration_contracts tests/validation/test_cht_integration_contracts.cpp)
+    add_cfdx_test(test_sst_integration_gaps tests/audit/test_sst_integration_gaps.cpp)
     add_cfdx_test(test_production_radiation_contract tests/unit/test_production_radiation_contract.cpp)
     add_cfdx_test(test_n11_sa_negative_nu_tilde_policy tests/unit/test_n11_sa_negative_nu_tilde_policy.cpp)
     add_cfdx_test(test_n11_bound_reporting_hardening tests/unit/test_n11_bound_reporting_hardening.cpp)
@@ -1124,6 +1135,33 @@ add_cfdx_test(test_n11_fw2_positivity_clamp_counter tests/unit/test_n11_fw2_posi
             WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
     endif()
     if(Python3_Interpreter_FOUND AND CFDX_HAS_HDF5)
+        add_test(NAME test_solid_conduction_contract_e2e
+            COMMAND ${Python3_EXECUTABLE}
+                    ${CMAKE_SOURCE_DIR}/tests/python/test_solid_conduction_contract_e2e.py)
+        set_tests_properties(test_solid_conduction_contract_e2e PROPERTIES
+            LABELS "thermal;integration;validation"
+            ENVIRONMENT "CFDX_THERMAL_SOLVER=$<TARGET_FILE:cfdx_thermal_solver>"
+            TIMEOUT 120)
+        add_test(NAME test_t2_campaign_verdict
+            COMMAND ${Python3_EXECUTABLE}
+                    ${CMAKE_SOURCE_DIR}/tests/python/test_t2_campaign_verdict.py)
+        set_tests_properties(test_t2_campaign_verdict PROPERTIES
+            LABELS "thermal;qualification;validation"
+            TIMEOUT 30)
+        add_test(NAME test_production_thermal_restart_e2e
+            COMMAND ${Python3_EXECUTABLE}
+                    ${CMAKE_SOURCE_DIR}/tests/python/test_production_thermal_restart_e2e.py)
+        set_tests_properties(test_production_thermal_restart_e2e PROPERTIES
+            LABELS "thermal;integration;restart"
+            ENVIRONMENT "CFDX_PRODUCTION_SOLVER=$<TARGET_FILE:cfdx_production_solver>"
+            TIMEOUT 120)
+        add_test(NAME test_imported_thermal_contracts
+            COMMAND ${Python3_EXECUTABLE}
+                    ${CMAKE_SOURCE_DIR}/tests/python/test_imported_thermal_contracts.py)
+        set_tests_properties(test_imported_thermal_contracts PROPERTIES
+            LABELS "thermal;integration;io"
+            ENVIRONMENT "PYTHONPATH=${CMAKE_SOURCE_DIR}/src/cfdx/python"
+            TIMEOUT 30)
         add_test(NAME test_n9_coupled_schur_case_selection
             COMMAND ${Python3_EXECUTABLE}
                     ${CMAKE_SOURCE_DIR}/tests/python/test_n9_coupled_schur_case_selection.py)
diff --git a/apps/production_solver/main.cpp b/apps/production_solver/main.cpp
index 73cf5a8c..cbd79d30 100644
--- a/apps/production_solver/main.cpp
+++ b/apps/production_solver/main.cpp
@@ -13,6 +13,7 @@
 #include "cfdx/io/runtime/convergence_history.h"
 #include "cfdx/io/runtime/execution_summary.h"
 #include "cfdx/io/vtu/vtu_writer.h"
+#include "cfdx/physics/production_time_contract.h"
 
 #include <algorithm>
 #include <cmath>
@@ -337,6 +338,10 @@ std::string energy_convergence_json(
         << (energy.global_energy_balance_passed ? "true" : "false") << ",\n";
     out << "  \"local_energy_balance_passed\": "
         << (energy.local_energy_balance_passed ? "true" : "false") << ",\n";
+    out << "  \"temperature_admissibility_passed\": "
+        << (energy.temperature_admissibility_passed ? "true" : "false") << ",\n";
+    out << "  \"nonpositive_temperature\": " << energy.nonpositive_temperature
+        << ",\n";
     out << "  \"iterations\": " << energy.iterations << ",\n";
     out << "  \"independent_conservation_residual\": "
         << json_number(energy.independent_conservation_residual) << ",\n";
@@ -396,6 +401,10 @@ std::string energy_summary_json(const ProductionEnergyStageResult& energy)
         << json_number(energy.boundary_flux_out) << ",\n";
     out << "  \"nonfinite_temperature\": " << energy.nonfinite_temperature
         << ",\n";
+    out << "  \"nonpositive_temperature\": " << energy.nonpositive_temperature
+        << ",\n";
+    out << "  \"temperature_admissibility_passed\": "
+        << (energy.temperature_admissibility_passed ? "true" : "false") << ",\n";
     out << "  \"volumetric_power\": "
         << json_number(energy.volumetric_power) << ",\n";
     out << "  \"patches\": [\n";
@@ -478,7 +487,9 @@ std::string radiation_s2s_json(const ProductionS2SOuterLoopResult& stage)
             << ", \"net_power\": "
             << json_number(radiation.net_power[index])
             << ", \"injected_heat_flux\": "
-            << json_number(radiation.injected_heat_flux.at(name)) << "}";
+            << json_number(radiation.injected_heat_flux.at(name))
+            << ", \"applied_heat_flux\": "
+            << json_number(radiation.applied_heat_flux.at(name)) << "}";
         out << (index + 1 == radiation.surface_names.size() ? "\n"
                                                             : ",\n");
     }
@@ -577,6 +588,7 @@ int main(int argc, char** argv)
         // load time; the flow-energy buoyancy outer loop below is the only
         // solver path that may run it.
         if (is_case_hdf5) {
+            validate_production_steady_time_selection(case_setup, "production solver");
             validate_production_turbulence_selection(case_setup);
             validate_production_energy_model(case_setup);
             validate_production_radiation_selection(case_setup);
@@ -879,6 +891,27 @@ int main(int argc, char** argv)
                 return g_stop_requested == 0;
             };
 
+        const bool energy_requested = is_case_hdf5 &&
+            production_energy_requested(case_setup);
+        Field<double, Location::CELL> temperature(
+            mesh.n_cells(), "temperature", "K", 1);
+        const bool thermal_restart = energy_requested && !options.restart.empty();
+        if (thermal_restart) {
+            // Validate the complete thermal checkpoint before running the
+            // carrier. The carrier's own initialization subsequently consumes
+            // U/p and its iteration metadata; no live fields are changed here.
+            auto candidate_velocity = U;
+            auto candidate_pressure = p;
+            DatRestartFields fields;
+            fields.temperature = &temperature;
+            read_dat_restart_fields(options.restart.string(), mesh,
+                candidate_velocity, candidate_pressure, fields);
+            for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
+                if (!std::isfinite(temperature(cell)) || !(temperature(cell) > 0.0))
+                    throw std::invalid_argument(
+                        "production solver: restarted temperature must be finite and positive (kelvin)");
+        }
+
         IncompressibleSolveResult result;
         bool buoyancy_ran = false;
         ProductionBuoyancyOuterLoopResult buoyancy_result;
@@ -893,7 +926,7 @@ int main(int argc, char** argv)
             // never an intermediate Picard pass.
             buoyancy_result = run_production_buoyant_outer_loop(
                 mesh, ubc, pbc, controls, case_setup, U, p,
-                options.restart.string());
+                options.restart.string(), thermal_restart ? &temperature : nullptr);
             result = buoyancy_result.flow;
             buoyancy_ran = true;
             std::cout << "Buoyant flow-energy outer loop "
@@ -938,15 +971,10 @@ int main(int argc, char** argv)
         // diverged flow is not evidence, and the global exit code below
         // is flow AND energy. The energy model was already validated
         // fail-closed right after case load.
-        bool energy_requested = false;
         bool energy_ran = false;
         ProductionEnergyStageResult energy_result;
         ProductionS2SOuterLoopResult s2s_outer_result;
         bool s2s_ran = false;
-        Field<double, Location::CELL> temperature(
-            mesh.n_cells(), "temperature", "K", 1);
-        if (is_case_hdf5)
-            energy_requested = production_energy_requested(case_setup);
         if (buoyancy_ran) {
             // The buoyant outer loop already produced the coupled energy
             // response (and, when declared, the energy-S2S fixed point)
@@ -955,7 +983,6 @@ int main(int argc, char** argv)
             // responses; they are gated on the loop having completed at
             // least one full outer pass, exactly like the frozen-flow
             // path.
-            energy_requested = true;
             energy_ran = buoyancy_result.energy_ran;
             s2s_ran = buoyancy_result.s2s_ran;
             energy_result = buoyancy_result.energy;
@@ -967,8 +994,9 @@ int main(int argc, char** argv)
             // Initial temperature: the declared initial condition, else
             // the case reference temperature. Zero means "not declared"
             // (stage-1 convention) and falls back to ref_temperature.
-            double initial_temperature = case_setup.ref_temperature;
-            if (case_setup.has_initial_condition &&
+            double initial_temperature = thermal_restart
+                ? temperature(0) : case_setup.ref_temperature;
+            if (!thermal_restart && case_setup.has_initial_condition &&
                 case_setup.initial_condition.temperature != 0.0)
                 initial_temperature =
                     case_setup.initial_condition.temperature;
@@ -980,8 +1008,9 @@ int main(int argc, char** argv)
                 throw std::invalid_argument(
                     "production solver: initial temperature must be finite "
                     "and positive (kelvin)");
-            for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
-                temperature(cell) = initial_temperature;
+            if (!thermal_restart)
+                for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
+                    temperature(cell) = initial_temperature;
             if (production_s2s_radiation_requested(case_setup)) {
                 // External energy-radiation iteration (issue #885 step
                 // 3): energy pass on the frozen flow, S2S radiosity from
@@ -1053,6 +1082,10 @@ int main(int argc, char** argv)
                       << (energy_result.global_energy_balance_passed ? 1 : 0)
                       << " local_balance="
                       << (energy_result.local_energy_balance_passed ? 1 : 0)
+                      << " temperature_admissibility="
+                      << (energy_result.temperature_admissibility_passed ? 1 : 0)
+                      << " nonpositive_temperature="
+                      << energy_result.nonpositive_temperature
                       << " audit_scale=" << energy_result.audit_power_scale
                       << "\n";
             std::cout << "Energy temperatures ["
@@ -1108,17 +1141,36 @@ int main(int argc, char** argv)
               (!s2s_ran || s2s_outer_result.converged)));
         const auto solver_exit_code =
             execution_exit_code_for_solver(overall_converged);
+        auto overall_status = result.convergence_status;
+        std::string overall_reason = result.convergence_reason;
+        if (result.converged && !overall_converged) {
+            overall_status = cfdx::core::ConvergenceStatus::DIVERGED;
+            if (energy_ran && !energy_result.temperature_admissibility_passed)
+                overall_reason = "energy temperature admissibility gate failed";
+            else if (energy_ran && !energy_result.converged)
+                overall_reason = "energy stage linear or conservation gate failed";
+            else if (s2s_ran && !s2s_outer_result.converged)
+                overall_reason = "energy-radiation fixed point did not converge";
+            else if (buoyancy_ran && !buoyancy_result.converged)
+                overall_reason = "flow-energy buoyancy fixed point did not converge";
+            else
+                overall_reason = "requested energy stage did not complete";
+        }
 
         const auto dat = options.output_dir / "restart.dat";
-        write_dat_restart(dat.string(), mesh, U, p, result.iterations, 0.0);
+        DatRestartFields checkpoint_fields;
+        if (energy_ran)
+            checkpoint_fields.temperature = &temperature;
+        write_dat_restart_fields(dat.string(), mesh, U, p, checkpoint_fields,
+                                 result.iterations, 0.0);
         const auto convergence = options.output_dir / "convergence.json";
         write_convergence_history(convergence.string(), result);
         write_execution_summary(execution.string(), ExecutionSummary{
             solver_exit_code,
             overall_converged,
             result.iterations,
-            static_cast<int>(result.convergence_status),
-            result.convergence_reason,
+            static_cast<int>(overall_status),
+            overall_reason,
             std::filesystem::is_regular_file(dat),
             std::filesystem::is_regular_file(convergence)});
         std::cout << "Checkpoint " << dat.string() << "\n";
@@ -1145,14 +1197,10 @@ int main(int argc, char** argv)
                     std::filesystem::is_regular_file(execution_output_dir / "convergence.json");
                 write_execution_summary(execution.string(), ExecutionSummary{
                     kExecutionExitError,
-                    execution_result.has_value() ? execution_result->converged : false,
+                    false,
                     execution_result.has_value() ? execution_result->iterations : 0,
-                    execution_result.has_value()
-                        ? static_cast<int>(execution_result->convergence_status)
-                        : 0,
-                    execution_result.has_value()
-                        ? execution_result->convergence_reason
-                        : exc.what(),
+                    static_cast<int>(cfdx::core::ConvergenceStatus::DIVERGED),
+                    exc.what(),
                     has_restart,
                     has_convergence});
             } catch (const std::exception& summary_exc) {
diff --git a/apps/thermal_solver/main.cpp b/apps/thermal_solver/main.cpp
index 85fa5d06..1ec55412 100644
--- a/apps/thermal_solver/main.cpp
+++ b/apps/thermal_solver/main.cpp
@@ -36,6 +36,7 @@
 //     otherwise ref_temperature. The linear solve result does not depend
 //     on this initial guess.
 #include "cfdx/physics/thermal_zones.h"
+#include "cfdx/physics/production_solid_conduction_contract.h"
 #include "cfdx/io/hdf5/case_hdf5_io.h"
 #include "cfdx/io/runtime/execution_summary.h"
 #include "cfdx/io/vtu/vtu_writer.h"
@@ -477,14 +478,7 @@ int main(int argc, char** argv)
         if (mesh.n_cells() == 0)
             throw std::runtime_error("thermal solver: case carries no cells");
 
-        // Energy-model gate: the thermal path is only meaningful for a case
-        // that actually declares a thermal problem.
-        if (case_setup.energy_model != "conjugate" &&
-            case_setup.energy_model != "solid_conduction")
-            throw std::invalid_argument(
-                "thermal solver: case energy_model is '" +
-                case_setup.energy_model + "'; the thermal path requires "
-                "'conjugate' or 'solid_conduction'");
+        validate_production_solid_conduction_selection(case_setup);
 
         const auto geometry = build_fv_geometry(mesh);
         const auto specs = build_zone_specs(case_setup);
diff --git a/docs/audits/09-turbulence-rans.md b/docs/audits/09-turbulence-rans.md
index e2a188c0..121d41b6 100644
--- a/docs/audits/09-turbulence-rans.md
+++ b/docs/audits/09-turbulence-rans.md
@@ -12,7 +12,8 @@
 - W. P. Jones, B. E. Launder, « The prediction of laminarization with a two-equation model of turbulence », *Int. J. Heat Mass Transfer* 15:301-314, 1972. (k-epsilon original ; DOI non confirmé — citation bibliographique complète exigée)
 - B. E. Launder, D. B. Spalding, « The numerical computation of turbulent flows », *Comput. Methods Appl. Mech. Engrg.* 3(2):269-289, 1974 — DOI : 10.1016/0045-7825(74)90029-2. (constantes k-epsilon standard et lois de paroi)
 - D. C. Wilcox, « Reassessment of the scale-determining equation for advanced turbulence models », *AIAA Journal* 26(11):1299-1310, 1988 — DOI : 10.2514/3.10041. (k-omega)
-- F. R. Menter, « Two-equation eddy-viscosity turbulence models for engineering applications », *AIAA Journal* 32(8):1598-1605, 1994 — DOI : 10.2514/3.12149. (SST)
+- F. R. Menter, « Two-equation eddy-viscosity turbulence models for engineering applications », *AIAA Journal* 32(8):1598-1605, 1994 — DOI : 10.2514/3.12149. (SST original formulation).
+- NASA Turbulence Modeling Resource, « Menter Shear Stress Transport (SST) Model », SST-2003 formulation and production limiter: https://tmbwg.github.io/turbmodels/sst.html. The cell-source convention documented below uses the SST-2003 limited-production form; full-model qualification remains separate.
 - P. R. Spalart, S. R. Allmaras, « A One-Equation Turbulence Model for Aerodynamic Flows », *AIAA Paper* 92-0439, 30th Aerospace Sciences Meeting, 1992 — DOI : 10.2514/6.1992-439.
 - D. C. Wilcox, *Turbulence Modeling for CFD*, 3e éd., DCW Industries, 2006. (référence de synthèse)
 
@@ -29,8 +30,11 @@ avec $P_k=2\nu_t\,\mathbf{S}:\mathbf{S}$, $\mathbf{S}=\tfrac{1}{2}(\nabla\mathbf
 k-omega (Wilcox 1988) : $\nu_t=k/\omega$, même équation k avec $\sigma_k=0.5$, et
 $$\frac{\partial \omega}{\partial t}+\nabla\cdot(\mathbf{u}\omega)=\alpha\frac{\omega}{k}P_k-\beta\omega^2+\nabla\cdot\left[\left(\nu+\frac{\nu_t}{\sigma_\omega}\right)\nabla \omega\right]$$
 
-SST (Menter 1994) : mélange de k-omega près de la paroi et k-epsilon transformé à distance via $F_1$ ; limiteur :
-$$\nu_t=\frac{a_1 k}{\max(a_1\omega,\; S\,F_2)},\qquad a_1=0.31$$
+SST (variante SST-2003, Menter) : mélange de k-omega près de la paroi et k-epsilon transformé à distance via $F_1$ ; limiteur de viscosité :
+$\nu_t=\frac{a_1 k}{\max(a_1\omega,\; S\,F_2)},\qquad a_1=0.31$
+Production turbulente limitée :
+$\widetilde P_k=\min(\rho\nu_t S^2,\;10\beta^*\rho k\omega),\qquad P_\omega=\gamma\,\frac{\widetilde P_k}{\nu_t}$
+Le même $\widetilde P_k$ alimente la production de $k$ et celle de $\omega$ ; si $\nu_t=0$, la contribution de production de $\omega$ est définie comme nulle. Les termes de diffusion croisée restent ajoutés séparément. Cette convention doit être vérifiée contre les références SST-2003 avant qualification physique.
 
 Spalart-Allmaras (variable de travail $\tilde{\nu}$, $\nu_t=\tilde{\nu} f_{v1}$, $f_{v1}=\chi^3/(\chi^3+C_{v1}^3)$, $\chi=\tilde{\nu}/\nu$) :
 $$\frac{\partial \tilde{\nu}}{\partial t}+\nabla\cdot(\mathbf{u}\tilde{\nu})=c_{b1}\tilde{S}\tilde{\nu}+\nabla\cdot\left[(\nu+\tilde{\nu})\nabla\tilde{\nu}\right]+c_{b2}|\nabla\tilde{\nu}|^2-c_{w1}f_w\left(\frac{\tilde{\nu}}{d}\right)^2$$
@@ -98,8 +102,119 @@ $$a_P\phi_P=\sum_{nb}a_{nb}\phi_{nb}+b,\qquad a_P=\sum_{nb}a_{nb}+\frac{\rho V_P
 
 
 ## 12. Production outer-loop SST runtime (PR pending)
-- `IncompressibleSolverControls::sst_runtime` is an explicit opt-in runtime seam. When absent, the legacy fixed-viscosity path is unchanged.
+- `IncompressibleSolverControls::sst_runtime` is an explicit opt-in runtime seam. When absent and no cellwise turbulent-viscosity field is supplied, the legacy uniform-viscosity momentum operator is unchanged.
 - When present, the caller must provide a validated `SSTProductionState` and SST `TurbulenceTransportControls`; initialization from declared case data or DAT restart remains explicit through the existing initialization/restart helpers. Carrier and SST density/kinematic viscosity must agree exactly. No implicit SST initial values are fabricated.
 - Before the first momentum assembly, the solver refreshes wall distance, strain rate and `nu_t` from the actual initialized/restarted carrier field. The outer iteration then advances SST using the pressure-corrected authoritative face mass flux and the explicitly supplied omega face values when configured. A failed transport step raises the existing retryable-iteration path; the enclosing retry restores the SST snapshot together with carrier rollback.
 - After the carrier/SST transaction commits, the refreshed `nu_t` is converted to `mu_eff = rho * (nu + nu_t)` for the next momentum assembly. SST convergence is an additional hard gate based on the transport solver's convergence status, finite equation residuals, zero non-finite bound events, and dimensionless relative `nu_t` change no larger than the configured carrier relative tolerance. This prevents declaring the carrier solution converged before momentum has been assembled with the accepted SST closure. Raw dimensional SST residuals are not mixed into the normalized carrier convergence metric; carrier, continuity, and QoI gates are not relaxed.
 - This is integration scaffolding, not permission to enable the general RANS selector or a claim of physical qualification. Restart/state lifecycle and production regression coverage must be demonstrated by exact-head CI before merge.
+
+## 13. Fermeture des écarts d'intégration SST — 2026-10-24
+
+Cette section décrit le code du lot local d'intégration ; les constats des sections
+8 à 11 restent des états historiques. Le sélecteur RANS général reste fermé.
+
+### Frontières, gradients et état retourné
+
+`SSTProductionRuntime::wall_omega_sst_patches` conserve les obligations physiques
+`WALL_OMEGA_SST`, indépendamment des parois utilisées pour la distance minimale.
+Le contrôle `require_supported_sst_runtime_boundaries(contract, mesh, geometry, values)`
+parcourt les `face_ids` du maillage : chaque face de chaque paroi déclarée doit
+avoir une valeur de face globale, finie et strictement positive, avec une
+appartenance de patch cohérente. Ce contrôle précède l'initialisation de `U/p`.
+L'ancienne surcharge sans maillage refuse désormais un contrat avec wall omega :
+la présence d'une seule entrée finie ne permet pas de prouver la couverture.
+
+Les gradients utilisés par `F1` et par la diffusion croisée de l'équation omega
+emploient les mêmes états de frontière que le transport scalaire : valeur
+imposée, valeur de face prioritaire, gradient normal imposé, ou extrapolation
+de la valeur propriétaire si le gradient est nul. Le noyau commun est
+`scalar_gradient_with_transport_boundaries`. Il n'ajoute pas de correction de
+skewness. `F2` ne dépend pas de la diffusion croisée ; les reconstructions de
+`nu_t` seules n'ont donc pas besoin de calculer ce gradient.
+
+Après relaxation puis planchers de positivité, les vecteurs servant aux résidus
+sont remis en cohérence avec les champs effectivement retournés. La convergence
+conserve les contrôles existants et ajoute
+
+$$e_\phi=\max_i\frac{|(A\phi-b)_i|}{\sum_j|a_{ij}|}\leq\mathrm{tolerance},
+\qquad\phi\in\{k,\omega\}.$$
+
+Ce défaut d'équation a les unités du champ et utilise le même seuil absolu que
+le changement d'itéré ; il reste défini pour un second membre nul. Les résidus
+bruts dimensionnels et le nombre/amplitude des corrections de plancher restent
+rapportés. Aucun plancher ne peut masquer un défaut matériel de l'équation.
+
+La tolérance interne relative du solveur linéaire est traduite à partir de ce
+seuil absolu, séparément pour k et omega :
+
+$$\mathrm{tol}_{lin}=\mathrm{tolerance}\,\min\left(1,
+\frac{\min_i\sum_j|a_{ij}|}{\max(\|b\|_2,10^{-15})}\right).$$
+
+Cette traduction peut seulement resserrer le critère linéaire. Elle évite qu'un
+itéré chaud soit rendu inchangé parce que son résidu relatif au RHS satisfait
+le solveur linéaire alors que son défaut en unités du champ échoue à la porte
+SST. Les seuils externes, budgets et relaxations sont inchangés. Une norme non
+finie, une ligne nulle ou une tolérance résultante sous-débordée échouent
+explicitement. Le test de reprise DAT a révélé ce désaccord d'échelles ; la
+continuation directe du même état sans sérialisation le reproduisait aussi.
+
+### Tenseur visqueux et pression
+
+Pour le runtime SST ou un champ de viscosité turbulente cellulaire explicite,
+la quantité de mouvement incompressible emploie
+
+$$\nabla\cdot(\rho\mathbf U\otimes\mathbf U)
+=-\nabla p^*+\nabla\cdot\left[\mu_{eff}
+(\nabla\mathbf U+\nabla\mathbf U^T)\right]+\rho\mathbf f,
+\quad\mu_{eff}=\rho(\nu+\nu_t),\quad p^*=p_{statique}+\tfrac23\rho k.$$
+
+La diffusion de chaque composante demeure implicite. Le terme transposé
+complémentaire est un flux conservatif explicite de l'itération externe,
+
+$$F^T_{i,f}=\mu_f\sum_j(\partial_i U_j)_f S_{f,j}.$$
+
+Un même flux orienté propriétaire est ajouté/soustrait aux cellules adjacentes.
+La viscosité interne suit l'interpolation harmonique pondérée par les distances
+normales de la diffusion scalaire ; en frontière, la viscosité et le gradient
+de face utilisent la cellule propriétaire. Les valeurs de vitesse imposées
+entrent dans le gradient de Gauss. Le terme est réévalué avec `U` final pour le
+résidu non linéaire, et reconstruit comme flux dans le bilan indépendant.
+Le bilan additionne flux sortant et `grad(p*) - rho*f` : son signe et le facteur
+densité sont testés avec une force volumique non nulle. Le parcours historique
+à viscosité uniforme, sans runtime SST ni champ cellulaire, conserve sa
+réduction incompressible en laplacien ; l'égalité des assembleurs scalaires à
+viscosité uniforme ne signifie pas l'identité de toute itération avec un terme
+transposé discret explicite.
+
+Dans ce chemin SST, les valeurs initiales et les conditions de pression portent
+sur `p*`. Le noyau ne convertit pas automatiquement une pression statique
+imposée et ne certifie pas des sorties `Cp`. Toute utilisation comme pression
+statique doit soustraire explicitement `2*rho*k/3` et déclarer le même choix aux
+frontières.
+
+### Preuves bornées et limites
+
+`test_sst_integration_gaps` contient les oracles de couverture complète et
+d'atomicité avant initialisation, un gradient affine avec BC réelles, le
+contre-exemple de faux résidu après plancher, un premier prédicteur à trois
+cellules avec flux non nul, un rejet après commit SST suivi d'un retry, et une
+réhydratation DAT de l'état primaire suivie d'un premier prédicteur. Le fixture
+à parois mobiles vérifie un raccordement logiciel ; il ne valide pas un canal
+turbulent ni la convergence RANS globale.
+
+La MMS du terme transposé impose `U=(y,0,0)` et
+`mu=1+x^2*(1-x)^2`, donc une force exacte
+`(0,2*x*(1-x)*(1-2*x),0)` ; elle compare les erreurs sur 16, 32 et 64 cellules,
+le bilan des flux et la limite constante. Sa dérivée de viscosité s'annule aux
+frontières : elle ne qualifie pas une loi de paroi, des variations arbitraires
+de viscosité en frontière ni les maillages non orthogonaux. Le traitement
+explicite doit encore être étudié pour les forts contrastes et les difficultés
+de convergence des cas physiques.
+
+DAT v4 conserve `U/p/k/omega`, mais pas le flux de face autoritaire. Le test de
+réhydratation ne prouve donc pas la reprise d'une trajectoire couplée générale.
+Restent notamment la MMS des deux équations SST sur au moins trois maillages,
+les qualifications DNS canal, plaque plane puis NACA0012, les lois de paroi et
+`y+`, les reprises complètes avec flux/historique, et les parcours MPI. La
+qualification du système turbulence-énergie-rayonnement reste distincte des
+tests de raccordement des champs `nu_t` et `lambda_eff`.
diff --git a/docs/audits/10-thermique-cht.md b/docs/audits/10-thermique-cht.md
index 51172162..7d1c57ca 100644
--- a/docs/audits/10-thermique-cht.md
+++ b/docs/audits/10-thermique-cht.md
@@ -6,7 +6,7 @@
 - Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json + docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json
 - Issues liées : #34 (roadmap, Phase 11) ; docs/validation/THERMOPHYSICAL_MODELS.md
 - Dernier commit vérifié : f5a728f205cb45ad6186c9fab92554ae4af06bdd (master, 2026-10-12 — présence des sources et tests confirmée)
-- Limitations connues : bilans d'énergie d'interface non exécutés depuis l'environnement d'audit ; statut déclaratif jusqu'à campagne sur HEAD exact.
+- Limites de promotion : les preuves locales ajoutées au §9 couvrent des contrats ciblés ; elles ne qualifient pas le pilote CHT multi-régions de production, encore absent, ni les campagnes physiques et la parité MPI.
 
 ## 1. Références sources
 - S. V. Patankar, *Numerical Heat Transfer and Fluid Flow*, Hemisphere, 1980 — DOI : 10.1201/9781482234213. (couplage conjugué, moyenne harmonique, linéarisation)
@@ -28,6 +28,7 @@ Forme discrète sur la cellule $P$ :
 $$a_P T_P=\sum_{nb}a_{nb}T_{nb}+b,\qquad \sum_f \dot m_f c_p (\rho)_f \ldots$$
 - **Flux convectif** : porte $\rho c_p$ — leçon d'audit #419 **B2** : quand $c_p$ varie (ou par prudence systématique), le produit $\rho c_p$ doit être interpolé et appliqué au flux convectif ; jamais de flux convectif en $\rho$ seul additionné d'un $c_p$ séparé après coup. Le bilan d'énergie doit être vérifié exact sur ce point.
 - **Flux diffusif** : conductivité de face $k_f$ par **moyenne harmonique** des valeurs de cellules (Patankar) pour les coefficients discontinus à l'interface fluide-solide ; coefficients assemblés dans la matrice $A$.
+- **Contrat des propriétés CHT** : les conductivités déclarées dans `ChtInterfaceControls` doivent être identiques (à l'arrondi relatif près) aux conductivités utilisées par les deux `EnergySolverControls`. Une configuration contradictoire est rejetée avant résolution ; le bilan d'interface ne peut pas certifier une propriété différente de celle de l'opérateur énergie.
 - **Interface CHT conservatif** : le flux d'interface $q_{\Gamma}$ est calculé **une seule fois** et débité/credité symétriquement aux deux régions ($-q_{\Gamma}$ d'un côté, $+q_{\Gamma}$ de l'autre) — conservation machine de l'énergie d'interface ; jamais deux flux calculés indépendamment.
 - **Terme source** $S_T$ : linéarisé en $S_T=S^{exp}+S^{lin}T_P$ quand il dépend de T (source radiatif, voir fiche 11) ; signe vérifié par bilan.
 
@@ -39,7 +40,7 @@ $$a_P T_P=\sum_{nb}a_{nb}T_{nb}+b,\qquad \sum_f \dot m_f c_p (\rho)_f \ldots$$
 
 ## 5. Pièges connus
 - **B2** : $\rho c_p$ dans la convection d'énergie — pas $\rho$ seul.
-- Bilan d'énergie d'interface : l'énergie qui quitte le fluide doit entrer dans le solide à la tolérance machine — sinon le couplage fuit.
+- Bilan d'énergie d'interface : l'énergie qui quitte le fluide doit entrer dans le solide à la tolérance machine — sinon le couplage fuit. Les flux doivent être reconstruits avec les conductivités réellement assemblées dans les équations régionales ; le seul équilibre du calcul d'interface ne suffit pas si les contrats de propriétés divergent.
 - Boussinesq : seulement pour $\Delta T$ petits ; hors validité, utiliser la flottabilité pleine (variable $\rho$) et documenter le choix.
 - $c_p$ variable : la formulation température peut perdre la conservativité — la formulation enthalpie est l'alternative documentée.
 - Moyenne harmonique obligatoire pour $k$ à l'interface (l'arithmétique sur-lisse les discontinuités et casse le flux exact).
@@ -66,3 +67,85 @@ $$a_P T_P=\sum_{nb}a_{nb}T_{nb}+b,\qquad \sum_f \dot m_f c_p (\rho)_f \ldots$$
 - 2026-10-12 — approbation du contrat (PR d'approbation #861) ; l'audit du code contre cette fiche reste requis avant toute promotion (Règle d'Or 1).
 - 2026-10-21 — étape 4 du plan #885 : stage de production flottabilité `production_buoyancy_stage.h` + routage `buoyant` dans `cfdx_production_solver` — le contrat Boussinesq de §2 (source a = −β(T−T_ref)·g_vec per-cellule, part constante absorbée dans la pression modifiée, valide βΔT ≪ 1) est implémenté via le canal cellwise du porteur (fiche 04 §7) dans une boucle externe flow-énergie (ω = 0,25 sur la température vue par la source, budget 200 fail-closed, critères bloquants séparés sur la réponse non relaxée ≤ 1e-10 : stationnarité T/U/p normalisée ; verdict = point fixe couplé exact, INCONCLUSIVE n'est jamais PASS) ; oracles : équivalence bit-exacte avec le `body_force` uniforme équivalent, identité T = T_ref, signes, replay du point fixe, catalogue fail-closed (test `test_production_buoyancy_outer_loop_vv`) ; aucune tolérance/gate/oracle/budget de la fiche modifiée ; la flottabilité pleine ρ variable reste hors validité documentée (§5) et exigerait une nouvelle fiche.
 - 2026-10-21 — correctif de la PR #1016 (échecs A/C/D de `test_production_buoyancy_outer_loop_vv`) : convergence-at-entry ajoutée à la boucle externe flottabilité — acceptation du porteur warm classé `STAGNATED` au plancher machine (résidus pression/continuité/flux corrigé ≤ 1e-14, momentum RHS-normalisé ≤ 1e-12), restauration U/p d'entrée, conservation du porteur précédemment accepté, comptabilité `flow_converged_at_entry` ; cas A à la 3e itération (itérations 1-2 sans at-entry, itération 3 at-entry, porteur inchangé vs itération 2) ; échecs C/D relèvent de la calibration de la gate locale de l'étage énergie (fiche 14 §7/§8) ; aucune tolérance/gate/oracle/budget de la fiche modifiée.
+
+
+## 9. Intégration locale du 24 octobre 2026 — contrats vérifiés et limites
+
+Cette révision complète les corrections CHT/énergie du lot construit sur
+`038e3aeff40f02af95b8845b931245ce23a1a094`. Les sources et les tests sont
+présents dans le patch local de la révision ; aucun statut du registre de
+qualification n'est promu par ces modifications.
+
+### Formulation, propriétés et source physique
+
+`energy_solver.h` conserve une équation en température avec `rho` et `cp`
+constants. La source physique est `Su + Sp*T_final`, y compris dans le bilan
+indépendant global, le bilan par cellule et les diagnostics après relaxation.
+Le diagnostic ne réutilise pas `A*T-b` : `reconstruct_energy_face_flux`
+reconstruit les watts avec `cp*mass_flux` et la conductivité physique.
+
+Un champ optionnel de conductivité par cellule et des conditions par face
+parcourent l'assemblage, la résolution et les deux reconstructions indépendantes.
+Le noyau scalaire utilise la résistance série géométrique
+`A/(d_owner/lambda_owner + d_neighbour/lambda_neighbour)`.
+La paroi à flux imposé emploie `q_flux/lambda_owner` face par face.
+
+`build_effective_thermal_conductivity` reçoit une viscosité turbulente
+**cinématique** acceptée et un `Pr_t` explicite pour calculer
+`lambda_eff = lambda + rho*cp*nu_t/Pr_t`. Les champs, propriétés, positivité,
+finitude et débordements sont contrôlés. Cette API ne constitue pas le
+raccordement du sélecteur RANS de l'application ; ce dernier reste fermé.
+
+### Conditions d'acceptation et interfaces
+
+Le stage de production impose une température absolue finale finie et
+strictement positive, en plus des portes linéaire et de conservation.
+`temperature_admissibility_passed` et `nonpositive_temperature` rendent ce
+refus explicite dans les sorties. Les noyaux mathématiques génériques ne
+bornent pas artificiellement une solution négative.
+
+Le noyau CHT à deux régions exige une interface non vide, des faces frontières
+uniques, des paramètres finis, des normales extérieures opposées, le même pas
+physique dans les deux régions et un flux massique nul à l'interface.
+Les deux régions sont vérifiées avant mutation. Ces contrôles ne prouvent pas
+une correspondance générale des polygones et ne fournissent pas un rollback
+pour toute erreur apparue pendant les itérations couplées.
+
+### Exécutables, import et reprise
+
+L'exécutable solide refuse les déclarations de rayonnement, RANS, multiphase,
+physique de fluide incompatible et frontières non traitées. Les deux
+exécutables stationnaires refusent une évolution temporelle déclarée, y
+compris une sélection numérique canonique de la famille `Temporal`.
+
+Une propriété matériau absente reste non renseignée (zéro ou chaîne vide)
+dans les schémas Python/C++ et le chargement HDF5 ; elle n'est plus remplacée
+par une propriété d'air. Les validateurs propres à chaque physique décident
+quelles propriétés sont obligatoires. L'import Fluent sérialise la gravité
+comme trois scalaires natifs séparés par des espaces.
+
+Le chemin énergie de `cfdx_production_solver` écrit et recharge T dans le DAT.
+T est obligatoire si une reprise thermique est demandée, contrôlée avant le
+porteur, prioritaire sur l'initialisation scalaire et transmise à la première
+source de flottabilité. Cette reprise stationnaire n'est pas une équivalence
+de trajectoire : le flux de faces autoritaire et les historiques physiques
+complets ne sont pas encore persistés. Une exception d'énergie ne peut plus
+publier le succès du seul porteur comme succès de la simulation complète.
+
+### Preuves discriminantes et travaux restants
+
+Les sources des oracles sont `test_energy_integration_contracts.cpp`,
+`test_cht_integration_contracts.cpp` et les tests Python d'application solide,
+d'import et de reprise thermique. Le rapport d'intégration conserve les
+commandes, journaux, nombres de cas, versions et l'arbre Git exact testé.
+Exemples d'oracles : équilibre `Su=600,Sp=-2,T=300` à bilan nul ; deux couches
+`lambda={1,4}` à T={360,310} K et flux160 W ; transmission d'un état SST
+accepté sur trois cellules à l'énergie avec solution tridiagonale indépendante ;
+champ repris T={310,320,330} K alimentant la première force de Boussinesq ;
+refus d'une solution absolue négative malgré des bilans fermés.
+
+Restent à implémenter ou qualifier : schéma et pilote multi-régions, contact
+par interface et conductivité variable dans CHT, choix `Pr_t` et construction
+applicative de l'état SST, transaction et reprise de tous les champs couplés,
+enthalpie et propriétés variables en T, temps physique, domaine Boussinesq,
+MMS/GCI et références de convection turbulente, CHT et rayonnement.
diff --git a/docs/audits/11-rayonnement.md b/docs/audits/11-rayonnement.md
index 7bae7500..ec571663 100644
--- a/docs/audits/11-rayonnement.md
+++ b/docs/audits/11-rayonnement.md
@@ -5,8 +5,8 @@
 - Statut : APPROUVÉE (2026-10-12, PR #861) — contrat mathématique approuvé ; le code Phase 12 existant doit être audité contre cette fiche avant toute promotion (Règle d'Or 1).
 - Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json + docs/validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json
 - Issues liées : #34 (roadmap, Phase 12) ; docs/validation/THERMAL_RADIATION_REGRESSION_MATRIX.md
-- Dernier commit vérifié : f5a728f205cb45ad6186c9fab92554ae4af06bdd (master, 2026-10-12 — présence des sources et tests confirmée)
-- Limitations connues : bilans radiatifs non exécutés depuis l'environnement d'audit ; statut déclaratif jusqu'à campagne sur HEAD exact.
+- Base de la revue du 2026-10-24 : `038e3aeff40f02af95b8845b931245ce23a1a094`, avec le lot d'intégration local testé ; les résultats de tests concernent l'arbre local exact et ne constituent pas une CI distante. Le relevé historique du 2026-10-12 portait sur `f5a728f205cb45ad6186c9fab92554ae4af06bdd`.
+- Limitations connues : seuls les échanges S2S fermés à facteurs de vue fournis par le cas sont raccordés au chemin de production. Les noyaux P1/DOM/Rosseland et les estimateurs géométriques ne sont pas une qualification industrielle ni un raccordement production implicite ; voir §7.
 
 ## 1. Références sources
 - M. F. Modest, *Radiative Heat Transfer*, 3e éd., Academic Press, 2013 — ISBN : 978-0-12-386944-9. (référence canonique pour tous les modèles ci-dessous)
@@ -25,11 +25,14 @@ La convention entre grandeurs intensives ($I_b=\sigma T^4/\pi$) et intégrales (
 $$\mathbf{q}_r=-k_R\nabla T,\qquad k_R=\frac{16\sigma n^2 T^3}{3\beta_R}$$
 
 **P1** (approximation différentielle) : $I\approx\frac{1}{4\pi}(G+3\,\mathbf{s}\cdot\mathbf{q}_r)$ :
-$$\nabla\cdot\left(\frac{1}{3\kappa}\nabla G\right)=\kappa G-4\kappa\sigma T^4$$
+$$-\nabla\cdot(D\nabla G)+\kappa_aG=4\kappa_a\sigma T^4,\qquad D=\frac{1}{3\beta},\quad \beta=\kappa_a+\kappa_s>0.$$
+Dans le milieu non diffusant de la formulation de base, $\kappa_s=0$ ; le coefficient ci-dessus couvre aussi la diffusion isotrope. Le noyau renvoie la divergence du flux $q_{rad}=\kappa_a(4\sigma T^4-G)$, positive lorsque la matière perd de l'énergie. Un helper utilisant l'intensité moyenne $J_{moy}=G/(4\pi)$ doit donc calculer $q_{rad}=4\pi\kappa_a(I_b-J_{moy})$ ; $J_{moy}$ n'est pas la radiosité surfacique $J_i$ du modèle S2S.
 **Leçon B13** : le facteur $\pi$ (irradiation $G=4\sigma T^4$ du corps noir intégral vs $I_b=\sigma T^4/\pi$) est le piège classique de la mise en œuvre P1 — la convention intégrale vs intensive doit être cohérente de la source à la frontière.
 
 **S2S** : échanges radiatifs entre surfaces grises diffuses par facteurs de forme $F_{ij}$ :
-$$q_i=\frac{\varepsilon_i A_i}{1-\varepsilon_i}\left(E_{bi}-J_i\right)+A_i\sum_j (J_i-J_j)F_{ij},\qquad \sum_j F_{ij}=1$$
+$$G_i=\sum_jF_{ij}J_j+G_{ext,i},\qquad J_i=\varepsilon_i E_{bi}+(1-\varepsilon_i)G_i,\qquad E_{bi}=\sigma T_{wi}^4,$$
+$$q_i''=J_i-G_i=\varepsilon_i(E_{bi}-G_i),\qquad Q_i=A_iq_i''.$$
+$J_i$, $G_i$ et $q_i''$ sont en W/m² ; $Q_i$ est en W. Le flux est positif lorsque la surface perd de l'énergie ; le flux entrant dans l'équation d'énergie vaut $-q_i''$. Ces expressions restent régulières pour $\varepsilon_i=1$. Pour une enceinte fermée sans irradiation extérieure, $\sum_jF_{ij}=1$ et $A_iF_{ij}=A_jF_{ji}$ impliquent $\sum_iQ_i=0$. Pour un système ouvert, $\sum_jF_{ij}\leq1$ ; le noyau peut affecter la fraction manquante à l'ambiance. Une irradiation extérieure explicite remplace ce terme ambiant et ne doit pas être comptée une seconde fois. Le contrat courant exige $F_{ii}=0$ ; cette restriction interdit de représenter en un seul patch une surface concave qui se voit elle-même.
 **DOM (Fiveland)** : RTE discrétisée en $N$ ordonnées $\mathbf{s}_k$ avec quadrature $w_k$ :
 $$\mu_k\frac{\partial I_k}{\partial x}+\eta_k\frac{\partial I_k}{\partial y}+\xi_k\frac{\partial I_k}{\partial z}=\kappa(I_b-I_k)$$
 
@@ -38,49 +41,56 @@ $$S_r=\kappa\left(G-4\sigma T^4\right)$$
 
 ## 3. Discrétisation Volume Fini
 - **Rosseland** : le terme $-k_R\nabla T$ entre la matrice de conduction/convection comme conductivité effective totale $k_{eff}=k+k_R$ — **leçon B18** : la conductivité **moléculaire k ne doit jamais être perdue** (aucun remplacement, uniquement l'addition) ; $k_R$ est assemblé dans la matrice, itéré de Picard externe sur T.
-- **P1** : équation de diffusion elliptique pour $G$ discrétisée comme le Laplacien (fiche 01) avec $\Gamma=1/(3\kappa)$ par face ; maillages non orthogonaux : correction de non-orthogonalité (fiche 01 §3).
-- **S2S** : facteurs de forme par ray tracing ou contours (leçon **B5** : réciprocité $A_iF_{ij}=A_jF_{ji}$ et fermeture $\sum_jF_{ij}=1$ vérifiées à la tolérance machine par test dédié — un facteur de forme erroné viole la conservation d'énergie radiative).
+- **P1** : équation de diffusion elliptique pour $G$ discrétisée comme le Laplacien (fiche 01) avec $\Gamma=D$ par face ; le contrôle d'ordre et de non-orthogonalité reste requis pour chaque chemin qui le revendique.
+- **S2S** : résoudre le système de radiosité construit avec les facteurs de forme. Leçon **B5** : réciprocité $A_iF_{ij}=A_jF_{ji}$ et fermeture d'une enceinte $\sum_jF_{ij}=1$ sont nécessaires à la conservation, mais ne prouvent pas l'exactitude géométrique des facteurs. L'estimateur par centroïdes est une approximation différentielle sans occlusion ; une normalisation conservatrice ne lui ajoute pas la visibilité manquante.
 - **DOM** : balayage par direction, schéma UPWIND spatial borné ; symétrie des quadratures conservée (aucune direction lâchée) ; verrouillage (ray effect) documenté par nombre d'ordres.
 - **Source couplé** $S_r=\kappa(G-4\sigma T^4)$ — **leçon B3** : linéarisation de Newton pour l'implicite :
 $$\frac{\partial S_r}{\partial T}=\kappa\left(\frac{\partial G}{\partial T}-16\sigma T^3\right)$$
-le terme $-16\kappa\sigma T^3$ entre au coefficient diagonal $a_P$ (stabilise), $\partial G/\partial T$ selon le modèle (0 pour S2S, implicite pour P1/Rosseland).
+Pour un couplage ségrégé qui fige $G$ à l'itération extérieure, $S_p=-16\kappa\sigma T_0^3\leq0$ et $S_u=S_r(T_0)-S_pT_0$ ; l'assemblage ajoute $-S_pV$ à la diagonale. Il s'agit d'une linéarisation partielle de la source d'émission, pas du Jacobien complet incluant $\partial G/\partial T$. Le bilan indépendant doit reconstruire $S_u+S_pT$ à la température effectivement résolue. S2S apporte un flux de frontière ; cette source volumique P1 ne décrit pas son couplage.
 
 ## 4. Traitement des frontières
 - Surfaces grises diffuses : condition Robin radiative sur $G$ (P1, Marshak) ou contribution directe au système S2S ; émissivités $\varepsilon_i$ déclarées par patch.
-- Frontières noires ($\varepsilon=1$) : $J=E_b$ ; réfléchissantes ($\varepsilon=0$) : symétrie radiative.
-- Ouvertures (vers l'extérieur à $T_\infty$) : flux net $\sigma(T_w^4-T_\infty^4)$ linéarisé de la même façon (B3).
-- **Gauge** : le rayonnement S2S n'a pas de null-space (bilan fermé), P1 a le null-space du Laplacien si toutes les frontières sont adiabatiques radiativement — point de référence à la P1 (fiche 04) si le cas se présente.
+- Frontières noires ($\varepsilon=1$) : $J=E_b$ ; surface opaque parfaitement réfléchissante ($\varepsilon=0$) : $J=G$ et flux net nul. Le stage S2S de production exige $\varepsilon>0$ ; une réflexion diffuse n'est pas une condition de symétrie géométrique générale.
+- Ouvertures : le flux $\sigma(T_w^4-T_\infty^4)$ correspond à une surface noire voyant uniquement une ambiance noire ; une surface grise isolée voit le facteur $\varepsilon$, et une ouverture partielle exige la fraction de vue correspondante.
+- **Null-space** : une absorption $\kappa_a>0$ ajoute un terme de réaction à P1 et supprime le mode constant, même avec des frontières Neumann. Pour $\kappa_a=0$, extinction positive et Neumann homogène, le Laplacien retrouve ce mode et exige une compatibilité et une référence. Une enceinte S2S entièrement réfléchissante peut également être singulière ; le contrat courant $\varepsilon_i>0$ exclut ce cas. La fermeture seule ne garantit pas l'inversibilité.
+- **Température S2S** : pour une face de flux imposé, $T_{w,f}=T_P+g_fd_{n,f}$, avec $g_f=(q''_{imposé}+q''_{rad,appliqué})/\lambda_{eff,P}$ et $d_{n,f}=|({\bf x}_f-{\bf x}_P)\cdot{\bf S}_f|/|{\bf S}_f|$. Utiliser les mêmes gradients par face et le même flux radiatif appliqué que le solve énergie. Une distance euclidienne ou la réponse radiative future change le problème résolu.
 
 ## 5. Pièges connus
 - **B13** : facteur pi du P1 (convention $G$ intégrale vs $I_b$ intensive) — cohérence de bout en bout exigée.
-- **B3** : signe du source radiatif dans l'équation d'énergie ($S_r$ chauffe < 0 ou refroidit > 0 selon $G$ vs $4\sigma T^4$) + linéarisation de Newton obligatoire pour la stabilité du couplage implicite.
+- **B3** : $S_r=\kappa(G-4\sigma T^4)$ chauffe si $S_r>0$ et refroidit si $S_r<0$ ; le noyau $q_{rad}=-S_r$ emploie le signe opposé. Pour S2S, injecter $-q_i''$ au domaine. Une convention cohérente doit couvrir noyau, assemblage et bilan reconstruit.
 - **B18** : Rosseland — conserver $k$ moléculaire en plus de $k_R$ ; le $c_p$ du bilan d'énergie reste celui de la fiche 10 (aucune confusion $k_R$/$c_p$).
 - **B5** : S2S — réciprocité et fermeture des facteurs de forme à la tolérance machine ; petit nombre de surfaces : conditionnement du système (facteurs de forme) à surveiller.
 - Domaine de validité : Rosseland uniquement $\tau\gg 1$ ; P1 pour géométries sans obstruction forte ; DOM/S2S sinon — le modèle utilisé doit être déclaré avec sa justification dans le rapport de cas (pas de Rosseland silencieux sur cas optiquement mince).
-- Couplage instable si le source radiatif est purement explicite (oscillation T) — toujours linéariser (B3).
+- Un couplage explicite fort peut osciller. Le couplage volumique peut stabiliser l'émission par $S_p\leq0$ ; le stage S2S utilise une itération de Picard sous-relaxée, avec résidu de point fixe évalué sur la réponse non relaxée. Ni une sous-relaxation ni un solve linéaire convergé ne remplacent ce critère de couplage.
+- Rejeter NaN/Inf et les valeurs non physiques avant les calculs. Les helpers P1 rejettent également les résultats non représentables ; les centroïdes, normales et aires S2S sont validés même pour un appel direct du noyau par paire ou une seule surface. Les températures reconstruites doivent être positives et finies face par face avant moyenne.
 - Anisotropie : hors scope ici (milieu non diffusant) ; toute extension devra avoir sa fiche.
 
 ## 6. Ordre de précision attendu et validation exigée
-- Ordre : P1/Rosseland = ordre du Laplacien ($O(h^2)$ orthogonales) ; DOM : $O(h)$ en espace par direction, erreur angulaire dominée par la quadrature (à documenter) ; S2S : exact par paires, l'erreur vient des facteurs de forme (tolérance machine exigée sur réciprocité/fermeture).
+- Ordre attendu à vérifier : P1/Rosseland = ordre du Laplacien ($O(h^2)$ sur maillages orthogonaux réguliers) ; DOM : $O(h)$ en espace pour le schéma amont, avec erreur angulaire à mesurer séparément ; S2S : le système de radiosité est résolu algébriquement pour des surfaces isothermes et des facteurs donnés. Les erreurs de facteurs, de visibilité et d'agrégation des températures subsistent même lorsque tous les bilans ferment.
 - MMS pour le P1 (éq. de diffusion sur $G$) sur au moins 3 maillages ; GCI sur les QoI ($q_r$ pariétal, $T_{max}$, fraction radiative du flux total).
 - Oracles quantitatifs : **cavité chauffée** (benchmark classique : parois à températures imposées, solution S2S exacte en géométrie simple), bilan d'énergie radiatif fermé à la tolérance machine (ce qui entre = ce qui sort + ce qui est absorbé), comparaison croisée entre modèles P1/S2S/DOM sur le même cas (doivent converger vers la même physique dans leur domaine de validité) ; matrice docs/validation/THERMAL_RADIATION_REGRESSION_MATRIX.md.
 - Conservation : bilan radiatif fermé machine ; bornitude de T avec télémétrie.
 - Parallèle : parité CPU/MPI à 1e-12 sur le même maillage (fiche 12).
 
 ## 7. Statut d'implémentation CFDX
-- Fichiers source concernés : src/cfdx/physics/radiation.h, radiation_models.h, radiation_advanced.h, radiation_s2s.h, radiation_solver.h, radiation_streaming.h.
-- Tests exécutables : tests/validation/test_phase12_radiation_vv.cpp, test_s2s_radiation_vv.cpp, test_radiation_vv.cpp, test_thermal_radiation_regression.cpp, test_thermal_radiation_model_matrix.cpp.
-- Écarts connus code <-> fiche : à auditer ligne à ligne lors de la Phase 3 — B3, B5, B13 et B18 sont les quatre premiers points de contrôle.
-- Gaps et prochaines actions : bilan radiatif fermé comme oracle machine systématique ; cavité chauffée benchmark avec superposition graphique dans le rapport ; comparaison croisée P1/S2S/DOM documentée ; promotion BLOCKED derrière Phase 1.
+- Fichiers source concernés : `radiation.h`, `radiation_models.h`, `radiation_advanced.h`, `radiation_s2s.h`, `radiation_solver.h`, `radiation_streaming.h`, `m1_m4_models.h`, `production_radiation_contract.h` et `production_radiation_s2s_stage.h`, sous `src/cfdx/physics/`.
+- Tests exécutables : `tests/validation/test_phase12_radiation_vv.cpp`, `test_s2s_radiation_vv.cpp`, `test_radiation_vv.cpp`, `test_thermal_radiation_regression.cpp`, `test_thermal_radiation_model_matrix.cpp`, `test_production_s2s_energy_stage_vv.cpp` et `test_production_s2s_energy_outer_loop_vv.cpp`.
+- Noyaux disponibles : P1, DOM, Rosseland et S2S. Le contrat `validate_production_radiation_selection` ne raccorde que S2S au chemin de production décrit ici ; il continue de rejeter explicitement les autres modèles. L'existence d'un test de noyau n'autorise pas à les annoncer comme modèles de production intégrés.
+- Contrôles de robustesse du 2026-10-24 : helpers P1 finis et facteur $4\pi$ ; aires, centres, normales et sommets triangulaires finis avant leurs calculs ; origines et directions finies pour le noyau d'intersection direct ; stabilité de la normalisation des normales ; rejet des résultats non représentables ; maintien de la réciprocité et de la diagonale nulle. Des tests négatifs reproduisent les acceptations silencieuses antérieures de NaN/Inf. Ces gardes ne résolvent ni l'occlusion générale ni la précision géométrique.
+- Températures S2S du chemin de production : moyenne surfacique des températures de face reconstruites avec les BC du solve énergie ; les coefficients locaux optionnels $\lambda_{eff,P}$ passent par `build_production_temperature_face_conditions`. Le résultat conserve séparément `applied_heat_flux` (flux appliqué au solve énergie et à la reconstruction) et `injected_heat_flux` (nouvelle réponse radiosité). Le stage deux-passes et la boucle externe propagent le même champ de conductivité. Un test avec décalage tangentiel impose l'usage de $d_n$, un test variable contrôle $q''/\lambda_{eff,P}$, et le replay indépendant contrôle la boucle couplée.
+- Admissibilité thermique : un second passage énergie rejeté, notamment pour $T\leq0$ malgré des bilans fermés, est retourné avec son verdict structuré et son champ non écrêté. Le stage deux-passes ne calcule pas ensuite le diagnostic de retard pariétal sur ce champ rejeté ; une exception de reconstruction ne doit pas masquer la cause thermique. Le diagnostic de retard n'est donc exploitable qu'après un passage énergie accepté.
+- Restriction d'agrégation : une température moyenne et une radiosité uniformes représentent chaque patch. En général $\overline{T_w}^{\,4}\neq\overline{T_w^4}$. Des parois fortement non isothermes nécessitent des surfaces plus fines ou un modèle par face accompagné de facteurs de vue cohérents ; la correction de reconstruction ne supprime pas cette erreur.
+- Critères restant à démontrer avant extension de qualification : précision des facteurs de vue sur géométries obstruées avec oracles analytiques ou références indépendantes ; études de maillage et d'agrégation ; convergence angulaire DOM ; MMS P1 à trois niveaux et GCI ; références physiques du domaine visé et parité MPI effective. Les bilans et tests présents constituent des vérifications ciblées, pas la clôture de ces critères.
 - Étape 2 du plan #885 (stage production S2S) : fichiers `src/cfdx/physics/production_radiation_s2s_stage.h` et test `tests/validation/test_production_s2s_energy_stage_vv.cpp` — un passage de Picard retardé sur écoulement figé : surfaces participantes = patches WALL déclarant `emissivity` (parsing strict, (0, 1]), facteurs de vue fournis par le cas avec fermeture |Σ_j F_ij − 1| ≤ 1e-10 et réciprocité B5 validées machine (rejets fail-closed sinon), gates bloquantes identité par surface et bilan radiatif fermé ≤ 1e-10, flux nets ré-injectés dans l'étage énergie via le canal `radiative_heat_flux_into_domain` (interface consolidée #999/#1004, nommée `extra_wall_heat_flux` lors de l'étape 2 ; signe −q_i, positif entrant, même convention que `heat_flux`).
 
-- Étape 3 du plan #885 (boucle externe énergie–rayonnement) : `run_production_s2s_energy_outer_loop` dans `src/cfdx/physics/production_radiation_s2s_stage.h` (test `tests/validation/test_production_s2s_energy_outer_loop_vv.cpp`) itère les deux primitives de l'étape 2 INCHANGÉES (passage énergie + solve radiosité, gates internes fiches 14/11 inchangées) : sous-relaxation fixe ω = 0,25 du flux radiatif appliqué (paramètre de chemin uniquement), budget 200 itérations fail-closed, critères bloquants SÉPARÉS évalués sur la réponse NON relaxée — résidu de point fixe radiatif max_i |R_i(T*) − F_i| / scale_i ≤ 1e-10 (scale = même échelle brute que la gate d'identité par surface) ET stationnarité des températures de paroi ≤ 1e-10 ; le verdict est le point fixe couplé exact (INCONCLUSIVE n'est jamais PASS) ; historique par itération dans `radiation_s2s.json`.
+- Étape 3 du plan #885 (boucle externe énergie–rayonnement) : `run_production_s2s_energy_outer_loop` dans `src/cfdx/physics/production_radiation_s2s_stage.h` (test `tests/validation/test_production_s2s_energy_outer_loop_vv.cpp`) itère le passage énergie et le solve radiosité, avec la reconstruction de frontière décrite ci-dessus : sous-relaxation fixe ω = 0,25 du flux radiatif appliqué, budget 200 itérations avec échec explicite si épuisé, critères bloquants séparés évalués sur la réponse non relaxée — résidu de point fixe radiatif max_i |R_i(T_w(T*,F_i)) − F_i| / scale_i ≤ 1e-10 (échelle brute de la gate d'identité par surface) et stationnarité des températures de paroi ≤ 1e-10. Le verdict est la convergence du point fixe discret à ces tolérances, et ne constitue ni une solution mathématiquement exacte ni une qualification physique ; historique par itération dans `radiation_s2s.json`.
 
 - Étape 4 du plan #885 (flottabilité) : les gates du contrat S2S (`production_radiation_contract.h`) et du stage acceptent désormais `energy_model = "buoyant"` en plus de `"conjugate"` (messages de rejet mis à jour, rejet fail-closed inchangé pour toute autre valeur) ; primitives, critères, tolérances et budget des étapes 2 et 3 INCHANGÉS ; la boucle externe énergie-rayonnement de l'étape 3 reste consommable telle quelle et devient l'opérateur de réponse interne (énergie + rayonnement) de la boucle externe flottabilité de l'étape 4 (fiche 10 §7).
 
 - Correctif de la PR #1016 : les échecs A/C/D de `test_production_buoyancy_outer_loop_vv` sont résolus en amont de l'étape 3 — la boucle S2S-énergie imbriquée (étape 3) reste INCHANGÉE : porteur warm accepté à l'entrée par la boucle flottabilité (fiche 10 §7), gate locale de l'étage énergie calibrée au plancher d'arrondi flottant (fiche 14 §7/§8) ; aucune tolérance/gate/oracle/budget de la fiche modifiée.
 
 ## 8. Journal des modifications de la fiche
+- 2026-10-24 — correction des équations S2S (radiosité, irradiation, flux et puissance), du signe de la source énergie, des conditions de null-space et de la linéarisation partielle de l'émission P1 ; distinction explicite noyaux/production/qualification. Gardes finies P1 et géométrie, reconstruction des températures de face selon la loi Fourier avec flux radiatif réellement appliqué et conductivité locale optionnelle. Oracles S2S adaptés à cette loi de frontière, sans modification des tolérances, de la sous-relaxation ni du budget de la boucle externe. Les entrées historiques ci-dessous décrivent leur état à leur date, avant cette correction de reconstruction.
 - 2026-10-12 — création (Phase 3 du Plan Directeur, PR empilée sur governance/constitution-technique).
 - 2026-10-12 — approbation du contrat (PR d'approbation #861) ; l'audit du code contre cette fiche reste requis avant toute promotion (Règle d'Or 1).
 - 2026-10-20 — étape 2 du plan #885 : le contrat `production_radiation_contract.h` accepte désormais « s2s » uniquement avec `energy_model = "conjugate"` (rejet fail-closed explicite sinon, message pointant « conjugate ») ; nouveau stage `production_radiation_s2s_stage.h` couplant la radiosité S2S à l'étage énergie de production (un passage de Picard retardé, écoulement figé, aucune convergence du point fixe couplé revendiquée — l'itération externe est l'étape 3) ; bilan radiatif fermé machine + réciprocité B5 = oracles bloquants (test `test_production_s2s_energy_stage_vv`) ; aucune tolérance/gate/oracle de la fiche modifiée.
diff --git a/scripts/t2_heated_channel_calibration.py b/scripts/t2_heated_channel_calibration.py
index b7321f68..b0c157d9 100644
--- a/scripts/t2_heated_channel_calibration.py
+++ b/scripts/t2_heated_channel_calibration.py
@@ -21,9 +21,10 @@ Locked case decisions (N9/T2 audit, heated-channel campaign #884):
   * external energy check: mdot*cp*(Tb_out - Tb_in) vs q''*2*L, with mdot
     measured from the u_x probes (not an oracle)
 
-The script never modifies the solver. Evidence (JSON + TXT) is always
-written to --output-dir before exiting; the exit code is 1 only when a
-mesh fails to run/converge (the err_rel measurement itself is not gated).
+The script never modifies the solver. Evidence (JSON + TXT) distinguishes
+execution, Poiseuille and Nusselt acceptance. By default the exit code remains
+diagnostic (failure to execute/converge); --require-qualification additionally
+gates the documented physical targets, without changing their thresholds.
 """
 
 from __future__ import annotations
@@ -608,6 +609,38 @@ def _fit_order(dxs: list[float], errs: list[float]) -> tuple[float, float]:
     return p, c
 
 
+def evaluate_campaign(meshes: list[dict], execution_ok: bool) -> dict:
+    """Apply the frozen campaign targets; a successful process is not a pass.
+
+    The complete prescribed mesh sequence is required, and every numerical
+    gate explicitly rejects non-finite values. These are the T2 campaign's
+    targets, not a qualification of unrelated models or coupled physics.
+    """
+    complete = ([(entry.get("nx"), entry.get("ny")) for entry in meshes] == MESHES)
+    execution_passed = bool(execution_ok and complete and
+                            all(entry.get("status") == "PASS" for entry in meshes))
+
+    def within(value, limit):
+        return (isinstance(value, (int, float)) and math.isfinite(value) and
+                0.0 <= value <= limit)
+
+    poiseuille_passed = execution_passed and all(
+        within(entry.get("l2_poiseuille"), L2_GATE) for entry in meshes)
+    nusselt_passed = (execution_passed and
+                      all(within(entry.get("err_rel"), float("inf")) for entry in meshes) and
+                      within(meshes[-1].get("err_rel"), ERR_TARGET))
+    qualification_passed = bool(poiseuille_passed and nusselt_passed)
+    return {
+        "execution_passed": execution_passed,
+        "poiseuille_passed": bool(poiseuille_passed),
+        "nusselt_passed": bool(nusselt_passed),
+        "qualification_passed": qualification_passed,
+        "verdict": ("physics_qualified" if qualification_passed else
+                    "execution_complete_physics_not_qualified" if execution_passed else
+                    "execution_failure"),
+    }
+
+
 def _write_evidence(out_dir: Path, campaign: dict) -> None:
     run_id = campaign["run_id"]
     json_path = out_dir / f"t2_nu_calibration_{run_id}.json"
@@ -650,6 +683,8 @@ def _write_evidence(out_dir: Path, campaign: dict) -> None:
         )
     lines.append("")
     lines.append(f"verdict: {campaign['verdict']}")
+    for key in ("execution_passed", "poiseuille_passed", "nusselt_passed", "qualification_passed"):
+        lines.append(f"{key}: {str(campaign[key]).lower()}")
     lines.append("interpretation (locked protocol):")
     lines.append(
         "- err_rel <= 0.05 on the finest mesh -> wall-treatment defect dominates"
@@ -672,6 +707,8 @@ def main() -> int:
     parser.add_argument("--solver", required=True, help="production solver binary")
     parser.add_argument("--output-dir", required=True, help="evidence directory")
     parser.add_argument("--run-id", default="local", help="campaign run id")
+    parser.add_argument("--require-qualification", action="store_true",
+                        help="also fail the process when the frozen physical targets are not met")
     args = parser.parse_args()
 
     solver = Path(args.solver)
@@ -732,27 +769,26 @@ def main() -> int:
                 print(f"T2NU_MESH nx={nx} ny={ny} FAILED error={message}")
             campaign["meshes"].append(entry)
 
-    if ok:
+    campaign.update(evaluate_campaign(campaign["meshes"], ok))
+    if campaign["execution_passed"]:
         dxs = [entry["dx"] for entry in campaign["meshes"]]
         errs = [entry["err_rel"] for entry in campaign["meshes"]]
         p, c = _fit_order(dxs, errs)
         campaign["order_fit"] = {"dx": dxs, "err": errs, "p": p, "C": c}
         finest = campaign["meshes"][-1]
-        within = finest["err_rel"] <= ERR_TARGET
-        campaign["verdict"] = (
-            "err_rel_within_target" if within else "err_rel_above_target"
-        )
         print(f"T2NU_ORDER p={p:.3f} C={c:.6f}")
         print(
-            f"T2NU_CAMPAIGN PASS finest_err_rel={finest['err_rel']:.6f} "
+            f"T2NU_DIAGNOSTIC COMPLETE execution_passed=true "
+            f"qualification_passed={str(campaign['qualification_passed']).lower()} "
+            f"finest_err_rel={finest['err_rel']:.6f} "
             f"target={ERR_TARGET} verdict={campaign['verdict']}"
         )
     else:
-        campaign["verdict"] = "execution_failure"
-        print("T2NU_CAMPAIGN FAIL execution_failure")
+        print("T2NU_DIAGNOSTIC FAIL execution_failure")
 
     _write_evidence(out_dir, campaign)
-    return 0 if ok else 1
+    passed = campaign["qualification_passed"] if args.require_qualification else campaign["execution_passed"]
+    return 0 if passed else 1
 
 
 if __name__ == "__main__":
diff --git a/src/cfdx/core/linalg/solver_status.h b/src/cfdx/core/linalg/solver_status.h
index eeb7cff2..fabf688a 100644
--- a/src/cfdx/core/linalg/solver_status.h
+++ b/src/cfdx/core/linalg/solver_status.h
@@ -61,6 +61,14 @@ struct SolverResult {
     double min_true_residual = std::numeric_limits<double>::infinity();
     std::size_t first_residual_replacement = 0;
     std::size_t last_residual_replacement = 0;
+    // Contract reference: docs/audits/11-solveur-energie-et-convergence-lineaire.md.
+    // For solvers that relax the returned field after solving the linear
+    // system, distinguish the linear candidate residual from the residual
+    // of the actual returned field. NaN means not measured/not applicable.
+    double residual_before_relaxation = std::numeric_limits<double>::quiet_NaN();
+    double residual_after_relaxation = std::numeric_limits<double>::quiet_NaN();
+    double residual_relative_before_relaxation = std::numeric_limits<double>::quiet_NaN();
+    double residual_relative_after_relaxation = std::numeric_limits<double>::quiet_NaN();
 };
 
 }  // namespace core
diff --git a/src/cfdx/io/hdf5/case_hdf5_io.cpp b/src/cfdx/io/hdf5/case_hdf5_io.cpp
index 7f8d87cb..8fbea134 100644
--- a/src/cfdx/io/hdf5/case_hdf5_io.cpp
+++ b/src/cfdx/io/hdf5/case_hdf5_io.cpp
@@ -646,12 +646,16 @@ static void case_setup_from_json(const mini_json::value& json, CaseSetup& setup)
             if (!mat_v.is_object()) continue;
             MaterialSpec mat;
             mat.name = json_get_string(mat_v, "name", "air");
-            mat.density = json_get_number(mat_v, "density", 1.225);
-            mat.dynamic_viscosity = json_get_number(mat_v, "dynamic_viscosity", 1.789e-5);
-            mat.thermal_conductivity = json_get_number(mat_v, "thermal_conductivity", 0.024);
-            mat.specific_heat = json_get_number(mat_v, "specific_heat", 1006.43);
-            mat.molecular_weight = json_get_number(mat_v, "molecular_weight", 28.97);
-            mat.eos_model = json_get_string(mat_v, "eos_model", "ideal_gas");
+            // Preserve the typed schema's unspecified values. Supplying air
+            // properties here lets incomplete material declarations pass the
+            // production validators and changes the physical problem silently.
+            // Each selected physics validates the properties it actually needs.
+            mat.density = json_get_number(mat_v, "density", 0.0);
+            mat.dynamic_viscosity = json_get_number(mat_v, "dynamic_viscosity", 0.0);
+            mat.thermal_conductivity = json_get_number(mat_v, "thermal_conductivity", 0.0);
+            mat.specific_heat = json_get_number(mat_v, "specific_heat", 0.0);
+            mat.molecular_weight = json_get_number(mat_v, "molecular_weight", 0.0);
+            mat.eos_model = json_get_string(mat_v, "eos_model");
             mat.turbulence_model = json_get_string(mat_v, "turbulence_model", "laminar");
 
             const mini_json::value* ep = mini_json::find(mat_v, "extra_properties");
diff --git a/src/cfdx/physics/cht_solver.h b/src/cfdx/physics/cht_solver.h
index c36cc659..e66bd902 100644
--- a/src/cfdx/physics/cht_solver.h
+++ b/src/cfdx/physics/cht_solver.h
@@ -23,6 +23,8 @@ struct ChtInterfaceControls {
     double matching_tolerance = 1e-8;
     double temperature_tolerance = 1e-8;
     double area_relative_tolerance = 1e-6;
+    // Bound |1+n1.n2| for unit outward normals at a conformal interface.
+    double normal_opposition_tolerance = 1e-8;
 };
 
 struct ChtInterfaceFacePair {
@@ -38,13 +40,32 @@ inline std::vector<ChtInterfaceFacePair> match_cht_interface(
     const ChtInterfaceControls& c)
 {
     if(c.region1_patch.empty() || c.region2_patch.empty() ||
+       !std::isfinite(c.conductivity1) || !std::isfinite(c.conductivity2) ||
+       !std::isfinite(c.tolerance) || !std::isfinite(c.relaxation) ||
+       !std::isfinite(c.matching_tolerance) ||
+       !std::isfinite(c.temperature_tolerance) ||
+       !std::isfinite(c.area_relative_tolerance) ||
+       !std::isfinite(c.normal_opposition_tolerance) ||
        c.conductivity1<=0.0 || c.conductivity2<=0.0 ||
        c.tolerance<=0.0 || c.max_iterations==0 ||
        !(c.relaxation>0.0 && c.relaxation<=1.0) ||
        c.matching_tolerance<=0.0 || c.temperature_tolerance<=0.0 ||
-       c.area_relative_tolerance<=0.0)
+       c.area_relative_tolerance<=0.0 ||
+       !(c.normal_opposition_tolerance>0.0 && c.normal_opposition_tolerance<1.0))
         throw std::invalid_argument("invalid CHT interface controls");
 
+    const auto validate_geometry = [](const cfdx::core::Mesh& mesh,
+                                      const FvGeometry& g) {
+        if (g.face_centres.size()!=mesh.n_faces() ||
+            g.face_area_vectors.size()!=mesh.n_faces() ||
+            g.face_patch.size()!=mesh.n_faces() ||
+            g.cell_centres.size()!=mesh.n_cells() ||
+            g.cell_volumes.size()!=mesh.n_cells())
+            throw std::invalid_argument("CHT geometry does not match the region mesh");
+    };
+    validate_geometry(mesh1,g1);
+    validate_geometry(mesh2,g2);
+
     std::size_t p1=mesh1.boundary().n_patches();
     std::size_t p2=mesh2.boundary().n_patches();
     for(std::size_t p=0;p<p1;++p)
@@ -56,9 +77,41 @@ inline std::vector<ChtInterfaceFacePair> match_cht_interface(
 
     const auto& f1=mesh1.boundary().patch(p1).face_ids;
     const auto& f2=mesh2.boundary().patch(p2).face_ids;
+    if(f1.empty() || f2.empty())
+        throw std::invalid_argument("CHT interface patches must be non-empty");
     if(f1.size()!=f2.size())
         throw std::invalid_argument("CHT interface face count mismatch");
 
+    const auto validate_faces = [](const cfdx::core::Mesh& mesh,
+                                   const FvGeometry& g, std::size_t patch) {
+        if (mesh.boundary().patch(patch).type==cfdx::core::PatchType::EMPTY)
+            throw std::invalid_argument("CHT interface cannot be an EMPTY patch");
+        std::vector<bool> seen(mesh.n_faces(),false);
+        for (const auto f : mesh.boundary().patch(patch).face_ids) {
+            if (f>=mesh.n_faces() || mesh.ownership().neighbour(f)>=0 ||
+                mesh.ownership().owner(f)>=mesh.n_cells() ||
+                g.face_patch[f]!=patch)
+                throw std::invalid_argument("CHT interface requires valid boundary faces");
+            if (seen[f])
+                throw std::invalid_argument("CHT interface has duplicate boundary faces");
+            seen[f]=true;
+            const auto& sf=g.face_area_vectors[f];
+            const auto& cf=g.face_centres[f];
+            const auto& cc=g.cell_centres[mesh.ownership().owner(f)];
+            const double area=sf.mag();
+            if (!(area>0.0) || !std::isfinite(area) ||
+                !std::isfinite(cf.x) || !std::isfinite(cf.y) || !std::isfinite(cf.z) ||
+                !std::isfinite(cc.x) || !std::isfinite(cc.y) || !std::isfinite(cc.z) ||
+                !(sf.dot(cf-cc)>0.0))
+                throw std::invalid_argument("CHT interface geometry must be finite and outward-oriented");
+            const double distance=boundary_normal_distance(g,f,mesh.ownership().owner(f));
+            if (!(distance>0.0) || !std::isfinite(distance))
+                throw std::invalid_argument("CHT interface normal distance must be finite and positive");
+        }
+    };
+    validate_faces(mesh1,g1,p1);
+    validate_faces(mesh2,g2,p2);
+
     std::vector<bool> used(f2.size(),false);
     std::vector<ChtInterfaceFacePair> pairs;
     pairs.reserve(f1.size());
@@ -78,6 +131,10 @@ inline std::vector<ChtInterfaceFacePair> match_cht_interface(
         const double area_scale=std::max({a1,a2,1e-30});
         if(std::abs(a1-a2)/area_scale>c.area_relative_tolerance)
             throw std::runtime_error("CHT interface face area mismatch");
+        const auto n1=g1.face_area_vectors[face1]*(1.0/a1);
+        const auto n2=g2.face_area_vectors[f2[bestj]]*(1.0/a2);
+        if (std::abs(1.0+n1.dot(n2))>c.normal_opposition_tolerance)
+            throw std::invalid_argument("CHT interface requires opposite outward normals");
         used[bestj]=true;
         pairs.push_back({
             face1,f2[bestj],
@@ -111,8 +168,61 @@ inline ChtSolveResult solve_two_region_cht(
     const ScalarBoundaryConditions& bcs1 = {},
     const ScalarBoundaryConditions& bcs2 = {})
 {
+    // Validate both regions before either solve can modify a temperature.
+    // Preserve the early checks from PR #1095 when consolidating with #1097.
+    validate_energy_controls(energy1);
+    validate_energy_controls(energy2);
+    if (energy1.dt!=energy2.dt)
+        throw std::invalid_argument("CHT regions must use the same physical time step");
+    const auto validate_region = [](const cfdx::core::Mesh& mesh,
+                                    const FvGeometry& geometry,
+                                    const auto& flux, const auto& temperature,
+                                    const auto& source) {
+        if (flux.size()!=mesh.n_faces() || flux.dimension()!=1 ||
+            temperature.size()!=mesh.n_cells() || temperature.dimension()!=1 ||
+            source.size()!=mesh.n_cells() || source.dimension()!=1 ||
+            geometry.cell_volumes.size()!=mesh.n_cells())
+            throw std::invalid_argument("CHT region field dimensions do not match the mesh");
+        for (std::size_t i=0;i<mesh.n_cells();++i) {
+            if (!std::isfinite(source(i)))
+                throw std::invalid_argument("CHT region source must be finite");
+            if (!std::isfinite(temperature(i)))
+                throw std::invalid_argument("CHT region temperature must be finite");
+            if (!(geometry.cell_volumes[i]>0.0) || !std::isfinite(geometry.cell_volumes[i]))
+                throw std::invalid_argument("CHT region volumes must be finite and positive");
+        }
+        for (std::size_t f=0;f<mesh.n_faces();++f)
+            if (!std::isfinite(flux(f)))
+                throw std::invalid_argument("CHT region mass flux must be finite");
+    };
+    validate_region(mesh1,g1,flux1,T1,source1);
+    validate_region(mesh2,g2,flux2,T2,source2);
+    // The interface update and its flux audit must use the same conductivities
+    // as the two regional energy operators. Otherwise the interface can report
+    // zero imbalance while the assembled energy equations carry different heat
+    // rates. Reject contradictory contracts before mutating either temperature.
+    const auto conductivity_matches = [](double interface_k, double energy_k) {
+        if (!std::isfinite(interface_k) || !std::isfinite(energy_k) ||
+            !(interface_k > 0.0) || !(energy_k > 0.0))
+            return false;
+        const double scale = std::max(std::abs(interface_k), std::abs(energy_k));
+        return std::abs(interface_k - energy_k) <= 1e-12 * scale;
+    };
+    if (!conductivity_matches(controls.conductivity1, energy1.conductivity) ||
+        !conductivity_matches(controls.conductivity2, energy2.conductivity))
+        throw std::invalid_argument(
+            "CHT interface conductivities must match the regional energy conductivities");
+
     const auto pairs=match_cht_interface(mesh1,g1,mesh2,g2,controls);
+    // This kernel couples conductive heat through an impermeable interface.
+    // Advected enthalpy crossing a fluid/solid interface is not represented
+    // by its two-resistance law and must not be silently omitted from it.
+    for (const auto& p : pairs)
+        if (flux1(p.face1)!=0.0 || flux2(p.face2)!=0.0)
+            throw std::invalid_argument("CHT conductive interface must be impermeable (zero mass flux)");
     ChtSolveResult result;
+    const auto T1_previous_time=T1;
+    const auto T2_previous_time=T2;
     std::vector<double> previous_interface_temperature(pairs.size(), std::numeric_limits<double>::quiet_NaN());
 
     for(std::size_t iter=1;iter<=controls.max_iterations;++iter) {
@@ -144,8 +254,8 @@ inline ChtSolveResult solve_two_region_cht(
             values2[p.face2]=Tint;
         }
 
-        auto r1=solve_energy(mesh1,g1,flux1,T1,source1,energy1,bcs1,&fv1);
-        auto r2=solve_energy(mesh2,g2,flux2,T2,source2,energy2,bcs2,&fv2);
+        auto r1=solve_energy(mesh1,g1,flux1,T1,source1,energy1,bcs1,&fv1,nullptr,&T1_previous_time);
+        auto r2=solve_energy(mesh2,g2,flux2,T2,source2,energy2,bcs2,&fv2,nullptr,&T2_previous_time);
         if(!r1.converged || !r2.converged)
             throw std::runtime_error("CHT region energy solve did not converge");
 
diff --git a/src/cfdx/physics/energy_solver.h b/src/cfdx/physics/energy_solver.h
index 43c4e524..36a3bbdd 100644
--- a/src/cfdx/physics/energy_solver.h
+++ b/src/cfdx/physics/energy_solver.h
@@ -18,7 +18,7 @@ namespace cfdx::physics {
 struct EnergySolverControls {
     double density = 1.0;
     double cp = 1000.0;
-    double conductivity = 1.0;
+    double conductivity = 1.0;       // fallback if no cell conductivity is supplied
     double dt = 0.0;                 // <=0 => steady
     double relaxation = 0.9;
     std::size_t max_iterations = 500;
@@ -89,7 +89,7 @@ inline void validate_energy_implicit_source(
     std::size_t n_cells)
 {
     if (!source_implicit) return;
-    if (source_implicit->size()!=n_cells)
+    if (source_implicit->size()!=n_cells || source_implicit->dimension()!=1)
         throw std::invalid_argument("energy implicit source field size mismatch");
     for(std::size_t i=0;i<n_cells;++i)
         if (!std::isfinite((*source_implicit)(i)) || (*source_implicit)(i)>0.0)
@@ -97,6 +97,20 @@ inline void validate_energy_implicit_source(
                 "energy implicit source must be finite and non-positive");
 }
 
+// Conductivity is in W/(m.K), including any frozen turbulent contribution.
+// Density and cp remain constant in this temperature-based energy equation.
+inline void validate_energy_conductivity(
+    const std::vector<double>* conductivity_values, std::size_t n_cells)
+{
+    if (!conductivity_values) return;
+    if (conductivity_values->size()!=n_cells)
+        throw std::invalid_argument("energy conductivity field size mismatch");
+    for (const double k : *conductivity_values)
+        if (!(k>=0.0) || !std::isfinite(k))
+            throw std::invalid_argument(
+                "energy conductivity must be finite and non-negative");
+}
+
 inline ScalarEquation assemble_energy_equation(
     const cfdx::core::Mesh& mesh,
     const FvGeometry& geometry,
@@ -106,12 +120,16 @@ inline ScalarEquation assemble_energy_equation(
     const EnergySolverControls& c,
     const ScalarBoundaryConditions& bcs = {},
     const ScalarBoundaryFaceValues* face_values = nullptr,
-    const cfdx::core::Field<double,cfdx::core::Location::CELL>* source_implicit = nullptr)
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>* source_implicit = nullptr,
+    const std::vector<double>* conductivity_values = nullptr,
+    const ScalarBoundaryFaceConditions* face_conditions = nullptr)
 {
     validate_energy_controls(c);
-    if (source.size()!=mesh.n_cells() || old_temperature.size()!=mesh.n_cells())
+    if (source.size()!=mesh.n_cells() || source.dimension()!=1 ||
+        old_temperature.size()!=mesh.n_cells() || old_temperature.dimension()!=1)
         throw std::invalid_argument("energy field size mismatch");
     validate_energy_implicit_source(source_implicit,mesh.n_cells());
+    validate_energy_conductivity(conductivity_values,mesh.n_cells());
 
     cfdx::core::Field<double,cfdx::core::Location::CELL> su(
         mesh.n_cells(),"energy_source","W/m3",1);
@@ -132,7 +150,38 @@ inline ScalarEquation assemble_energy_equation(
     const auto enthalpy_flux=enthalpy_face_flux(mass_flux,c.cp);
     return assemble_scalar_equation(
         mesh,geometry,enthalpy_flux,c.conductivity,su,sp,bcs,true,
-        face_values,&transient_diag,&transient_rhs);
+        face_values,&transient_diag,&transient_rhs,conductivity_values,
+        ConvectionScheme::UPWIND,nullptr,face_conditions);
+}
+
+// Rebuild physical watts directly from the final state, independently of the
+// matrix: m_dot*cp*T_upwind - lambda_eff*grad(T). The same cell conductivities
+// and resolved face conditions must be supplied to assembly and all audits.
+inline cfdx::core::Field<double,cfdx::core::Location::FACE>
+reconstruct_energy_face_flux(
+    const cfdx::core::Mesh& mesh,
+    const FvGeometry& geometry,
+    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>& temperature,
+    const EnergySolverControls& controls,
+    const ScalarBoundaryConditions& bcs = {},
+    const ScalarBoundaryFaceValues* face_values = nullptr,
+    const std::vector<double>* conductivity_values = nullptr,
+    const ScalarBoundaryFaceConditions* face_conditions = nullptr)
+{
+    validate_energy_controls(controls);
+    validate_energy_conductivity(conductivity_values,mesh.n_cells());
+    if (mass_flux.size()!=mesh.n_faces() || mass_flux.dimension()!=1 ||
+        temperature.size()!=mesh.n_cells() || temperature.dimension()!=1)
+        throw std::invalid_argument("energy flux reconstruction field size mismatch");
+    if (face_conditions &&
+        (face_conditions->conditions.size()!=mesh.n_faces() ||
+         face_conditions->specified.size()!=mesh.n_faces()))
+        throw std::invalid_argument("energy boundary face conditions size mismatch");
+    const auto enthalpy_flux=enthalpy_face_flux(mass_flux,controls.cp);
+    return reconstruct_scalar_transport_flux(
+        mesh,geometry,enthalpy_flux,temperature,controls.conductivity,bcs,true,
+        ConvectionScheme::UPWIND,face_conditions,conductivity_values,face_values);
 }
 
 inline double energy_balance_relative_from_equation(
@@ -170,72 +219,32 @@ inline double energy_balance_relative(
     const EnergySolverControls& controls,
     const ScalarBoundaryConditions& bcs,
     const ScalarBoundaryFaceValues* face_values = nullptr,
-    const std::vector<double>* conductivity_values = nullptr)
+    const std::vector<double>* conductivity_values = nullptr,
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>* source_implicit = nullptr,
+    const ScalarBoundaryFaceConditions* face_conditions = nullptr)
 {
-    if (conductivity_values && conductivity_values->size()!=mesh.n_cells())
-        throw std::invalid_argument("energy conductivity field size mismatch");
-    double net_flux = 0.0;
-    double source_total = 0.0;
-    double accumulation = 0.0;
-    const auto& own = mesh.ownership();
-
-    for(std::size_t c=0;c<mesh.n_cells();++c) {
-        source_total += source(c) * geometry.cell_volumes[c];
-        if(controls.dt > 0.0)
-            accumulation += controls.density * controls.cp *
-                geometry.cell_volumes[c] *
+    if (temperature.size()!=mesh.n_cells() || temperature.dimension()!=1 ||
+        old_temperature.size()!=mesh.n_cells() || old_temperature.dimension()!=1 ||
+        source.size()!=mesh.n_cells() || source.dimension()!=1)
+        throw std::invalid_argument("energy balance field size mismatch");
+    validate_energy_implicit_source(source_implicit,mesh.n_cells());
+    const auto flux=reconstruct_energy_face_flux(
+        mesh,geometry,mass_flux,temperature,controls,bcs,face_values,
+        conductivity_values,face_conditions);
+    double net_flux=0.0, source_total=0.0, accumulation=0.0;
+    for (std::size_t c=0;c<mesh.n_cells();++c) {
+        const double physical_source=source(c)+
+            (source_implicit ? (*source_implicit)(c)*temperature(c) : 0.0);
+        source_total+=physical_source*geometry.cell_volumes[c];
+        if (controls.dt>0.0)
+            accumulation+=controls.density*controls.cp*geometry.cell_volumes[c]*
                 (temperature(c)-old_temperature(c))/controls.dt;
     }
-
-    for(std::size_t f=0;f<mesh.n_faces();++f) {
-        const std::size_t o=own.owner(f);
-        const auto nraw=own.neighbour(f);
-        if(nraw>=0) continue;
-
-        const std::size_t patch=geometry.face_patch[f];
-        ScalarBoundaryCondition bc;
-        if(patch<mesh.boundary().n_patches()) {
-            const auto& name=mesh.boundary().patch(patch).name;
-            auto it=bcs.find(name);
-            if(it!=bcs.end()) bc=it->second;
-        }
-
-        const auto& Sf0=geometry.face_area_vectors[f];
-        const auto dvec=geometry.face_centres[f]-geometry.cell_centres[o];
-        const double sign=(Sf0.dot(dvec)>=0.0)?1.0:-1.0;
-        const auto Sf=Sf0*sign;
-        const double area=Sf.mag();
-        const double d=dvec.mag();
-        const double F=mass_flux(f);
-        double Tf=temperature(o);
-        bool face_override_valid=false;
-        if(face_values && patch<mesh.boundary().n_patches()) {
-            const auto it=face_values->values.find(mesh.boundary().patch(patch).name);
-            if(it!=face_values->values.end() && f<it->second.size() &&
-               std::isfinite(it->second[f])) {
-                Tf=it->second[f];
-                face_override_valid=true;
-            }
-        }
-        if(!face_override_valid) {
-            if(bc.type==ScalarBoundaryType::FIXED_VALUE)
-                Tf=bc.value;
-            else if(bc.type==ScalarBoundaryType::FIXED_GRADIENT)
-                Tf=temperature(o)+bc.gradient*d;
-        }
-
-        const double k=conductivity_values
-            ? (*conductivity_values)[o] : controls.conductivity;
-        double conductive=-k*bc.gradient*area;
-        if(d>0.0 && (bc.type==ScalarBoundaryType::FIXED_VALUE || face_override_valid))
-            conductive=-k*(Tf-temperature(o))/d*area;
-        const double convective=controls.cp*F*(F>=0.0 ? temperature(o) : Tf);
-        net_flux += convective + conductive;
-    }
-
-    const double imbalance=accumulation + net_flux - source_total;
+    for (std::size_t f=0;f<mesh.n_faces();++f)
+        if (mesh.ownership().neighbour(f)<0) net_flux+=flux(f);
+    const double imbalance=accumulation+net_flux-source_total;
     const double scale=std::max({1.0,std::abs(accumulation),
-                                 std::abs(net_flux),std::abs(source_total)});
+                                std::abs(net_flux),std::abs(source_total)});
     return std::abs(imbalance)/scale;
 }
 
@@ -269,24 +278,20 @@ inline IndependentEnergyBalance reconstruct_energy_balance(
     const cfdx::core::Field<double,cfdx::core::Location::CELL>& source,
     const EnergySolverControls& controls,
     const ScalarBoundaryConditions& bcs = {},
-    const ScalarBoundaryFaceValues* face_values = nullptr)
+    const ScalarBoundaryFaceValues* face_values = nullptr,
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>* source_implicit = nullptr,
+    const std::vector<double>* conductivity_values = nullptr,
+    const ScalarBoundaryFaceConditions* face_conditions = nullptr)
 {
-    if (temperature.size()!=mesh.n_cells() || old_temperature.size()!=mesh.n_cells() ||
-        source.size()!=mesh.n_cells())
+    if (temperature.size()!=mesh.n_cells() || temperature.dimension()!=1 ||
+        old_temperature.size()!=mesh.n_cells() || old_temperature.dimension()!=1 ||
+        source.size()!=mesh.n_cells() || source.dimension()!=1)
         throw std::invalid_argument("reconstruct_energy_balance: field size mismatch");
 
-    // Rebuild the energy face flux from the final temperature and the
-    // authoritative mass flux. This does not use the assembled matrix.
-    // The scalar reconstruction is performed in temperature units first:
-    // F*T - (k/cp)*grad(T). Multiplication by cp then recovers the
-    // physical energy flux F*cp*T - k*grad(T), including boundary terms.
-    auto temperature_flux = reconstruct_scalar_transport_flux(
-        mesh, geometry, mass_flux, temperature,
-        controls.conductivity / controls.cp, bcs, true,
-        ConvectionScheme::UPWIND, nullptr, nullptr, face_values);
-    for (std::size_t f = 0; f < temperature_flux.size(); ++f)
-        temperature_flux(f) *= controls.cp;
-    const auto& flux = temperature_flux;
+    validate_energy_implicit_source(source_implicit,mesh.n_cells());
+    const auto flux=reconstruct_energy_face_flux(
+        mesh,geometry,mass_flux,temperature,controls,bcs,face_values,
+        conductivity_values,face_conditions);
 
     IndependentEnergyBalance r;
     const auto& own=mesh.ownership();
@@ -309,7 +314,11 @@ inline IndependentEnergyBalance reconstruct_energy_balance(
 
     double sum2=0.0;
     for(std::size_t c=0;c<mesh.n_cells();++c) {
-        const double s=source(c);
+        // The implicit coefficient belongs to the physical volumetric
+        // source Su+Sp*T_final, not just the matrix diagonal. Using Su alone
+        // reports a spurious source even at exact radiative equilibrium.
+        const double s=source(c)+
+            (source_implicit ? (*source_implicit)(c)*temperature(c) : 0.0);
         const double accumulation=controls.dt>0.0
             ? controls.density*controls.cp*geometry.cell_volumes[c]*
               (temperature(c)-old_temperature(c))/controls.dt : 0.0;
@@ -354,15 +363,20 @@ inline EnergySolveResult solve_energy(
     const ScalarBoundaryConditions& bcs = {},
     const ScalarBoundaryFaceValues* face_values = nullptr,
     const cfdx::core::Field<double,cfdx::core::Location::CELL>* source_implicit = nullptr,
-    const cfdx::core::Field<double,cfdx::core::Location::CELL>* previous_time_temperature = nullptr)
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>* previous_time_temperature = nullptr,
+    const std::vector<double>* conductivity_values = nullptr,
+    const ScalarBoundaryFaceConditions* face_conditions = nullptr)
 {
     validate_energy_controls(controls);
-    if(temperature.size()!=mesh.n_cells() || source.size()!=mesh.n_cells())
+    if(temperature.size()!=mesh.n_cells() || temperature.dimension()!=1 ||
+       source.size()!=mesh.n_cells() || source.dimension()!=1)
         throw std::invalid_argument("energy field size mismatch");
 
     validate_energy_implicit_source(source_implicit,mesh.n_cells());
+    validate_energy_conductivity(conductivity_values,mesh.n_cells());
     if (previous_time_temperature &&
-        previous_time_temperature->size()!=mesh.n_cells())
+        (previous_time_temperature->size()!=mesh.n_cells() ||
+         previous_time_temperature->dimension()!=1))
         throw std::invalid_argument("energy previous-time field size mismatch");
 
     EnergySolveResult result;
@@ -380,7 +394,8 @@ inline EnergySolveResult solve_energy(
         controls.energy_balance_tolerance > 0.0 ? controls.energy_balance_tolerance : controls.tolerance;
     for(std::size_t iter=1;iter<=controls.max_iterations;++iter) {
         auto eq=assemble_energy_equation(
-            mesh,geometry,mass_flux,source,T_previous_time,controls,bcs,face_values,source_implicit);
+            mesh,geometry,mass_flux,source,T_previous_time,controls,bcs,face_values,
+            source_implicit,conductivity_values,face_conditions);
         cfdx::core::Vector candidate(temperature.size(),0.0);
         for(std::size_t i=0;i<temperature.size();++i)
             candidate(i)=temperature(i);
@@ -466,7 +481,8 @@ inline EnergySolveResult solve_energy(
             relative_temperature_change <= temperature_tolerance;
         const auto independent_energy = reconstruct_energy_balance(
             mesh, geometry, mass_flux, temperature, T_previous_time,
-            source, controls, bcs, face_values);
+            source, controls, bcs, face_values,source_implicit,
+            conductivity_values,face_conditions);
         EnergyIteration iteration_report;
         iteration_report.iteration = iter;
         iteration_report.residual = res;
diff --git a/src/cfdx/physics/finite_volume_transport.h b/src/cfdx/physics/finite_volume_transport.h
index 1f9f632d..9bd43d24 100644
--- a/src/cfdx/physics/finite_volume_transport.h
+++ b/src/cfdx/physics/finite_volume_transport.h
@@ -151,11 +151,36 @@ inline FvGeometry build_fv_geometry(const cfdx::core::Mesh& mesh)
     return g;
 }
 
+// Method contract: docs/audits/01-diffusion-laplacien.md
+struct FluxImbalanceMetrics {
+    double signed_global = 0.0;
+    double global_l1 = 0.0;
+    double local_linf = 0.0;
+};
+
+inline FluxImbalanceMetrics flux_imbalance_metrics(
+    const std::vector<double>& cell_flux_residuals)
+{
+    FluxImbalanceMetrics metrics;
+    for (const double residual : cell_flux_residuals) {
+        if (!std::isfinite(residual))
+            throw std::invalid_argument("flux_imbalance_metrics: non-finite residual");
+        metrics.signed_global += residual;
+        metrics.global_l1 += std::abs(residual);
+        metrics.local_linf = std::max(metrics.local_linf, std::abs(residual));
+    }
+    return metrics;
+}
+
 struct ScalarEquation {
     cfdx::core::SparseMatrix matrix;
     cfdx::core::Vector rhs;
     std::vector<double> diagonal;
+    // Compatibility alias: now represents the maximum cell-local imbalance.
     double max_imbalance = 0.0;
+    double signed_global_imbalance = 0.0;
+    double global_l1_imbalance = 0.0;
+    double max_local_imbalance = 0.0;
     // Assembly diagnostics retained temporarily for Issue #377 root-cause analysis.
     double max_abs_face_flux = 0.0;
     double min_internal_distance = std::numeric_limits<double>::infinity();
@@ -546,9 +571,11 @@ inline ScalarEquation assemble_scalar_equation(
     eq.min_face_area = min_face_area;
     eq.max_face_area = max_face_area;
 
-    double total_flux = 0.0;
-    for (std::size_t c = 0; c < nc; ++c) total_flux += div_phi[c];
-    eq.max_imbalance = std::abs(total_flux);
+    const auto imbalance = flux_imbalance_metrics(div_phi);
+    eq.signed_global_imbalance = imbalance.signed_global;
+    eq.global_l1_imbalance = imbalance.global_l1;
+    eq.max_local_imbalance = imbalance.local_linf;
+    eq.max_imbalance = imbalance.local_linf;
     return eq;
 }
 
@@ -746,6 +773,7 @@ reconstruct_scalar_transport_flux(
     return result;
 }
 
+// Method contract: docs/audits/01-diffusion-laplacien.md
 inline cfdx::core::SolverResult solve_scalar_equation(
     const ScalarEquation& equation,
     cfdx::core::Vector& solution,
@@ -768,24 +796,41 @@ inline cfdx::core::SolverResult solve_scalar_equation(
         const double candidate_value = equation.rhs(0) / diagonal;
         const double candidate_residual =
             std::abs(diagonal * candidate_value - equation.rhs(0));
-        if (candidate_residual <= controls.tolerance) {
-            // A scalar 1x1 system has an exact closed-form solution. Do not
-            // under-relax an already converged linear solve: doing so makes
-            // the returned state fail the requested residual even though the
-            // reported linear residual is zero.
+        // H04: the tolerance is relative, as on the Krylov paths. Keep the
+        // candidate/returned-state diagnostics introduced by PR #1088.
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
+            cfdx::core::SolverResult result{
                 cfdx::core::SolverStatus::CONVERGED,
-                1, candidate_residual, candidate_residual
+                1, candidate_residual, relative
             };
+            result.residual_before_relaxation = candidate_residual;
+            result.residual_after_relaxation = candidate_residual;
+            result.residual_relative_before_relaxation = relative;
+            result.residual_relative_after_relaxation = result.residual_relative_before_relaxation;
+            return result;
         }
-        solution(0) += controls.relaxation * (candidate_value - solution(0));
+        // A rejected predictor must not replace the accepted iterate.
         const double residual =
             std::abs(diagonal * solution(0) - equation.rhs(0));
-        return {
+        cfdx::core::SolverResult result{
             cfdx::core::SolverStatus::MAX_ITER_REACHED,
-            controls.max_iterations, residual, residual
+            1, residual, rhs_scale > 0.0 ? residual / rhs_scale : residual
         };
+        result.residual_before_relaxation = candidate_residual;
+        result.residual_after_relaxation = residual;
+        result.residual_relative_before_relaxation = relative;
+        result.residual_relative_after_relaxation = result.residual_relative;
+        return result;
     }
 
     cfdx::core::Vector candidate = solution;
@@ -940,8 +985,36 @@ inline cfdx::core::SolverResult solve_scalar_equation(
             result.status = cfdx::core::SolverStatus::DIVERGED;
             return result;
         }
+        // These paired diagnostics use the same physical operator, L-infinity
+        // norm, and RHS normalization. Solver-native Krylov diagnostics remain
+        // available in residual/residual_relative and are intentionally separate.
+        const auto physical_residual_inf = [&](const cfdx::core::Vector& state) {
+            double norm = 0.0;
+            for (std::size_t i = 0; i < state.size(); ++i) {
+                double ri = -equation.rhs(i);
+                const auto row_begin = equation.matrix.row_offsets_data()[i];
+                const auto row_end = equation.matrix.row_offsets_data()[i + 1];
+                for (std::uint32_t k = row_begin; k < row_end; ++k)
+                    ri += equation.matrix.values_data()[k] *
+                          state(equation.matrix.columns_data()[k]);
+                norm = std::max(norm, std::abs(ri));
+            }
+            return norm;
+        };
+        double rhs_norm = 0.0;
+        for (std::size_t i = 0; i < solution.size(); ++i)
+            rhs_norm = std::max(rhs_norm, std::abs(equation.rhs(i)));
+        const double rhs_scale = std::max(rhs_norm, 1e-300);
+        result.residual_before_relaxation = physical_residual_inf(candidate);
+        result.residual_relative_before_relaxation =
+            result.residual_before_relaxation / rhs_scale;
+
         for (std::size_t i = 0; i < solution.size(); ++i)
             solution(i) += controls.relaxation * (candidate(i) - solution(i));
+
+        result.residual_after_relaxation = physical_residual_inf(solution);
+        result.residual_relative_after_relaxation =
+            result.residual_after_relaxation / rhs_scale;
     }
     // A MAX_ITER_REACHED result is not an acceptable predictor solution.
     // Do not inject an unconverged Krylov iterate into the nonlinear solver:
diff --git a/src/cfdx/physics/m1_m4_models.h b/src/cfdx/physics/m1_m4_models.h
index 3da07cf4..5b757386 100644
--- a/src/cfdx/physics/m1_m4_models.h
+++ b/src/cfdx/physics/m1_m4_models.h
@@ -269,16 +269,32 @@ inline void validate_reciprocity(const std::vector<double>& F,
 
 inline double p1_blackbody_intensity(double T)
 {
-    if (T < 0.0) throw std::invalid_argument("temperature must be non-negative");
-    return blackbody(T) / M_PI;
+    if (!std::isfinite(T) || T < 0.0)
+        throw std::invalid_argument("P1 temperature must be finite and non-negative");
+    // Apply the small prefactor before T^4 to avoid overflowing pow(T,4)
+    // when the physical intensity itself is still representable.
+    const double intensity = (sigma_sb / M_PI) * T * T * T * T;
+    if (!std::isfinite(intensity))
+        throw std::overflow_error("P1 blackbody intensity is not representable");
+    return intensity;
 }
 
+// Method contract: docs/audits/11-rayonnement.md
 inline double p1_source(double absorption, double mean_intensity, double T)
 {
-    if (absorption < 0.0 || mean_intensity < 0.0 || T < 0.0)
-        throw std::invalid_argument("invalid P1 source input");
-    // P1 source convention: S_r = 4*kappa*(I_b - J).
-    return 4.0 * absorption * (p1_blackbody_intensity(T) - mean_intensity);
+    if (!std::isfinite(absorption) || !std::isfinite(mean_intensity) ||
+        !std::isfinite(T) || absorption < 0.0 || mean_intensity < 0.0 || T < 0.0)
+        throw std::invalid_argument("P1 source inputs must be finite and non-negative");
+    if (absorption == 0.0) return 0.0;
+    // P1 source convention: div(q_r) = 4*pi*kappa*(I_b - J),
+    // where I_b = sigma*T^4/pi and J is mean intensity (G = 4*pi*J).
+    // Form the intensity difference first: equilibrium remains exactly zero
+    // even for a very large, finite absorption coefficient.
+    const double source = 4.0 * M_PI *
+        (absorption * (p1_blackbody_intensity(T) - mean_intensity));
+    if (!std::isfinite(source))
+        throw std::overflow_error("P1 radiative source is not representable");
+    return source;
 }
 
 struct Direction {
diff --git a/src/cfdx/physics/momentum_viscous_stress.h b/src/cfdx/physics/momentum_viscous_stress.h
new file mode 100644
index 00000000..fadcc1a2
--- /dev/null
+++ b/src/cfdx/physics/momentum_viscous_stress.h
@@ -0,0 +1,111 @@
+#pragma once
+
+#include "cfdx/physics/scalar_transport_gradient.h"
+#include <array>
+
+namespace cfdx::physics {
+
+// Owner-oriented flux of mu_eff * transpose(grad(U)). The scalar diffusion
+// operator already supplies mu_eff * grad(U). This complementary conservative
+// flux is required when mu_eff varies, even for divergence-free velocity.
+// Internal mu uses the same normal-distance harmonic interpolation as scalar
+// diffusion; boundary mu and face gradients use their owner-cell values.
+// Boundary U values enter the Green--Gauss reconstruction explicitly. No
+// non-orthogonal/skew correction or higher-order wall extrapolation is implied.
+inline cfdx::core::Field<double,cfdx::core::Location::FACE>
+transpose_viscous_stress_flux(
+    const cfdx::core::Mesh& mesh,
+    const FvGeometry& geometry,
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>& velocity,
+    const std::vector<double>& cell_mu,
+    const std::array<ScalarBoundaryConditions,3>& velocity_bcs)
+{
+    using namespace cfdx::core;
+    const auto nc=mesh.n_cells(), nf=mesh.n_faces();
+    if (velocity.size()!=nc || velocity.dimension()!=3 || cell_mu.size()!=nc)
+        throw std::invalid_argument("transpose viscous stress: invalid field shape");
+    for (const double mu : cell_mu)
+        if (!std::isfinite(mu) || mu<0.0)
+            throw std::invalid_argument("transpose viscous stress: invalid dynamic viscosity");
+    std::array<Field<double,Location::CELL>,3> gradients;
+    for (std::size_t component=0; component<3; ++component) {
+        Field<double,Location::CELL> scalar(nc,"U_stress_component","m/s",1);
+        for (std::size_t cell=0; cell<nc; ++cell) scalar(cell)=velocity(cell,component);
+        gradients[component]=scalar_gradient_with_transport_boundaries(
+            scalar,mesh,geometry,velocity_bcs[component]);
+    }
+    Field<double,Location::FACE> flux(nf,"transpose_viscous_stress_flux","N",3);
+    flux.fill(0.0);
+    for (std::size_t face=0; face<nf; ++face) {
+        const auto owner=mesh.ownership().owner(face);
+        const auto neighbour=mesh.ownership().neighbour(face);
+        const auto patch=geometry.face_patch[face];
+        if (neighbour<0 && patch<mesh.boundary().n_patches() &&
+            mesh.boundary().patch(patch).type==PatchType::EMPTY)
+            continue;
+        double mu=cell_mu[owner], weight=0.0;
+        const auto& Sf=geometry.face_area_vectors[face];
+        if (neighbour>=0) {
+            const auto nb=static_cast<std::size_t>(neighbour);
+            const double area=Sf.mag();
+            const double dp=std::abs((geometry.face_centres[face]-geometry.cell_centres[owner]).dot(Sf))/area;
+            const double dn=std::abs((geometry.cell_centres[nb]-geometry.face_centres[face]).dot(Sf))/area;
+            if (!std::isfinite(dp) || !std::isfinite(dn) || !(dp>0.0) || !(dn>0.0))
+                throw std::invalid_argument("transpose viscous stress: invalid normal half-distance");
+            mu=mu>0.0 && cell_mu[nb]>0.0 ? (dp+dn)/(dp/mu+dn/cell_mu[nb]) : 0.0;
+            weight=geometric_interpolation_weight(
+                geometry.face_centres[face],geometry.cell_centres[owner],&geometry.cell_centres[nb]);
+        }
+        const std::array<double,3> areas{Sf.x,Sf.y,Sf.z};
+        for (std::size_t component=0; component<3; ++component) {
+            double contraction=0.0;
+            for (std::size_t direction=0; direction<3; ++direction) {
+                double gradient=gradients[direction](owner,component);
+                if (neighbour>=0)
+                    gradient=(1.0-weight)*gradient+
+                        weight*gradients[direction](static_cast<std::size_t>(neighbour),component);
+                contraction+=gradient*areas[direction];
+            }
+            flux(face,component)=mu*contraction;
+            if (!std::isfinite(flux(face,component)))
+                throw std::runtime_error("transpose viscous stress: non-finite face flux");
+        }
+    }
+    return flux;
+}
+
+inline cfdx::core::Field<double,cfdx::core::Location::CELL>
+viscous_face_flux_divergence(
+    const cfdx::core::Mesh& mesh,
+    const FvGeometry& geometry,
+    const cfdx::core::Field<double,cfdx::core::Location::FACE>& flux)
+{
+    using namespace cfdx::core;
+    if (flux.size()!=mesh.n_faces() || flux.dimension()!=3 ||
+        geometry.cell_volumes.size()!=mesh.n_cells())
+        throw std::invalid_argument("viscous face divergence: invalid field shape");
+    Field<double,Location::CELL> result(mesh.n_cells(),"transpose_viscous_force","N/m3",3);
+    result.fill(0.0);
+    for (std::size_t face=0; face<mesh.n_faces(); ++face) {
+        const auto owner=mesh.ownership().owner(face);
+        const auto neighbour=mesh.ownership().neighbour(face);
+        if (owner>=mesh.n_cells() ||
+            (neighbour>=0 && static_cast<std::size_t>(neighbour)>=mesh.n_cells()))
+            throw std::invalid_argument("viscous face divergence: invalid face ownership");
+        for (std::size_t component=0; component<3; ++component) {
+            if (!std::isfinite(flux(face,component)))
+                throw std::invalid_argument("viscous face divergence: non-finite flux");
+            result(owner,component)+=flux(face,component);
+            if (neighbour>=0) result(static_cast<std::size_t>(neighbour),component)-=flux(face,component);
+        }
+    }
+    for (std::size_t cell=0; cell<mesh.n_cells(); ++cell) {
+        if (!(geometry.cell_volumes[cell]>0.0) || !std::isfinite(geometry.cell_volumes[cell]))
+            throw std::invalid_argument("viscous face divergence: invalid cell volume");
+        for (std::size_t component=0; component<3; ++component)
+            result(cell,component)/=geometry.cell_volumes[cell];
+    }
+    return result;
+}
+
+} // namespace cfdx::physics
diff --git a/src/cfdx/physics/production_buoyancy_stage.h b/src/cfdx/physics/production_buoyancy_stage.h
index e07cf05b..0666b897 100644
--- a/src/cfdx/physics/production_buoyancy_stage.h
+++ b/src/cfdx/physics/production_buoyancy_stage.h
@@ -372,7 +372,8 @@ inline ProductionBuoyancyOuterLoopResult run_production_buoyant_outer_loop(
     const cfdx::io::CaseSetup& setup,
     cfdx::core::Field<double, cfdx::core::Location::CELL>& velocity,
     cfdx::core::Field<double, cfdx::core::Location::CELL>& pressure,
-    const std::string& restart_path = {})
+    const std::string& restart_path = {},
+    const cfdx::core::Field<double, cfdx::core::Location::CELL>* initial_temperature_field = nullptr)
 {
     if (!production_buoyancy_requested(setup))
         throw std::invalid_argument(
@@ -385,8 +386,23 @@ inline ProductionBuoyancyOuterLoopResult run_production_buoyant_outer_loop(
     const FvGeometry geometry = build_fv_geometry(mesh);
     const std::size_t n_cells = mesh.n_cells();
 
+    // An explicitly restored thermal field owns the first buoyancy source.
+    // Validate it completely before the carrier can modify velocity/pressure.
+    if (initial_temperature_field) {
+        if (initial_temperature_field->size()!=n_cells ||
+            initial_temperature_field->dimension()!=1)
+            throw std::invalid_argument(
+                "production buoyancy stage: initial temperature field must be a scalar matching the mesh");
+        for (std::size_t cell=0;cell<n_cells;++cell)
+            if (!std::isfinite((*initial_temperature_field)(cell)) ||
+                !((*initial_temperature_field)(cell)>0.0))
+                throw std::invalid_argument(
+                    "production buoyancy stage: initial temperature field must be finite and positive (kelvin)");
+    }
+
     // The temperature seen by the buoyancy source. Seeded with the
-    // declared initial temperature, else the case reference temperature
+    // restored field when supplied; otherwise the declared initial
+    // temperature, else the case reference temperature
     // (the stage-1 convention of the frozen-flow path; at the reference
     // temperature the buoyancy source is identically zero, a consistent
     // cold start).
@@ -394,15 +410,16 @@ inline ProductionBuoyancyOuterLoopResult run_production_buoyant_outer_loop(
     if (setup.has_initial_condition &&
         setup.initial_condition.temperature != 0.0)
         initial_temperature = setup.initial_condition.temperature;
-    if (!std::isfinite(initial_temperature) ||
-        !(initial_temperature > 0.0))
+    if (!initial_temperature_field &&
+        (!std::isfinite(initial_temperature) || !(initial_temperature > 0.0)))
         throw std::invalid_argument(
             "production buoyancy stage: initial temperature must be finite "
             "and positive (kelvin)");
     cfdx::core::Field<double, cfdx::core::Location::CELL> relaxed_temperature(
         n_cells, "temperature", "K", 1);
     for (std::size_t cell = 0; cell < n_cells; ++cell)
-        relaxed_temperature(cell) = initial_temperature;
+        relaxed_temperature(cell) = initial_temperature_field
+            ? (*initial_temperature_field)(cell) : initial_temperature;
 
     // Previous UNRELAXED responses for the stationarity gates. Empty
     // until the first response exists.
diff --git a/src/cfdx/physics/production_energy_stage.h b/src/cfdx/physics/production_energy_stage.h
index 57825f1a..1f2b8d8b 100644
--- a/src/cfdx/physics/production_energy_stage.h
+++ b/src/cfdx/physics/production_energy_stage.h
@@ -48,6 +48,11 @@
 //     mass flux, never a flux reconstructed from the cell velocity, and the
 //     finite-volume geometry must match the mesh (sizes, positive finite
 //     volumes, finite centres and area vectors, valid patch indices).
+//   * An optional frozen cell conductivity (W/(m.K), finite and positive)
+//     overrides the molecular conductivity in assembly and every flux audit.
+//     Wall gradients are resolved per face as q''/lambda_eff_owner; density
+//     and cp remain the declared constants. The caller owns any nonlinear
+//     property update and the accepted turbulence state used to build lambda_t.
 //   * Conservation audits are independent of the assembled matrix: the
 //     energy flux is rebuilt from the final temperature and the same mass
 //     flux (reconstruct_energy_balance), and per-patch fluxes are rebuilt
@@ -61,6 +66,10 @@
 //     1 W), so a local defect cannot hide behind a large domain scale -
 //     AND every temperature/flux/cell residual to be finite. A solution
 //     whose independent audit fails cannot be declared converged.
+//   * Production temperatures are absolute kelvins. The final field must
+//     also be strictly positive; an excessive heat sink can satisfy the
+//     linear equation and conservation while producing an inadmissible
+//     temperature. This separate physical gate never clips the solution.
 
 #include "cfdx/core/field/field.h"
 #include "cfdx/io/cfdx_io/case_schema.h"
@@ -213,6 +222,37 @@ inline ProductionEnergyCase build_production_energy_case(
     return energy_case;
 }
 
+// Gradient-diffusion heat transport for an accepted, frozen turbulent state:
+// lambda_eff = lambda + rho*cp*nu_t/Pr_t. This is an opt-in coupling seam;
+// it does not select a turbulence model or run a flow/thermal outer loop.
+inline std::vector<double> build_effective_thermal_conductivity(
+    const ProductionEnergyCase& energy_case,
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>& turbulent_viscosity,
+    double turbulent_prandtl)
+{
+    if (!(turbulent_prandtl>0.0) || !std::isfinite(turbulent_prandtl))
+        throw std::invalid_argument("turbulent Prandtl number must be finite and positive");
+    if (!(energy_case.conductivity>0.0) || !std::isfinite(energy_case.conductivity) ||
+        !(energy_case.density>0.0) || !std::isfinite(energy_case.density) ||
+        !(energy_case.specific_heat>0.0) || !std::isfinite(energy_case.specific_heat))
+        throw std::invalid_argument("effective conductivity requires positive finite material properties");
+    if (turbulent_viscosity.dimension()!=1)
+        throw std::invalid_argument("turbulent viscosity must be a scalar cell field");
+    std::vector<double> conductivity(turbulent_viscosity.size());
+    for (std::size_t c=0;c<turbulent_viscosity.size();++c) {
+        const double nu_t=turbulent_viscosity(c);
+        if (!(nu_t>=0.0) || !std::isfinite(nu_t))
+            throw std::invalid_argument("turbulent viscosity must be finite and non-negative");
+        const long double value=static_cast<long double>(energy_case.conductivity)+
+            static_cast<long double>(energy_case.density)*energy_case.specific_heat*
+            nu_t/turbulent_prandtl;
+        conductivity[c]=static_cast<double>(value);
+        if (!(conductivity[c]>0.0) || !std::isfinite(conductivity[c]))
+            throw std::invalid_argument("effective conductivity must be finite and positive");
+    }
+    return conductivity;
+}
+
 inline ScalarBoundaryConditions build_production_temperature_bcs(
     const cfdx::core::Mesh& mesh,
     const cfdx::io::CaseSetup& setup,
@@ -367,16 +407,68 @@ inline ScalarBoundaryConditions build_production_temperature_bcs(
     return bcs;
 }
 
+// Resolve the canonical patch conditions against each boundary owner's
+// conductivity. A uniform prescribed heat flux must remain q'' on every
+// face, even when one wall patch touches cells with different lambda_eff.
+inline ScalarBoundaryFaceConditions build_production_temperature_face_conditions(
+    const cfdx::core::Mesh& mesh,
+    const ScalarBoundaryConditions& bcs,
+    const ProductionEnergyCase& energy_case,
+    const std::vector<double>& conductivity_values)
+{
+    validate_energy_conductivity(&conductivity_values,mesh.n_cells());
+    if (!(energy_case.conductivity>0.0) || !std::isfinite(energy_case.conductivity))
+        throw std::invalid_argument("production boundary conductivity must be finite and positive");
+    for (const double k : conductivity_values)
+        if (!(k>0.0))
+            throw std::invalid_argument("production effective conductivity must be positive");
+    ScalarBoundaryFaceConditions faces(mesh.n_faces());
+    for (std::size_t p=0;p<mesh.boundary().n_patches();++p) {
+        const auto& patch=mesh.boundary().patch(p);
+        if (patch.type==cfdx::core::PatchType::EMPTY) continue;
+        const auto it=bcs.find(patch.name);
+        if (it==bcs.end())
+            throw std::invalid_argument("missing production thermal patch condition");
+        for (const auto f : patch.face_ids) {
+            if (f>=mesh.n_faces() || mesh.ownership().neighbour(f)>=0 ||
+                mesh.ownership().owner(f)>=mesh.n_cells())
+                throw std::invalid_argument("invalid production thermal boundary face");
+            auto& resolved=faces.conditions[f];
+            const auto& bc=it->second;
+            if (bc.type==ScalarBoundaryType::FIXED_VALUE) {
+                resolved.type=ScalarBoundaryFaceCondition::Type::FIXED_VALUE;
+                resolved.value=bc.value;
+            } else {
+                resolved.type=ScalarBoundaryFaceCondition::Type::FIXED_GRADIENT;
+                const double gradient=bc.type==ScalarBoundaryType::FIXED_GRADIENT
+                    ? bc.gradient : 0.0;
+                // The patch builder encodes q'' as gradient*lambda. Use
+                // wide intermediate arithmetic to avoid a spurious overflow
+                // when the final owner-specific quotient is representable.
+                resolved.value=static_cast<double>(
+                    static_cast<long double>(gradient)*energy_case.conductivity/
+                    conductivity_values[mesh.ownership().owner(f)]);
+            }
+            if (!std::isfinite(resolved.value))
+                throw std::invalid_argument("production thermal face condition must be finite");
+            faces.specified[f]=1;
+        }
+    }
+    return faces;
+}
+
 struct ProductionEnergyStageResult {
     // Blocking verdict, see run_production_energy_stage: linear convergence
     // AND independent global balance AND independent local balance AND all
-    // results finite. A red independent audit can never yield converged.
+    // results finite AND absolute temperatures positive. Linear convergence,
+    // conservation and physical admissibility are reported separately.
     bool converged = false;
     bool linear_solver_converged = false;
     int linear_solver_status = -1;
     double linear_residual_relative = std::numeric_limits<double>::infinity();
     bool global_energy_balance_passed = false;
     bool local_energy_balance_passed = false;
+    bool temperature_admissibility_passed = false;
     std::size_t iterations = 0;
     double min_temperature = std::numeric_limits<double>::infinity();
     double max_temperature = -std::numeric_limits<double>::infinity();
@@ -402,6 +494,7 @@ struct ProductionEnergyStageResult {
     double max_cell_residual_normalized = 0.0;
     std::size_t worst_cell_normalized = 0;
     std::size_t nonfinite_temperature = 0;
+    std::size_t nonpositive_temperature = 0; // finite T <= 0 K
     std::size_t nonfinite_faces = 0;
     std::size_t nonfinite_cells = 0;
     std::vector<EnergyIteration> history;
@@ -415,7 +508,8 @@ inline ProductionEnergyStageResult run_production_energy_stage(
     const cfdx::io::CaseSetup& setup,
     const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
     cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature,
-    const std::map<std::string, double>& radiative_heat_flux_into_domain = {})
+    const std::map<std::string, double>& radiative_heat_flux_into_domain = {},
+    const std::vector<double>* conductivity_values = nullptr)
 {
     validate_production_energy_model(setup);
     if (!production_energy_requested(setup))
@@ -494,6 +588,10 @@ inline ProductionEnergyStageResult run_production_energy_stage(
     const auto energy_case = build_production_energy_case(setup);
     const auto bcs = build_production_temperature_bcs(
         mesh, setup, energy_case, radiative_heat_flux_into_domain);
+    const auto resolved_faces=conductivity_values
+        ? build_production_temperature_face_conditions(mesh,bcs,energy_case,*conductivity_values)
+        : ScalarBoundaryFaceConditions{};
+    const auto* face_conditions=conductivity_values ? &resolved_faces : nullptr;
 
     cfdx::core::Field<double, cfdx::core::Location::CELL> source(
         mesh.n_cells(), "volumetric_power", "W/m3", 1);
@@ -503,7 +601,8 @@ inline ProductionEnergyStageResult run_production_energy_stage(
     const auto previous_temperature = temperature;
     const auto solve = solve_energy(
         mesh, geometry, mass_flux, temperature, source,
-        energy_case.controls, bcs);
+        energy_case.controls, bcs,nullptr,nullptr,nullptr,
+        conductivity_values,face_conditions);
 
     ProductionEnergyStageResult result;
     result.linear_solver_converged = solve.linear_solver_converged;
@@ -519,6 +618,7 @@ inline ProductionEnergyStageResult run_production_energy_stage(
             ++result.nonfinite_temperature;
             continue;
         }
+        if (!(t>0.0)) ++result.nonpositive_temperature;
         result.min_temperature = std::min(result.min_temperature, t);
         result.max_temperature = std::max(result.max_temperature, t);
     }
@@ -530,7 +630,8 @@ inline ProductionEnergyStageResult run_production_energy_stage(
     // the assembled matrix.
     const auto audit = reconstruct_energy_balance(
         mesh, geometry, mass_flux, temperature, previous_temperature,
-        source, energy_case.controls, bcs);
+        source, energy_case.controls, bcs,nullptr,nullptr,
+        conductivity_values,face_conditions);
     result.independent_conservation_residual = audit.residual;
     result.independent_conservation_normalized = audit.normalized_residual;
     result.boundary_flux_out = audit.boundary_flux;
@@ -541,14 +642,11 @@ inline ProductionEnergyStageResult run_production_energy_stage(
     result.nonfinite_faces = audit.nonfinite_faces;
     result.nonfinite_cells = audit.nonfinite_cells;
 
-    // Per-patch fluxes use the same independent reconstruction (temperature
-    // units first, then scaled back by cp to physical watts).
-    auto energy_flux = reconstruct_scalar_transport_flux(
-        mesh, geometry, mass_flux, temperature,
-        energy_case.conductivity / energy_case.specific_heat,
-        bcs, true, ConvectionScheme::UPWIND, nullptr, nullptr, nullptr);
-    for (std::size_t f = 0; f < energy_flux.size(); ++f)
-        energy_flux(f) *= energy_case.specific_heat;
+    // Every balance uses the same physical face flux, including variable
+    // conductivity and the actual owner-specific applied wall gradients.
+    const auto energy_flux=reconstruct_energy_face_flux(
+        mesh,geometry,mass_flux,temperature,energy_case.controls,bcs,nullptr,
+        conductivity_values,face_conditions);
     const auto& own = mesh.ownership();
     for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
         if (own.neighbour(f) >= 0)
@@ -593,15 +691,23 @@ inline ProductionEnergyStageResult run_production_energy_stage(
         std::isfinite(audit.max_cell_residual_normalized) &&
         audit.max_cell_residual_normalized <= balance_tolerance;
 
+    // Conservation alone does not establish physical admissibility. Keep
+    // the raw temperatures and successful algebraic/balance diagnostics
+    // visible when this independent kelvin-domain gate rejects the result.
+    result.temperature_admissibility_passed =
+        result.nonfinite_temperature == 0 &&
+        result.nonpositive_temperature == 0 &&
+        std::isfinite(result.min_temperature) &&
+        std::isfinite(result.max_temperature) &&
+        result.min_temperature > 0.0;
+
     result.converged =
         result.linear_solver_converged &&
         result.global_energy_balance_passed &&
         result.local_energy_balance_passed &&
-        result.nonfinite_temperature == 0 &&
+        result.temperature_admissibility_passed &&
         result.nonfinite_faces == 0 &&
-        result.nonfinite_cells == 0 &&
-        std::isfinite(result.min_temperature) &&
-        std::isfinite(result.max_temperature);
+        result.nonfinite_cells == 0;
     return result;
 }
 
diff --git a/src/cfdx/physics/production_radiation_s2s_stage.h b/src/cfdx/physics/production_radiation_s2s_stage.h
index 323dee4b..932d9328 100644
--- a/src/cfdx/physics/production_radiation_s2s_stage.h
+++ b/src/cfdx/physics/production_radiation_s2s_stage.h
@@ -31,9 +31,11 @@
 //   * Surface areas are the sum of |face_area_vectors| over the patch's
 //     boundary faces and must be finite and positive.
 //   * The wall temperature of each surface is the |Sf|-weighted average of
-//     the adjacent (owner) cell temperatures of the first-pass energy
-//     field: this is the link between the radiation and the resolved
-//     thermal field. The radiosity system is the production
+//     the boundary-face temperatures reconstructed from the energy BC:
+//     T_w = T_owner + (q_imposed + q_radiative_applied)*d_normal/lambda_owner
+//     for prescribed flux. The applied flux is the one used by the last
+//     energy pass, not the as-yet unapplied radiosity response. One patch
+//     remains one isothermal radiosity surface. The radiosity system is the production
 //     solve_s2s_radiosity with machine controls (tolerance 1e-10, direct
 //     dense solve).
 //   * Sign convention (fiche 11, B3): the net radiative flux per surface is
@@ -287,20 +289,32 @@ inline ProductionS2SRadiationCase build_production_s2s_radiation_case(
     return radiation_case;
 }
 
-// |Sf|-weighted average of the adjacent (owner) cell temperatures over
-// each participating patch: the wall temperature the radiosity system
-// sees. This is the only coupling from the resolved thermal field to the
-// radiation.
+// |Sf|-weighted boundary-face temperature for each participating patch.
+// Boundary states use the same outward normal gradient and normal distance
+// as the finite-volume energy operator. Omitted BCs mean zero gradient (the
+// legacy four-argument helper remains usable for that specific condition).
 inline std::vector<double> production_s2s_wall_temperatures(
     const cfdx::core::Mesh& mesh,
     const FvGeometry& geometry,
     const ProductionS2SRadiationCase& radiation_case,
-    const cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature)
+    const cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature,
+    const ScalarBoundaryConditions& temperature_bcs = {},
+    const ScalarBoundaryFaceConditions* face_conditions = nullptr)
 {
     if (temperature.size() != mesh.n_cells() || temperature.dimension() != 1)
         throw std::invalid_argument(
             "production S2S radiation stage: temperature field does not "
             "match the mesh");
+    if (geometry.face_area_vectors.size() != mesh.n_faces() ||
+        geometry.face_centres.size() != mesh.n_faces() ||
+        geometry.cell_centres.size() != mesh.n_cells() ||
+        geometry.cell_volumes.size() != mesh.n_cells() ||
+        (face_conditions &&
+         (face_conditions->conditions.size() != mesh.n_faces() ||
+          face_conditions->specified.size() != mesh.n_faces())))
+        throw std::invalid_argument(
+            "production S2S radiation stage: wall reconstruction geometry "
+            "or face conditions do not match the mesh");
     const std::size_t n = radiation_case.surface_names.size();
     std::vector<double> wall_temperatures(n, 0.0);
     const auto& own = mesh.ownership();
@@ -319,6 +333,9 @@ inline std::vector<double> production_s2s_wall_temperatures(
                 "' does not match a mesh patch");
         double weighted = 0.0;
         double total_weight = 0.0;
+        ScalarBoundaryCondition patch_bc;
+        const auto bc_it = temperature_bcs.find(patch->name);
+        if (bc_it != temperature_bcs.end()) patch_bc = bc_it->second;
         for (const auto face : patch->face_ids) {
             if (face >= mesh.n_faces() || own.neighbour(face) >= 0)
                 throw std::invalid_argument(
@@ -333,10 +350,43 @@ inline std::vector<double> production_s2s_wall_temperatures(
                 throw std::invalid_argument(
                     "production S2S radiation stage: boundary face owner "
                     "is invalid");
-            weighted += weight * temperature(cell);
+            const auto finite_vector = [](const cfdx::core::Vec3& value) {
+                return std::isfinite(value.x) && std::isfinite(value.y) &&
+                       std::isfinite(value.z);
+            };
+            if (!(weight > 0.0) || !std::isfinite(weight) ||
+                !finite_vector(geometry.face_centres[face]) ||
+                !finite_vector(geometry.cell_centres[cell]) ||
+                !std::isfinite(temperature(cell)) || !(temperature(cell) > 0.0))
+                throw std::invalid_argument(
+                    "production S2S radiation stage: invalid boundary geometry "
+                    "or owner temperature");
+            auto bc = patch_bc;
+            if (face_conditions && face_conditions->has(face)) {
+                const auto& resolved = face_conditions->conditions[face];
+                bc.type = resolved.type == ScalarBoundaryFaceCondition::Type::FIXED_VALUE
+                    ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::FIXED_GRADIENT;
+                bc.value = resolved.value;
+                bc.gradient = resolved.value;
+            }
+            double face_temperature = temperature(cell);
+            if (bc.type == ScalarBoundaryType::FIXED_VALUE) {
+                face_temperature = bc.value;
+            } else if (bc.type == ScalarBoundaryType::FIXED_GRADIENT) {
+                face_temperature += bc.gradient *
+                    boundary_normal_distance(geometry, face, cell);
+            }
+            // Reject each face before averaging: a nonphysical cold face
+            // must not be masked by a hotter face of the same patch.
+            if (!std::isfinite(face_temperature) || !(face_temperature > 0.0))
+                throw std::invalid_argument(
+                    "production S2S radiation stage: reconstructed wall "
+                    "temperature must be finite and positive (kelvin)");
+            weighted += weight * face_temperature;
             total_weight += weight;
         }
-        if (!(total_weight > 0.0) || !std::isfinite(weighted))
+        if (!(total_weight > 0.0) || !std::isfinite(total_weight) ||
+            !std::isfinite(weighted))
             throw std::invalid_argument(
                 "production S2S radiation stage: surface '" +
                 radiation_case.surface_names[i] +
@@ -362,7 +412,7 @@ struct ProductionS2SRadiationResult {
     std::vector<std::string> surface_names;
     std::vector<double> areas;              // m^2
     std::vector<double> emissivities;
-    std::vector<double> wall_temperatures;   // K (T0-based)
+    std::vector<double> wall_temperatures;   // K, reconstructed boundary states
     std::vector<double> radiosity;           // W/m^2, J_i
     std::vector<double> irradiation;         // W/m^2, G_i
     std::vector<double> net_flux;            // W/m^2, q_i = J_i - G_i,
@@ -371,6 +421,9 @@ struct ProductionS2SRadiationResult {
     // Flux injected into the energy equation per patch (W/m^2, positive
     // into the domain): -q_i, the same convention as declared heat_flux.
     std::map<std::string, double> injected_heat_flux;
+    // Radiative map used by the preceding energy solve and wall reconstruction.
+    // Distinct from the new response above until the outer fixed point closes.
+    std::map<std::string, double> applied_heat_flux;
     double total_net_power = 0.0;     // W, sum A_i q_i (closes to ~0)
     double balance_scale = 0.0;       // W, emitted + absorbed, floor 1 W
     double balance_normalized = 0.0;  // |total| / scale
@@ -384,7 +437,9 @@ inline ProductionS2SRadiationResult solve_production_s2s_radiation(
     const cfdx::core::Mesh& mesh,
     const FvGeometry& geometry,
     const cfdx::io::CaseSetup& setup,
-    const cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature)
+    const cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature,
+    const std::map<std::string, double>& applied_radiative_heat_flux_into_domain = {},
+    const std::vector<double>* conductivity_values = nullptr)
 {
     if (!production_s2s_radiation_requested(setup))
         throw std::invalid_argument(
@@ -409,8 +464,16 @@ inline ProductionS2SRadiationResult solve_production_s2s_radiation(
 
     const auto radiation_case =
         build_production_s2s_radiation_case(mesh, geometry, setup);
+    const auto energy_case = build_production_energy_case(setup);
+    const auto bcs = build_production_temperature_bcs(
+        mesh, setup, energy_case, applied_radiative_heat_flux_into_domain);
+    const auto face_conditions = conductivity_values
+        ? build_production_temperature_face_conditions(
+              mesh, bcs, energy_case, *conductivity_values)
+        : ScalarBoundaryFaceConditions{};
     const auto wall_temperatures = production_s2s_wall_temperatures(
-        mesh, geometry, radiation_case, temperature);
+        mesh, geometry, radiation_case, temperature, bcs,
+        conductivity_values ? &face_conditions : nullptr);
 
     S2SControls controls;
     controls.tolerance = 1e-10;
@@ -425,6 +488,11 @@ inline ProductionS2SRadiationResult solve_production_s2s_radiation(
     result.areas = radiation_case.areas;
     result.emissivities = radiation_case.emissivities;
     result.wall_temperatures = wall_temperatures;
+    result.applied_heat_flux = applied_radiative_heat_flux_into_domain;
+    // The first energy pass has no radiative map. Keep the returned surface
+    // contract complete even if an isothermal case converges on that pass.
+    for (const auto& name : radiation_case.surface_names)
+        result.applied_heat_flux.try_emplace(name, 0.0);
     result.radiosity = solution.radiosity;
     result.irradiation = solution.irradiation;
     result.net_flux.resize(radiation_case.surface_names.size(), 0.0);
@@ -510,7 +578,8 @@ inline ProductionS2SEnergyStageResult run_production_s2s_energy_stage(
     const FvGeometry& geometry,
     const cfdx::io::CaseSetup& setup,
     const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
-    cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature)
+    cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature,
+    const std::vector<double>* conductivity_values = nullptr)
 {
     if (!production_s2s_radiation_requested(setup))
         throw std::invalid_argument(
@@ -527,14 +596,15 @@ inline ProductionS2SEnergyStageResult run_production_s2s_energy_stage(
     ProductionS2SEnergyStageResult result;
     const auto pass1 =
         run_production_energy_stage(mesh, geometry, setup, mass_flux,
-                                    temperature);
+                                    temperature, {}, conductivity_values);
     result.energy = pass1;
     result.energy_pass1_converged = pass1.converged;
     if (!pass1.converged)
         return result;
 
     const auto radiation =
-        solve_production_s2s_radiation(mesh, geometry, setup, temperature);
+        solve_production_s2s_radiation(mesh, geometry, setup, temperature,
+                                      {}, conductivity_values);
     result.radiation = radiation;
     result.s2s_radiation_converged = radiation.converged;
     if (!radiation.converged)
@@ -545,18 +615,29 @@ inline ProductionS2SEnergyStageResult run_production_s2s_energy_stage(
     // declared heat_flux by the energy BC builder.
     const auto extra_flux = radiation.injected_heat_flux;
     const auto pass2 = run_production_energy_stage(
-        mesh, geometry, setup, mass_flux, temperature, extra_flux);
+        mesh, geometry, setup, mass_flux, temperature, extra_flux,
+        conductivity_values);
     result.energy = pass2;
     result.energy_pass2_converged = pass2.converged;
     result.converged = pass2.converged;
+    if (!pass2.converged)
+        return result; // preserve the structured energy verdict and rejected field
 
     // Delayed-Picard diagnostic: how much the wall temperatures moved
     // between the radiation input (T0) and the re-solved field (T1). Not
     // a gate; the external energy-radiation iteration is step 3.
     const auto radiation_case =
         build_production_s2s_radiation_case(mesh, geometry, setup);
+    const auto energy_case = build_production_energy_case(setup);
+    const auto bcs = build_production_temperature_bcs(
+        mesh, setup, energy_case, extra_flux);
+    const auto face_conditions = conductivity_values
+        ? build_production_temperature_face_conditions(
+              mesh, bcs, energy_case, *conductivity_values)
+        : ScalarBoundaryFaceConditions{};
     const auto wall_t1 = production_s2s_wall_temperatures(
-        mesh, geometry, radiation_case, temperature);
+        mesh, geometry, radiation_case, temperature, bcs,
+        conductivity_values ? &face_conditions : nullptr);
     for (std::size_t i = 0; i < radiation.wall_temperatures.size(); ++i)
         result.max_wall_temperature_lag =
             std::max(result.max_wall_temperature_lag,
@@ -657,7 +738,8 @@ inline ProductionS2SOuterLoopResult run_production_s2s_energy_outer_loop(
     const FvGeometry& geometry,
     const cfdx::io::CaseSetup& setup,
     const cfdx::core::Field<double, cfdx::core::Location::FACE>& mass_flux,
-    cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature)
+    cfdx::core::Field<double, cfdx::core::Location::CELL>& temperature,
+    const std::vector<double>* conductivity_values = nullptr)
 {
     if (!production_s2s_radiation_requested(setup))
         throw std::invalid_argument(
@@ -677,10 +759,18 @@ inline ProductionS2SOuterLoopResult run_production_s2s_energy_outer_loop(
     // builder once, before any solve).
     const auto radiation_case =
         build_production_s2s_radiation_case(mesh, geometry, setup);
+    const auto energy_case = build_production_energy_case(setup);
+    const auto initial_bcs = build_production_temperature_bcs(
+        mesh, setup, energy_case);
+    const auto initial_face_conditions = conductivity_values
+        ? build_production_temperature_face_conditions(
+              mesh, initial_bcs, energy_case, *conductivity_values)
+        : ScalarBoundaryFaceConditions{};
     // Wall temperatures of the incoming field: the k = 1 stationarity
     // reference.
     auto previous_wall_temperatures = production_s2s_wall_temperatures(
-        mesh, geometry, radiation_case, temperature);
+        mesh, geometry, radiation_case, temperature, initial_bcs,
+        conductivity_values ? &initial_face_conditions : nullptr);
 
     // F(k): the radiative flux map applied to the energy pass. Empty on
     // the first pass (no radiative boundary contribution at all), then the
@@ -697,7 +787,8 @@ inline ProductionS2SOuterLoopResult run_production_s2s_energy_outer_loop(
         record.iteration = iteration;
 
         const auto energy = run_production_energy_stage(
-            mesh, geometry, setup, mass_flux, temperature, applied_flux);
+            mesh, geometry, setup, mass_flux, temperature, applied_flux,
+            conductivity_values);
         result.energy = energy;
         record.energy_converged = energy.converged;
         result.outer_iterations = iteration;
@@ -708,7 +799,7 @@ inline ProductionS2SOuterLoopResult run_production_s2s_energy_outer_loop(
         }
 
         const auto radiation = solve_production_s2s_radiation(
-            mesh, geometry, setup, temperature);
+            mesh, geometry, setup, temperature, applied_flux, conductivity_values);
         result.radiation = radiation;
         record.radiation_converged = radiation.converged;
         if (!radiation.converged) {
@@ -716,8 +807,7 @@ inline ProductionS2SOuterLoopResult run_production_s2s_energy_outer_loop(
             return result;  // fail-closed: red radiosity machine gate
         }
 
-        const auto wall_temperatures = production_s2s_wall_temperatures(
-            mesh, geometry, radiation_case, temperature);
+        const auto& wall_temperatures = radiation.wall_temperatures;
 
         // (3) Radiative fixed-point residual on the UNRELAXED response:
         // the map the final field was solved with must be the map the
@@ -767,7 +857,7 @@ inline ProductionS2SOuterLoopResult run_production_s2s_energy_outer_loop(
             std::isfinite(max_wall_change) && max_wall_change <= tolerance;
         if (result.fixed_point_flux_passed &&
             result.fixed_point_wall_temperature_passed) {
-            result.converged = true;  // exact coupled fixed point
+            result.converged = true;  // coupled fixed point within both tolerances
             return result;
         }
 
diff --git a/src/cfdx/physics/production_solid_conduction_contract.h b/src/cfdx/physics/production_solid_conduction_contract.h
new file mode 100644
index 00000000..712eacbb
--- /dev/null
+++ b/src/cfdx/physics/production_solid_conduction_contract.h
@@ -0,0 +1,58 @@
+#pragma once
+
+#include "cfdx/io/cfdx_io/case_schema.h"
+#include "cfdx/physics/production_turbulence_contract.h"
+#include "cfdx/physics/production_time_contract.h"
+
+#include <stdexcept>
+#include <string>
+
+namespace cfdx::physics {
+
+// The solid executable solves steady conduction only. A declaration must not
+// be silently discarded merely because the conduction subsystem converges.
+// "conjugate" is retained as the documented legacy alias for solid-only cases;
+// it does not authorize a fluid/solid or radiation coupling in this executable.
+inline void validate_production_solid_conduction_selection(
+    const cfdx::io::CaseSetup& setup)
+{
+    if (setup.energy_model != "solid_conduction" &&
+        setup.energy_model != "conjugate")
+        throw std::invalid_argument(
+            "thermal solver: energy_model must be 'solid_conduction' or the "
+            "solid-only legacy alias 'conjugate'");
+
+    const auto inactive = [](const std::string& value) {
+        const auto key = normalize_turbulence_model_key(value);
+        return value.find_first_not_of(" \t\r\n") == std::string::npos ||
+            key == "none" || key == "off" || key == "disabled";
+    };
+    if (!inactive(setup.radiation_model))
+        throw std::invalid_argument(
+            "thermal solver: declared radiation is not coupled to the solid "
+            "conduction executable; use a supported radiation-energy path");
+    if (!inactive(setup.multiphase_model))
+        throw std::invalid_argument(
+            "thermal solver: multiphase physics is not supported by steady solid conduction");
+    validate_production_turbulence_selection(setup);
+
+    validate_production_steady_time_selection(setup, "thermal solver");
+
+    const auto physics = normalize_turbulence_model_key(setup.physics_model);
+    if (physics == "buoyant" || physics.rfind("compressible", 0) == 0 ||
+        physics == "euler")
+        throw std::invalid_argument(
+            "thermal solver: declared flow physics is not solved by steady solid conduction");
+
+    for (const auto& boundary : setup.boundary_conditions) {
+        using cfdx::io::BCType;
+        if (boundary.type != BCType::WALL &&
+            boundary.type != BCType::SYMMETRY && boundary.type != BCType::EMPTY)
+            throw std::invalid_argument(
+                "thermal solver: patch '" + boundary.patch_name +
+                "' declares flow, coupling, periodic or unknown boundary physics "
+                "that steady solid conduction does not implement");
+    }
+}
+
+}  // namespace cfdx::physics
diff --git a/src/cfdx/physics/production_time_contract.h b/src/cfdx/physics/production_time_contract.h
new file mode 100644
index 00000000..157c73a0
--- /dev/null
+++ b/src/cfdx/physics/production_time_contract.h
@@ -0,0 +1,34 @@
+#pragma once
+
+#include "cfdx/io/cfdx_io/case_schema.h"
+#include "cfdx/physics/production_turbulence_contract.h"
+
+#include <stdexcept>
+#include <string>
+
+namespace cfdx::physics {
+
+inline void validate_production_steady_time_selection(
+    const cfdx::io::CaseSetup& setup, const std::string& executable)
+{
+    const auto mode = normalize_turbulence_model_key(setup.solver_mode);
+    const auto temporal = normalize_turbulence_model_key(setup.numerics.transient_scheme);
+    const auto declared = [](const std::string& value) {
+        return value.find_first_not_of(" \t\r\n") != std::string::npos;
+    };
+    bool explicit_temporal = false;
+    for (const auto& entry : setup.numerical_config.entries)
+        explicit_temporal = explicit_temporal ||
+            entry.family == cfdx::core::NumericalMethodFamily::Temporal;
+    for (const auto& entry : setup.numerical_report.resolved)
+        explicit_temporal = explicit_temporal ||
+            entry.family == cfdx::core::NumericalMethodFamily::Temporal;
+    if (setup.transient || explicit_temporal ||
+        (declared(setup.solver_mode) && mode != "steady") ||
+        (declared(setup.numerics.transient_scheme) && temporal != "steady"))
+        throw std::invalid_argument(
+            executable + ": transient evolution is not implemented by this steady executable; "
+            "the declared time integration cannot be silently ignored");
+}
+
+}  // namespace cfdx::physics
diff --git a/src/cfdx/physics/radiation_advanced.h b/src/cfdx/physics/radiation_advanced.h
index 942a92ba..4f8b84d7 100644
--- a/src/cfdx/physics/radiation_advanced.h
+++ b/src/cfdx/physics/radiation_advanced.h
@@ -202,45 +202,70 @@ struct ViewFactorPatch {
     double area = 0.0;
 };
 
+inline void validate_view_factor_patch(const ViewFactorPatch& patch)
+{
+    for (std::size_t axis = 0; axis < 3; ++axis) {
+        if (!std::isfinite(patch.center[axis]) ||
+            !std::isfinite(patch.normal[axis]))
+            throw std::invalid_argument(
+                "view-factor patch centre and normal must be finite");
+    }
+    if (!std::isfinite(patch.area) || patch.area <= 0.0)
+        throw std::invalid_argument(
+            "view-factor patch area must be finite and positive");
+    if (patch.normal[0] == 0.0 && patch.normal[1] == 0.0 &&
+        patch.normal[2] == 0.0)
+        throw std::invalid_argument("view-factor patch normal must be non-zero");
+}
+
 inline double patch_pair_view_factor(
     const ViewFactorPatch& source,
     const ViewFactorPatch& target)
 {
+    // Validate before coincident-centroid or back-facing early returns, so
+    // malformed geometry cannot silently become a zero exchange.
+    validate_view_factor_patch(source);
+    validate_view_factor_patch(target);
     const double rx = target.center[0]-source.center[0];
     const double ry = target.center[1]-source.center[1];
     const double rz = target.center[2]-source.center[2];
-    const double r2 = rx*rx+ry*ry+rz*rz;
-    if (!std::isfinite(r2) || r2 <= 0.0)
+    const double r = std::hypot(rx, ry, rz);
+    if (!std::isfinite(r))
+        throw std::overflow_error("view-factor patch separation is not representable");
+    if (r == 0.0)
         return 0.0;
-    const double r = std::sqrt(r2);
-    const double ns = std::sqrt(source.normal[0]*source.normal[0] +
-                                 source.normal[1]*source.normal[1] +
-                                 source.normal[2]*source.normal[2]);
-    const double nt = std::sqrt(target.normal[0]*target.normal[0] +
-                                 target.normal[1]*target.normal[1] +
-                                 target.normal[2]*target.normal[2]);
-    if (!(source.area > 0.0) || !(target.area > 0.0) ||
-        !(ns > 0.0) || !(nt > 0.0))
-        throw std::invalid_argument("invalid view-factor patch");
-    const double cs = (source.normal[0]*rx + source.normal[1]*ry +
-                       source.normal[2]*rz)/(ns*r);
-    const double ct = -(target.normal[0]*rx + target.normal[1]*ry +
-                        target.normal[2]*rz)/(nt*r);
+    const auto unit_normal = [](const std::array<double, 3>& normal) {
+        // A normal's arbitrary scale must not cause squared-norm overflow.
+        const double scale = std::max({std::abs(normal[0]),
+                                       std::abs(normal[1]),
+                                       std::abs(normal[2])});
+        std::array<double, 3> unit{
+            normal[0] / scale, normal[1] / scale, normal[2] / scale};
+        const double norm = std::hypot(unit[0], unit[1], unit[2]);
+        for (auto& component : unit) component /= norm;
+        return unit;
+    };
+    const auto ns = unit_normal(source.normal);
+    const auto nt = unit_normal(target.normal);
+    const double cs = ns[0]*(rx/r) + ns[1]*(ry/r) + ns[2]*(rz/r);
+    const double ct = -(nt[0]*(rx/r) + nt[1]*(ry/r) + nt[2]*(rz/r));
     if (cs <= 0.0 || ct <= 0.0)
         return 0.0;
-    return std::max(0.0, std::min(1.0, target.area * cs * ct /
-                                        (M_PI * r2)));
+    // The differential-area estimate is capped at one. Logarithms avoid
+    // overflow/underflow in A/r^2 before this cap is applied.
+    const double log_factor = std::log(target.area) + std::log(cs) +
+                              std::log(ct) - std::log(M_PI) - 2.0*std::log(r);
+    return log_factor >= 0.0 ? 1.0 : std::exp(log_factor);
 }
 
+// Method contract: docs/audits/11-rayonnement.md
 inline std::vector<double> estimate_view_factor_matrix(
     const std::vector<ViewFactorPatch>& patches)
 {
     if (patches.empty())
         throw std::invalid_argument("view-factor patch set is empty");
-    for (const auto& patch : patches) {
-        if (!std::isfinite(patch.area) || patch.area <= 0.0)
-            throw std::invalid_argument("view-factor patch area must be finite and positive");
-    }
+    // Also checks a singleton patch set, for which no pair kernel is called.
+    for (const auto& patch : patches) validate_view_factor_patch(patch);
     const std::size_t n = patches.size();
     std::vector<double> F(n*n, 0.0);
     for (std::size_t i=0; i<n; ++i)
@@ -248,10 +273,12 @@ inline std::vector<double> estimate_view_factor_matrix(
             if (i != j)
                 F[i*n+j] = patch_pair_view_factor(patches[i], patches[j]);
 
-    // The raw centroid estimator is not reciprocal on unequal-area patches.
+    // Capping the centroid estimate at one can break reciprocity on
+    // unequal-area patches.
     // First symmetrize the area-weighted exchange Q_ij = A_i F_ij. Then apply
     // one global closure scale; a row-by-row normalization would destroy
-    // reciprocity again. Finally, assign the remaining row fraction to F_ii.
+    // reciprocity again. The S2S consumer contract requires F_ii = 0; any
+    // remaining row fraction represents exchange with unmodeled/open space.
     // Positive patch areas are required explicitly: silently flooring a
     // degenerate geometric area would hide an invalid radiation surface and
     // could corrupt reciprocity by changing the physical area.
@@ -259,13 +286,18 @@ inline std::vector<double> estimate_view_factor_matrix(
     for (std::size_t i=0; i<n; ++i) {
         for (std::size_t j=i+1; j<n; ++j) {
             const double exchange =
-                0.5*(patches[i].area*F[i*n+j] +
-                     patches[j].area*F[j*n+i]);
+                0.5*patches[i].area*F[i*n+j] +
+                0.5*patches[j].area*F[j*n+i];
             F[i*n+j] = exchange/patches[i].area;
             F[j*n+i] = exchange/patches[j].area;
+            if (!std::isfinite(F[i*n+j]) || !std::isfinite(F[j*n+i]))
+                throw std::overflow_error(
+                    "reciprocal view-factor exchange is not representable");
         }
         double row_sum = 0.0;
         for (std::size_t j=0; j<n; ++j) row_sum += F[i*n+j];
+        if (!std::isfinite(row_sum))
+            throw std::overflow_error("view-factor row sum is not representable");
         max_row_sum = std::max(max_row_sum,row_sum);
     }
 
@@ -275,12 +307,9 @@ inline std::vector<double> estimate_view_factor_matrix(
     // established above.
     const double closure_scale = 1.0/std::max(1.0,max_row_sum);
     for (std::size_t i=0; i<n; ++i) {
-        double row_sum = 0.0;
-        for (std::size_t j=0; j<n; ++j) {
-            if (i!=j) F[i*n+j] *= closure_scale;
-            row_sum += F[i*n+j];
-        }
-        F[i*n+i] = std::max(0.0,1.0-row_sum);
+        F[i*n+i] = 0.0;
+        for (std::size_t j=0; j<n; ++j)
+            if (i != j) F[i*n+j] *= closure_scale;
     }
     return F;
 }
@@ -385,6 +414,15 @@ struct RadiationTriangle {
     std::array<double,3> c{0,0,0};
 };
 
+inline void validate_radiation_triangle_coordinates(const RadiationTriangle& triangle)
+{
+    for (std::size_t axis = 0; axis < 3; ++axis)
+        if (!std::isfinite(triangle.a[axis]) ||
+            !std::isfinite(triangle.b[axis]) ||
+            !std::isfinite(triangle.c[axis]))
+            throw std::invalid_argument("radiation triangle coordinates must be finite");
+}
+
 inline std::array<double,3> radiation_cross(
     const std::array<double,3>& a,const std::array<double,3>& b)
 {
@@ -419,19 +457,26 @@ inline std::array<double,3> radiation_scale(
 
 inline double radiation_norm(const std::array<double,3>& a)
 {
-    return std::sqrt(radiation_dot(a,a));
+    return std::hypot(a[0], a[1], a[2]);
 }
 
 inline double radiation_triangle_area(const RadiationTriangle& t)
 {
-    return 0.5*radiation_norm(radiation_cross(
+    validate_radiation_triangle_coordinates(t);
+    const double area = 0.5*radiation_norm(radiation_cross(
         radiation_sub(t.b,t.a),radiation_sub(t.c,t.a)));
+    if (!std::isfinite(area))
+        throw std::overflow_error("radiation triangle area is not representable");
+    return area;
 }
 
 inline std::array<double,3> radiation_triangle_normal(const RadiationTriangle& t)
 {
+    validate_radiation_triangle_coordinates(t);
     const auto n=radiation_cross(radiation_sub(t.b,t.a),radiation_sub(t.c,t.a));
     const double l=radiation_norm(n);
+    if (!std::isfinite(l))
+        throw std::overflow_error("radiation triangle normal is not representable");
     if(!(l>0.0)) throw std::invalid_argument("degenerate radiation triangle");
     return radiation_scale(n,1.0/l);
 }
@@ -442,20 +487,32 @@ inline bool radiation_ray_triangle_hit(
     const RadiationTriangle& tri,
     double& distance)
 {
+    validate_radiation_triangle_coordinates(tri);
+    for (std::size_t axis = 0; axis < 3; ++axis)
+        if (!std::isfinite(origin[axis]) || !std::isfinite(direction[axis]))
+            throw std::invalid_argument("radiation ray origin and direction must be finite");
     constexpr double eps=1e-11;
     const auto e1=radiation_sub(tri.b,tri.a);
     const auto e2=radiation_sub(tri.c,tri.a);
     const auto h=radiation_cross(direction,e2);
     const double det=radiation_dot(e1,h);
+    if (!std::isfinite(det))
+        throw std::overflow_error("radiation ray determinant is not representable");
     if(std::abs(det)<eps) return false;
     const double inv=1.0/det;
     const auto s=radiation_sub(origin,tri.a);
     const double u=inv*radiation_dot(s,h);
+    if (!std::isfinite(u))
+        throw std::overflow_error("radiation ray barycentric coordinate is not representable");
     if(u<0.0 || u>1.0) return false;
     const auto q=radiation_cross(s,e1);
     const double v=inv*radiation_dot(direction,q);
+    if (!std::isfinite(v))
+        throw std::overflow_error("radiation ray barycentric coordinate is not representable");
     if(v<0.0 || u+v>1.0) return false;
     const double t=inv*radiation_dot(e2,q);
+    if (!std::isfinite(t))
+        throw std::overflow_error("radiation ray distance is not representable");
     if(t<=eps) return false;
     distance=t;
     return true;
@@ -479,11 +536,9 @@ struct RadiationAabb {
 
 inline RadiationAabb radiation_triangle_aabb(const RadiationTriangle& t)
 {
+    validate_radiation_triangle_coordinates(t);
     RadiationAabb b;
     for (int k=0;k<3;++k) {
-        if (!std::isfinite(t.a[k]) || !std::isfinite(t.b[k]) ||
-            !std::isfinite(t.c[k]))
-            throw std::invalid_argument("radiation triangle coordinates must be finite");
         b.lo[k]=std::min({t.a[k],t.b[k],t.c[k]});
         b.hi[k]=std::max({t.a[k],t.b[k],t.c[k]});
     }
@@ -664,6 +719,8 @@ inline double estimate_view_factor_ray_traced(
     RadiationBvh blocker_bvh(blockers);
     double source_area=0.0;
     for(const auto& t:source) source_area+=radiation_triangle_area(t);
+    if (!std::isfinite(source_area))
+        throw std::overflow_error("radiation source surface area is not representable");
     if(!(source_area>0.0)) throw std::invalid_argument("source surface has zero area");
 
     std::vector<double> cumulative;
diff --git a/src/cfdx/physics/scalar_transport_gradient.h b/src/cfdx/physics/scalar_transport_gradient.h
new file mode 100644
index 00000000..4c32fcae
--- /dev/null
+++ b/src/cfdx/physics/scalar_transport_gradient.h
@@ -0,0 +1,114 @@
+#pragma once
+
+#include "cfdx/physics/finite_volume_transport.h"
+
+namespace cfdx::physics {
+
+// Explicit face data is an all-faces contract for each named patch. Entries
+// outside that patch may be NaN; no value on an actual patch face may be absent.
+// Inspect mesh face_ids, not just the supplied vectors or geometry lookup.
+inline void validate_transport_boundary_face_values(
+    const cfdx::core::Mesh& mesh,
+    const FvGeometry& geometry,
+    const ScalarBoundaryFaceValues& values,
+    bool require_positive = false)
+{
+    if (geometry.face_patch.size() != mesh.n_faces())
+        throw std::invalid_argument("boundary face values: geometry shape mismatch");
+    for (const auto& [name, data] : values.values) {
+        const auto index = mesh.boundary().find(name);
+        if (index == mesh.boundary().n_patches() || data.size() != mesh.n_faces())
+            throw std::invalid_argument("boundary face values: unknown patch or invalid global-face shape: " + name);
+        const auto& patch = mesh.boundary().patch(index);
+        if (patch.face_ids.empty() || patch.type == cfdx::core::PatchType::EMPTY)
+            throw std::invalid_argument("boundary face values: empty patch: " + name);
+        for (const auto face : patch.face_ids) {
+            if (face >= mesh.n_faces() || geometry.face_patch[face] != index ||
+                mesh.ownership().neighbour(face) >= 0 ||
+                mesh.ownership().owner(face) >= mesh.n_cells())
+                throw std::invalid_argument("boundary face values: invalid mesh membership: " + name);
+            if (!std::isfinite(data[face]) || (require_positive && data[face] <= 0.0))
+                throw std::invalid_argument("boundary face values: missing or invalid value at face " +
+                                            std::to_string(face) + " on patch " + name);
+        }
+    }
+}
+
+// Green--Gauss reconstruction with exactly the scalar transport boundary
+// states: prescribed value, prescribed outward normal gradient, or owner value
+// for zero gradient. Per-face values override patch values. Internal faces use
+// the same longitudinal geometric interpolation as the core Gauss operator.
+// This does not add a skewness correction or claim affine exactness on skew meshes.
+inline cfdx::core::Field<double, cfdx::core::Location::CELL>
+scalar_gradient_with_transport_boundaries(
+    const cfdx::core::Field<double, cfdx::core::Location::CELL>& field,
+    const cfdx::core::Mesh& mesh,
+    const FvGeometry& geometry,
+    const ScalarBoundaryConditions& bcs = {},
+    const ScalarBoundaryFaceValues* face_values = nullptr)
+{
+    using namespace cfdx::core;
+    const std::size_t nc = mesh.n_cells(), nf = mesh.n_faces();
+    if (field.size() != nc || field.dimension() != 1 ||
+        geometry.cell_centres.size() != nc || geometry.cell_volumes.size() != nc ||
+        geometry.face_centres.size() != nf || geometry.face_area_vectors.size() != nf ||
+        geometry.face_patch.size() != nf)
+        throw std::invalid_argument("scalar transport gradient: field or geometry shape mismatch");
+    if (face_values)
+        validate_transport_boundary_face_values(mesh, geometry, *face_values);
+    std::vector<Vec3> sums(nc);
+    for (std::size_t cell = 0; cell < nc; ++cell) {
+        if (!std::isfinite(field(cell)) || !(geometry.cell_volumes[cell] > 0.0) ||
+            !std::isfinite(geometry.cell_volumes[cell]))
+            throw std::invalid_argument("scalar transport gradient: invalid field or cell volume");
+    }
+    for (std::size_t face = 0; face < nf; ++face) {
+        const std::size_t owner = mesh.ownership().owner(face);
+        const auto neighbour = mesh.ownership().neighbour(face);
+        if (owner >= nc || (neighbour >= 0 && static_cast<std::size_t>(neighbour) >= nc))
+            throw std::invalid_argument("scalar transport gradient: invalid face ownership");
+        double value = field(owner);
+        if (neighbour >= 0) {
+            const auto nb = static_cast<std::size_t>(neighbour);
+            const double weight = geometric_interpolation_weight(
+                geometry.face_centres[face], geometry.cell_centres[owner],
+                &geometry.cell_centres[nb]);
+            value = (1.0 - weight) * value + weight * field(nb);
+        } else {
+            const auto patch = geometry.face_patch[face];
+            if (patch < mesh.boundary().n_patches()) {
+                const auto& boundary = mesh.boundary().patch(patch);
+                if (boundary.type == PatchType::EMPTY)
+                    continue;
+                ScalarBoundaryCondition bc;
+                const auto entry = bcs.find(boundary.name);
+                if (entry != bcs.end()) bc = entry->second;
+                if (face_values) {
+                    const auto data = face_values->values.find(boundary.name);
+                    if (data != face_values->values.end())
+                        bc = {ScalarBoundaryType::FIXED_VALUE, data->second[face], 0.0};
+                }
+                if (bc.type == ScalarBoundaryType::FIXED_VALUE)
+                    value = bc.value;
+                else if (bc.type == ScalarBoundaryType::FIXED_GRADIENT)
+                    value += bc.gradient * boundary_normal_distance(geometry, face, owner);
+            }
+        }
+        const auto contribution = geometry.face_area_vectors[face] * value;
+        if (!std::isfinite(contribution.x) || !std::isfinite(contribution.y) ||
+            !std::isfinite(contribution.z))
+            throw std::invalid_argument("scalar transport gradient: non-finite face contribution");
+        sums[owner] += contribution;
+        if (neighbour >= 0)
+            sums[static_cast<std::size_t>(neighbour)] += contribution * -1.0;
+    }
+    Field<double, Location::CELL> result(nc, field.name() + "_grad_bc", "", 3);
+    for (std::size_t cell = 0; cell < nc; ++cell) {
+        result(cell, 0) = sums[cell].x / geometry.cell_volumes[cell];
+        result(cell, 1) = sums[cell].y / geometry.cell_volumes[cell];
+        result(cell, 2) = sums[cell].z / geometry.cell_volumes[cell];
+    }
+    return result;
+}
+
+} // namespace cfdx::physics
diff --git a/src/cfdx/physics/sst_production_coupling.h b/src/cfdx/physics/sst_production_coupling.h
index d396ef36..15bcd2ec 100644
--- a/src/cfdx/physics/sst_production_coupling.h
+++ b/src/cfdx/physics/sst_production_coupling.h
@@ -293,13 +293,14 @@ inline SSTProductionStepResult advance_sst_transport_transactionally(
             return step;
         }
 
-        const auto cross = cell_gradient_dot(k, omega, mesh);
+        // F2 contains no cross diffusion; no scalar-boundary gradient is
+        // needed to refresh nu_t. F1/transport use the true BCs in sst_solver.
         cfdx::core::Field<double, cfdx::core::Location::CELL> F2(
             n, "sst_F2", "1", 1);
         for (std::size_t i = 0; i < n; ++i) {
             const auto blend = compute_sst_blending(
                 k(i), omega(i), wall_distance(i), controls.molecular_viscosity,
-                controls.beta_star, cross[i], controls.density,
+                controls.beta_star, 0.0, controls.density,
                 controls.sst_sigma_w2);
             F2(i) = blend.second;
         }
@@ -375,12 +376,11 @@ inline cfdx::io::DatRestartState restore_sst_production_state_from_dat(
                 "SST restart requires finite positive current-mesh wall distance");
     }
 
-    const auto cross = cell_gradient_dot(candidate.k, candidate.omega, mesh);
     Field<double, Location::CELL> F2(n, "sst_F2", "1", 1);
     for (std::size_t cell = 0; cell < n; ++cell) {
         const auto blending = compute_sst_blending(
             candidate.k(cell), candidate.omega(cell), candidate.wall_distance(cell),
-            controls.molecular_viscosity, controls.beta_star, cross[cell],
+            controls.molecular_viscosity, controls.beta_star, 0.0,
             controls.density, controls.sst_sigma_w2);
         F2(cell) = blending.second;
     }
@@ -426,13 +426,12 @@ inline SSTProductionState initialize_sst_production_state(
         state.strain_rate(cell) = strain_rate[cell];
     }
 
-    const auto cross = cell_gradient_dot(state.k, state.omega, mesh);
     cfdx::core::Field<double, cfdx::core::Location::CELL> F2(
         n, "sst_F2", "1", 1);
     for (std::size_t cell = 0; cell < n; ++cell) {
         const auto blending = compute_sst_blending(
             state.k(cell), state.omega(cell), state.wall_distance(cell),
-            controls.molecular_viscosity, controls.beta_star, cross[cell],
+            controls.molecular_viscosity, controls.beta_star, 0.0,
             controls.density, controls.sst_sigma_w2);
         F2(cell) = blending.second;
     }
@@ -494,10 +493,9 @@ inline void refresh_sst_production_derived_state_from_mesh(
         candidate.wall_distance(cell) = wall_distance.distance[cell];
         candidate.strain_rate(cell) = strain_field(cell);
     }
-    const auto cross = cell_gradient_dot(candidate.k, candidate.omega, mesh);
     cfdx::core::Field<double, cfdx::core::Location::CELL> F2(n, "sst_F2", "1", 1);
     for (std::size_t cell = 0; cell < n; ++cell) {
-        const auto blending = compute_sst_blending(candidate.k(cell), candidate.omega(cell), candidate.wall_distance(cell), controls.molecular_viscosity, controls.beta_star, cross[cell], controls.density, controls.sst_sigma_w2);
+        const auto blending = compute_sst_blending(candidate.k(cell), candidate.omega(cell), candidate.wall_distance(cell), controls.molecular_viscosity, controls.beta_star, 0.0, controls.density, controls.sst_sigma_w2);
         F2(cell) = blending.second;
     }
     refresh_sst_production_eddy_viscosity(candidate.k, candidate.omega, candidate.strain_rate, F2, controls, candidate.nu_t);
@@ -706,29 +704,69 @@ inline ScalarBoundaryFaceValues make_sst_wall_omega_face_values(
         if (face_counts[name] == 0)
             throw std::invalid_argument("SST wall omega patch has no boundary faces in mesh: " + name);
     }
+    validate_transport_boundary_face_values(mesh, geometry, values, true);
     return values;
 }
 
-// Fail closed before transport assembly until every declared wall-omega
-// condition has a corresponding production-qualified wall operator.
+// A geometry-free call cannot certify global-face coverage. Retain the old
+// overload for contracts without wall omega, but require the mesh-aware API
+// whenever WALL_OMEGA_SST is declared.
 inline void require_supported_sst_runtime_boundaries(
     const SSTRuntimeBoundaryContract& runtime,
     const ScalarBoundaryFaceValues* omega_face_values = nullptr)
 {
+    (void)omega_face_values;
     if (runtime.wall_omega_sst_patches.empty())
         return;
+    throw std::invalid_argument(
+        "SST production transport is blocked: WALL_OMEGA_SST requires mesh-aware "
+        "validation of every boundary face");
+}
+
+// Fail closed before carrier initialization or any SST mutation. Every mesh
+// face on every declared wall must have its own finite positive omega value.
+inline void require_supported_sst_runtime_boundaries(
+    const SSTRuntimeBoundaryContract& runtime,
+    const cfdx::core::Mesh& mesh,
+    const FvGeometry& geometry,
+    const ScalarBoundaryFaceValues* omega_face_values = nullptr)
+{
+    if (omega_face_values)
+        validate_transport_boundary_face_values(mesh, geometry, *omega_face_values, true);
+    if (runtime.wall_omega_sst_patches.empty())
+        return;
+    if (geometry.face_centres.size() != mesh.n_faces() ||
+        geometry.face_area_vectors.size() != mesh.n_faces() ||
+        geometry.cell_centres.size() != mesh.n_cells())
+        throw std::invalid_argument("SST wall omega geometry shape mismatch");
     if (omega_face_values == nullptr)
         throw std::invalid_argument(
             "SST production transport is blocked: WALL_OMEGA_SST requires explicit "
             "per-face wall values");
+    std::map<std::string, bool> seen;
     for (const auto& patch : runtime.wall_omega_sst_patches) {
-        const auto it = omega_face_values->values.find(patch);
-        if (it == omega_face_values->values.end() ||
-            std::none_of(it->second.begin(), it->second.end(),
-                         [](double value) { return std::isfinite(value); }))
+        const auto index = mesh.boundary().find(patch);
+        if (patch.empty() || !seen.emplace(patch, true).second ||
+            index == mesh.boundary().n_patches() ||
+            mesh.boundary().patch(index).type != cfdx::core::PatchType::WALL)
+            throw std::invalid_argument("SST wall omega requires a unique mesh WALL patch: " + patch);
+        if (omega_face_values->values.find(patch) == omega_face_values->values.end())
             throw std::invalid_argument(
                 "SST production transport is blocked: wall omega values are missing "
                 "for patch '" + patch + "'");
+        for (const auto face : mesh.boundary().patch(index).face_ids) {
+            const auto& area_vector = geometry.face_area_vectors[face];
+            const double area = area_vector.mag();
+            const auto owner = mesh.ownership().owner(face);
+            const double distance = area > 0.0
+                ? std::abs((geometry.face_centres[face] - geometry.cell_centres[owner]).dot(area_vector)) / area
+                : 0.0;
+            // WALL_OMEGA_SST has a resolved-wall geometric contract: never
+            // accept a characteristic-length fallback for a degenerate face.
+            if (!(area > 0.0) || !std::isfinite(area) ||
+                !(distance > 0.0) || !std::isfinite(distance))
+                throw std::invalid_argument("SST wall omega requires a positive geometric normal distance: " + patch);
+        }
     }
 }
 
diff --git a/src/cfdx/physics/sst_solver.h b/src/cfdx/physics/sst_solver.h
index e56ad2b5..87f26b0f 100644
--- a/src/cfdx/physics/sst_solver.h
+++ b/src/cfdx/physics/sst_solver.h
@@ -1,6 +1,7 @@
 #pragma once
 
 #include "cfdx/physics/turbulence_solver.h"
+#include "cfdx/physics/scalar_transport_gradient.h"
 #include <algorithm>
 #include <cmath>
 #include <stdexcept>
@@ -8,6 +9,89 @@
 
 namespace cfdx::physics {
 
+inline std::vector<double> sst_boundary_gradient_dot(
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>& k,
+    const cfdx::core::Field<double,cfdx::core::Location::CELL>& omega,
+    const cfdx::core::Mesh& mesh,
+    const FvGeometry& geometry,
+    const ScalarBoundaryConditions& k_bcs,
+    const ScalarBoundaryConditions& omega_bcs,
+    const ScalarBoundaryFaceValues* omega_face_values = nullptr)
+{
+    const auto grad_k = scalar_gradient_with_transport_boundaries(k, mesh, geometry, k_bcs);
+    const auto grad_omega = scalar_gradient_with_transport_boundaries(
+        omega, mesh, geometry, omega_bcs, omega_face_values);
+    std::vector<double> result(mesh.n_cells());
+    for (std::size_t cell = 0; cell < result.size(); ++cell) {
+        result[cell] = 0.0;
+        for (std::size_t component = 0; component < 3; ++component)
+            result[cell] += grad_k(cell, component) * grad_omega(cell, component);
+        if (!std::isfinite(result[cell]))
+            throw std::invalid_argument("SST boundary gradient contraction is non-finite");
+    }
+    return result;
+}
+
+// Residual in solution units. Scaling each row by sum_j |a_ij| makes this
+// comparable with the existing absolute iterate-change tolerance. Unlike a
+// RHS-relative residual, it remains meaningful for a homogeneous equation.
+// The unscaled physical residual is still reported in the transport result.
+inline double sst_equation_solution_defect(
+    const ScalarEquation& equation, const cfdx::core::Vector& solution)
+{
+    if (solution.size() != equation.rhs.size())
+        throw std::invalid_argument("SST residual: field size mismatch");
+    double maximum = 0.0;
+    for (std::size_t cell = 0; cell < solution.size(); ++cell) {
+        double residual = -equation.rhs(cell), scale = 0.0;
+        for (auto entry = equation.matrix.row_offsets_data()[cell];
+             entry < equation.matrix.row_offsets_data()[cell + 1]; ++entry) {
+            const double coefficient = equation.matrix.values_data()[entry];
+            residual += coefficient * solution(equation.matrix.columns_data()[entry]);
+            scale += std::abs(coefficient);
+        }
+        if (!std::isfinite(residual) || !std::isfinite(scale) || !(scale > 0.0))
+            return std::numeric_limits<double>::infinity();
+        maximum = std::max(maximum, std::abs(residual) / scale);
+    }
+    return maximum;
+}
+
+// Translate the absolute solution-unit gate into a sufficient relative
+// Krylov tolerance. Otherwise a warm start can satisfy ||r||2/||b||2 while
+// failing the row-scaled SST gate, so every solve returns the same rejected
+// iterate forever. This only tightens the requested linear tolerance; the
+// external tolerance, iteration budgets and relaxation are unchanged.
+inline double sst_linear_tolerance_for_solution_defect(
+    const ScalarEquation& equation, double tolerance)
+{
+    if (!std::isfinite(tolerance) || tolerance<=0.0 || equation.rhs.size()==0)
+        throw std::invalid_argument("SST linear tolerance: invalid requested tolerance or empty equation");
+    double rhs_norm=0.0;
+    double minimum_row_norm=std::numeric_limits<double>::infinity();
+    for (std::size_t cell=0; cell<equation.rhs.size(); ++cell) {
+        if (!std::isfinite(equation.rhs(cell)))
+            throw std::invalid_argument("SST linear tolerance: non-finite right-hand side");
+        rhs_norm=std::hypot(rhs_norm,equation.rhs(cell));
+        double row_norm=0.0;
+        for (auto entry=equation.matrix.row_offsets_data()[cell];
+             entry<equation.matrix.row_offsets_data()[cell+1]; ++entry)
+            row_norm+=std::abs(equation.matrix.values_data()[entry]);
+        if (!std::isfinite(row_norm) || !(row_norm>0.0))
+            throw std::invalid_argument("SST linear tolerance: zero or non-finite matrix row norm");
+        minimum_row_norm=std::min(minimum_row_norm,row_norm);
+    }
+    if (!std::isfinite(rhs_norm))
+        throw std::invalid_argument("SST linear tolerance: right-hand side norm overflow");
+    // Native CG/BiCGStab use max(||b||2,1e-15). GMRES/dense fallbacks use
+    // smaller RHS scales, so this is sufficient for those paths as well.
+    const double ratio=std::min(1.0,minimum_row_norm/std::max(rhs_norm,1e-15));
+    const double translated=tolerance*ratio;
+    if (!std::isfinite(translated) || !(translated>0.0))
+        throw std::domain_error("SST linear tolerance: requested absolute defect cannot be represented");
+    return translated;
+}
+
 inline std::pair<double,double> compute_sst_blending(
     double k, double omega, double wall_distance,
     double molecular_viscosity, double beta_star,
@@ -58,7 +142,10 @@ inline SSTCellSources sst_cell_sources(
     s.sk=Pk;
     s.spk=-c.density*c.beta_star*wi;
     const double cross=2.0*(1.0-f1)*c.density*c.sst_sigma_w2/wi*grad_k_dot_grad_omega;
-    s.sw=blend(c.gamma1,c.gamma2)*c.density*S*S+std::max(cross,0.0);
+    // SST-2003 uses the limited k-production in the omega production term.
+    // If nut is zero, the limited production is also zero; avoid a 0/0 path.
+    const double omega_production = s.nut > 0.0 ? Pk / s.nut : 0.0;
+    s.sw=blend(c.gamma1,c.gamma2)*omega_production+std::max(cross,0.0);
     s.spw=-c.density*blend(c.beta1,c.beta2)*wi+std::min(cross,0.0)/wi;
     s.gamma_k=c.density*(c.molecular_viscosity+
         blend(c.sst_sigma_k1,c.sst_sigma_k2)*s.nut);
@@ -87,12 +174,20 @@ inline TurbulenceTransportResult solve_sst_transport(
     if(controls.model!=TurbulenceModel::SST)
         throw std::invalid_argument("solve_sst_transport requires SST model");
     const std::size_t n=mesh.n_cells();
-    if(F1.size()!=n || F2.size()!=n || strain_rate.size()!=n)
+    if(F1.size()!=n || F2.size()!=n || strain_rate.size()!=n ||
+       k.size()!=n || omega.size()!=n || mass_flux.size()!=mesh.n_faces() ||
+       k.dimension()!=1 || omega.dimension()!=1 || strain_rate.dimension()!=1 ||
+       F1.dimension()!=1 || F2.dimension()!=1 || mass_flux.dimension()!=1)
         throw std::invalid_argument("SST blending field size mismatch");
+    if (max_iterations == 0 || !std::isfinite(tolerance) || tolerance <= 0.0)
+        throw std::invalid_argument("SST transport requires valid iteration controls");
+    if (omega_face_values)
+        validate_transport_boundary_face_values(mesh, geometry, *omega_face_values, true);
 
     TurbulenceTransportResult result;
     for(std::size_t iter=1;iter<=max_iterations;++iter) {
-        const auto cross = cell_gradient_dot(k,omega,mesh);
+        const auto cross = sst_boundary_gradient_dot(
+            k, omega, mesh, geometry, k_bcs, omega_bcs, omega_face_values);
         auto oldk=k;
         auto oldw=omega;
         cfdx::core::Field<double,cfdx::core::Location::CELL> sk(n,"Sk","W/m3",1),sw(n,"Sw","W/m3",1);
@@ -111,15 +206,18 @@ inline TurbulenceTransportResult solve_sst_transport(
             mesh,geometry,mass_flux,0.0,sw,spw,omega_bcs,true,
             omega_face_values,nullptr,nullptr,&gamma_w);
 
-        ScalarSolveControls sc{2000,tolerance,0.7};
+        const ScalarSolveControls k_sc{
+            2000,sst_linear_tolerance_for_solution_defect(eqk,tolerance),0.7};
+        const ScalarSolveControls omega_sc{
+            2000,sst_linear_tolerance_for_solution_defect(eqw,tolerance),0.7};
         cfdx::core::Vector k_solution(n, 0.0);
         cfdx::core::Vector omega_solution(n, 0.0);
         for (std::size_t i = 0; i < n; ++i) {
             k_solution(i) = k(i);
             omega_solution(i) = omega(i);
         }
-        const auto rk = solve_scalar_equation(eqk, k_solution, sc);
-        const auto rw = solve_scalar_equation(eqw, omega_solution, sc);
+        const auto rk = solve_scalar_equation(eqk, k_solution, k_sc);
+        const auto rw = solve_scalar_equation(eqw, omega_solution, omega_sc);
         for (std::size_t i = 0; i < n; ++i) {
             k(i) = k_solution(i);
             omega(i) = omega_solution(i);
@@ -130,13 +228,20 @@ inline TurbulenceTransportResult solve_sst_transport(
         for(std::size_t i=0;i<n;++i) {
             dk=std::max(dk,std::abs(k(i)-oldk(i)));
             dw=std::max(dw,std::abs(omega(i)-oldw(i)));
+            // Residuals must describe the state returned to the caller, after
+            // both linear relaxation and turbulence positivity enforcement.
+            k_solution(i)=k(i);
+            omega_solution(i)=omega(i);
         }
         result.k_residual=scalar_equation_residual_inf(eqk,k_solution);
         result.second_residual=scalar_equation_residual_inf(eqw,omega_solution);
         result.iterations=iter;
         if (result.nonfinite_bound_events==0 && rk.status==cfdx::core::SolverStatus::CONVERGED &&
            rw.status==cfdx::core::SolverStatus::CONVERGED &&
-           std::max(dk,dw)<=tolerance) {
+           std::isfinite(result.k_residual) && std::isfinite(result.second_residual) &&
+           std::max(dk,dw)<=tolerance &&
+           sst_equation_solution_defect(eqk,k_solution)<=tolerance &&
+           sst_equation_solution_defect(eqw,omega_solution)<=tolerance) {
             result.converged=true;
             break;
         }
@@ -157,15 +262,20 @@ inline TurbulenceTransportResult solve_sst_transport_dynamic_blending(
 {
     validate_turbulence_controls(controls);
     const std::size_t n=mesh.n_cells();
-    if(wall_distance.size()!=n||k.size()!=n||omega.size()!=n||strain_rate.size()!=n||mass_flux.size()!=mesh.n_faces())
+    if(n==0||wall_distance.size()!=n||k.size()!=n||omega.size()!=n||strain_rate.size()!=n||mass_flux.size()!=mesh.n_faces()||
+       wall_distance.dimension()!=1||k.dimension()!=1||omega.dimension()!=1||strain_rate.dimension()!=1||mass_flux.dimension()!=1)
         throw std::invalid_argument("SST dynamic-blending field size mismatch");
+    if (controls.model!=TurbulenceModel::SST || max_iterations==0 ||
+        !std::isfinite(tolerance) || tolerance<=0.0)
+        throw std::invalid_argument("SST dynamic blending requires SST and valid iteration controls");
     TurbulenceTransportResult result;
     for(std::size_t iter=1;iter<=max_iterations;++iter){
         auto old_k=k,old_w=omega;
         cfdx::core::Field<double,cfdx::core::Location::CELL> F1(n,"F1","1",1),F2(n,"F2","1",1);
         // Compute the full-cell gradient contraction once per outer iteration.
         // Recomputing this n-cell vector inside the n-cell loop made this O(n^2).
-        const auto cross=cell_gradient_dot(k,omega,mesh);
+        const auto cross=sst_boundary_gradient_dot(
+            k, omega, mesh, geometry, k_bcs, omega_bcs, omega_face_values);
         for(std::size_t i=0;i<n;++i){
             const auto b=compute_sst_blending(k(i),omega(i),wall_distance(i),controls.molecular_viscosity,controls.beta_star,cross[i],controls.density,controls.sst_sigma_w2);
             F1(i)=b.first; F2(i)=b.second;
diff --git a/src/cfdx/physics/steady_incompressible_solver.h b/src/cfdx/physics/steady_incompressible_solver.h
index 2bfa725d..46537766 100644
--- a/src/cfdx/physics/steady_incompressible_solver.h
+++ b/src/cfdx/physics/steady_incompressible_solver.h
@@ -15,6 +15,7 @@
 #include "cfdx/io/restart/dat_restart.h"
 #include "cfdx/physics/advanced_convergence.h"
 #include "cfdx/physics/finite_volume_transport.h"
+#include "cfdx/physics/momentum_viscous_stress.h"
 #include "cfdx/physics/pressure_velocity_algorithms.h"
 #include "cfdx/physics/pressure_velocity_system.h"
 #include "cfdx/physics/solver_control.h"
@@ -177,6 +178,10 @@ struct SSTProductionRuntime {
     // State must be explicitly initialized from case data or restored from a
     // DAT checkpoint with the helpers in sst_production_coupling.h. The
     // production solver advances this state transactionally when non-null.
+    // In this incompressible SST path the carrier pressure and pressure BCs
+    // represent p*=p_static+2*rho*k/3. The isotropic Reynolds stress is thus
+    // absorbed in the pressure gradient. Convert explicitly before reporting
+    // static pressure/Cp or prescribing static-pressure boundary data.
     SSTProductionState* state = nullptr;
     TurbulenceTransportControls transport;
     ScalarBoundaryConditions k_boundary_conditions;
@@ -186,6 +191,9 @@ struct SSTProductionRuntime {
     // remains valid for the duration of the update.
     ScalarBoundaryFaceValues omega_face_values;
     std::vector<std::string> wall_patch_names;
+    // Physical WALL_OMEGA_SST obligations from the case contract. This is
+    // distinct from the list of walls used to compute nearest wall distance.
+    std::vector<std::string> wall_omega_sst_patches;
     std::size_t max_transport_iterations = 100;
     double transport_tolerance = 1e-8;
 };
@@ -279,6 +287,15 @@ struct IncompressibleSolverControls {
         const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
         const cfdx::core::Field<double, cfdx::core::Location::CELL>&)>;
     InitializedStateCallback initialized_state_callback;
+    // Evidence hook at every actual momentum predictor (including retries).
+    // The viscosity and flux are the exact arrays consumed by assembly; U is
+    // the current carrier iterate. All references are read-only and transient.
+    using MomentumPredictorCallback = std::function<void(
+        std::size_t,
+        const std::vector<double>&,
+        const cfdx::core::Field<double, cfdx::core::Location::CELL>&,
+        const cfdx::core::Field<double, cfdx::core::Location::FACE>&)>;
+    MomentumPredictorCallback momentum_predictor_callback;
     // Evidence-only hook: invoked at the exact Krylov call point of the
     // coupled solve - after assembly, gauge treatment, finite checks and
     // preconditioner setup - so the snapshot matrix, RHS and initial guess
@@ -2519,8 +2536,16 @@ inline IncompressibleSolveResult solve_steady_incompressible(
     // Validate the caller-owned SST bundle before initialization can mutate
     // the carrier fields. Derived SST quantities must be refreshed only after
     // initialization/restart has established the actual velocity field.
-    if (controls.sst_runtime)
+    const FvGeometry geometry = build_fv_geometry(mesh);
+    if (controls.sst_runtime) {
         validate_sst_production_state(*controls.sst_runtime->state, mesh.n_cells());
+        const auto& runtime = *controls.sst_runtime;
+        const SSTRuntimeBoundaryContract boundary_contract{
+            runtime.k_boundary_conditions, runtime.omega_boundary_conditions,
+            runtime.wall_omega_sst_patches};
+        require_supported_sst_runtime_boundaries(
+            boundary_contract, mesh, geometry, &runtime.omega_face_values);
+    }
 
     apply_initialization(
         controls.initialization, mesh, U, p, restart_fields, restart_path);
@@ -2530,6 +2555,11 @@ inline IncompressibleSolveResult solve_steady_incompressible(
     // Authoritative dynamic viscosity used by predictor assembly, nonlinear
     // acceptance residuals, and the independent momentum conservation audit.
     std::vector<double> mu_eff(mesh.n_cells());
+    // The legacy uniform-viscosity incompressible path uses mu*Laplacian(U).
+    // Cellwise/SST viscosity requires the full symmetric stress; its transpose
+    // term is conservative and explicit in each nonlinear predictor.
+    const bool include_transpose_stress = controls.sst_runtime != nullptr ||
+        !controls.cell_turbulent_viscosity.empty();
     if (controls.sst_runtime) {
         std::map<std::string, std::array<double, 3>> fixed_velocity_values;
         for (const auto& [patch, bc] : velocity_bcs) {
@@ -2556,7 +2586,6 @@ inline IncompressibleSolveResult solve_steady_incompressible(
                 "effective cellwise dynamic viscosity must be finite and non-negative");
     }
 
-    const FvGeometry geometry = build_fv_geometry(mesh);
     cfdx::core::ConvergenceMonitorControls monitor_controls;
     monitor_controls.residual.absolute_tolerance = controls.convergence.absolute_tolerance;
     monitor_controls.residual.relative_tolerance = controls.convergence.relative_tolerance;
@@ -2881,6 +2910,21 @@ inline IncompressibleSolveResult solve_steady_incompressible(
         // Momentum predictor uses the conservative corrected flux from the
         // previous outer iteration. This preserves the SIMPLE phi state.
 
+        Field<double, Location::CELL> transpose_force;
+        if (include_transpose_stress) {
+            const auto stress_flux = transpose_viscous_stress_flux(
+                mesh, geometry, U, mu_eff, {ubc_x, ubc_y, ubc_z});
+            transpose_force = viscous_face_flux_divergence(mesh, geometry, stress_flux);
+            for (std::size_t cell=0; cell<mesh.n_cells(); ++cell) {
+                body_x(cell) += transpose_force(cell,0);
+                body_y(cell) += transpose_force(cell,1);
+                body_z(cell) += transpose_force(cell,2);
+            }
+        }
+
+        if (controls.momentum_predictor_callback)
+            controls.momentum_predictor_callback(iter, mu_eff, U, mass_flux);
+
         auto ex = assemble_momentum_component(
             mesh, geometry, mass_flux, grad_p, body_x, mu_eff, ubc_x, 0,
             controls.use_bounded_convection, controls.convection_scheme, &u_x_field);
@@ -3551,6 +3595,17 @@ inline IncompressibleSolveResult solve_steady_incompressible(
             final_uz_field(c) = U.component_data(2)[c];
         }
         auto final_grad_p = gauss_gradient_with_boundary(p, mesh, geometry, pressure_bcs);
+        if (include_transpose_stress) {
+            const auto final_stress_flux = transpose_viscous_stress_flux(
+                mesh, geometry, U, mu_eff, {ubc_x, ubc_y, ubc_z});
+            const auto final_transpose_force = viscous_face_flux_divergence(
+                mesh, geometry, final_stress_flux);
+            for (std::size_t cell=0; cell<mesh.n_cells(); ++cell) {
+                body_x(cell) += final_transpose_force(cell,0)-transpose_force(cell,0);
+                body_y(cell) += final_transpose_force(cell,1)-transpose_force(cell,1);
+                body_z(cell) += final_transpose_force(cell,2)-transpose_force(cell,2);
+            }
+        }
         auto final_ex = assemble_momentum_component(
             mesh, geometry, mass_flux, final_grad_p, body_x, mu_eff, ubc_x, 0,
             controls.use_bounded_convection, controls.convection_scheme, &final_ux_field);
@@ -4538,6 +4593,10 @@ inline IncompressibleSolveResult solve_steady_incompressible(
             // The flux is rebuilt from the accepted final field, independently
             // of the assembled matrix residual. Internal faces therefore have
             // one owner-oriented flux and are audited separately for balance.
+            Field<double, Location::FACE> audit_transpose_flux;
+            if (include_transpose_stress)
+                audit_transpose_flux = transpose_viscous_stress_flux(
+                    mesh, geometry, U, mu_eff, {ubc_x, ubc_y, ubc_z});
             for (std::size_t d = 0; d < 3; ++d) {
                 cfdx::core::Field<double, cfdx::core::Location::CELL> component(
                     mesh.n_cells(), "U_reconstructed", "m/s", 1);
@@ -4574,12 +4633,22 @@ inline IncompressibleSolveResult solve_steady_incompressible(
                                    : (d == 1
                                           ? controls.body_force.y
                                           : controls.body_force.z));
-                    source(cell) = body_component - final_grad_p.component_data(d)[cell];
+                    // audit_integrated_balance sums outward flux + source.
+                    // The physical equation is div(convective-viscous flux)
+                    // + grad(p*) - rho*f = 0.
+                    source(cell) = final_grad_p.component_data(d)[cell] -
+                                   controls.density * body_component;
                 }
-                const auto flux = reconstruct_scalar_transport_flux(
+                auto flux = reconstruct_scalar_transport_flux(
                     mesh, geometry, mass_flux, component, 0.0, scalar_bcs,
                     controls.use_bounded_convection, controls.convection_scheme,
                     nullptr, &mu_eff);
+                // Transport flux is convection minus viscous stress. Rebuild
+                // the transpose contribution from final U independently of
+                // the predictor source and assembled matrix residual.
+                if (include_transpose_stress)
+                    for (std::size_t face=0; face<mesh.n_faces(); ++face)
+                        flux(face) -= audit_transpose_flux(face,d);
                 const auto balance = cfdx::core::audit_integrated_balance(
                     mesh, flux, source, geometry.cell_volumes);
                 h.momentum_conservation_residual[d] = balance.residual;
diff --git a/src/cfdx/physics/thermal_zones.h b/src/cfdx/physics/thermal_zones.h
index cbb6f2fb..287176d7 100644
--- a/src/cfdx/physics/thermal_zones.h
+++ b/src/cfdx/physics/thermal_zones.h
@@ -429,7 +429,12 @@ struct ThermalConductionResult {
     bool converged = false;
     cfdx::core::SolverStatus status = cfdx::core::SolverStatus::MAX_ITER_REACHED;
     std::size_t iterations = 0;
+    // Status and candidate_converged describe the linear candidate;
+    // converged describes the actual returned temperature field.
     double residual_inf = 0.0;
+    double residual_relative_inf = 0.0;
+    bool candidate_converged = false;
+    bool returned_state_converged = false;
     double min_temperature = 0.0;             // K
     double max_temperature = 0.0;             // K
     double total_source_power = 0.0;         // W
@@ -461,9 +466,18 @@ inline ThermalConductionResult solve_thermal_conduction(
     ThermalConductionResult result;
     result.status = solve.status;
     result.iterations = solve.iterations;
+    result.candidate_converged = solve.status == SolverStatus::CONVERGED;
     result.residual_inf = scalar_equation_residual_inf(equation, temperature);
-    result.converged = (solve.status == SolverStatus::CONVERGED) &&
-                       std::isfinite(result.residual_inf);
+    double rhs_norm_inf = 0.0;
+    for (std::size_t i = 0; i < equation.rhs.size(); ++i)
+        rhs_norm_inf = std::max(rhs_norm_inf, std::abs(equation.rhs(i)));
+    result.residual_relative_inf =
+        result.residual_inf / std::max(rhs_norm_inf, 1e-300);
+    result.returned_state_converged =
+        std::isfinite(result.residual_inf) &&
+        std::isfinite(result.residual_relative_inf) &&
+        result.residual_relative_inf <= controls.tolerance;
+    result.converged = result.candidate_converged && result.returned_state_converged;
 
     for (std::size_t c = 0; c < nc; ++c) {
         if (c == 0) {
diff --git a/src/cfdx/physics/turbulence.h b/src/cfdx/physics/turbulence.h
index b9813e78..85b947c0 100644
--- a/src/cfdx/physics/turbulence.h
+++ b/src/cfdx/physics/turbulence.h
@@ -37,15 +37,13 @@ inline double turbulent_kinematic_viscosity_komega_sst(
     return a1 * k / std::max(a1 * omega, strain_rate * F2);
 }
 
-// Backward-compatible overload retained for existing callers using the
-// historical algebraic four-argument form.
+// Method contract: docs/audits/09-turbulence-rans.md
+// The historical overload omitted strain_rate and silently substituted
+// omega in its place. That is only correct for the special case S == omega.
+// Reject it at compile time so callers must provide the physical invariant.
+[[deprecated("Use the five-argument SST viscosity overload and provide strain_rate")]]
 inline double turbulent_kinematic_viscosity_komega_sst(
-    double k, double omega, double F2, double a1 = 0.31)
-{
-    if (k < 0.0 || omega <= 0.0 || F2 < 0.0 || F2 > 1.0 || a1 <= 0.0)
-        throw std::invalid_argument("k-omega SST: invalid k, omega, F2 or a1");
-    return a1 * k / std::max(a1 * omega, omega * F2);
-}
+    double k, double omega, double F2, double a1 = 0.31) = delete;
 
 inline double strain_rate_magnitude(const cfdx::core::Vec3& s)
 {
diff --git a/src/cfdx/python/cfdx/io/adapters/fluent.py b/src/cfdx/python/cfdx/io/adapters/fluent.py
index c68817ed..4e3a6498 100644
--- a/src/cfdx/python/cfdx/io/adapters/fluent.py
+++ b/src/cfdx/python/cfdx/io/adapters/fluent.py
@@ -1030,9 +1030,9 @@ def _extract_cas_h5_reference_values(rp_vars: dict, setup: CaseSetup) -> None:
         gx = float(rp_vars.get("gravity/x", 0.0))
         gy = float(rp_vars.get("gravity/y", 0.0))
         gz = float(rp_vars.get("gravity/z", 0.0))
-        setup.gravity_vector = f"[{gx}, {gy}, {gz}]"
+        setup.gravity_vector = f"{gx:.17g} {gy:.17g} {gz:.17g}"
     else:
-        setup.gravity_vector = "[0, 0, 0]"
+        setup.gravity_vector = "0 0 0"
 
     setup.energy_model = "energy" if rp_vars.get("energy?", False) else "isothermal"
 
diff --git a/src/cfdx/python/cfdx/io/schema.py b/src/cfdx/python/cfdx/io/schema.py
index da91922d..89b16e9d 100644
--- a/src/cfdx/python/cfdx/io/schema.py
+++ b/src/cfdx/python/cfdx/io/schema.py
@@ -59,12 +59,14 @@ class SourceInfo(BaseModel):
 
 class MaterialSpec(BaseModel):
     name: str = "air"
-    density: float = 1.225
-    dynamic_viscosity: float = 1.789e-5
-    thermal_conductivity: float = 0.024
-    specific_heat: float = 1006.43
-    molecular_weight: float = 28.97
-    eos_model: str = "ideal_gas"
+    # Match the native typed schema: zero/empty means unspecified. The selected
+    # physics must validate its required properties instead of inventing air.
+    density: float = 0.0
+    dynamic_viscosity: float = 0.0
+    thermal_conductivity: float = 0.0
+    specific_heat: float = 0.0
+    molecular_weight: float = 0.0
+    eos_model: str = ""
     turbulence_model: str = "laminar"
     extra_properties: dict[str, str] = Field(default_factory=dict)
 
diff --git a/tests/audit/chain_fixture.h b/tests/audit/chain_fixture.h
new file mode 100644
index 00000000..7d348b32
--- /dev/null
+++ b/tests/audit/chain_fixture.h
@@ -0,0 +1,73 @@
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
diff --git a/tests/audit/sst_channel_fixture.h b/tests/audit/sst_channel_fixture.h
new file mode 100644
index 00000000..dba8ddb0
--- /dev/null
+++ b/tests/audit/sst_channel_fixture.h
@@ -0,0 +1,60 @@
+#pragma once
+
+#include "chain_fixture.h"
+#include "cfdx/physics/steady_incompressible_solver.h"
+
+using namespace cfdx::physics;
+
+// Three-cell production lifecycle fixture with an accepted SST transport step.
+// It is not a physically validated turbulent channel benchmark.
+struct SSTChannelCase {
+    Mesh mesh = make_1d_chain(3);
+    Field<double,Location::CELL> U{3,"U","m/s",3};
+    Field<double,Location::CELL> p{3,"p","Pa",1};
+    SSTProductionState state{3};
+    SSTProductionRuntime runtime;
+    IncompressibleSolverControls controls;
+    VelocityBoundaryConditions ubc;
+    ScalarBoundaryConditions pbc;
+
+    SSTChannelCase() {
+        mesh.boundary().patch(0).type = PatchType::INLET;
+        mesh.boundary().patch(1).type = PatchType::OUTLET;
+        U.fill(-7.0); p.fill(4.0);
+        state.k.fill(0.01); state.omega.fill(2.0); state.nu_t.fill(9.0);
+        state.wall_distance.fill(17.0); state.strain_rate.fill(11.0);
+        controls.density = 1.2;
+        controls.kinematic_viscosity = 0.01;
+        controls.convergence.max_iterations = 1;
+        controls.linear_tolerance = 1e-12;
+        controls.convergence.continuity_tolerance = 1e-12;
+        controls.initialization.mode = InitializationMode::Uniform;
+        controls.initialization.uniform_velocity = {1,0,0};
+        controls.initialization.uniform_pressure = 0.0;
+        ubc["inlet"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}};
+        ubc["outlet"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT,{0,0,0}};
+        // Moving walls make the carrier a uniform exact solution; this is a
+        // production lifecycle test, not a channel-flow validation benchmark.
+        ubc["wall"] = {VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}};
+        pbc["inlet"] = {ScalarBoundaryType::ZERO_GRADIENT,0,0};
+        pbc["outlet"] = {ScalarBoundaryType::FIXED_VALUE,0,0};
+        pbc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT,0,0};
+        runtime.state = &state;
+        runtime.transport.model = TurbulenceModel::SST;
+        runtime.transport.density = controls.density;
+        runtime.transport.molecular_viscosity = controls.kinematic_viscosity;
+        runtime.wall_patch_names = {"wall"};
+        runtime.wall_omega_sst_patches = {"wall"};
+        runtime.k_boundary_conditions = {
+            {"inlet",{ScalarBoundaryType::FIXED_VALUE,0.01,0}},
+            {"outlet",{ScalarBoundaryType::ZERO_GRADIENT,0,0}},
+            {"wall",{ScalarBoundaryType::FIXED_VALUE,0,0}}};
+        runtime.omega_boundary_conditions = {
+            {"inlet",{ScalarBoundaryType::FIXED_VALUE,2,0}},
+            {"outlet",{ScalarBoundaryType::ZERO_GRADIENT,0,0}}};
+        runtime.omega_face_values = make_sst_wall_omega_face_values(
+            mesh, build_fv_geometry(mesh), runtime.wall_omega_sst_patches, runtime.transport);
+        controls.sst_runtime = &runtime;
+    }
+    SSTChannelCase(const SSTChannelCase&) = delete;
+};
diff --git a/tests/audit/test_audit_cht_time.cpp b/tests/audit/test_audit_cht_time.cpp
new file mode 100644
index 00000000..06eb713a
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
diff --git a/tests/audit/test_audit_energy_boundary_metric.cpp b/tests/audit/test_audit_energy_boundary_metric.cpp
new file mode 100644
index 00000000..42580ac5
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
diff --git a/tests/audit/test_sst_integration_gaps.cpp b/tests/audit/test_sst_integration_gaps.cpp
new file mode 100644
index 00000000..887a2527
--- /dev/null
+++ b/tests/audit/test_sst_integration_gaps.cpp
@@ -0,0 +1,382 @@
+#include "sst_channel_fixture.h"
+#include "cfdx/physics/sst_production_coupling.h"
+#include "cfdx/physics/steady_incompressible_solver.h"
+#include "common/test_harness.h"
+#include <limits>
+
+using namespace cfdx::physics;
+using namespace cfdx::testing;
+
+static Mesh make_shear_chain(std::size_t n)
+{
+    auto mesh=make_1d_chain(n);
+    const auto faces=mesh.boundary().patch(2).face_ids;
+    mesh.boundary().remove_patch(2);
+    const std::array<std::string,4> names{"y0","y1","z0","z1"};
+    for (std::size_t direction=0; direction<4; ++direction) {
+        Patch patch; patch.name=names[direction]; patch.type=PatchType::WALL;
+        for (std::size_t cell=0; cell<n; ++cell)
+            patch.face_ids.push_back(faces[direction*n+cell]);
+        mesh.boundary().add_patch(patch);
+    }
+    return mesh;
+}
+
+static std::array<ScalarBoundaryConditions,3> shear_scalar_bcs(const Mesh& mesh)
+{
+    std::array<ScalarBoundaryConditions,3> bcs;
+    for (std::size_t component=0; component<3; ++component)
+        for (const auto& patch : mesh.boundary().patches()) {
+            const double value=component==0 ? (patch.name=="y0" ? 0.0 :
+                                               (patch.name=="y1" ? 1.0 : 0.5)) : 0.0;
+            bcs[component][patch.name]={ScalarBoundaryType::FIXED_VALUE,value,0.0};
+        }
+    return bcs;
+}
+
+int main()
+{
+    run_case("sst_wall_omega_requires_every_mesh_face_before_carrier_mutation", [] {
+        SSTChannelCase c;
+        const auto geometry = build_fv_geometry(c.mesh);
+        const SSTRuntimeBoundaryContract contract{
+            c.runtime.k_boundary_conditions, c.runtime.omega_boundary_conditions, {"wall"}};
+        require_supported_sst_runtime_boundaries(contract, c.mesh, geometry, &c.runtime.omega_face_values);
+        const auto wall_faces = c.mesh.boundary().patch(2).face_ids;
+        EXPECT_TRUE(wall_faces.size() == 12);
+        const auto good = c.runtime.omega_face_values;
+        for (const double invalid : {std::numeric_limits<double>::quiet_NaN(),
+                                     std::numeric_limits<double>::infinity(), 0.0, -1.0}) {
+            c.runtime.omega_face_values = good;
+            c.runtime.omega_face_values.values.at("wall")[wall_faces.back()] = invalid;
+            EXPECT_THROW_WITH(solve_steady_incompressible(c.mesh,c.U,c.p,c.ubc,c.pbc,c.controls),
+                              std::invalid_argument, "missing or invalid value");
+            for (std::size_t cell=0; cell<c.mesh.n_cells(); ++cell) {
+                for (std::size_t d=0; d<3; ++d) EXPECT_NEAR(c.U(cell,d),-7.0,0.0);
+                EXPECT_NEAR(c.p(cell),4.0,0.0);
+                EXPECT_NEAR(c.state.k(cell),0.01,0.0);
+                EXPECT_NEAR(c.state.omega(cell),2.0,0.0);
+                EXPECT_NEAR(c.state.nu_t(cell),9.0,0.0);
+                EXPECT_NEAR(c.state.wall_distance(cell),17.0,0.0);
+                EXPECT_NEAR(c.state.strain_rate(cell),11.0,0.0);
+            }
+        }
+        auto short_values = good;
+        short_values.values.at("wall").resize(1);
+        EXPECT_THROW(require_supported_sst_runtime_boundaries(contract,c.mesh,geometry,&short_values),
+                     std::invalid_argument);
+        auto wrong_geometry = geometry;
+        wrong_geometry.face_patch[wall_faces.front()] = 0;
+        EXPECT_THROW(require_supported_sst_runtime_boundaries(contract,c.mesh,wrong_geometry,&good),
+                     std::invalid_argument);
+        auto flat_geometry = geometry;
+        const auto face = wall_faces.front();
+        flat_geometry.face_centres[face] = geometry.cell_centres[c.mesh.ownership().owner(face)];
+        EXPECT_THROW_WITH(require_supported_sst_runtime_boundaries(contract,c.mesh,flat_geometry,&good),
+                          std::invalid_argument,"positive geometric normal distance");
+    });
+
+    run_case("sst_cross_diffusion_uses_real_scalar_and_face_boundaries", [] {
+        const auto mesh = make_1d_chain(3);
+        const auto geometry = build_fv_geometry(mesh);
+        Field<double,Location::CELL> k(3,"k","",1), omega(3,"omega","",1);
+        for (std::size_t cell=0; cell<3; ++cell) {
+            const auto x = geometry.cell_centres[cell];
+            k(cell)=2+x.x;
+            omega(cell)=3+2*x.x+3*x.y;
+        }
+        ScalarBoundaryConditions kbc{
+            {"inlet",{ScalarBoundaryType::FIXED_VALUE,2,0}},
+            {"outlet",{ScalarBoundaryType::FIXED_VALUE,3,0}}};
+        ScalarBoundaryConditions wbc{
+            {"inlet",{ScalarBoundaryType::FIXED_VALUE,4.5,0}},
+            {"outlet",{ScalarBoundaryType::FIXED_VALUE,6.5,0}},
+            {"wall",{ScalarBoundaryType::FIXED_VALUE,999,0}}};
+        ScalarBoundaryFaceValues faces;
+        faces.values["wall"].assign(mesh.n_faces(),std::numeric_limits<double>::quiet_NaN());
+        for (const auto face : mesh.boundary().patch(2).face_ids) {
+            const auto x=geometry.face_centres[face];
+            faces.values["wall"][face]=3+2*x.x+3*x.y;
+        }
+        const auto contraction=sst_boundary_gradient_dot(k,omega,mesh,geometry,kbc,wbc,&faces);
+        const auto gradw=scalar_gradient_with_transport_boundaries(omega,mesh,geometry,wbc,&faces);
+        for (std::size_t cell=0; cell<3; ++cell) {
+            EXPECT_NEAR(contraction[cell],2.0,2e-14);
+            EXPECT_NEAR(gradw(cell,0),2.0,2e-14);
+            EXPECT_NEAR(gradw(cell,1),3.0,2e-14);
+            EXPECT_NEAR(gradw(cell,2),0.0,2e-14);
+        }
+        const auto old_contraction=cell_gradient_dot(k,omega,mesh);
+        EXPECT_NEAR(old_contraction.front(),0.5,2e-14);
+        EXPECT_NEAR(old_contraction.back(),0.5,2e-14);
+        kbc["inlet"]={ScalarBoundaryType::FIXED_GRADIENT,0,-1};
+        kbc["outlet"]={ScalarBoundaryType::FIXED_GRADIENT,0,1};
+        const auto gradk=scalar_gradient_with_transport_boundaries(k,mesh,geometry,kbc);
+        for (std::size_t cell=0; cell<3; ++cell) EXPECT_NEAR(gradk(cell,0),1.0,2e-14);
+    });
+
+    run_case("sst_residuals_and_convergence_describe_bounded_returned_fields", [] {
+        const auto mesh=make_1d_chain(1);
+        const auto geometry=build_fv_geometry(mesh);
+        Field<double,Location::FACE> flux(mesh.n_faces(),"phi","",1); flux.fill(0);
+        Field<double,Location::CELL> k(1,"k","",1),w(1,"omega","",1),S(1,"S","",1),F(1,"F","",1);
+        k.fill(1); w.fill(1); S.fill(0); F.fill(1);
+        TurbulenceTransportControls controls;
+        controls.model=TurbulenceModel::SST;
+        controls.k_min=1; controls.omega_min=1;
+        const auto result=solve_sst_transport(mesh,geometry,flux,k,w,S,F,F,controls,{}, {},3,1e-8);
+        EXPECT_FALSE(result.converged);
+        EXPECT_NEAR(k(0),1.0,0.0); EXPECT_NEAR(w(0),1.0,0.0);
+        EXPECT_NEAR(result.k_residual,controls.density*controls.beta_star,1e-14);
+        EXPECT_NEAR(result.second_residual,controls.density*controls.beta1,1e-14);
+        EXPECT_TRUE(result.k_floor_events==3 && result.second_floor_events==3);
+        std::cout << "    bounded-state residuals: k=" << result.k_residual
+                  << " omega=" << result.second_residual << " converged=" << result.converged << '\n';
+        ScalarEquation finite_equation{SparseMatrix(1,1),Vector(1,2.0),{1.0},0.0};
+        finite_equation.matrix.push_back(0,0,1.0);
+        finite_equation.matrix.finalize();
+        EXPECT_NEAR(sst_linear_tolerance_for_solution_defect(finite_equation,1e-8),5e-9,1e-24);
+        EXPECT_THROW_WITH(sst_linear_tolerance_for_solution_defect(
+                              finite_equation,std::numeric_limits<double>::denorm_min()),
+                          std::domain_error,"cannot be represented");
+        finite_equation.rhs(0)=std::numeric_limits<double>::infinity();
+        EXPECT_THROW(sst_linear_tolerance_for_solution_defect(finite_equation,1e-8),std::invalid_argument);
+        ScalarEquation zero_row{SparseMatrix(1,1),Vector(1,1.0),{0.0},0.0};
+        zero_row.matrix.push_back(0,0,0.0);
+        zero_row.matrix.finalize();
+        EXPECT_THROW_WITH(sst_linear_tolerance_for_solution_defect(zero_row,1e-8),
+                          std::invalid_argument,"zero or non-finite matrix row norm");
+    });
+
+    run_case("sst_multicell_first_predictor_uses_refreshed_viscosity_and_nonzero_flux", [] {
+        SSTChannelCase c;
+        std::size_t calls=0;
+        c.controls.momentum_predictor_callback = [&](std::size_t iteration,const auto& mu,const auto& U,const auto& flux) {
+            ++calls;
+            EXPECT_TRUE(iteration==1);
+            for (std::size_t cell=0; cell<c.mesh.n_cells(); ++cell) {
+                EXPECT_NEAR(U(cell,0),1.0,1e-14);
+                EXPECT_NEAR(c.state.strain_rate(cell),0.0,1e-14);
+                EXPECT_NEAR(c.state.wall_distance(cell),0.5,1e-14);
+                EXPECT_NEAR(mu[cell],1.2*(0.01+0.01/2.0),1e-14);
+            }
+            EXPECT_NEAR(flux(0),-1.2,1e-14);
+            EXPECT_NEAR(flux(1),1.2,1e-14);
+            for (std::size_t face=0; face<c.mesh.n_faces(); ++face)
+                if (c.mesh.ownership().neighbour(face)>=0) EXPECT_NEAR(flux(face),1.2,1e-14);
+        };
+        const auto result=solve_steady_incompressible(c.mesh,c.U,c.p,c.ubc,c.pbc,c.controls);
+        EXPECT_TRUE(calls==1 && result.history.size()==1);
+        EXPECT_TRUE(result.authoritative_mass_flux.size()==c.mesh.n_faces());
+        EXPECT_NEAR(result.authoritative_mass_flux(0),-1.2,1e-12);
+        EXPECT_NEAR(result.authoritative_mass_flux(1),1.2,1e-12);
+        EXPECT_TRUE(c.state.k(0)>c.state.k(2));
+        EXPECT_TRUE(c.state.k(0)<0.01);
+        std::cout << "    nonzero-flux SST: k=[" << c.state.k(0) << ',' << c.state.k(1) << ',' << c.state.k(2)
+                  << "] omega=[" << c.state.omega(0) << ',' << c.state.omega(1) << ',' << c.state.omega(2) << "]\n";
+    });
+
+    run_case("sst_multicell_retry_restores_every_field_and_authoritative_flux", [] {
+        SSTChannelCase c;
+        c.controls.nonlinear_retry.max_retries=1;
+        std::size_t predictor_calls=0,output_calls=0;
+        Field<double,Location::CELL> first_U;
+        Field<double,Location::FACE> first_flux;
+        SSTProductionState first_state(3);
+        std::vector<double> first_mu;
+        c.controls.momentum_predictor_callback = [&](std::size_t,const auto& mu,const auto& U,const auto& flux) {
+            ++predictor_calls;
+            if (predictor_calls==1) { first_U=U; first_flux=flux; first_state=c.state; first_mu=mu; }
+            else {
+                for (std::size_t cell=0; cell<3; ++cell) {
+                    for (std::size_t d=0; d<3; ++d) EXPECT_NEAR(U(cell,d),first_U(cell,d),0.0);
+                    EXPECT_NEAR(c.state.k(cell),first_state.k(cell),0.0);
+                    EXPECT_NEAR(c.state.omega(cell),first_state.omega(cell),0.0);
+                    EXPECT_NEAR(c.state.nu_t(cell),first_state.nu_t(cell),0.0);
+                    EXPECT_NEAR(c.state.wall_distance(cell),first_state.wall_distance(cell),0.0);
+                    EXPECT_NEAR(c.state.strain_rate(cell),first_state.strain_rate(cell),0.0);
+                    EXPECT_NEAR(mu[cell],first_mu[cell],0.0);
+                }
+                for (std::size_t face=0; face<flux.size(); ++face) EXPECT_NEAR(flux(face),first_flux(face),0.0);
+            }
+        };
+        c.controls.iteration_output_callback = [&](std::size_t,double,const auto&,const auto&,const auto&) {
+            if (++output_calls==1) {
+                EXPECT_TRUE(c.state.k(0)!=first_state.k(0));
+                throw NonlinearRetryableFailure("test rejection after SST commit with nonzero flux");
+            }
+            return false;
+        };
+        const auto result=solve_steady_incompressible(c.mesh,c.U,c.p,c.ubc,c.pbc,c.controls);
+        EXPECT_TRUE(predictor_calls==2 && output_calls==2 && result.history.size()==1);
+        SSTChannelCase reference;
+        reference.controls.coupling.alpha_u=result.history[0].effective_alpha_u;
+        reference.controls.coupling.alpha_p=result.history[0].effective_alpha_p;
+        const auto clean=solve_steady_incompressible(reference.mesh,reference.U,reference.p,
+                                                   reference.ubc,reference.pbc,reference.controls);
+        for (std::size_t cell=0; cell<3; ++cell) {
+            for (std::size_t d=0; d<3; ++d) EXPECT_NEAR(c.U(cell,d),reference.U(cell,d),1e-12);
+            EXPECT_NEAR(c.p(cell),reference.p(cell),1e-12);
+            EXPECT_NEAR(c.state.k(cell),reference.state.k(cell),1e-12);
+            EXPECT_NEAR(c.state.omega(cell),reference.state.omega(cell),1e-12);
+            EXPECT_NEAR(c.state.nu_t(cell),reference.state.nu_t(cell),1e-12);
+        }
+        for (std::size_t face=0; face<c.mesh.n_faces(); ++face)
+            EXPECT_NEAR(result.authoritative_mass_flux(face),clean.authoritative_mass_flux(face),1e-12);
+    });
+
+    run_case("sst_multicell_dat_rehydrates_primary_fields_before_first_predictor", [] {
+        SSTChannelCase saved;
+        const auto accepted=solve_steady_incompressible(
+            saved.mesh,saved.U,saved.p,saved.ubc,saved.pbc,saved.controls);
+        EXPECT_NEAR(accepted.authoritative_mass_flux(0),-1.2,1e-12);
+        const auto path=unique_temp_path("_sst_multicell_restore.dat");
+        cfdx::io::DatRestartFields fields;
+        fields.k=&saved.state.k;
+        fields.omega=&saved.state.omega;
+        cfdx::io::write_dat_restart_fields(path.string(),saved.mesh,saved.U,saved.p,fields,1,0.0);
+        SSTChannelCase restored;
+        const std::map<std::string,std::array<double,3>> fixed_velocity{
+            {"inlet",{1,0,0}},{"wall",{1,0,0}}};
+        const auto metadata=restore_sst_production_state_from_dat(
+            path.string(),restored.mesh,restored.U,restored.p,restored.state,
+            {0.5,0.5,0.5},restored.runtime.transport,fixed_velocity);
+        EXPECT_TRUE(metadata.iteration==1);
+        for (std::size_t cell=0; cell<3; ++cell) {
+            for (std::size_t d=0; d<3; ++d) EXPECT_NEAR(restored.U(cell,d),saved.U(cell,d),0.0);
+            EXPECT_NEAR(restored.p(cell),saved.p(cell),0.0);
+            EXPECT_NEAR(restored.state.k(cell),saved.state.k(cell),0.0);
+            EXPECT_NEAR(restored.state.omega(cell),saved.state.omega(cell),0.0);
+            EXPECT_NEAR(restored.state.nu_t(cell),saved.state.nu_t(cell),1e-15);
+        }
+        restored.controls.initialization.mode=InitializationMode::Provided;
+        std::size_t predictor_calls=0;
+        restored.controls.momentum_predictor_callback=[&](std::size_t,const auto& mu,const auto& U,const auto&) {
+            ++predictor_calls;
+            for (std::size_t cell=0; cell<3; ++cell) {
+                EXPECT_NEAR(mu[cell],1.2*(0.01+saved.state.nu_t(cell)),1e-15);
+                for (std::size_t d=0; d<3; ++d) EXPECT_NEAR(U(cell,d),saved.U(cell,d),0.0);
+            }
+        };
+        solve_steady_incompressible(restored.mesh,restored.U,restored.p,
+                                   restored.ubc,restored.pbc,restored.controls);
+        EXPECT_TRUE(predictor_calls==1);
+        std::filesystem::remove(path);
+        // DAT v4 does not store authoritative face flux. This test certifies
+        // primary-field hydration and derived viscosity only, not equivalence
+        // of a general interrupted coupled trajectory.
+    });
+
+    run_case("variable_viscosity_transpose_stress_mms_converges_and_is_conservative", [] {
+        // Divergence-free U=(y,0,0), mu=1+x^2(1-x)^2. The scalar diffusion
+        // alone gives zero, whereas div(mu*grad(U)^T)=(0,mu'(x),0).
+        // mu'(0)=mu'(1)=0 makes the declared owner boundary-mu convention
+        // consistent to second order on this orthogonal family.
+        double previous_error=0.0;
+        for (const std::size_t n : {16u,32u,64u}) {
+            const auto mesh=make_shear_chain(n);
+            const auto geometry=build_fv_geometry(mesh);
+            Field<double,Location::CELL> U(n,"U","m/s",3); U.fill(0.0);
+            std::vector<double> mu(n);
+            for (std::size_t cell=0; cell<n; ++cell) {
+                U(cell,0)=geometry.cell_centres[cell].y;
+                const double x=geometry.cell_centres[cell].x;
+                mu[cell]=1+x*x*(1-x)*(1-x);
+            }
+            const auto flux=transpose_viscous_stress_flux(mesh,geometry,U,mu,shear_scalar_bcs(mesh));
+            const auto force=viscous_face_flux_divergence(mesh,geometry,flux);
+            double error_sq=0.0,omitted_sq=0.0,volume_sum=0.0,boundary_sum=0.0;
+            for (std::size_t cell=0; cell<n; ++cell) {
+                const double x=geometry.cell_centres[cell].x;
+                const double exact=2*x*(1-x)*(1-2*x);
+                error_sq+=std::pow(force(cell,1)-exact,2)*geometry.cell_volumes[cell];
+                omitted_sq+=exact*exact*geometry.cell_volumes[cell];
+                volume_sum+=force(cell,1)*geometry.cell_volumes[cell];
+                EXPECT_NEAR(force(cell,0),0.0,2e-12);
+                EXPECT_NEAR(force(cell,2),0.0,2e-12);
+            }
+            for (std::size_t face=0; face<mesh.n_faces(); ++face)
+                if (mesh.ownership().neighbour(face)<0) boundary_sum+=flux(face,1);
+            EXPECT_NEAR(volume_sum,boundary_sum,2e-13);
+            const double error=std::sqrt(error_sq);
+            EXPECT_TRUE(error<std::sqrt(omitted_sq)/20.0);
+            double order=0.0;
+            if (previous_error>0.0) {
+                order=std::log(previous_error/error)/std::log(2.0);
+                EXPECT_TRUE(order>1.8);
+            }
+            std::cout << "    transpose-stress MMS N=" << n << " L2=" << error
+                      << " omitted=" << std::sqrt(omitted_sq) << " order=" << order << '\n';
+            previous_error=error;
+            std::fill(mu.begin(),mu.end(),1.0);
+            const auto constant_flux=transpose_viscous_stress_flux(mesh,geometry,U,mu,shear_scalar_bcs(mesh));
+            const auto constant_force=viscous_face_flux_divergence(mesh,geometry,constant_flux);
+            for (std::size_t cell=0; cell<n; ++cell)
+                for (std::size_t component=0; component<3; ++component)
+                    EXPECT_NEAR(constant_force(cell,component),0.0,2e-12);
+        }
+    });
+
+    run_case("variable_viscosity_production_predictor_consumes_transpose_stress", [] {
+        const auto mesh=make_shear_chain(8);
+        const auto geometry=build_fv_geometry(mesh);
+        Field<double,Location::CELL> U(8,"U","m/s",3),p(8,"p","Pa",1);
+        U.fill(0.0); p.fill(0.0);
+        VelocityBoundaryConditions ubc;
+        ScalarBoundaryConditions pbc;
+        const auto scalar_bcs=shear_scalar_bcs(mesh);
+        for (const auto& patch : mesh.boundary().patches()) {
+            ubc[patch.name]={VelocityBoundaryCondition::Type::FIXED_VALUE,
+                            {scalar_bcs[0].at(patch.name).value,0,0}};
+            pbc[patch.name]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
+        }
+        IncompressibleSolverControls controls;
+        controls.kinematic_viscosity=1.0;
+        controls.convergence.max_iterations=1;
+        controls.linear_tolerance=1e-12;
+        controls.cell_turbulent_viscosity.resize(8);
+        for (std::size_t cell=0; cell<8; ++cell) {
+            const double x=geometry.cell_centres[cell].x;
+            U(cell,0)=0.5;
+            controls.cell_turbulent_viscosity[cell]=x*x*(1-x)*(1-x);
+        }
+        const auto result=solve_steady_incompressible(mesh,U,p,ubc,pbc,controls);
+        EXPECT_TRUE(result.history.size()==1);
+        // With the old scalar-only operator Uy stays exactly zero. The new
+        // transverse motion must have the sign of mu'(x), not its opposite.
+        EXPECT_TRUE(U(0,1)>1e-5);
+        EXPECT_TRUE(U(7,1)<-1e-5);
+        std::cout << "    transverse predictor Uy=" << U(0,1) << ',' << U(7,1) << '\n';
+    });
+
+    run_case("momentum_balance_uses_physical_source_sign_and_density_in_laminar_limit", [] {
+        const auto mesh=make_1d_chain(1);
+        Field<double,Location::CELL> U(1,"U","m/s",3),p(1,"p","Pa",1);
+        U.fill(0.0); p.fill(0.0);
+        VelocityBoundaryConditions ubc;
+        ScalarBoundaryConditions pbc;
+        for (const auto& patch : mesh.boundary().patches()) {
+            ubc[patch.name]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
+            pbc[patch.name]={ScalarBoundaryType::ZERO_GRADIENT,0,0};
+        }
+        IncompressibleSolverControls controls;
+        controls.density=2.5;
+        controls.kinematic_viscosity=0.4;
+        controls.cell_turbulent_viscosity={0.0};
+        controls.body_force={3,0,0};
+        controls.convergence.max_iterations=1;
+        controls.coupling.alpha_u=0.7;
+        const auto result=solve_steady_incompressible(mesh,U,p,ubc,pbc,controls);
+        // Six unit faces at distance 1/2 give conductance 12*mu=12.
+        // Outward diffusive flux is 12*Ux, source integral is -rho*f=-7.5.
+        EXPECT_NEAR(U(0,0),0.7*7.5/12.0,1e-12);
+        EXPECT_NEAR(result.history.at(0).momentum_conservation_residual[0],-2.25,1e-12);
+        const auto cellwise_U=U;
+        U.fill(0.0); p.fill(0.0);
+        controls.cell_turbulent_viscosity.clear();
+        const auto uniform=solve_steady_incompressible(mesh,U,p,ubc,pbc,controls);
+        EXPECT_NEAR(U(0,0),cellwise_U(0,0),1e-12);
+        EXPECT_NEAR(uniform.history.at(0).momentum_conservation_residual[0],-2.25,1e-12);
+    });
+    return run_all();
+}
diff --git a/tests/common/test_harness.h b/tests/common/test_harness.h
index fadc0bfa..80e1b92c 100644
--- a/tests/common/test_harness.h
+++ b/tests/common/test_harness.h
@@ -55,8 +55,16 @@ void run_case(const std::string& name, std::function<void()> func) {
 
 #define EXPECT_NEAR(a, b, eps) \
     do { \
-        if (std::abs((a) - (b)) > (eps)) { \
-            std::cerr << "  FAILED: " << #a " ≈ " #b " (diff=" << std::abs((a) - (b)) << " > " << eps << ")" << " in " << __FUNCTION__ << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
+        const double cfdx_expect_near_a_ = static_cast<double>(a); \
+        const double cfdx_expect_near_b_ = static_cast<double>(b); \
+        const double cfdx_expect_near_eps_ = static_cast<double>(eps); \
+        const double cfdx_expect_near_diff_ = std::abs(cfdx_expect_near_a_ - cfdx_expect_near_b_); \
+        if (!std::isfinite(cfdx_expect_near_a_) || !std::isfinite(cfdx_expect_near_b_) || \
+            !std::isfinite(cfdx_expect_near_eps_) || cfdx_expect_near_eps_ < 0.0 || \
+            !std::isfinite(cfdx_expect_near_diff_) || cfdx_expect_near_diff_ > cfdx_expect_near_eps_) { \
+            std::cerr << "  FAILED: " << #a " ≈ " #b " (diff=" << cfdx_expect_near_diff_ \
+                      << ", tolerance=" << cfdx_expect_near_eps_ << ")" << " in " << __FUNCTION__ \
+                      << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
             throw std::runtime_error("Assertion failed: " #a " ≈ " #b); \
         } \
     } while (0)
diff --git a/tests/python/test_imported_thermal_contracts.py b/tests/python/test_imported_thermal_contracts.py
new file mode 100644
index 00000000..199d03be
--- /dev/null
+++ b/tests/python/test_imported_thermal_contracts.py
@@ -0,0 +1,36 @@
+"""Canonical import data must preserve unknown properties and parseable gravity."""
+from __future__ import annotations
+
+import unittest
+
+from cfdx.io.schema import CaseSetup, MaterialSpec
+from cfdx.io.adapters.fluent import _extract_cas_h5_reference_values
+
+
+class ImportedThermalContracts(unittest.TestCase):
+    def test_unspecified_material_does_not_become_air(self):
+        material = MaterialSpec(name="unknown")
+        for name in ("density", "dynamic_viscosity", "thermal_conductivity",
+                     "specific_heat", "molecular_weight"):
+            self.assertEqual(getattr(material, name), 0.0)
+        self.assertEqual(material.eos_model, "")
+        self.assertEqual(MaterialSpec.model_validate_json(material.model_dump_json()), material)
+
+    def test_explicit_material_values_survive_serialization(self):
+        material = MaterialSpec(name="declared", density=1.5, dynamic_viscosity=0.002,
+                                thermal_conductivity=2.0, specific_heat=900.0)
+        loaded = MaterialSpec.model_validate_json(material.model_dump_json())
+        self.assertEqual(loaded, material)
+
+    def test_fluent_gravity_uses_native_three_scalar_contract(self):
+        setup = CaseSetup()
+        _extract_cas_h5_reference_values(
+            {"gravity?": True, "gravity/x": 1.5, "gravity/y": -9.81, "gravity/z": 0.25}, setup)
+        self.assertEqual([float(value) for value in setup.gravity_vector.split()],
+                         [1.5, -9.81, 0.25])
+        _extract_cas_h5_reference_values({"gravity?": False}, setup)
+        self.assertEqual([float(value) for value in setup.gravity_vector.split()], [0.0, 0.0, 0.0])
+
+
+if __name__ == "__main__":
+    unittest.main(verbosity=2)
diff --git a/tests/python/test_production_thermal_restart_e2e.py b/tests/python/test_production_thermal_restart_e2e.py
new file mode 100644
index 00000000..9b21d700
--- /dev/null
+++ b/tests/python/test_production_thermal_restart_e2e.py
@@ -0,0 +1,157 @@
+"""Thermal state must survive the real production executable's checkpoint."""
+from __future__ import annotations
+
+import json
+import os
+from pathlib import Path
+import subprocess
+import tempfile
+import unittest
+
+from test_solid_conduction_contract_e2e import write_cube
+
+
+class ProductionThermalRestart(unittest.TestCase):
+    @classmethod
+    def setUpClass(cls):
+        cls.solver = Path(os.environ["CFDX_PRODUCTION_SOLVER"])
+        if not cls.solver.is_file():
+            raise RuntimeError(f"production executable is missing: {cls.solver}")
+
+    def setUp(self):
+        self.tmp = tempfile.TemporaryDirectory(prefix="cfdx_temperature_restart_")
+        self.root = Path(self.tmp.name)
+        self.case = self.root / "case.cfdx.h5"
+        def configure(setup):
+            setup.update(energy_model="conjugate", turbulence_model="laminar")
+            setup["materials"][0].update(density=1.0, dynamic_viscosity=0.1,
+                                         specific_heat=1.0)
+            setup["numerics"] = {"selection": {
+                "entries": [
+                    {"family": "gradient", "configuration_key": "numerics.gradient.gauss"},
+                    {"family": "convection", "configuration_key": "numerics.convection.upwind"},
+                    {"family": "pressure_velocity", "configuration_key": "pressure_velocity.simple"},
+                    {"family": "linear_solver", "configuration_key": "linear.fgmres"},
+                ],
+                "required_families": ["gradient", "convection", "pressure_velocity"],
+            }}
+        write_cube(self.case, configure)
+
+    def tearDown(self):
+        self.tmp.cleanup()
+
+    def run_solver(self, name, restart=None):
+        output = self.root / name
+        command = [str(self.solver), "--mesh", str(self.case),
+                   "--output-dir", str(output), "--iterations", "4"]
+        if restart is not None:
+            command += ["--restart", str(restart)]
+        result = subprocess.run(command, capture_output=True, text=True,
+                                timeout=30, check=False)
+        return result, output
+
+    def test_checkpoint_preserves_temperature_and_warm_energy_initialization(self):
+        first, first_dir = self.run_solver("cold")
+        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
+        checkpoint = first_dir / "restart.dat"
+        text = checkpoint.read_text()
+        self.assertIn("field T 1\n", text)
+        self.assertEqual(float(text.split("field T 1\n", 1)[1].split()[0]), 350.0)
+        # The checkpoint is the thermal initializer. A cold solve would reject
+        # this absent positive reference; a restart must consume its saved T.
+        import h5py
+        with h5py.File(self.case, "r+") as h5:
+            setup = json.loads(h5.attrs["case_setup_json"])
+            setup["ref_temperature"] = 0.0
+            h5.attrs["case_setup_json"] = json.dumps(setup)
+        second, second_dir = self.run_solver("warm", checkpoint)
+        self.assertEqual(second.returncode, 0, second.stdout + second.stderr)
+        warm_energy = json.loads((second_dir / "energy_convergence.json").read_text())
+        self.assertTrue(warm_energy["converged"])
+        self.assertTrue(warm_energy["temperature_admissibility_passed"])
+        self.assertEqual(warm_energy["nonpositive_temperature"], 0)
+        self.assertEqual(warm_energy["iterations"], 1)
+        self.assertEqual(float((second_dir / "restart.dat").read_text().split(
+            "field T 1\n", 1)[1].split()[0]), 350.0)
+
+    def test_thermal_restart_requires_temperature_instead_of_cold_fallback(self):
+        first, output = self.run_solver("cold")
+        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
+        checkpoint = self.root / "missing_temperature.dat"
+        text = (output / "restart.dat").read_text()
+        checkpoint.write_text(text.split("field T 1\n", 1)[0])
+        result, rejected_dir = self.run_solver("rejected", checkpoint)
+        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
+        self.assertIn("requested T but DAT has none", result.stderr)
+        self.assertFalse((rejected_dir / "restart.dat").exists())
+        execution = json.loads((rejected_dir / "execution.json").read_text())
+        self.assertFalse(execution["converged"])
+
+    def test_negative_restart_temperature_is_rejected_before_flow(self):
+        first, output = self.run_solver("cold")
+        self.assertEqual(first.returncode, 0, first.stdout + first.stderr)
+        checkpoint = self.root / "negative_temperature.dat"
+        text = (output / "restart.dat").read_text()
+        checkpoint.write_text(text.split("field T 1\n", 1)[0] + "field T 1\n-1\n")
+        result, rejected_dir = self.run_solver("rejected", checkpoint)
+        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
+        self.assertIn("restarted temperature", result.stderr)
+        self.assertFalse((rejected_dir / "restart.dat").exists())
+
+    def test_physical_transient_flag_cannot_be_solved_as_steady(self):
+        import h5py
+        with h5py.File(self.case, "r+") as h5:
+            setup = json.loads(h5.attrs["case_setup_json"])
+            setup.update(transient=True, time_step=0.1, end_time=1.0)
+            h5.attrs["case_setup_json"] = json.dumps(setup)
+        result, output = self.run_solver("rejected_time")
+        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
+        self.assertIn("transient evolution", result.stderr)
+        self.assertFalse((output / "restart.dat").exists())
+
+    def test_energy_exception_cannot_report_successful_carrier_as_overall_convergence(self):
+        import h5py
+        for missing in ("ref_temperature", "specific_heat"):
+            with self.subTest(missing=missing):
+                with h5py.File(self.case, "r+") as h5:
+                    setup = json.loads(h5.attrs["case_setup_json"])
+                    setup["ref_temperature"] = 300.0
+                    setup["materials"][0]["specific_heat"] = 1.0
+                    if missing == "ref_temperature":
+                        setup[missing] = 0.0
+                    else:
+                        setup["materials"][0].pop(missing)
+                    h5.attrs["case_setup_json"] = json.dumps(setup)
+                result, output = self.run_solver("missing_" + missing)
+                self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
+                execution = json.loads((output / "execution.json").read_text())
+                self.assertFalse(execution["converged"])
+                self.assertEqual(execution["process_exit_code"], 2)
+                self.assertIn(execution["convergence_reason"], result.stderr)
+
+    def test_conservative_negative_temperature_is_not_physically_accepted(self):
+        import h5py
+        with h5py.File(self.case, "r+") as h5:
+            setup = json.loads(h5.attrs["case_setup_json"])
+            # Two unit-area walls at 300/400 K, lambda=2 and d=0.5:
+            # 8*T = 2800 + Su. Su=-3600 yields exactly T=-100 K.
+            setup["materials"][0]["extra_properties"] = {"volumetric_power": "-3600"}
+            h5.attrs["case_setup_json"] = json.dumps(setup)
+        result, output = self.run_solver("negative_final")
+        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
+        energy = json.loads((output / "energy_convergence.json").read_text())
+        summary = json.loads((output / "energy_summary.json").read_text())
+        execution = json.loads((output / "execution.json").read_text())
+        self.assertTrue(energy["linear_solver_converged"])
+        self.assertTrue(energy["global_energy_balance_passed"])
+        self.assertTrue(energy["local_energy_balance_passed"])
+        self.assertFalse(energy["temperature_admissibility_passed"])
+        self.assertEqual(energy["nonpositive_temperature"], 1)
+        self.assertEqual(summary["min_temperature"], -100.0)
+        self.assertFalse(energy["converged"])
+        self.assertFalse(execution["converged"])
+        self.assertIn("temperature admissibility", execution["convergence_reason"])
+
+
+if __name__ == "__main__":
+    unittest.main(verbosity=2)
diff --git a/tests/python/test_solid_conduction_contract_e2e.py b/tests/python/test_solid_conduction_contract_e2e.py
new file mode 100644
index 00000000..c4b402f2
--- /dev/null
+++ b/tests/python/test_solid_conduction_contract_e2e.py
@@ -0,0 +1,115 @@
+"""Exercise physics selection through the real solid-conduction executable."""
+from __future__ import annotations
+
+import json
+import os
+from pathlib import Path
+import subprocess
+import tempfile
+import unittest
+
+import h5py
+import numpy as np
+
+
+def write_cube(path: Path, change=None) -> None:
+    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
+                       [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
+    faces = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4],
+             [3, 7, 6, 2], [0, 4, 7, 3], [1, 2, 6, 5]]
+    names = ["bottom", "top", "front", "back", "left", "right"]
+    setup = {
+        "energy_model": "solid_conduction", "solver_mode": "steady",
+        "ref_temperature": 300.0,
+        "materials": [{"name": "solid", "thermal_conductivity": 2.0}],
+        "boundary_conditions": [
+            {"patch_name": name, "type": "wall", "value_type": "wall_thermal",
+             "wall_temperature": 300.0 if name == "left" else 400.0 if name == "right" else 0.0}
+            for name in names],
+    }
+    if change is not None:
+        change(setup)
+    with h5py.File(path, "w") as h5:
+        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
+        h5.attrs["format_version"] = "1"
+        h5.attrs["schema_version"] = "1"
+        h5.attrs["mesh_topology"] = "cfdx-csr-v1"
+        h5.attrs["case_setup_json"] = json.dumps(setup)
+        h5.attrs["boundary_patches"] = ";".join(
+            f"{name}:{index}:1:0" for index, name in enumerate(names))
+        h5.create_dataset("points", data=points)
+        h5.create_dataset("face_vertices", data=np.array(faces, dtype=np.uint64).ravel())
+        h5.create_dataset("face_offsets", data=np.arange(0, 25, 4, dtype=np.uint64))
+        h5.create_dataset("owner", data=np.zeros(6, dtype=np.uint64))
+        h5.create_dataset("neighbour", data=np.full(6, -1, dtype=np.int64))
+        h5.create_dataset("cell_faces", data=np.arange(6, dtype=np.uint64))
+        h5.create_dataset("cell_offsets", data=np.array([0, 6], dtype=np.uint64))
+        h5.create_dataset("patch_face_ids", data=np.arange(6, dtype=np.uint64))
+        h5.create_dataset("patch_face_offsets", data=np.arange(7, dtype=np.uint64))
+
+
+class SolidConductionApplicationContract(unittest.TestCase):
+    @classmethod
+    def setUpClass(cls):
+        cls.solver = Path(os.environ["CFDX_THERMAL_SOLVER"])
+        if not cls.solver.is_file():
+            raise RuntimeError(f"thermal executable is missing: {cls.solver}")
+
+    def run_case(self, change=None):
+        with tempfile.TemporaryDirectory(prefix="cfdx_solid_contract_") as work:
+            root = Path(work)
+            write_cube(root / "case.cfdx.h5", change)
+            result = subprocess.run(
+                [str(self.solver), "--mesh", str(root / "case.cfdx.h5"),
+                 "--output-dir", str(root / "output")],
+                capture_output=True, text=True, timeout=30, check=False)
+            execution = json.loads((root / "output/execution.json").read_text())
+            summary_path = root / "output/thermal_summary.json"
+            summary = json.loads(summary_path.read_text()) if summary_path.exists() else None
+            return result, execution, summary
+
+    def assert_rejected(self, change, message):
+        result, execution, summary = self.run_case(change)
+        self.assertEqual(result.returncode, 2, result.stdout + result.stderr)
+        self.assertEqual(execution["process_exit_code"], 2)
+        self.assertFalse(execution["converged"])
+        self.assertIsNone(summary)
+        self.assertIn(message, result.stderr)
+
+    def test_supported_solid_case_solves(self):
+        result, execution, summary = self.run_case()
+        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
+        self.assertTrue(execution["converged"])
+        self.assertIsNotNone(summary)
+        self.assertIn("Temperature min 350 max 350", result.stdout)
+
+    def test_radiation_is_not_silently_ignored(self):
+        self.assert_rejected(lambda case: case.update(radiation_model="S2S"), "radiation")
+
+    def test_transient_evolution_is_not_silently_ignored(self):
+        self.assert_rejected(lambda case: case.update(transient=True, time_step=0.1), "transient")
+
+    def test_turbulence_is_not_silently_ignored(self):
+        self.assert_rejected(lambda case: case.update(turbulence_model="k_omega_sst"), "turbulence")
+
+    def test_canonical_temporal_selection_is_not_silently_ignored(self):
+        self.assert_rejected(lambda case: case.update(numerics={"selection": {
+            "entries": [{"family": "temporal", "configuration_key": "numerics.temporal.bdf2"}],
+            "required_families": ["temporal"],
+        }}), "transient evolution")
+
+    def test_multiphase_is_not_silently_ignored(self):
+        self.assert_rejected(lambda case: case.update(multiphase_model="VOF"), "multiphase")
+
+    def test_fluid_inlet_is_not_replaced_by_adiabatic_wall(self):
+        def change(case):
+            case["boundary_conditions"][0].update(type="inlet", velocity_magnitude=1.0)
+        self.assert_rejected(change, "boundary physics")
+
+    def test_absent_conductivity_is_not_replaced_by_air(self):
+        self.assert_rejected(lambda case: case["materials"][0].pop("thermal_conductivity"),
+                             "conductivity")
+
+
+if __name__ == "__main__":
+    unittest.main(verbosity=2)
diff --git a/tests/python/test_t2_campaign_verdict.py b/tests/python/test_t2_campaign_verdict.py
new file mode 100644
index 00000000..46d47005
--- /dev/null
+++ b/tests/python/test_t2_campaign_verdict.py
@@ -0,0 +1,69 @@
+"""Independent acceptance examples for the T2 campaign's frozen criteria."""
+from __future__ import annotations
+
+import importlib.util
+from pathlib import Path
+import unittest
+
+path = Path(__file__).resolve().parents[2] / "scripts/t2_heated_channel_calibration.py"
+spec = importlib.util.spec_from_file_location("t2_campaign", path)
+campaign = importlib.util.module_from_spec(spec)
+spec.loader.exec_module(campaign)
+
+
+def observations():
+    return [{"nx": nx, "ny": ny, "status": "PASS", "l2_poiseuille": 0.001,
+             "err_rel": error} for (nx, ny), error in
+            zip(campaign.MESHES, [0.16, 0.08, 0.04])]
+
+
+class CampaignVerdictTests(unittest.TestCase):
+    def test_complete_campaign_meets_physical_targets(self):
+        result = campaign.evaluate_campaign(observations(), True)
+        self.assertTrue(result["qualification_passed"])
+        self.assertEqual(result["verdict"], "physics_qualified")
+
+    def test_execution_success_does_not_qualify_wrong_nusselt(self):
+        data = observations()
+        data[-1]["err_rel"] = 0.0993
+        result = campaign.evaluate_campaign(data, True)
+        self.assertTrue(result["execution_passed"])
+        self.assertTrue(result["poiseuille_passed"])
+        self.assertFalse(result["nusselt_passed"])
+        self.assertFalse(result["qualification_passed"])
+        self.assertEqual(result["verdict"], "execution_complete_physics_not_qualified")
+
+    def test_poiseuille_gate_applies_to_every_mesh(self):
+        data = observations()
+        data[0]["l2_poiseuille"] = 0.006
+        self.assertFalse(campaign.evaluate_campaign(data, True)["qualification_passed"])
+
+    def test_nan_infinity_negative_and_missing_metrics_cannot_pass(self):
+        for metric in ("l2_poiseuille", "err_rel"):
+            for invalid in (float("nan"), float("inf"), -1.0, None):
+                with self.subTest(metric=metric, invalid=invalid):
+                    for index in (0, -1):
+                        data = observations()
+                        data[index][metric] = invalid
+                        self.assertFalse(campaign.evaluate_campaign(data, True)["qualification_passed"])
+
+    def test_missing_or_failed_mesh_does_not_qualify(self):
+        for data in ([], observations()[:-1], observations()[1:]):
+            self.assertFalse(campaign.evaluate_campaign(data, True)["execution_passed"])
+        data = observations()
+        data[1]["status"] = "FAIL"
+        self.assertFalse(campaign.evaluate_campaign(data, True)["qualification_passed"])
+        self.assertFalse(campaign.evaluate_campaign(observations(), False)["qualification_passed"])
+
+    def test_original_thresholds_are_inclusive_and_unchanged(self):
+        self.assertEqual(campaign.L2_GATE, 0.005)
+        self.assertEqual(campaign.ERR_TARGET, 0.05)
+        data = observations()
+        for entry in data:
+            entry["l2_poiseuille"] = 0.005
+        data[-1]["err_rel"] = 0.05
+        self.assertTrue(campaign.evaluate_campaign(data, True)["qualification_passed"])
+
+
+if __name__ == "__main__":
+    unittest.main(verbosity=2)
diff --git a/tests/unit/test_case_hdf5_io.cpp b/tests/unit/test_case_hdf5_io.cpp
index 62145c7f..67abfe36 100644
--- a/tests/unit/test_case_hdf5_io.cpp
+++ b/tests/unit/test_case_hdf5_io.cpp
@@ -124,6 +124,43 @@ static GapAnalysis make_test_gap() {
 
 int main() {
 
+    run_case("case_hdf5_missing_material_properties_remain_unspecified", []() {
+        const std::string filename = "/tmp/cfdx_case_unspecified_material.h5";
+        std::remove(filename.c_str());
+        EXPECT_TRUE(write_case_cfdx_h5(filename, make_unit_cube(),
+            make_test_source(), make_default_setup(), GapAnalysis{}));
+        const std::string json = R"({"materials":[{"name":"unspecified"}]})";
+        const hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
+        EXPECT_TRUE(file >= 0);
+        EXPECT_TRUE(H5Adelete(file, "case_setup_json") >= 0);
+        const hid_t type = H5Tcopy(H5T_C_S1);
+        EXPECT_TRUE(H5Tset_size(type, json.size() + 1) >= 0);
+        const hid_t space = H5Screate(H5S_SCALAR);
+        const hid_t attr = H5Acreate2(file, "case_setup_json", type, space,
+            H5P_DEFAULT, H5P_DEFAULT);
+        EXPECT_TRUE(attr >= 0);
+        EXPECT_TRUE(H5Awrite(attr, type, json.c_str()) >= 0);
+        H5Aclose(attr);
+        H5Sclose(space);
+        H5Tclose(type);
+        H5Fclose(file);
+        Mesh mesh;
+        SourceInfo source;
+        CaseSetup setup;
+        GapAnalysis gap;
+        EXPECT_TRUE(read_case_cfdx_h5(filename, mesh, source, setup, gap));
+        EXPECT_TRUE(setup.materials.size() == 1);
+        const auto& material = setup.materials.front();
+        EXPECT_TRUE(material.name == "unspecified");
+        EXPECT_TRUE(material.density == 0.0);
+        EXPECT_TRUE(material.dynamic_viscosity == 0.0);
+        EXPECT_TRUE(material.thermal_conductivity == 0.0);
+        EXPECT_TRUE(material.specific_heat == 0.0);
+        EXPECT_TRUE(material.molecular_weight == 0.0);
+        EXPECT_TRUE(material.eos_model.empty());
+        std::remove(filename.c_str());
+    });
+
     run_case("case_hdf5_roundtrip_mesh_and_metadata", []() {
         const std::string filename = "/tmp/cfdx_case_roundtrip.h5";
         std::remove(filename.c_str());
diff --git a/tests/unit/test_finite_volume_transport.cpp b/tests/unit/test_finite_volume_transport.cpp
index fac8dbcb..3b95dba5 100644
--- a/tests/unit/test_finite_volume_transport.cpp
+++ b/tests/unit/test_finite_volume_transport.cpp
@@ -168,6 +168,37 @@ static Mesh make_unit_cube()
 
 int main()
 {
+    run_case("flux_imbalance_metrics_do_not_hide_local_cancellation", [] {
+        const auto metrics = flux_imbalance_metrics({10.0, -10.0});
+        EXPECT_NEAR(metrics.signed_global, 0.0, 1e-14);
+        EXPECT_NEAR(metrics.global_l1, 20.0, 1e-14);
+        EXPECT_NEAR(metrics.local_linf, 10.0, 1e-14);
+    });
+
+    run_case("assembled_flux_imbalance_exposes_local_cancellation", [] {
+        const Mesh m = make_1d_chain(2);
+        const FvGeometry g = build_fv_geometry(m);
+        Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1);
+        Field<double,Location::CELL> su(m.n_cells(),"Su","K/s",1);
+        Field<double,Location::CELL> sp(m.n_cells(),"Sp","1/s",1);
+        flux.fill(0.0); su.fill(0.0); sp.fill(0.0);
+        std::size_t internal = m.n_faces();
+        for (std::size_t face = 0; face < m.n_faces(); ++face) {
+            if (m.ownership().neighbour(face) >= 0) {
+                internal = face;
+                break;
+            }
+        }
+        EXPECT_TRUE(internal < m.n_faces());
+        flux(internal) = 1.0;
+        ScalarBoundaryConditions bc;
+        const auto eq = assemble_scalar_equation(m,g,flux,1.0,su,sp,bc);
+        EXPECT_NEAR(eq.signed_global_imbalance, 0.0, 1e-14);
+        EXPECT_NEAR(eq.global_l1_imbalance, 2.0, 1e-14);
+        EXPECT_NEAR(eq.max_local_imbalance, 1.0, 1e-14);
+        EXPECT_NEAR(eq.max_imbalance, 1.0, 1e-14);
+    });
+
     run_case("geometric_interpolation_weight_uses_longitudinal_projection", [] {
         const Vec3 owner(0.0,0.0,0.0);
         const Vec3 neighbour(1.0,0.0,0.0);
@@ -255,6 +286,90 @@ int main()
         EXPECT_NEAR(scalar_equation_residual_inf(eq,x),0.0,1e-12);
     });
 
+    run_case("scalar_residual_pair_uses_same_physical_norm_and_scale", [] {
+        const auto check_identity_case = [](double a0, double a1, double b0,
+                                            double b1, double x0, double x1,
+                                            double expected_inf, double expected_relative,
+                                            cfdx::core::LinearScalingModel scaling) {
+            ScalarEquation eq;
+            eq.matrix = SparseMatrix(2,2);
+            eq.matrix.push_back(0,0,a0);
+            eq.matrix.push_back(1,1,a1);
+            eq.matrix.finalize();
+            eq.rhs = Vector(2);
+            eq.rhs(0)=b0; eq.rhs(1)=b1;
+            Vector x(2);
+            x(0)=x0; x(1)=x1;
+            ScalarSolveControls controls{20,1e-3,1.0};
+            controls.problem=cfdx::core::LinearProblemKind::Diffusion;
+            controls.linear_solver.krylov=cfdx::core::KrylovModel::CG;
+            controls.linear_solver.preconditioner=cfdx::core::PreconditionerModel::Jacobi;
+            controls.linear_solver.scaling=scaling;
+            const auto result=solve_scalar_equation(eq,x,controls);
+            EXPECT_TRUE(result.status==SolverStatus::CONVERGED);
+            EXPECT_NEAR(result.residual_before_relaxation,expected_inf,1e-10);
+            EXPECT_NEAR(result.residual_after_relaxation,expected_inf,1e-10);
+            EXPECT_NEAR(result.residual_relative_before_relaxation,expected_relative,1e-12);
+            EXPECT_NEAR(result.residual_relative_after_relaxation,expected_relative,1e-12);
+        };
+
+        // Initial candidate is already admissible under the solver tolerance,
+        // so alpha=1 leaves it unchanged. Both diagnostics must report the
+        // same physical L-infinity residual, not a Krylov L2 residual.
+        check_identity_case(1.0,1.0,3.0,4.0,2.9996,3.9996,
+                            0.0004,0.0001,cfdx::core::LinearScalingModel::None);
+
+        // Row scaling makes the algebraic system well-scaled but does not
+        // change the physical residual contract of the returned field.
+        check_identity_case(1.0e6,1.0,3.0e6,4.0,2.9996,3.9996,
+                            400.0,400.0/3.0e6,cfdx::core::LinearScalingModel::Row);
+    });
+
+    run_case("scalar_relaxation_residual_matches_closed_form_oracle", [] {
+        ScalarEquation eq;
+        eq.matrix=SparseMatrix(2,2);
+        eq.matrix.push_back(0,0,6.0); eq.matrix.push_back(0,1,-2.0);
+        eq.matrix.push_back(1,0,-2.0); eq.matrix.push_back(1,1,6.0);
+        eq.matrix.finalize();
+        eq.rhs=Vector(2); eq.rhs(0)=0.5; eq.rhs(1)=0.5;
+        Vector x(2,0.0);
+        ScalarSolveControls controls{20,1e-12,0.25};
+        controls.problem=cfdx::core::LinearProblemKind::Diffusion;
+        controls.linear_solver.krylov=cfdx::core::KrylovModel::CG;
+        controls.linear_solver.preconditioner=cfdx::core::PreconditionerModel::Jacobi;
+        const auto result=solve_scalar_equation(eq,x,controls);
+        EXPECT_TRUE(result.status==SolverStatus::CONVERGED);
+        EXPECT_NEAR(x(0),1.0/32.0,1e-14);
+        EXPECT_NEAR(x(1),1.0/32.0,1e-14);
+        EXPECT_NEAR(result.residual_after_relaxation,3.0/8.0,1e-14);
+        EXPECT_NEAR(result.residual_relative_after_relaxation,3.0/4.0,1e-14);
+    });
+
+    run_case("relaxed_scalar_solve_reports_returned_field_residual", [] {
+        const Mesh m = make_1d_chain(2);
+        const FvGeometry g = build_fv_geometry(m);
+        Field<double,Location::FACE> flux(m.n_faces(),"phi","kg/s",1);
+        Field<double,Location::CELL> su(m.n_cells(),"Su","K/s",1);
+        Field<double,Location::CELL> sp(m.n_cells(),"Sp","1/s",1);
+        flux.fill(0.0);
+        su.fill(1.0);
+        sp.fill(0.0);
+        ScalarBoundaryConditions bc;
+        bc["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
+        bc["outlet"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
+        const auto eq = assemble_scalar_equation(m,g,flux,1.0,su,sp,bc);
+        Vector x(m.n_cells(),0.0);
+        const auto r = solve_scalar_equation(eq,x,{100,1e-12,0.25});
+        EXPECT_TRUE(r.status == SolverStatus::CONVERGED);
+        EXPECT_NEAR(r.residual_before_relaxation, 0.0, 1e-12);
+        EXPECT_NEAR(r.residual_after_relaxation, scalar_equation_residual_inf(eq,x), 1e-14);
+        EXPECT_TRUE(r.residual_after_relaxation > r.residual_before_relaxation);
+        EXPECT_NEAR(r.residual_relative_before_relaxation, r.residual_relative, 1e-14);
+        EXPECT_NEAR(r.residual_relative_after_relaxation,
+                    r.residual_after_relaxation /
+                        std::max(std::abs(eq.rhs(0)), std::abs(eq.rhs(1))), 1e-14);
+    });
+
     run_case("exact_one_cell_solution_is_not_under_relaxed", [] {
         const Mesh m = make_unit_cube();
         const FvGeometry g = build_fv_geometry(m);
@@ -463,6 +578,44 @@ int main()
         }
     });
 
+    run_case("high_resolution_flux_boundedness_is_independent_of_source_correction", [] {
+        const Mesh m = make_1d_chain(3);
+        const FvGeometry g = build_fv_geometry(m);
+        Field<double,Location::FACE> flux(m.n_faces(),"mass_flux","kg/s",1);
+        Field<double,Location::CELL> su(m.n_cells(),"su","unit/s",1);
+        Field<double,Location::CELL> sp(m.n_cells(),"sp","1/s",1);
+        Field<double,Location::CELL> phi(m.n_cells(),"phi","unit",1);
+        flux.fill(0.0); su.fill(0.0); sp.fill(0.0);
+        phi(0)=0.0; phi(1)=1.0; phi(2)=1.01;
+        std::size_t target = m.n_faces();
+        for (std::size_t f = 0; f < m.n_faces(); ++f) {
+            if (m.ownership().owner(f) == 1 && m.ownership().neighbour(f) == 2) {
+                flux(f) = 1.0;
+                target = f;
+                break;
+            }
+        }
+        EXPECT_TRUE(target < m.n_faces());
+        const std::vector<double> extra_diagonal{1.0, 1.0, 1.0};
+        ScalarBoundaryConditions bc;
+        bc["inlet"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
+        bc["outlet"] = {ScalarBoundaryType::FIXED_VALUE, 1.01, 0.0};
+        bc["wall"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
+
+        for (const auto scheme : {ConvectionScheme::SECOND_ORDER_UPWIND, ConvectionScheme::TVD}) {
+            const auto equation = assemble_scalar_equation(
+                m, g, flux, 0.0, su, sp, bc, false, nullptr, &extra_diagonal,
+                nullptr, nullptr, scheme, &phi);
+            const auto reconstructed = reconstruct_scalar_transport_flux(
+                m, g, flux, phi, 0.0, bc, false, scheme);
+            EXPECT_NEAR(reconstructed(target), 1.01, 1e-12);
+            const auto upwind = assemble_scalar_equation(
+                m, g, flux, 0.0, su, sp, bc, false, nullptr, &extra_diagonal,
+                nullptr, nullptr, ConvectionScheme::UPWIND, &phi);
+            EXPECT_NEAR(equation.rhs(1) - upwind.rhs(1), -0.01, 1e-10);
+        }
+    });
+
     run_case("reconstructed_fixed_value_boundary_flux_matches_upwind_matrix", [] {
         const Mesh m = make_unit_cube();
         const FvGeometry g = build_fv_geometry(m);
diff --git a/tests/unit/test_m1_m4_physics.cpp b/tests/unit/test_m1_m4_physics.cpp
index 3c90c98a..d597a37b 100644
--- a/tests/unit/test_m1_m4_physics.cpp
+++ b/tests/unit/test_m1_m4_physics.cpp
@@ -4,6 +4,8 @@
 #include "cfdx/physics/thermal.h"
 #include "cfdx/physics/radiation.h"
 #include "cfdx/physics/radiation_advanced.h"
+#include "cfdx/physics/m1_m4_models.h"
+#include "cfdx/physics/radiation_s2s.h"
 #include "common/test_harness.h"
 #include <cmath>
 #include <limits>
@@ -14,6 +16,20 @@ using namespace cfdx::testing;
 
 int main()
 {
+    run_case("m1_m4_p1_source_uses_mean_intensity_normalization", [] {
+        // Independent dimensional reference: 4*pi*kappa*(Ib-J), with
+        // Ib=sigma*T^4/pi and J=0. The P1 production helper is checked
+        // separately below by the same physical convention.
+        constexpr double sigma = 5.670374419e-8;
+        constexpr double temperature = 1000.0;
+        const double ib = sigma * std::pow(temperature, 4) / M_PI;
+        // With zero irradiation, 4*pi*kappa*Ib simplifies independently to
+        // 4*sigma*T^4; do not reuse the implementation's source expression.
+        const double expected = 4.0 * sigma * std::pow(temperature, 4);
+        EXPECT_NEAR(cfdx::physics::m1m4::p1_source(1.0, 0.0, temperature), expected, 1e-8);
+        EXPECT_NEAR(cfdx::physics::m1m4::p1_source(1.0, ib, temperature), 0.0, 1e-8);
+    });
+
     run_case("pressure_velocity_relaxation_and_rhie_chow", [] {
         CouplingControls c;
         validate_coupling_controls(c);
@@ -23,7 +39,7 @@ int main()
 
     run_case("rans_and_les_models_are_positive", [] {
         EXPECT_NEAR(turbulent_kinematic_viscosity_kepsilon(4.0, 2.0), 0.72, 1e-12);
-        EXPECT_TRUE(turbulent_kinematic_viscosity_komega_sst(1.0, 2.0, 0.5) > 0.0);
+        EXPECT_TRUE(turbulent_kinematic_viscosity_komega_sst(1.0, 2.0, 0.5, 2.0, 0.31) > 0.0);
         EXPECT_NEAR(smagorinsky_eddy_viscosity(0.1, 20.0),
                     0.00578, 1e-12);
         EXPECT_TRUE(des_eddy_viscosity(0.1, 0.02, 20.0) >= 0.0);
@@ -138,6 +154,9 @@ int main()
             {{0,0,0},{0,0,1},1.0},
             {{0,0,1},{0,0,-1},1.0}};
         auto F=estimate_view_factor_matrix(patches);
+        EXPECT_NEAR(F[0], 0.0, 1e-14);
+        EXPECT_NEAR(F[3], 0.0, 1e-14);
+        validate_s2s_view_factors(F, std::vector<double>{1.0, 1.0}, 1e-10);
         EXPECT_TRUE(F[1] >= 0.0 && F[1] <= 1.0);
         EXPECT_TRUE(F[0] >= 0.0 && F[0] <= 1.0);
         EXPECT_TRUE(F[1] >= 0.0 && F[1] <= 1.0);
diff --git a/tests/unit/test_solid_conduction_selection.cpp b/tests/unit/test_solid_conduction_selection.cpp
new file mode 100644
index 00000000..02706dc4
--- /dev/null
+++ b/tests/unit/test_solid_conduction_selection.cpp
@@ -0,0 +1,104 @@
+#include "cfdx/physics/production_solid_conduction_contract.h"
+#include "common/test_harness.h"
+
+#include <stdexcept>
+
+using namespace cfdx::io;
+using namespace cfdx::physics;
+using namespace cfdx::testing;
+
+static CaseSetup solid_case()
+{
+    CaseSetup setup;
+    setup.energy_model = "solid_conduction";
+    setup.solver_mode = "steady";
+    BoundarySpec wall;
+    wall.patch_name = "wall";
+    wall.type = BCType::WALL;
+    wall.value_type = BCValueType::WALL_THERMAL;
+    wall.wall_temperature = 300.0;
+    setup.boundary_conditions.push_back(wall);
+    return setup;
+}
+
+static bool rejected(const CaseSetup& setup)
+{
+    try {
+        validate_production_solid_conduction_selection(setup);
+        return false;
+    } catch (const std::invalid_argument&) {
+        return true;
+    }
+}
+
+int main()
+{
+    run_case("solid_conduction_accepts_supported_physics_and_legacy_alias", [] {
+        auto setup = solid_case();
+        EXPECT_TRUE(!rejected(setup));
+        setup.energy_model = "conjugate";
+        setup.radiation_model = "off";
+        setup.turbulence_model = "laminar";
+        EXPECT_TRUE(!rejected(setup));
+    });
+    run_case("solid_conduction_does_not_silently_drop_radiation", [] {
+        auto setup = solid_case();
+        for (const auto* model : {"S2S", "P1", "DO", "solar", "???"}) {
+            setup.radiation_model = model;
+            EXPECT_TRUE(rejected(setup));
+        }
+    });
+    run_case("solid_conduction_does_not_silently_drop_transient_evolution", [] {
+        auto setup = solid_case();
+        setup.transient = true;
+        EXPECT_TRUE(rejected(setup));
+        setup = solid_case();
+        setup.solver_mode = "unsteady";
+        EXPECT_TRUE(rejected(setup));
+        setup.solver_mode = "???";
+        EXPECT_TRUE(rejected(setup));
+        setup = solid_case();
+        setup.numerics.transient_scheme = "unsteady_first_order";
+        EXPECT_TRUE(rejected(setup));
+        setup.numerics.transient_scheme = "???";
+        EXPECT_TRUE(rejected(setup));
+        setup = solid_case();
+        setup.numerical_config.entries.push_back({
+            cfdx::core::NumericalMethodFamily::Temporal, "numerics.temporal.bdf2"});
+        EXPECT_TRUE(rejected(setup));
+        setup.numerical_report = cfdx::core::resolve_case_numerics(
+            "solid_transient", setup.numerical_config);
+        setup.numerical_config.entries.clear();
+        EXPECT_TRUE(setup.numerical_report.valid());
+        EXPECT_TRUE(rejected(setup));
+    });
+    run_case("solid_conduction_rejects_unimplemented_physics", [] {
+        auto setup = solid_case();
+        setup.turbulence_model = "k_omega_sst";
+        EXPECT_TRUE(rejected(setup));
+        setup = solid_case();
+        setup.multiphase_model = "VOF";
+        EXPECT_TRUE(rejected(setup));
+        setup = solid_case();
+        setup.physics_model = "buoyant";
+        EXPECT_TRUE(rejected(setup));
+        for (const auto* physics : {"compressible_euler", "compressible_navier_stokes"}) {
+            setup = solid_case();
+            setup.physics_model = physics;
+            EXPECT_TRUE(rejected(setup));
+        }
+        setup = solid_case();
+        setup.energy_model = "isothermal";
+        EXPECT_TRUE(rejected(setup));
+    });
+    run_case("solid_conduction_does_not_replace_fluid_boundaries_with_adiabatic_walls", [] {
+        auto setup = solid_case();
+        for (const auto type : {BCType::INLET, BCType::OUTLET,
+                               BCType::PRESSURE_OUTLET, BCType::INTERFACE,
+                               BCType::PERIODIC, BCType::UNKNOWN}) {
+            setup.boundary_conditions.front().type = type;
+            EXPECT_TRUE(rejected(setup));
+        }
+    });
+    return run_all();
+}
diff --git a/tests/unit/test_sst_blending.cpp b/tests/unit/test_sst_blending.cpp
index 40a10725..b7f989d0 100644
--- a/tests/unit/test_sst_blending.cpp
+++ b/tests/unit/test_sst_blending.cpp
@@ -273,22 +273,20 @@ int main()
         EXPECT_NEAR(omega(0), 8.0, 0.0);
         EXPECT_NEAR(nu_t(0), 0.02, 0.0);
 
-        // The transport solver currently accepts this field's cell count and
-        // reads its first component; the production-closure contract rejects
-        // the non-scalar dimension after transport has converged. This forces
-        // the transactional wrapper's exception rollback path without a
-        // production-only fault-injection hook.
+        // Malformed dimensions must now be rejected before direct transport
+        // mutates k/omega, as well as by the transactional wrapper. Exception
+        // rollback after a real inner commit is exercised separately below.
         cfdx::core::Field<double, cfdx::core::Location::CELL> malformed_strain(
             mesh.n_cells(), "strain_bad_dimension", "1/s", 2);
         malformed_strain.fill(1.0);
         auto exception_direct_k = k;
         auto exception_direct_omega = omega;
-        const auto exception_direct = solve_sst_transport_dynamic_blending(
+        EXPECT_THROW(solve_sst_transport_dynamic_blending(
             mesh, geometry, flux, exception_direct_k, exception_direct_omega,
-            malformed_strain, wall_distance, controls, k_bcs, omega_bcs, 100, 1e-4);
-        EXPECT_TRUE(exception_direct.converged);
-        EXPECT_TRUE(exception_direct_k(0) != k(0) ||
-                    exception_direct_omega(0) != omega(0));
+            malformed_strain, wall_distance, controls, k_bcs, omega_bcs, 100, 1e-4),
+            std::invalid_argument);
+        EXPECT_NEAR(exception_direct_k(0), k(0), 0.0);
+        EXPECT_NEAR(exception_direct_omega(0), omega(0), 0.0);
 
         EXPECT_THROW(advance_sst_transport_transactionally(
             mesh, geometry, flux, k, omega, malformed_strain, wall_distance,
@@ -735,7 +733,9 @@ int main()
         const auto runtime = SSTRuntimeBoundaryContract{{}, {}, {"wall"}};
         EXPECT_THROW(require_supported_sst_runtime_boundaries(runtime),
                      std::invalid_argument);
-        require_supported_sst_runtime_boundaries(runtime, &face_values);
+        EXPECT_THROW(require_supported_sst_runtime_boundaries(runtime, &face_values),
+                     std::invalid_argument);
+        require_supported_sst_runtime_boundaries(runtime, mesh, geometry, &face_values);
     });
 
     run_case("sst_wall_omega_value_uses_resolved_wall_formula", [] {
diff --git a/tests/unit/test_turbulence_qualification_equations.cpp b/tests/unit/test_turbulence_qualification_equations.cpp
index 08becab7..4f95d1eb 100644
--- a/tests/unit/test_turbulence_qualification_equations.cpp
+++ b/tests/unit/test_turbulence_qualification_equations.cpp
@@ -3,13 +3,27 @@
 #include "cfdx/physics/turbulence_solver.h"
 #include "cfdx/physics/sst_solver.h"
 #include "cfdx/physics/spalart_allmaras.h"
+#include "common/test_harness.h"
 #include <algorithm>
 #include <cmath>
 #include <iostream>
+#include <limits>
 #include <stdexcept>
+#include <sstream>
 
 namespace {
-void expect_close(const char* label,double a,double b,double tol=1e-12) {    if(std::abs(a-b)>tol*std::max({1.0,std::abs(a),std::abs(b)})) {        std::cerr << "FAIL " << label << ": actual=" << a << " reference=" << b << " abs_error=" << std::abs(a-b) << "\n";
+void expect_close(const char* label, double a, double b, double tol = 1e-12) {
+    if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(tol) || tol < 0.0) {
+        std::cerr << "FAIL " << label
+                  << ": non-finite operand or invalid tolerance; actual=" << a
+                  << " reference=" << b << " tolerance=" << tol << "\n";
+        throw std::runtime_error("turbulence equation reference is non-finite or tolerance invalid");
+    }
+    const double error = std::abs(a - b);
+    const double limit = tol * std::max({1.0, std::abs(a), std::abs(b)});
+    if (!std::isfinite(error) || !std::isfinite(limit) || error > limit) {
+        std::cerr << "FAIL " << label << ": actual=" << a << " reference=" << b
+                  << " abs_error=" << error << " tolerance_limit=" << limit << "\n";
         throw std::runtime_error("turbulence equation reference mismatch");
     }
 }
@@ -17,6 +31,42 @@ void expect_close(const char* label,double a,double b,double tol=1e-12) {    if(
 }
 
 int main() {
+    // Prove the qualification comparator rejects non-finite values and invalid tolerances.
+    const auto expect_comparator_rejects = [](double actual, double reference, double tolerance) {
+        bool rejected = false;
+        try {
+            expect_close("comparator self-test", actual, reference, tolerance);
+        } catch (const std::runtime_error&) {
+            rejected = true;
+        }
+        if (!rejected)
+            throw std::runtime_error("numeric comparator accepted invalid input");
+    };
+    expect_comparator_rejects(std::numeric_limits<double>::quiet_NaN(), 1.0, 1e-12);
+    expect_comparator_rejects(std::numeric_limits<double>::infinity(), 1.0, 1e-12);
+    expect_comparator_rejects(1.0, 1.0, std::numeric_limits<double>::quiet_NaN());
+    expect_comparator_rejects(1.0, 1.0, -1.0);
+    expect_comparator_rejects(
+        2.0, 0.0, std::numeric_limits<double>::max());
+
+    const auto expect_near_rejects = [](double actual, double reference, double tolerance) {
+        bool rejected = false;
+        std::ostringstream suppressed_diagnostic;
+        auto* previous_buffer = std::cerr.rdbuf(suppressed_diagnostic.rdbuf());
+        try {
+            EXPECT_NEAR(actual, reference, tolerance);
+        } catch (const std::runtime_error&) {
+            rejected = true;
+        }
+        std::cerr.rdbuf(previous_buffer);
+        if (!rejected)
+            throw std::runtime_error("EXPECT_NEAR accepted invalid numeric input");
+    };
+    expect_near_rejects(std::numeric_limits<double>::quiet_NaN(), 1.0, 1e-12);
+    expect_near_rejects(1.0, std::numeric_limits<double>::infinity(), 1e-12);
+    expect_near_rejects(1.0, 1.0, std::numeric_limits<double>::quiet_NaN());
+    expect_near_rejects(1.0, 1.0, -1.0);
+
     using namespace cfdx::physics;
     TurbulenceTransportControls c;
     const double k=2.0, eps=4.0, omega=5.0, S=3.0, delta=0.2, volume=delta*delta*delta, y=0.05;
@@ -83,6 +133,16 @@ int main() {
         const double P=std::min(c.density*nut*S*S,c.sst_production_limiter*c.beta_star*c.density*k*omega);
         const double cross_term=2.0*(1.0-F1)*c.density*c.sst_sigma_w2/omega*cross;
         EXPECT_CLOSE(q.nut,nut); EXPECT_CLOSE(q.sk,P);        EXPECT_CLOSE(q.spk,-c.density*c.beta_star*omega);        EXPECT_CLOSE(q.sw,blend(c.gamma1,c.gamma2)*c.density*S*S+std::max(cross_term,0.0));        EXPECT_CLOSE(q.spw,-c.density*blend(c.beta1,c.beta2)*omega+std::min(cross_term,0.0)/omega);        EXPECT_CLOSE(q.gamma_k,c.density*(c.molecular_viscosity+blend(c.sst_sigma_k1,c.sst_sigma_k2)*nut));        EXPECT_CLOSE(q.gamma_w,c.density*(c.molecular_viscosity+blend(c.sst_sigma_w1,c.sst_sigma_w2)*nut));
+
+        // Independent SST-2003 production-limiter oracle. At S=100, raw
+        // Pk=77.5 is capped at 9; omega production must use Pk_limited/nut.
+        const auto limited=sst_cell_sources(k,omega,100.0,F1,F2,0.0,c);
+        if (!std::isfinite(limited.nut) || !std::isfinite(limited.sk) ||
+            !std::isfinite(limited.sw))
+            throw std::runtime_error("SST active-limiter oracle received a non-finite source");
+        EXPECT_CLOSE(limited.nut,0.00775);
+        EXPECT_CLOSE(limited.sk,9.0);
+        EXPECT_CLOSE(limited.sw,551.2258064516129);
     }
 
     c.model=TurbulenceModel::SPALART_ALLMARAS;
diff --git a/tests/validation/test_benchmark_matrix.cpp b/tests/validation/test_benchmark_matrix.cpp
index fc8ad3b9..235c4ab2 100644
--- a/tests/validation/test_benchmark_matrix.cpp
+++ b/tests/validation/test_benchmark_matrix.cpp
@@ -88,8 +88,8 @@ int main()
         EXPECT_NEAR(turbulent_kinematic_viscosity_kepsilon(1.0,2.0),0.045,1e-14);
     });
     run_case("M2_SST_eddy_viscosity_limits", [] {
-        EXPECT_NEAR(turbulent_kinematic_viscosity_komega_sst(0.5,2.0,0.0),0.25,1e-14);
-        EXPECT_TRUE(turbulent_kinematic_viscosity_komega_sst(0.5,2.0,1.0)>0.0);
+        EXPECT_NEAR(turbulent_kinematic_viscosity_komega_sst(0.5,2.0,0.0,0.0,0.31),0.25,1e-14);
+        EXPECT_TRUE(turbulent_kinematic_viscosity_komega_sst(0.5,2.0,1.0,2.0,0.31)>0.0);
     });
     run_case("M2_Smagorinsky_scaling", [] {
         const double n1=smagorinsky_eddy_viscosity(0.1,4.0);
diff --git a/tests/validation/test_cht_integration_contracts.cpp b/tests/validation/test_cht_integration_contracts.cpp
new file mode 100644
index 00000000..b9bfac7c
--- /dev/null
+++ b/tests/validation/test_cht_integration_contracts.cpp
@@ -0,0 +1,103 @@
+// A conductive CHT interface is conformal and impermeable, and both regions
+// advance the same physical step. Invalid contracts fail before any mutation.
+#include "audit/chain_fixture.h"
+#include "cfdx/physics/cht_solver.h"
+#include "common/test_harness.h"
+#include <limits>
+
+using namespace cfdx::core;
+using namespace cfdx::physics;
+using namespace cfdx::testing;
+
+struct PairFixture {
+    Mesh a=make_1d_chain(1), b=make_1d_chain(1);
+    FvGeometry ga,gb;
+    Field<double,Location::FACE> fa{a.n_faces(),"F","kg/s",1};
+    Field<double,Location::FACE> fb{b.n_faces(),"F","kg/s",1};
+    Field<double,Location::CELL> ta{1,"T","K",1},tb{1,"T","K",1};
+    Field<double,Location::CELL> qa{1,"Q","W/m3",1},qb{1,"Q","W/m3",1};
+    EnergySolverControls ea,eb;
+    ChtInterfaceControls c;
+    PairFixture()
+    {
+        for (std::size_t i=0;i<b.n_points();++i)
+            b.points().set(i,b.points().x(i)+1.0,b.points().y(i),b.points().z(i));
+        ga=build_fv_geometry(a); gb=build_fv_geometry(b);
+        fa.fill(0.0); fb.fill(0.0); qa.fill(0.0); qb.fill(0.0);
+        ta.fill(300.0); tb.fill(300.0);
+        ea.density=eb.density=1.0; ea.cp=eb.cp=1.0;
+        ea.conductivity=eb.conductivity=1.0; ea.dt=eb.dt=1.0;
+        ea.relaxation=eb.relaxation=1.0;
+        c.region1_patch="outlet"; c.region2_patch="inlet";
+        c.max_iterations=200;
+    }
+    ChtSolveResult solve()
+    {
+        return solve_two_region_cht(a,ga,b,gb,fa,fb,ta,tb,qa,qb,ea,eb,c);
+    }
+    void unchanged() const
+    {
+        EXPECT_NEAR(ta(0),300.0,0.0); EXPECT_NEAR(tb(0),300.0,0.0);
+    }
+};
+
+int main()
+{
+    run_case("cht_accepts_adjacent_opposed_interface", [] {
+        PairFixture f;
+        const auto pairs=match_cht_interface(f.a,f.ga,f.b,f.gb,f.c);
+        EXPECT_TRUE(pairs.size()==1);
+        EXPECT_TRUE(f.solve().converged);
+        f.unchanged();
+    });
+    run_case("cht_rejects_coincident_domains_with_parallel_outward_normals", [] {
+        PairFixture f;
+        f.c.region2_patch="outlet";
+        EXPECT_THROW_WITH(match_cht_interface(f.a,f.ga,f.a,f.ga,f.c),
+                          std::invalid_argument,"opposite outward normals");
+    });
+    run_case("cht_rejects_empty_interface", [] {
+        PairFixture f;
+        f.a.boundary().patch(1).face_ids.clear();
+        f.b.boundary().patch(0).face_ids.clear();
+        f.ga=build_fv_geometry(f.a); f.gb=build_fv_geometry(f.b);
+        EXPECT_THROW_WITH(f.solve(),std::invalid_argument,"non-empty");
+        f.unchanged();
+    });
+    run_case("cht_rejects_nonfinite_interface_controls", [] {
+        PairFixture f;
+        f.c.matching_tolerance=std::numeric_limits<double>::quiet_NaN();
+        EXPECT_THROW_WITH(f.solve(),std::invalid_argument,"interface controls");
+        f.unchanged();
+    });
+    run_case("cht_rejects_different_physical_steps", [] {
+        for (const double dt : {0.0,0.5,2.0}) {
+            PairFixture f;
+            f.eb.dt=dt;
+            EXPECT_THROW_WITH(f.solve(),std::invalid_argument,"same physical time step");
+            f.unchanged();
+        }
+    });
+    run_case("cht_rejects_mass_transport_at_impermeable_interface", [] {
+        for (const bool first : {true,false}) {
+            PairFixture f;
+            if(first) f.fa(1)=1e-12;
+            else f.fb(0)=-1e-12;
+            EXPECT_THROW_WITH(f.solve(),std::invalid_argument,"impermeable");
+            f.unchanged();
+        }
+    });
+    run_case("cht_rejects_malformed_second_region_before_first_solve", [] {
+        PairFixture f;
+        f.qb(0)=std::numeric_limits<double>::quiet_NaN();
+        EXPECT_THROW_WITH(f.solve(),std::invalid_argument,"region source");
+        f.unchanged();
+    });
+    run_case("cht_rejects_nonboundary_interface_face", [] {
+        PairFixture f;
+        f.a.ownership().set_neighbour(1,0);
+        EXPECT_THROW_WITH(f.solve(),std::invalid_argument,"boundary faces");
+        f.unchanged();
+    });
+    return run_all();
+}
diff --git a/tests/validation/test_cht_validation.cpp b/tests/validation/test_cht_validation.cpp
index 6e6342e4..32594fe3 100644
--- a/tests/validation/test_cht_validation.cpp
+++ b/tests/validation/test_cht_validation.cpp
@@ -6,13 +6,13 @@ using namespace cfdx::core;
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
@@ -20,7 +20,13 @@ static Mesh cube(const char* patch_name)
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
 
@@ -42,11 +48,54 @@ int main()
         EXPECT_TRUE(euclidean_distance > normal_distance);
     });
 
-    run_case("two_region_cht_independent_energy_conservation",[] {
+    run_case("cht_rejects_inconsistent_interface_conductivity", [] {
         Mesh m1=cube("interface1"),m2=cube("interface2");
         auto g1=build_fv_geometry(m1),g2=build_fv_geometry(m2);
         Field<double,Location::FACE> phi1(m1.n_faces(),"phi1","kg/s",1);
         Field<double,Location::FACE> phi2(m2.n_faces(),"phi2","kg/s",1);
+        phi1.fill(0.0); phi2.fill(0.0);
+        Field<double,Location::CELL> T1(1,"T1","K",1),T2(1,"T2","K",1);
+        Field<double,Location::CELL> q1(1,"q1","W/m3",1),q2(1,"q2","W/m3",1);
+        T1(0)=400.0; T2(0)=300.0; q1.fill(0.0); q2.fill(0.0);
+        EnergySolverControls e1; e1.conductivity=1.0;
+        EnergySolverControls e2=e1; e2.conductivity=2.0;
+        ChtInterfaceControls c;
+        c.region1_patch="interface1"; c.region2_patch="interface2";
+        c.conductivity1=1.0; c.conductivity2=1.0;
+        bool rejected=false;
+        try {
+            (void)solve_two_region_cht(
+                m1,g1,m2,g2,phi1,phi2,T1,T2,q1,q2,e1,e2,c);
+        } catch(const std::invalid_argument&) {
+            rejected=true;
+        }
+        EXPECT_TRUE(rejected);
+        EXPECT_NEAR(T1(0),400.0,0.0);
+        EXPECT_NEAR(T2(0),300.0,0.0);
+
+        // A malformed second region must be rejected before the first region
+        // is solved. Conductivities now match, isolating the early validation.
+        e2.conductivity=1.0;
+        e2.density=0.0;
+        bool invalid_region_rejected=false;
+        try {
+            (void)solve_two_region_cht(
+                m1,g1,m2,g2,phi1,phi2,T1,T2,q1,q2,e1,e2,c);
+        } catch(const std::invalid_argument&) {
+            invalid_region_rejected=true;
+        }
+        EXPECT_TRUE(invalid_region_rejected);
+        EXPECT_NEAR(T1(0),400.0,0.0);
+        EXPECT_NEAR(T2(0),300.0,0.0);
+    });
+
+    run_case("two_region_cht_independent_energy_conservation",[] {
+        // Adjacent regions, one common face and opposite outward normals.
+        // The old fixture overlapped both cubes and coupled all six faces.
+        Mesh m1=cube("interface1",true,false),m2=cube("interface2",true,true);
+        auto g1=build_fv_geometry(m1),g2=build_fv_geometry(m2);
+        Field<double,Location::FACE> phi1(m1.n_faces(),"phi1","kg/s",1);
+        Field<double,Location::FACE> phi2(m2.n_faces(),"phi2","kg/s",1);
         phi1.fill(0);phi2.fill(0);
         Field<double,Location::CELL> T1(1,"T1","K",1),T2(1,"T2","K",1);
         Field<double,Location::CELL> q1(1,"q1","W/m3",1),q2(1,"q2","W/m3",1);
@@ -65,8 +114,10 @@ int main()
 
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
@@ -91,7 +142,7 @@ int main()
 
         // Independent interface heat-flow equality: each side must carry the
         // same physical heat rate, with opposite outward signs.
-        const double interface_area=g1.face_area_vectors[0].mag()*static_cast<double>(m1.boundary().patch(0).face_ids.size());
+        const double interface_area=g1.face_area_vectors[face1].mag();
         const double qflux1=c.conductivity1*(T1(0)-tint)/d1;
         const double qflux2=c.conductivity2*(tint-T2(0))/d2;
         EXPECT_TRUE(std::isfinite(qflux1) && std::isfinite(qflux2));
diff --git a/tests/validation/test_energy_integration_contracts.cpp b/tests/validation/test_energy_integration_contracts.cpp
new file mode 100644
index 00000000..88a0e68d
--- /dev/null
+++ b/tests/validation/test_energy_integration_contracts.cpp
@@ -0,0 +1,399 @@
+// Independent energy oracles for implicit sources and variable conductivity.
+// The baseline macro runs only contracts expressible by the previous API.
+#include "audit/chain_fixture.h"
+#include "cfdx/physics/production_energy_stage.h"
+#include "common/test_harness.h"
+#include <limits>
+#ifndef CFDX_THERMAL_CONTRACTS_BASELINE
+#include "audit/sst_channel_fixture.h"
+#include "cfdx/physics/production_buoyancy_stage.h"
+#endif
+
+using namespace cfdx::core;
+using namespace cfdx::physics;
+using namespace cfdx::testing;
+
+static EnergySolverControls energy_controls()
+{
+    EnergySolverControls c;
+    c.density=1.0; c.cp=1.0; c.conductivity=1.0;
+    c.relaxation=1.0; c.tolerance=1e-12;
+    return c;
+}
+
+static cfdx::io::CaseSetup thermal_case()
+{
+    cfdx::io::CaseSetup setup;
+    setup.energy_model="conjugate";
+    cfdx::io::MaterialSpec material;
+    material.name="test_fluid"; material.density=2.0;
+    material.dynamic_viscosity=0.001;
+    material.thermal_conductivity=1.0; material.specific_heat=10.0;
+    setup.materials.push_back(material);
+    cfdx::io::BoundarySpec inlet;
+    inlet.patch_name="inlet"; inlet.type=cfdx::io::BCType::WALL;
+    inlet.value_type=cfdx::io::BCValueType::WALL_THERMAL;
+    inlet.wall_temperature=300.0;
+    setup.boundary_conditions.push_back(inlet);
+    cfdx::io::BoundarySpec outlet;
+    outlet.patch_name="outlet"; outlet.type=cfdx::io::BCType::WALL;
+    setup.boundary_conditions.push_back(outlet);
+    cfdx::io::BoundarySpec wall;
+    wall.patch_name="wall"; wall.type=cfdx::io::BCType::WALL;
+    wall.raw_params["heat_flux"]="10";
+    setup.boundary_conditions.push_back(wall);
+    return setup;
+}
+
+int main()
+{
+    run_case("implicit_source_equilibrium_uses_final_temperature", [] {
+        const auto mesh=make_1d_chain(1);
+        const auto g=build_fv_geometry(mesh);
+        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
+        Field<double,Location::CELL> T(1,"T","K",1),Su(1,"Su","W/m3",1),Sp(1,"Sp","W/m3/K",1);
+        phi.fill(0.0); T(0)=200.0; Su(0)=600.0; Sp(0)=-2.0;
+        auto c=energy_controls(); c.conductivity=0.0;
+        const auto r=solve_energy(mesh,g,phi,T,Su,c,{},nullptr,&Sp);
+        EXPECT_TRUE(r.converged);
+        EXPECT_NEAR(T(0),300.0,1e-12);
+        EXPECT_NEAR(r.history.back().independent_conservation_residual,0.0,1e-12);
+        EXPECT_NEAR(r.history.back().independent_conservation_normalized,0.0,1e-12);
+    });
+
+    run_case("implicit_source_underrelaxed_state_remains_out_of_balance", [] {
+        const auto mesh=make_1d_chain(1);
+        const auto g=build_fv_geometry(mesh);
+        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
+        Field<double,Location::CELL> T(1,"T","K",1),Su(1,"Su","W/m3",1),Sp(1,"Sp","W/m3/K",1);
+        phi.fill(0.0); T(0)=200.0; Su(0)=600.0; Sp(0)=-2.0;
+        auto c=energy_controls(); c.conductivity=0.0;
+        c.relaxation=0.5; c.max_iterations=1;
+        const auto r=solve_energy(mesh,g,phi,T,Su,c,{},nullptr,&Sp);
+        EXPECT_FALSE(r.converged);
+        EXPECT_NEAR(T(0),250.0,1e-12);
+        // No boundary flux or accumulation: residual = -(600 - 2*250).
+        EXPECT_NEAR(r.history.back().independent_conservation_residual,-100.0,1e-12);
+        EXPECT_NEAR(r.history.back().independent_conservation_normalized,1.0,1e-12);
+    });
+
+#ifndef CFDX_THERMAL_CONTRACTS_BASELINE
+    run_case("production_temperature_admissibility_is_separate_from_conservation", [] {
+        const auto mesh=make_1d_chain(1);
+        const auto g=build_fv_geometry(mesh);
+        auto setup=thermal_case();
+        for (auto& bc : setup.boundary_conditions) {
+            bc.type=cfdx::io::BCType::WALL;
+            bc.value_type=cfdx::io::BCValueType::WALL_THERMAL;
+            bc.wall_temperature=300.0;
+            bc.raw_params.clear();
+        }
+        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1); phi.fill(0.0);
+        Field<double,Location::CELL> T(1,"T","K",1);
+        // Six unit-area walls at distance .5: 12*T = 12*300 + Su*V.
+        // Do not clip the mathematical solution: report it and fail the
+        // physical stage if that absolute temperature is not positive.
+        for (const int power : {-4800,-3600,-3588}) {
+            setup.materials[0].extra_properties["volumetric_power"]=std::to_string(power);
+            T(0)=300.0;
+            const auto r=run_production_energy_stage(mesh,g,setup,phi,T);
+            const double exact=300.0+power/12.0;
+            std::cout << "    admissibility oracle: Su=" << power << " T=" << T(0)
+                      << " linear=" << r.linear_solver_converged
+                      << " global=" << r.global_energy_balance_passed
+                      << " local=" << r.local_energy_balance_passed
+                      << " residual=" << r.independent_conservation_residual
+                      << " converged=" << r.converged << '\n';
+            EXPECT_NEAR(T(0),exact,1e-12);
+            EXPECT_NEAR(r.min_temperature,exact,1e-12);
+            EXPECT_TRUE(r.linear_solver_converged);
+            EXPECT_TRUE(r.global_energy_balance_passed && r.local_energy_balance_passed);
+            EXPECT_NEAR(r.independent_conservation_residual,0.0,1e-12);
+            EXPECT_TRUE(r.converged==(exact>0.0));
+#ifndef CFDX_THERMAL_ADMISSIBILITY_BASELINE
+            EXPECT_TRUE(r.temperature_admissibility_passed==(exact>0.0));
+            EXPECT_TRUE(r.nonpositive_temperature==(exact>0.0 ? 0u : 1u));
+#endif
+        }
+        // The lower-level mathematical solver still supports relative or
+        // manufactured temperature fields, including non-positive values.
+        auto controls=energy_controls();
+        ScalarBoundaryConditions bcs;
+        for (const auto* name : {"inlet","outlet","wall"})
+            bcs[name]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
+        Field<double,Location::CELL> source(1,"Su","W/m3",1); source(0)=-4800.0;
+        T(0)=300.0;
+        const auto kernel=solve_energy(mesh,g,phi,T,source,controls,bcs);
+        EXPECT_TRUE(kernel.converged);
+        EXPECT_NEAR(T(0),-100.0,1e-12);
+    });
+
+    run_case("independent_source_audit_detects_perturbation_and_rejects_invalid_sp", [] {
+        const auto mesh=make_1d_chain(1);
+        const auto g=build_fv_geometry(mesh);
+        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
+        Field<double,Location::CELL> T(1,"T","K",1),Su(1,"Su","W/m3",1),Sp(1,"Sp","W/m3/K",1);
+        phi.fill(0.0); T(0)=310.0; Su(0)=600.0; Sp(0)=-2.0;
+        const auto c=energy_controls();
+        const auto audit=reconstruct_energy_balance(mesh,g,phi,T,T,Su,c,{},nullptr,&Sp);
+        EXPECT_NEAR(audit.source_integral,-20.0,1e-12);
+        EXPECT_NEAR(audit.residual,20.0,1e-12);
+        EXPECT_NEAR(audit.max_cell_residual,20.0,1e-12);
+        EXPECT_NEAR(audit.normalized_residual,1.0,1e-12);
+        EXPECT_NEAR(energy_balance_relative(mesh,g,phi,T,T,Su,c,{},nullptr,nullptr,&Sp),1.0,1e-12);
+        for (const double invalid : {1.0, std::numeric_limits<double>::quiet_NaN(),
+                                     std::numeric_limits<double>::infinity()}) {
+            Sp(0)=invalid;
+            EXPECT_THROW_WITH(reconstruct_energy_balance(mesh,g,phi,T,T,Su,c,{},nullptr,&Sp),
+                              std::invalid_argument,"implicit source");
+        }
+        Field<double,Location::CELL> short_sp;
+        EXPECT_THROW_WITH(reconstruct_energy_balance(mesh,g,phi,T,T,Su,c,{},nullptr,&short_sp),
+                          std::invalid_argument,"implicit source field");
+    });
+
+    run_case("two_material_conduction_matches_series_resistance_and_both_audits", [] {
+        const auto mesh=make_1d_chain(2);
+        const auto g=build_fv_geometry(mesh);
+        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
+        Field<double,Location::CELL> T(2,"T","K",1),Su(2,"Su","W/m3",1);
+        phi.fill(0.0); T.fill(300.0); Su.fill(0.0);
+        const auto old=T;
+        const std::vector<double> k{1.0,4.0};
+        auto c=energy_controls(); c.cp=17.0;
+        ScalarBoundaryConditions bcs;
+        bcs["inlet"]={ScalarBoundaryType::FIXED_VALUE,400.0,0.0};
+        bcs["outlet"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
+        const auto r=solve_energy(mesh,g,phi,T,Su,c,bcs,nullptr,nullptr,nullptr,&k);
+        EXPECT_TRUE(r.converged);
+        // R=0.5/1+0.5/4, q=100/R=160; cell centres lie at x=.25,.75.
+        EXPECT_NEAR(T(0),360.0,1e-10);
+        EXPECT_NEAR(T(1),310.0,1e-10);
+        const auto flux=reconstruct_energy_face_flux(mesh,g,phi,T,c,bcs,nullptr,&k);
+        EXPECT_NEAR(flux(0),-160.0,1e-10);
+        EXPECT_NEAR(flux(1),160.0,1e-10);
+        const auto audit=reconstruct_energy_balance(mesh,g,phi,T,old,Su,c,bcs,nullptr,nullptr,&k);
+        EXPECT_NEAR(audit.residual,0.0,1e-9);
+        EXPECT_TRUE(audit.max_cell_residual<1e-9);
+        EXPECT_TRUE(energy_balance_relative(mesh,g,phi,T,old,Su,c,bcs,nullptr,&k)<1e-9);
+        T(1)+=1.0;
+        const auto disturbed=reconstruct_energy_balance(mesh,g,phi,T,old,Su,c,bcs,nullptr,nullptr,&k);
+        EXPECT_TRUE(disturbed.max_cell_residual>10.0);
+    });
+
+    run_case("turbulent_conductivity_is_laminar_limit_and_prandtl_scaled", [] {
+        const auto ec=build_production_energy_case(thermal_case());
+        Field<double,Location::CELL> nu(2,"nu_t","m2/s",1);
+        nu(0)=0.0; nu(1)=0.15;
+        const auto k=build_effective_thermal_conductivity(ec,nu,1.0);
+        EXPECT_NEAR(k[0],1.0,0.0);
+        EXPECT_NEAR(k[1],4.0,1e-14);
+        const auto k_half_pr=build_effective_thermal_conductivity(ec,nu,0.5);
+        EXPECT_NEAR(k_half_pr[1],7.0,1e-14);
+        for (const double invalid : {0.0,-1.0,std::numeric_limits<double>::quiet_NaN(),
+                                     std::numeric_limits<double>::infinity()}) {
+            EXPECT_THROW_WITH(build_effective_thermal_conductivity(ec,nu,invalid),
+                              std::invalid_argument,"turbulent Prandtl");
+        }
+        for (const double invalid : {-1.0,std::numeric_limits<double>::quiet_NaN(),
+                                     std::numeric_limits<double>::infinity()}) {
+            nu(1)=invalid;
+            EXPECT_THROW_WITH(build_effective_thermal_conductivity(ec,nu,1.0),
+                              std::invalid_argument,"turbulent viscosity");
+        }
+        nu(1)=std::numeric_limits<double>::max();
+        EXPECT_THROW_WITH(build_effective_thermal_conductivity(ec,nu,1.0),
+                          std::invalid_argument,"effective conductivity");
+    });
+
+    run_case("production_flux_bc_uses_each_owner_conductivity", [] {
+        const auto mesh=make_1d_chain(2);
+        const auto g=build_fv_geometry(mesh);
+        const auto setup=thermal_case();
+        const auto ec=build_production_energy_case(setup);
+        const std::vector<double> k{1.0,4.0};
+        const auto bcs=build_production_temperature_bcs(mesh,setup,ec);
+        const auto faces=build_production_temperature_face_conditions(mesh,bcs,ec,k);
+        for (const auto f : mesh.boundary().patch(2).face_ids) {
+            const auto owner=mesh.ownership().owner(f);
+            EXPECT_TRUE(faces.has(f));
+            EXPECT_NEAR(faces.conditions[f].value,10.0/k[owner],1e-14);
+        }
+        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
+        Field<double,Location::CELL> T(2,"T","K",1);
+        phi.fill(0.0); T.fill(300.0);
+        const auto r=run_production_energy_stage(mesh,g,setup,phi,T,{},&k);
+        EXPECT_TRUE(r.converged);
+        // Each cell receives 20 W from four side faces. The left boundary
+        // removes 40 W: T0=300+40/4, T1=T0+20/(1/(.25/1+.25/4)).
+        EXPECT_NEAR(T(0),310.0,1e-10);
+        EXPECT_NEAR(T(1),316.25,1e-10);
+        EXPECT_NEAR(r.patch_flux_out.at("wall"),-40.0,1e-10);
+        EXPECT_NEAR(r.patch_flux_out.at("inlet"),40.0,1e-9);
+        EXPECT_TRUE(r.global_energy_balance_passed && r.local_energy_balance_passed);
+    });
+
+    run_case("variable_conductivity_rejected_before_temperature_mutation", [] {
+        const auto mesh=make_1d_chain(2);
+        const auto g=build_fv_geometry(mesh);
+        Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1);
+        Field<double,Location::CELL> T(2,"T","K",1),Su(2,"Su","W/m3",1);
+        phi.fill(0.0); T.fill(300.0); Su.fill(0.0);
+        auto c=energy_controls(); c.dt=1.0;
+        for (const std::vector<double>& invalid : {
+                 std::vector<double>{1.0},std::vector<double>{1.0,-1.0},
+                 std::vector<double>{1.0,std::numeric_limits<double>::quiet_NaN()},
+                 std::vector<double>{1.0,std::numeric_limits<double>::infinity()}}) {
+            EXPECT_THROW_WITH(solve_energy(mesh,g,phi,T,Su,c,{},nullptr,nullptr,nullptr,&invalid),
+                              std::invalid_argument,"conductivity");
+            EXPECT_NEAR(T(0),300.0,0.0); EXPECT_NEAR(T(1),300.0,0.0);
+        }
+        const std::vector<double> zero_boundary_k{1.0,0.0};
+        EXPECT_THROW_WITH(run_production_energy_stage(mesh,g,thermal_case(),phi,T,{},&zero_boundary_k),
+                          std::invalid_argument,"positive");
+        EXPECT_NEAR(T(0),300.0,0.0); EXPECT_NEAR(T(1),300.0,0.0);
+        // The flux audit also requires a scalar temperature. Reject an
+        // incompatible input before a successful solve can change component 0.
+        Field<double,Location::CELL> vector_T(2,"T","K",2); vector_T.fill(300.0);
+        Su.fill(10.0);
+        EXPECT_THROW_WITH(solve_energy(mesh,g,phi,vector_T,Su,c),
+                          std::invalid_argument,"energy field size");
+        for (std::size_t cell=0;cell<2;++cell)
+            for (std::size_t d=0;d<2;++d) EXPECT_NEAR(vector_T(cell,d),300.0,0.0);
+    });
+
+    run_case("accepted_sst_state_drives_energy_with_authoritative_flux", [] {
+        SSTChannelCase carrier;
+        const auto flow=solve_steady_incompressible(
+            carrier.mesh,carrier.U,carrier.p,carrier.ubc,carrier.pbc,carrier.controls);
+        // One accepted coupled iteration proves this software seam. It is
+        // deliberately not a claim of a fully converged turbulent channel.
+        EXPECT_TRUE(flow.history.size()==1);
+        EXPECT_TRUE(carrier.state.k(0)<0.01 && carrier.state.k(0)>carrier.state.k(2));
+        EXPECT_NEAR(flow.authoritative_mass_flux(0),-1.2,1e-12);
+        EXPECT_NEAR(flow.authoritative_mass_flux(1),1.2,1e-12);
+        auto setup=thermal_case();
+        setup.materials[0].density=1.2;
+        setup.materials[0].dynamic_viscosity=0.012;
+        setup.materials[0].extra_properties["volumetric_power"]="30";
+        setup.boundary_conditions[0].type=cfdx::io::BCType::INLET;
+        setup.boundary_conditions[0].value_type=cfdx::io::BCValueType::FIXED;
+        setup.boundary_conditions[0].wall_temperature=0.0;
+        setup.boundary_conditions[0].temperature=300.0;
+        setup.boundary_conditions[1].type=cfdx::io::BCType::OUTLET;
+        setup.boundary_conditions[2].raw_params.clear();
+        const auto ec=build_production_energy_case(setup);
+        const auto k=build_effective_thermal_conductivity(ec,carrier.state.nu_t,0.9);
+        for (std::size_t cell=0;cell<3;++cell) {
+            EXPECT_NEAR(k[cell],1.0+1.2*10.0*carrier.state.nu_t(cell)/0.9,1e-14);
+            EXPECT_TRUE(k[cell]>1.0);
+        }
+        Field<double,Location::CELL> T(3,"T","K",1); T.fill(300.0);
+        auto laminar=T;
+        const auto g=build_fv_geometry(carrier.mesh);
+        const auto r=run_production_energy_stage(carrier.mesh,g,setup,
+            flow.authoritative_mass_flux,T,{},&k);
+        const auto r_laminar=run_production_energy_stage(carrier.mesh,g,setup,
+            flow.authoritative_mass_flux,laminar);
+        EXPECT_TRUE(r.converged && r_laminar.converged);
+        // Independent 3x3 upwind/series-resistance system: dx=1/3,
+        // m_dot*cp=12 W/K and Su*V=10 W in every cell. No solver assembly
+        // or energy face-reconstruction helper is used by this oracle.
+        const double d_left=6.0*k[0];
+        const double d01=1.0/(1.0/(6.0*k[0])+1.0/(6.0*k[1]));
+        const double d12=1.0/(1.0/(6.0*k[1])+1.0/(6.0*k[2]));
+        const double lower[3]={0.0,-(12.0+d01),-(12.0+d12)};
+        double diagonal[3]={d_left+d01+12.0,d01+d12+12.0,d12+12.0};
+        const double upper[3]={-d01,-d12,0.0};
+        double exact[3]={10.0+(d_left+12.0)*300.0,10.0,10.0};
+        for (std::size_t cell=1;cell<3;++cell) {
+            const double factor=lower[cell]/diagonal[cell-1];
+            diagonal[cell]-=factor*upper[cell-1];
+            exact[cell]-=factor*exact[cell-1];
+        }
+        exact[2]/=diagonal[2];
+        for (int cell=1;cell>=0;--cell)
+            exact[cell]=(exact[cell]-upper[cell]*exact[cell+1])/diagonal[cell];
+        double max_effect=0.0;
+        for (std::size_t cell=0;cell<3;++cell) {
+            EXPECT_NEAR(T(cell),exact[cell],1e-9);
+            max_effect=std::max(max_effect,std::abs(T(cell)-laminar(cell)));
+        }
+        EXPECT_TRUE(max_effect>1e-4);
+        EXPECT_NEAR(r.boundary_flux_out,30.0,1e-8);
+        EXPECT_TRUE(r.local_energy_balance_passed && r.global_energy_balance_passed);
+        std::cout << "    SST-energy seam: max thermal effect=" << max_effect
+                  << " K; balance=" << r.independent_conservation_residual << " W\n";
+    });
+
+    run_case("restored_temperature_drives_the_first_buoyancy_source", [] {
+        SSTChannelCase carrier;
+        carrier.controls.sst_runtime=nullptr;
+        auto setup=thermal_case();
+        setup.energy_model="buoyant";
+        setup.ref_temperature=300.0;
+        setup.gravity_vector="0 -9.81 0";
+        setup.materials[0].density=1.2;
+        setup.materials[0].dynamic_viscosity=0.012;
+        setup.materials[0].extra_properties["thermal_expansion_coefficient"]="0.001";
+        Field<double,Location::CELL> loaded_T(3,"T","K",1);
+        for (std::size_t cell=0;cell<3;++cell) loaded_T(cell)=310.0+10.0*cell;
+        auto expected_U=carrier.U; auto expected_p=carrier.p;
+        auto reference_controls=carrier.controls;
+        // Build the oracle directly, independently of the thermal seeding
+        // and of build_cell_buoyancy_acceleration.
+        reference_controls.cell_body_force={
+            {0.0,0.001*10.0*9.81,0.0},
+            {0.0,0.001*20.0*9.81,0.0},
+            {0.0,0.001*30.0*9.81,0.0}};
+        const auto reference=solve_steady_incompressible(carrier.mesh,expected_U,expected_p,
+            carrier.ubc,carrier.pbc,reference_controls);
+        const auto result=run_production_buoyant_outer_loop(carrier.mesh,carrier.ubc,carrier.pbc,
+            carrier.controls,setup,carrier.U,carrier.p,{},&loaded_T);
+        EXPECT_TRUE(result.outer_iterations==1 && reference.history.size()==1);
+        double vertical_motion=0.0;
+        for (std::size_t cell=0;cell<3;++cell) {
+            for (std::size_t d=0;d<3;++d)
+                EXPECT_NEAR(carrier.U(cell,d),expected_U(cell,d),1e-13);
+            EXPECT_NEAR(carrier.p(cell),expected_p(cell),1e-13);
+            vertical_motion=std::max(vertical_motion,std::abs(carrier.U(cell,1)));
+            EXPECT_NEAR(loaded_T(cell),310.0+10.0*cell,0.0);
+        }
+        EXPECT_TRUE(vertical_motion>1e-4);
+        std::cout << "    restored buoyancy: first-iterate max Uy=" << vertical_motion << " m/s\n";
+    });
+
+    run_case("invalid_restored_temperature_is_rejected_before_carrier_mutation", [] {
+        SSTChannelCase carrier;
+        carrier.controls.sst_runtime=nullptr;
+        auto setup=thermal_case();
+        setup.energy_model="buoyant"; setup.ref_temperature=300.0;
+        setup.gravity_vector="0 -9.81 0";
+        setup.materials[0].extra_properties["thermal_expansion_coefficient"]="0.001";
+        Field<double,Location::CELL> loaded_T(3,"T","K",1); loaded_T.fill(300.0);
+        std::size_t predictor_calls=0;
+        carrier.controls.momentum_predictor_callback=[&](std::size_t,const auto&,const auto&,const auto&) {
+            ++predictor_calls;
+        };
+        for (const double invalid : {0.0,-1.0,std::numeric_limits<double>::quiet_NaN(),
+                                     std::numeric_limits<double>::infinity()}) {
+            loaded_T(2)=invalid;
+            EXPECT_THROW_WITH(run_production_buoyant_outer_loop(carrier.mesh,carrier.ubc,carrier.pbc,
+                carrier.controls,setup,carrier.U,carrier.p,{},&loaded_T),
+                std::invalid_argument,"initial temperature field must be finite and positive");
+        }
+        Field<double,Location::CELL> short_T(2,"T","K",1); short_T.fill(300.0);
+        Field<double,Location::CELL> vector_T(3,"T","K",2); vector_T.fill(300.0);
+        for (const auto* invalid : {&short_T,&vector_T})
+            EXPECT_THROW_WITH(run_production_buoyant_outer_loop(carrier.mesh,carrier.ubc,carrier.pbc,
+                carrier.controls,setup,carrier.U,carrier.p,{},invalid),
+                std::invalid_argument,"initial temperature field must be a scalar matching the mesh");
+        EXPECT_TRUE(predictor_calls==0);
+        for (std::size_t cell=0;cell<3;++cell) {
+            for (std::size_t d=0;d<3;++d) EXPECT_NEAR(carrier.U(cell,d),-7.0,0.0);
+            EXPECT_NEAR(carrier.p(cell),4.0,0.0);
+        }
+    });
+#endif
+    return run_all();
+}
diff --git a/tests/validation/test_level_b_reference_benchmarks.cpp b/tests/validation/test_level_b_reference_benchmarks.cpp
index c79bd903..4dd47041 100644
--- a/tests/validation/test_level_b_reference_benchmarks.cpp
+++ b/tests/validation/test_level_b_reference_benchmarks.cpp
@@ -16,7 +16,7 @@ double flat_plate_turbulent_cf(double Re_x){if(Re_x<=0)throw std::invalid_argume
 int main(){try{
  const double k=0.81,eps=0.27;
  if(std::abs(turbulent_kinematic_viscosity_kepsilon(k,eps)-0.09*k*k/eps)>1e-14)throw std::runtime_error("k-epsilon closure");
- if(std::abs(turbulent_kinematic_viscosity_komega_sst(k,2.0,0.4)-0.31*k/std::max(0.31*2.0,2.0*0.4))>1e-14)throw std::runtime_error("SST closure");
+ if(std::abs(turbulent_kinematic_viscosity_komega_sst(k,2.0,0.4,2.0,0.31)-0.313875)>1e-14)throw std::runtime_error("SST closure reference value");
  const Vec3 s{2,-1,0.5}; if(!(smagorinsky_eddy_viscosity(0.2,strain_rate_magnitude(s))>0))throw std::runtime_error("Smagorinsky");
  if(!(des_eddy_viscosity(0.2,0.01,1.0)>0))throw std::runtime_error("DES");
  for(const auto&p:channel_dns)if(!std::isfinite(p.y_plus)||!std::isfinite(p.u_plus))throw std::runtime_error("DNS reference");
diff --git a/tests/validation/test_level_c_coupled_verification.cpp b/tests/validation/test_level_c_coupled_verification.cpp
index 0f040cce..2979aa89 100644
--- a/tests/validation/test_level_c_coupled_verification.cpp
+++ b/tests/validation/test_level_c_coupled_verification.cpp
@@ -168,6 +168,12 @@ int main()
         {
             Mesh m1=cube_with_patches();
             Mesh m2=cube_with_patches();
+            // Cold region lies above the hot one, not on top of its volume.
+            // Its lower y face is the common interface.
+            for (std::size_t i=0;i<m2.n_points();++i)
+                m2.points().set(i,m2.points().x(i),m2.points().y(i)+1.0,m2.points().z(i));
+            m2.boundary().patch(2).name="interface";
+            m2.boundary().patch(3).name="y_max";
             auto g1=build_fv_geometry(m1);
             auto g2=build_fv_geometry(m2);
 
@@ -195,6 +201,8 @@ int main()
 
             auto bc1=fixed_wall_conditions(500.0,"interface");
             auto bc2=fixed_wall_conditions(300.0,"interface");
+            bc2.erase("y_min");
+            bc2["y_max"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
 
             const auto r=solve_two_region_cht(
                 m1,g1,m2,g2,f1,f2,T1,T2,s1,s2,ec,ec,cc,bc1,bc2);
diff --git a/tests/validation/test_phase10_sst_limiter.cpp b/tests/validation/test_phase10_sst_limiter.cpp
index 0a15c625..d79c1b6d 100644
--- a/tests/validation/test_phase10_sst_limiter.cpp
+++ b/tests/validation/test_phase10_sst_limiter.cpp
@@ -9,7 +9,9 @@ int main() {
         const double k=0.1, omega=10.0, F2=0.8, strain=5.0, a1=0.31;
         const double nut=cfdx::physics::turbulent_kinematic_viscosity_komega_sst(
             k,omega,F2,strain,a1);
-        const double expected=a1*k/std::max(a1*omega,strain*F2);
+        // Independent numeric oracle for these inputs: max(3.1, 4.0) = 4.0,
+        // so nu_t = 0.31*0.1/4.0 = 0.00775. The old omega*F2 path gives 0.003875.
+        constexpr double expected = 0.00775;
         if(std::abs(nut-expected)>1e-14) throw std::runtime_error("SST limiter mismatch");
         std::cout<<"PHASE10_SST_LIMITER: PASS\n";
         return 0;
diff --git a/tests/validation/test_production_s2s_energy_outer_loop_vv.cpp b/tests/validation/test_production_s2s_energy_outer_loop_vv.cpp
index d2c1c453..395c6919 100644
--- a/tests/validation/test_production_s2s_energy_outer_loop_vv.cpp
+++ b/tests/validation/test_production_s2s_energy_outer_loop_vv.cpp
@@ -246,12 +246,12 @@ static CaseSetup make_s2s_setup()
     return setup;
 }
 
-// Independent |Sf|-weighted wall temperature of a patch: the only link
-// between the resolved thermal field and the radiosity system, recomputed
-// here without any production helper.
+// Independent boundary-law reconstruction, without the production helper.
 static double patch_wall_temperature(
     const Mesh& mesh, const FvGeometry& geometry, const std::string& name,
-    const Field<double, Location::CELL>& temperature)
+    const Field<double, Location::CELL>& temperature,
+    double heat_flux_into_domain,
+    const std::vector<double>* conductivity_values = nullptr)
 {
     double weighted = 0.0;
     double total = 0.0;
@@ -265,7 +265,14 @@ static double patch_wall_temperature(
                 area.x * area.x + area.y * area.y + area.z * area.z);
             const auto owner = static_cast<std::size_t>(
                 mesh.ownership().owner(face));
-            weighted += weight * temperature(owner);
+            const auto delta = geometry.face_centres[face] -
+                               geometry.cell_centres[owner];
+            const double normal_distance = std::abs(
+                delta.x * area.x + delta.y * area.y + delta.z * area.z) / weight;
+            const double conductivity = conductivity_values
+                ? conductivity_values->at(owner) : 1.0;
+            weighted += weight * (temperature(owner) +
+                heat_flux_into_domain * normal_distance / conductivity);
             total += weight;
         }
     }
@@ -283,6 +290,30 @@ static std::size_t surface_index(const ProductionS2SRadiationResult& result,
 
 int main()
 {
+    run_case("production_s2s_outer_loop_isothermal_zero_applied_map", [] {
+        auto mesh = box_mesh(2, 2, 0.25, 0.25);
+        const auto geometry = build_fv_geometry(mesh);
+        auto setup = make_s2s_setup();
+        for (auto& spec : setup.boundary_conditions) {
+            if (spec.patch_name == "left" || spec.patch_name == "right")
+                spec.wall_temperature = 350.0;
+            if (spec.patch_name == "bottom") spec.raw_params["heat_flux"] = "0";
+        }
+        Field<double, Location::FACE> mass_flux(mesh.n_faces(), "mass_flux", "kg/s", 1);
+        mass_flux.fill(0.0);
+        Field<double, Location::CELL> temperature(mesh.n_cells(), "T", "K", 1);
+        temperature.fill(350.0);
+        const auto result = run_production_s2s_energy_outer_loop(
+            mesh, geometry, setup, mass_flux, temperature);
+        EXPECT_TRUE(result.converged);
+        EXPECT_TRUE(result.outer_iterations == 1);
+        EXPECT_TRUE(result.radiation.applied_heat_flux.size() == 2);
+        for (const auto& name : result.radiation.surface_names)
+            EXPECT_NEAR(result.radiation.applied_heat_flux.at(name), 0.0, 0.0);
+        for (const double flux : result.radiation.net_flux)
+            EXPECT_NEAR(flux, 0.0, 1e-10);
+    });
+
     // -----------------------------------------------------------------
     // Happy path: 4x4x1 box, two-surface parallel-plate enclosure, outer
     // energy-radiation iteration to the coupled fixed point.
@@ -341,9 +372,11 @@ int main()
 
         // Independent wall temperatures of the FINAL field.
         const double tb =
-            patch_wall_temperature(mesh, geometry, "bottom", temperature);
+            patch_wall_temperature(mesh, geometry, "bottom", temperature,
+                50.0 + stage.radiation.applied_heat_flux.at("bottom"));
         const double tt =
-            patch_wall_temperature(mesh, geometry, "top", temperature);
+            patch_wall_temperature(mesh, geometry, "top", temperature,
+                stage.radiation.applied_heat_flux.at("top"));
         EXPECT_NEAR(stage.radiation.wall_temperatures[bi], tb, 1e-12);
         EXPECT_NEAR(stage.radiation.wall_temperatures[ti], tt, 1e-12);
 
@@ -422,7 +455,8 @@ int main()
             mesh, geometry, setup, mass_flux, two_pass_field);
         EXPECT_TRUE(two_pass.converged);  // its own stage-2 gates pass
         const auto response = solve_production_s2s_radiation(
-            mesh, geometry, setup, two_pass_field);
+            mesh, geometry, setup, two_pass_field,
+            two_pass.radiation.injected_heat_flux);
         double two_pass_residual = 0.0;
         for (std::size_t i = 0; i < response.surface_count; ++i) {
             const std::string& name = response.surface_names[i];
@@ -451,6 +485,52 @@ int main()
         EXPECT_TRUE(stage.max_wall_temperature_lag >= 0.0);
     });
 
+    run_case("production_s2s_outer_loop_local_conductivity_replay", [] {
+        auto mesh = box_mesh(4, 4, 0.25, 0.25);
+        const auto geometry = build_fv_geometry(mesh);
+        const auto setup = make_s2s_setup();
+        Field<double, Location::FACE> mass_flux(mesh.n_faces(), "mass_flux", "kg/s", 1);
+        mass_flux.fill(0.0);
+        Field<double, Location::CELL> temperature(mesh.n_cells(), "T", "K", 1);
+        temperature.fill(350.0);
+        std::vector<double> conductivity(mesh.n_cells());
+        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
+            conductivity[c] = 2.0 + geometry.cell_centres[c].x;
+        const auto result = run_production_s2s_energy_outer_loop(
+            mesh, geometry, setup, mass_flux, temperature, &conductivity);
+        EXPECT_TRUE(result.converged);
+        EXPECT_TRUE(result.fixed_point_flux_passed);
+        EXPECT_TRUE(result.fixed_point_wall_temperature_passed);
+        EXPECT_TRUE(result.energy.global_balance_normalized <= 1e-12);
+        EXPECT_TRUE(result.energy.local_balance_normalized <= 1e-12);
+        const double tb = patch_wall_temperature(
+            mesh, geometry, "bottom", temperature,
+            50.0 + result.radiation.applied_heat_flux.at("bottom"), &conductivity);
+        const double tt = patch_wall_temperature(
+            mesh, geometry, "top", temperature,
+            result.radiation.applied_heat_flux.at("top"), &conductivity);
+        const auto bi = surface_index(result.radiation, "bottom");
+        const auto ti = surface_index(result.radiation, "top");
+        EXPECT_NEAR(result.radiation.wall_temperatures[bi], tb, 1e-12);
+        EXPECT_NEAR(result.radiation.wall_temperatures[ti], tt, 1e-12);
+        const double oracle = STEFAN_BOLTZMANN *
+            (std::pow(tb, 4) - std::pow(tt, 4)) /
+            (1.0 / 0.8 + 1.0 / 0.9 - 1.0);
+        EXPECT_NEAR(result.radiation.net_flux[bi], oracle, 1e-10);
+        EXPECT_NEAR(result.radiation.net_flux[ti], -oracle, 1e-10);
+        Field<double, Location::CELL> replay(mesh.n_cells(), "replay", "K", 1);
+        replay.fill(350.0);
+        const std::map<std::string, double> oracle_flux{
+            {"bottom", -oracle}, {"top", oracle}};
+        const auto replay_result = run_production_energy_stage(
+            mesh, geometry, setup, mass_flux, replay, oracle_flux, &conductivity);
+        EXPECT_TRUE(replay_result.converged);
+        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
+            EXPECT_NEAR(replay(c), temperature(c), 1e-6);
+        EXPECT_NEAR(result.energy.patch_flux_out.at("left") +
+                    result.energy.patch_flux_out.at("right"), 50.0, 1e-6);
+    });
+
     // -----------------------------------------------------------------
     // Fail-closed catalogue on a 2x2x1 box.
     // -----------------------------------------------------------------
diff --git a/tests/validation/test_production_s2s_energy_stage_vv.cpp b/tests/validation/test_production_s2s_energy_stage_vv.cpp
index bfecb2c0..63de3d41 100644
--- a/tests/validation/test_production_s2s_energy_stage_vv.cpp
+++ b/tests/validation/test_production_s2s_energy_stage_vv.cpp
@@ -10,8 +10,8 @@
 //     hand with Cramer's rule, and whose net flux is checked against the
 //     parallel-plate identity q = sigma (T_b^4 - T_t^4) / (1/eps_b +
 //     1/eps_t - 1). Wall temperatures are recomputed INDEPENDENTLY by the
-//     test (|Sf|-weighted average of adjacent cell temperatures over the
-//     patch faces) from the pass-1 production energy field.
+//     test (|Sf|-weighted boundary temperature, extrapolated from owner
+//     cells with the imposed thermal gradient) from the pass-1 energy field.
 //   * Blocking machine gates: the global enclosure balance closes
 //     (balance_passed, balance_normalized <= 1e-10) and the per-surface
 //     identity holds (surface_identity_passed).
@@ -260,12 +260,12 @@ static CaseSetup make_s2s_setup()
     return setup;
 }
 
-// Independent |Sf|-weighted wall temperature of a patch: the only link
-// between the resolved thermal field and the radiosity system, recomputed
-// here without any production helper.
+// Independent boundary-law reconstruction, without the production helper.
 static double patch_wall_temperature(
     const Mesh& mesh, const FvGeometry& geometry, const std::string& name,
-    const Field<double, Location::CELL>& temperature)
+    const Field<double, Location::CELL>& temperature,
+    double heat_flux_into_domain = 0.0,
+    const std::vector<double>* conductivity_values = nullptr)
 {
     double weighted = 0.0;
     double total = 0.0;
@@ -279,7 +279,14 @@ static double patch_wall_temperature(
                 area.x * area.x + area.y * area.y + area.z * area.z);
             const auto owner = static_cast<std::size_t>(
                 mesh.ownership().owner(face));
-            weighted += weight * temperature(owner);
+            const auto delta = geometry.face_centres[face] -
+                               geometry.cell_centres[owner];
+            const double normal_distance = std::abs(
+                delta.x * area.x + delta.y * area.y + delta.z * area.z) / weight;
+            const double conductivity = conductivity_values
+                ? conductivity_values->at(owner) : 1.0;
+            weighted += weight * (temperature(owner) +
+                heat_flux_into_domain * normal_distance / conductivity);
             total += weight;
         }
     }
@@ -297,6 +304,129 @@ static std::size_t surface_index(const ProductionS2SRadiationResult& result,
 
 int main()
 {
+    run_case("production_s2s_wall_temperature_obeys_thermal_gradient", [] {
+        auto mesh = box_mesh(2, 2, 0.25, 0.25);
+        auto geometry = build_fv_geometry(mesh);
+        const auto setup = make_s2s_setup();
+        Field<double, Location::CELL> temperature(mesh.n_cells(), "T", "K", 1);
+        temperature.fill(350.0);
+        // Tangential offset must not increase the normal extrapolation length.
+        for (std::size_t p = 0; p < mesh.boundary().n_patches(); ++p) {
+            const auto& patch = mesh.boundary().patch(p);
+            if (patch.name == "bottom")
+                for (const auto face : patch.face_ids)
+                    geometry.face_centres[face].x += 0.2;
+        }
+        const auto result = solve_production_s2s_radiation(
+            mesh, geometry, setup, temperature);
+        const auto bi = surface_index(result, "bottom");
+        const auto ti = surface_index(result, "top");
+        // Fourier boundary law: q_into = lambda*dT/dn, d_normal=dy/2.
+        const double expected_bottom = 350.0 + 50.0 * 0.125 / 1.0;
+        EXPECT_NEAR(result.wall_temperatures[bi], expected_bottom, 1e-12);
+        EXPECT_NEAR(result.wall_temperatures[ti], 350.0, 1e-12);
+        EXPECT_NEAR(result.applied_heat_flux.at("bottom"), 0.0, 0.0);
+        EXPECT_NEAR(result.applied_heat_flux.at("top"), 0.0, 0.0);
+        const double expected_flux = STEFAN_BOLTZMANN *
+            (std::pow(expected_bottom, 4) - std::pow(350.0, 4)) /
+            (1.0 / 0.8 + 1.0 / 0.9 - 1.0);
+        EXPECT_NEAR(result.net_flux[bi], expected_flux, 1e-10);
+        EXPECT_NEAR(result.net_flux[ti], -expected_flux, 1e-10);
+        EXPECT_TRUE(result.converged);
+    });
+
+    run_case("production_s2s_rejects_nonphysical_reconstructed_wall_temperature", [] {
+        auto mesh = box_mesh(2, 2, 0.25, 0.25);
+        const auto geometry = build_fv_geometry(mesh);
+        auto setup = make_s2s_setup();
+        for (auto& spec : setup.boundary_conditions)
+            if (spec.patch_name == "bottom")
+                spec.raw_params["heat_flux"] = "-10000";
+        Field<double, Location::CELL> temperature(mesh.n_cells(), "T", "K", 1);
+        temperature.fill(350.0);
+        EXPECT_THROW(solve_production_s2s_radiation(
+                         mesh, geometry, setup, temperature),
+                     std::invalid_argument);
+        // A positive patch mean must not hide one negative reconstructed face.
+        for (auto& spec : setup.boundary_conditions)
+            if (spec.patch_name == "bottom")
+                spec.raw_params["heat_flux"] = "-1000";
+        for (std::size_t p = 0; p < mesh.boundary().n_patches(); ++p) {
+            const auto& patch = mesh.boundary().patch(p);
+            if (patch.name != "bottom") continue;
+            temperature(mesh.ownership().owner(patch.face_ids[0])) = 100.0;
+            temperature(mesh.ownership().owner(patch.face_ids[1])) = 2000.0;
+        }
+        EXPECT_THROW(solve_production_s2s_radiation(
+                         mesh, geometry, setup, temperature),
+                     std::invalid_argument);
+    });
+
+    run_case("production_s2s_wall_temperature_uses_applied_flux_and_local_conductivity", [] {
+        auto mesh = box_mesh(2, 2, 0.25, 0.25);
+        const auto geometry = build_fv_geometry(mesh);
+        const auto setup = make_s2s_setup();
+        Field<double, Location::CELL> temperature(mesh.n_cells(), "T", "K", 1);
+        temperature.fill(350.0);
+        std::vector<double> conductivity(mesh.n_cells());
+        for (std::size_t c = 0; c < mesh.n_cells(); ++c)
+            conductivity[c] = 2.0 + static_cast<double>(c);
+        const std::map<std::string, double> applied{{"bottom", -10.0}, {"top", 10.0}};
+        const auto result = solve_production_s2s_radiation(
+            mesh, geometry, setup, temperature, applied, &conductivity);
+        const auto bi = surface_index(result, "bottom");
+        const auto ti = surface_index(result, "top");
+        const double tb = patch_wall_temperature(
+            mesh, geometry, "bottom", temperature, 40.0, &conductivity);
+        const double tt = patch_wall_temperature(
+            mesh, geometry, "top", temperature, 10.0, &conductivity);
+        EXPECT_NEAR(result.wall_temperatures[bi], tb, 1e-12);
+        EXPECT_NEAR(result.wall_temperatures[ti], tt, 1e-12);
+        EXPECT_NEAR(result.applied_heat_flux.at("bottom"), -10.0, 0.0);
+        const double oracle = STEFAN_BOLTZMANN *
+            (std::pow(tb, 4) - std::pow(tt, 4)) /
+            (1.0 / 0.8 + 1.0 / 0.9 - 1.0);
+        EXPECT_NEAR(result.net_flux[bi], oracle, 1e-10);
+        EXPECT_NEAR(result.net_flux[ti], -oracle, 1e-10);
+        EXPECT_TRUE(result.converged);
+        conductivity.pop_back();
+        EXPECT_THROW(solve_production_s2s_radiation(
+                         mesh, geometry, setup, temperature, applied, &conductivity),
+                     std::invalid_argument);
+    });
+
+    run_case("production_s2s_second_energy_pass_preserves_temperature_rejection", [] {
+        auto mesh = box_mesh(4, 4, 0.25, 0.25);
+        const auto geometry = build_fv_geometry(mesh);
+        auto setup = make_s2s_setup();
+        setup.materials.front().thermal_conductivity = 0.01;
+        Field<double, Location::FACE> mass_flux(mesh.n_faces(), "mass_flux", "kg/s", 1);
+        mass_flux.fill(0.0);
+        Field<double, Location::CELL> temperature(mesh.n_cells(), "T", "K", 1);
+        temperature.fill(350.0);
+
+        // A large one-pass radiative feedback can drive the second linear
+        // energy solution outside the kelvin domain despite perfect balances.
+        // Return that structured verdict; do not attempt a wall-temperature
+        // diagnostic on the rejected field and replace it with an exception.
+        const auto result = run_production_s2s_energy_stage(
+            mesh, geometry, setup, mass_flux, temperature);
+        EXPECT_TRUE(result.energy_pass1_converged);
+        EXPECT_TRUE(result.s2s_radiation_converged);
+        EXPECT_FALSE(result.energy_pass2_converged);
+        EXPECT_FALSE(result.converged);
+        EXPECT_TRUE(result.energy.linear_solver_converged);
+        EXPECT_TRUE(result.energy.global_energy_balance_passed);
+        EXPECT_TRUE(result.energy.local_energy_balance_passed);
+        EXPECT_FALSE(result.energy.temperature_admissibility_passed);
+        EXPECT_TRUE(result.energy.nonpositive_temperature > 0);
+        EXPECT_TRUE(result.energy.min_temperature < 0.0);
+        bool contains_negative = false;
+        for (std::size_t cell = 0; cell < mesh.n_cells(); ++cell)
+            contains_negative = contains_negative || temperature(cell) < 0.0;
+        EXPECT_TRUE(contains_negative); // evidence is retained, never clipped
+    });
+
     // -----------------------------------------------------------------
     // Happy path: 4x4x1 box, two-surface parallel-plate enclosure.
     // -----------------------------------------------------------------
@@ -323,7 +453,7 @@ int main()
 
         // Independent wall temperatures from the pass-1 field.
         const double tb0 =
-            patch_wall_temperature(mesh, geometry, "bottom", t0);
+            patch_wall_temperature(mesh, geometry, "bottom", t0, 50.0);
         const double tt0 =
             patch_wall_temperature(mesh, geometry, "top", t0);
         EXPECT_TRUE(tb0 > 0.0 && tt0 > 0.0);
diff --git a/tests/validation/test_radiation_vv.cpp b/tests/validation/test_radiation_vv.cpp
index 85c7b463..88469c16 100644
--- a/tests/validation/test_radiation_vv.cpp
+++ b/tests/validation/test_radiation_vv.cpp
@@ -1,9 +1,11 @@
 #include "cfdx/physics/radiation.h"
 #include "cfdx/physics/radiation_solver.h"
+#include "cfdx/physics/m1_m4_models.h"
 #include "common/test_harness.h"
 
 #include <cmath>
 #include <iostream>
+#include <limits>
 #include <string>
 
 using namespace cfdx::core;
@@ -93,5 +95,38 @@ int main()
             expected,1e-12*expected);
     });
 
+    run_case("radiation_p1_helper_rejects_nonfinite_inputs", [] {
+        const double invalid[] = {
+            std::numeric_limits<double>::quiet_NaN(),
+            std::numeric_limits<double>::infinity(),
+            -std::numeric_limits<double>::infinity(), -1.0};
+        for (double value : invalid) {
+            EXPECT_THROW(m1m4::p1_blackbody_intensity(value),
+                         std::invalid_argument);
+            EXPECT_THROW(m1m4::p1_source(value, 0.0, 800.0),
+                         std::invalid_argument);
+            EXPECT_THROW(m1m4::p1_source(1.0, value, 800.0),
+                         std::invalid_argument);
+            EXPECT_THROW(m1m4::p1_source(1.0, 0.0, value),
+                         std::invalid_argument);
+        }
+        // Valid transparent/cold limits and heating must not be rejected.
+        EXPECT_NEAR(m1m4::p1_blackbody_intensity(0.0), 0.0, 0.0);
+        EXPECT_NEAR(m1m4::p1_source(0.0, 12.0, 800.0), 0.0, 0.0);
+        EXPECT_NEAR(m1m4::p1_source(2.0, 3.0, 0.0),
+                    -24.0 * M_PI, 1e-12);
+    });
+
+    run_case("radiation_p1_helper_rejects_unrepresentable_output", [] {
+        const double largest = std::numeric_limits<double>::max();
+        EXPECT_THROW(m1m4::p1_blackbody_intensity(largest),
+                     std::overflow_error);
+        EXPECT_THROW(m1m4::p1_source(largest, 0.0, 800.0),
+                     std::overflow_error);
+        // A large finite absorption must still yield zero at equilibrium.
+        const double Ib = m1m4::p1_blackbody_intensity(800.0);
+        EXPECT_NEAR(m1m4::p1_source(largest, Ib, 800.0), 0.0, 0.0);
+    });
+
     return run_all();
 }
diff --git a/tests/validation/test_s2s_radiation_vv.cpp b/tests/validation/test_s2s_radiation_vv.cpp
index 0ea2d66f..87b1abd1 100644
--- a/tests/validation/test_s2s_radiation_vv.cpp
+++ b/tests/validation/test_s2s_radiation_vv.cpp
@@ -3,6 +3,7 @@
 
 #include <cmath>
 #include <iostream>
+#include <limits>
 #include <vector>
 
 using namespace cfdx::physics;
@@ -86,6 +87,100 @@ int main()
         EXPECT_NEAR(A[0]*F[1],A[1]*F[2],1e-12);
     });
 
+    run_case("s2s_centroid_geometry_rejects_nonfinite_inputs", [] {
+        const ViewFactorPatch source{{0, 0, 0}, {0, 0, 1}, 2.0};
+        const ViewFactorPatch target{{0, 0, 2}, {0, 0, -1}, 1.0};
+        const double invalid[] = {
+            std::numeric_limits<double>::quiet_NaN(),
+            std::numeric_limits<double>::infinity(),
+            -std::numeric_limits<double>::infinity()};
+        for (double value : invalid) {
+            for (std::size_t axis = 0; axis < 3; ++axis) {
+                auto bad = source;
+                bad.center[axis] = value;
+                EXPECT_THROW(patch_pair_view_factor(bad, target),
+                             std::invalid_argument);
+                EXPECT_THROW(patch_pair_view_factor(target, bad),
+                             std::invalid_argument);
+                EXPECT_THROW(estimate_view_factor_matrix({bad}),
+                             std::invalid_argument);
+                bad = source;
+                bad.normal[axis] = value;
+                EXPECT_THROW(patch_pair_view_factor(bad, target),
+                             std::invalid_argument);
+                EXPECT_THROW(patch_pair_view_factor(target, bad),
+                             std::invalid_argument);
+                EXPECT_THROW(estimate_view_factor_matrix({bad}),
+                             std::invalid_argument);
+            }
+            auto bad = source;
+            bad.area = value;
+            EXPECT_THROW(patch_pair_view_factor(bad, target),
+                         std::invalid_argument);
+            EXPECT_THROW(patch_pair_view_factor(target, bad),
+                         std::invalid_argument);
+            EXPECT_THROW(estimate_view_factor_matrix({bad}),
+                         std::invalid_argument);
+        }
+        auto bad = source;
+        bad.normal = {0, 0, 0};
+        EXPECT_THROW(estimate_view_factor_matrix({bad}),
+                     std::invalid_argument);
+        // Reject geometry before even the coincident-centroid early return.
+        EXPECT_THROW(patch_pair_view_factor(bad, source),
+                     std::invalid_argument);
+        for (double area : {0.0, -1.0}) {
+            bad = source;
+            bad.area = area;
+            EXPECT_THROW(patch_pair_view_factor(bad, source),
+                         std::invalid_argument);
+            EXPECT_THROW(estimate_view_factor_matrix({bad}),
+                         std::invalid_argument);
+        }
+        auto distant = target;
+        distant.center[2] = std::numeric_limits<double>::max();
+        auto opposite = source;
+        opposite.center[2] = -std::numeric_limits<double>::max();
+        EXPECT_THROW(patch_pair_view_factor(opposite, distant),
+                     std::overflow_error);
+        // Finite geometry may still exceed the representable exchange range.
+        // Reject it before a zero closure scale could multiply infinity.
+        auto tiny = source;
+        tiny.area = 1e-300;
+        auto huge = target;
+        huge.area = 1e308;
+        huge.center[2] = 0.1;
+        EXPECT_THROW(estimate_view_factor_matrix({tiny, huge}),
+                     std::overflow_error);
+    });
+
+    run_case("s2s_centroid_geometry_preserves_reciprocity_and_limits", [] {
+        const ViewFactorPatch source{{0, 0, 0}, {0, 0, 2}, 2.0};
+        const ViewFactorPatch target{{0, 0, 2}, {0, 0, -3}, 1.0};
+        const auto F = estimate_view_factor_matrix({source, target});
+        const std::vector<double> areas{source.area, target.area};
+        validate_s2s_view_factors(F, areas, 1e-12);
+        EXPECT_NEAR(F[1], 1.0 / (4.0 * M_PI), 1e-14);
+        EXPECT_NEAR(F[2], 2.0 / (4.0 * M_PI), 1e-14);
+        EXPECT_NEAR(areas[0] * F[1], areas[1] * F[2], 1e-14);
+        EXPECT_NEAR(F[0], 0.0, 0.0);
+        EXPECT_NEAR(F[3], 0.0, 0.0);
+        EXPECT_NEAR(patch_pair_view_factor(source, source), 0.0, 0.0);
+        auto facing_away = target;
+        facing_away.normal = {0, 0, 1};
+        EXPECT_NEAR(patch_pair_view_factor(source, facing_away), 0.0, 0.0);
+
+        // Scaling normals cannot change the projected solid angle.
+        auto large_normal = source;
+        large_normal.normal[2] = 1e200;
+        EXPECT_NEAR(patch_pair_view_factor(large_normal, target), F[1], 1e-14);
+        const auto enclosure = solve_s2s_radiosity(
+            {1.0, 1.0}, {0.7, 0.5}, {800.0, 400.0},
+            {0.0, 1.0, 1.0, 0.0});
+        EXPECT_TRUE(enclosure.converged);
+        EXPECT_NEAR(enclosure.total_power, 0.0, 1e-9);
+    });
+
     run_case("s2s_monte_carlo_visibility_and_blocking", [] {
         const std::vector<RadiationTriangle> s1{{
             {0,0,0},{1,0,0},{0,1,0}}}; // +z
@@ -116,6 +211,42 @@ int main()
         EXPECT_NEAR(areas[0] * matrix[1], areas[1] * matrix[2], 1e-12);
     });
 
+    run_case("s2s_triangle_kernel_rejects_nonfinite_geometry", [] {
+        const RadiationTriangle triangle{{0, 0, 1}, {1, 0, 1}, {0, 1, 1}};
+        const double invalid[] = {
+            std::numeric_limits<double>::quiet_NaN(),
+            std::numeric_limits<double>::infinity(),
+            -std::numeric_limits<double>::infinity()};
+        for (double value : invalid) {
+            for (std::size_t axis = 0; axis < 3; ++axis) {
+                for (std::size_t vertex = 0; vertex < 3; ++vertex) {
+                    auto bad = triangle;
+                    auto& coordinate = vertex == 0 ? bad.a : vertex == 1 ? bad.b : bad.c;
+                    coordinate[axis] = value;
+                    double distance = 0.0;
+                    EXPECT_THROW(radiation_ray_triangle_hit(
+                        {0.1, 0.1, 0}, {0, 0, 1}, bad, distance), std::invalid_argument);
+                    EXPECT_THROW(radiation_triangle_area(bad), std::invalid_argument);
+                    EXPECT_THROW(radiation_triangle_normal(bad), std::invalid_argument);
+                    EXPECT_THROW(estimate_view_factor_ray_traced(
+                        {bad}, {triangle}, {}, 1), std::invalid_argument);
+                }
+                std::array<double, 3> bad_ray{0, 0, 1};
+                bad_ray[axis] = value;
+                double distance = 0.0;
+                EXPECT_THROW(radiation_ray_triangle_hit(
+                    bad_ray, {0, 0, 1}, triangle, distance), std::invalid_argument);
+                EXPECT_THROW(radiation_ray_triangle_hit(
+                    {0.1, 0.1, 0}, bad_ray, triangle, distance), std::invalid_argument);
+            }
+        }
+        double distance = 0.0;
+        EXPECT_TRUE(radiation_ray_triangle_hit(
+            {0.1, 0.1, 0}, {0, 0, 1}, triangle, distance));
+        EXPECT_NEAR(distance, 1.0, 1e-14);
+        EXPECT_NEAR(radiation_triangle_area(triangle), 0.5, 1e-14);
+    });
+
     run_case("s2s_model_selector", [] {
         RadiationModelSelector selector;
         EXPECT_TRUE(selector.select_surface_to_surface(false)==RadiationApproximation::S2S);
diff --git a/tests/validation/test_thermal_challenge_vv.cpp b/tests/validation/test_thermal_challenge_vv.cpp
index c0ddc478..54075352 100644
--- a/tests/validation/test_thermal_challenge_vv.cpp
+++ b/tests/validation/test_thermal_challenge_vv.cpp
@@ -326,7 +326,7 @@ int main()
         T1(0)=400; T2(0)=300;
         Field<double,Location::CELL> s1(1,"s1","W/m3",1), s2(1,"s2","W/m3",1); s1.fill(0); s2.fill(0);
         EnergySolverControls e1; e1.conductivity=2; e1.relaxation=1; e1.max_iterations=20; e1.tolerance=1e-10;
-        EnergySolverControls e2=e1;
+        EnergySolverControls e2=e1; e2.conductivity=1;
         ChtInterfaceControls c; c.region1_patch="interface"; c.region2_patch="interface";
         c.conductivity1=2; c.conductivity2=1; c.tolerance=1e-10; c.temperature_tolerance=1e-10;
         ScalarBoundaryConditions b1{{"outer_hot",{ScalarBoundaryType::FIXED_VALUE,400,0}},
diff --git a/tests/validation/test_thermal_radiation_regression.cpp b/tests/validation/test_thermal_radiation_regression.cpp
index 0f740ba1..1c81a013 100644
--- a/tests/validation/test_thermal_radiation_regression.cpp
+++ b/tests/validation/test_thermal_radiation_regression.cpp
@@ -168,10 +168,12 @@ int main()
             {{0,0,1},{0,0,-1},1.0}
         };
         const auto F=estimate_view_factor_matrix(patches);
-        validate_view_factor_matrix(F,2,{2.0,1.0},1e-12);
+        validate_s2s_view_factors(F,{2.0,1.0},1e-12);
+        EXPECT_NEAR(F[0],0.0,1e-14);
+        EXPECT_NEAR(F[3],0.0,1e-14);
         EXPECT_NEAR(2.0*F[1],1.0*F[2],1e-12);
-        EXPECT_NEAR(F[0]+F[1],1.0,1e-12);
-        EXPECT_NEAR(F[2]+F[3],1.0,1e-12);
+        EXPECT_TRUE(F[0]+F[1] <= 1.0+1e-12);
+        EXPECT_TRUE(F[2]+F[3] <= 1.0+1e-12);
 
         auto degenerate=patches;
         degenerate[1].area=0.0;
diff --git a/tests/validation/test_thermal_zones_vv.cpp b/tests/validation/test_thermal_zones_vv.cpp
index c4a3307e..624a72a2 100644
--- a/tests/validation/test_thermal_zones_vv.cpp
+++ b/tests/validation/test_thermal_zones_vv.cpp
@@ -287,6 +287,28 @@ int main()
         EXPECT_TRUE(result.energy_imbalance_relative<1e-12);
     });
 
+    run_case("thermal_conduction_rejects_underrelaxed_unconverged_state", [] {
+        const std::size_t n = 2;
+        const double dx = 0.5;
+        Mesh m = chain_mesh(n, dx);
+        auto g = build_fv_geometry(m);
+        const std::vector<ThermalZoneSpec> specs = {
+            {"solid", 1, 1.0, 0.0, 0.0, true},
+        };
+        const auto zones = build_thermal_zones(m, std::vector<int>{1, 1}, specs);
+        const auto bc = chain_bcs(400.0, 300.0);
+        Vector temperature(n, 300.0);
+        const auto result = solve_thermal_conduction(
+            m, g, zones, bc, temperature, ScalarSolveControls{5000, 1e-12, 0.25});
+
+        EXPECT_TRUE(result.candidate_converged);
+        EXPECT_TRUE(!result.returned_state_converged);
+        EXPECT_TRUE(!result.converged);
+        EXPECT_NEAR(temperature(0), 318.75, 1e-12);
+        EXPECT_NEAR(temperature(1), 306.25, 1e-12);
+        EXPECT_NEAR(result.residual_inf, 300.0, 1e-10);
+    });
+
     run_case("thermal_zones_contact_resistance_exact", [] {
         const std::size_t n = 6;
         const double dx = 0.5;
Annexe B — Preuves consolidées de la v2
Commandes et codes de retour
{
  "release": {
    "variant": "integration_v2",
    "started_at_utc": "2026-10-24T17:06:00.168174+00:00",
    "base_sha": "038e3aeff40f02af95b8845b931245ce23a1a094",
    "stages": [
      {
        "stage": "build",
        "command": [
          "/workspace/scratch/cae91c5c59b2/audit_tools/cmake/data/bin/cmake",
          "--build",
          "/workspace/scratch/cae91c5c59b2/build_integration_v2",
          "--parallel",
          "4",
          "--target",
          "test_m1_m4_physics",
          "test_m2_m4_solver_validation",
          "test_finite_volume_transport",
          "test_turbulence_qualification_equations",
          "test_production_turbulence_contract",
          "test_production_radiation_contract",
          "test_sst_blending",
          "test_phase10_sst_limiter",
          "test_phase10_turbulence_hardening",
          "test_steady_incompressible_solver",
          "test_cht_validation",
          "test_thermal_vv",
          "test_thermal_zones_vv",
          "test_production_energy_stage_vv",
          "test_production_s2s_energy_stage_vv",
          "test_production_s2s_energy_outer_loop_vv",
          "test_production_buoyancy_outer_loop_vv",
          "test_radiation_vv",
          "test_s2s_radiation_vv",
          "test_thermal_radiation_model_matrix",
          "test_thermal_radiation_regression",
          "test_thermal_challenge_vv",
          "test_thermophysical_models_vv",
          "test_phase12_radiation_vv",
          "test_thermo_transport_contracts",
          "test_thermophysical_models",
          "test_dat_restart_identity",
          "test_case_hdf5_io",
          "test_benchmark_matrix",
          "test_level_b_reference_benchmarks",
          "test_level_c_coupled_verification",
          "test_audit_cht_time",
          "test_audit_energy_boundary_metric",
          "test_solid_conduction_selection",
          "test_energy_integration_contracts",
          "test_cht_integration_contracts",
          "test_sst_integration_gaps",
          "test_hdf5_roundtrip",
          "test_field_roundtrip",
          "cfdx_thermal_solver",
          "cfdx_production_solver"
        ],
        "exit_code": 0,
        "elapsed_seconds": 28.031921026999953
      },
      {
        "stage": "ctest",
        "command": [
          "/workspace/scratch/cae91c5c59b2/audit_tools/cmake/data/bin/ctest",
          "--test-dir",
          "/workspace/scratch/cae91c5c59b2/build_integration_v2",
          "--output-on-failure",
          "--no-tests=error",
          "--parallel",
          "3",
          "--test-output-size-passed",
          "65536",
          "--test-output-size-failed",
          "65536",
          "--timeout",
          "180",
          "-R",
          "^(test_m1_m4_physics|test_m2_m4_solver_validation|test_finite_volume_transport|test_turbulence_qualification_equations|test_production_turbulence_contract|test_production_radiation_contract|test_sst_blending|test_phase10_sst_limiter|test_phase10_turbulence_hardening|test_steady_incompressible_solver|test_cht_validation|test_thermal_vv|test_thermal_zones_vv|test_production_energy_stage_vv|test_production_s2s_energy_stage_vv|test_production_s2s_energy_outer_loop_vv|test_production_buoyancy_outer_loop_vv|test_radiation_vv|test_s2s_radiation_vv|test_thermal_radiation_model_matrix|test_thermal_radiation_regression|test_thermal_challenge_vv|test_thermophysical_models_vv|test_phase12_radiation_vv|test_thermo_transport_contracts|test_thermophysical_models|test_dat_restart_identity|test_case_hdf5_io|test_benchmark_matrix|test_level_b_reference_benchmarks|test_level_c_coupled_verification|test_audit_cht_time|test_audit_energy_boundary_metric|test_solid_conduction_selection|test_energy_integration_contracts|test_cht_integration_contracts|test_sst_integration_gaps|test_solid_conduction_contract_e2e|test_t2_campaign_verdict|test_hdf5_roundtrip|test_field_roundtrip|test_production_thermal_restart_e2e|test_imported_thermal_contracts)$",
          "--output-junit",
          "/workspace/scratch/cae91c5c59b2/audit_evidence/integration_v2/results.xml"
        ],
        "exit_code": 0,
        "elapsed_seconds": 1.132395851000183
      }
    ]
  },
  "sanitizers": {
    "variant": "integration_v2_san",
    "started_at_utc": "2026-10-24T17:07:14.226891+00:00",
    "base_sha": "038e3aeff40f02af95b8845b931245ce23a1a094",
    "stages": [
      {
        "stage": "build",
        "command": [
          "/workspace/scratch/cae91c5c59b2/audit_tools/cmake/data/bin/cmake",
          "--build",
          "/workspace/scratch/cae91c5c59b2/build_integration_v2_san",
          "--parallel",
          "2",
          "--target",
          "test_sst_integration_gaps",
          "test_energy_integration_contracts",
          "test_cht_integration_contracts",
          "test_solid_conduction_selection",
          "test_radiation_vv",
          "test_s2s_radiation_vv",
          "test_production_s2s_energy_stage_vv",
          "test_production_s2s_energy_outer_loop_vv",
          "test_case_hdf5_io"
        ],
        "exit_code": 0,
        "elapsed_seconds": 109.55191645800005
      },
      {
        "stage": "ctest",
        "command": [
          "/workspace/scratch/cae91c5c59b2/audit_tools/cmake/data/bin/ctest",
          "--test-dir",
          "/workspace/scratch/cae91c5c59b2/build_integration_v2_san",
          "--output-on-failure",
          "--no-tests=error",
          "--parallel",
          "2",
          "--test-output-size-passed",
          "65536",
          "--test-output-size-failed",
          "65536",
          "--timeout",
          "180",
          "-R",
          "^(test_sst_integration_gaps|test_energy_integration_contracts|test_cht_integration_contracts|test_solid_conduction_selection|test_radiation_vv|test_s2s_radiation_vv|test_production_s2s_energy_stage_vv|test_production_s2s_energy_outer_loop_vv|test_case_hdf5_io)$",
          "--output-junit",
          "/workspace/scratch/cae91c5c59b2/audit_evidence/integration_v2_san/results.xml"
        ],
        "exit_code": 0,
        "elapsed_seconds": 0.21730161699997552
      }
    ],
    "sanitizer_environment": {
      "ASAN_OPTIONS": "detect_leaks=1:halt_on_error=1",
      "UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1"
    }
  }
}
Journal CTest Release
Test project /workspace/scratch/cae91c5c59b2/build_integration_v2
      Start 169: test_production_buoyancy_outer_loop_vv
      Start 235: test_imported_thermal_contracts
      Start 232: test_solid_conduction_contract_e2e
 1/43 Test #232: test_solid_conduction_contract_e2e .........   Passed    0.31 sec
      Start 234: test_production_thermal_restart_e2e
 2/43 Test #235: test_imported_thermal_contracts ............   Passed    0.32 sec
      Start 233: test_t2_campaign_verdict
 3/43 Test #233: test_t2_campaign_verdict ...................   Passed    0.19 sec
      Start 166: test_production_energy_stage_vv
 4/43 Test #166: test_production_energy_stage_vv ............   Passed    0.04 sec
      Start 225: test_case_hdf5_io
 5/43 Test #225: test_case_hdf5_io ..........................   Passed    0.02 sec
      Start 223: test_hdf5_roundtrip
 6/43 Test #223: test_hdf5_roundtrip ........................   Passed    0.01 sec
      Start 174: test_thermal_challenge_vv
 7/43 Test #174: test_thermal_challenge_vv ..................   Passed    0.02 sec
      Start   1: test_audit_cht_time
 8/43 Test   #1: test_audit_cht_time ........................   Passed    0.01 sec
      Start 171: test_s2s_radiation_vv
 9/43 Test #171: test_s2s_radiation_vv ......................   Passed    0.01 sec
      Start 224: test_field_roundtrip
10/43 Test #224: test_field_roundtrip .......................   Passed    0.01 sec
      Start  98: test_steady_incompressible_solver
11/43 Test  #98: test_steady_incompressible_solver ..........   Passed    0.01 sec
      Start 162: test_level_c_coupled_verification
12/43 Test #234: test_production_thermal_restart_e2e ........   Passed    0.34 sec
      Start 168: test_production_s2s_energy_outer_loop_vv
13/43 Test #162: test_level_c_coupled_verification ..........   Passed    0.01 sec
      Start 105: test_m2_m4_solver_validation
14/43 Test #105: test_m2_m4_solver_validation ...............   Passed    0.01 sec
      Start 112: test_sst_integration_gaps
15/43 Test #168: test_production_s2s_energy_outer_loop_vv ...   Passed    0.01 sec
      Start 110: test_energy_integration_contracts
16/43 Test #110: test_energy_integration_contracts ..........   Passed    0.00 sec
      Start 213: test_sst_blending
17/43 Test #112: test_sst_integration_gaps ..................   Passed    0.01 sec
      Start 165: test_thermal_zones_vv
18/43 Test #165: test_thermal_zones_vv ......................   Passed    0.00 sec
      Start 167: test_production_s2s_energy_stage_vv
19/43 Test #213: test_sst_blending ..........................   Passed    0.01 sec
      Start  97: test_finite_volume_transport
20/43 Test #167: test_production_s2s_energy_stage_vv ........   Passed    0.00 sec
      Start  68: test_dat_restart_identity
21/43 Test  #97: test_finite_volume_transport ...............   Passed    0.00 sec
      Start  88: test_m1_m4_physics
22/43 Test  #88: test_m1_m4_physics .........................   Passed    0.00 sec
      Start 111: test_cht_integration_contracts
23/43 Test  #68: test_dat_restart_identity ..................   Passed    0.00 sec
      Start 107: test_turbulence_qualification_equations
24/43 Test #107: test_turbulence_qualification_equations ....   Passed    0.00 sec
      Start 160: test_benchmark_matrix
25/43 Test #111: test_cht_integration_contracts .............   Passed    0.00 sec
      Start 170: test_radiation_vv
26/43 Test #160: test_benchmark_matrix ......................   Passed    0.00 sec
      Start 106: test_phase10_turbulence_hardening
27/43 Test #170: test_radiation_vv ..........................   Passed    0.00 sec
      Start 164: test_thermal_vv
28/43 Test #106: test_phase10_turbulence_hardening ..........   Passed    0.00 sec
      Start 163: test_cht_validation
29/43 Test #164: test_thermal_vv ............................   Passed    0.00 sec
      Start 113: test_production_radiation_contract
30/43 Test #113: test_production_radiation_contract .........   Passed    0.00 sec
      Start 215: test_thermophysical_models
31/43 Test #163: test_cht_validation ........................   Passed    0.00 sec
      Start 173: test_thermal_radiation_regression
32/43 Test #215: test_thermophysical_models .................   Passed    0.00 sec
      Start 175: test_thermophysical_models_vv
33/43 Test #175: test_thermophysical_models_vv ..............   Passed    0.00 sec
      Start   2: test_audit_energy_boundary_metric
34/43 Test #173: test_thermal_radiation_regression ..........   Passed    0.00 sec
      Start  60: test_phase12_radiation_vv
35/43 Test  #60: test_phase12_radiation_vv ..................   Passed    0.00 sec
      Start 124: test_phase10_sst_limiter
36/43 Test   #2: test_audit_energy_boundary_metric ..........   Passed    0.00 sec
      Start 109: test_solid_conduction_selection
37/43 Test #124: test_phase10_sst_limiter ...................   Passed    0.00 sec
      Start 214: test_thermo_transport_contracts
38/43 Test #109: test_solid_conduction_selection ............   Passed    0.00 sec
      Start 108: test_production_turbulence_contract
39/43 Test #214: test_thermo_transport_contracts ............   Passed    0.00 sec
      Start 172: test_thermal_radiation_model_matrix
40/43 Test #108: test_production_turbulence_contract ........   Passed    0.00 sec
      Start 161: test_level_b_reference_benchmarks
41/43 Test #172: test_thermal_radiation_model_matrix ........   Passed    0.00 sec
42/43 Test #161: test_level_b_reference_benchmarks ..........   Passed    0.00 sec
43/43 Test #169: test_production_buoyancy_outer_loop_vv .....   Passed    1.12 sec

100% tests passed out of 43

Label Time Summary:
audit            =   0.01 sec*proc (2 tests)
integration      =   0.97 sec*proc (3 tests)
io               =   0.32 sec*proc (1 test)
phase12          =   0.00 sec*proc (1 test)
phase13          =   1.23 sec*proc (18 tests)
qualification    =   1.36 sec*proc (6 tests)
regression       =   0.01 sec*proc (2 tests)
restart          =   0.34 sec*proc (1 test)
thermal          =   1.17 sec*proc (6 tests)
validation       =   1.73 sec*proc (21 tests)

Total Test time (real) =   1.12 sec
Journal CTest ASan/UBSan
Test project /workspace/scratch/cae91c5c59b2/build_integration_v2_san
    Start 225: test_case_hdf5_io
    Start 171: test_s2s_radiation_vv
1/9 Test #171: test_s2s_radiation_vv ......................   Passed    0.06 sec
    Start 112: test_sst_integration_gaps
2/9 Test #225: test_case_hdf5_io ..........................   Passed    0.08 sec
    Start 168: test_production_s2s_energy_outer_loop_vv
3/9 Test #168: test_production_s2s_energy_outer_loop_vv ...   Passed    0.04 sec
    Start 110: test_energy_integration_contracts
4/9 Test #112: test_sst_integration_gaps ..................   Passed    0.06 sec
    Start 167: test_production_s2s_energy_stage_vv
5/9 Test #167: test_production_s2s_energy_stage_vv ........   Passed    0.03 sec
    Start 111: test_cht_integration_contracts
6/9 Test #110: test_energy_integration_contracts ..........   Passed    0.04 sec
    Start 170: test_radiation_vv
7/9 Test #111: test_cht_integration_contracts .............   Passed    0.03 sec
    Start 109: test_solid_conduction_selection
8/9 Test #170: test_radiation_vv ..........................   Passed    0.02 sec
9/9 Test #109: test_solid_conduction_selection ............   Passed    0.02 sec

100% tests passed out of 9

Label Time Summary:
phase13          =   0.16 sec*proc (4 tests)
qualification    =   0.08 sec*proc (2 tests)
validation       =   0.16 sec*proc (4 tests)

Total Test time (real) =   0.21 sec
Régressions Python existantes et quatre skips PyFluent
................................ssss.................................... [ 79%]
...................                                                      [100%]
=========================== short test summary info ============================
SKIPPED [1] tests/python/test_fluent_adapter.py:106: pyfluent not installed
SKIPPED [1] tests/python/test_fluent_adapter.py:116: pyfluent not installed
SKIPPED [1] tests/python/test_fluent_adapter.py:129: pyfluent not installed
SKIPPED [1] tests/python/test_fluent_adapter.py:154: pyfluent not installed
87 passed, 4 skipped in 3.06s
Sorties ciblées de la campagne Release finale
test_energy_integration_contracts
Running 11 test cases...
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0 backward_error=0 linear_relative=0 dT_relative=0.333333 energy_balance=0 accepted_state=1
  PASS: implicit_source_equilibrium_uses_final_temperature
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=100 backward_error=0.0909091 linear_relative=0 dT_relative=0.2 energy_balance=0.0909091 accepted_state=0
  PASS: implicit_source_underrelaxed_state_remains_out_of_balance
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0 backward_error=0 linear_relative=0 dT_relative=4 energy_balance=0 accepted_state=1
    admissibility oracle: Su=-4800 T=-100 linear=1 global=1 local=1 residual=0 converged=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0 backward_error=0 linear_relative=0 dT_relative=300 energy_balance=0 accepted_state=1
    admissibility oracle: Su=-3600 T=0 linear=1 global=1 local=1 residual=0 converged=0
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0 backward_error=0 linear_relative=0 dT_relative=299 energy_balance=0 accepted_state=1
    admissibility oracle: Su=-3588 T=1 linear=1 global=1 local=1 residual=0 converged=1
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0 backward_error=0 linear_relative=0 dT_relative=4 energy_balance=0 accepted_state=1
  PASS: production_temperature_admissibility_is_separate_from_conservation
  PASS: independent_source_audit_detects_perturbation_and_rejects_invalid_sp
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=9.09495e-13 backward_error=7.64024e-17 linear_relative=1.85287e-16 dT_relative=0.166667 energy_balance=5.73018e-17 accepted_state=1
  PASS: two_material_conduction_matches_series_resistance_and_both_audits
  PASS: turbulent_conductivity_is_laminar_limit_and_prandtl_scaled
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=4.54747e-13 backward_error=1.12339e-16 linear_relative=4.16684e-16 dT_relative=0.0513834 energy_balance=5.0935e-17 accepted_state=1
  PASS: production_flux_bc_uses_each_owner_conductivity
  PASS: variable_conductivity_rejected_before_temperature_mutation
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=9.09495e-12 backward_error=9.92088e-16 linear_relative=2.08943e-15 dT_relative=0.00707485 energy_balance=9.64606e-17 accepted_state=1
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=4.54747e-12 backward_error=4.18983e-16 linear_relative=1.08909e-15 dT_relative=0.00713393 energy_balance=9.00209e-18 accepted_state=1
    SST-energy seam: max thermal effect=0.0179786 K; balance=9.05942e-13 W
  PASS: accepted_sst_state_drives_energy_with_authoritative_flux
    restored buoyancy: first-iterate max Uy=0.104895 m/s
  PASS: restored_temperature_drives_the_first_buoyancy_source
  PASS: invalid_restored_temperature_is_rejected_before_carrier_mutation

Results: 11 passed, 0 failed
test_cht_integration_contracts
Running 8 test cases...
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0 backward_error=0 linear_relative=0 dT_relative=0 energy_balance=0 accepted_state=1
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0 backward_error=0 linear_relative=0 dT_relative=0 energy_balance=0 accepted_state=1
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0 backward_error=0 linear_relative=0 dT_relative=0 energy_balance=0 accepted_state=1
THERMAL_RESIDUAL: iteration=1 linear_status=0 residual=0 backward_error=0 linear_relative=0 dT_relative=0 energy_balance=0 accepted_state=1
  PASS: cht_accepts_adjacent_opposed_interface
  PASS: cht_rejects_coincident_domains_with_parallel_outward_normals
  PASS: cht_rejects_empty_interface
  PASS: cht_rejects_nonfinite_interface_controls
  PASS: cht_rejects_different_physical_steps
  PASS: cht_rejects_mass_transport_at_impermeable_interface
  PASS: cht_rejects_malformed_second_region_before_first_solve
  PASS: cht_rejects_nonboundary_interface_face

Results: 8 passed, 0 failed
test_sst_integration_gaps
Running 9 test cases...
  PASS: sst_wall_omega_requires_every_mesh_face_before_carrier_mutation
  PASS: sst_cross_diffusion_uses_real_scalar_and_face_boundaries
    bounded-state residuals: k=0.09 omega=0.075 converged=0
  PASS: sst_residuals_and_convergence_describe_bounded_returned_fields
    nonzero-flux SST: k=[0.00916065,0.00838028,0.00770093] omega=[1.94838,1.90061,1.85891]
  PASS: sst_multicell_first_predictor_uses_refreshed_viscosity_and_nonzero_flux
CFDX nonlinear iteration 1 rejected; retry=1 alpha_u=0.35 alpha_p=0.15 reason=test rejection after SST commit with nonzero flux
  PASS: sst_multicell_retry_restores_every_field_and_authoritative_flux
  PASS: sst_multicell_dat_rehydrates_primary_fields_before_first_predictor
    transpose-stress MMS N=16 L2=0.00437466 omitted=0.138032 order=0
    transpose-stress MMS N=32 L2=0.00111903 omitted=0.138014 order=1.96692
    transpose-stress MMS N=64 L2=0.000282884 omitted=0.138013 order=1.98397
  PASS: variable_viscosity_transpose_stress_mms_converges_and_is_conservative
    transverse predictor Uy=0.000596783,-0.000619322
  PASS: variable_viscosity_production_predictor_consumes_transpose_stress
  PASS: momentum_balance_uses_physical_source_sign_and_density_in_laminar_limit

Results: 9 passed, 0 failed
test_solid_conduction_contract_e2e
test_absent_conductivity_is_not_replaced_by_air (__main__.SolidConductionApplicationContract.test_absent_conductivity_is_not_replaced_by_air) ... ok
test_canonical_temporal_selection_is_not_silently_ignored (__main__.SolidConductionApplicationContract.test_canonical_temporal_selection_is_not_silently_ignored) ... ok
test_fluid_inlet_is_not_replaced_by_adiabatic_wall (__main__.SolidConductionApplicationContract.test_fluid_inlet_is_not_replaced_by_adiabatic_wall) ... ok
test_multiphase_is_not_silently_ignored (__main__.SolidConductionApplicationContract.test_multiphase_is_not_silently_ignored) ... ok
test_radiation_is_not_silently_ignored (__main__.SolidConductionApplicationContract.test_radiation_is_not_silently_ignored) ... ok
test_supported_solid_case_solves (__main__.SolidConductionApplicationContract.test_supported_solid_case_solves) ... ok
test_transient_evolution_is_not_silently_ignored (__main__.SolidConductionApplicationContract.test_transient_evolution_is_not_silently_ignored) ... ok
test_turbulence_is_not_silently_ignored (__main__.SolidConductionApplicationContract.test_turbulence_is_not_silently_ignored) ... ok

----------------------------------------------------------------------
Ran 8 tests in 0.116s

OK
test_t2_campaign_verdict
test_complete_campaign_meets_physical_targets (__main__.CampaignVerdictTests.test_complete_campaign_meets_physical_targets) ... ok
test_execution_success_does_not_qualify_wrong_nusselt (__main__.CampaignVerdictTests.test_execution_success_does_not_qualify_wrong_nusselt) ... ok
test_missing_or_failed_mesh_does_not_qualify (__main__.CampaignVerdictTests.test_missing_or_failed_mesh_does_not_qualify) ... ok
test_nan_infinity_negative_and_missing_metrics_cannot_pass (__main__.CampaignVerdictTests.test_nan_infinity_negative_and_missing_metrics_cannot_pass) ... ok
test_original_thresholds_are_inclusive_and_unchanged (__main__.CampaignVerdictTests.test_original_thresholds_are_inclusive_and_unchanged) ... ok
test_poiseuille_gate_applies_to_every_mesh (__main__.CampaignVerdictTests.test_poiseuille_gate_applies_to_every_mesh) ... ok

----------------------------------------------------------------------
Ran 6 tests in 0.000s

OK
test_production_thermal_restart_e2e
test_checkpoint_preserves_temperature_and_warm_energy_initialization (__main__.ProductionThermalRestart.test_checkpoint_preserves_temperature_and_warm_energy_initialization) ... ok
test_conservative_negative_temperature_is_not_physically_accepted (__main__.ProductionThermalRestart.test_conservative_negative_temperature_is_not_physically_accepted) ... ok
test_energy_exception_cannot_report_successful_carrier_as_overall_convergence (__main__.ProductionThermalRestart.test_energy_exception_cannot_report_successful_carrier_as_overall_convergence) ... ok
test_negative_restart_temperature_is_rejected_before_flow (__main__.ProductionThermalRestart.test_negative_restart_temperature_is_rejected_before_flow) ... ok
test_physical_transient_flag_cannot_be_solved_as_steady (__main__.ProductionThermalRestart.test_physical_transient_flag_cannot_be_solved_as_steady) ... ok
test_thermal_restart_requires_temperature_instead_of_cold_fallback (__main__.ProductionThermalRestart.test_thermal_restart_requires_temperature_instead_of_cold_fallback) ... ok

----------------------------------------------------------------------
Ran 6 tests in 0.116s

OK
test_imported_thermal_contracts
test_explicit_material_values_survive_serialization (__main__.ImportedThermalContracts.test_explicit_material_values_survive_serialization) ... ok
test_fluent_gravity_uses_native_three_scalar_contract (__main__.ImportedThermalContracts.test_fluent_gravity_uses_native_three_scalar_contract) ... ok
test_unspecified_material_does_not_become_air (__main__.ImportedThermalContracts.test_unspecified_material_does_not_become_air) ... ok

----------------------------------------------------------------------
Ran 3 tests in 0.001s

OK
Annexe C — Historique de la révision 1 et provenance des PR
Cette annexe est une photographie historique du premier audit. Les versions anciennes, résultats et défauts qui y figurent ne décrivent pas automatiquement l’état v2. En particulier, #1096 a depuis réparé CMake, les défauts de source implicite, de température de paroi, de parois SST et de stress transposé sont corrigés dans le patch actuel. Le chapitre 2 et les chapitres 4–7 font foi pour l’état courant.
La demande correspond à deux familles de travaux. Le socle thermique/rayonnement est suivi par #885 et s’est construit notamment avec #887/#888, #999/#1004, #1011 et #1016. Le socle SST transactionnel a été fusionné via #1047. Il faut compléter et corriger ces chemins, plutôt que réintroduire une ancienne branche d’intégration.
Socles déjà présents dans la base
PR
Rôle
Preuve / décision
#887 / #888
T1 conduction multimatériaux / T2 énergie fluide
Présents dans le code et le registre ; statuts machine READY, pas PASS
#994
Refus des modèles radiatifs non raccordés
Garde de production à conserver
#999 / #1004 / #1005
S2S dans l’énergie, canal de flux et irradiation externe
Primitives présentes ; préserver signes et priorité de l’override
#1011
Itération externe énergie–S2S
Déjà fusionnée ; tests de point fixe présents
#1016
Boucle écoulement–énergie de Boussinesq
Déjà fusionnée ; domaine de validité à conserver
#1047
SST opt-in transactionnel dans le noyau incompressible
Fusion 3f7a764799ba4db5fbb5e1495b85e8fc5c31c6bb, ascendance vérifiée
#1059
Diagnostics bruts et rejet des échecs énergie
Déjà fusionnée ; ne pas réintroduire l’ancien solveur de #1055
#1085
Distance normale aux interfaces CHT
Déjà fusionnée ; ne ferme pas à elle seule le contrat CHT
#1091
Oracles de fermetures turbulence
Dans la base initiale ; durcissement des comparateurs encore nécessaire
#1084
Contact thermique asymétrique
Fusionnée pendant l’audit ; incluse dans la base d’application du patch final
PR examinées et choix technique
Les mises à jour des versions figurent au chapitre 2 de la révision courante. Les liens de PR servent à retrouver les discussions ; les SHA servent à identifier le code auquel s’appliquent les conclusions.
PR
Version examinée
Décision dans ce dossier
#1045
7024a014
Ancienne intégration SST ; ne pas la fusionner comme un socle neuf après #1047
#1055
1471f9f0
Porter sélectivement les verdicts de campagne T2 ; ne pas remplacer les fichiers énergie actuels par leur ancienne version
#1068
f95d0ac9
Tracker documentaire uniquement ; pas d’implémentation des trois lacunes E2E/restart annoncées
#1079
bc7891d4
Retenue : source auxiliaire P1 normalisée par 4*pi
#1080
755d7ac0
Retenue : distinguer compensation globale et défaut local de continuité
#1081
3ada8a04
Retenue : suppression de la surcharge SST sans invariant de déformation explicite
#1084
d2f1c364
Déjà fusionnée ; formule et oracle intégré inclus dans les essais
#1086
4d51f934
Retenue : cohérence producteur/consommateur des facteurs de vue, avec limites géométriques maintenues
#1088
dd71176d, puis fefc7bb9
Retenue dans sa forme corrigée ; première version échouait sur deux cas configurés CG sans profil diffusion
#1089
13dbbd71
Retenue : régression de bornitude convective et reconstruction des flux
#1092
b010a6fa
Hors lot ciblé testé ici ; impacts gradients/Krylov à revalider séparément avant combinaison
#1093
0ca90b94
Remplacée techniquement par #1100 ; un contre-exemple de débordement est encore accepté
#1094
636c8f4a
Retenue : même production limitée dans les équations k et omega
#1095
f5b5b75c
Ne pas cumuler en bloc avec #1097 ; seuls les deux appels précoces aux validateurs énergie sont repris
#1096
314f3603
Non intégrable en bloc : CMake tronqué ; F04/H04/H05 et fixtures/tests repris sélectivement dans le lot testé
#1097
a74c75cf
Retenue pour la cohérence réellement relative des conductivités CHT
#1098
26018bc2
Recouvre les comparateurs ; ne pas cumuler mécaniquement avec #1100, plus complet sur le débordement du comparateur F34
#1099
7ad610f9
Non intégrable en l’état : erreur de compilation reproduite ; retenir #1094 pour cette correction
#1100
e0b4ca94
Retenue : comparateurs finis et harnais partagé renforcé
#1087 est fermée sans fusion dans les métadonnées consultées ; ne pas la compter comme une correction intégrée. Les autres modèles de turbulence et les autres familles numériques du dépôt ne sont pas promus par cette sélection SST.
Résultats discriminants du premier audit
Sujet
Observation sur la base ou le candidat incorrect
Observation après correction / verdict
P1, κ=1, J=0, T=1000 K
72 197,449437 W/m³ au lieu de 226 814,97676
#1079 : 226 814,97675999996, erreur relative 1,28×10⁻¹⁶
SST omega, k=.5, omega=2, S=10, F1=F2=1
55,5555556 au lieu de 32,2580645
Source limitée de #1094 vérifiée ; oracle actif supplémentaire dans le test officiel
Contact asymétrique
Mauvaise résistance quand les deux demi-distances sont supposées égales
#1084, désormais dans la base d’application : Q=8000/39 W ; températures centrales et bilan vérifiés
Producteur → consommateur S2S
Diagonale non nulle et rejet du validateur
#1086 : diagonale nulle, F12=1/π, F21=2/π pour aires 2 et 1 ; réciprocité respectée
CHT F04 sans H04
Faux rejets du solveur 1×1 à grand second membre ; corrigé F04 seul insuffisant
F04 + H04 + H05 : 24 contrôles, 0 échec, six combinaisons dt × relaxation
CHT avec H04 seul
24 contrôles physiques, 24 échecs
Le gel des états physiques F04 est indispensable en plus du bon résidu scalaire
H05, cellule cisaillée
Bilan matriciel nul mais diagnostic relatif 0,292893 / 0,683772
Correction distance normale : maximum 7,96×10⁻¹⁶ ; 12 contrôles réussis
Comparateur #1093
Opérandes finis ±DBL_MAX avec seuil débordant acceptés
#1100 rejette le débordement ; harnais partagé durci et utilisateurs ciblés retestés
Compilation #1099
Variables SST hors portée et membres inexistants d’un objet RNG k-epsilon
Échec reproduit ; #1094 retenue à sa place
Configuration #1096
Flow control statements are not properly nested à CMakeLists.txt:1080
Extraits thermiques enregistrés dans le CMake complet ; PR entière non approuvée
Pour F04, les six configurations dt={0.01,0.1,1} × relaxation {0.5,1} convergent en 9 à 152 itérations. Sur l’overlay F04+H04+H05, les erreurs maximales sont 2,38×10⁻¹¹ K et 3,47×10⁻¹¹ K, avec un défaut capacitif maximal de 1,52×10⁻¹⁰ J. Le même test est ensuite exécuté avec succès dans le lot final.
L’isolement des causes a été conservé : base F04, 12 contrôles atteints / 12 échecs avec quatre exceptions ; F04+H05 sans H04, 12 contrôles / 4 échecs ; H04 seul, 24/24 échecs physiques ; F04+H04+H05, 24/0. Les nombres de contrôles atteints diffèrent parce que les exceptions interrompent certaines configurations : ils ne doivent pas être additionnés comme un total de validations indépendantes.
La conduction avec puissance volumique donne un ordre mesuré de 2 sur les trois maillages du test. Ce résultat concerne cette famille orthogonale et cette équation ; il n’est pas transférable au couplage SST, au rayonnement ou à tous les maillages du dépôt.
Annexe D — Sondes indépendantes complémentaires
Les tests intégrés figurent déjà dans le patch. Les sondes ci-dessous documentent deux comparaisons supplémentaires ; leurs anciens résultats doivent être obtenus sur l’état ancien spécifié, pas en compilant deux fois les headers corrigés. Elles ne constituent pas un nouveau pilote de production.
Mesures du point fixe radiatif et commandes exactes
{
  "recorded_at": "2026-10-24T16:54:10.445746+00:00",
  "scope": "Same current thermal implementation; previous S2S header vs corrected S2S header. Actual applied boundary heat flux recovered from independent energy balance.",
  "source_sha256": "92ed5f2d067b81b9ae5d71980a9b977ef39b7b28f8ef2d49e26cd334d90fd53d",
  "results": [
    {
      "name": "fixed_point_previous",
      "compile_command": [
        "mpic++",
        "-std=c++17",
        "-O2",
        "-DNDEBUG",
        "-DCFDX_HAS_MPI=1",
        "-DCFDX_HAS_HDF5=1",
        "-I/usr/include/hdf5/serial",
        "-I/workspace/scratch/cae91c5c59b2/audit_evidence/radiation_wall_previous",
        "-I/workspace/scratch/cae91c5c59b2/cfdx_integration_v2/src",
        "-I/workspace/scratch/cae91c5c59b2/cfdx_integration_v2/tests",
        "/workspace/scratch/cae91c5c59b2/audit_evidence/radiation_wall_fixed_point_probe.cpp",
        "/workspace/scratch/cae91c5c59b2/build_baseline/libcfdx_core.a",
        "-L/usr/lib/x86_64-linux-gnu/hdf5/serial",
        "-lhdf5",
        "-lcrypto",
        "-lcurl",
        "-lsz",
        "-lz",
        "-ldl",
        "-lm",
        "-o",
        "/workspace/scratch/cae91c5c59b2/audit_evidence/fixed_point_previous_v2_verified"
      ],
      "compile_returncode": 0,
      "compile_output": "[32e80b77a0e4:00006] opal_ifinit: socket() failed with errno=1\n",
      "returncode": 0,
      "metrics": {
        "converged": true,
        "outer_iterations": 11,
        "reported_wall_bottom_K": 357.0037421510402,
        "reported_wall_top_K": 354.6034007061026,
        "fourier_wall_bottom_K": 361.001593585563,
        "fourier_wall_top_K": 356.85554927157983,
        "reported_radiative_flux_W_m2": 18.017188534245065,
        "fourier_based_radiative_flux_W_m2": 31.948514245302267,
        "radiative_constitutive_error_W_m2": 13.931325711057202,
        "energy_global_balance_normalized": 4.547473508864639e-16,
        "radiative_balance_normalized": 0,
        "fixed_point_flux_normalized": 6.375709443131633e-12
      }
    },
    {
      "name": "fixed_point_corrected",
      "compile_command": [
        "mpic++",
        "-std=c++17",
        "-O2",
        "-DNDEBUG",
        "-DCFDX_HAS_MPI=1",
        "-DCFDX_HAS_HDF5=1",
        "-I/usr/include/hdf5/serial",
        "-I/workspace/scratch/cae91c5c59b2/cfdx_integration_v2/src",
        "-I/workspace/scratch/cae91c5c59b2/cfdx_integration_v2/tests",
        "/workspace/scratch/cae91c5c59b2/audit_evidence/radiation_wall_fixed_point_probe.cpp",
        "/workspace/scratch/cae91c5c59b2/build_baseline/libcfdx_core.a",
        "-L/usr/lib/x86_64-linux-gnu/hdf5/serial",
        "-lhdf5",
        "-lcrypto",
        "-lcurl",
        "-lsz",
        "-lz",
        "-ldl",
        "-lm",
        "-o",
        "/workspace/scratch/cae91c5c59b2/audit_evidence/fixed_point_corrected_v2_verified"
      ],
      "compile_returncode": 0,
      "compile_output": "[32e80b77a0e4:00013] opal_ifinit: socket() failed with errno=1\n",
      "returncode": 0,
      "metrics": {
        "converged": true,
        "outer_iterations": 23,
        "reported_wall_bottom_K": 360.259797296162,
        "reported_wall_top_K": 357.5973455609808,
        "fourier_wall_bottom_K": 360.259797296162,
        "fourier_wall_top_K": 357.5973455609808,
        "reported_radiative_flux_W_m2": 20.515870904989242,
        "fourier_based_radiative_flux_W_m2": 20.515870904989228,
        "radiative_constitutive_error_W_m2": 1.4210854715202004e-14,
        "energy_global_balance_normalized": 0,
        "radiative_balance_normalized": 0,
        "fixed_point_flux_normalized": 8.454522104912623e-11
      }
    }
  ]
}
Comparaison du point fixe radiatif avec la température réelle de frontière
#include "cfdx/physics/production_radiation_s2s_stage.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::physics;
static Mesh box_mesh(std::size_t nx, std::size_t ny, double dx, double dy,
                     PatchType front_type = PatchType::WALL)
{
    const auto pid = [&](std::size_t i, std::size_t j, std::size_t l) {
        return (i * (ny + 1) + j) * 2 + l;
    };
    const std::size_t n_points = 2 * (nx + 1) * (ny + 1);
    const std::size_t base_fy = (nx + 1) * ny;
    const std::size_t base_fz = base_fy + nx * (ny + 1);
    const std::size_t n_faces = base_fz + 2 * nx * ny;
    const auto fx = [&](std::size_t i, std::size_t j) { return i * ny + j; };
    const auto fy = [&](std::size_t i, std::size_t j) {
        return base_fy + i * (ny + 1) + j;
    };
    const auto fz = [&](std::size_t i, std::size_t j, std::size_t l) {
        return base_fz + 2 * (i * ny + j) + l;
    };
    const auto cell = [&](std::size_t i, std::size_t j) {
        return i * ny + j;
    };

    Mesh m;
    m.points().resize(n_points);
    for (std::size_t i = 0; i <= nx; ++i)
        for (std::size_t j = 0; j <= ny; ++j)
            for (std::size_t l = 0; l < 2; ++l)
                m.points().set(pid(i, j, l),
                               static_cast<double>(i) * dx,
                               static_cast<double>(j) * dy,
                               static_cast<double>(l));

    for (std::size_t i = 0; i <= nx; ++i)
        for (std::size_t j = 0; j < ny; ++j) {
            // Face normals point owner-to-neighbour for internal faces and
            // outward from the owner for boundary faces.
            if (i == 0)
                m.faces().push_face({pid(i, j, 0), pid(i, j, 1),
                                     pid(i, j + 1, 1), pid(i, j + 1, 0)});
            else
                m.faces().push_face({pid(i, j, 0), pid(i, j + 1, 0),
                                     pid(i, j + 1, 1), pid(i, j, 1)});
        }
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j <= ny; ++j) {
            if (j == 0)
                m.faces().push_face({pid(i, j, 0), pid(i + 1, j, 0),
                                     pid(i + 1, j, 1), pid(i, j, 1)});
            else
                m.faces().push_face({pid(i, j, 0), pid(i, j, 1),
                                     pid(i + 1, j, 1), pid(i + 1, j, 0)});
        }
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t l = 0; l < 2; ++l) {
                if (l == 0)
                    m.faces().push_face({pid(i, j, l), pid(i, j + 1, l),
                                         pid(i + 1, j + 1, l), pid(i + 1, j, l)});
                else
                    m.faces().push_face({pid(i, j, l), pid(i + 1, j, l),
                                         pid(i + 1, j + 1, l), pid(i, j + 1, l)});
            }

    m.ownership().resize(n_faces);
    for (std::size_t j = 0; j < ny; ++j) {
        m.ownership().set_owner(fx(0, j), cell(0, j));
        m.ownership().set_neighbour(fx(0, j), FaceOwnership::BOUNDARY);
        m.ownership().set_owner(fx(nx, j), cell(nx - 1, j));
        m.ownership().set_neighbour(fx(nx, j), FaceOwnership::BOUNDARY);
    }
    for (std::size_t i = 1; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j) {
            m.ownership().set_owner(fx(i, j), cell(i - 1, j));
            m.ownership().set_neighbour(fx(i, j),
                                        static_cast<int>(cell(i, j)));
        }
    for (std::size_t i = 0; i < nx; ++i) {
        m.ownership().set_owner(fy(i, 0), cell(i, 0));
        m.ownership().set_neighbour(fy(i, 0), FaceOwnership::BOUNDARY);
        m.ownership().set_owner(fy(i, ny), cell(i, ny - 1));
        m.ownership().set_neighbour(fy(i, ny), FaceOwnership::BOUNDARY);
    }
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 1; j < ny; ++j) {
            m.ownership().set_owner(fy(i, j), cell(i, j - 1));
            m.ownership().set_neighbour(fy(i, j),
                                        static_cast<int>(cell(i, j)));
        }
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t l = 0; l < 2; ++l) {
            m.ownership().set_owner(fz(i, j, l), cell(i, j));
            m.ownership().set_neighbour(fz(i, j, l),
                                        FaceOwnership::BOUNDARY);
        }

    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
            m.cells().push_cell({fx(i, j), fx(i + 1, j),
                                 fy(i, j), fy(i, j + 1),
                                 fz(i, j, 0), fz(i, j, 1)});

    Patch left;
    left.name = "left";
    left.type = PatchType::WALL;
    for (std::size_t j = 0; j < ny; ++j) left.face_ids.push_back(fx(0, j));
    Patch right;
    right.name = "right";
    right.type = PatchType::WALL;
    for (std::size_t j = 0; j < ny; ++j) right.face_ids.push_back(fx(nx, j));
    Patch bottom;
    bottom.name = "bottom";
    bottom.type = PatchType::WALL;
    for (std::size_t i = 0; i < nx; ++i) bottom.face_ids.push_back(fy(i, 0));
    Patch top;
    top.name = "top";
    top.type = PatchType::WALL;
    for (std::size_t i = 0; i < nx; ++i) top.face_ids.push_back(fy(i, ny));
    Patch front;
    front.name = "front";
    front.type = front_type;
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
            front.face_ids.push_back(fz(i, j, 0));
    Patch back;
    back.name = "back";
    back.type = PatchType::WALL;
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
            back.face_ids.push_back(fz(i, j, 1));
    m.boundary().add_patch(left);
    m.boundary().add_patch(right);
    m.boundary().add_patch(bottom);
    m.boundary().add_patch(top);
    m.boundary().add_patch(front);
    m.boundary().add_patch(back);
    return m;
}

// ---------------------------------------------------------------------------
// Case: conjugate energy + S2S radiation. Dirichlet walls left (400 K) and
// right (300 K); radiative participants bottom (eps 0.8, declared heat
// flux 50 W/m^2, view factor 1 to top) and top (eps 0.9, view factor 1 to
// bottom); front / back without spec (adiabatic).
// ---------------------------------------------------------------------------
static CaseSetup make_s2s_setup()
{
    CaseSetup setup;
    setup.energy_model = "conjugate";
    setup.radiation_model = "S2S";

    MaterialSpec water;
    water.name = "water";
    water.density = 1000.0;
    water.dynamic_viscosity = 1e-3;
    water.thermal_conductivity = 1.0;
    water.specific_heat = 1000.0;
    setup.materials.push_back(water);

    BoundarySpec left;
    left.patch_name = "left";
    left.type = BCType::WALL;
    left.value_type = BCValueType::WALL_THERMAL;
    left.wall_temperature = 400.0;
    setup.boundary_conditions.push_back(left);

    BoundarySpec right;
    right.patch_name = "right";
    right.type = BCType::WALL;
    right.value_type = BCValueType::WALL_THERMAL;
    right.wall_temperature = 300.0;
    setup.boundary_conditions.push_back(right);

    BoundarySpec bottom;
    bottom.patch_name = "bottom";
    bottom.type = BCType::WALL;
    bottom.raw_params["emissivity"] = "0.8";
    bottom.raw_params["heat_flux"] = "50";
    bottom.raw_params["s2s_view_factor_to_top"] = "1";
    setup.boundary_conditions.push_back(bottom);

    BoundarySpec top;
    top.patch_name = "top";
    top.type = BCType::WALL;
    top.raw_params["emissivity"] = "0.9";
    top.raw_params["s2s_view_factor_to_bottom"] = "1";
    setup.boundary_conditions.push_back(top);

    // front / back: deliberately no spec (adiabatic walls).
    return setup;
}


int main() {
    auto mesh=box_mesh(4,4,0.25,0.25);
    const auto geometry=build_fv_geometry(mesh);
    const auto setup=make_s2s_setup();
    Field<double,Location::FACE> phi(mesh.n_faces(),"phi","kg/s",1); phi.fill(0.0);
    Field<double,Location::CELL> T(mesh.n_cells(),"T","K",1); T.fill(350.0);
    const auto r=run_production_s2s_energy_outer_loop(mesh,geometry,setup,phi,T);
    std::size_t b=0,t=0;
    for(std::size_t i=0;i<r.radiation.surface_count;++i){
        if(r.radiation.surface_names[i]=="bottom") b=i;
        if(r.radiation.surface_names[i]=="top") t=i;
    }
    const auto wall_from_energy=[&](const std::string& name,double qinto){
        double sum=0.0,area=0.0;
        for(std::size_t p=0;p<mesh.boundary().n_patches();++p){
            const auto& patch=mesh.boundary().patch(p);if(patch.name!=name)continue;
            for(const auto f:patch.face_ids){
                const auto c=mesh.ownership().owner(f);
                const auto sf=geometry.face_area_vectors[f];
                const auto d=geometry.face_centres[f]-geometry.cell_centres[c];
                const double a=sf.mag(), dn=std::abs(d.dot(sf))/a;
                sum+=a*(T(c)+qinto*dn);area+=a;
            }
        } return sum/area;
    };
    const double qb=r.radiation.net_flux[b];
    // Pure conduction: the independent energy flux audit recovers the
    // actually applied inward flux, before any fixed-point response update.
    const double bottom_q_into=-r.energy.patch_flux_out.at("bottom")/r.radiation.areas[b];
    const double top_q_into=-r.energy.patch_flux_out.at("top")/r.radiation.areas[t];
    const double tb=wall_from_energy("bottom",bottom_q_into);
    const double tt=wall_from_energy("top",top_q_into);
    const double physical_q=STEFAN_BOLTZMANN*(std::pow(tb,4)-std::pow(tt,4))/(1.0/0.8+1.0/0.9-1.0);
    std::cout<<std::setprecision(17)
        <<"{\"converged\":"<<(r.converged?"true":"false")
        <<",\"outer_iterations\":"<<r.outer_iterations
        <<",\"reported_wall_bottom_K\":"<<r.radiation.wall_temperatures[b]
        <<",\"reported_wall_top_K\":"<<r.radiation.wall_temperatures[t]
        <<",\"fourier_wall_bottom_K\":"<<tb
        <<",\"fourier_wall_top_K\":"<<tt
        <<",\"reported_radiative_flux_W_m2\":"<<qb
        <<",\"fourier_based_radiative_flux_W_m2\":"<<physical_q
        <<",\"radiative_constitutive_error_W_m2\":"<<std::abs(physical_q-qb)
        <<",\"energy_global_balance_normalized\":"<<r.energy.global_balance_normalized
        <<",\"radiative_balance_normalized\":"<<r.radiation.balance_normalized
        <<",\"fixed_point_flux_normalized\":"<<r.max_flux_fixed_point_normalized<<"}\n";
    return r.converged?0:1;
}
Contre-exemple du résidu SST après bornage
#include "audit/chain_fixture.h"
#include "cfdx/physics/sst_solver.h"
#include <iostream>
#include <iomanip>
using namespace cfdx::physics;
int main() {
 const auto mesh=make_1d_chain(1); const auto geometry=build_fv_geometry(mesh);
 Field<double,Location::FACE> flux(mesh.n_faces(),"phi","",1);flux.fill(0);
 Field<double,Location::CELL> k(1,"k","",1),w(1,"omega","",1),S(1,"S","",1),F(1,"F","",1);
 k.fill(1);w.fill(1);S.fill(0);F.fill(1);
 TurbulenceTransportControls c;c.model=TurbulenceModel::SST;c.k_min=1;c.omega_min=1;
 const auto result=solve_sst_transport(mesh,geometry,flux,k,w,S,F,F,c,{}, {},3,1e-8);
 std::cout<<std::setprecision(17)<<"converged="<<result.converged<<" k="<<k(0)<<" omega="<<w(0)
 <<" k_residual="<<result.k_residual<<" omega_residual="<<result.second_residual
 <<" iterations="<<result.iterations<<" k_floor_events="<<result.k_floor_events
 <<" omega_floor_events="<<result.second_floor_events<<"\n";
}
Décision de revue
Le lot exact ci-dessus est prêt pour une revue de code et les contrôles obligatoires du dépôt, avec les preuves locales documentées. Les modifications corrigeant les opérateurs, les diagnostics et les contrats peuvent être examinées sans ouvrir prématurément le sélecteur RANS. L’intégration multiphysique complète sera acquise lorsque les pilotes, le format des cas, la reprise globale et les campagnes de référence du chapitre 7 seront effectivement réalisés et satisfaits.