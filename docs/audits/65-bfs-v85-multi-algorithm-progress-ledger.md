# FICHE MÉTHODE — Progression non monotone bornée V85 pour les solveurs pression-vitesse

**En-tête d'évidence (gouvernance §16.5)**
- Scope : politique de progression externe V85 appliquée aux résidus physiques réassemblés V8 ; SIMPLE, SIMPLEC et COUPLED avec ses six variantes de Schur. Hors scope : modification des équations FV, des préconditionneurs linéaires, des conditions aux limites, des relaxations et des critères finaux.
- Statut : DRAFT (2026-11-20) — revue et qualification numérique requises avant APPROVED.
- Registre autoritaire : docs/validation/CFDX_QUALIFICATION_REGISTRY.json + NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json.
- Issue / PR liées : #1196.
- Dernier commit vérifié : e6d1d86f0e11ac5607c091ff37ae8c89381c488c (2026-11-20).
- Limitations connues : aucune preuve de convergence BFS ni de reproduction Driver–Seegmiller n'est apportée par cette politique. Les campagnes par algorithme et par Schur doivent être rejouées.

## 1. Références sources (obligatoire : DOI ou ISBN)

- C. T. Kelley, *Iterative Methods for Linear and Nonlinear Equations*, SIAM, 1995, ISBN 978-0-89871-352-7. Référence pour les méthodes de résolution non linéaire, les critères de résidu et les stratégies de globalisation.
- M. A. Elman, D. J. Silvester, A. J. Wathen, *Finite Elements and Fast Iterative Solvers*, 2e éd., Oxford University Press, 2014, ISBN 978-0-19-967879-2. Référence de contexte pour les systèmes blocs et les approximations de Schur ; V85 ne modifie pas ces opérateurs.

## 2. Formulation continue

V85 ne définit aucune nouvelle équation physique. Il gouverne uniquement l'acceptation provisoire d'états successifs pour le système stationnaire incompressible avec SST :

\[
\mathcal{R}(U,p,\phi,k,\omega)=
\begin{bmatrix}
\mathcal{R}_{U}(U,p,\phi,\nu_t)\\
\mathcal{R}_{c}(U,\phi)\\
\mathcal{R}_{k}(U,\phi,k,\omega)\\
\mathcal{R}_{\omega}(U,\phi,k,\omega)
\end{bmatrix}=0.
\]

La solution finale reste celle qui satisfait les critères physiques originaux V8, calculés par réassemblage indépendant sur l'état retourné. V85 n'est ni un critère d'arrêt ni une tolérance alternative.

## 3. Discrétisation Volume Fini

Les composantes utilisées par le ledger sont les résidus physiques réassemblés V8 : trois équations momentum, continuité conservative basée sur le flux de face autoritaire, transport de k, transport de omega et fermeture de viscosité turbulente. Les résidus de Krylov, les estimations AMG et les indicateurs de moniteur ne remplacent aucune de ces quantités.

Pour un état candidat \(s_n\), le mérite est le maximum des résidus normalisés indépendants :

\[
M_n=\max\left(
r_{U,global},r_{U,max},r_c,r_k,r_\omega,r_{\nu_t}
\right).
\]

La politique compare les états uniquement si l'opérateur de convection reste identique. Chaque historique est isolé par couple (algorithme pression-vitesse, modèle de Schur) ; un changement d'identité dans un ledger est rejeté.

Paramètres V85 du code :
- fenêtre historique : 6 états ;
- croissance bornée du mérite par rapport au meilleur état de la fenêtre : 1,01 par défaut ;
- garde momentum : croissance maximale de 2 % par rapport à l'état précédent ;
- garde continuité : au plus 5 fois la tolérance finale, avec contrôle relatif au précédent ;
- dette SST : excursions bornées par les multiplicateurs 8 (meilleur état) et 3,5 (état précédent), avec récupération exigée avant épuisement de la fenêtre ;
- fermeture de viscosité turbulente : garde dédiée ;
- au moins une famille de résidus physiques doit progresser.

Les constantes ci-dessus sont des limites d'acceptation provisoire, pas des tolérances de convergence.

## 4. Traitement des frontières

V85 n'ajoute ni ne modifie de condition de Dirichlet, Neumann ou Robin. Les contributions de frontière restent celles de l'assemblage physique pression-vitesse et SST utilisé par le solveur.

Le flux de masse aux faces demeure l'état conservatif autoritaire. Il n'est jamais remplacé par une interpolation reconstruite de la vitesse pour faire paraître le résidu de continuité plus faible. Les contrôles de conservation, de positivité et de bornitude du solveur restent actifs.

## 5. Pièges connus (Gotchas)

- Ne pas confondre progression provisoire et convergence : seul le certificat final V8 peut déclarer la convergence.
- Ne pas partager un historique entre SIMPLE, SIMPLEC ou deux modèles de Schur différents.
- Ne pas mélanger deux opérateurs de convection dans la même comparaison.
- Le transport SST provisoire V83 reste qualifié uniquement pour SIMPLE. Pour SIMPLEC et COUPLED, V85 ne peut être activé qu'avec un état SST entièrement commis ; aucune extension de V83 n'est revendiquée.
- Une simulation peut satisfaire le ledger V85 tout en restant physiquement non convergée. Les profils BFS, le recollement, Cf, Cp, les contraintes de Reynolds et la convergence de maillage doivent être évalués séparément.
- Un succès de compilation ou de CI ne constitue pas une qualification physique.

## 6. Ordre de précision attendu et validation exigée

V85 ne change pas l'ordre spatial ou temporel du schéma FV : il ne modifie aucune équation ni aucun coefficient d'assemblage. Il ne peut donc pas être crédité d'un gain d'ordre ou d'une amélioration de précision physique.

Validation logicielle exigée :
- tests de l'acceptation, de la dette SST, de la bornitude, de la continuité et du refus des faux certificats ;
- tests d'historiques indépendants SIMPLEC et COUPLED ;
- test des six configurations COUPLED : BlockLocal, PCD, LSC, BFBT, SIMPLE-Schur et SIMPLEC-Schur ;
- rejet explicite d'un changement d'algorithme ou de Schur dans un historique ;
- CI DebugSanitizers et Release.

Validation numérique/physique exigée séparément :
- même maillage, conditions aux limites, initialisation et budgets pour chaque modèle ;
- résidus vrais des systèmes linéaires et résidus FV non linéaires réassemblés ;
- conservation de masse et cohérence flux-vitesse ;
- profils U/Uref et V/Uref aux stations expérimentales, Cf, Cp, contraintes de Reynolds, y+ et position de recollement ;
- au moins trois maillages pour une étude de sensibilité/GCI lorsque les résultats convergent.

Aucune tolérance physique ou linéaire ne doit être relâchée pour faire passer cette validation.

## 7. Statut d'implémentation CFDX

- Fichiers source concernés : src/cfdx/physics/global_outer_progress_v85.h ; src/cfdx/physics/steady_incompressible_solver.h.
- Tests exécutables : tests/unit/test_global_outer_progress_v85.cpp ; cible CMake test_global_outer_progress_v85.
- État : extension du contrat V85 pour les algorithmes avec SST commis ; V83 provisoire reste SIMPLE-only. CI déclenchée sur la PR #1196 ; résultats à vérifier avant approbation.
- Écarts connus code ↔ fiche : la fiche est DRAFT ; aucune campagne comparative BFS n'a encore démontré la convergence des modèles ni leur précision physique.
- Gaps et prochaines actions : valider les tests et les builds ; préparer une matrice de campagnes strictement comparable pour SIMPLE, SIMPLEC et les six Schur ; analyser ensuite les profils physiques.

## 8. Journal des modifications de la fiche

- 2026-11-20 — PR #1196 — création de la fiche pour l'extension multi-algorithmes V85.
