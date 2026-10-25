# CFDX - Etat d'integration et reste a faire

<a id="integration-2026-10-25"></a>
## 13. Point d'integration du 25 octobre 2026 : code fusionne, candidats locaux et travaux restants

### 13.1 Verdict et limite temporelle

**L'integration a avance, mais les preuves disponibles ne portent pas sur un seul commit reunissant tous les correctifs.** Il faut distinguer le master, les branches de PR, le patch R4 et le lot thermique/rayonnement/turbulence v2. Le nombre de correctifs locaux testes ne constitue pas un pourcentage d'integration GitHub.

Cette section est une **revision documentaire et une revue croisee des patches**, pas un nouveau patch de production. Le patch embarque et les preuves R4 restent inchanges. Les resultats v2 cites ci-dessous sont des executions anterieures rapportees, non des executions refaites aujourd'hui. La seule nouvelle execution est la sonde isolee I01 en section 13.4.

Consultation publique tentee le **25 octobre 2026**. Les pages accessibles sont des captures incompatibles entre elles : 19 PR dans la liste, 27 dans d'autres pages ; #1084 est encore affichee ouverte alors que le dossier v2 la documente deja fusionnee dans sa base. L'API du dernier master, les diffs et les checks actuels n'ont pas ete recuperes ; les telechargements directs ont echoue. **Aucun HEAD courant n'a ete telecharge, compile ou teste dans cette revision.** Un echec de lecture ne prouve ni absence de correctif ni absence de fusion.

La derniere inspection Git precise retrouvee dans les dossiers date du **24 octobre 2026 a 16:46 UTC (18:46 a Paris)** [D2, chapitre 2]. Les identifiants ci-dessous sont donc dates, et non presentes comme des HEAD garantis actuels.

| Source | Revision / identite | Ce qu'elle prouve |
|---|---|---|
| D1 - Audit consolide R4 | ZIP `3f7a764799ba4db5fbb5e1495b85e8fc5c31c6bb` + patch local 54 fichiers | Correctifs locaux et trois prototypes, sans integration generale de F03/F06/F11 |
| D2 - Integration thermique v2 | Base `038e3aeff40f02af95b8845b931245ce23a1a094` ; arbre local `89adbb733827bdc932ff546ecf701eec11fe65f6` | Lot local de 58 fichiers, dont de vrais raccordements applicatifs ; aucun push/merge par ce dossier |
| D3 - Audit de #1092 | `7ccb31d8e2fb5882acf0fe1beb9ecfbb616bf227` | Defauts et regressions sur CETTE ancienne revision uniquement |
| D2 - Releve plus recent de #1092 | `15bb69c7917518e2bacf1017a4002bbddc65d57b` | La branche a change ; le verdict de D3 doit etre reexecute, pas recopie |
| Pages publiques reconsultees | #1047, #1067, #1068, #1084, #1094 | Etat affiche et declarations ; pas un arbre source synchronise du master |

### 13.2 Avancement par groupe de PR

| Groupe / PR | Etat suffisamment etaye | Reste a faire |
|---|---|---|
| #1047 | Fusion confirmee par la page publique, commit `3f7a764...` [W1]. Runtime SST opt-in present. | Ne pas annoncer RANS general : creation depuis un vrai cas et reprise complete restent des exigences distinctes. |
| #887/#888, #999/#1004/#1005, #1011/#1016, #1059 | D2 les identifie deja dans sa base : conduction/energie, S2S, boucle energie-S2S, Boussinesq, diagnostics energie. | Garder les capacites existantes ; READY ne signifie pas PASS de qualification physique. |
| #1085 / #1091 | D2 les identifie dans sa base : distance normale CHT et oracles turbulence. | Ni le contrat CHT complet ni la robustesse generale des comparateurs ne sont clos par ces seules fusions. |
| #1084 | D2 documente la fusion dans `038e3ae...`, HEAD de PR `d2f1c36426669623ce3c9d23706bbdfc3e593078`. La page publique est plus ancienne/contradictoire [W4]. | Ne pas reintroduire l'ancienne formule. Retester l'oracle integre Q=8000/39 W sur la cible. |
| #1079/#1080/#1081/#1086/#1088/#1089/#1094/#1097/#1100 | Changements selectionnes et testes dans D2, pas preuve de fusion distante. | Integrer ces versions sans perdre leurs corrections ni leurs tests ; CI du commit combine. |
| #1096 | HEAD date `265ad02b432fcf45c63ded89926eb75d422bfb10`. D2 confirme que la configuration CMake complete passe. | L'ancien defaut CMake tronque n'est plus a presenter comme courant. Compiler/tester toute la PR ; v2 n'en prend que F04/H04/H05 et certains tests. |
| #1092 | HEAD date `15bb69c7917518e2bacf1017a4002bbddc65d57b`, absent du lot D2. | Reauditer G01-G09 puis combiner avec D2 ; les 43 CTest thermiques ne qualifient pas cette PR. |
| #1078 | HEAD date `8c00dfd20e90ec648ff200a3acd40d918ff99a74`, absent du lot D2. | Verifier simultanement flux Rhie-Chow, continuite et Schur. |
| #1095 / #1097 | D2 retient #1097 et deux validations precoces de #1095. | Integration selective du contrat CHT, pas deux merges aveugles de versions du meme fichier. |
| #1093/#1098/#1100 | D2 privilegie #1100 `e0b4ca94` pour les comparateurs, dont l'overflow. | Conserver les tests discriminants utiles, ne pas cumuler les implementations concurrentes. |
| #1099 | `7ad610f9e3d189c6db8912c0098d12318982523d` exclu de D2 apres echec de compilation reproduit. | Ne pas l'integrer en l'etat documente ; retenir #1094 pour sa correction SST, reevaluer seulement un nouveau SHA. |
| #1045 / #1055 | Ancienne integration SST / ancien etage energie. | Eviter de restaurer des fichiers anciens : #1047 et #1059 existent ; porter uniquement les verdicts utiles de campagne T2. |
| #1068 | La capture publique est un draft documentaire pour #1067 [W2,W3]. D2 apporte certaines preuves locales supplementaires. | Publier/raccorder ces preuves et finir fabrique de runtime et reprise de trajectoire ; ne pas confondre tracker et implementation. |
| #1083 | PR de courbure reperee, hors lot thermique ; dernier diff non recupere. | Revue et tests VOF separes avant promotion. |

**Progression a retenir :** F11 n'est plus seulement un prototype dans le travail LOCAL le plus recent. D2 ajoute `momentum_viscous_stress.h`, utilise sa contribution dans le predicteur, le residu final et le bilan de quantite de mouvement. La MMS du terme donne des ordres 1,96692 et 1,98397 et un test observe un effet dans le predicteur reel. Cela ne prouve ni fusion ni qualification sur tous les murs/maillages [D2, chapitre 6].

Les essais historiques de D2 sont **43/43 CTest Release, 9/9 CTest ASan/UBSan, 87 tests Python reussis et 4 ignores faute de PyFluent**. Ils concernent son arbre local, sont cibles, et ne doivent pas etre additionnes aux compteurs R4. Compiler avec MPI ne prouve pas une execution multirang.

### 13.3 Le probleme d'integration transversal : 15 fichiers communs

Comparaison effectuee aujourd'hui sur les blocs `diff --git` des deux Markdown : **54 fichiers dans R4, 58 dans v2, 15 chemins communs**. Il s'agit de chevauchements de fichiers, pas de 15 conflits Git prouves. Les bases sont differentes ; l'absence d'un fichier dans un patch ne prouve pas que son correctif manque dans la base.

| Groupe | Fichiers communs / enjeu | Regle de resolution proposee |
|---|---|---|
| Construction | `CMakeLists.txt` | Preserver la configuration complete et l'union des tests pertinents ; detecter tests manquants et doublons dans CTest, sans convertir une campagne omise en succes. |
| Transport/energie | `finite_volume_transport.h`, `energy_solver.h` | Conserver flux/matrice, sources Su+Sp*T, conductivite variable, metrique normale, diagnostics avant/apres et gardes non finies. |
| Interfaces | `cht_solver.h`, `thermal_zones.h` | Conserver T^n fixe, validation des deux regions avant mutation, k coherents, distance normale, contact asymetrique et impermeabilite. |
| Fermetures | `sst_solver.h`, `m1_m4_models.h` | Source limitee, BC reelles, residuals apres bornage et P1 normalise ; pas de retour aux versions anterieures. |
| Essais communs | Huit chemins de fixtures/tests | Comparer les oracles et le perimetre ; un test dupliquant une formule n'a pas la valeur d'un oracle independant. |

Chemins communs complets :

- `CMakeLists.txt`
- `src/cfdx/physics/cht_solver.h`
- `src/cfdx/physics/energy_solver.h`
- `src/cfdx/physics/finite_volume_transport.h`
- `src/cfdx/physics/m1_m4_models.h`
- `src/cfdx/physics/sst_solver.h`
- `src/cfdx/physics/thermal_zones.h`
- `tests/audit/chain_fixture.h`
- `tests/audit/test_audit_cht_time.cpp`
- `tests/audit/test_audit_energy_boundary_metric.cpp`
- `tests/unit/test_finite_volume_transport.cpp`
- `tests/unit/test_turbulence_qualification_equations.cpp`
- `tests/validation/test_cht_validation.cpp`
- `tests/validation/test_level_c_coupled_verification.cpp`
- `tests/validation/test_thermal_challenge_vv.cpp`

Les fichiers de solveurs Krylov modifies par R4 ne sont pas dans le patch thermique v2. Les applications et le nouveau tenseur visqueux de v2 ne sont pas dans R4. **Aucun de ces deux ensembles ne remplace automatiquement l'autre.** L'objectif doit etre un nouvel arbre combine identifie et teste, pas l'application consecutive de deux patches sur une base non controlee.

### 13.4 I01 - Garde de residu a ne pas perdre lors de l'integration

**Constat nouveau de cette revue croisee :** la lambda `physical_residual_inf` ajoutee par le patch v2, autour des lignes 2508-2525 de D2, contient encore `norm = std::max(norm, std::abs(ri))`, sans verifier la finitude de `ri`. Sa reduction du second membre suit le meme motif. Le correctif H08 du helper R4 ne protege donc pas automatiquement ce second diagnostic.

Une sonde C++ a ete compilee aujourd'hui avec `g++ -std=c++17 -O2 -Wall -Wextra -pedantic`. Le corps de la lambda est **extrait automatiquement du patch v2**, avec de petits conteneurs compatibles pour exposer sa reduction. Il ne s'agit PAS d'une execution de `solve_scalar_equation`, ni d'une preuve que ses validations amont laisseraient atteindre tous ces cas. C'est une preuve du contrat insuffisant de ce reducteur lorsqu'il recoit un residu non fini. Un cas utilise uniquement des entrees finies dont les produits debordent avant de se compenser.

| Cas de la sonde | Reduction v2 extraite | Proposition avec gardes | Verdict de la proposition |
|---|---:|---:|---|
| Solution finie exacte | 0 | 0 | Conforme |
| Residu fini non nul | 3 | 3 | Conforme |
| NaN dans le champ | **0** | +Inf | Rejet explicite |
| NaN dans un coefficient | **0** | +Inf | Rejet explicite |
| NaN dans le RHS | **0** | +Inf | Rejet explicite |
| Entrees finies, overflow puis Inf-Inf | **0** | +Inf | Rejet explicite |

**6 cas de la proposition passent ; 4 faux zeros du reducteur extrait sont reproduits.** Ces essais ne sont pas ajoutes aux 43/43 CTest historiques ni aux 56 constats historiques. I01 est un point de raccordement F32/H08 ; la proposition de garde n'est pas appliquee au depot.

Correction requise : valider dimensions/indices et finitude des donnees, des produits et de l'accumulation ; retourner un echec ou +Inf explicite sur calcul non fini ; proteger aussi la norme et la reference du RHS. Recalculer les diagnostics du champ apres relaxation. Enfin, ajouter les cas au vrai solveur et aux consommateurs energie/SST, afin de prouver le chemin complet au lieu de s'arreter a cette sonde. Une strategie de norme robuste ou de rescaling peut eviter certains overflows ; accepter un faux zero ne le peut pas.

### 13.5 Ordre d'integration conseille et criteres de cloture

Les lots suivants sont un **plan d'implementation propose**, pas des travaux annonces comme deja effectues. Les priorites indiquent les dependances de livraison ; elles ne remplacent pas la severite scientifique de chaque constat historique.

#### Lot 0 - Figer les sources et rendre les preuves de fusion non ambigues

Recuperer le master et chaque HEAD de PR, enregistrer date UTC, SHA, base, arbre, et statut des checks. Distinguer commit de PR, eventuel merge ref de CI et commit reellement fusionne. En cas de squash, la non-ascendance d'un ancien HEAD ne prouve pas l'absence du correctif : utiliser la PR, son commit de fusion, le diff semantique et les tests. Ne pas transposer les erreurs observees sur `7ccb31d...` au candidat `15bb69c...` sans les reexecuter.

Produire une matrice `constat -> fichier/symbole -> PR -> commit integre -> test enregistre -> execution -> campagne physique`. Un statut inconnu reste inconnu. Ni fermeture de PR ni fichier de test present ne valent execution reussie.

**Cloture :** arbre combine propre, provenance conservee, tous les tests obligatoires effectivement decouverts par CTest ; aucune base ou campagne implicite.

#### Lot 1 - Solveurs, gradients et diagnostics : #1092 + #1096 + gardes I01

Dans `gmres_solver.h`, mutualiser la construction des operateurs de travail et de certification. Tous les retours CONVERGED doivent etre compatibles avec le residu de l'operateur FP64 original. Tester GMRES et FGMRES separement, y compris avec budget insuffisant. Si b est non nul, garder la reference annoncee au RHS ; pour b=0 definir une reference explicite et stable, commune au critere et au diagnostic. Un grand x0 ne doit pas relacher silencieusement le seuil.

Dans `gradient.h`, `cell_gradient` et les laplaciens, propager une politique de frontiere et les BC, pas seulement modifier une valeur par defaut. Un fantome de gradient nul n'est pas une valeur analytique generique. Tester tetraedres de bord, champs affines avec BC coherentes, stencils deficits et cas 2D declares. Conserver le repli QUICK avec les vrais voisins et la normalisation geometrique. Rejouer les produits scalaires owned/ghost en MPI, y compris partitions de tailles differentes.

Combiner les changements de #1096 (thermodynamique, temps, normes et cas 1x1) sans perdre la fin de CMake ni H06-H08 de R4. Enregistrer `test_audit_numerics` apres resolution de sa dependance #1092 et verifier sa presence dans CTest/CI. Le commentaire public annonce qu'il avait ete ajoute sans enregistrement [W5] : l'etat final doit etre controle, pas suppose.

**Cloture :** tous les contre-exemples G01-G09, NaN/Inf/overflow, RHS nul/non nul, plusieurs echelles et appels consommateurs passent au SHA combine. L'ancien audit #1092 montre quoi retester, pas un verdict de la nouvelle version.

#### Lot 2 - F03 : fermer le contrat pression-vitesse, pas seulement son flux

Cibles : `steady_incompressible_solver.h`, assemblage pression-vitesse et Schur/preconditionneurs ; comparer #1078 au prototype R4.

Avec U physique, la reconstruction doit respecter `H/A=U+R*G*p` et `phi=rho*(I(U+R*G*p).S - T(p))`. Par divergence, la meme compensation doit figurer dans le bloc de continuite. Ne pas remplacer le champ d'entree sans modifier son contrat ; ne pas corriger seulement le post-traitement ; ne pas supprimer le traitement anti-damier.

**Cloture :** pression affine/vitesse constante, invariance de jauge, pression en damier, identite bloc/flux avec oracle independant, tests des frontieres et MMS vitesse-pression. Puis cas Couette/Poiseuille/cavite et comparaison des variantes effectivement proposees. Inclure viscosite variable et combinaison avec le stress v2. Le PCD qui emploie la viscosite moleculaire alors que l'operateur emploie mu_eff est a qualifier en robustesse ; ce n'est pas une preuve que l'operateur omet mu_eff.

#### Lot 3 - F06/F08 : diffusion et gradients communs aux equations reelles

Cibles : `finite_volume_transport.h`, `laplacian.h`, gradients du noyau et gradients internes du solveur. La correction doit etre consommee par momentum, energie et turbulence, pas uniquement par un test du noyau.

Definir une decomposition coherente `S=alpha*d+T`, une contribution a deux points implicite et une correction transverse reconstruite. Employer la meme orientation, les memes coefficients et les memes gradients dans l'assemblage et le flux restitue ; verifier le bilan avec une reconstruction independante. Definir Dirichlet, Neumann et Robin et les interfaces a coefficients discontinus. L'interpolation geometrique corrigee ne traite pas a elle seule la non-orthogonalite.

**Cloture :** constantes/affines, champ lisse non polynomial, raffinement, rotations et changement d'unites, maillages etires/cisailles/polyedriques, bilans locaux. Mesurer interieurs et frontieres separement. R4 conserve deux ordres grossiers d'environ 1,77 sous la gate 1,8 ; les ordres fins 1,92/1,97 ne rendent pas ces gates retroactivement vertes. Examiner le regime asymptotique et le traitement des frontieres ; toute modification de protocole doit etre justifiee avant d'annoncer la qualification.

#### Lot 4 - Preserver et publier les ameliorations du lot v2

Fichiers centraux : `energy_solver.h`, `production_energy_stage.h`, `production_radiation_s2s_stage.h`, `cht_solver.h`, `sst_solver.h`, `scalar_transport_gradient.h`, `momentum_viscous_stress.h`, `steady_incompressible_solver.h`, applications et HDF5.

Conserver `Su+Sp*T_final` dans les bilans ; lambda par cellule jusqu'aux BC et aux flux ; T absolue admissible ; parois radiatives calculees aux faces ; couverture complete des parois omega ; residus apres bornage k/omega ; premier predicteur et retry observes sur les vrais tableaux consommes. Pour F11 conserver la contribution transposee dans le predicteur, le residu final et les efforts/bilans, sans double comptage.

La pression SST de v2 est `p*=p_statique+2*rho*k/3`. Declarer ce contrat dans le cas, les BC, les sorties et le restart ; convertir explicitement avant de publier pression statique/Cp. Tester forts gradients de viscosite aux murs, variations/contrastes, traitement explicite, maillages skew et criteres couples. La MMS actuelle du stress annule la derivee de mu au bord ; elle ne qualifie pas les murs generaux.

**Cloture :** 43 cibles historiques et nouveaux tests de combinaison rejoues au commit integre, puis qualification physique propre a chaque methode. Aucun compteur de l'ancien arbre n'est transfere automatiquement.

#### Lot 5 - Rendre SST utilisable depuis un cas et alimenter l'energie

Cibles : `apps/production_solver/main.cpp`, `production_turbulence_contract.h`, `sst_production_coupling.h`, schema C++/Python et HDF5.

Construire l'etat SST depuis les metadonnees reelles et les champs initialises/repris : k, omega, choix de modele et constantes, BC, parois, valeurs omega par face et convention de pression. Recalculer distance, deformation, nu_t et mu_eff au bon moment avant le premier predicteur. Le petit fixture C++ v2 apporte des preuves du cycle de vie ; il ne remplace pas la creation depuis un fichier utilisateur.

Declarer et persister `Pr_t`, lire uniquement `nu_t` accepte, puis transmettre `lambda_eff=lambda+rho*cp*nu_t/Pr_t` a l'energie et au S2S. La formule recoit une viscosite CINEMATIQUE ; ne pas remultiplier par rho une viscosite dynamique. Commencer par energie passive sur ecoulement accepte ; les retroactions viennent dans une boucle explicitement convergee.

**Cloture :** cas CLI non fabrique a la main, flux autoritaire non nul, effet thermique mesurable d'une variation de Pr_t/nu_t, limite laminaire, bilans et rejet sans mutation d'un etat invalide. Le selecteur RANS general reste refuse tant que ses criteres de promotion ne sont pas satisfaits.

#### Lot 6 - CHT multi-regions, contact et pilotage global

Definir regions, maillages, materiaux, etats fluides/solides et interfaces versionnes dans `case_schema.h`, le schema Python et HDF5. Le solveur solide multizone n'est pas le pilote fluide-solide general. D2 ne trouve pas d'appel applicatif a `solve_two_region_cht` dans `apps/` [D2, TH-02].

Creer la boucle regionale : porteur sur fluides, energie sur fluides/solides, correspondance des interfaces et powers opposes. Conserver les controles v2 : faces adjacentes, normales opposees, meme pas physique, sources valides, impermeabilite et etat ancien fige. Etendre le noyau CHT aux conductivites par cellule : son interface ne doit pas continuer d'utiliser deux scalaires lorsque l'energie utilise lambda_eff.

Une resistance de contact surfacique appartient a une INTERFACE. Employer `q''=(T_P1-T_P2)/(d1/lambda1+R''+d2/lambda2)`, et deux temperatures de surface telles que `T_gamma1-T_gamma2=R''*q''`. Une temperature de surface unique des deux cotes est incompatible avec un saut non nul. Un materiau central touche par deux voisins doit pouvoir avoir deux resistances distinctes.

**Cloture :** aller-retour C++/Python/HDF5, vrai cas CLI a deux regions, parite des flux, bilan individuel de chaque region et bilan global, deux contacts differents a trois regions, limites Rc=0 et nu_t=0, erreurs regionales provoquant restauration atomique. Les interfaces non conformes ne sont pas couvertes implicitement : les refuser jusqu'a un traitement declare.

#### Lot 7 - Transaction multiphysique, reprise et temps physique

Une iteration ne doit publier l'etat qu'apres satisfaction des gates du porteur, SST, energie, CHT et rayonnement. Rejeter/restaurer U,p,phi,k,omega,nu_t,T par region, coefficients derives, radiosites/flux appliques, interfaces, moniteur et historiques. Injecter un echec apres chaque sous-etape et comparer a un essai propre utilisant les memes relaxations acceptees.

D2 confirme que DAT v4 ne sauvegarde pas le flux de face autoritaire. La sauvegarde de U/p/T/k/omega ne prouve donc pas la reprise d'une trajectoire generale. Completer version, maillage, conventions, champs de face et historiques temporels ; valider dans un etat temporaire avant publication. Autoriser un champ derive reconstruit seulement avec un contrat deterministe et une preuve d'equivalence.

Introduire un pilote de temps physique distinct des iterations non lineaires. Garder T^n et, pour BDF2, les etats et pas precedents pendant les corrections internes ; avancer le temps seulement apres acceptation. Les executables v2 refusent explicitement les regimes temporels non implementes : ne pas enlever cette garde avant le pilote et sa preuve.

**Cloture :** calcul continu vs interrompu/repris sur plusieurs pas, flux/coefficients non uniformes, erreurs de checkpoint sans mutation, bilan accumulatif et etude d'ordre temporel.

#### Lot 8 - Rayonnement general, autres fermetures et campagnes physiques

S2S : conserver la temperature de surface de la BC et distinguer flux applique de flux recalcule. Tester reciprocite `Ai*Fij=Aj*Fji`, visibilite/occultation, auto-vue d'une face plane, fermeture uniquement en enceinte fermee et echanges avec l'environnement en enceinte ouverte. Subdiviser les parois non isothermes ; une renormalisation independante des lignes peut detruire la reciprocite.

P1/DOM/Rosseland : un noyau present ne vaut pas raccordement production. Ajouter variables radiatives, coefficients d'absorption/diffusion, BC, sources volumiques avec signe/coherence energie, et iterations couples. Qualifier selon domaine optique : plaque/enclosure, convergence spatiale et angulaire pour DOM, bilan source-volume/flux-frontiere. La correction auxiliaire 4*pi de #1079 n'effectue pas ces raccordements.

Clore separement F02 (Wilcox complet), F18 (EOS), F22-F24 (WALE/compressibilite/Kato-Launder/QCR), F27 (VOF), A03 (meshio) et A05 (qualification). Ne pas ouvrir un modele non qualifie sous pretexte que la variante SST locale fonctionne.

Matrice minimale : MMS des DEUX transports k/omega sur au moins trois maillages, canal turbulent avec reference figee, y+ et frottement, plaque plane puis cas aerodynamique ; canal chauffe avec Nu et Pr/Re traces ; resistance composite/CHT ; bilan radiatif et references adaptees. Publier normes, ordres, maillages, proprietes, BC, versions et SHA. Ajouter MPI a partitions inegales ; GPU seulement apres une campagne propre. Marquer PASS, FAIL, SKIPPED et NOT_RUN distinctement.

### 13.6 Registre de raccordement des 56 constats historiques

Cette table complete les statuts R4 sans les requalifier. Elle indique la preuve la plus recente retrouvee, et le critere qui manque pour l'arbre d'integration. **Aucune ligne n'est certifiee aujourd'hui sur le dernier master.**

| ID | Sujet | Etat documentaire/local le plus recent | Action de cloture |
|---|---|---|---|
| F01 | Production de omega non limitee dans SST-2003 | v2 : source limitee reprise de #1094. Page publique Draft. | Rejouer les sources limites puis MMS k/omega sur le commit combine. |
| F02 | La fermeture nommee komega2006 omet le limiteur de Wilcox | R4 partiel ; #1096 le vise ; hors sous-ensemble thermique v2. | Completer beta 3D/Wilcox et ses invariants ; ne pas confondre viscosite limitee et diffusivite. |
| F03 | La branche COUPLED utilise une vitesse physique comme H/A dans le flux de pression | #1078 8c00dfd... repertorie, exclu du lot v2 ; prototype R4. | Integrer flux U/HbyA et bloc de continuite ensemble, puis Schur et cas physiques. |
| F04 | Le couplage CHT transitoire avance le temps a chaque iteration d'interface | v2 : extrait thermique #1096 combine avec H04/H05. | Conserver T^n fixe et reexecuter les six cas dt/relaxation apres integration. |
| F05 | Moyenne harmonique non ponderee sur cellules non uniformes | Annonce present dans la base du lot v2 ; non relu au master courant. | Verifier les resistances ponderees avec contrastes et interfaces asymetriques. |
| F06 | La diffusion du solveur n'utilise pas le laplacien non orthogonal corrige | Prototype R4 seulement ; deux gates de convergence echouent. | Raccorder diffusion non orthogonale matrice/flux/bilan/BC dans le transport reel. |
| F07 | Interpolation LINEAR centrale 1/2 non exacte sur grille etiree | Annonce present dans la base ; aucune nouvelle preuve du master. | Rejouer champs affines et grilles etirees sur chaque consommateur. |
| F08 | Le poids du gradient de Gauss confond distance et projection | R4 partiel ; le gradient SST v2 respecte ses BC sans correction generale de skewness. | Unifier politiques/BC ; conserver les erreurs de bord et les etudes non affines. |
| F09 | Le flux convectif sortant reconstruit ne correspond pas a la matrice | Annonce present dans la base du lot v2. | Verifier entree/sortie et identite des flux assembles/reconstruits. |
| F10 | Lecture d'un gradient vide pendant la reconstruction de flux | Annonce present dans la base du lot v2. | Rejouer tous les schemas sous ASan/UBSan sur le commit combine. |
| F11 | Diffusion de vitesse scalaire incomplete lorsque la viscosite varie | Nouveau statut LOCAL : v2 raccorde stress transpose au predicteur, residu et bilan. | Publier/revoir cette integration, tester murs, mu discontinu, skewness et convergence externe. |
| F12 | Helper momentum: convection integree melangee a des termes volumiques | Correctif helper R4 ; pas de nouvelle trace v2 de ce fichier. | Verifier sa presence dans la cible puis tests dimensionnels et de raccordement. |
| F13 | Helper thermique: convection en W et conduction en W/kg | Correctif helper R4 ; pas de nouvelle trace v2 de ce fichier. | Verifier les taux K/s et comparaison independante avec equation energie. |
| F14 | Integrale d'enthalpie tabulee incorrecte hors table | R4 corrige ; vise par #1096, non couvert par le lot thermique selectionne. | Rejouer integrale tabulee, extrapolation, antisymetrie et h(T,T)=0. |
| F15 | L'enthalpie n'integre pas le Cp effectivement borne | R4 corrige ; vise par #1096. | Verifier dh/dT=cp effectivement borne, notamment aux points de changement de borne. |
| F16 | Cp par defaut et conductivite de l'air errones dans une surcharge transport | R4 corrige ; vise par #1096. | Conserver la correction des proprietes air et les nouvelles exigences de cas explicite. |
| F17 | Deux IdealGasEOS renvoient des derivees incompatibles | R4 corrige ; vise par #1096. | Verifier les deux EOS et leurs derivees sur l arbre combine. |
| F18 | Identites incompatibles dans IncompressibleEOS pour beta non nul | Toujours ouvert dans les preuves disponibles. | Separer dilatation thermique, compressibilite et compressibilite artificielle ; tests derives/unites. |
| F19 | BDF2 a pas variable incorrect dans le pilote generique double temps | R4 corrige ; vise par #1096 ; pas de pilote transitoire applicatif qualifie. | BDF2 variable : trajectoire polynomiale, changements de dt et reprise des historiques. |
| F20 | Une iteration implicite peut accepter un NaN comme convergence | R4 corrige ; vise par #1096. | Rejets NaN/Inf a toutes les sorties implicites et budgets epuises. |
| F21 | Les integrateurs low-storage ne mettent a jour que la composante 0 | R4 corrige ; vise par #1096. | Verifier toutes les composantes, dont vecteurs et champs multi-composantes. |
| F22 | Une fonction nommee WALE calcule la loi de Smagorinsky | Ouvert : nom WALE auxiliaire ambigu. | Retirer/renommer le helper scalaire ; utiliser la vraie fermeture tensorielle et migrer les appelants. |
| F23 | Correction de compressibilite turbulente non adimensionnelle | Ouvert : correction compressible dimensionnellement incoherente. | Choisir une variante documentee et un Mach turbulent adimensionnel ; invariance d unites. |
| F24 | Correction Kato-Launder ajoutee avec des dimensions incompatibles | Ouvert : Kato-Launder/QCR a distinguer. | Production dimensionnellement correcte et fermeture QCR tensorielle separee ; stagnation/rotation. |
| F25 | Surcharge historique SST sans invariant de deformation | v2 reprend #1081 pour retirer la surcharge SST ambigue. | Rejouer les appelants, interdire silencieusement l ancienne signature. |
| F26 | Facteur pi manquant dans une source P1 auxiliaire | v2 reprend #1079 : normalisation de la source P1 auxiliaire. | Conserver 4*pi et verifier signes/unites lors du futur raccordement volumique. |
| F27 | Courbure VOF: derivees de composantes sur le meme axe | #1083 reperee ; non incluse dans v2 ; HEAD courant non verifie. | Plan/cercle/sphere, signe, rotations, convergence et force capillaire du chemin VOF. |
| F28 | Interface CHT: distances et conductivites de deux contrats differents | v2 combine #1097, deux validations de #1095 et corrections thermiques. | Ne pas perdre validation avant mutation, distance normale, coherence k, impermeabilite et dt commun. |
| F29 | Resistance thermique de zone suppose des demi-distances egales | #1084 documentee fusionnee dans base 038e3a... ; page web contradictoire. | Ne pas reappliquer en bloc ; reexecuter Q=8000/39 et distinguer contact par interface encore absent. |
| F30 | Le producteur de facteurs de vue viole le contrat du solveur S2S | v2 reprend #1086 : contrat producteur/consommateur corrige localement. | Visibilite, surfaces concaves, enceinte ouverte et quadrature convergente restent a qualifier. |
| F31 | max_imbalance est un bilan global, pas le maximum local | v2 reprend #1080 : global signe/L1 et maximum local separes. | Auditer les reductions distribuees, debordements et normalisations. |
| F32 | Statut lineaire rapporte avant la sous-relaxation du champ retourne | v2 reprend #1088 corrigee ; diagnostics avant/apres presents ; garde I01 manquante. | Garde finitude commune et emploi du residu du champ accepte dans la convergence externe. |
| F33 | Bornage haute resolution different entre assemblage et reconstruction | Correction annoncee dans la base, regression #1089 reprise v2. | Rejouer identite assemblage/flux, schemas et bornage active/desactive. |
| F34 | Plusieurs tests recodent la formule fautive au lieu d'un oracle independant | v2 retient #1100 et durcit les comparateurs ; revue globale des oracles encore partielle. | Ne pas empiler #1093/#1098 ; tests de mutation, NaN/Inf/overflow, references independantes. |
| G01 | Politiques de frontière incohérentes entre LS et WLS | R4 corrige localement ; #1092 15bb69c... repertoriee mais hors lot v2 ; audit ancien sur 7ccb31d... seulement. | Propager politiques et BC dans WLS, cell_gradient, laplaciens ; tetraedres de rang insuffisant. |
| G02 | CG interne : convergence dépendante de l’échelle et résidu mal normalisé | R4 corrige localement ; #1092 15bb69c... repertoriee mais hors lot v2 ; audit ancien sur 7ccb31d... seulement. | Certification vraie FP64 sur Laplace multi-echelle, avec/sans preconditionneur. |
| G03 | Précision mixte : le vrai résidu est calculé avec l’opérateur FP32 | R4 corrige localement ; #1092 15bb69c... repertoriee mais hors lot v2 ; audit ancien sur 7ccb31d... seulement. | Tester GMRES et FGMRES separement : travail FP32, certification FP64 et absence de faux succes. |
| G04 | BiCGStab : seuils de rupture absolus et diagnostic nul sur échec | R4 corrige localement ; #1092 15bb69c... repertoriee mais hors lot v2 ; audit ancien sur 7ccb31d... seulement. | Seuils de rupture relatifs et diagnostic correct sur vraie divergence BiCGStab. |
| G05 | GMRES : second membre nul refusé avec un itéré initial non nul | R4 corrige localement ; #1092 15bb69c... repertoriee mais hors lot v2 ; audit ancien sur 7ccb31d... seulement. | b=0 et b non nul avec grand x0 ; reference relative constante et diagnostic coherent. |
| G06 | GMRES : produits scalaires Arnoldi incompatibles avec la politique de réduction | R4 corrige localement ; #1092 15bb69c... repertoriee mais hors lot v2 ; audit ancien sur 7ccb31d... seulement. | Produits scalaires sur owned seulement et reductions MPI avec tailles locales inegales. |
| G07 | QUICK : repli linéaire sur un stencil vide | R4 corrige localement ; #1092 15bb69c... repertoriee mais hors lot v2 ; audit ancien sur 7ccb31d... seulement. | Repli QUICK avec voisins disponibles, exactitude affine et domaine 2D declare. |
| G08 | QUICK : base quadratique non normalisée et rang dépendant des unités | R4 corrige localement ; #1092 15bb69c... repertoriee mais hors lot v2 ; audit ancien sur 7ccb31d... seulement. | Normalisation geometrique, rang/rotation, plusieurs echelles et stencils difficiles. |
| G09 | Oracle ExactSchur : NaN traité comme un second membre nul | R4 corrige localement ; #1092 15bb69c... repertoriee mais hors lot v2 ; audit ancien sur 7ccb31d... seulement. | Schur exact : RHS zero vs NaN/Inf, oracle independant. |
| A01 | Bilan energetique global susceptible de compenser des erreurs locales | v2 ameliore bilans energie avec Su+Sp*T et rejets locaux ; pas de cloture generale. | Etendre a accumulation, interfaces, rayonnement et chaque equation sans compensation locale. |
| A02 | Compilation CPU sans MPI interrompue par une inclusion MPI non conditionnelle | R4 corrige ; #1096 vise la garde MPI. | Compiler une configuration sans MPI complete ; un build avec MPI ne la valide pas. |
| A03 | Dependance meshio non declaree pour les chemins de conversion qui l'exigent | Integration de meshio dans packaging non reverifiee. | Installation propre plus conversions OpenFOAM ; pas de skip silencieux. |
| A04 | Une diffusivite NaN peut etre transformee en zero dans le flux | Correctif local R4 ; presence dans cible finale a verifier. | Rejets NaN sur coefficients, geometrie et flux avant mutation. |
| A05 | Qualification inachevee de familles numeriques et physiques | Qualification globale toujours ouverte. | Campagnes exact-SHA avec erreurs, ordre, bilans, et statut explicite des tests ignores. |
| H01 | Normes instables et fausses convergences par overflow/underflow | R4 corrige ; norme robuste visee par #1096, pas certifiee par lot v2. | Rejouer changements d echelle et produits/reductions potentiellement debordants. |
| H02 | Newton retourne un echec apres son dernier pas pourtant convergent | R4 corrige ; vise par #1096. | Newton au dernier pas convergent, dernier pas non convergent et NaN. |
| H03 | Diagonale matrix-free artificielle et coefficient NaN accepte | R4 corrige ; vise par #1096. | Matrix-free : diagonale effective, coefficient NaN et equivalence avec matrice assemblee. |
| H04 | Tolerance relative traitee comme absolue dans le solveur scalaire 1x1 | v2 reprend le correctif relatif 1x1 combine a #1088. | Rejouer grands/petits RHS, zero RHS et contrat predicteur/champ accepte. |
| H05 | Diagnostic energetique de frontiere calcule avec la mauvaise distance | v2 reprend la metrique normale du bilan thermique. | Verifier les memes distances dans matrice, flux et audit, sans reutiliser le RHS comme oracle. |
| H06 | Poids Gauss non fini ou dependant des unites | Correctif R4 local ; pas de preuve d integration cible. | Porter/rejouer geometries extremes, non finies et degenerees. |
| H07 | Flux entrant zero-gradient omis sans ancien champ convecte | Correctif R4 local ; presence dans cible non verifiee. | Rejouer entree ZERO_GRADIENT sans ancien champ, convection bornee/non bornee. |
| H08 | Diagnostic de residu nul sur des donnees NaN | Helper R4 corrige, mais un autre reducteur v2 masque encore NaN (I01). | Mutualiser et appliquer les gardes a tous les diagnostics et a la normalisation du RHS. |

### 13.7 Definition de fini pour le prochain commit d'integration

Un lot est pret seulement si sa formule est documentee, son appelant reel l'utilise, ses tests negatifs et references independantes sont enregistres ET executes, ses diagnostics portent sur les champs acceptes, et sa preuve correspond au commit final. Une CI d'une branche ne remplace pas la CI de la combinaison. Les campagnes physiques sautées ne sont pas des succes.

Ordre recommande : **provenance et harnais -> #1092 / diagnostics -> #1096 et petits correctifs -> comparaison/integration du lot v2 -> #1078 et verification combinee -> diffusion/gradients -> SST depuis cas et chaleur turbulente -> CHT global -> transaction/reprise/temps -> qualification elargie**. Certains travaux de schema peuvent avancer independamment, mais pas promouvoir un solveur avant ses gates.

### 13.8 Sources et reproductibilite de cette revue

- **D1** : `CFDX_audit_consolide_correctifs_testes_2026-10-24.md`, revision R4, conservee avant ajout de cette section. Patch de 54 fichiers ; empreinte du document original : `afebd987fa1c23b194fd85636209fcd4aeec855a5db308c19b8439371da91ae9`.
- **D2** : `CFDX_integration_thermique_rayonnement_turbulence_2026-10-24.md`, revision 2, relevé Git 24/10 16:46 UTC. Sections 1-2 (SHA/selection), 4-7 (contrats/limites), 8 (tests), annexe A (patch). Empreinte du document relu : `8e42c4e52224b2ab71615bedfab6e4e17d81775e480533575f45f0089446347e` ; empreinte du patch : `1cef46e6ad105f0b1818bceeffeff313f5af1bacb8af4546b6a6fa9dbde1ebf3`.
- **D3** : `CFDX_PR1092_audit_7ccb31d8_2026-10-24.md`, verdict limite a `7ccb31d8e2fb5882acf0fe1beb9ecfbb616bf227`, anterieur a `15bb69c...`.
- **W1** : https://github.com/stevendaix/cfdx_dev/pull/1047 (fusion affichee).
- **W2** : https://github.com/stevendaix/cfdx_dev/issues/1067 (gates SST et campagnes sautees).
- **W3** : https://github.com/stevendaix/cfdx_dev/pull/1068 (draft documentaire dans la capture accessible).
- **W4** : https://github.com/stevendaix/cfdx_dev/pull/1084 (capture ouverte, contradictoire au releve Git plus recent de D2).
- **W5** : https://github.com/stevendaix/cfdx_dev/pull/1094 (source SST locale et commentaire de #1096 ; test numerique annonce non enregistre dans ce commentaire historique).

Sources web reconsultees le 25/10/2026 ; pas de certification de leur fraicheur par rapport au Git courant. Aucun push, merge ou nouveau test complet CFDX n'a ete effectue dans cette revision. La sonde suivante utilise des conteneurs minimaux et ne doit pas etre presentee comme un solveur CFDX execute.

#### Sonde I01 - source compilee (lambda extraite du patch v2 et garde proposee)

```cpp
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
namespace cfdx { namespace core {
struct Vector {
    std::vector<double> data;
    std::size_t size() const { return data.size(); }
    double operator()(std::size_t i) const { return data.at(i); }
};
}}
struct Matrix {
    std::vector<std::uint32_t> offsets, columns;
    std::vector<double> values;
    const std::uint32_t* row_offsets_data() const {return offsets.data();}
    const std::uint32_t* columns_data() const {return columns.data();}
    const double* values_data() const {return values.data();}
};
struct Equation { Matrix matrix; cfdx::core::Vector rhs; };
double original(const Equation& equation,const cfdx::core::Vector& state) {
        const auto physical_residual_inf = [&](const cfdx::core::Vector& state) {
            double norm = 0.0;
            for (std::size_t i = 0; i < state.size(); ++i) {
                double ri = -equation.rhs(i);
                const auto row_begin = equation.matrix.row_offsets_data()[i];
                const auto row_end = equation.matrix.row_offsets_data()[i + 1];
                for (std::uint32_t k = row_begin; k < row_end; ++k)
                    ri += equation.matrix.values_data()[k] *
                          state(equation.matrix.columns_data()[k]);
                norm = std::max(norm, std::abs(ri));
            }
            return norm;
        };
    return physical_residual_inf(state);
}
// Guard prototype only; not applied to CFDX. Invalid arithmetic cannot certify a residual.
double guarded(const Equation& eq,const cfdx::core::Vector& state) {
    const auto bad=std::numeric_limits<double>::infinity();
    if (state.size()!=eq.rhs.size() || eq.matrix.offsets.size()!=state.size()+1 ||
        eq.matrix.columns.size()!=eq.matrix.values.size()) return bad;
    for (double x:state.data) if(!std::isfinite(x)) return bad;
    double result=0;
    for(std::size_t i=0;i<state.size();++i) {
        if(!std::isfinite(eq.rhs(i))) return bad;
        auto lo=eq.matrix.offsets[i],hi=eq.matrix.offsets[i+1];
        if(lo>hi || hi>eq.matrix.values.size()) return bad;
        double ri=-eq.rhs(i);
        for(auto k=lo;k<hi;++k) {
            if(eq.matrix.columns[k]>=state.size() || !std::isfinite(eq.matrix.values[k])) return bad;
            double term=eq.matrix.values[k]*state(eq.matrix.columns[k]);
            if(!std::isfinite(term)) return bad;
            ri+=term;
            if(!std::isfinite(ri)) return bad;
        }
        result=std::max(result,std::abs(ri));
    }
    return result;
}
int main() {
    const double nan=std::numeric_limits<double>::quiet_NaN();
    const double huge=std::numeric_limits<double>::max();
    struct Test {std::string name;Equation eq;cfdx::core::Vector x;bool invalid;double expected;};
    const std::vector<Test> cases={
        {"finite_exact",{{{0,1},{0},{2}},{{4}}},{{2}},false,0},
        {"finite_nonzero",{{{0,1},{0},{2}},{{1}}},{{2}},false,3},
        {"nan_field",{{{0,1},{0},{1}},{{1}}},{{nan}},true,0},
        {"nan_coefficient",{{{0,1},{0},{nan}},{{1}}},{{1}},true,0},
        {"nan_rhs",{{{0,1},{0},{1}},{{nan}}},{{1}},true,0},
        {"finite_inputs_overflow_cancel",{{{0,2,3},{0,1,1},{huge,-huge,1}},{{0,2}}},{{2,2}},true,0}
    };
    int failures=0;
    std::cout<<"case\toriginal\tguarded\tguard_pass\n"<<std::setprecision(17);
    for(const auto& c:cases) {
        const double a=original(c.eq,c.x),b=guarded(c.eq,c.x);
        const bool ok=c.invalid?std::isinf(b):(b==c.expected);
        failures+=!ok;
        std::cout<<c.name<<'\t'<<a<<'\t'<<b<<'\t'<<ok<<'\n';
    }
    return failures?1:0;
}
```

#### Sonde I01 - sortie obtenue

```text
case	original	guarded	guard_pass
finite_exact	0	0	1
finite_nonzero	3	3	1
nan_field	0	inf	1
nan_coefficient	0	inf	1
nan_rhs	0	inf	1
finite_inputs_overflow_cancel	0	inf	1
```
