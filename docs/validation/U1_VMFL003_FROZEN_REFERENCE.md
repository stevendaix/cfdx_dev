# U1 — VMFL003 : référence gelée du canal turbulent RANS

**Statut : RÉFÉRENCE GELÉE AVANT EXÉCUTION** (gouvernance #118 : la référence, les quantités d'intérêt et les tolérances sont fixées avant toute interprétation ; aucun relâchement a posteriori).

- Issue de cas : #882 (U1) — Parent : #118 (autorité de qualification) — V&V équations par modèle : #473 — Wall-distance : #466.
- Fiche de cas : `TUR-CHANNEL_QUALIFICATION_CASE.md` (PLANNED/BLOCKED — ce gel est le préalable documentaire du déblocage).
- Références numériques : `LEVEL_B_REFERENCES.md` (Moser–Kim–Mansour DNS, Re_tau = 178.12, sous-ensemble déterministe déjà versionné dans le repo).

## 1. Configuration de référence gelée (étage 1 : canal plan 2D)

| Grandeur | Valeur gelée | Justification |
|---|---|---|
| Géométrie | Canal plan, demi-hauteur h = 0.01 m, L/h = 100 (section développée L/h ≥ 50 mesurée, amont ignoré) | canal plan = première étape de #882 ; 3D/axisymétrique (#445) en variante ultérieure |
| Fluide | Air, rho = 1.2 kg/m^3, nu = 1.5e-5 m^2/s, propriétés constantes | gel explicite |
| Régime | Développé, RANS à deux équations (k-epsilon standard primaire ; SST en variante) | #882 : « RANS à deux équations » |
| Vitesse débitante | U_b = 4.2 m/s | donne Re_b(h) = U_b·h/nu = 2800 |
| Re_bulk (base h) | 2800 | cohérent avec le point DNS MKM Re_tau = 178.12 |
| Nombre de friction visé | Re_tau ≈ 178 | point DNS disponible dans le repo |

## 2. Quantités d'intérêt (QoI) et valeurs de référence

Toutes les valeurs DNS sont **dérivées du sous-ensemble MKM gelé du repo** (jamais codées en dur dans les assertions ; le driver les relit depuis le fichier de référence).

| QoI | Valeur de référence gelée | Source |
|---|---|---|
| u_tau | 0.26718 m/s | Re_tau·nu/h, Re_tau = 178.12 (MKM) |
| U_b+ (= U_b/u_tau) | 15.72 | Re_b/Re_tau |
| C_f = tau_w / (rho·U_b^2/2) | 8.09e-3 (= 2/U_b+^2) | DNS MKM dérivée |
| tau_w | 8.57e-2 Pa | rho·u_tau^2 |
| dp/dx (développé) | 8.57 Pa/m | tau_w/h (bilan de quantité de mouvement exact du canal plan) |
| Delta_p sur L/h = 50 mesuré | 4.28 Pa | (tau_w/h)·L |
| Profil U+(y+) | sous-ensemble MKM chan180 | fichier de référence du repo |

**Référence secondaire (recoupement ingénierie, gelée avec sa définition)** : corrélation de Dean C_f = 0.073·Re_H^(-1/4) avec Re_H = U_b·(2h)/nu = 5600 → C_f = 8.44e-3 (écart DNS ≈ 4 %). La définition de Re (hauteur totale 2h) fait partie du gel. Variante conduit strict VMFL003 (3D-axisym., #445) : Colebrook–White lisse (forme de Haaland) sur Re_D — à geler dans un addendum au moment de la variante, pas rétroactivement.

## 3. Gates de réussite (gelées avant interprétation — méthode #118)

1. **Accord quantitatif** : C_f (donc tau_w et dp/dx, cohérents entre eux) dans **±5 % en relatif** de la valeur DNS dérivée. Le recoupement Dean doit être rapporté mais n'est pas un gate (écart de définition documenté).
2. **Conservation (gate indépendante)** : bilan de masse relatif < 1e-10 (fermeture indépendante du solveur) ; bilan de quantité de mouvement : |-dp/dx - tau_w/h| / (tau_w/h) < 1 % sur la section développée.
3. **Stabilité et positivité** : k > 0, epsilon > 0 (omega > 0, nu_tilde > 0 selon le modèle) à toutes les itérations ; les compteurs N11 (#480 : fallbacks nearest_normal, clamps de positivité) doivent rester à zéro sur la config gelée, sinon rapportés explicitement — un clamp silencieux invalide le PASS.
4. **Profil** : écart relatif maximal au plateau logarithmique (30 < y+ < 100, cellules disponibles) < 10 % vs MKM.
5. **Étude de maillage** : 3 niveaux de raffinement pariétal (p.ex. y+ ≈ 60 / 30 / < 1 avec traitement cohérent), QoI dans la tolérance au niveau de production ; ordre observé rapporté, pas gating sans justification.
6. **Régime pariétal** : y+ de la première cellule et régime (wall functions vs low-Re) documentés dans les artefacts ; changement de régime = nouvelle configuration, nouvelle exécution, pas de mélange.
7. **Reproductibilité** : driver + configuration enregistrés, logs et QoI conservés (artefacts CI), exécution sur le HEAD exact de master (A0 #879).

## 4. Dépendances et blocages

- A0 (#879) : cartographie du driver turbulent exécutable — doit rester satisfaite au HEAD de lancement.
- CI validation (#870) : la promotion U1 exige la CI verte au HEAD exact ; le gel de référence peut être mergé indépendamment (documentaire).
- Wall-distance (#466) : requis pour SST uniquement ; le timeout du job `wall-distance` (#870) ne doit pas masquer un échec de convergence du service.
- Gaps de capability rencontrés à l'exécution → tracés vers #422 ; V&V équation par équation → #473.

## 5. Ce que ce document n'est pas

Ce fichier ne déclare **aucun PASS**. Il ne contient aucune valeur mesurée par CFDX. Toute valeur CFDX obtenue après exécution sera comparée à ces gates gelés, sans ajustement de tolérance (constitution #118 : aucune tolérance élargie, aucun oracle substitué, aucun PASS à itérations fixes).
