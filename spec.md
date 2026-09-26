# CFDX — Spécification technique complète

**Version : 0.7**
**Statut : Architecture cible / spécification de développement**
**Projet : CFDX — High-Performance General-Purpose CFD Framework**

---

# 1. Vision

CFDX est un framework CFD généraliste destiné au calcul scientifique haute performance sur :

* CPU ;
* GPU ;
* systèmes hybrides CPU/GPU ;
* clusters MPI ;
* environnements HPC ;
* problèmes de grande taille dépassant éventuellement la mémoire GPU.

L'objectif n'est pas de reproduire OpenFOAM ligne par ligne.

CFDX doit reprendre les concepts numériques éprouvés de la méthode des volumes finis tout en séparant radicalement :

1. représentation des données ;
2. maillage ;
3. géométrie ;
4. champs ;
5. opérateurs FVM ;
6. algèbre linéaire ;
7. physique ;
8. exécution CPU/GPU ;
9. gestion mémoire ;
10. parallélisme ;
11. stockage ;
12. orchestration Python.

Le principe fondamental est :

> **La physique ne doit pas connaître le backend d'exécution.**

Une équation physique doit pouvoir être exécutée :

* sur CPU ;
* sur GPU ;
* en mode GPU out-of-core ;
* en MPI ;

sans modifier sa formulation.

---

# 2. Principes architecturaux

## 2.1 Séparation des responsabilités

CFDX est organisé en couches indépendantes :

```text
Python API / CLI
        │
        ▼
Case / Runtime
        │
        ├──────── Configuration
        ├──────── Memory Planner
        ├──────── Execution Policy
        └──────── Checkpoint / I/O
        │
        ▼
Physics
        │
        ▼
FVM Operators
        │
        ▼
Linear Algebra
        │
        ▼
Execution Engine
        │
        ├──────── CPU backend
        ├──────── GPU backend
        └──────── GPU-OOC backend
        │
        ▼
Mesh / Geometry / Fields
        │
        ▼
Storage
        │
        └──────── HDF5
```

Aucune couche supérieure ne doit dépendre inutilement d'une couche d'implémentation inférieure.

---

# 3. Philosophie de développement

## 3.1 Module 0 avant la physique

Aucun solveur physique complet ne doit être développé avant validation du noyau de calcul.

Le premier objectif est :

> **Construire un noyau FVM générique capable de manipuler n'importe quel maillage polyédrique statique valide.**

Le Module 0 doit être :

* indépendant de la physique ;
* indépendant d'OpenFOAM ;
* indépendant du GPU ;
* utilisable en CPU ;
* compatible MPI ;
* capable de lire/écrire un cas complet HDF5.

---

# 4. Objectifs principaux

CFDX doit à terme fournir :

* maillage polyédrique arbitraire ;
* maillages structurés et non structurés ;
* import OpenFOAM ;
* import Gmsh ;
* import meshio ;
* HDF5 natif ;
* champs cellule/face/point/boundary ;
* opérateurs FVM ;
* algèbre linéaire sparse ;
* solveurs linéaires ;
* MPI ;
* CPU ;
* CUDA/GPU ;
* GPU out-of-core ;
* Navier-Stokes ;
* turbulence RANS ;
* thermique/CHT ;
* rayonnement ;
* VOF ;
* maillage dynamique ;
* FSI ;
* Python API ;
* checkpoints/restart ;
* post-processing ;
* validation quantitative.

---

# 5. Format de cas CFDX

## 5.1 Principe

Un cas CFDX doit être **autoportant**.

Un fichier :

```text
case.cfdx.h5
```

doit contenir toutes les informations nécessaires pour :

* reconstruire le maillage ;
* reconstruire le setup ;
* reconstruire les champs ;
* reconstruire les conditions limites ;
* reconstruire la physique ;
* reconstruire les paramètres numériques ;
* reconstruire le solveur ;
* reprendre un calcul ;
* identifier précisément l'état du cas.

Le but est d'éviter la dépendance à une arborescence de fichiers dispersée.

---

# 6. Structure HDF5

Structure cible :

```text
/
├── meta/
├── case/
├── mesh/
│   ├── points/
│   ├── faces/
│   ├── cells/
│   ├── boundary/
│   └── geometry/
├── fields/
├── boundary_conditions/
├── physics/
├── numerics/
├── solver/
├── runtime/
├── decomposition/
├── output/
└── results/
```

HDF5 fournit naturellement une hiérarchie de groupes et datasets permettant ce type d'organisation.

---

# 7. `/meta`

Contient les informations d'identification du fichier.

Exemple :

```text
/meta/
    format_version
    cfdx_version
    schema_version
    creation_date
    modification_date
    dimension
    precision
    endian
    mesh_hash
    topology_hash
    geometry_hash
    case_hash
```

## 7.1 Hash

Le hash doit permettre d'identifier précisément :

* la topologie ;
* le maillage ;
* le setup ;
* éventuellement l'état du calcul.

Le `topology_hash` doit être construit à partir des données fondamentales du maillage et non des données géométriques dérivées.

---

# 8. `/case`

Configuration générale :

```text
/case/
    name
    description
    author
    units
    dimension
    start_time
    end_time
    current_time
    time_step
```

---

# 9. `/mesh`

## 9.1 Données fondamentales

Les données fondamentales sont :

```text
points
face_vertices
face_vertices_offsets
face_owner
face_neighbour
cell_faces
cell_faces_offsets
boundary patches
```

Elles constituent la représentation topologique minimale.

---

# 10. Points

```text
/mesh/points
```

Format :

```text
float64[nPoints, 3]
```

Chaque point possède une position :

```text
x y z
```

---

# 11. Faces

## 11.1 Connectivité CSR

```text
/mesh/faces/vertices
/mesh/faces/vertices_offsets
```

Exemple :

```text
vertices:
[0,1,2,3, 4,5,6,7, ... ]

offsets:
[0,4,8,...]
```

Cela permet de représenter des faces de nombre de sommets arbitraire sans allocation par face.

---

# 12. Owner / neighbour

```text
/mesh/faces/owner
/mesh/faces/neighbour
```

Types :

```text
uint32 owner
int32 neighbour
```

Une face interne possède :

```text
owner >= 0
neighbour >= 0
```

Une face frontière :

```text
neighbour = -1
```

Cette représentation reprend le principe général de `polyMesh`, sans imposer son format de fichier.

---

# 13. Cell connectivity

```text
/mesh/cells/faces
/mesh/cells/faces_offsets
```

Format CSR.

Une cellule peut donc posséder un nombre quelconque de faces.

---

# 14. Boundary patches

```text
/mesh/boundary/patches/
```

Chaque patch contient notamment :

```text
name
type
face_ids
metadata
```

Exemple :

```text
wall
inlet
outlet
symmetry
periodic
interface
```

Le concept de patch reste indépendant des conditions physiques.

---

# 15. Géométrie

La géométrie est **dérivée de la topologie**.

```text
/mesh/geometry/
    cell_centres
    cell_volumes
    face_centres
    face_area_vectors
    face_areas
    delta_coeffs
    non_orthogonal_correction
    skewness
    non_orthogonality
```

Les données géométriques peuvent être stockées pour éviter leur recalcul.

Mais :

> **La topologie reste la source de vérité.**

Si la géométrie stockée est absente ou invalide, CFDX doit pouvoir la reconstruire.

---

# 16. Faces non planaires

Une face peut être non planaire.

CFDX ne doit pas modifier la topologie pour résoudre ce problème.

La topologie originale est conservée.

Une triangulation éventuelle est considérée comme une :

```text
geometric integration view
```

et non comme une nouvelle topologie.

Cela permet :

* calcul de surface ;
* centre de face ;
* intégration ;
* calcul des normales ;

sans perdre l'information originale.

---

# 17. Calcul géométrique

Le moteur géométrique doit calculer :

### Face

* centre ;
* aire ;
* vecteur surface `Sf` ;
* normale ;
* triangulation d'intégration si nécessaire.

### Cellule

* centre ;
* volume ;
* métriques géométriques.

### Qualité

* skewness ;
* non-orthogonality ;
* aspect ratio ;
* volume ratio ;
* distance owner/neighbour ;
* qualité des faces.

---

# 18. Invariants géométriques

Pour une cellule fermée :

```text
Σ Sf ≈ 0
```

Pour un volume valide :

```text
V > 0
```

Les critères doivent utiliser des tolérances absolues et relatives adaptées à l'échelle géométrique.

Il ne faut pas imposer arbitrairement :

```text
error < 1e-9
```

à tous les calculs.

---

# 19. Mesh Validator

Le validator doit être indépendant d'OpenFOAM.

Il vérifie :

### Topologie

* indices valides ;
* absence de références invalides ;
* absence de faces orphelines ;
* cohérence owner/neighbour ;
* cohérence cellules/faces ;
* doublons ;
* patches valides.

### Géométrie

* volumes positifs ;
* surfaces non nulles ;
* NaN ;
* Inf ;
* centres valides ;
* normales cohérentes.

### Qualité

* skewness ;
* non-orthogonality ;
* aspect ratio ;
* volume ratio ;
* cellules dégénérées.

### Conservation

```text
Σ Sf
```

et autres invariants géométriques.

---

# 20. Import / Export

Architecture :

```text
external format
      │
      ▼
parser
      │
      ▼
external representation
      │
      ▼
normalization
      │
      ▼
topology validation
      │
      ▼
boundary identification
      │
      ▼
geometry
      │
      ▼
quality
      │
      ▼
CFDX Mesh
      │
      ▼
case.cfdx.h5
```

Les importeurs ne doivent jamais contaminer le modèle interne CFDX.

---

# 21. OpenFOAM importer

Support :

```text
constant/polyMesh/
    points
    faces
    owner
    neighbour
    boundary
```

Le format `polyMesh` est basé précisément sur cette séparation points/faces/owner/neighbour/boundary.

Usage :

```python
mesh = cfdx.io.openfoam.read("constant/polyMesh")
```

---

# 22. Gmsh importer

Support prévu :

```text
.msh
```

avec conversion vers la représentation polyédrique CFDX.

---

# 23. meshio importer

`meshio` doit être supporté comme **adaptateur externe**, et non comme dépendance du cœur.

`meshio` fournit une interface Python pour lire de nombreux formats de maillage et expose notamment `points`, `cells`, etc.

Architecture :

```text
cfdx.io.meshio
```

Exemple :

```python
mesh = cfdx.io.meshio.read("mesh.vtu")
```

ou :

```python
mesh = cfdx.io.read(
    "mesh.vtu",
    format="meshio"
)
```

Le cœur CFDX ne doit jamais dépendre de l'objet `meshio.Mesh`.

---

# 23. Fluent importer

Support de deux formats Fluent :

**ASCII legacy** (`.cas` / `.dat`) — parseur S-expression pur, aucune dépendance externe :

```text
.cas   — maillage (S-expression Lisp style)
.dat   — données de résultats (format binaire/texte)
```

**HDF5** (`.cas.h5` + `.dat.h5`) — via les lecteurs *offline* de pyfluent (`CaseFile` + `FileSession`), aucune licence Fluent requise :

```text
.cas.h5 — maillage HDF5 (via CaseFile, h5py direct)
.dat.h5 — données de résultats (via FileSession, h5py direct)
```

L'adaptateur détecte automatiquement le format via `detect_source()` et route vers le parseur approprié :
- Si le fichier se termine par `.cas.h5` → utilise `CaseFile` de pyfluent (optionnel)
- Sinon → parseur ASCII legacy

Installation optionnelle :

```bash
pip install ansys-fluent-core  # fournit CaseFile + FileSession
```

---

# 23b. Adaptateurs solveurs externes (issue #425)

Adaptateurs additionnels implémentés :

| Solveur | Format | Parseur | Licence requise |
|---------|--------|---------|-----------------|
| STAR-CCM+ | `.sim` | Extraction de chaînes ASCII (scan binaire sûr) | Non |
| SU2 | `.su2` / `.cfg` | Parseur texte pur | Non |
| Code_Saturne | `.xml` / `.py` | Parseur XML + Python (regex) | Non |

STAR-CCM+ : le format `.sim` est binaire propriétaire. L'adaptateur effectue un *scan binaire sûr* pour extraire les métadonnées (noms de régions, matériaux, paramètres physiques) mais ne reconstruit pas la topologie du maillage. Rapport d'écart généré automatiquement.

Code_Saturne : l'adaptateur parse les fichiers de configuration XML et Python (`case.set()`, `case.set_*`) pour extraire BCs, matériaux, paramètres physiques. L'import du maillage nécessite une étape d'export vers un format neutre (MED/CGNS/VTK) — signalé dans le rapport d'écart.

---

# 24. Fields

Les champs sont génériques.

Concept :

```cpp
Field<T, Location>
```

Locations :

```text
CELL
FACE
POINT
BOUNDARY
```

Exemples :

```cpp
Field<double, Cell> p;
Field<Vec3, Cell> U;
Field<double, Face> phi;
```

---

# 25. Field metadata

Chaque champ possède :

```text
name
location
size
dimension
unit
precision
storage
```

Exemple :

```text
p:
    location = CELL
    dimension = pressure
    unit = Pa
```

Les unités doivent être résolues avant les kernels HPC.

Aucun objet unité lourd ne doit être présent dans les boucles numériques critiques.

---

# 26. BoundaryField

Le champ volumique et sa représentation frontière sont séparés.

```text
Field
BoundaryField
PatchField
```

Exemple :

```text
p.internal
p.boundary["inlet"]
p.boundary["wall"]
```

Les conditions limites sont donc attachées aux patches et non au maillage lui-même.

---

# 27. FVM operators

Module 0 expose :

```text
gradient
divergence
laplacian
interpolate
flux
surface_integrate
volume_integrate
```

---

# 28. Gradient

Première implémentation :

```text
Gauss linear
```

Puis :

```text
least_squares
weighted_least_squares
limited_gradient
```

Exemple :

```python
grad_p = cfdx.fvm.gradient(p)
```

---

# 29. Interpolation

Support :

```text
cell → face
```

avec plusieurs schémas :

```text
linear
upwind
limited
```

L'interpolation doit être générique et indépendante de la physique.

---

# 30. Divergence

Exemple :

```text
div(phi)
```

avec :

```text
Σ_f phi_f
```

et traitement cohérent des faces frontières.

---

# 31. Laplacien

Support :

```text
orthogonal
non-orthogonal corrected
non-orthogonal limited
uncorrected
```

La décomposition générale est :

```text
Laplacian =
orthogonal contribution
+
non-orthogonal correction
```

---

# 32. Rhie-Chow

Rhie-Chow **n'appartient pas au Module 0**.

Il appartient à la formulation pression-vitesse du Module 1.

Même principe pour :

* SIMPLE ;
* SIMPLEC ;
* PISO ;
* PIMPLE.

Module 0 fournit les opérateurs nécessaires, pas l'algorithme de couplage.

---

# 33. Algèbre linéaire

L'algèbre linéaire est indépendante de la physique.

Objets :

```text
Vector
SparseMatrix
LinearSystem
Preconditioner
LinearSolver
```

---

# 34. Sparse matrix

La représentation interne initiale peut être :

```text
CSR
```

mais l'API ne doit pas exposer directement cette représentation.

Exemple :

```cpp
SparseMatrix A;
Vector b;
Vector x;

LinearSystem system(A, b, x);
```

Cela permet de changer ultérieurement :

* CSR ;
* SELL-C-σ ;
* ELLPACK ;
* BSR ;
* formats GPU spécifiques ;

sans modifier la physique.

---

# 35. Solveurs linéaires

Premiers solveurs :

```text
CG
BiCGStab
GMRES
```

Puis éventuellement :

```text
FGMRES
MINRES
```

---

# 36. Préconditionneurs

Première génération :

```text
Jacobi
Gauss-Seidel
ILU
```

Puis :

```text
AMG
Block-AMG
domain decomposition
```

Le préconditionneur est séparé du solveur.

---

# 37. Module 0 — contenu définitif

Le Module 0 contient :

```text
M0.1  Mesh topology
M0.2  Geometry
M0.3  Mesh quality
M0.4  Fields
M0.5  Boundary fields
M0.6  Interpolation
M0.7  FVM operators
M0.8  Linear algebra
M0.9  HDF5
M0.10 Import/export
M0.11 MPI infrastructure
M0.12 Execution abstraction
```

Il ne contient pas :

```text
SIMPLE
PISO
PIMPLE
Rhie-Chow
turbulence
VOF
thermal models
radiation
FSI
```

---

# 38. CPU memory model

La représentation mémoire doit être optimisée pour :

* cache ;
* SIMD ;
* vectorisation ;
* faible overhead ;
* accès séquentiels.

La représentation logique ne doit pas imposer un layout physique unique.

Le backend peut choisir :

```text
AoS
SoA
AoSoA
```

selon le kernel.

---

# 39. Device abstraction

Un `Field` ne doit pas contenir directement :

```cpp
device_ptr
```

ou une logique CUDA spécifique.

Architecture :

```text
Field
  │
  ▼
StorageHandle
  │
  ├── HostStorage
  ├── DeviceStorage
  └── WorkingSetStorage
```

---

# 40. Host storage

Le stockage CPU est la représentation persistante principale.

Mais il ne faut pas considérer que le CPU est toujours synchronisé avec le GPU.

Pendant un calcul GPU :

```text
Host : iteration n
GPU  : iteration n+1
```

est un état valide.

La synchronisation est explicite.

---

# 41. Storage states

Le système doit pouvoir représenter :

```text
HOST_ONLY
DEVICE_ONLY
SYNCHRONIZED
HOST_DIRTY
DEVICE_DIRTY
```

ou un mécanisme équivalent plus robuste.

Un simple :

```cpp
bool is_device_synced;
```

est insuffisant.

---

# 42. Execution Policy

CFDX définit :

```cpp
enum class ExecutionPolicy
{
    CPU,
    GPU,
    GPU_OUT_OF_CORE,
    AUTO
};
```

`AUTO` choisit la stratégie **avant le lancement du calcul**.

Il n'y a pas de changement silencieux de backend pendant une simulation.

Cela garantit :

* reproductibilité ;
* benchmarking ;
* prédictibilité ;
* diagnostic des performances.

---

# 43. Runtime

Le Runtime devient responsable de :

```text
ExecutionPolicy
MemoryManager
MemoryPlanner
StorageManager
WorkingSetManager
TileManager
HaloManager
AsyncTransferManager
CPU backend
GPU backend
```

---

# 44. Memory Planner

Le Memory Planner doit estimer :

```text
mesh
geometry
fields
linear algebra
preconditioner
temporary buffers
MPI buffers
CUDA workspace
halo buffers
output buffers
```

Il ne doit pas simplement regarder :

```text
cudaMemGetInfo()
```

---

# 45. Memory budget

Objet :

```text
MemoryBudget
```

Contient :

```text
VRAM total
runtime reserve
CFDX available budget
safety margin
estimated usage
peak usage
```

Exemple conceptuel :

```text
VRAM = 24 GB

runtime reserve = 2 GB
safety margin = 2 GB

CFDX budget = 20 GB
```

---

# 46. Full GPU

Si :

```text
estimated_memory <= CFDX_GPU_budget
```

alors :

```text
GPU mode
```

Le maillage, les champs et les structures numériques nécessaires restent en VRAM.

Flux :

```text
H2D initial
      ↓
GPU iterations
      ↓
GPU iterations
      ↓
GPU iterations
      ↓
D2H checkpoint/output
```

Il ne faut surtout pas transférer les champs CPU↔GPU à chaque itération.

---

# 47. GPU Out-of-Core

Si le problème ne tient pas dans la VRAM :

```text
GPU_OUT_OF_CORE
```

Ce mode n'est **pas** une résolution indépendante de chunks.

Il repose sur :

```text
domain decomposition
+
tiles
+
halo exchange
+
working sets
+
asynchronous transfers
```

---

# 48. OOC architecture

```text
Global mesh
    │
    ▼
Domain decomposition
    │
    ├── Tile 0 + halo
    ├── Tile 1 + halo
    ├── Tile 2 + halo
    └── ...
            │
            ▼
      GPU working set
            │
            ▼
       compute kernel
            │
            ▼
       halo update
```

Les interactions entre tuiles doivent être conservées.

---

# 49. Halo management

Chaque tile possède :

```text
owned cells
halo cells
neighbor relations
```

Le `HaloManager` gère :

```text
pack
transfer
unpack
synchronization
```

Cela doit fonctionner :

* CPU↔CPU ;
* MPI↔MPI ;
* CPU↔GPU ;
* GPU↔GPU ;
* GPU↔MPI.

---

# 50. Async transfers

Le système GPU doit utiliser lorsque possible :

```text
CUDA streams
asynchronous H2D
asynchronous D2H
double buffering
```

Architecture :

```text
Transfer tile A
       │
       ├── compute tile B
       │
       └── transfer tile C
```

Le but est de masquer une partie du coût PCIe/NVLink.

---

# 51. Pinned memory

La mémoire pinned n'est pas utilisée comme mémoire générale.

CFDX utilise :

```text
Host RAM
    +
PinnedBufferPool
```

Le `PinnedBufferPool` est borné.

Il sert principalement aux transferts fréquents :

```text
H2D
D2H
MPI staging
```

---

# 52. Performance policy

Trois modes ont des objectifs différents :

| Mode    | Objectif                        |
| ------- | ------------------------------- |
| CPU     | capacité + simplicité           |
| GPU     | performance                     |
| GPU-OOC | capacité GPU au-delà de la VRAM |

Le GPU-OOC ne doit pas être présenté comme nécessairement plus rapide que le CPU.

Pour certaines applications :

```text
CPU RAM bandwidth
>
PCIe transfer + GPU compute
```

peut être vrai.

---

# 53. Memory-aware execution

Le runtime choisit :

```text
CPU
GPU
GPU-OOC
```

selon :

* taille du problème ;
* mémoire disponible ;
* backend disponible ;
* nombre de GPU ;
* MPI ;
* coût estimé des transferts ;
* configuration utilisateur.

---

# 54. MPI

MPI appartient à l'infrastructure d'exécution.

Le modèle initial est :

```text
MPI rank
    ↓
local mesh
    ↓
local fields
    ↓
local operators
```

---

# 55. Domain decomposition

La décomposition doit séparer :

```text
global topology
```

et :

```text
local topology
```

Chaque rank possède :

```text
owned cells
ghost cells
local faces
interface faces
```

---

# 56. MPI + GPU

Architecture cible :

```text
MPI Rank
    │
    ├── CPU memory
    ├── GPU memory
    │
    └── local execution
```

Les communications inter-ranks doivent pouvoir être réalisées :

```text
CPU staging
```

ou ultérieurement :

```text
GPU-aware MPI
```

---

# 57. Checkpoint / Restart

Un checkpoint doit permettre de reconstruire :

* maillage ;
* champs ;
* état physique ;
* temps ;
* timestep ;
* solveur ;
* configuration numérique ;
* décomposition si nécessaire ;
* état nécessaire au restart.

---

# 58. Restart indépendant du nombre de MPI ranks

Objectif important :

```text
restart generated with 64 MPI ranks
```

doit pouvoir être repris avec :

```text
32 MPI ranks
128 MPI ranks
```

si les données globales nécessaires sont disponibles.

La décomposition doit donc être considérée comme une donnée d'exécution et non comme une partie immuable du cas physique.

---

# 59. Output

Configuration :

```text
/output/
```

Exemple :

```yaml
fields:
    enabled: true
    frequency: 200
    variables:
        - U
        - p
        - alpha

checkpoints:
    frequency: 1000
    keep_last: 3

residuals:
    frequency: 1

statistics:
    frequency: 10
```

---

# 60. Results

Les résultats peuvent être stockés dans :

```text
/results/
```

avec :

```text
time
iteration
fields
residuals
statistics
forces
fluxes
probes
```

Le format doit permettre une lecture indépendante du solveur.

---

# 61. `data.h5`

À côté de :

```text
case.cfdx.h5
```

CFDX peut utiliser :

```text
data.h5
```

pour les données scientifiques générales.

Différence :

```text
case.cfdx.h5
```

= cas CFD autoportant.

```text
data.h5
```

= conteneur scientifique générique.

---

# 62. Python API

Python est la couche d'orchestration.

Exemple :

```python
case = cfdx.Case("cavity.cfdx.h5")

case.mesh.check()

case.fields["U"].set(...)
case.fields["p"].set(...)

case.physics.enable("incompressible")

case.run()
```

---

# 63. Python ne doit pas être le moteur HPC

Python gère :

* configuration ;
* construction ;
* inspection ;
* automation ;
* post-processing ;
* workflow ;
* scripting ;
* orchestration.

C++ gère :

* mesh kernels ;
* geometry ;
* FVM kernels ;
* sparse algebra ;
* MPI ;
* CUDA ;
* memory management ;
* solveurs.

---

# 64. C++ API

Le C++ constitue le cœur numérique.

Architecture cible :

```text
cfdx/
├── core/
│   ├── mesh/
│   ├── geometry/
│   ├── field/
│   ├── boundary/
│   ├── operators/
│   ├── linalg/
│   ├── parallel/
│   └── units/
│
├── io/
│   ├── hdf5/
│   ├── openfoam/
│   ├── gmsh/
│   └── meshio/
│
├── numerics/
│   ├── interpolation/
│   ├── gradients/
│   ├── laplacian/
│   └── convection/
│
├── runtime/
│   ├── execution/
│   ├── memory/
│   ├── cpu/
│   ├── gpu/
│   └── ooc/
│
├── solvers/
│
├── physics/
│   ├── incompressible/
│   ├── turbulence/
│   ├── thermal/
│   ├── radiation/
│   ├── vof/
│   ├── dynamic_mesh/
│   └── fsi/
│
├── python/
└── tests/
```

---

# 65. Physics modules

## M1 — Incompressible laminar flow

Contient :

* Navier-Stokes ;
* continuity ;
* pressure-velocity coupling ;
* SIMPLE ;
* SIMPLEC ;
* PISO ;
* PIMPLE ;
* Rhie-Chow.

---

# 66. M2 — Turbulence

Support initial :

```text
RANS
```

Modèles :

```text
k-epsilon
k-omega
SST
```

Puis :

```text
LES
DES
```

---

# 67. M3 — Thermal / CHT

Support :

```text
energy
conduction
convection
conjugate heat transfer
solid regions
fluid regions
```

Couplage multi-régions.

---

# 68. M4 — Radiation

Support futur :

```text
surface radiation
participating media
view factors
DOM
P1
```

Le modèle exact sera défini après validation du noyau thermique.

---

# 69. M5 — VOF

Support :

```text
volume fraction
interface reconstruction
surface tension
contact angle
compressive schemes
```

---

# 70. M6 — Dynamic mesh

Support :

```text
mesh motion
deformation
remeshing
topology changes
```

La géométrie devra alors devenir recalculable à chaque modification topologique/géométrique.

---

# 71. M7 — FSI

Support :

```text
fluid
structure
interface coupling
```

avec éventuellement :

```text
partitioned coupling
monolithic coupling
```

---

# 72. Configuration

Le setup utilisateur doit être représentable dans HDF5 mais également manipulable en Python.

Exemple :

```python
case.numerics.gradient.scheme = "gauss_linear"

case.numerics.laplacian.scheme = "corrected"

case.solver.pressure.algorithm = "GAMG"

case.runtime.execution = "GPU"
```

Le format interne reste indépendant de la syntaxe Python.

---

# 73. Reproductibilité

Un cas doit pouvoir enregistrer :

```text
CFDX version
schema version
compiler
compiler version
CUDA version
MPI implementation
GPU
CPU
precision
numerical schemes
solver settings
mesh hash
case hash
```

---

# 74. Déterminisme

CFDX doit distinguer :

```text
deterministic mode
performance mode
```

Le mode déterministe peut imposer :

* ordre de réduction ;
* ordre de somme ;
* synchronisations supplémentaires.

Le mode performance permet :

* reductions parallèles ;
* atomics ;
* ordre non déterministe contrôlé.

---

# 75. Precision

Support cible :

```text
float32
float64
```

Le modèle doit permettre à terme :

```text
mixed precision
```

par exemple :

```text
FP32 compute
FP64 residual accumulation
```

sans modifier la physique.

---

# 76. Validation analytique

Le Module 0 doit être validé sur :

```text
cube
tetrahedron
pyramid
prism
arbitrary polyhedron
```

---

# 77. Tests gradient

Champs :

```text
constant
linear
quadratic
```

Vérification :

```text
∇constant = 0
```

et pour un champ linéaire :

```text
∇φ = exact
```

dans les limites de la discrétisation.

---

# 78. Tests laplacien

Tester :

```text
constant
linear
quadratic
```

avec :

```text
Δconstant = 0
Δlinear = 0
```

et cas quadratique analytique.

---

# 79. Tests conservation

Tester :

```text
Σ flux = 0
```

pour des volumes fermés.

Tester également la conservation globale sur :

* maillage orthogonal ;
* maillage non orthogonal ;
* maillage skewed.

---

# 80. Comparaison OpenFOAM

OpenFOAM est une **référence de comparaison**, pas la définition mathématique de la vérité.

Comparaison :

```text
same mesh
same BC
same discretization
same physical parameters
```

puis :

```text
geometry
operators
linear system
solution
conservation
performance
memory
```

---

# 81. Niveau de comparaison

La validation doit être hiérarchique :

### Niveau 1

Géométrie.

### Niveau 2

Opérateurs.

### Niveau 3

Matrice.

### Niveau 4

Solution analytique.

### Niveau 5

Comparaison OpenFOAM.

### Niveau 6

Performance.

Cela permet de savoir où apparaît une divergence.

---

# 82. Performance benchmark

Mesures :

```text
wall time
CPU time
memory
VRAM
bandwidth
FLOP/s
MPI communication
GPU utilization
PCIe traffic
```

Pour GPU-OOC :

```text
H2D bandwidth
D2H bandwidth
halo traffic
compute/transfer overlap
```

---

# 83. Objectifs de performance

Les objectifs du type :

```text
-20 % memory
-15 % wall time
```

ne doivent pas être des critères arbitraires avant mesure.

La procédure est :

```text
baseline
→ benchmark
→ profiling
→ optimization
→ benchmark
```

Les gains sont mesurés par rapport à une baseline documentée.

---

# 84. Memory benchmark

Pour chaque cas :

```text
mesh memory
field memory
matrix memory
preconditioner memory
temporary memory
peak memory
```

Comparer :

```text
CFDX CPU
CFDX GPU
CFDX GPU-OOC
OpenFOAM
```

---

# 85. Profiling

CFDX doit prévoir des hooks de profiling :

```text
mesh
geometry
assembly
gradient
divergence
laplacian
matrix-vector
preconditioner
solver
MPI
H2D
D2H
halo
I/O
```

---

# 86. Logging

Niveaux :

```text
ERROR
WARNING
INFO
DEBUG
TRACE
```

Le runtime doit pouvoir produire un résumé :

```text
Execution policy : GPU
GPU memory        : 18.4 / 24 GB
Peak memory       : 19.1 GB
Iterations        : 1200
Solver time       : ...
H2D time          : ...
D2H time          : ...
MPI time          : ...
```

---

# 87. CLI

Objectif :

```bash
cfdx check case.cfdx.h5
cfdx info case.cfdx.h5
cfdx run case.cfdx.h5
cfdx convert mesh.msh case.cfdx.h5
cfdx inspect case.cfdx.h5
cfdx benchmark case.cfdx.h5
```

---

# 88. Inspection HDF5

Un cas doit être inspectable sans CFDX complet.

Exemple :

```bash
h5ls case.cfdx.h5
h5dump case.cfdx.h5
```

HDF5 est conçu précisément pour fournir des datasets, groupes et métadonnées auto-descriptifs.

---

# 89. Robustesse du format

Le format HDF5 doit être versionné :

```text
schema_version
```

Exemple :

```text
1.0
1.1
2.0
```

CFDX doit pouvoir détecter :

```text
unsupported schema
missing mandatory dataset
invalid dtype
invalid dimension
```

---

# 90. Compatibilité

Le format doit privilégier :

* stabilité ;
* auto-description ;
* migration explicite ;
* backward compatibility raisonnable.

Les changements incompatibles doivent incrémenter la version majeure du schema.

---

# 91. Séparation case / runtime

Les paramètres physiques et numériques doivent être séparés des paramètres d'exécution.

Exemple :

```text
case:
    physics
    numerics
    mesh
```

versus :

```text
runtime:
    CPU/GPU
    MPI ranks
    tile size
    memory policy
    streams
```

Ainsi le même cas peut être exécuté :

```text
CPU
GPU
GPU-OOC
32 MPI
128 MPI
```

sans modifier sa définition physique.

---

# 92. Execution graph

Architecture finale :

```text
                 CFDX CASE
                     │
                     ▼
              Configuration
                     │
                     ▼
              Memory Planner
                     │
        ┌────────────┼────────────┐
        ▼            ▼            ▼
       CPU          GPU        GPU-OOC
        │            │            │
        │            │       Tile Manager
        │            │            │
        │            │       Halo Manager
        │            │            │
        └────────────┴────────────┘
                     │
                     ▼
              Execution Engine
                     │
                     ▼
                 FVM kernels
                     │
                     ▼
               Linear Algebra
                     │
                     ▼
                  Physics
```

La relation physique réelle est :

```text
Physics
   ↓
FVM operators
   ↓
Execution Engine
   ↓
backend
```

et non :

```text
Physics
   ↓
CUDA
```

---

# 93. Principe essentiel du GPU

Le GPU est un **backend d'exécution**, pas une architecture imposée au modèle scientifique.

Le même objet :

```text
Field<double, Cell>
```

doit pouvoir être exécuté :

```text
CPU
GPU
GPU-OOC
```

---

# 94. Principe essentiel du OOC

Le mode out-of-core est un mode de **capacité**, pas simplement un mode de performance.

Il permet :

```text
Problem size > VRAM
```

mais peut être plus lent qu'une exécution CPU.

Il doit donc être explicitement visible dans les métriques et le reporting.

---

# 95. Roadmap

## Phase 0 — Data model

Construire :

```text
Case
Mesh
Field
BoundaryField
HDF5
Hash
```

Objectif :

```text
case.cfdx.h5
```

autoportant.

---

## Phase 1 — Import

Implémenter :

```text
OpenFOAM
Gmsh
meshio
```

avec pipeline commun.

---

## Phase 2 — Geometry

Implémenter :

```text
face geometry
cell geometry
quality
validation
```

---

## Phase 3 — Operators

Implémenter :

```text
gradient
interpolation
divergence
laplacian
flux
```

---

## Phase 4 — Scalar solver

Construire :

```text
Poisson
Laplace
```

sur maillage arbitraire.

C'est le premier véritable solveur CFDX.

---

## Phase 5 — Linear algebra

Implémenter :

```text
CSR
CG
BiCGStab
GMRES
Jacobi
ILU
AMG
```

progressivement.

---

## Phase 6 — MPI

Implémenter :

```text
domain decomposition
ghost cells
halo exchange
parallel HDF5
```

---

## Phase 7 — GPU

Implémenter :

```text
CUDA backend
device storage
memory planner
GPU kernels
```

---

## Phase 8 — GPU OOC

Implémenter :

```text
tile manager
halo manager
working sets
pinned buffer pool
async transfers
double buffering
```

---

## Phase 9 — Incompressible

Module 1 :

```text
Navier-Stokes
continuity
SIMPLE
PISO
PIMPLE
Rhie-Chow
```

---

## Phase 10 — Physics

Puis :

```text
M2 turbulence
M3 thermal / CHT
M4 radiation
M5 VOF
M6 dynamic mesh
M7 FSI
```

---

# 96. Critères de sortie du Module 0

Le Module 0 est considéré comme terminé uniquement lorsque CFDX peut :

### Data

* [ ] créer un `.cfdx.h5` ;
* [ ] recharger le fichier ;
* [ ] reconstruire le mesh ;
* [ ] reconstruire les fields ;
* [ ] calculer les hashes ;
* [ ] vérifier l'intégrité.

### Mesh

* [ ] cube ;
* [ ] tetrahedron ;
* [ ] prism ;
* [ ] pyramid ;
* [ ] polyhedron arbitraire ;
* [ ] mesh OpenFOAM ;
* [ ] mesh Gmsh ;
* [ ] mesh meshio.

### Geometry

* [ ] centres ;
* [ ] volumes ;
* [ ] surfaces ;
* [ ] normales ;
* [ ] skewness ;
* [ ] non-orthogonality.

### Operators

* [ ] gradient ;
* [ ] interpolation ;
* [ ] divergence ;
* [ ] laplacian ;
* [ ] conservation.

### Linear algebra

* [ ] sparse matrix ;
* [ ] vector ;
* [ ] CG ;
* [ ] BiCGStab ;
* [ ] GMRES ;
* [ ] au moins un préconditionneur robuste.

### Solver

* [ ] Poisson analytique ;
* [ ] Laplace analytique.

### Parallel

* [ ] MPI ;
* [ ] decomposition ;
* [ ] ghost cells ;
* [ ] halo exchange.

### GPU

* [ ] CPU backend ;
* [ ] GPU backend ;
* [ ] memory planner ;
* [ ] full GPU ;
* [ ] GPU-OOC prototype.

### Validation

* [ ] tests analytiques ;
* [ ] comparaison OpenFOAM ;
* [ ] benchmark ;
* [ ] profiling ;
* [ ] memory benchmark.

---

# 97. Architecture finale visée

La vision complète est :

```text
                         PYTHON
                           │
                 ┌─────────┴─────────┐
                 │       CASE        │
                 │ configuration     │
                 │ workflow          │
                 └─────────┬─────────┘
                           │
                           ▼
                    ┌─────────────┐
                    │   RUNTIME   │
                    │             │
                    │ execution   │
                    │ memory      │
                    │ MPI         │
                    │ I/O         │
                    └──────┬──────┘
                           │
              ┌────────────┼────────────┐
              ▼            ▼            ▼
             CPU          GPU        GPU-OOC
              │            │            │
              └────────────┼────────────┘
                           ▼
                    ┌─────────────┐
                    │ FVM ENGINE  │
                    └──────┬──────┘
                           │
             ┌─────────────┼─────────────┐
             ▼             ▼             ▼
          gradient      divergence    laplacian
             │             │             │
             └─────────────┼─────────────┘
                           ▼
                  ┌─────────────────┐
                  │ LINEAR ALGEBRA  │
                  └────────┬────────┘
                           │
                           ▼
                    ┌─────────────┐
                    │   PHYSICS   │
                    ├─────────────┤
                    │ NavierStokes│
                    │ Turbulence  │
                    │ Thermal     │
                    │ Radiation   │
                    │ VOF         │
                    │ DynamicMesh │
                    │ FSI         │
                    └─────────────┘

                 STORAGE
                    │
                    ▼
              case.cfdx.h5
                    │
             ┌──────┴──────┐
             ▼             ▼
          /mesh          /fields
             │             │
             ├─────────────┤
             │             │
             ▼             ▼
         /physics       /solver
             │             │
             └──────┬──────┘
                    ▼
                /results
```

---

# 98. Décisions d'architecture figées

Les décisions suivantes sont considérées comme fondamentales pour CFDX v0.7 :

1. **HDF5 est le format natif du cas.**
2. **`case.cfdx.h5` est autoportant.**
3. **La topologie est la source de vérité du maillage.**
4. **La géométrie est dérivée de la topologie.**
5. **OpenFOAM est une référence, pas la définition de CFDX.**
6. **meshio est un adaptateur, pas une dépendance du cœur.**
7. **Les Fields sont indépendants du backend matériel.**
8. **L'algèbre linéaire est indépendante de la physique.**
9. **Les opérateurs FVM sont indépendants du backend.**
10. **SIMPLE/PISO/PIMPLE/Rhie-Chow appartiennent au Module 1.**
11. **CPU/GPU/GPU-OOC sont des Execution Policies.**
12. **Aucun fallback GPU→CPU silencieux pendant un calcul.**
13. **Le choix CPU/GPU/GPU-OOC est fait par le Runtime avant exécution.**
14. **Le GPU-OOC utilise domain decomposition + tiles + halos.**
15. **Les transferts CPU↔GPU ne doivent pas être présents dans la boucle GPU normale.**
16. **Pinned memory est limitée à des buffers de staging.**
17. **La physique ne connaît jamais CUDA/MPI/mémoire GPU directement.**
18. **Python orchestre ; C++ calcule.**
19. **Le parallélisme est une propriété du Runtime, pas de la physique.**
20. **Les benchmarks sont mesurés et non définis arbitrairement à l'avance.**

---

# 99. Principe directeur final

CFDX ne doit pas être conçu comme :

> « OpenFOAM réécrit en C++ moderne avec un wrapper Python et un backend GPU ».

Il doit être conçu comme :

> **un moteur FVM généraliste où le modèle scientifique, les opérateurs numériques, les données et l'exécution sont explicitement découplés.**

La conséquence architecturale essentielle est :

```text
                 SCIENTIFIC MODEL
                        │
                        ▼
                   FVM OPERATORS
                        │
                        ▼
                  LINEAR ALGEBRA
                        │
                        ▼
                EXECUTION ENGINE
                        │
          ┌─────────────┼─────────────┐
          ▼             ▼             ▼
         CPU           GPU          GPU-OOC
```

et non :

```text
Physics → CPU
Physics → CUDA
Physics → MPI
```

C'est cette séparation qui doit permettre à CFDX d'évoluer progressivement d'un noyau FVM simple vers un solveur HPC complet sans devoir réécrire les modules physiques à chaque changement d'architecture matérielle.
