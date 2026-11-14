# CFDX V23 — preconditionnement local 2x2 SST et sauvegarde du JFNK V22

Date : 2026-11-10. **Statut experimental, non qualifie physiquement.**

## Constat et methode

Le BFS 712 cellules V22 montre un pas Newton-Jacobi reel accepte a la dixieme
iteration interne, mais pas de convergence SST ou U-p-k-omega globale.
Un simple equilibrage des lignes (V22) n'approxime pas les derivees croisees
entre production, destruction et diffusion croisee SST. V23 implemente un
preconditionneur a droite par blocs k-omega, sans changer l'operateur FV ni
les criteres de convergence.

Avec z = [log(k/k_ref), log(omega/omega_ref)], le Newton vrai reste
J(z) delta_z = -F(z), ou F est le vecteur des residus volumes finis
reassembles. Une matrice P locale par cellule approxime J et GMRES resout
(J P^{-1})y=-F, delta_z=P^{-1}y. Le controle Krylov est toujours effectue
sur le produit directionnel vrai J*delta_z+F, et l'acceptation est effectuee
sur le residu non lineaire reassemble apres changement de paroi.

Pour chaque cellule, A_ii provient des diagonales FV existantes, incluant
la destruction implicitement gelee. Avec S_k=sk+spk*k, S_w=sw+spw*omega,
les approximations locales des derivees physiques sont :

- J_kk = A_k,ii - V*(d(sk)/dk + k*d(spk)/dk)
- J_kw = -V*(d(sk)/domega + k*d(spk)/domega)
- J_wk = -V*(d(sw)/dk + omega*d(spw)/dk)
- J_ww = A_w,ii - V*(d(sw)/domega + omega*d(spw)/domega)

Les derivees de sk, sw, spk, spw sont evaluees localement avec F1/F2
recalcules sur la cellule et la diffusion croisee de base gelee. Les blocs
sont transformes en variables logarithmiques et avec les echelles de lignes
figees du solveur Newton. Cette construction est O(N) pour N cellules. Elle
n'inclut PAS les derivees spatiales des coefficients de diffusion, les
gradients de diffusion croisee, la sensibilite au rafraichissement de G ni
les derivees de la cible omega de paroi. Toutes ces dependances restent
presentes dans le produit Jacobien-vecteur calcule par reassemblage FV.
La ligne omega proprietaire contrainte n'est PAS traitee comme une source
volumique. En cas de bloc singulier ou non fini, seul le preconditionneur
revient a une action identite; la physique et les criteres demeurent inchanges.

La premiere tentative utilisant seulement les blocs a echoue : residu
lineaire directionnel GMRES ~0.824 a la dixieme iteration SST, donc direction
rejete; merite du pas Picard ~0.956, moins bon que V22. V23 applique alors
une **comparaison de propositions**: l'ancienne direction Jacobi est toujours
calculee et conservee; une direction bloc supplementaire est essayee a la
premiere iteration interne puis toutes les dix iterations. Picard, Newton-
Jacobi, Newton-bloc sont tries par merit physique apres reassemblage, sans
changer les tolerances ni le nombre nominal d'iterations SST.

## Preuves sur le BFS reduit

La meme image HDF5 de 712 cellules, les memes options physiques et le meme
budget diagnostic de dix iterations SST sont utilises. L'execution longue
complete n'a pas ete validee : les appels ont ete interrompus pendant le
rejeu des retries exterieurs. Le traceur prend le premier essai exterieur.

| A l'iteration SST 10 | V22 Jacobi | V23 blocs seuls | V23 sauvegarde |
|---|---:|---:|---:|
| Merite local reassemble | 0.8899074514 | 0.9557896839 | 0.8899074514 |
| R_k brut | 3.960597672e-6 | 0.001211639826 | 3.960597672e-6 |
| R_omega brut | 1.769146644 | 0.6352146560 | 1.769146644 |
| GMRES vrai, direction bloc | sans objet | 0.82377488678 | 0.74693614522 |
| GMRES vrai, direction Jacobi | 0.24793702144 | sans objet | 0.24793702144 |
| Direction selectionnee | Newton Jacobi | Picard | Newton Jacobi |
| Pas SST physiquement converge | NON | NON | NON |

Ce tableau interdit de decrire le bloc local comme une amelioration de
convergence. Il donne une preuve d'absence de regression **a cet instant**
seulement pour la version sauvegardee, pas une garantie globale.

## Limitations et campagnes requises

- La tolerance physique SST demeure 1e-8; le forcing Krylov 0.25 est une
  regle interne, non un certificat de convergence.
- Le portage exact des poids OpenFOAM Foundation `omegaWallFunction::updateCoeffs`
  et `manipulateMatrix` n'est pas fait; les differences de contraintes sont
  une source plausible de mauvaise qualite des blocs pres des parois.
- Ni premier ni deuxieme etat exterieur V23 converges, ni preuve de 3 etats
  acceptes; BFS 20540, six Schur, comparaison NASA non executes.
- Le preconditionneur est CPU serial local. Aucune preuve de performance GPU,
  MPI ni multi-GPU. Le calcul O(N) est adaptable a CPU parallel/GPU si les
  sorties sources et les echelles sont organisees en tableaux de structures
  explicites, sans promesse de performance.
- Le parcours de qualification doit confronter le Jacobien local approximate
  aux blocs vrais de differences finies sur des cellules temoins 357,363,426,
  puis evaluer une factorisation bloc des matrices de convection-diffusion.

## Sources publiques

- PETSc SNES : https://petsc.org/release/manual/snes/
- OpenFOAM Foundation v11 : https://cpp.openfoam.org/v11/omegaWallFunctionFvPatchScalarField_8C_source.html

Le mode active par `numerics.raw_settings["sst_block_preconditioner_v23_experimental"]`
necessite egalement `sst_jfnk_v22_experimental=true` et reste desactive par
defaut. Aucun push GitHub ni changement de seuil physique.
