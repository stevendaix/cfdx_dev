# %% [markdown]
# # Thermohydraulique turbulente — base théorique et V&V
#
# Ce chapitre fixe les équations, conventions et limites d'interprétation communes
# aux campagnes TUR-CHANNEL, SEP-BFS-DRIVER, TH-MIXED-CONVECTION-PIPE et TH-TJUNCTION.
# TH-07 est une autre reference : DNS de canal plan vertical Wagner-Wetzel.
# Ne pas assimiler l'experience 1998 de conduite a cette DNS.
# Il constitue une base scientifique, pas une preuve d'exécution ni de qualification.
#
# ## 1. Équations moyennées de quantité de mouvement
#
# Pour un fluide incompressible newtonien à propriétés constantes, la décomposition
# de Reynolds est :
#
# \\[
# u_i=U_i+u_i', \\qquad p=P+p', \\qquad
# \\overline{u_i'}=0.
# \\]
#
# La continuité et les équations RANS s'écrivent :
#
# \\[
# \\frac{\\partial U_i}{\\partial x_i}=0,
# \\]
# \\[
# \\rho U_j\\frac{\\partial U_i}{\\partial x_j}
# =-\\frac{\\partial P}{\\partial x_i}
# +\\frac{\\partial}{\\partial x_j}
# \\left(2\\mu S_{ij}-\\rho\\overline{u_i'u_j'}\\right)
# +f_i,
# \\qquad
# S_{ij}=\\frac12\\left(
# \\frac{\\partial U_i}{\\partial x_j}
# +\\frac{\\partial U_j}{\\partial x_i}
# \\right).
# \\]
#
# Le tenseur de contraintes de Reynolds
# \\(R_{ij}=\\overline{u_i'u_j'}\\) est une inconnue de fermeture.
# L'hypothèse de Boussinesq des modèles à viscosité turbulente pose :
#
# \\[
# -\\rho R_{ij}=2\\mu_t S_{ij}
# -\\frac23\\rho k\\delta_{ij},
# \\qquad k=\\frac12R_{ii}.
# \\]
#
# Cette hypothèse doit être explicitement prise en compte lors de la comparaison
# avec les six composantes DNS : un modèle linéaire à viscosité turbulente ne
# représente pas librement l'anisotropie du tenseur. Les modèles RSM doivent être
# évalués séparément.
#
# ## 2. Équation d'énergie moyenne et flux thermique turbulent
#
# Pour \\(c_p\\) constant et sans source volumique :
#
# \\[
# \\rho c_p U_j\\frac{\\partial\\overline T}{\\partial x_j}
# =\\frac{\\partial}{\\partial x_j}
# \\left(k\\frac{\\partial\\overline T}{\\partial x_j}
# -\\rho c_p\\overline{u_j'T'}\\right)+\\dot q_v.
# \\]
#
# Le flux turbulent \\(\\overline{u_j'T'}\\) est une grandeur distincte du
# flux conductif. Une fermeture gradient-diffusion classique est :
#
# \\[
# -\\rho c_p\\overline{u_j'T'}
# =\\frac{\\mu_t c_p}{Pr_t}
# \\frac{\\partial\\overline T}{\\partial x_j},
# \\qquad
# \\alpha_t=\\frac{\\nu_t}{Pr_t}.
# \\]
#
# Le nombre de Prandtl turbulent \\(Pr_t\\) n'est pas universel. Sa valeur,
# sa corrélation et son éventuelle dépendance à \\(Re_t\\), \\(Pr\\) ou
# \\(y^+\\) doivent être gelées dans chaque campagne. Une validation de la
# vitesse ne valide pas automatiquement le flux thermique turbulent.
#
# ## 3. Flottabilité et convection mixte
#
# L'approximation de Boussinesq suppose des variations de densité faibles,
# conservées uniquement dans le terme de gravité :
#
# \\[
# \\rho(T)\\simeq\\rho_0[1-\\beta(T-T_0)],
# \\qquad
# \\boldsymbol f_b=-\\rho_0\\beta(T-T_0)\\boldsymbol g.
# \\]
#
# Les propriétés \\(\\rho_0,\\mu,k,c_p,\\beta\\) sont évaluées à une
# température de référence. Si les variations de propriétés sont significatives,
# il faut résoudre avec les lois \\(\\rho(T),\\mu(T),k(T),c_p(T)\\) du
# benchmark plutôt que conserver Boussinesq par commodité.
#
# Les nombres sans dimension usuels sont :
#
# \\[
# Re=\\frac{UL}{\\nu},\\quad
# Pr=\\frac{\\nu}{\\alpha}=\\frac{\\mu c_p}{k},\\quad
# Gr=\\frac{g\\beta\\Delta T L^3}{\\nu^2},\\quad
# Ri=\\frac{Gr}{Re^2}=\\frac{g\\beta\\Delta T L}{U^2}.
# \\]
#
# \\(Ri\\) mesure le rapport des effets de flottabilité aux effets inertiels.
# Son signe dépend de l'orientation de la gravité, du gradient thermique et du
# sens de l'écoulement : il faut publier la convention de signe, pas seulement
# \\(|Ri|\\). Les régimes aiding/opposing peuvent modifier la turbulence,
# le frottement et le transfert thermique, voire provoquer une laminarisation
# partielle dans certaines conduites verticales chauffées.
#
# ## 4. Grandeurs pariétales et transfert thermique
#
# \\[
# u_\\tau=\\sqrt{\\frac{|\\tau_w|}{\\rho_w}},\\qquad
# y^+=\\frac{u_\\tau y}{\\nu_w},\\qquad
# U^+=\\frac{U}{u_\\tau}.
# \\]
#
# \\[
# C_f=\\frac{\\tau_w}{\\tfrac12\\rho U_{ref}^2},\\qquad
# Nu=\\frac{hL_c}{k_f},\\qquad
# h=\\frac{q_w''}{T_w-T_b}.
# \\]
#
# Le facteur de frottement de Darcy est \\(f_D=8\\tau_w/(\\rho U_b^2)\\)
# pour une conduite circulaire pleinement développée. Il ne faut pas confondre
# \\(C_f\\), le facteur de Fanning et le facteur de Darcy.
# Pour les écoulements séparés, le coefficient de pression est
# \\(C_p=(p-p_{ref})/(\\tfrac12\\rho U_{ref}^2)\\).
#
# ## 5. Statistiques et comparaison
#
# Pour une grandeur instantanée \\(q\\), la moyenne temporelle est
# \\(\\overline q=(1/T_s)\\int_{t_0}^{t_0+T_s}q(t)dt\\), et la fluctuation
# \\(q'=q-\\overline q\\). Les RMS et corrélations doivent être calculés
# après retrait du transitoire initial et sur une fenêtre statistiquement
# stationnaire. Les moyennes, RMS, spectres et contraintes ne sont pas
# interchangeables.
#
# Comparer les mêmes stations, coordonnées, normalisations, filtres et fenêtres
# d'échantillonnage que la référence. Rapporter l'erreur normalisée, le biais,
# l'incertitude expérimentale et la sensibilité à la fenêtre temporelle.
#
# ## 6. Vérification numérique avant validation physique
#
# Chaque campagne doit distinguer :
#
# - convergence itérative : résidus, bilans et stabilisation des QoI ;
# - erreur spatiale : au moins trois maillages systématiquement raffinés ;
# - erreur temporelle : au moins trois pas de temps pour les campagnes instationnaires ;
# - conservation : masse, quantité de mouvement et énergie selon le cas ;
# - incertitude de modèle : comparaison RANS/URANS/LES uniquement lorsque ces
#   formulations sont effectivement disponibles et correctement résolues.
#
# Un maillage grossier qui converge n'est pas une preuve de précision. Aucun
# seuil d'acceptation ne doit être choisi après observation des résultats.
#
# ## 7. Domaine d'applicabilité des modèles
#
# - k-omega SST : premier candidat RANS pour couches limites avec gradient de
#   pression adverse et séparation ; tester la sensibilité à y+ et au traitement
#   pariétal.
# - Realizable k-epsilon : comparaison utile pour le backward-facing step,
#   notamment avec traitement de paroi cohérent avec sa formulation.
# - RSM : nécessaire pour étudier explicitement les composantes anisotropes
#   des contraintes de Reynolds.
# - LES/URANS : requis pour prétendre comparer des fluctuations temporelles,
#   RMS thermiques ou spectres ; ne pas les déduire d'un RANS stationnaire.
#
# ## Références scientifiques principales
#
# - Moser, R. D., Kim, J. & Mansour, N. N. (1999), Physics of Fluids 11(4),
#   943–945, DOI: 10.1063/1.869966.
# - Driver, D. M. & Seegmiller, H. L. (1985), AIAA Journal 23, 163–171.
# - OECD/NEA (2011), Vattenfall T-junction benchmark, NEA/CSNI/R(2011)5.
# - Upflow turbulent mixed convection heat transfer in vertical pipes (1998),
#   International Journal of Heat and Mass Transfer 41, 4037–4054,
#   DOI: 10.1016/S0017-9310(98)00177-X.
# - Pope, S. B. (2000), Turbulent Flows, Cambridge University Press.
# - Wilcox, D. C. (2006), Turbulence Modeling for CFD, 3rd ed.
