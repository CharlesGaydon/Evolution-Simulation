# Plan d'expérience

## Hypothèses

- les gènes restent intouchés par les mutations.
- relire le sujet.
- Unité d'indel < taille maximale d'un gène

## Critère 
Macro :
- ratio fitness de fin / fitness de début à temps de simulation égal
- (fitness de début - fitness de fin) / temps de simulation à temps de simulation non égal

Micro :
- taille génome
- espacement moyen/médiane/variance de distance inter-gènes
- nombre de ++, --, +- (succession de gènes)

Figures :

- Une exploration de l'espace des probabilités retranscrite sur un triangle de mélange en deux dimension par une heatmap. 

### Tousgenesidentiques

#### Plans
- TASK1 : pIns = 1 : potentiellement intéressant.
- TASK2 : pDel = 1 : probablement vite limité
- TASK3 : pInv = 1¨: le meilleur potentiel.

- INDEL1 : pIns=pDel=0.5 : ordre et orientation des gènes ne change pas ; seul leur espacement relatif peut varier, ce qui jouera sur leur interactions de manière limitée.

- HOMOG1 : pIns=pDel=pInv=1/3

- design d'expérimentation sur un triangle.

#### Variance de la mesure de fitness

Avec quelle précision la simulation nous donne-t-elle la fitness ?
Sur cent répétitions nous obtenons : sd = 0.000749.
Ramenée à la valeur initiale de la fitness, de moyenne sur cent répétition  mean = -0.069414, on trouve une erreur d'environ 1.1%. Evidemment, la dispersion des simulation peut être dépendante de la structure du génome, mais en première approximation on la considère comme négligeable face aux variations induites par les mutations.
L'intervalle de crédibilité autour de chaque valeur mesurée correspond à un écart-type de part et d'autre. Dans le cas où la nouvelle fitness est inférieure à la précédente, si l'écart delta entre elles est inférieur à 2 x sigma, la différence peut être vue comme due à la dispersion car les intervalles se chevauchent. On choisit de faire en sorte que pour un delta négatif = 4 x sigma, 90% des valeurs soient refusées. Cela nous conduit à choisir : p = exp(-720 x delta)
