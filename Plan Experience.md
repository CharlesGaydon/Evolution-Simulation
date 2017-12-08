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

- pIns = 1 : grand potentiel.
- pDel = 1 : vite limité
- pInv = 1¨: possibilité

- pIns=pDel=0.5 : ordre des gènes ne change pas ; seul leur espacement relatif peut varier.

- pIns=pDel=pInv=1/3

- design d'expérimentation sur un triangle.


