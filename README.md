# Projet CGAL : Octree et simplification

## Compilation

Pour compiler le programme, il faut se rendre dans le répertoire `build` et exécuter la commande `cmake ..` pour préparer la compilation. Ensuite, il faut compiler le programme avec la commande `make`.

## Utilisation

Le programme implémente la segmentation d'un fichier .OFF. Pour tester le programme avec un fichier d'entrée `.off`, on peut utiliser la commande suivante à partir du répertoire `build` :

```
./octree_base ../data/cube.off
```

On peut également ajouter deux paramètres optionnels :

Le premier paramètre est la profondeur max que peu prendre l'Octree. Par defaut elle est a `10`.
Le deuxième paramètre est le nombre maximum de sommet par feuilles. Par defaut elle est a `50`.

Voici un exemple d'utilisation avec les deux paramètres optionnels :

```
./color ../data/pig.off 5 20
```

Ce qui produit la sortie suivante:

```
Le résultat a été exporté dans colorMesh.off !
Le résultat a été exporté dans Octree.off !
Le résultat a été exporté dans simple.off !
```

Cette commande permet d'obtenir trois fichiers situés dans le dossier `build`. Le premier fichier est `colorMesh.off `qui contient le maillage initial avec une coloration des nœuds de l'octree. Le deuxième fichier est `Octree.off` qui représente la boîte englobante de l'octree. 


En plus des trois fichiers mentionnés précédemment, la commande produit également le fichier `Octree.json` qui représente l'arborescence de l'octree. Ce fichier est également situé dans le dossier `build`.
