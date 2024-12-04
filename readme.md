# TP Assemblage HAU902I

Dans le cadre de l'UE HAU902I, ce projet développe un assembleur très simple en python.



## Execution
Nécessite python 3.10. Lancer le main en ligne de commande avec deux arguments dans l'ordre :
1) path : chemin vers le fichier de read (format fastq)
2) k : longueur des kmers

Les contigs générés sont écrits dans un fichier "contigs_<k>.fa"


## Exemple

```console
python TP_assemblage/main.py reads.fastq.fq 13
```
Sortie :
```console
Assembly from file reads.fastq.fq with k=13.
Executed in 1.4860658249999688s.
```


## Auteurs :

Margaux Imbert et Florent Marchal


