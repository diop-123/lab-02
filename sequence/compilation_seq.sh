#compilatiom
gcc -O2 -Wall -Wextra -o determinant_seq determinant_seq.c -lm

#Exemple practique
# Matrice 6x6 aléatoire (graine fixe = 42)
./determinant_seq --n 6

# Affichage de la matrice + 3 répétitions (on garde le meilleur temps)
./determinant_seq --n 7 --verbose --repeat 3

# Depuis un fichier
./determinant_seq --input matrix.txt