#compilation

gcc -O2 -Wall -Wextra -pthread -o determinant_par determinant_par.c -lm

#Exemple pratique
# Version parallèle avec 4 threads
./determinant_par --n 9 --threads 4

# Avec 8 threads, 5 répétitions, vérification séquentielle
./determinant_par --n 10 --threads 8 --repeat 5 --verify

# Depuis un fichier
./determinant_par --input matrix.txt --threads 8 --verify

# Inspection des threads avec ps -L / top -H
./determinant_par --n 12 --threads 8 --hold 20
# dans un autre terminal :
#   ps -L -p <PID> | tail -n +2 | wc -l
#   top -H -p <PID>