#90	90	1.0	1.0	1.0
#1.0 0.05
#1
#3 F 1e-3
#3 0.001 1000

#1ere ligne  ->  paramètres géométriques: Nx, Ny, Lx, Ly, D
#2eme ligne  ->  paramètres temporels : Tmax, Deltat
#3eme ligne  ->  numéro du cas test à résoudre : nb_probleme (1: stationnaire, 2:idem, 3: instationnaire)
#4eme ligne  ->  décomposition : recouvrement et choix de l'algo : T pour l'algo multiplicatif, F pour l'algo additif, tolérance des algo de Schwarz
#5eme ligne  ->  paramètres solveur (Jocobi=1 / GaussSeidel=2 / GradConj=3), tolérance, itérations max


#To do : 

#+ Comparer le temps d'execution et/ou la précision en faisant varier le recouvrement pour un algo donné. Entre 0 et 4
#(1er arg ligne 4)
#+ Comparer Additif et Multiplicatif pour un meme recouvrement 
#+ Speed Up pour un recouvrement donné


#Comparaison Additif, Multiplicatif
echo "Debut comparaison additif/multiplicatif"
for i in {T,F}
do
    echo "Comparaison $i"
    printf "" > test_comp_$i.dat

    for j in {0..4}
    do
	echo "Recouvrement $j"
	(head -3 "parametres"; echo "$j $i 1e-3"; tail -n +5 "parametres") > foobar && mv foobar parametres;
	printf "$j " >> test_comp_$i.dat
	mpirun -n 4 ./run >> test_comp_$i.dat
    done
done

#Speedup


echo "Debut courbe de speedup"
rec='2'
for i in {T,F}
do
    echo "Type $i"
    (head -3 "parametres"; echo "$rec $i 1e-3"; tail -n +5 "parametres") > foobar && mv foobar parametres;
    #ICI FAUT METTRE LE TEMPS DE LALGO SEQUENTIEL
    printf "" > test_seq_$i.dat
    ./run_seq >> test_seq_$i.dat

    printf "" > test_speedup_$i.dat
    for j in {1..24}
    do
	echo "nb threads: $j"
	printf "$j " >> test_speedup_$i.dat
	mpirun -n $j ./run >> test_speedup_$i.dat
    done
done
