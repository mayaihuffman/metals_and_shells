for pdb in `cat $1`
    do
	python atom_walk.py $pdb $2 $3
    done
