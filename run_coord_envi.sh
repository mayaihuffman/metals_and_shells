for pdb in `cat $1`
    do
	python coord_envi.py $pdb $2 $3
    done
