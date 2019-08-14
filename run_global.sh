for pdb in `cat $1`
    do
	python metal_site_global.py $pdb $2 $3
    done
