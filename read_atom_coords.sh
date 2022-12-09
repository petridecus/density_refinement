cat ../Desktop/5nrodel.pdb | grep ATOM | cut -c 33- | sed 's/^/[/' | sed 's/./,/8' | sed 's/./,/15' | sed 's/./], #/24' | sed 's/^/                /' > array_vals.txt
