
#Afivo results
#set key autotitle columnhead
set datafile separator ','
set logscale y 10
set yrange [10:10e21]
plot "Afivo_data.dat" u ($1/1e-9):2 with lines lt 1 lc 1 title 'e',\
     "Afivo_data.dat" u ($1/1e-9):5 with lines lt 1 lc 2 title 'N2+',\
     "Afivo_data.dat" u ($1/1e-9):6 with lines lt 1 lc 3 title 'O',\
     "Afivo_data.dat" u ($1/1e-9):7 with lines lt 1 lc 4 title 'O2-',\
     "Afivo_data.dat" u ($1/1e-9):8 with lines lt 1 lc 5 title 'O2+',\
     "Afivo_data.dat" u ($1/1e-9):9 with lines lt 1 lc 6 title 'O-'
     
#ZDPlaskin results
set datafile separator ' '
plot "ZDPlaskin_simplerAir.dat" u ($1/1e-9):2 with points lt 2 lc 1 title 'e',\
     "ZDPlaskin_simplerAir.dat" u ($1/1e-9):5 with points lt 2 lc 2 title 'N2+',\
     "ZDPlaskin_simplerAir.dat" u ($1/1e-9):9 with points lt 2 lc 3 title 'O',\
     "ZDPlaskin_simplerAir.dat" u ($1/1e-9):7 with points lt 2 lc 4 title 'O2-',\
     "ZDPlaskin_simplerAir.dat" u ($1/1e-9):6 with points lt 2 lc 5 title 'O2+',\
     "ZDPlaskin_simplerAir.dat" u ($1/1e-9):8 with points lt 2 lc 6 title 'O-'
