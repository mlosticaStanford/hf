#make rhf curve

for i in $(seq 0.01 0.01 10.0)
  do  
 #awk -v dist=$i '{if (NR==4) print $1,$2,$3,$dist; else print $0}' input > r_"$i".in;
  sed "s/dist/$i/" inTemp.txt > input.txt;
  ./main ;
  mv output.txt energies/"$i".out ;
  cd ./energies;
  awk '{if ($1=="Total" && $2=="ground") print $6}' "$i".out >> energy.txt;
  cd ../ ;
  done
