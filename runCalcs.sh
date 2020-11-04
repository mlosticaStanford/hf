python writeIn.py;

for i in $(seq 0.1 0.1 10.0)
  do 
  mv $i.in input.txt;
  ./main > text.txt;
  mv output.txt outputFiles/$i.out;
  mv text.txt outputFiles/$i.info;
  done

cd outputFiles;
cp ../writeOut.py ./;
python writeOut.py;
cd ..;

