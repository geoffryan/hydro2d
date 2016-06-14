#N="0008 0016 0032 0064 0128 0256 0512 1024 2048 4096 8192"
N="0008 0016 0032 0064 0128 0256 0512 1024 2048 4096"
parfile="parfiles/isentrope.par"
name=${parfile##*/}
mkdir -p temp

echo $name
cp $parfile temp/$name

for n in $N
do
    filen=temp/$name.$n
    cp -f temp/$name $filen
    sed -e "s/\(Nx1\)[ \t]*\([0-9]*\)\(.*\)/\1 $n \3/" $filen > $filen.temp
    sed -e "s/\(NumCheckpoints\)[ \t]*\([0-9]*\)\(.*\)/\1 1 \3/" $filen.temp > $filen 
    rm $filen.temp
    echo "Running Nx=$n"
    bin/hydro2d $filen > temp/out.$n
    cp "checkpoint_00000.h5" "temp/chk0.$n.h5"
    cp "checkpoint_00001.h5" "temp/chk1.$n.h5"
done
