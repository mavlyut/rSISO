if [[ -d img ]]; then
    rm -r img
fi
mkdir img

cat $2 > input.txt
echo "
Simulate 8 1000 1
" >> input.txt
g++ -DLOG -DVISUALIZE $1 src/*.cpp -o sol.out
./sol.out
rm sol.out
for file in `ls img`; do
    inputfile="img/$file"
    outputfile="img/${file/dot/svg}"
    dot -Tsvg $inputfile -o $outputfile
    rm $inputfile
done