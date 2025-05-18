cat $2 > input.txt
echo "
Simulate 100 1000 10
" >> input.txt
g++ -DLOG -DTIMELOG -DPRINT_SECTION_TREE $1 src/*.cpp -o sol.out
./sol.out
rm sol.out