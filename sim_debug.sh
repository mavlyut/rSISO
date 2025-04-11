cat $2 > input.txt
echo "
Simulate 8 10000 100
" >> input.txt
g++ -DLOG -DTIMELOG $1 src/*.cpp -o sol.out
./sol.out
rm sol.out