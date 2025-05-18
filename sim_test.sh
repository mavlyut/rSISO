cat $2 > input.txt
echo "
Simulate 0.0 100000 100
Simulate 0.5 100000 100
Simulate 1.0 100000 100
Simulate 1.5 100000 100
Simulate 2.0 100000 100
Simulate 2.5 100000 100
Simulate 3.0 100000 100
Simulate 3.5 100000 100
Simulate 4.0 100000 100
Simulate 4.5 100000 100
Simulate 5.0 100000 100
" >> input.txt
g++ -DTEST -DTIMELOG $1 src/*.cpp -o sol.out
./sol.out
rm sol.out
while read -r line; do
    cnt=`echo $line | sed 's/\./,/'`
    echo $cnt
done < output.txt