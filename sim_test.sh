cat $2 > input.txt
echo "Simulate 3 10 10" >> input.txt
g++ -DCNTLOG -DLOG_SECTIONS $1 -o sol.out
./sol.out
rm sol.out
echo "Operations count:" `cat output.txt`
cat $2 > input.txt
echo "
Simulate  0.0 100000 100
Simulate  0.5 100000 100
Simulate  1.0 100000 100
Simulate  1.5 100000 100
Simulate  2.0 100000 100
Simulate  2.5 100000 100
Simulate  3.0 100000 100
Simulate  4.0 100000 100
--Simulate  5.0 100000 100
" >> input.txt
g++ -DTEST -DTIMELOG $1 -o sol.out
./sol.out
rm sol.out
while read -r line; do
    errr=`echo $line | grep -Eo "^\S+" | sed 's/\./,/'`
    echo $errr
done < output.txt
echo
while read -r line; do
    cnt=`echo $line | grep -Eo "\S+$" | sed 's/\./,/'`
    echo $cnt
done < output.txt