cat $1 > input.txt
echo "
Simulate  0.0 100000 100
Simulate  0.5 100000 100
Simulate  1.0 100000 100
Simulate  1.5 100000 100
Simulate  2.0 100000 100
Simulate  2.5 100000 100
Simulate  3.0 100000 100
Simulate  4.0 100000 100
Simulate  5.0 100000 100" >> input.txt
g++ main.cpp -o sol.out
./sol.out
rm sol.out
p=""
while read -r line; do
    if [[ -z $p ]]; then
        p="1"
        continue
    else
        errr=`echo $line | grep -Eo "^\S+" | sed 's/\./,/'`
        echo $errr
    fi
done < output.txt
echo
p=""
while read -r line; do
    if [[ -z $p ]]; then
        p="1"
        continue
    else
        time=`echo $line | grep -Eo " \S+ " | sed -e 's/\./,/; s/ //'`
        echo $time
    fi
done < output.txt
echo
p=""
while read -r line; do
    if [[ -z $p ]]; then
        p="1"
        continue
    else
        cnt=`echo $line | grep -Eo "\S+$" | sed 's/\./,/'`
        echo $cnt
    fi
done < output.txt