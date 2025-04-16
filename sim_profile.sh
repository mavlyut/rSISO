cat $2 > input.txt
echo "
Simulate 0.5 100000 100
Simulate 1.0 100000 100
Simulate 1.5 100000 100
Simulate 2.0 100000 100
Simulate 2.5 100000 100
--Simulate 3.0 100000 100
--Simulate 3.5 100000 100
--Simulate 4.0 100000 100
--Simulate 4.5 100000 100
" >> input.txt
g++ -DTEST -DTIMELOG -pg $1 src/*.cpp -o sol.out
./sol.out
cat gmon.out | gprof ./sol.out > _log_profiler.txt && rm gmon.out && rm sol.out
gprof2dot ./_log_profiler.txt | dot -Tsvg -o output.svg