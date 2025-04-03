After the make has been compiled by `make` , run program with the following command (using huge.fa as an example)

`mpirun -np 8 . /cdhit -i huge.fa -o test -T 4 -M 60000`

For this -M parameter, for huge probably about 1200 is enough, this place 60000 is what I wrote blindly. The value x written here means xM.

If there are other parameters you can keep them the same as cdhit.

Working Logs:

QUESTIONS:

有些地方的代码会在某次运行时会报错，但不是经常出现。