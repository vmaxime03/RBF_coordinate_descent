```shell
make build CXX=g++ CC=gcc

# run the optimization
make run

# requires python3, numpy, matplotlib, pandas 
make anim           # show coordinate descent animation
make gif            # generate the gif in ./output/anim.gif 
make correlation    # print correlation matrix

# requires octave 
make octave         # optimization result plot
```
