# Ridgidbody Simulation

Ridgidbody simmulation project I made. The final product is essentially just from the papper:

```
KENNY ERLEBEN : Velocity-Based Shock Propagation for Multibody Dynamics Animation : ACM Transactions on Graphics, Vol. 26, 2007
```

It also comes with with a realtime vissualization app based on the [polyscope](https://polyscope.run/) library

## Compilation instructions

### Preraquisites

* Eigen
* Polyscope(gets downloaded automaticly through git submodules)
* Make
* Cmake

Get the code:
```
> git clone --recurse-submodules https://github.com/Dr-McFish/Ridgidbody-Simulation.git
```

Build:
```
> cd Ridgidbody-Simulation
> cmake -DCMAKE_BUILD_TYPE=DEBUG -S . -B build
> cd build
> make
```

Run the vissualiser app:
```
> cd ..
> ./build/bin/main.out
```

## Written material(in French)

[Slides](https://github.com/user-attachments/files/16102983/Anonymized_TIPE_slides-5.pdf)
[MCOT](https://github.com/user-attachments/files/16102995/Mcot_Anonymized.pdf)


### Future plans:
I would definitly like to experiment with an acceleration based moddel, and test more integration methods.



