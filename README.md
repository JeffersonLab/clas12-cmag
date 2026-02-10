# cmag
C++ version of the CLAS12 magnetic field package


# Installation (meson)

Pre-requisites: [meson](https://mesonbuild.com/index.html)

Assuming `/path/to/cmag` is the desired installation location:

```shell
meson setup build --prefix=/path/to/cmag
meson install -C build
```


# Installation (Makefile)

Pre-requisites: `make`


```shell
cd src
make
```