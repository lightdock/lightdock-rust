# LightDock-Rust

A Rust implementation of the [LightDock](https://lightdock.org) macromolecular software with the DFIRE scoring function.

## Installation
1. Clone this repository:

 ```
 git clone https://github.com/lightdock/lightdock-rust.git
 ```

2. Compile it with Rust (you may install Rust using [rustup](https://rustup.rs/)):

 ```
 cd lightdock-rust
 cargo build --release
 ```
 
## Examples

Several examples can be found in the `example` folder.

| :information_source: Data Path          |
|:---------------------------|
| You may set an environment variable `LIGHTDOCK_DATA` to point to the data folder included in this repository to avoid copying it: `export LIGHTDOCK_DATA=/path/to/lightdock-rust/data`  |

### 1k4c

```
cd example/1k4c
cp -R ../../data .
time ../../target/release/lightdock-rust setup.json initial_positions_0.dat 100 dfire

...

real    3m53.851s
user    3m52.550s
sys     0m0.717s
```

### 1ppe

```
cd example/1ppe
cp -R ../../data .
time ../../target/release/lightdock-rust setup.json initial_positions_0.dat 100 dfire

...

real    0m9.968s
user    0m9.640s
sys     0m0.271s
```

### 2uuy

```
cd example/2uuy
cp -R ../../data .
time ../../target/release/lightdock-rust setup.json initial_positions_0.dat 100 dfire

...

real    0m18.042s
user    0m17.123s
sys     0m0.621s
```

### 1czy

```
cd example/1czy
cp -R ../../data .
time ../../target/release/lightdock-rust setup.json initial_positions_0.dat 100 dfire

...

real    0m20.868s
user    0m19.692s
sys     0m1.070s
```

