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

### 1k4c

```
cd example/1k4c
cp -R ../../data .
time ../../target/release/lightdock-rust setup.json initial_positions_0.dat 100 dfire

...

real    4m27.412s
user    4m21.112s
sys     0m1.595s
```

### 1ppe

```
cd example/1ppe
cp -R ../../data .
time ../../target/release/lightdock-rust setup_1ppe.json initial_positions_0.dat 100 dfire

...

real    0m11.923s
user    0m11.151s
sys     0m0.484s
```

### 2uuy

```
cd example/2uuy
cp -R ../../data .
time ../../target/release/lightdock-rust setup.json initial_positions_0.dat 100 dfire

...

real    0m20.868s
user    0m19.692s
sys     0m1.070s
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

