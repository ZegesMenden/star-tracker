# Star Tracking

## Pipeline

The general process for star tracking is:

1. Centroid identification — finding the stars in the image
2. Star ID classification — identifying the centroids and matching them to known stars
3. Attitude estimation — turning the known map of stars in frame to a rotation of the camera

This implementation uses:

* Center-of-gravity algorithm for centroid identification
* Pyramid algorithm for star classification
* QUEST algorithm for attitude estimation

---

## Building

The project uses CMake >= 3.14 and requires a C++17 compiler.

### Basic Build

```bash
mkdir build && cd build
cmake ..
make
```

### CMake Options

| Option                      | Description                                  | Default |
|-----------------------------|----------------------------------------------|---------|
| `STARTRACKER_USE_FLOAT`     | Use `float` instead of `double`              | `ON`    |
| `STARTRACKER_BUILD_TESTS`   | Build the unit-test executable               | `ON`    |
| `STARTRACKER_BUILD_TOOLS`   | Build helper tools (generate-db, benchmark)  | `ON`    |

Options are passed with `-D` at configure time:

```bash
# Release build, double precision, no tools
cmake -DCMAKE_BUILD_TYPE=Release \
      -DSTARTRACKER_USE_FLOAT=OFF \
      -DSTARTRACKER_BUILD_TOOLS=OFF \
      ..
```

### Common CMake Configurations

```bash
# Debug build (default) — all features enabled
cmake ..

# Release build with optimizations
cmake -DCMAKE_BUILD_TYPE=Release ..

# Build only the library (no tests, no tools)
cmake -DSTARTRACKER_BUILD_TESTS=OFF -DSTARTRACKER_BUILD_TOOLS=OFF ..

# Double precision build
cmake -DSTARTRACKER_USE_FLOAT=OFF ..
```

### Build Targets

```bash
# Build everything
make

# Build only the library
make startracker

# Build only the test executable
make startracker_tests

# Build individual tools
make generate-db
make generate-stars
make generate-catalog
make benchmark
```

### Install

```bash
cmake --build build --target install
# or
cd build && make install
```

This installs headers to `<prefix>/include/startracker/` and the library to the appropriate lib directory. A CMake package config is also installed so downstream projects can use:

```cmake
find_package(startracker REQUIRED)
target_link_libraries(myapp PRIVATE startracker::startracker)
```

---

## Testing

### Unit Tests (Catch2)

The project uses [Catch2](https://github.com/catchorg/Catch2) (single-header, vendored in `vendor/`) for unit testing.

**Run all tests via CTest:**

```bash
cd build
ctest
```

Or with verbose output:

```bash
ctest --verbose
# or
ctest -V
```

**Run the test executable directly** for more control:

```bash
./startracker_tests
```

**Catch2 CLI options** (passed directly to the test executable):

```bash
# List all test cases
./startracker_tests --list-tests

# Run a specific test case by name
./startracker_tests "test case name"

# Run tests matching a tag
./startracker_tests [tag]

# Run with verbose output
./startracker_tests -s

# Run a specific test file's tests (by tag or substring)
./startracker_tests "attitude*"
./startracker_tests "camera*"
```

**Test source files:**

| File                  | Coverage                |
|-----------------------|-------------------------|
| `test/test_attitude.cc` | Attitude estimation   |
| `test/test_camera.cc`   | Camera model          |
| `test/test_centroid.cc`  | Centroid detection    |
| `test/test_database.cc`  | Star database         |
| `test/test_serialize.cc` | Serialization         |
| `test/test_star.cc`      | Star utilities        |
| `test/test_starid.cc`    | Star identification   |
| `test/test_vec.cc`       | Vector math           |

### Benchmark

The `benchmark` tool runs the full star-tracker pipeline on a single image and prints the resulting attitude for comparison against a known solution.

```bash
./build/benchmark \
    --image <path>           \
    --database <db-file>     \
    --catalog <bsc.tsv>      \
    --focal-length <px>
```

Or using physical camera parameters:

```bash
./build/benchmark \
    --image <path>           \
    --database <db-file>     \
    --catalog <bsc.tsv>      \
    --pixel-size <µm>        \
    --focal-length-mm <mm>
```

**Output:** `RA_deg DEC_deg ROLL_deg NUM_IDENTIFIED NUM_CENTROIDS`

### Automated Benchmark Script

The `benchmark/generate-and-test.sh` script generates synthetic star field images at random orientations and tests the pipeline against known ground-truth attitudes.

```bash
./benchmark/generate-and-test.sh [OPTIONS]
```

**Options:**

| Flag               | Description                              | Default    |
|--------------------|------------------------------------------|------------|
| `-n <count>`       | Number of images to generate             | `10`       |
| `-c <catalog>`     | Path to bright-star-catalog TSV          | auto-detect|
| `-d <database>`    | Path to serialized MultiDatabase file    | auto-detect|
| `-b <benchmark>`   | Path to benchmark executable             | `build/benchmark` |
| `-g <generate>`    | Path to generate-stars executable        | `build/generate-stars` |
| `-f <fov-deg>`     | Field of view in degrees                 | `15`       |
| `-t <tol-deg>`     | Angular tolerance for pass/fail          | `1.0`      |
| `--width <px>`     | Image width                              | `512`      |
| `--height <px>`    | Image height                             | `512`      |
| `--spread <σ>`     | Star spread std-dev                      | `2.0`      |
| `--photons <N>`    | Zero-mag total photons                   | `200000`   |
| `--saturation <N>` | Saturation photons                       | `1000`     |
| `--dark-current <f>` | Dark current (0–1)                     | `0.0`      |
| `--read-noise <f>` | Read noise std-dev                       | `0.01`     |
| `--shot-noise`     | Enable shot noise                        | off        |
| `--seed <int>`     | Random seed                              | time-based |

---

## Tools

Built when `STARTRACKER_BUILD_TOOLS=ON` (default).

| Executable         | Description                                              |
|--------------------|----------------------------------------------------------|
| `generate-db`      | Build a serialized star database from a catalog          |
| `generate-stars`   | Generate synthetic star field images                     |
| `generate-catalog` | Generate/filter a star catalog                           |
| `benchmark`        | Run the full pipeline on a real or synthetic image       |